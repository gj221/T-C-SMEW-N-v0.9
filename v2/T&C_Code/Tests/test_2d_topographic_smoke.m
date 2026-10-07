function test_2d_topographic_smoke(n_days)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% test_2d_topographic_smoke
%%
%% Integration smoke test of the topographic 2-D driver: build a small
%% multi-cell topology on a synthetic tilted plane, attach it, and run
%% MAIN_FRAME_BG_ERW_2D over a few days.  Checks that every new topographic
%% code path executes and that the domain water / solute ledgers close.
%%
%% Reuses the forcing bundled in the equivalence baseline
%%   v2/SIM/EQUIV_1D_TMP/TC_ERW_results_Kantola_crop_control.mat
%% (generate it as test_2d_driver_equivalence documents).
%%
%% Usage
%%   test_2d_topographic_smoke        % 5 simulated days
%%   test_2d_topographic_smoke(10)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 1 || isempty(n_days); n_days = 5; end
fprintf('\n=== test_2d_topographic_smoke (%d days) ===\n', n_days);

here = fileparts(mfilename('fullpath'));
v2   = fileparts(fileparts(here));
code = fullfile(v2,'T&C_Code');
addpath(code, fullfile(code,'Inputs'), fullfile(code,'Parameters'), here);

base_file = fullfile(v2,'SIM','EQUIV_1D_TMP', ...
    'TC_ERW_results_Kantola_crop_control.mat');
if ~isfile(base_file)
    error('test_2d_topographic_smoke:NoForcing', ...
        ['Forcing bundle not found: %s\nGenerate it first with:\n' ...
         '  TC_CASE=5 TC_SCENARIO=EQUIV_1D_TMP TC_NDAYS=%d matlab -batch ' ...
         '"cd(''v2/run''); ERW_run"'], base_file, max(n_days,10));
end
R1 = load(base_file);
n_hours = n_days*24;
if R1.NN < n_hours
    error('test_2d_topographic_smoke:ShortForcing', ...
        'Forcing bundle has %d h; need %d.', R1.NN, n_hours);
end

n_pass = 0; n_fail = 0;
    function check(name,cond,detail)
        if nargin < 3; detail = ''; end
        if cond; n_pass = n_pass+1; fprintf('  PASS  %s\n',name);
        else;     n_fail = n_fail+1; fprintf('  FAIL  %s  %s\n',name,detail); end
    end

%% --- topology: 3x3 main grid + buffer on a 6% tilted plane -----------
topo = Build_Cell_Topology(60,5,3,3);
z = 100 + 0.06*(max(topo.cell_yc) - topo.cell_yc) ...      % slope toward -y
      + 0.5*sin(topo.cell_xc/15);                          % gentle cross-roll
topo = Attach_Topography_To_Topology(topo, struct('cell_z_m',z));
check('topography attaches to the sub-field graph', ...
    isfield(topo,'topography_attached') && topo.topography_attached);
check('a downslope outlet is identified', any(topo.cell_is_outlet));
check('edge elevation drops are non-trivial', max(abs(topo.edge_dz_m)) > 0.1);

[Dom, runoff_total_m3] = run_topographic_driver(R1, code, topo, n_hours, 1);

%% --- ledger + state checks -----------------------------------------
check('driver completes without error', isstruct(Dom));
check('surface solute export vector is finite and sized right', ...
    numel(Dom.exit_mass_surface_g) == numel(Dom.exit_mass_surface_species) && ...
    all(isfinite(Dom.exit_mass_surface_g)));

%% Subsurface lateral flow must have redistributed water on the slope.
lat_water = sum(abs(Dom.lateral.water_buf_to_main),'omitnan');
check('subsurface lateral flow redistributes water down the slope', ...
    lat_water > 0, sprintf('|net buffer<->main| = %.4g m3', lat_water));
check('a dissolved species is carried with that lateral water', ...
    sum(abs(Dom.lateral.NO3_buf_to_main),'omitnan') > 0);
check('subsurface lateral sub-steps stayed bounded', ...
    max(Dom.lateral.sub_steps_max) < 5000);

%% Surface routing: assert conservation whenever runoff was generated.
if runoff_total_m3 > 1e-9
    surf_out = sum(Dom.lateral.surface_water_out_m3);
    surf_reinf = sum(Dom.lateral.surface_reinfilt_m3);
    check('overland runoff is routed (exit + re-infiltration = runoff)', ...
        abs((surf_out + surf_reinf) - runoff_total_m3) < 1e-6*max(1,runoff_total_m3), ...
        sprintf('out %.4g + reinf %.4g vs runoff %.4g', ...
        surf_out, surf_reinf, runoff_total_m3));
    check('the topographic outlet exported some water', surf_out > 0);
else
    fprintf('  SKIP  surface routing checks (no runoff in this %d-day window)\n', ...
        n_hours/24);
end

check('total domain water export is non-negative and finite', ...
    all(isfinite(Dom.exit_water_total_m3)) && all(Dom.exit_water_total_m3 >= -1e-9));
check('no NaN in the area-weighted pore pH series', ...
    all(isfinite(Dom.PH_pore_d(Dom.PH_pore_d ~= 0))));

%% --- WET pass: force overland runoff so the surface path is covered ---
%% The bundled window is dry, so without amplified rain the router, the
%% ponded-solute store and the ERW overland export never execute.
%% One day is enough to drive the routing, the ponded-solute store and
%% the ERW overland export.  Amplifying the whole 5-day window instead
%% keeps the column saturated, which sends the adaptive Darcy sub-stepping
%% into a grind, and 20x rainfall is not a physical test in any case.
%% The bundled window can be dry for days, so scaling its rain proves
%% nothing; inject a storm instead -- 25 mm/h for four hours (100 mm), a
%% heavy but physical Virginia convective event -- which exceeds any
%% loam's infiltration capacity and forces overland flow on the slope.
wet_hours = min(n_hours,24);
fprintf('  --- wet pass (25 mm/h storm, hours 6-9 of %d) ---\n', wet_hours);
[DomW, runoffW_m3] = run_topographic_driver(R1, code, topo, wet_hours, 1, 25);
check('the wet pass actually generates overland runoff', runoffW_m3 > 1e-6, ...
    sprintf('runoff %.4g m3', runoffW_m3));
surfW_out = sum(DomW.lateral.surface_water_out_m3);
surfW_reinf = sum(DomW.lateral.surface_reinfilt_m3);
check('routed runoff is conserved (outlet export + re-infiltration = runoff)', ...
    abs((surfW_out + surfW_reinf) - runoffW_m3) < 1e-6*max(1,runoffW_m3), ...
    sprintf('out %.6g + reinf %.6g vs runoff %.6g', ...
    surfW_out, surfW_reinf, runoffW_m3));
check('the topographic outlet exported water', surfW_out > 0, ...
    sprintf('exported %.4g m3',surfW_out));
check('solute left with it, and only where water did', ...
    all(isfinite(DomW.exit_mass_surface_g)) && all(DomW.exit_mass_surface_g >= -1e-9) && ...
    any(DomW.exit_mass_surface_g > 0));
check('the field-origin export never exceeds the total export', ...
    DomW.exit_water_field_surface_m3 <= sum(DomW.exit_water_surface_m3) + 1e-9 && ...
    DomW.exit_water_field_surface_m3 >= -1e-9, ...
    sprintf('field %.6g of %.6g m3', DomW.exit_water_field_surface_m3, ...
    sum(DomW.exit_water_surface_m3)));
check('field-origin solute export is bounded by the total likewise', ...
    all(DomW.exit_mass_surface_field_g <= DomW.exit_mass_surface_g + 1e-9) && ...
    all(DomW.exit_mass_surface_field_g >= -1e-9));
check('ERW aqueous species are exported overland too (Na, Al, An, DIC)', ...
    numel(DomW.exit_mass_surface_erw_umol) == 4 && ...
    all(isfinite(DomW.exit_mass_surface_erw_umol)) && ...
    all(DomW.exit_mass_surface_erw_umol >= -1e-6), ...
    sprintf('%s umol', mat2str(DomW.exit_mass_surface_erw_umol.',3)));
check('DIC leaves with the exported water, so alkalinity export is not lost', ...
    DomW.exit_mass_surface_erw_umol(4) > 0, ...
    sprintf('DIC export %.6g umol', DomW.exit_mass_surface_erw_umol(4)));
check('the wet pass keeps pore pH finite', ...
    all(isfinite(DomW.PH_pore_d(DomW.PH_pore_d ~= 0))));

%% --- this tree's switch: OPT_DEM and the topology must agree ----------
flat = Build_Cell_Topology(60,5,3,3);
check('elevation in the topology with OPT_DEM = 0 is refused', ...
    raises(@() run_topographic_driver(R1,code,topo,24,1,0,0), ...
    'MAIN_FRAME_BG_ERW_2D:TopographyBehindSwitch'));
check('OPT_DEM = 1 without topography is refused', ...
    raises(@() run_topographic_driver(R1,code,flat,24,1,0,1), ...
    'MAIN_FRAME_BG_ERW_2D:NoTopographyAttached'));
check('a switch value other than 0 or 1 is refused', ...
    raises(@() run_topographic_driver(R1,code,flat,24,1,0,2), ...
    'MAIN_FRAME_BG_ERW_2D:BadSwitch'));

%% --- a polygon grid (v1.0's majority rule) under a storm --------------
%% Built from polyshape alone, as a saved field grid is loaded on the
%% cluster, with its outer ring of squares labelled margin. Partial
%% perimeter squares make the polygon weights and the field-origin
%% ledgers differ from the plain totals, which is what this checks.
field = polyshape([0 31 34 25 4 -2],[0 -1 14 24 22 9]);
[grid_polygon,~,offset] = Place_Field_Polygon_Grid(field,5);
ptopo = Build_Polygon_Grid_Topology(grid_polygon,5,'majority',offset);
ptopo.cell_region(ptopo.outer_face_width_m > 0) = {'buffer'};
zp = 100 + 0.05*(max(ptopo.cell_yc) - ptopo.cell_yc);
ptopo = Attach_Topography_To_Topology(ptopo, struct('cell_z_m',zp));
fprintf('  --- polygon grid: %d cells of 5 m, storm as above ---\n', ptopo.N_cells);
[DomP, runoffP_m3] = run_topographic_driver(R1, code, ptopo, wet_hours, 1, 25);
check('the polygon grid weighs the domain by its overlap areas', ...
    strcmp(DomP.field_aggregation,'polygon_overlap_area') && ...
    isequal(DomP.cell_overlap_area_m2,ptopo.cell_overlap_area_m2) && ...
    abs(DomP.area_m2-area(field)) < 1e-6*area(field), ...
    sprintf('domain %.6g m2 for a %.6g m2 field', DomP.area_m2, area(field)));
check('main and margin areas add up to the captured polygon', ...
    abs(DomP.area_main_m2+DomP.area_buffer_m2-sum(ptopo.cell_overlap_area_m2)) < 1e-9);
surfP_out = sum(DomP.lateral.surface_water_out_m3);
surfP_reinf = sum(DomP.lateral.surface_reinfilt_m3);
check('routed runoff on the polygon grid is conserved', ...
    abs((surfP_out + surfP_reinf) - runoffP_m3) < 1e-6*max(1,runoffP_m3), ...
    sprintf('out %.6g + reinf %.6g vs runoff %.6g', surfP_out, surfP_reinf, runoffP_m3));
check('field-origin exports are bounded by the totals', ...
    DomP.exit_water_field_surface_m3 >= -1e-9 && ...
    DomP.exit_water_field_surface_m3 <= sum(DomP.exit_water_surface_m3)+1e-9 && ...
    DomP.exit_water_field_vert_m3 <= sum(DomP.exit_water_vert_m3)+1e-9 && ...
    all(DomP.exit_mass_surface_erw_field_umol <= DomP.exit_mass_surface_erw_umol+1e-6));
check('the polygon grid keeps pore pH finite', ...
    all(isfinite(DomP.PH_pore_d(DomP.PH_pore_d ~= 0))));

fprintf('  ---------------------------------------\n');
fprintf('  %d passed, %d failed\n\n', n_pass, n_fail);
if n_fail > 0
    error('test_2d_topographic_smoke: %d check(s) failed', n_fail);
end
end

%% =======================================================================
function [Domain, runoff_total_m3] = run_topographic_driver(R1, code, topo, n_hours, rain_factor, storm_mm_h, opt_dem)
if nargin < 6; storm_mm_h = 0; end
if nargin < 7; opt_dem = 1; end
%% This tree gates topography on an explicit switch (OPT_DEM), which the
%% frame checks against the topology; per-cell soils stay off.
OPT_DEM = opt_dem;                                            %#ok<NASGU>
OPT_CELL_INPUTS = 0;                                          %#ok<NASGU>
%% Rebuild the caller context MAIN_FRAME_BG_ERW_2D expects (mirrors
%% test_2d_driver_equivalence/run_driver_single) but for N_col > 1 with a
%% topography-attached topology, and truncate the forcing to n_hours.

NN = n_hours; NNd = ceil(NN/24)+1;                             %#ok<NASGU>
ms = R1.ms; cc = R1.cc; Zs = R1.Zs; dt = R1.dt; dth = R1.dth;  %#ok<NASGU>
%% The bundled series are [1 x NN_full] rows; Datam / IrD are [NN_full x k].
cr = @(a) a(:, 1:min(size(a,2),n_hours));
cm = @(a) a(1:min(size(a,1),n_hours), :);
Datam = cm(R1.Datam); Date = cr(R1.Date);                     %#ok<NASGU>
Pr=cr(R1.Pr)*rain_factor; Ta=cr(R1.Ta); Ds=cr(R1.Ds); Ws=cr(R1.Ws); %#ok<NASGU>
if storm_mm_h > 0
    hrs = 6:min(9,numel(Pr));
    Pr(hrs) = Pr(hrs) + storm_mm_h;          % [mm/h] synthetic storm
    Ta(hrs) = max(Ta(hrs),4);                % keep it rain, not snow
end
ea=cr(R1.ea); N=cr(R1.N); Pre=cr(R1.Pre); Tdew=cr(R1.Tdew);   %#ok<NASGU>
SAB1=cr(R1.SAB1); SAB2=cr(R1.SAB2); SAD1=cr(R1.SAD1); SAD2=cr(R1.SAD2); %#ok<NASGU>
PARB=cr(R1.PARB); PARD=cr(R1.PARD); Ca=cr(R1.Ca); Oa=R1.Oa;   %#ok<NASGU>
Lat=R1.Lat; Lon=R1.Lon; DeltaGMT=R1.DeltaGMT;                 %#ok<NASGU>
t_bef=R1.t_bef; t_aft=R1.t_aft; id_location=R1.id_location;   %#ok<NASGU>
NN_case_study_ERW_MOD = R1.NN_case_study_ERW_MOD;             %#ok<NASGU>
MOD_PARAM_input = R1.MOD_PARAM_input;
ERW_PARAM_input = R1.ERW_PARAM_input;
OPT_ERW = R1.OPT_ERW;                                         %#ok<NASGU>
OPT_SoilBiogeochemistry = R1.OPT_SoilBiogeochemistry;         %#ok<NASGU>
if isfield(R1,'B_out');     B_out = R1.B_out;         end     %#ok<NASGU>
if isfield(R1,'O_out');     O_out = R1.O_out;         end     %#ok<NASGU>
if isfield(R1,'NavlI_out'); NavlI_out = R1.NavlI_out; end     %#ok<NASGU>
if isfield(R1,'RexmyI_out');RexmyI_out = R1.RexmyI_out;end    %#ok<NASGU>
if isfield(R1,'NupI_L_out');NupI_L_out = R1.NupI_L_out;end    %#ok<NASGU>
if isfield(R1,'IrD');       IrD = cm(R1.IrD);         end     %#ok<NASGU>

Parameter_Directory = fullfile(code,'Parameters');            %#ok<NASGU>
apply_rock = logical(R1.OPT_ERW && ERW_PARAM_input.diss_f > 0);
%% The column settings v0.9's equivalence test gives the frame: the field
%% record's last application, no replaced literals, native saturation.
ColCfg = repmat(struct( ...
    'PARAM_IC', fullfile(Parameter_Directory, [char(MOD_PARAM_input.MOD_PARAM) '.m']), ...
    'PFT_OVERRIDE','', 'last_application_year',2019, 'biogeo_overrides',[], ...
    'saturation_exponent',1, 'apply_rock',apply_rock), topo.N_cells, 1); %#ok<NASGU>

old = cd(code); cleanup = onCleanup(@() cd(old));
MAIN_FRAME_BG_ERW_2D;
runoff_total_m3 = 0;
for cc_ = 1:topo.N_cells
    runoff_total_m3 = runoff_total_m3 + ...
        sum((Rh(:,1,cc_) + Rd(:,1,cc_))/1000) * topo.cell_area_m2(cc_);
end
end

function ok = raises(fcn,identifier)
ok = false;
try
    fcn();
catch ME
    ok = strcmp(ME.identifier,identifier);
    if ~ok
        fprintf('    raised %s: %s\n',ME.identifier,ME.message);
    end
end
end
