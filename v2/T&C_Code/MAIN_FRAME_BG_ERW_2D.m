%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%  MULTI-COLUMN MAIN FRAME  --  T&C-SMEW-N sub-plot domains %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Drives one soil-vegetation column per sub-plot cell of a
%% Build_Cell_Topology partition, coupling the columns laterally every
%% hydrologic timestep.
%%
%% RELATION TO THE SINGLE-COLUMN MAIN FRAME
%%   The science is unchanged: this script calls exactly the same unit
%%   functions as MAIN_FRAME_BG_ERW_10mins (HYDROLOGIC_UNIT, BIOGEO_UNIT,
%%   VEGGIE_UNIT, Plant_Exports, CropHeightType, ERW_function) with the
%%   same arguments.  What changes is that every state array carries a
%%   trailing column dimension and the columns exchange water and
%%   dissolved solute and heat after each hourly step.  With N_col = 1 and a
%%   single-cell topology it reproduces the single-column result, which
%%   is how it is verified (test_2d_driver_equivalence).
%%
%% WHAT IS RESOLVED IN 2D, AND WHAT IS NOT
%%   Water     lateral Buckingham-Darcy exchange between adjacent cells,
%%             per soil layer, every hour (Darcy_Lateral_Topology).
%%   Mass      advective transport of the aqueous phase with that water,
%%             upwind (Solute_Lateral_Transport).  Sorbed and solid
%%             phases do not move laterally. ERW-only Na, Al, residual
%%             anions and DIC use the same water flux in the active layer.
%%   Energy    lateral Fourier conduction and liquid-water sensible heat
%%             advection couple the soil thermal profiles.
%%   Gas       soil-air CO2 diffuses across cell faces with the existing
%%             Millington-Quirk effective diffusivity.
%%
%% METEOROLOGY
%%   All columns share one meteorological forcing series.  Without
%%   topography sub-domains diverge through vegetation, soil state and
%%   treatment, not through their driving weather; with OPT_DEM each
%%   cell's shortwave is projected onto its own slope, aspect and sky view
%%   (Topographic_Radiation_Projection).
%%
%% TOPOGRAPHY AND CELL GEOMETRY (from T&C-SMEW-N v1.0)
%%   OPT_DEM = 1 with a DEM attached to the topology
%%   (Attach_Topography_To_Topology) drives the lateral solvers on total
%%   head, routes runoff down the DEM to the outlet cell(s), and projects
%%   radiation per cell.  A polygon grid (Build_Polygon_Grid_Topology)
%%   weights every domain result by each cell's overlap with its polygon.
%%   Neither is used unless asked for, and without them the frame is the
%%   flat sub-plot frame bit for bit.
%%
%% EXPECTED IN THE CALLER WORKSPACE
%%   topo                Build_Cell_Topology or Build_Polygon_Grid_Topology
%%                       struct (N_col = topo.N_cells), with topography
%%                       attached when OPT_DEM = 1
%%   ColCfg              [N_col x 1] struct: .PARAM_IC .PFT_OVERRIDE
%%                       .apply_rock (logical)
%%   Parameter_Directory path holding the parameter scripts
%%   MOD_PARAM_input, ERW_PARAM_input, NN, NNd, ms, Zs, cc, dt, dth,
%%   Datam, Date, and the meteorological series, as for the 1-D frame.
%%
%% RECORDS AND PARALLEL SQUARES (package y)
%%   A one-column run keeps every hourly and ten-minute row, as before. A
%%   multi-square run keeps the last HOURLY_WINDOW hours of each hourly
%%   array, one day's block of ten-minute chemistry, and what the readers
%%   use per square and day (DAILY STORES); its result holds those stores,
%%   hourly field means of the tower fluxes and the final state, not the
%%   whole workspace. The daily, hourly and ERW work of the squares runs
%%   on the thread pool, square by square; nothing a square computes reads
%%   another square, so every number is the serial loop's, bit for bit.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N_col = topo.N_cells;
if ~exist('parallel_workers','var') || isempty(parallel_workers)
    parallel_workers = 0;
end
fprintf('\n=== 2D driver: %d columns over %.0f m2 ===\n', ...
    N_col, sum(topo.cell_area_m2));

if ~exist('OPT_SoilBiogeochemistry','var'); OPT_SoilBiogeochemistry = 1; end
if ~exist('OPT_ERW','var');                 OPT_ERW = 1;                 end
if ~exist('IrD','var') || isempty(IrD);     IrD = zeros(NN,1);           end

NCP = 8;    %% number of carbon pools

%% ========================================================================
%% PER-COLUMN PARAMETERS
%%   Each column's parameter script runs in its own private workspace, so
%%   columns may differ in PFT, management, rooting depth and soil.
%% ========================================================================
%% Build field by field: struct('f',M) with a non-scalar M would create a
%% struct ARRAY the size of M, not a scalar struct holding M.
ctx = struct();
ctx.MOD_PARAM_input = MOD_PARAM_input;
ctx.ERW_PARAM_input = ERW_PARAM_input;
ctx.NN = NN;   ctx.NNd = NNd;   ctx.ms = ms;
ctx.Datam = Datam;
ctx.Lat = Lat; ctx.Lon = Lon;   ctx.DeltaGMT = DeltaGMT;
ctx.t_bef = t_bef; ctx.t_aft = t_aft;
ctx.Ta = Ta;   ctx.Ca = Ca;
ctx.OPT_SoilBiogeochemistry = OPT_SoilBiogeochemistry;
ctx.NN_case_study_ERW_MOD = NN_case_study_ERW_MOD;
if exist('B_out','var');      ctx.B_out = B_out;           end
if exist('O_out','var');      ctx.O_out = O_out;           end
if exist('NavlI_out','var');  ctx.NavlI_out = NavlI_out;   end
if exist('RexmyI_out','var'); ctx.RexmyI_out = RexmyI_out; end
if exist('NupI_L_out','var'); ctx.NupI_L_out = NupI_L_out; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FIELD-SCALE HETEROGENEITY, two independent switches, both default OFF.
%%%   OPT_DEM         = 1  topography from a DEM: total-head Darcy,
%%%                        slope-normal faces, overland routing of runoff
%%%                        and its dissolved load, per-cell sky-view and
%%%                        beam-angle radiation.
%%%   OPT_CELL_INPUTS = 1  per-cell soil properties (pH, sand, clay, OM,
%%%                        CEC) from a per-cell agronomic table
%%%                        (Configure_Cell_Agronomic_Inputs) instead of one
%%%                        shared row.
%%% They are orthogonal at runtime -- the per-cell path reads only cell ids
%%% and plan coordinates, the DEM path adds only elevation-derived fields --
%%% so all four combinations are meaningful. Switched OFF, every branch
%%% below is skipped and the frame is bit-identical to the flat, uniform
%%% case. Ported from T&C-SMEW-N v1.0, which gated both implicitly by file
%%% presence; the explicit switches are this tree's addition.
%%%
%%% The polygon grid (Build_Polygon_Grid_Topology) is neither switch. It is
%%% the cell geometry, recognised by its overlap areas, and it changes only
%%% how the results are weighted (DOMAIN AGGREGATION) and what the overland
%%% router credits to the field.
if ~exist('OPT_DEM','var');         OPT_DEM         = 0; end
if ~exist('OPT_CELL_INPUTS','var'); OPT_CELL_INPUTS = 0; end
if ~(isequal(OPT_DEM,0) || isequal(OPT_DEM,1)) || ...
        ~(isequal(OPT_CELL_INPUTS,0) || isequal(OPT_CELL_INPUTS,1))
    error('MAIN_FRAME_BG_ERW_2D:BadSwitch', ...
        'OPT_DEM and OPT_CELL_INPUTS must each be 0 or 1.');
end

PAR = cell(N_col,1);
biogeo_exchange_scheme = repmat({'TC_GT'},N_col,1);
biogeo_KNmaxd = cell(N_col,1);
biogeo_overrides = cell(N_col,1);
for c = 1:N_col
    ovr = '';
    if isfield(ColCfg(c),'PFT_OVERRIDE') && ~isempty(ColCfg(c).PFT_OVERRIDE)
        ovr = fullfile(Parameter_Directory, [char(ColCfg(c).PFT_OVERRIDE) '.m']);
    end
    %%% Per-cell parameter context. Only under OPT_CELL_INPUTS does a column
    %%% get its own MOD/ERW parameter row; otherwise every column reads the
    %%% one shared context exactly as before. The ERW row in ColCfg already
    %%% exists on every run (control vs. treated), so that branch MUST be
    %%% gated on the switch -- taken unconditionally, as v1.0 has it, it
    %%% silently re-parameterises every margin column of every existing run.
    cell_ctx = ctx;
    if OPT_CELL_INPUTS == 1
        if isfield(ColCfg(c),'MOD_PARAM') && istable(ColCfg(c).MOD_PARAM) && ...
                height(ColCfg(c).MOD_PARAM) == 1
            cell_ctx.MOD_PARAM_input = ColCfg(c).MOD_PARAM;
        end
        if isfield(ColCfg(c),'ERW_PARAM') && istable(ColCfg(c).ERW_PARAM) && ...
                height(ColCfg(c).ERW_PARAM) == 1
            cell_ctx.ERW_PARAM_input = ColCfg(c).ERW_PARAM;
        end
    end
    PAR{c} = Load_Cell_Parameters(char(ColCfg(c).PARAM_IC), ovr, cell_ctx);
    if isfield(PAR{c},'bulk_pH_params') && ...
            isfield(PAR{c}.bulk_pH_params,'scheme')
        biogeo_exchange_scheme{c} = PAR{c}.bulk_pH_params.scheme;
    end
    if isfield(ColCfg(c),'KNmaxd')
        biogeo_KNmaxd{c} = ColCfg(c).KNmaxd;
    end
    if isfield(ColCfg(c),'biogeo_overrides')
        biogeo_overrides{c} = ColCfg(c).biogeo_overrides;
    end
    if isfield(ColCfg(c),'appelo_alpha') && ~isempty(ColCfg(c).appelo_alpha)
        if ~isfield(PAR{c},'bulk_pH_params')
            error('MAIN_FRAME_BG_ERW_2D:MissingBulkPHParameters', ...
                'AppeloAlpha requires bulk_pH_params in column %d.',c);
        end
        PAR{c}.bulk_pH_params.appelo.alpha = ColCfg(c).appelo_alpha;
    end
    if isfield(ColCfg(c),'cec_dpH') && ~isempty(ColCfg(c).cec_dpH)
        %%% Variable charge is a property of the soil, not of the exchange
        %%% scheme, but the sensitivity that reproduces a given titration
        %%% curve depends on the selectivity constants it is fitted with.
        %%% A scalar supplied here applies to whichever scheme is active;
        %%% per-scheme values are declared in the parameter file instead.
        if ~isfield(PAR{c},'bulk_pH_params')
            error('MAIN_FRAME_BG_ERW_2D:MissingBulkPHParameters', ...
                'CECpHSensitivity requires bulk_pH_params in column %d.',c);
        end
        PAR{c}.bulk_pH_params.dCEC_dpH = ColCfg(c).cec_dpH;
    end
    if isfield(ColCfg(c),'appelo_logKH') && ~isempty(ColCfg(c).appelo_logKH)
        if ~isfield(PAR{c},'bulk_pH_params')
            error('MAIN_FRAME_BG_ERW_2D:MissingBulkPHParameters', ...
                'AppeloLogKH requires bulk_pH_params in column %d.',c);
        end
        PAR{c}.bulk_pH_params.appelo.logKH = ColCfg(c).appelo_logKH;
    end
    fprintf('  column %2d  %-14s  %s\n', c, topo.cell_name{c}, ...
        ternary(isempty(ovr),'base parameters','PFT override'));
end
if exist('run_provenance','var') && isstruct(run_provenance)
    appelo_alpha_effective = nan(N_col,1);
    appelo_logKH_effective = nan(N_col,5);
    for c = 1:N_col
        if strcmpi(biogeo_exchange_scheme{c},'appelo')
            appelo_params = [];
            if isfield(PAR{c},'bulk_pH_params') && ...
                    isfield(PAR{c}.bulk_pH_params,'appelo')
                appelo_params = PAR{c}.bulk_pH_params.appelo;
            end
            [~,~,~,~,~,appelo_alpha_effective(c), ...
                appelo_logKH_effective(c,:)] = ...
                Resolve_Appelo_Constants(appelo_params);
        end
    end
    cec_dpH_effective = nan(N_col,1);
    cec_ref_pH_effective = nan(N_col,1);
    for c = 1:N_col
        if isfield(PAR{c},'bulk_pH_params')
            [cec_dpH_effective(c),cec_ref_pH_effective(c)] = ...
                ERW_Resolve_CEC_pH_Sensitivity(PAR{c}.bulk_pH_params, ...
                    char(biogeo_exchange_scheme{c}));
        end
    end
    run_provenance.exchange_scheme = biogeo_exchange_scheme;
    run_provenance.AppeloAlpha = appelo_alpha_effective;
    run_provenance.AppeloLogKH = appelo_logKH_effective;
    run_provenance.CECpHSensitivity = cec_dpH_effective;
    run_provenance.CECrefpH = cec_ref_pH_effective;
end

%% Soil layering must agree across columns: the lateral solver exchanges
%% water layer by layer, which is only defined on a shared discretisation.
for c = 2:N_col
    if ~isequal(PAR{c}.Zs, PAR{1}.Zs)
        error('MAIN_FRAME_BG_ERW_2D:LayerMismatch', ...
            ['Column %d uses a different soil layer mesh from column 1. ' ...
             'Lateral exchange is computed layer by layer and is undefined ' ...
             'across differing meshes.'], c);
    end
end
if any(abs(cellfun(@(p) p.Zbio,PAR)-PAR{1}.Zbio) > 1e-9)
    error('MAIN_FRAME_BG_ERW_2D:ActiveDepthMismatch', ...
        'All columns must use the same active biogeochemical depth.');
end
Zlayer_m = diff(PAR{1}.Zs)/1000;      %% [m] layer thickness, shared

%% TOPOGRAPHY LIVES IN THE TOPOLOGY, NEVER IN THE PARAMETER FILE. With
%% OPT_DEM the lateral solvers run on total head H = psi + z and the surface
%% router moves ponded runoff down the DEM; the vertical column stays 1-D
%% (Slo_top = 0 in the HYDROLOGIC_UNIT call), so its own kinematic Qi_out
%% cannot double-count the graph exchange. A non-zero parameter-file Slo_top
%% would re-introduce that double count and is rejected on every run.
%%
%% The flag is derived from BOTH the switch and the topology, and the two
%% must agree. A topology that carries elevation while the switch is off is
%% refused rather than silently engaging: the lateral solvers and the
%% router key on field presence, so a deserialised checkpoint or a reused
%% struct could otherwise turn topography on behind the switch.
topography_attached = OPT_DEM == 1 && ...
    isfield(topo,'topography_attached') && topo.topography_attached;
if OPT_DEM == 0 && isfield(topo,'cell_z_m') && ~isempty(topo.cell_z_m)
    error('MAIN_FRAME_BG_ERW_2D:TopographyBehindSwitch', ...
        ['The topology carries elevation (cell_z_m) but OPT_DEM is 0. ' ...
         'Either enable OPT_DEM or build the topology without topography.']);
end
if OPT_DEM == 1 && ~topography_attached
    error('MAIN_FRAME_BG_ERW_2D:NoTopographyAttached', ...
        ['OPT_DEM is 1 but no topography is attached to the topology. ' ...
         'Attach a DEM with Attach_Topography_To_Topology in the driver.']);
end
if N_col > 1
    for c = 1:N_col
        if abs(PAR{c}.Slo_top) > 1e-12 || any(abs(PAR{c}.Slo_pot) > 1e-12)
            error('MAIN_FRAME_BG_ERW_2D:SlopedDomainUnsupported', ...
                ['Column %d carries a non-zero parameter-file Slo_top/Slo_pot. ' ...
                 'Put field topography in the cell topology via ' ...
                 'Attach_Topography_To_Topology (topo.cell_z_m), not in the ' ...
                 'parameter file.'], c);
        end
    end
    if topography_attached
        fprintf(['  topography attached: %s; slope %.1f-%.1f%%, ' ...
            'relief %.2f m, %d outlet cell(s)\n'], ...
            char(string(topo.topography_source)), ...
            100*min(topo.cell_slope), 100*max(topo.cell_slope), ...
            max(topo.cell_z_m)-min(topo.cell_z_m), sum(topo.cell_is_outlet));
    end
end
%% Solar geometry window for the per-cell beam projection (v1.0).
%% SetSunVariables averages the sun over [HR-t_bef, HR+t_aft]; the run
%% script's window is [HR+0.67, HR+1.67], while an hour-ending forcing
%% (ERA5-Land) labels record HR with [HR-1, HR]. The two differ by ~1.7 h,
%% ~25 deg of azimuth at midday, enough to reverse the east/west aspect
%% contrast the projection exists to resolve. The global t_bef/t_aft is
%% left alone; only the per-cell BEAM PROJECTION uses this window, which
%% defaults to hour-ending [1 0]. Override with TOPO_SOLAR_WINDOW =
%% [t_bef t_aft]. Inert without topography.
if ~exist('TOPO_SOLAR_WINDOW','var') || isempty(TOPO_SOLAR_WINDOW)
    TOPO_SOLAR_WINDOW = [1 0];
end
topo_t_bef = TOPO_SOLAR_WINDOW(1); topo_t_aft = TOPO_SOLAR_WINDOW(2);
if topography_attached
    %% Report how well each window tracks the observed shortwave, on a
    %% subsample, so a forcing with another time convention shows in the
    %% run log instead of silently rotating every aspect.
    probe = 1:max(1,floor(NN/500)):NN;
    Rsw_probe = SAB1(probe)+SAB2(probe)+SAD1(probe)+SAD2(probe);
    sh_topo = zeros(numel(probe),1); sh_run = zeros(numel(probe),1);
    for pk = 1:numel(probe)
        sh_topo(pk) = max(0,sin(SetSunVariables(Datam(probe(pk),:), ...
            DeltaGMT,Lon,Lat,topo_t_bef,topo_t_aft)));
        sh_run(pk) = max(0,sin(SetSunVariables(Datam(probe(pk),:), ...
            DeltaGMT,Lon,Lat,t_bef,t_aft)));
    end
    if std(Rsw_probe(:)) > 0 && std(sh_topo) > 0 && std(sh_run) > 0
        r_topo = corrcoef(sh_topo,Rsw_probe(:)); c_topo = r_topo(1,2);
        r_run  = corrcoef(sh_run,Rsw_probe(:));  c_run  = r_run(1,2);
        fprintf(['  topographic solar window [HR%+.2f, HR%+.2f]: sun-vs-' ...
            'shortwave correlation %.3f (run-script window %.3f)\n'], ...
            -topo_t_bef,topo_t_aft,c_topo,c_run);
        if c_topo < c_run
            warning('MAIN_FRAME_BG_ERW_2D:SolarWindowMismatch', ...
                ['The forcing tracks the run-script solar window (%.3f) ' ...
                 'better than the topographic one (%.3f); set ' ...
                 'TOPO_SOLAR_WINDOW to match this forcing time convention ' ...
                 'or the per-cell beam projection will be rotated.'], ...
                c_run,c_topo);
        end
    end
end

bio_active_fraction = zeros(N_col,ms);
for c = 1:N_col
    active_thickness_m = PAR{c}.Bio_Zs*(PAR{c}.Zbio/1000);
    bio_active_fraction(c,:) = active_thickness_m./Zlayer_m;
end
if any(~isfinite(bio_active_fraction(:))) || ...
        any(bio_active_fraction(:) < -1e-12) || ...
        any(bio_active_fraction(:) > 1+1e-12)
    error('MAIN_FRAME_BG_ERW_2D:BadActiveLayer', ...
        'Biogeochemical active-layer fractions must lie within each soil layer.');
end
bio_active_fraction = min(1,max(0,bio_active_fraction));

%% Per-column soil parameters for the lateral solver
soilpar_cell = Soil_Param_To_Routing(PAR{1}.Soil_Param, PAR{1}.SPAR, PAR{1}.aR);
soilpar_cell = repmat(soilpar_cell, N_col, 1);
for c = 1:N_col
    soilpar_cell(c) = Soil_Param_To_Routing(PAR{c}.Soil_Param, PAR{c}.SPAR, PAR{c}.aR);
end
A_cell = topo.cell_area_m2(:);             %% [m2]
A_tot  = sum(A_cell);

%% ========================================================================
%% RECORDS: THE WHOLE RUN OR A ROLLING WINDOW (package y)
%% ========================================================================
%% A one-column run keeps every hourly and ten-minute row, as before, so the
%% 1-D read-outs, validations and the equivalence test see the same arrays.
%% A multi-square run keeps only what the model reads again: the last
%% HOURLY_WINDOW hours of each hourly array (hour i-1, the 24 hours the
%% daily step averages, and the ERW water bins, which reach back at most to
%% hour i-49) and one day's block of ten-minute chemistry; what the readers
%% use is stored per square and day (DAILY STORES below). TC_FULL_RECORDS
%% = 1 or 0, or FULL_RECORDS in the caller's workspace, forces either.
HOURLY_WINDOW = 50;
if ~exist('FULL_RECORDS','var') || isempty(FULL_RECORDS)
    FULL_RECORDS = N_col == 1;
    records_choice = strtrim(getenv('TC_FULL_RECORDS'));
    if strcmp(records_choice,'1')
        FULL_RECORDS = true;
    elseif strcmp(records_choice,'0')
        FULL_RECORDS = false;
    elseif ~isempty(records_choice)
        error('MAIN_FRAME_BG_ERW_2D:BadRecords','TC_FULL_RECORDS must be 0 or 1.');
    end
end
FULL_RECORDS = logical(FULL_RECORDS);
%% Hour h of an hourly array sits in row mod(h-1,NH)+1.
if FULL_RECORDS
    NH = NN;
else
    NH = min(NN,HOURLY_WINDOW);
end

%% ========================================================================
%% STATE ALLOCATION   (trailing dimension = column; hourly arrays hold
%% the last NH hours, daily arrays every day)
%% ========================================================================
z2 = @(a,b) zeros(a,b,N_col);
z3 = @(a,b,d) zeros(a,b,d,N_col);

V=z2(NH,ms); O=z2(NH,ms); Qi_out=z2(NH,ms); WTR=z2(NH,ms); POT=z2(NH,ms);
Tdp=z2(NH,ms); Vice=z2(NH,ms); Oice=z2(NH,ms);
OF=z2(NH,1); OS=z2(NH,1); ZWT=z2(NH,1);
Pr_sno=z2(NH,1); Pr_liq=z2(NH,1); er=z2(NH,1);
Ts=z2(NH,1); ra=z2(NH,1); r_soil=z2(NH,1); b_soil=z2(NH,1); Tice=z2(NH,1);
Rn=z2(NH,1); H=z2(NH,1); QE=z2(NH,1); Qv=z2(NH,1); Lpho=z2(NH,1);
EG=z2(NH,1); ESN=z2(NH,1); ESN_In=z2(NH,1); EWAT=z2(NH,1);
EIn_urb=z2(NH,1); EIn_rock=z2(NH,1); dw_SNO=z2(NH,1);
G=z2(NH,1); SWE=z2(NH,1); SND=z2(NH,1); ros=z2(NH,1);
In_SWE=z2(NH,1); SP_wc=z2(NH,1); WR_SP=z2(NH,1); U_SWE=z2(NH,1); NIn_SWE=z2(NH,1);
dQ=z2(NH,1); Qfm=z2(NH,1); t_sls=z2(NH,1); DQ=z2(NH,1); DT=z2(NH,1);
In_urb=z2(NH,1); In_rock=z2(NH,1); SE_rock=z2(NH,1); SE_urb=z2(NH,1);
Csno=z2(NH,1); WIS=z2(NH,1); Lk=z2(NH,1); f=z2(NH,1); Rh=z2(NH,1); Rd=z2(NH,1);
NDVI=z2(NH,1); alp_soil=z2(NH,1); tau_sno=z2(NH,1); e_sno=z2(NH,1); ALB=z2(NH,1);
EK=z2(NH,1); dQVEG=z2(NH,1); TsV=z2(NH,1); HV=z2(NH,1); QEV=z2(NH,1);
Cice=z2(NH,1); Lk_wat=z2(NH,1); Lk_rock=z2(NH,1); EICE=z2(NH,1);
WAT=z2(NH,1); ICE=z2(NH,1); ICE_D=z2(NH,1); WR_IP=z2(NH,1); NIce=z2(NH,1);
Cicew=z2(NH,1); IP_wc=z2(NH,1); Csnow=z2(NH,1); FROCK=z2(NH,1);
Imelt=z2(NH,1); Smelt=z2(NH,1); Tdamp=z2(NH,1); Gfin=z2(NH,1);
In_Litter=z2(NH,1); ELitter=z2(NH,1); Ws_under=z2(NH,1);
r_litter=z2(NH,cc);
In_H=z2(NH,cc); In_L=z2(NH,cc); OH=z2(NH,cc); OL=z2(NH,cc);
Psi_s_H=z2(NH,cc); Psi_s_L=z2(NH,cc); Tdp_H=z2(NH,cc); Tdp_L=z2(NH,cc);
rap_H=z2(NH,cc); rap_L=z2(NH,cc); Dr_H=z2(NH,cc); Dr_L=z2(NH,cc);
T_H=z2(NH,cc); T_L=z2(NH,cc); EIn_H=z2(NH,cc); EIn_L=z2(NH,cc);
rb_H=z2(NH,cc); rb_L=z2(NH,cc);
An_H=z2(NH,cc); Rdark_H=z2(NH,cc); An_L=z2(NH,cc); Rdark_L=z2(NH,cc);
Ci_sunH=z2(NH,cc); Ci_sunL=z2(NH,cc); Ci_shdH=z2(NH,cc); Ci_shdL=z2(NH,cc);
rs_sunH=z2(NH,cc); rs_sunL=z2(NH,cc); rs_shdH=z2(NH,cc); rs_shdL=z2(NH,cc);
Vx_H=z2(NH,cc); Vl_H=z2(NH,cc); Vx_L=z2(NH,cc); Vl_L=z2(NH,cc);
Psi_x_H=z2(NH,cc); Psi_l_H=z2(NH,cc); Psi_x_L=z2(NH,cc); Psi_l_L=z2(NH,cc);
gsr_H=z2(NH,cc); Jsx_H=z2(NH,cc); Jxl_H=z2(NH,cc); Kleaf_H=z2(NH,cc); Kx_H=z2(NH,cc);
gsr_L=z2(NH,cc); Jsx_L=z2(NH,cc); Jxl_L=z2(NH,cc); Kleaf_L=z2(NH,cc); Kx_L=z2(NH,cc);
fapar_H=z2(NH,cc); fapar_L=z2(NH,cc); SIF_H=z2(NH,cc); SIF_L=z2(NH,cc);
Tdeb = zeros(NH, max(1,length(PAR{1}.Zs_deb)-1), N_col);

%% Daily vegetation and biogeochemistry state
LAI_L=z2(NNd,cc); B_L=z3(NNd,cc,NCP); NPP_L=z2(NNd,cc); Rg_L=z2(NNd,cc);
RA_L=z2(NNd,cc); Rmc_L=z2(NNd,cc); ANPP_L=z2(NNd,cc);
Rms_L=z2(NNd,cc); Rmr_L=z2(NNd,cc); PHE_S_L=z2(NNd,cc); dflo_L=z2(NNd,cc);
AgeL_L=z2(NNd,cc); e_rel_L=ones(NNd,cc,N_col); SAI_L=z2(NNd,cc); hc_L=z2(NNd,cc);
e_relN_L=ones(NNd,cc,N_col); BA_L=z2(NNd,cc);
LAI_H=z2(NNd,cc); B_H=z3(NNd,cc,NCP); NPP_H=z2(NNd,cc); Rg_H=z2(NNd,cc);
RA_H=z2(NNd,cc); Rms_H=z2(NNd,cc); Rmr_H=z2(NNd,cc); ANPP_H=z2(NNd,cc);
PHE_S_H=z2(NNd,cc); dflo_H=z2(NNd,cc); Rmc_H=z2(NNd,cc);
AgeL_H=z2(NNd,cc); e_rel_H=ones(NNd,cc,N_col); SAI_H=z2(NNd,cc); hc_H=z2(NNd,cc);
e_relN_H=ones(NNd,cc,N_col); BA_H=z2(NNd,cc);
Rrootl_H=z2(NNd,cc); Rrootl_L=z2(NNd,cc);
Bfac_dayH=ones(NNd,cc,N_col); Bfac_weekH=ones(NNd,cc,N_col);
NPPI_H=z2(NNd,cc); TdpI_H=z2(NNd,cc);
Bfac_dayL=ones(NNd,cc,N_col); Bfac_weekL=ones(NNd,cc,N_col);
NPPI_L=z2(NNd,cc); TdpI_L=z2(NNd,cc);
Sr_H=z2(NNd,cc); Slf_H=z2(NNd,cc); Sfr_H=z2(NNd,cc); Swm_H=z2(NNd,cc); Sll_H=z2(NNd,cc);
Sr_L=z2(NNd,cc); Slf_L=z2(NNd,cc); Sfr_L=z2(NNd,cc); Swm_L=z2(NNd,cc); Sll_L=z2(NNd,cc);
LAIdead_H=z2(NNd,cc); LAIdead_L=z2(NNd,cc);
Rexmy_H=z3(NNd,cc,3); AgeDL_H=z2(NNd,cc);
Nreserve_H=z2(NNd,cc); Preserve_H=z2(NNd,cc); Kreserve_H=z2(NNd,cc);
Careserve_H=z2(NNd,cc); Mgreserve_H=z2(NNd,cc); Sireserve_H=z2(NNd,cc);
rNc_H=z2(NNd,cc); rPc_H=z2(NNd,cc); rKc_H=z2(NNd,cc);
rCac_H=z2(NNd,cc); rMgc_H=z2(NNd,cc); rSic_H=z2(NNd,cc);
Rexmy_L=z3(NNd,cc,3); AgeDL_L=z2(NNd,cc);
Nreserve_L=z2(NNd,cc); Preserve_L=z2(NNd,cc); Kreserve_L=z2(NNd,cc);
Careserve_L=z2(NNd,cc); Mgreserve_L=z2(NNd,cc); Sireserve_L=z2(NNd,cc);
rNc_L=z2(NNd,cc); rPc_L=z2(NNd,cc); rKc_L=z2(NNd,cc);
rCac_L=z2(NNd,cc); rMgc_L=z2(NNd,cc); rSic_L=z2(NNd,cc);
NBLeaf_H=z2(NNd,cc); PARI_H=z3(NNd,cc,3); NBLI_H=z2(NNd,cc);
NBLeaf_L=z2(NNd,cc); PARI_L=z3(NNd,cc,3); NBLI_L=z2(NNd,cc);
NupI_H=z3(NNd,cc,6); NupI_L=z3(NNd,cc,6); NavlI=z2(NNd,6);
RB_L=z3(NNd,cc,7); RB_H=z3(NNd,cc,7);
NuLit_H=z3(NNd,cc,6); NuLit_L=z3(NNd,cc,6);
BLit=z2(NNd,cc);
Nuptake_H=z2(NNd,cc); NH4uptake_H=z2(NNd,cc); NO3uptake_H=z2(NNd,cc);
Puptake_H=z2(NNd,cc); Kuptake_H=z2(NNd,cc);
Cauptake_H=z2(NNd,cc); Mguptake_H=z2(NNd,cc); Siuptake_H=z2(NNd,cc);
FNC_H=z2(NNd,cc);
TexC_H=z2(NNd,cc); TexN_H=z2(NNd,cc); TexP_H=z2(NNd,cc); TexK_H=z2(NNd,cc);
TNIT_H=z2(NNd,cc); TPHO_H=z2(NNd,cc); TPOT_H=z2(NNd,cc);
SupN_H=z2(NNd,cc); SupP_H=z2(NNd,cc); SupK_H=z2(NNd,cc);
TexCa_H=z2(NNd,cc); TexMg_H=z2(NNd,cc); TexSi_H=z2(NNd,cc);
TCAL_H=z2(NNd,cc); TMAG_H=z2(NNd,cc); TSIL_H=z2(NNd,cc);
SupCa_H=z2(NNd,cc); SupMg_H=z2(NNd,cc); SupSi_H=z2(NNd,cc);
Nuptake_L=z2(NNd,cc); NH4uptake_L=z2(NNd,cc); NO3uptake_L=z2(NNd,cc);
Puptake_L=z2(NNd,cc); Kuptake_L=z2(NNd,cc);
Cauptake_L=z2(NNd,cc); Mguptake_L=z2(NNd,cc); Siuptake_L=z2(NNd,cc);
FNC_L=z2(NNd,cc);
TexC_L=z2(NNd,cc); TexN_L=z2(NNd,cc); TexP_L=z2(NNd,cc); TexK_L=z2(NNd,cc);
TNIT_L=z2(NNd,cc); TPHO_L=z2(NNd,cc); TPOT_L=z2(NNd,cc);
SupN_L=z2(NNd,cc); SupP_L=z2(NNd,cc); SupK_L=z2(NNd,cc);
TexCa_L=z2(NNd,cc); TexMg_L=z2(NNd,cc); TexSi_L=z2(NNd,cc);
TCAL_L=z2(NNd,cc); TMAG_L=z2(NNd,cc); TSIL_L=z2(NNd,cc);
SupCa_L=z2(NNd,cc); SupMg_L=z2(NNd,cc); SupSi_L=z2(NNd,cc);
ISOIL_H=z3(NNd,cc,27); ISOIL_L=z3(NNd,cc,27);
P=z2(NNd,79);
R_litter=z2(NNd,1); R_microbe=z2(NNd,1); R_litter_sur=z2(NNd,1); R_ew=z2(NNd,1);
VOL=z2(NNd,1); N2flx=z2(NNd,1); N2Oflx=z2(NNd,1); ...
    NITflx=z2(NNd,1); DENITflx=z2(NNd,1); N2Oflx_nit=z2(NNd,1); ...
    FertN_d=z2(NNd,1); DepN_d=z2(NNd,1); PH_d=z2(NNd,1); ...
    PH_pore_d=z2(NNd,1); PH_bulk_d=z2(NNd,1);
Min_N=z2(NNd,1); Min_P=z2(NNd,1); R_bacteria=z2(NNd,1);
RmycAM=z2(NNd,1); RmycEM=z2(NNd,1); Prod_B=z2(NNd,1); Prod_F=z2(NNd,1);
BfixN=z2(NNd,1); LitFirEmi=z2(NNd,2);
LEAK_NH4=z2(NNd,1); LEAK_NO3=z2(NNd,1); LEAK_P=z2(NNd,1); LEAK_K=z2(NNd,1);
LEAK_DOC=z2(NNd,1); LEAK_DON=z2(NNd,1); LEAK_DOP=z2(NNd,1);
LEAK_Ca=z2(NNd,1); LEAK_Mg=z2(NNd,1); LEAK_Si=z2(NNd,1);
RexmyI=z2(NNd,3);
AgrHarNut=z2(NNd,6);
SeedInput=z2(NNd,7); %% [C N P K Ca Mg Si], [g element m-2 ground]
Ccrown_t=ones(NNd,cc,N_col);
AgePl_H=z2(NNd,cc); AgePl_L=z2(NNd,cc); Tden_H=z2(NNd,cc); Tden_L=z2(NNd,cc);
TBio_Ht=z2(NNd,cc); TBio_Lt=z2(NNd,cc);
q_runon=z2(NH,1); Qi_in=z2(NH,ms);
ManIH=zeros(cc,N_col); ManIL=zeros(cc,N_col);

%% Per-column scalars carried between steps
Ccrown_c = zeros(N_col,cc);  Cbare_c = zeros(N_col,1);
ZR95_H_c = zeros(N_col,cc);  ZR95_L_c = zeros(N_col,cc);
RfH_Zs_c = cell(N_col,1);    RfL_Zs_c = cell(N_col,1);
PH_c     = zeros(N_col,1);   PH_bulk_c = zeros(N_col,1); Tstm0_c = zeros(N_col,1);
Bam_c    = zeros(N_col,1);   Bem_c    = zeros(N_col,1);
snow_alb_c = cell(N_col,1);
B_IO_c   = cell(N_col,1);
TBio_H_c = zeros(N_col,cc);  TBio_L_c = zeros(N_col,cc);

for c = 1:N_col
    p = PAR{c};
    Ccrown_c(c,:) = p.Ccrown;   Cbare_c(c) = p.Cbare;
    ZR95_H_c(c,:) = p.ZR95_H;   ZR95_L_c(c,:) = p.ZR95_L;
    RfH_Zs_c{c} = p.RfH_Zs;     RfL_Zs_c{c} = p.RfL_Zs;
    %% Cold-start pH seed. p.PH is the OBSERVED BULK-SOIL pH from the
    %% parameter file, the quantity BIOGEO_UNIT consumes, so seeding
    %% PH_bulk_c from it is correct. The pore-water pH the ERW chemistry
    %% starts from is only known once ERW_MOD_PARAM_excel has inverted that
    %% observation (bulk_pH_params.initial_pH_mode = 'bulk_observed',
    %% Invert_Bulk_pH), which happens inside the first ERW_Daily_Step_2D
    %% call; the daily records written before that call are backfilled
    %% from the ERW initial condition there ("cold-start backfill" below).
    %% Before the first ERW call PH_c is only recorded (PH_pore_d), never
    %% consumed, so seeding it from the bulk value is harmless.
    PH_c(c) = p.PH(1);          PH_bulk_c(c) = p.PH(1); snow_alb_c{c} = p.snow_alb;
    if isfield(p,'ERW_biogeochem') && isfield(p.ERW_biogeochem,'pH_bulk') && ...
            isfinite(p.ERW_biogeochem.pH_bulk(1))
        PH_bulk_c(c) = p.ERW_biogeochem.pH_bulk(1);
    end
    B_IO_c{c} = p.B_IO;
    TBio_H_c(c,:) = p.TBio_H;   TBio_L_c(c,:) = p.TBio_L;

    %% Seed row 1 of every state array from this column's parameter file
    O(1,:,c) = p.O(1,:);  V(1,:,c) = p.V(1,:);
    Ts(1,c) = p.Ts(1);  Tstm0_c(c) = p.Ts(1);
    Tdp(1,:,c) = p.Tdp(1,:);  Tdamp(1,c) = p.Tdamp(1);
    Tdp_H(1,:,c) = p.Tdp_H(1,:);  Tdp_L(1,:,c) = p.Tdp_L(1,:);
    TdpI_H(1,:,c) = p.TdpI_H(1,:); TdpI_L(1,:,c) = p.TdpI_L(1,:);
    SWE(1,c)=p.SWE(1); SND(1,c)=p.SND(1); ros(1,c)=p.ros(1);
    In_SWE(1,c)=p.In_SWE(1); SP_wc(1,c)=p.SP_wc(1); t_sls(1,c)=p.t_sls(1);
    e_sno(1,c)=p.e_sno(1); tau_sno(1,c)=p.tau_sno(1); EK(1,c)=p.EK(1);
    WAT(1,c)=p.WAT(1); ICE(1,c)=p.ICE(1); IP_wc(1,c)=p.IP_wc(1);
    ICE_D(1,c)=p.ICE_D(1); FROCK(1,c)=p.FROCK(1); Ws_under(1,c)=p.Ws_under(1);
    In_L(1,:,c)=p.In_L(1,:); In_H(1,:,c)=p.In_H(1,:);
    In_urb(1,c)=p.In_urb(1); In_rock(1,c)=p.In_rock(1); In_Litter(1,c)=p.In_Litter(1);
    Ci_sunL(1,:,c)=p.Ci_sunL(1,:); Ci_sunH(1,:,c)=p.Ci_sunH(1,:);
    Ci_shdL(1,:,c)=p.Ci_shdL(1,:); Ci_shdH(1,:,c)=p.Ci_shdH(1,:);
    LAI_H(1,:,c)=p.LAI_H(1,:); LAI_L(1,:,c)=p.LAI_L(1,:);
    B_H(1,:,:,c)=p.B_H(1,:,:);  B_L(1,:,:,c)=p.B_L(1,:,:);
    Rrootl_H(1,:,c)=p.Rrootl_H(1,:); Rrootl_L(1,:,c)=p.Rrootl_L(1,:);
    PHE_S_H(1,:,c)=p.PHE_S_H(1,:); PHE_S_L(1,:,c)=p.PHE_S_L(1,:);
    dflo_H(1,:,c)=p.dflo_H(1,:); dflo_L(1,:,c)=p.dflo_L(1,:);
    AgeL_H(1,:,c)=p.AgeL_H(1,:); AgeL_L(1,:,c)=p.AgeL_L(1,:);
    e_rel_H(1,:,c)=p.e_rel_H(1,:); e_rel_L(1,:,c)=p.e_rel_L(1,:);
    hc_H(1,:,c)=p.hc_H(1,:); hc_L(1,:,c)=p.hc_L(1,:);
    SAI_H(1,:,c)=p.SAI_H(1,:); SAI_L(1,:,c)=p.SAI_L(1,:);
    Nreserve_H(1,:,c)=p.Nreserve_H(1,:); Preserve_H(1,:,c)=p.Preserve_H(1,:);
    Kreserve_H(1,:,c)=p.Kreserve_H(1,:);
    Careserve_H(1,:,c)=p.Careserve_H(1,:); Mgreserve_H(1,:,c)=p.Mgreserve_H(1,:);
    Sireserve_H(1,:,c)=p.Sireserve_H(1,:);
    Nreserve_L(1,:,c)=p.Nreserve_L(1,:); Preserve_L(1,:,c)=p.Preserve_L(1,:);
    Kreserve_L(1,:,c)=p.Kreserve_L(1,:);
    Careserve_L(1,:,c)=p.Careserve_L(1,:); Mgreserve_L(1,:,c)=p.Mgreserve_L(1,:);
    Sireserve_L(1,:,c)=p.Sireserve_L(1,:);
    FNC_H(1,:,c)=p.FNC_H(1,:); FNC_L(1,:,c)=p.FNC_L(1,:);
    NupI_H(1,:,:,c)=p.NupI_H(1,:,:); NupI_L(1,:,:,c)=p.NupI_L(1,:,:);
    RexmyI(1,:,c)=p.RexmyI(1,:); NavlI(1,:,c)=p.NavlI(1,:);
    Vx_H(1,:,c)=p.Vx_H; Vl_H(1,:,c)=p.Vl_H;
    Vx_L(1,:,c)=p.Vx_L; Vl_L(1,:,c)=p.Vl_L;
    BLit(1,:,c)=p.BLit(1);
    if OPT_SoilBiogeochemistry == 1
        P(1,:,c) = p.P(1,:);
        PH_d(1,c) = PH_bulk_c(c);
        PH_pore_d(1,c) = p.PH(1);
        PH_bulk_d(1,c) = PH_bulk_c(c);
    end
    Check_Land_Cover_Fractions(p.Crock,p.Curb,p.Cwat,p.Cbare,p.Ccrown);
end
CcrownFIX = Ccrown_c;
[~,~,lateral_species,lateral_pool_idx,lateral_fraction_fields] = Biogeo_Lateral_Pools( ...
    reshape(P(1,:,:),79,N_col).');
N_lateral_species = numel(lateral_species);
i_lateral_NO3 = find(strcmp(lateral_species,'NO3'),1);
lateral_fraction = zeros(N_col,N_lateral_species);
BiogeoPar_lateral = Biogeochemistry_Parameter(1);
for s = 1:N_lateral_species
    lateral_fraction(:,s) = BiogeoPar_lateral.(lateral_fraction_fields{s});
end
if any(~isfinite(lateral_fraction(:))) || any(lateral_fraction(:) < 0) || ...
        any(lateral_fraction(:) > 1)
    error('MAIN_FRAME_BG_ERW_2D:BadMobileFraction', ...
        'Biogeo_Leakage mobility fractions must be finite and within [0,1].');
end

%% ========================================================================
%% NUMERICAL OPTIONS  (identical to the single-column frame)
%% ========================================================================
Opt_CR = optimset('TolFun',1);
Opt_ST = optimset('TolFun',0.1);
OPT_SM  = odeset('AbsTol',0.05,'MaxStep',dth);
OPT_VD  = odeset('AbsTol',0.05);
OPT_STh = odeset('AbsTol',0.02);
OPT_PH  = odeset('AbsTol',0.01);
OPT_FR_SOIL = 1; OPT_VegSnow = 1; OPT_SoilTemp = 1;
if ~exist('OPT_VCA','var');  OPT_VCA = 0; OPT_ALLOME = 0; end
if ~exist('OPT_PlantHydr','var');      OPT_PlantHydr = 0;      end
if ~exist('OPT_EnvLimitGrowth','var'); OPT_EnvLimitGrowth = 1; end
OPT_WET = 0;

%% ERW per-column state
%% Cell array, not a struct array: each column's ERW state gains fields
%% after its first call, which a struct array cannot absorb elementwise.
ERW = cell(N_col,1);
for c = 1:N_col; ERW{c} = struct('init',false,'has_update',false); end
conv_mol = 1e6; conv_Al = 1e3;
Tdp_bio_c = zeros(N_col,1); VSUM_c = zeros(N_col,1); VTSUM_c = zeros(N_col,1);
lateral_diag = struct('water_buf_to_main', zeros(NNd,1), ...
                      'mass_buf_to_main',  zeros(NNd,N_lateral_species), ...
                      'NO3_buf_to_main',   zeros(NNd,1), ...
                      'heat_buf_to_main_J',zeros(NNd,1), ...
                      'heat_residual_J',   zeros(NNd,1), ...
                      'heat_residual_relative',zeros(NNd,1), ...
                      'heat_sub_steps_max',zeros(NNd,1), ...
                      'erw_mass_buf_to_main_umol',zeros(NNd,4), ...
                      'erw_CO2_buf_to_main_umol',zeros(NNd,1), ...
                      'erw_mass_residual_umol',zeros(NNd,4), ...
                      'erw_mass_abs_residual_umol',zeros(NNd,4), ...
                      'erw_mass_throughput_umol',zeros(NNd,4), ...
                      'sub_steps_max',     zeros(NNd,1), ...
                      'species_names',     {lateral_species}, ...
                      'erw_species_names', {{'Na','Al','An','DIC'}});

%% Domain-boundary accounting. The subsurface lateral perimeter is closed;
%% bottom leakage and unrouted surface runoff leave the modeled plot.
exit_water_vert_m3 = zeros(N_col,1);
exit_water_surface_m3 = zeros(N_col,1);
%% Overland routing (OPT_DEM). Runoff produced this step is routed down the
%% DEM and offered to the NEXT step as run-on; T&C's operator split is one
%% step lagged by design, and overland travel over a farm field is far
%% below one hydrological step. Solute that arrived with run-on and has not
%% yet infiltrated sits in a POND store: crediting it straight to the soil
%% would strand it in the first absorbing cell it reached. All zero and
%% never touched when the switch is off.
q_runon_next        = zeros(N_col,1);                  %% [mm/h]
pond_mass_area      = zeros(N_col,N_lateral_species);  %% [gX/m^2]
exit_mass_surface_g = zeros(N_lateral_species,1);      %% [g] at outlets
erw_pond_umol       = zeros(N_col,4);                  %% Na, Al, An, DIC
exit_mass_surface_erw_umol = zeros(4,1);
%% Field-origin provenance (v1.0). On a polygon grid a perimeter cell is a
%% whole square only partly inside the polygon, so the routed water and
%% solute that came from inside it are carried through the router and the
%% pond store as passive tracers; a stop-over on a partial cell then cannot
%% re-label the load with that cell's own overlap fraction. On a
%% rectangular grid every overlap fraction is 1 and the tracers equal the
%% totals.
field_grid = isfield(topo,'cell_overlap_area_m2');
[~,~,field_fraction] = Field_Aggregation_Areas(topo);
pond_field_water_mm  = zeros(N_col,1);                 %% [mm]
pond_field_mass_area = zeros(N_col,N_lateral_species); %% [gX/m^2]
erw_pond_field_umol  = zeros(N_col,4);
exit_water_field_surface_m3      = 0;
exit_mass_surface_field_g        = zeros(N_lateral_species,1);
exit_mass_surface_erw_field_umol = zeros(4,1);
if OPT_DEM == 1
    %% Added only when on, so Domain.lateral keeps its field set otherwise.
    lateral_diag.surface_water_out_m3   = zeros(NNd,1);
    lateral_diag.surface_mass_out_g     = zeros(NNd,N_lateral_species);
    lateral_diag.surface_reinfilt_m3    = zeros(NNd,1);
    %% Subsurface lateral throughput per day, the sum over faces, layers
    %% and hours of |Q_edge| [m3], and the net redistribution, half the
    %% sum over cells of each cell's net gain or loss [m3]: how much water
    %% the topography actually moved between cells (v1.0).
    lateral_diag.water_edge_abs_m3      = zeros(NNd,1);
    lateral_diag.water_net_m3           = zeros(NNd,1);
end
%% Each cell's net subsurface lateral gain per day (positive = received
%% from its neighbours), water [mm] and NO3 [g N m-2], for the field maps
%% (Plot_Field_Maps 'lateral_water', 'lateral_NO3'). Kept on a field grid
%% or with topography only, so the plot runs keep Domain.lateral's fields.
if field_grid || OPT_DEM == 1
    lateral_diag.water_cell_net_mm      = zeros(NNd,N_col);
    lateral_diag.NO3_cell_net_g_m2      = zeros(NNd,N_col);
end

%% ========================================================================
%% DAILY STORES AND HOURLY FIELD MEANS (package y)
%% ========================================================================
%% What the readers use, kept for every square and day (the default list
%% Greg approved on 8 Oct 2026, PACKAGE_Y_OUTPUT_LISTS.md):
%%   - the daily arrays in DAILY_STORE_NAMES, every row, as they are;
%%   - from the hourly records, per daily row (the hours its daily step
%%     averaged: row 1 none, row 2 the first hour alone): the mean water
%%     content of the active layers O_d and water-table depth ZWT_d, and
%%     the totals of drainage Lk_d, infiltration f_d, run-on q_runon_d and
%%     runoff Rh_d, Rd_d [mm]; the hours after the last daily step go to
%%     Hourly_tail at the end;
%%   - from the ten-minute chemistry, one entry per ERW day and square
%%     (ERW_daily{c}, ERW_Daily_Record);
%%   - every hour, the field means of the tower-comparison fluxes
%%     (Field_hourly, Field_Hourly_Means), weighted as Domain weighs.
%% Everything else a daily-store result keeps is the final state.
DAILY_NAMES = { ...
    'LAI_L','B_L','NPP_L','Rg_L','RA_L','Rmc_L', ...
    'ANPP_L','Rms_L','Rmr_L','PHE_S_L','dflo_L','AgeL_L', ...
    'e_rel_L','SAI_L','hc_L','e_relN_L','BA_L','LAI_H', ...
    'B_H','NPP_H','Rg_H','RA_H','Rms_H','Rmr_H', ...
    'ANPP_H','PHE_S_H','dflo_H','Rmc_H','AgeL_H','e_rel_H', ...
    'SAI_H','hc_H','e_relN_H','BA_H','Rrootl_H','Rrootl_L', ...
    'Bfac_dayH','Bfac_weekH','NPPI_H','TdpI_H','Bfac_dayL','Bfac_weekL', ...
    'NPPI_L','TdpI_L','Sr_H','Slf_H','Sfr_H','Swm_H', ...
    'Sll_H','Sr_L','Slf_L','Sfr_L','Swm_L','Sll_L', ...
    'LAIdead_H','LAIdead_L','Rexmy_H','AgeDL_H','Nreserve_H','Preserve_H', ...
    'Kreserve_H','Careserve_H','Mgreserve_H','Sireserve_H','rNc_H','rPc_H', ...
    'rKc_H','rCac_H','rMgc_H','rSic_H','Rexmy_L','AgeDL_L', ...
    'Nreserve_L','Preserve_L','Kreserve_L','Careserve_L','Mgreserve_L','Sireserve_L', ...
    'rNc_L','rPc_L','rKc_L','rCac_L','rMgc_L','rSic_L', ...
    'NBLeaf_H','PARI_H','NBLI_H','NBLeaf_L','PARI_L','NBLI_L', ...
    'NupI_H','NupI_L','NavlI','RB_L','RB_H','NuLit_H', ...
    'NuLit_L','BLit','Nuptake_H','NH4uptake_H','NO3uptake_H','Puptake_H', ...
    'Kuptake_H','Cauptake_H','Mguptake_H','Siuptake_H','FNC_H','TexC_H', ...
    'TexN_H','TexP_H','TexK_H','TNIT_H','TPHO_H','TPOT_H', ...
    'SupN_H','SupP_H','SupK_H','TexCa_H','TexMg_H','TexSi_H', ...
    'TCAL_H','TMAG_H','TSIL_H','SupCa_H','SupMg_H','SupSi_H', ...
    'Nuptake_L','NH4uptake_L','NO3uptake_L','Puptake_L','Kuptake_L','Cauptake_L', ...
    'Mguptake_L','Siuptake_L','FNC_L','TexC_L','TexN_L','TexP_L', ...
    'TexK_L','TNIT_L','TPHO_L','TPOT_L','SupN_L','SupP_L', ...
    'SupK_L','TexCa_L','TexMg_L','TexSi_L','TCAL_L','TMAG_L', ...
    'TSIL_L','SupCa_L','SupMg_L','SupSi_L','ISOIL_H','ISOIL_L', ...
    'P','R_litter','R_microbe','R_litter_sur','R_ew','VOL', ...
    'N2flx','N2Oflx','NITflx','DENITflx','N2Oflx_nit','FertN_d', ...
    'DepN_d','PH_d','PH_pore_d','PH_bulk_d','Min_N','Min_P', ...
    'R_bacteria','RmycAM','RmycEM','Prod_B','Prod_F','BfixN', ...
    'LitFirEmi','LEAK_NH4','LEAK_NO3','LEAK_P','LEAK_K','LEAK_DOC', ...
    'LEAK_DON','LEAK_DOP','LEAK_Ca','LEAK_Mg','LEAK_Si','RexmyI', ...
    'AgrHarNut','SeedInput','Ccrown_t','AgePl_H','AgePl_L','Tden_H', ...
    'Tden_L','TBio_Ht','TBio_Lt'};
DAILY_STORE_NAMES = { ...
    'P','N2Oflx','N2Oflx_nit','N2flx','VOL','NITflx', ...
    'DENITflx','Min_N','BfixN','FertN_d','DepN_d','LEAK_NH4', ...
    'LEAK_NO3','LEAK_P','LEAK_K','LEAK_DOC','LEAK_DON','LEAK_DOP', ...
    'LEAK_Ca','LEAK_Mg','LEAK_Si','PH_bulk_d','PH_pore_d','R_litter', ...
    'R_microbe','R_ew','NPP_H','NPP_L','LAI_H','LAI_L', ...
    'Ccrown_t','RB_L'};
DAILY_STATE_ONLY = setdiff(DAILY_NAMES,DAILY_STORE_NAMES,'stable');
HOURLY_NAMES = { ...
    'V','O','Qi_out','WTR','POT','Tdp', ...
    'Vice','Oice','OF','OS','ZWT','Pr_sno', ...
    'Pr_liq','er','Ts','ra','r_soil','b_soil', ...
    'Tice','Rn','H','QE','Qv','Lpho', ...
    'EG','ESN','ESN_In','EWAT','EIn_urb','EIn_rock', ...
    'dw_SNO','G','SWE','SND','ros','In_SWE', ...
    'SP_wc','WR_SP','U_SWE','NIn_SWE','dQ','Qfm', ...
    't_sls','DQ','DT','In_urb','In_rock','SE_rock', ...
    'SE_urb','Csno','WIS','Lk','f','Rh', ...
    'Rd','NDVI','alp_soil','tau_sno','e_sno','ALB', ...
    'EK','dQVEG','TsV','HV','QEV','Cice', ...
    'Lk_wat','Lk_rock','EICE','WAT','ICE','ICE_D', ...
    'WR_IP','NIce','Cicew','IP_wc','Csnow','FROCK', ...
    'Imelt','Smelt','Tdamp','Gfin','In_Litter','ELitter', ...
    'Ws_under','r_litter','In_H','In_L','OH','OL', ...
    'Psi_s_H','Psi_s_L','Tdp_H','Tdp_L','rap_H','rap_L', ...
    'Dr_H','Dr_L','T_H','T_L','EIn_H','EIn_L', ...
    'rb_H','rb_L','An_H','Rdark_H','An_L','Rdark_L', ...
    'Ci_sunH','Ci_sunL','Ci_shdH','Ci_shdL','rs_sunH','rs_sunL', ...
    'rs_shdH','rs_shdL','Vx_H','Vl_H','Vx_L','Vl_L', ...
    'Psi_x_H','Psi_l_H','Psi_x_L','Psi_l_L','gsr_H','Jsx_H', ...
    'Jxl_H','Kleaf_H','Kx_H','gsr_L','Jsx_L','Jxl_L', ...
    'Kleaf_L','Kx_L','fapar_H','fapar_L','SIF_H','SIF_L', ...
    'Tdeb','q_runon','Qi_in'};
n_bio_layers = 0;
for c = 1:N_col
    n_bio_layers = max(n_bio_layers,find(PAR{c}.Bio_Zs > 0,1,'last'));
end
O_d = nan(NNd,n_bio_layers,N_col);
ZWT_d = nan(NNd,1,N_col);
Lk_d = nan(NNd,1,N_col);
f_d = nan(NNd,1,N_col);
q_runon_d = nan(NNd,1,N_col);
Rh_d = nan(NNd,1,N_col);
Rd_d = nan(NNd,1,N_col);
if field_grid
    [~,w_field] = Field_Aggregation_Areas(topo);
else
    w_field = reshape(A_cell/A_tot, 1, N_col);
end
FIELD_HOURLY_NAMES = {'Rn','H','QE','G','Ts','ET','T','EG','GPP','O','Tdp','Lk','Runoff'};
Field_hourly = struct('weights',w_field,'layers',n_bio_layers);
for k_fh = 1:numel(FIELD_HOURLY_NAMES)
    if any(strcmp(FIELD_HOURLY_NAMES{k_fh},{'O','Tdp'}))
        Field_hourly.(FIELD_HOURLY_NAMES{k_fh}) = nan(NN,n_bio_layers);
    else
        Field_hourly.(FIELD_HOURLY_NAMES{k_fh}) = nan(NN,1);
    end
end
%% Hour 1 is the initial state.
fh = Field_Hour_Values(1,N_col,cc,n_bio_layers,Ccrown_c, ...
    Rn,H,QE,G,Ts,EG,ELitter,ESN,ESN_In,EWAT,EICE,EIn_urb,EIn_rock, ...
    Lk,Rh,Rd,T_H,T_L,EIn_H,EIn_L,An_H,An_L,Rdark_H,Rdark_L,O,Tdp);
fhm = Field_Hourly_Means(fh,w_field);
for k_fh = 1:numel(FIELD_HOURLY_NAMES)
    Field_hourly.(FIELD_HOURLY_NAMES{k_fh})(1,:) = fhm.(FIELD_HOURLY_NAMES{k_fh});
end
ERW_daily = cell(N_col,1);
%% Per-cell radiation and sky view exist only with topography; the parallel
%% hourly loop slices these placeholders otherwise and never reads them.
no_topo_rad = zeros(N_col,1);
if topography_attached
    cell_svf = topo.cell_svf;
else
    cell_svf = no_topo_rad;
end
%% Today's row of every daily array before the daily step writes it: the
%% values the arrays were allocated with.
daily_row_template = struct();
daily_row_template.LAI_L = zeros(1,cc);
daily_row_template.B_L = zeros(1,cc,NCP);
daily_row_template.NPP_L = zeros(1,cc);
daily_row_template.Rg_L = zeros(1,cc);
daily_row_template.RA_L = zeros(1,cc);
daily_row_template.Rmc_L = zeros(1,cc);
daily_row_template.ANPP_L = zeros(1,cc);
daily_row_template.Rms_L = zeros(1,cc);
daily_row_template.Rmr_L = zeros(1,cc);
daily_row_template.PHE_S_L = zeros(1,cc);
daily_row_template.dflo_L = zeros(1,cc);
daily_row_template.AgeL_L = zeros(1,cc);
daily_row_template.e_rel_L = ones(1,cc);
daily_row_template.SAI_L = zeros(1,cc);
daily_row_template.hc_L = zeros(1,cc);
daily_row_template.e_relN_L = ones(1,cc);
daily_row_template.BA_L = zeros(1,cc);
daily_row_template.LAI_H = zeros(1,cc);
daily_row_template.B_H = zeros(1,cc,NCP);
daily_row_template.NPP_H = zeros(1,cc);
daily_row_template.Rg_H = zeros(1,cc);
daily_row_template.RA_H = zeros(1,cc);
daily_row_template.Rms_H = zeros(1,cc);
daily_row_template.Rmr_H = zeros(1,cc);
daily_row_template.ANPP_H = zeros(1,cc);
daily_row_template.PHE_S_H = zeros(1,cc);
daily_row_template.dflo_H = zeros(1,cc);
daily_row_template.Rmc_H = zeros(1,cc);
daily_row_template.AgeL_H = zeros(1,cc);
daily_row_template.e_rel_H = ones(1,cc);
daily_row_template.SAI_H = zeros(1,cc);
daily_row_template.hc_H = zeros(1,cc);
daily_row_template.e_relN_H = ones(1,cc);
daily_row_template.BA_H = zeros(1,cc);
daily_row_template.Rrootl_H = zeros(1,cc);
daily_row_template.Rrootl_L = zeros(1,cc);
daily_row_template.Bfac_dayH = ones(1,cc);
daily_row_template.Bfac_weekH = ones(1,cc);
daily_row_template.NPPI_H = zeros(1,cc);
daily_row_template.TdpI_H = zeros(1,cc);
daily_row_template.Bfac_dayL = ones(1,cc);
daily_row_template.Bfac_weekL = ones(1,cc);
daily_row_template.NPPI_L = zeros(1,cc);
daily_row_template.TdpI_L = zeros(1,cc);
daily_row_template.Sr_H = zeros(1,cc);
daily_row_template.Slf_H = zeros(1,cc);
daily_row_template.Sfr_H = zeros(1,cc);
daily_row_template.Swm_H = zeros(1,cc);
daily_row_template.Sll_H = zeros(1,cc);
daily_row_template.Sr_L = zeros(1,cc);
daily_row_template.Slf_L = zeros(1,cc);
daily_row_template.Sfr_L = zeros(1,cc);
daily_row_template.Swm_L = zeros(1,cc);
daily_row_template.Sll_L = zeros(1,cc);
daily_row_template.LAIdead_H = zeros(1,cc);
daily_row_template.LAIdead_L = zeros(1,cc);
daily_row_template.Rexmy_H = zeros(1,cc,3);
daily_row_template.AgeDL_H = zeros(1,cc);
daily_row_template.Nreserve_H = zeros(1,cc);
daily_row_template.Preserve_H = zeros(1,cc);
daily_row_template.Kreserve_H = zeros(1,cc);
daily_row_template.Careserve_H = zeros(1,cc);
daily_row_template.Mgreserve_H = zeros(1,cc);
daily_row_template.Sireserve_H = zeros(1,cc);
daily_row_template.rNc_H = zeros(1,cc);
daily_row_template.rPc_H = zeros(1,cc);
daily_row_template.rKc_H = zeros(1,cc);
daily_row_template.rCac_H = zeros(1,cc);
daily_row_template.rMgc_H = zeros(1,cc);
daily_row_template.rSic_H = zeros(1,cc);
daily_row_template.Rexmy_L = zeros(1,cc,3);
daily_row_template.AgeDL_L = zeros(1,cc);
daily_row_template.Nreserve_L = zeros(1,cc);
daily_row_template.Preserve_L = zeros(1,cc);
daily_row_template.Kreserve_L = zeros(1,cc);
daily_row_template.Careserve_L = zeros(1,cc);
daily_row_template.Mgreserve_L = zeros(1,cc);
daily_row_template.Sireserve_L = zeros(1,cc);
daily_row_template.rNc_L = zeros(1,cc);
daily_row_template.rPc_L = zeros(1,cc);
daily_row_template.rKc_L = zeros(1,cc);
daily_row_template.rCac_L = zeros(1,cc);
daily_row_template.rMgc_L = zeros(1,cc);
daily_row_template.rSic_L = zeros(1,cc);
daily_row_template.NBLeaf_H = zeros(1,cc);
daily_row_template.PARI_H = zeros(1,cc,3);
daily_row_template.NBLI_H = zeros(1,cc);
daily_row_template.NBLeaf_L = zeros(1,cc);
daily_row_template.PARI_L = zeros(1,cc,3);
daily_row_template.NBLI_L = zeros(1,cc);
daily_row_template.NupI_H = zeros(1,cc,6);
daily_row_template.NupI_L = zeros(1,cc,6);
daily_row_template.NavlI = zeros(1,6);
daily_row_template.RB_L = zeros(1,cc,7);
daily_row_template.RB_H = zeros(1,cc,7);
daily_row_template.NuLit_H = zeros(1,cc,6);
daily_row_template.NuLit_L = zeros(1,cc,6);
daily_row_template.BLit = zeros(1,cc);
daily_row_template.Nuptake_H = zeros(1,cc);
daily_row_template.NH4uptake_H = zeros(1,cc);
daily_row_template.NO3uptake_H = zeros(1,cc);
daily_row_template.Puptake_H = zeros(1,cc);
daily_row_template.Kuptake_H = zeros(1,cc);
daily_row_template.Cauptake_H = zeros(1,cc);
daily_row_template.Mguptake_H = zeros(1,cc);
daily_row_template.Siuptake_H = zeros(1,cc);
daily_row_template.FNC_H = zeros(1,cc);
daily_row_template.TexC_H = zeros(1,cc);
daily_row_template.TexN_H = zeros(1,cc);
daily_row_template.TexP_H = zeros(1,cc);
daily_row_template.TexK_H = zeros(1,cc);
daily_row_template.TNIT_H = zeros(1,cc);
daily_row_template.TPHO_H = zeros(1,cc);
daily_row_template.TPOT_H = zeros(1,cc);
daily_row_template.SupN_H = zeros(1,cc);
daily_row_template.SupP_H = zeros(1,cc);
daily_row_template.SupK_H = zeros(1,cc);
daily_row_template.TexCa_H = zeros(1,cc);
daily_row_template.TexMg_H = zeros(1,cc);
daily_row_template.TexSi_H = zeros(1,cc);
daily_row_template.TCAL_H = zeros(1,cc);
daily_row_template.TMAG_H = zeros(1,cc);
daily_row_template.TSIL_H = zeros(1,cc);
daily_row_template.SupCa_H = zeros(1,cc);
daily_row_template.SupMg_H = zeros(1,cc);
daily_row_template.SupSi_H = zeros(1,cc);
daily_row_template.Nuptake_L = zeros(1,cc);
daily_row_template.NH4uptake_L = zeros(1,cc);
daily_row_template.NO3uptake_L = zeros(1,cc);
daily_row_template.Puptake_L = zeros(1,cc);
daily_row_template.Kuptake_L = zeros(1,cc);
daily_row_template.Cauptake_L = zeros(1,cc);
daily_row_template.Mguptake_L = zeros(1,cc);
daily_row_template.Siuptake_L = zeros(1,cc);
daily_row_template.FNC_L = zeros(1,cc);
daily_row_template.TexC_L = zeros(1,cc);
daily_row_template.TexN_L = zeros(1,cc);
daily_row_template.TexP_L = zeros(1,cc);
daily_row_template.TexK_L = zeros(1,cc);
daily_row_template.TNIT_L = zeros(1,cc);
daily_row_template.TPHO_L = zeros(1,cc);
daily_row_template.TPOT_L = zeros(1,cc);
daily_row_template.SupN_L = zeros(1,cc);
daily_row_template.SupP_L = zeros(1,cc);
daily_row_template.SupK_L = zeros(1,cc);
daily_row_template.TexCa_L = zeros(1,cc);
daily_row_template.TexMg_L = zeros(1,cc);
daily_row_template.TexSi_L = zeros(1,cc);
daily_row_template.TCAL_L = zeros(1,cc);
daily_row_template.TMAG_L = zeros(1,cc);
daily_row_template.TSIL_L = zeros(1,cc);
daily_row_template.SupCa_L = zeros(1,cc);
daily_row_template.SupMg_L = zeros(1,cc);
daily_row_template.SupSi_L = zeros(1,cc);
daily_row_template.ISOIL_H = zeros(1,cc,27);
daily_row_template.ISOIL_L = zeros(1,cc,27);
daily_row_template.P = zeros(1,79);
daily_row_template.R_litter = zeros(1,1);
daily_row_template.R_microbe = zeros(1,1);
daily_row_template.R_litter_sur = zeros(1,1);
daily_row_template.R_ew = zeros(1,1);
daily_row_template.VOL = zeros(1,1);
daily_row_template.N2flx = zeros(1,1);
daily_row_template.N2Oflx = zeros(1,1);
daily_row_template.NITflx = zeros(1,1);
daily_row_template.DENITflx = zeros(1,1);
daily_row_template.N2Oflx_nit = zeros(1,1);
daily_row_template.FertN_d = zeros(1,1);
daily_row_template.DepN_d = zeros(1,1);
daily_row_template.PH_d = zeros(1,1);
daily_row_template.PH_pore_d = zeros(1,1);
daily_row_template.PH_bulk_d = zeros(1,1);
daily_row_template.Min_N = zeros(1,1);
daily_row_template.Min_P = zeros(1,1);
daily_row_template.R_bacteria = zeros(1,1);
daily_row_template.RmycAM = zeros(1,1);
daily_row_template.RmycEM = zeros(1,1);
daily_row_template.Prod_B = zeros(1,1);
daily_row_template.Prod_F = zeros(1,1);
daily_row_template.BfixN = zeros(1,1);
daily_row_template.LitFirEmi = zeros(1,2);
daily_row_template.LEAK_NH4 = zeros(1,1);
daily_row_template.LEAK_NO3 = zeros(1,1);
daily_row_template.LEAK_P = zeros(1,1);
daily_row_template.LEAK_K = zeros(1,1);
daily_row_template.LEAK_DOC = zeros(1,1);
daily_row_template.LEAK_DON = zeros(1,1);
daily_row_template.LEAK_DOP = zeros(1,1);
daily_row_template.LEAK_Ca = zeros(1,1);
daily_row_template.LEAK_Mg = zeros(1,1);
daily_row_template.LEAK_Si = zeros(1,1);
daily_row_template.RexmyI = zeros(1,3);
daily_row_template.AgrHarNut = zeros(1,6);
daily_row_template.SeedInput = zeros(1,7);
daily_row_template.Ccrown_t = ones(1,cc);
daily_row_template.AgePl_H = zeros(1,cc);
daily_row_template.AgePl_L = zeros(1,cc);
daily_row_template.Tden_H = zeros(1,cc);
daily_row_template.Tden_L = zeros(1,cc);
daily_row_template.TBio_Ht = zeros(1,cc);
daily_row_template.TBio_Lt = zeros(1,cc);

j = 1;
tic;

%% ========================================================================
%% CHECKPOINT / RESUME
%% ========================================================================
%% A full column-year costs tens of minutes, so a run that aborts late
%% (or a laptop that hibernates) should not have to start from day 1.
%% Every ERW_ckpt_minutes of wall time the whole workspace is written to
%% one rolling file, atomically (temp + rename), so an interrupted save
%% can never leave a half-written checkpoint behind.
%%
%% Resuming is only safe when the checkpoint was produced by the SAME
%% code and the same run configuration: a checkpoint carries the
%% trajectory that the old solver produced, and replaying it under
%% changed chemistry would silently splice two model versions together.
%% The fingerprint below covers the run identity, the horizon, and the
%% mtime+size of the files that decide the chemistry.  On a mismatch the
%% run starts clean and says so.  Set ERW_ckpt_force_resume = true to
%% override deliberately (useful when debugging a late-timestep failure
%% where the earlier trajectory is known to be unaffected).
if ~exist('ERW_ckpt_minutes','var') || isempty(ERW_ckpt_minutes)
    ERW_ckpt_minutes = 10;
end
if ~exist('ERW_ckpt_file','var') || isempty(ERW_ckpt_file)
    ERW_ckpt_file = '';
end
if ~exist('ERW_ckpt_force_resume','var') || isempty(ERW_ckpt_force_resume)
    ERW_ckpt_force_resume = false;
end
ERW_runtime_parts = cell(1,N_col);
for c = 1:N_col
    appelo_alpha_c = [];
    appelo_logKH_c = [];
    knmaxd_c = [];
    cec_dpH_c = [];
    ERW_runtime_param_c = ERW_PARAM_input;
    if isfield(ColCfg(c),'ERW_PARAM')
        ERW_runtime_param_c = ColCfg(c).ERW_PARAM;
    end
    if isfield(ColCfg(c),'appelo_alpha')
        appelo_alpha_c = ColCfg(c).appelo_alpha;
    end
    if isfield(ColCfg(c),'appelo_logKH')
        appelo_logKH_c = ColCfg(c).appelo_logKH;
    end
    if isfield(ColCfg(c),'KNmaxd')
        knmaxd_c = ColCfg(c).KNmaxd;
    end
    if isfield(ColCfg(c),'cec_dpH')
        cec_dpH_c = ColCfg(c).cec_dpH;
    end
    last_app_c = [];
    if isfield(ColCfg(c),'last_application_year')
        last_app_c = ColCfg(c).last_application_year;
    end
    %% Replaced biogeochemistry literals and the saturation exponent are
    %% part of the run's identity: a checkpoint of one member must not
    %% resume another.
    bgovr_c = '';
    if isfield(ColCfg(c),'biogeo_overrides') && isstruct(ColCfg(c).biogeo_overrides)
        bgovr_names = sort(fieldnames(ColCfg(c).biogeo_overrides));
        bgovr_pairs = cellfun(@(n) sprintf('%s=%s',n, ...
            mat2str(ColCfg(c).biogeo_overrides.(n),17)),bgovr_names, ...
            'UniformOutput',false);
        bgovr_c = strjoin(bgovr_pairs,',');
    end
    sexp_c = [];
    if isfield(ColCfg(c),'saturation_exponent')
        sexp_c = ColCfg(c).saturation_exponent;
    end
    pft_override_c = '';
    if isfield(ColCfg(c),'PFT_OVERRIDE') && ~isempty(ColCfg(c).PFT_OVERRIDE)
        pft_override_c = char(ColCfg(c).PFT_OVERRIDE);
    end
    ERW_runtime_parts{c} = sprintf( ...
        ['c%d:column_param=%s;pft=%s;apply_rock=%d;erw_param=%s;' ...
         'alpha=%s;logKH=%s;pcla=%s;psan=%s;knmaxd=%s;kbot=%s;cecdpH=%s;' ...
         'lastapp=%s;bgovr=%s;sexp=%s'], ...
        c,char(ColCfg(c).PARAM_IC),pft_override_c,logical(ColCfg(c).apply_rock), ...
        ERW_Runtime_Parameter_Identity(ERW_runtime_param_c), ...
        mat2str(appelo_alpha_c,17), ...
        mat2str(appelo_logKH_c,17),mat2str(PAR{c}.Pcla,17), ...
        mat2str(PAR{c}.Psan,17),mat2str(knmaxd_c,17), ...
        mat2str(PAR{c}.Kbot,17), ...
        mat2str(cec_dpH_c,17),mat2str(last_app_c,17),bgovr_c, ...
        mat2str(sexp_c,17));
    %% A per-cell MOD row (OPT_CELL_INPUTS) is the column's soil; the ERW
    %% row above is already per column. Off, the identity is unchanged.
    if OPT_CELL_INPUTS == 1 && isfield(ColCfg(c),'MOD_PARAM') && ...
            istable(ColCfg(c).MOD_PARAM) && height(ColCfg(c).MOD_PARAM) == 1
        ERW_runtime_parts{c} = sprintf('%s;mod_param=%s',ERW_runtime_parts{c}, ...
            ERW_Runtime_Parameter_Identity(ColCfg(c).MOD_PARAM));
    end
end
ERW_forcing_identity = sprintf('date_start=%.17g;date_end=%.17g',Date(1),Date(end));
ERW_topology_identity = sprintf([ ...
    'topology_areas=%s;topology_edge_a=%s;topology_edge_b=%s;' ...
    'topology_edge_width=%s;topology_edge_distance=%s;' ...
    'topology_edge_distance_a=%s;topology_edge_distance_b=%s;' ...
    'topology_edge_active=%s;topology_outer_width=%s'], ...
    mat2str(topo.cell_area_m2,17),mat2str(topo.edge_a,17), ...
    mat2str(topo.edge_b,17),mat2str(topo.edge_face_width_m,17), ...
    mat2str(topo.edge_face_distance_m,17), ...
    mat2str(topo.edge_face_distance_a_m,17), ...
    mat2str(topo.edge_face_distance_b_m,17), ...
    mat2str(topo.edge_active,17),mat2str(topo.outer_face_width_m,17));
%% Topography enters the checkpoint identity ONLY when attached, so every
%% existing flat checkpoint keeps its fingerprint and stays resumable.
if topography_attached
    %% Attach_Topography_To_Topology writes all of these together, so a
    %% missing field here means a malformed topology and should error.
    topo_field = @(name) mat2str(reshape(double(topo.(name)),[],1),17);
    ERW_topology_identity = sprintf(['%s;' ...
        'topology_source=%s;topology_cell_z=%s;topology_cell_slope=%s;' ...
        'topology_cell_aspect=%s;topology_cell_svf=%s;topology_cell_twi=%s;' ...
        'topology_cell_curv_prof=%s;topology_cell_curv_plan=%s;' ...
        'topology_edge_dz=%s;topology_edge_surface_distance=%s;' ...
        'topology_overland_ab=%s;topology_overland_ba=%s;topology_is_outlet=%s'], ...
        ERW_topology_identity, ...
        char(string(topo.topography_source)), ...
        topo_field('cell_z_m'),topo_field('cell_slope'), ...
        topo_field('cell_aspect_deg'),topo_field('cell_svf'), ...
        topo_field('cell_twi'),topo_field('cell_curv_prof'), ...
        topo_field('cell_curv_plan'),topo_field('edge_dz_m'), ...
        topo_field('edge_face_distance_surface_m'), ...
        topo_field('edge_overland_weight_ab'), ...
        topo_field('edge_overland_weight_ba'),topo_field('cell_is_outlet'));
end
%% Likewise the polygon overlap, which weights the field-origin tracers.
if field_grid
    ERW_topology_identity = sprintf('%s;topology_overlap=%s', ...
        ERW_topology_identity, ...
        mat2str(reshape(double(topo.cell_overlap_area_m2),[],1),17));
end
ERW_numerical_identity = strtrim(getenv('TC_NUMERICAL_IDENTITY'));
if isempty(ERW_numerical_identity); ERW_numerical_identity = 'unset'; end
%% The static TC_NUMERICAL_IDENTITY string is what the submission script
%% asked for, not what MATLAB actually loaded.  A new module build can pass
%% campaign validation while an older checkpoint still resumes, splicing
%% two numerical runtimes into one trajectory, so the identity carries the
%% interpreter, platform, BLAS and thread count that are really in force.
ERW_execution_identity = sprintf( ...
    ['parallel_workers=%d;numerical_identity=%s;matlab=%s;computer=%s;' ...
     'blas=%s;max_threads=%d'], ...
    parallel_workers,ERW_numerical_identity,version,computer, ...
    regexprep(strtrim(evalc('version -blas')),'\s+',' '),maxNumCompThreads);
%% A checkpoint holds one record layout; the other cannot resume it.
ERW_execution_identity = sprintf('%s;records=%s',ERW_execution_identity, ...
    ternary(FULL_RECORDS,'full','daily'));
ERW_ckpt_runtime_identity = strjoin( ...
    [{ERW_forcing_identity,ERW_topology_identity,ERW_execution_identity} ...
     ERW_runtime_parts],'|');
ERW_ckpt_fingerprint = ERW_Checkpoint_Fingerprint(NN, N_col, ...
    ERW_ckpt_runtime_identity, ...
    fullfile(cd, '*.m'), ...
    fullfile(cd, 'Parameters', '*.m'), ...
    fullfile(cd, 'BIOGEOCHEMISTRY_DYNAMIC3.m'), ...
    fullfile(cd, 'BIOGEO_UNIT.m'), ...
    fullfile(cd, 'Biogeochemistry_Parameter.m'), ...
    fullfile(cd, 'Biogeo_Cation_Partition_Rates.m'), ...
    fullfile(cd, 'Biogeo_Mobile_Water_Volume.m'), ...
    fullfile(cd, 'Biogeo_Nitrification_PH_Response.m'), ...
    fullfile(cd, 'CropHeightType.m'), ...
    fullfile(cd, 'Plant_Exports.m'), ...
    fullfile(cd, 'Plant_Reserve_Fraction.m'), ...
    fullfile(cd, 'Plant_Nutrient_Suppression.m'), ...
    fullfile(cd, 'Plant_Harvest_Reserve_Export.m'), ...
	fullfile(cd, 'Crop_Planting_Harvest.m'), ...
	fullfile(cd, 'Plant_Establishment_Element_Boundary.m'), ...
	fullfile(cd, 'Plant_Seed_Input_Accumulator.m'), ...
    fullfile(cd, 'ERW_function.m'), fullfile(cd, 'ERW_Daily_Step_2D.m'), ...
    fullfile(cd, 'ERW_Particle_History.m'), ...
    fullfile(cd, 'ERW_Daily_Record.m'), ...
    fullfile(cd, 'Field_Daily_From_Hourly.m'), ...
    fullfile(cd, 'Field_Hourly_Means.m'), ...
    fullfile(cd, 'Field_Hour_Values.m'), ...
    fullfile(cd, 'ERW_Initialize_TC_Cation_Pools.m'), ...
    fullfile(cd, 'ERW_Normalize_TC_Daily_Deltas.m'), ...
    fullfile(cd, 'ERW_Should_Apply_TC_Daily_Delta.m'), ...
    fullfile(cd, 'ERW_Return_Weathered_Phosphorus.m'), ...
    fullfile(cd, 'ERW_Appelo_Mass_Conserving_State.m'), ...
    fullfile(cd, 'ERW_Appelo_State_At_H.m'), ...
    fullfile(cd, 'ERW_CEC_At_H.m'), ...
    fullfile(cd, 'ERW_AlOH3_LogK.m'), ...
    fullfile(cd, 'ERW_Calcite_Kinetics.m'), ...
    fullfile(cd, 'ERW_Resolve_CEC.m'), ...
    fullfile(cd, 'ERW_Resolve_CEC_pH_Sensitivity.m'), ...
    fullfile(cd, 'ERW_Seed_Water_Loss.m'), ...
    fullfile(cd, 'ERW_Calcite_Reaction_Limit.m'), ...
    fullfile(cd, 'ERW_Feedstock_Pathway.m'), ...
    fullfile(cd, 'ERW_Particle_Mass_Fractions.m'), ...
    fullfile(cd, 'ERW_Mineral_DIC_Stoichiometry.m'), ...
    fullfile(cd, 'ERW_Mineral_Molar_Mass_Row.m'), ...
    fullfile(cd, 'ERW_Mineral_Saturation.m'), ...
    fullfile(cd, 'ERW_Total_Calcite_Solid.m'), ...
    fullfile(cd, 'ERW_pre_processing.m'), fullfile(cd, 'MAIN_FRAME_BG_ERW_2D.m'), ...
    fullfile(cd, 'Parameters', 'ERW_MOD_PARAM_excel.m'), ...
    fullfile(cd, 'Parameters', 'MOD_PARAM_Kantola_Maize_Soy.m'), ...
    fullfile(cd, 'Inputs', 'ERW_MOD_PARAM_parameters.xlsx'), ...
    fullfile(cd, 'Inputs', 'Data_Kantola.mat'), ...
    fullfile(cd, 'Inputs', 'Data_Kantola_BG.mat'), ...
    fullfile(cd, 'Inputs', 'Ca_Data.mat'), ...
    fullfile(cd, 'Inputs', 'All_deposition_data.mat'), ...
    fullfile(cd, 'Inputs', 'blueridge_CDF_data.json'), ...
    fullfile(cd, '..', 'run', 'Run_One_Subfield.m'), ...
    fullfile(cd, '..', 'run', 'Run_Matrix_12.m'), ...
    fullfile(cd, 'Resolve_Appelo_Constants.m'), ...
    fullfile(cd, 'Invert_Bulk_pH.m'), ...
    fullfile(cd, 'Porewater_to_Bulk_pH.m'));
i_resume = 2;
if ~isempty(ERW_ckpt_file) && isfile(ERW_ckpt_file)
    ckpt_probe = load(ERW_ckpt_file, 'ERW_ckpt_fingerprint', 'i');
    ckpt_ok = isfield(ckpt_probe,'ERW_ckpt_fingerprint') && ...
        strcmp(ckpt_probe.ERW_ckpt_fingerprint, ERW_ckpt_fingerprint);
    if ckpt_ok || ERW_ckpt_force_resume
        keep_ckpt_file = ERW_ckpt_file;
        keep_ckpt_minutes = ERW_ckpt_minutes;
        keep_ckpt_fp = ERW_ckpt_fingerprint;
        %% The caller's options describe this invocation, not the one that
        %% wrote the checkpoint.  The fingerprint pins the scientific
        %% configuration, but SaveResults and ResultDirectory are not in
        %% it, so restoring them wholesale would make a resumed run ignore
        %% where it was told to write and whether to write at all.
        if exist('opt','var'); keep_opt = opt; else; keep_opt = []; end
        %% The checkpoint carries the provenance of the segment that wrote
        %% it.  Restoring it wholesale would make the final MAT report the
        %% first Slurm job and node as if the run had never been
        %% interrupted, so the segments are accumulated instead.
        if exist('run_provenance','var'); keep_provenance = run_provenance;
        else;                             keep_provenance = [];
        end
        load(ERW_ckpt_file);
        %% A daily-store checkpoint carries the daily arrays outside the
        %% stores as their last two rows, all the run reads of them again.
        if exist('ckpt_daily_rows','var') && isstruct(ckpt_daily_rows)
            ckpt_names = fieldnames(ckpt_daily_rows);
            for k_ckpt = 1:numel(ckpt_names)
                ckpt_rows = ckpt_daily_rows.(ckpt_names{k_ckpt}); %#ok<NASGU>
                eval([ckpt_names{k_ckpt} '(ckpt_daily_rows_first:' ...
                    'ckpt_daily_rows_first+size(ckpt_rows,1)-1,:,:,:) = ckpt_rows;']);
            end
            clear ckpt_daily_rows ckpt_names ckpt_rows
        end
        ERW_ckpt_file = keep_ckpt_file;
        ERW_ckpt_minutes = keep_ckpt_minutes;
        ERW_ckpt_fingerprint = keep_ckpt_fp;
        if ~isempty(keep_opt); opt = keep_opt; end
        if ~isempty(keep_provenance)
            prior_segments = {};
            if exist('run_provenance','var') && isstruct(run_provenance)
                if isfield(run_provenance,'execution_segments')
                    prior_segments = run_provenance.execution_segments;
                    run_provenance = rmfield(run_provenance, ...
                        'execution_segments');
                end
                prior_segments{end+1} = run_provenance;
            end
            run_provenance = keep_provenance;
            run_provenance.execution_segments = prior_segments;
            run_provenance.resumed_from_iteration = i;
        end
        i_resume = i + 1;
        if ckpt_ok
            fprintf(['Resuming from checkpoint at iter %d / %d (%.1f%%)\n'], ...
                i, NN, 100*i/NN);
        else
            warning('ERW:CheckpointCodeMismatch', ...
                ['Resuming from a checkpoint written by DIFFERENT code ' ...
                 '(ERW_ckpt_force_resume = true). The trajectory before ' ...
                 'iter %d came from the previous version.'], i);
        end
        tic;
    else
        fprintf(['Checkpoint found but code/config changed; starting clean.\n']);
    end
end
ERW_ckpt_last_toc = toc;

%% ========================================================================
%% TIME LOOP
%% ========================================================================
%% The hourly arrays hold the last NH hours (all NN with full records):
%% hour h sits in row mod(h-1,NH)+1, so si is hour i, sp hour i-1 and pd the
%% 24 hours before i. The daily arrays keep every day.
for i = i_resume:NN
    if (mod(i,1000) == 0) || (i == 2)
        fprintf('Iter %d / %d   (%.1f%%)\n', i, NN, 100*i/NN);
        %% Pool progress, printed to the run log only. Read-only: it
        %% touches no state and cannot alter the simulation. Column 1 is
        %% enough to spot a pool heading to zero or running away.
        if OPT_SoilBiogeochemistry == 1 && j > 1
            fprintf(['  POOLS  SOMC %.4g  POC %.4g  MOC %.4g  DOC %.4g  ' ...
                     'Bact %.4g  Fung %.4g  AM %.4g  SOMN %.4g  ' ...
                     'NH4 %.4g  NO3 %.4g  Pmin %.4g  Ca %.4g\n'], ...
                P(j-1,9,1), P(j-1,10,1), P(j-1,11,1), ...
                P(j-1,12,1)+P(j-1,13,1), P(j-1,18,1), P(j-1,19,1), ...
                P(j-1,20,1), P(j-1,26,1), P(j-1,31,1), P(j-1,32,1), ...
                P(j-1,43,1), P(j-1,60,1));
        end
    end
    pdind = max(1,i-24):i-1;
    si = mod(i-1,NH)+1;
    sp = mod(i-2,NH)+1;
    pd = mod(pdind-1,NH)+1;
    h_S = NaN;

    %% ------------------------------------------------------------------
    %% DAILY BLOCK  --  one square per worker
    %% ------------------------------------------------------------------
    %% Every square's day is independent of the others'. Each worker reads
    %% yesterday's rows (D1) and the 24 hours its square averaged (W1),
    %% builds today's row in o, and the rows are written back after the
    %% loop. Carried per-square state (Ccrown_c, ZR95_*_c, Rf*_Zs_c,
    %% Cbare_c, B_IO_c, Bam_c, Bem_c, Tdp_bio_c, VSUM_c, VTSUM_c) is read
    %% and written by its own square only.
    if (Datam(i,4) == 1)
        j = j+1;
        [jDay_i] = JULIAN_DAY(Datam(i,:));
        [~,~,~,~,~,L_day_j] = SetSunVariables(Datam(i,:),DeltaGMT,Lon,Lat,t_bef,t_aft);
        jDay(j) = jDay_i; L_day(j) = L_day_j; %#ok<SAGROW>

        %% Daily stores of the hourly records: the hours this daily row
        %% averages (row 2 is the first hour alone).
        O_d(j,:,:)       = Field_Daily_From_Hourly(O(pd,1:n_bio_layers,:),'mean',dth);
        ZWT_d(j,:,:)     = Field_Daily_From_Hourly(ZWT(pd,:,:),'mean',dth);
        Lk_d(j,:,:)      = Field_Daily_From_Hourly(Lk(pd,:,:),'rate_total',dth);
        f_d(j,:,:)       = Field_Daily_From_Hourly(f(pd,:,:),'rate_total',dth);
        q_runon_d(j,:,:) = Field_Daily_From_Hourly(q_runon(pd,:,:),'rate_total',dth);
        Rh_d(j,:,:)      = Field_Daily_From_Hourly(Rh(pd,:,:),'total',dth);
        Rd_d(j,:,:)      = Field_Daily_From_Hourly(Rd(pd,:,:),'total',dth);

        %% Yesterday's rows and the 24 hours each square averages.
        D1 = struct();
        D1.AgeDL_H = AgeDL_H(j-1,:,:);
        D1.AgeDL_L = AgeDL_L(j-1,:,:);
        D1.AgeL_H = AgeL_H(j-1,:,:);
        D1.AgeL_L = AgeL_L(j-1,:,:);
        D1.B_H = B_H(j-1,:,:,:);
        D1.B_L = B_L(j-1,:,:,:);
        D1.Bfac_weekH = Bfac_weekH(j-1,:,:);
        D1.Bfac_weekL = Bfac_weekL(j-1,:,:);
        D1.Careserve_H = Careserve_H(j-1,:,:);
        D1.Careserve_L = Careserve_L(j-1,:,:);
        D1.FNC_H = FNC_H(j-1,:,:);
        D1.FNC_L = FNC_L(j-1,:,:);
        D1.ISOIL_H = ISOIL_H(j-1,:,:,:);
        D1.ISOIL_L = ISOIL_L(j-1,:,:,:);
        D1.Kreserve_H = Kreserve_H(j-1,:,:);
        D1.Kreserve_L = Kreserve_L(j-1,:,:);
        D1.LAI_H = LAI_H(j-1,:,:);
        D1.LAI_L = LAI_L(j-1,:,:);
        D1.LAIdead_H = LAIdead_H(j-1,:,:);
        D1.LAIdead_L = LAIdead_L(j-1,:,:);
        D1.Mgreserve_H = Mgreserve_H(j-1,:,:);
        D1.Mgreserve_L = Mgreserve_L(j-1,:,:);
        D1.NBLI_H = NBLI_H(j-1,:,:);
        D1.NBLI_L = NBLI_L(j-1,:,:);
        D1.NBLeaf_H = NBLeaf_H(j-1,:,:);
        D1.NBLeaf_L = NBLeaf_L(j-1,:,:);
        D1.NPPI_H = NPPI_H(j-1,:,:);
        D1.NPPI_L = NPPI_L(j-1,:,:);
        D1.NPP_H = NPP_H(j-1,:,:);
        D1.NPP_L = NPP_L(j-1,:,:);
        D1.NavlI = NavlI(j-1,:,:);
        D1.Nreserve_H = Nreserve_H(j-1,:,:);
        D1.Nreserve_L = Nreserve_L(j-1,:,:);
        D1.NuLit_H = NuLit_H(j-1,:,:,:);
        D1.NuLit_L = NuLit_L(j-1,:,:,:);
        D1.NupI_H = NupI_H(j-1,:,:,:);
        D1.NupI_L = NupI_L(j-1,:,:,:);
        D1.P = P(j-1,:,:);
        D1.PARI_H = PARI_H(j-1,:,:,:);
        D1.PARI_L = PARI_L(j-1,:,:,:);
        D1.PHE_S_H = PHE_S_H(j-1,:,:);
        D1.PHE_S_L = PHE_S_L(j-1,:,:);
        D1.Preserve_H = Preserve_H(j-1,:,:);
        D1.Preserve_L = Preserve_L(j-1,:,:);
        D1.RexmyI = RexmyI(j-1,:,:);
        D1.Rexmy_H = Rexmy_H(j-1,:,:,:);
        D1.Rexmy_L = Rexmy_L(j-1,:,:,:);
        D1.SAI_H = SAI_H(j-1,:,:);
        D1.SAI_L = SAI_L(j-1,:,:);
        D1.Sireserve_H = Sireserve_H(j-1,:,:);
        D1.Sireserve_L = Sireserve_L(j-1,:,:);
        D1.SupCa_H = SupCa_H(j-1,:,:);
        D1.SupCa_L = SupCa_L(j-1,:,:);
        D1.SupK_H = SupK_H(j-1,:,:);
        D1.SupK_L = SupK_L(j-1,:,:);
        D1.SupMg_H = SupMg_H(j-1,:,:);
        D1.SupMg_L = SupMg_L(j-1,:,:);
        D1.SupN_H = SupN_H(j-1,:,:);
        D1.SupN_L = SupN_L(j-1,:,:);
        D1.SupP_H = SupP_H(j-1,:,:);
        D1.SupP_L = SupP_L(j-1,:,:);
        D1.SupSi_H = SupSi_H(j-1,:,:);
        D1.SupSi_L = SupSi_L(j-1,:,:);
        D1.TBio_Ht = TBio_Ht(j-1,:,:);
        D1.TBio_Lt = TBio_Lt(j-1,:,:);
        D1.TdpI_H = TdpI_H(j-1,:,:);
        D1.TdpI_L = TdpI_L(j-1,:,:);
        D1.dflo_H = dflo_H(j-1,:,:);
        D1.dflo_L = dflo_L(j-1,:,:);
        D1.hc_H = hc_H(j-1,:,:);
        D1.hc_L = hc_L(j-1,:,:);
        W1 = struct();
        W1.An_H = An_H(pd,:,:);
        W1.An_L = An_L(pd,:,:);
        W1.Lk = Lk(pd,:,:);
        W1.O = O(pd,:,:);
        W1.Psi_l_H = Psi_l_H(pd,:,:);
        W1.Psi_l_L = Psi_l_L(pd,:,:);
        W1.Psi_x_H = Psi_x_H(pd,:,:);
        W1.Psi_x_L = Psi_x_L(pd,:,:);
        W1.Rdark_H = Rdark_H(pd,:,:);
        W1.Rdark_L = Rdark_L(pd,:,:);
        W1.T_H = T_H(pd,:,:);
        W1.T_L = T_L(pd,:,:);
        W1.Tdp = Tdp(pd,:,:);
        W1.Tdp_H = Tdp_H(pd,:,:);
        W1.Tdp_L = Tdp_L(pd,:,:);
        W1.V = V(pd,:,:);
        %% The hour's and day's forcing as the workers read it.
        Datam_i = Datam(i,:);
        Ta_pd = Ta(pdind);
        PAR_pd = PARB(pdind)+PARD(pdind);
        dout = cell(N_col,1);
        carbon_n = cell(N_col,1);
        parfor (c = 1:N_col,parallel_workers)
            p = PAR{c};
            o = daily_row_template;
            ManIH_c = ManIH(:,c); ManIL_c = ManIL(:,c);
            Ccrown = Ccrown_c(c,:); Cbare = Cbare_c(c);
            ZR95_H = ZR95_H_c(c,:); ZR95_L = ZR95_L_c(c,:);
            RfH_Zs = RfH_Zs_c{c};   RfL_Zs = RfL_Zs_c{c};
            B_IO_c_now = B_IO_c{c};

            %% Rotation fertiliser rule for cropped cells
            if isempty(char(ColCfg(c).PFT_OVERRIDE))
                if ismember(Datam_i(1),B_IO_c_now.FertN_years)
                    B_IO_c_now.FertN(81:100) = 1.01; %% 202 kg N ha-1 (Kantola 2022; 2011/2014 were 168, outside 2016-2020)
                else
                    B_IO_c_now.FertN = zeros(1,366);
                end
                B_IO_c{c} = B_IO_c_now;
            end

            [Se_bio,Se_fc,Psi_bio,Tdp_bio,VSUM,VTSUM,f_anaer_bio,WFPS_bio] = Biogeo_environment( ...
                W1.Tdp(:,:,c),W1.O(:,:,c),W1.V(:,:,c),p.Soil_Param,p.Phy,p.SPAR,p.Bio_Zs,p.Zbio);
            Tdp_bio_c(c) = Tdp_bio; VSUM_c(c) = VSUM; VTSUM_c(c) = VTSUM;

            if OPT_SoilBiogeochemistry == 1
                IS    = Ccrown*squeeze(D1.ISOIL_L(1,:,:,c)) + Ccrown*squeeze(D1.ISOIL_H(1,:,:,c));
                REXMY = Ccrown*squeeze(D1.Rexmy_L(1,:,:,c)) + Ccrown*squeeze(D1.Rexmy_H(1,:,:,c));
                FireA = 1*((sum(ManIH_c==-5) + sum(ManIL_c==-5)) > 0);
                if sum(Ccrown) > 0
                    OPT_BfixN = p.ParEx_L(Ccrown>0).bfix;
                else
                    OPT_BfixN = 0;
                end
                if ~(isfinite(PH_c(c)) && PH_c(c) > 0 && PH_c(c) < 14)
                    error('BIOGEO:InvalidPorePH', ...
                        'Column %d pore-water pH must be finite and within (0,14).',c);
                end
                if ~(isfinite(PH_bulk_c(c)) && PH_bulk_c(c) > 0 && PH_bulk_c(c) < 14)
                    error('BIOGEO:InvalidBulkPH', ...
                        'Column %d microbial bulk-soil pH must be finite and within (0,14).',c);
                end
                o.PH_d(1) = PH_bulk_c(c);
                o.PH_pore_d(1) = PH_c(c);
                o.PH_bulk_d(1) = PH_bulk_c(c);

                [o.P(1,:),o.LEAK_NH4(1),o.LEAK_NO3(1),o.LEAK_P(1),o.LEAK_K(1), ...
                 o.LEAK_DOC(1),o.LEAK_DON(1),o.LEAK_DOP(1), ...
                 o.LEAK_Ca(1),o.LEAK_Mg(1),o.LEAK_Si(1), ...
                 o.Nuptake_H(1,:),o.NH4uptake_H(1,:),o.NO3uptake_H(1,:), ...
                 o.Puptake_H(1,:),o.Kuptake_H(1,:), ...
                 o.Nuptake_L(1,:),o.NH4uptake_L(1,:),o.NO3uptake_L(1,:), ...
                 o.Puptake_L(1,:),o.Kuptake_L(1,:), ...
                 o.Cauptake_H(1,:),o.Mguptake_H(1,:),o.Siuptake_H(1,:), ...
                 o.Cauptake_L(1,:),o.Mguptake_L(1,:),o.Siuptake_L(1,:), ...
                 o.RexmyI(1,:), ...
                 o.R_litter(1),o.R_microbe(1),o.R_litter_sur(1),o.R_ew(1), ...
                 o.VOL(1),o.N2flx(1),o.Min_N(1),o.Min_P(1), ...
                 o.R_bacteria(1),o.RmycAM(1),o.RmycEM(1),o.Prod_B(1),o.Prod_F(1), ...
                 o.BfixN(1),o.NavlI(1,:),o.LitFirEmi(1,:),o.N2Oflx(1), ...
                 o.FertN_d(1),o.DepN_d(1),o.NITflx(1),o.DENITflx(1),o.N2Oflx_nit(1)] = ...
                 BIOGEO_UNIT(D1.P(1,:,c),IS,p.Zbio,sum(p.Bio_Zs.*p.rsd), ...
                    PH_bulk_c(c),Tdp_bio,mean(Ta_pd),Psi_bio,Se_bio,Se_fc,VSUM,VTSUM, ...
                    Ccrown,p.Bio_Zs,RfH_Zs,RfL_Zs, ...
                    sum(W1.Lk(:,1,c)),sum(W1.T_H(:,:,c),1),sum(W1.T_L(:,:,c),1), ...
                    D1.B_H(1,:,3,c),D1.B_L(1,:,3,c),D1.LAI_H(1,:,c),D1.LAI_L(1,:,c), ...
                    D1.SupN_H(1,:,c),D1.SupP_H(1,:,c),D1.SupK_H(1,:,c), ...
                    D1.SupN_L(1,:,c),D1.SupP_L(1,:,c),D1.SupK_L(1,:,c), ...
                    D1.SupCa_H(1,:,c),D1.SupMg_H(1,:,c),D1.SupSi_H(1,:,c), ...
                    D1.SupCa_L(1,:,c),D1.SupMg_L(1,:,c),D1.SupSi_L(1,:,c), ...
                    REXMY,D1.RexmyI(1,:,c),p.ExEM,D1.NavlI(1,:,c),p.Pcla,p.Psan, ...
                    B_IO_c_now,Datam_i,FireA,0, ...
                    OPT_ERW,OPT_BfixN,f_anaer_bio,WFPS_bio, ...
                    biogeo_exchange_scheme{c},biogeo_KNmaxd{c}, ...
                    biogeo_overrides{c});

                o.BLit(1,:) = 0.002*sum(o.P(1,1:5))*Ccrown';
                Bam_c(c) = o.P(1,20);
                Bem_c(c) = o.P(1,21);
            else
                o.BLit(1,:) = 0;
                o.NavlI(1,:) = 0; o.P(1,:) = 0;
                Bam_c(c) = 0; Bem_c(c) = 0;
            end

            %% ---- Vegetation, per PFT slot -----------------------------
            for k = 1:length(Ccrown)
                if ZR95_H(k) > 0
                    [o.LAI_H(1,k),o.B_H(1,k,:),o.NPP_H(1,k),o.ANPP_H(1,k),o.Rg_H(1,k), ...
                     o.RA_H(1,k),o.Rms_H(1,k),o.Rmr_H(1,k),o.Rmc_H(1,k),o.PHE_S_H(1,k), ...
                     o.dflo_H(1,k),o.AgeL_H(1,k),o.e_rel_H(1,k),o.e_relN_H(1,k), ...
                     o.LAIdead_H(1,k),o.NBLeaf_H(1,k), ...
                     o.Sr_H(1,k),o.Slf_H(1,k),o.Sfr_H(1,k),o.Swm_H(1,k),o.Sll_H(1,k), ...
                     o.Rexmy_H(1,k,:),o.Rrootl_H(1,k), ...
                     o.AgeDL_H(1,k),o.Bfac_dayH(1,k),o.Bfac_weekH(1,k),o.NPPI_H(1,k), ...
                     o.TdpI_H(1,k),o.NupI_H(1,k,:),o.PARI_H(1,k,:),o.NBLI_H(1,k), ...
                     o.RB_H(1,k,:),o.FNC_H(1,k), ...
                     o.Nreserve_H(1,k),o.Preserve_H(1,k),o.Kreserve_H(1,k), ...
                     o.Careserve_H(1,k),o.Mgreserve_H(1,k),o.Sireserve_H(1,k), ...
                     o.rNc_H(1,k),o.rPc_H(1,k),o.rKc_H(1,k), ...
                     o.rCac_H(1,k),o.rMgc_H(1,k),o.rSic_H(1,k),ManIH_c(k),SeedInput_cc] = ...
                     VEGGIE_UNIT( ...
                        D1.B_H(1,k,:,c),D1.PHE_S_H(1,k,c),D1.dflo_H(1,k,c), ...
                        D1.AgeL_H(1,k,c),D1.AgeDL_H(1,k,c), ...
                        Ta_pd,W1.Tdp_H(:,k,c),PAR_pd, ...
                        W1.Psi_x_H(:,k,c),W1.Psi_l_H(:,k,c), ...
                        W1.An_H(:,k,c),W1.Rdark_H(:,k,c),D1.NPP_H(1,k,c),jDay_i,Datam_i, ...
                        D1.NPPI_H(1,k,c),D1.TdpI_H(1,k,c),D1.Bfac_weekH(1,k,c), ...
                        D1.NupI_H(1,k,:,c),o.NavlI(1,:),D1.PARI_H(1,k,:,c), ...
                        D1.NBLI_H(1,k,c),D1.NBLeaf_H(1,k,c), ...
                        L_day_j,p.Lmax_day,p.VegH_Param_Dyn,k, ...
                        D1.Nreserve_H(1,k,c),D1.Preserve_H(1,k,c),D1.Kreserve_H(1,k,c), ...
                        D1.Careserve_H(1,k,c),D1.Mgreserve_H(1,k,c),D1.Sireserve_H(1,k,c), ...
                        o.Nuptake_H(1,k),o.Puptake_H(1,k),o.Kuptake_H(1,k), ...
                        o.Cauptake_H(1,k),o.Mguptake_H(1,k),o.Siuptake_H(1,k), ...
                        D1.FNC_H(1,k,c),Se_bio,Tdp_bio, ...
                        p.ParEx_H(k),p.ExEM,Bam_c(c),Bem_c(c),p.Mpar_H(k), ...
                        D1.TBio_Ht(1,k,c),OPT_EnvLimitGrowth,OPT_VCA,OPT_VD,OPT_SoilBiogeochemistry);

                    [o.TexC_H(1,k),o.TexN_H(1,k),o.TexP_H(1,k),o.TexK_H(1,k), ...
                     o.TexCa_H(1,k),o.TexMg_H(1,k),o.TexSi_H(1,k), ...
                     o.TNIT_H(1,k),o.TPHO_H(1,k),o.TPOT_H(1,k), ...
                     o.TCAL_H(1,k),o.TMAG_H(1,k),o.TSIL_H(1,k), ...
                     o.NuLit_H(1,k,:),o.Nreserve_H(1,k),o.Preserve_H(1,k),o.Kreserve_H(1,k), ...
                     o.Careserve_H(1,k),o.Mgreserve_H(1,k),o.Sireserve_H(1,k), ...
                     o.SupN_H(1,k),o.SupP_H(1,k),o.SupK_H(1,k), ...
                     o.SupCa_H(1,k),o.SupMg_H(1,k),o.SupSi_H(1,k), ...
                     o.ISOIL_H(1,k,:)] = Plant_Exports(o.B_H(1,k,:),D1.B_H(1,k,:,c), ...
                        D1.NuLit_H(1,k,:,c),o.Slf_H(1,k),o.Sfr_H(1,k),o.Swm_H(1,k), ...
                        o.Sll_H(1,k),o.Sr_H(1,k), ...
                        o.Rexmy_H(1,k,:),p.Stoich_H(k),p.Mpar_H(k),p.fab_H(k),p.fbe_H(k), ...
                        o.RB_H(1,k,:),o.Nreserve_H(1,k),o.Preserve_H(1,k),o.Kreserve_H(1,k), ...
                        o.Careserve_H(1,k),o.Mgreserve_H(1,k),o.Sireserve_H(1,k), ...
                        o.rNc_H(1,k),o.rPc_H(1,k),o.rKc_H(1,k), ...
                        o.rCac_H(1,k),o.rMgc_H(1,k),o.rSic_H(1,k), ...
                        ManIH_c(k),OPT_SoilBiogeochemistry,id_location);

                    o.SAI_H(1,k) = D1.SAI_H(1,k,c);
                    if p.aSE_H(k) == 2
                        o.hc_H(1,k) = GrassHeight(o.LAI_H(1,k),o.LAIdead_H(1,k));
                    elseif p.aSE_H(k) == 5
                        [o.hc_H(1,k),o.SAI_H(1,k),o.B_H(1,:,:),Ccrown,ZR95_H, ...
                         o.Nreserve_H(1,:),o.Preserve_H(1,:),o.Kreserve_H(1,:), ...
                         o.Careserve_H(1,:),o.Mgreserve_H(1,:),o.Sireserve_H(1,:), ...
                         AgrHarNut_cc,RfH_Zs] = CropHeightType( ...
                            o.LAI_H(1,k),o.LAIdead_H(1,k),k,ZR95_H,o.B_H(1,:,:), ...
                            p.Zs,p.CASE_ROOT,Ccrown,Cbare+p.Cwat+p.Curb+p.Crock, ...
                            o.Nreserve_H(1,:),o.Preserve_H(1,:),o.Kreserve_H(1,:), ...
                            o.Careserve_H(1,:),o.Mgreserve_H(1,:),o.Sireserve_H(1,:), ...
                            ManIH_c,p.Mpar_H,p.VegH_Param_Dyn,OPT_SoilBiogeochemistry);
                        o.AgrHarNut(1,:) = Plant_Harvest_Reserve_Export( ...
                            o.AgrHarNut(1,:),AgrHarNut_cc,Ccrown(k));
                        o.SeedInput(1,:) = Plant_Seed_Input_Accumulator( ...
                            o.SeedInput(1,:),SeedInput_cc,Ccrown(k));
                    else
                        o.hc_H(1,k) = D1.hc_H(1,k,c);
                    end
                else
                    o.LAI_H(1,k)=0; o.B_H(1,k,:)=0; o.NPP_H(1,k)=0; o.ANPP_H(1,k)=0;
                    o.Rg_H(1,k)=0; o.RA_H(1,k)=0; o.Rms_H(1,k)=0; o.Rmr_H(1,k)=0;
                    o.Rmc_H(1,k)=0; o.PHE_S_H(1,k)=0; o.dflo_H(1,k)=0; o.AgeL_H(1,k)=0;
                    o.e_rel_H(1,k)=0; o.e_relN_H(1,k)=0; o.SAI_H(1,k)=0; o.hc_H(1,k)=0;
                    o.LAIdead_H(1,k)=D1.LAIdead_H(1,k,c);
                    o.Sr_H(1,k)=0; o.Slf_H(1,k)=0; o.Sfr_H(1,k)=0; o.Swm_H(1,k)=0;
                    o.Sll_H(1,k)=0; o.Rrootl_H(1,k)=0; o.Rexmy_H(1,k,:)=0;
                    o.AgeDL_H(1,k)=0; o.FNC_H(1,k)=1;
                    o.Nreserve_H(1,k)=0; o.Preserve_H(1,k)=0; o.Kreserve_H(1,k)=0;
                    o.Careserve_H(1,k)=0; o.Mgreserve_H(1,k)=0; o.Sireserve_H(1,k)=0;
                    o.rNc_H(1,k)=1; o.rPc_H(1,k)=1; o.rKc_H(1,k)=1;
                    o.rCac_H(1,k)=1; o.rMgc_H(1,k)=1; o.rSic_H(1,k)=1;
                    o.Bfac_dayH(1,k)=0; o.Bfac_weekH(1,k)=0; o.NPPI_H(1,k)=0;
                    o.TdpI_H(1,k)=0; o.NupI_H(1,k,:)=0; o.RB_H(1,k,:)=0;
                    o.TexC_H(1,k)=0; o.TexN_H(1,k)=0; o.TexP_H(1,k)=0; o.TexK_H(1,k)=0;
                    o.TexCa_H(1,k)=0; o.TexMg_H(1,k)=0; o.TexSi_H(1,k)=0;
                    o.TNIT_H(1,k)=0; o.TPHO_H(1,k)=0; o.TPOT_H(1,k)=0;
                    o.TCAL_H(1,k)=0; o.TMAG_H(1,k)=0; o.TSIL_H(1,k)=0;
                    o.ISOIL_H(1,k,:)=0; o.NuLit_H(1,k,:)=0;
                end

                if ZR95_L(k) > 0
                    [o.LAI_L(1,k),o.B_L(1,k,:),o.NPP_L(1,k),o.ANPP_L(1,k),o.Rg_L(1,k), ...
                     o.RA_L(1,k),o.Rms_L(1,k),o.Rmr_L(1,k),o.Rmc_L(1,k),o.PHE_S_L(1,k), ...
                     o.dflo_L(1,k),o.AgeL_L(1,k),o.e_rel_L(1,k),o.e_relN_L(1,k), ...
                     o.LAIdead_L(1,k),o.NBLeaf_L(1,k), ...
                     o.Sr_L(1,k),o.Slf_L(1,k),o.Sfr_L(1,k),o.Swm_L(1,k),o.Sll_L(1,k), ...
                     o.Rexmy_L(1,k,:),o.Rrootl_L(1,k), ...
                     o.AgeDL_L(1,k),o.Bfac_dayL(1,k),o.Bfac_weekL(1,k),o.NPPI_L(1,k), ...
                     o.TdpI_L(1,k),o.NupI_L(1,k,:),o.PARI_L(1,k,:),o.NBLI_L(1,k), ...
                     o.RB_L(1,k,:),o.FNC_L(1,k), ...
                     o.Nreserve_L(1,k),o.Preserve_L(1,k),o.Kreserve_L(1,k), ...
                     o.Careserve_L(1,k),o.Mgreserve_L(1,k),o.Sireserve_L(1,k), ...
                     o.rNc_L(1,k),o.rPc_L(1,k),o.rKc_L(1,k), ...
                     o.rCac_L(1,k),o.rMgc_L(1,k),o.rSic_L(1,k),ManIL_c(k),SeedInput_cc] = ...
                     VEGGIE_UNIT( ...
                        D1.B_L(1,k,:,c),D1.PHE_S_L(1,k,c),D1.dflo_L(1,k,c), ...
                        D1.AgeL_L(1,k,c),D1.AgeDL_L(1,k,c), ...
                        Ta_pd,W1.Tdp_L(:,k,c),PAR_pd, ...
                        W1.Psi_x_L(:,k,c),W1.Psi_l_L(:,k,c), ...
                        W1.An_L(:,k,c),W1.Rdark_L(:,k,c),D1.NPP_L(1,k,c),jDay_i,Datam_i, ...
                        D1.NPPI_L(1,k,c),D1.TdpI_L(1,k,c),D1.Bfac_weekL(1,k,c), ...
                        D1.NupI_L(1,k,:,c),o.NavlI(1,:),D1.PARI_L(1,k,:,c), ...
                        D1.NBLI_L(1,k,c),D1.NBLeaf_L(1,k,c), ...
                        L_day_j,p.Lmax_day,p.VegL_Param_Dyn,k, ...
                        D1.Nreserve_L(1,k,c),D1.Preserve_L(1,k,c),D1.Kreserve_L(1,k,c), ...
                        D1.Careserve_L(1,k,c),D1.Mgreserve_L(1,k,c),D1.Sireserve_L(1,k,c), ...
                        o.Nuptake_L(1,k),o.Puptake_L(1,k),o.Kuptake_L(1,k), ...
                        o.Cauptake_L(1,k),o.Mguptake_L(1,k),o.Siuptake_L(1,k), ...
                        D1.FNC_L(1,k,c),Se_bio,Tdp_bio, ...
                        p.ParEx_L(k),p.ExEM,Bam_c(c),Bem_c(c),p.Mpar_L(k), ...
                        D1.TBio_Lt(1,k,c),OPT_EnvLimitGrowth,OPT_VCA,OPT_VD,OPT_SoilBiogeochemistry);

                    [o.TexC_L(1,k),o.TexN_L(1,k),o.TexP_L(1,k),o.TexK_L(1,k), ...
                     o.TexCa_L(1,k),o.TexMg_L(1,k),o.TexSi_L(1,k), ...
                     o.TNIT_L(1,k),o.TPHO_L(1,k),o.TPOT_L(1,k), ...
                     o.TCAL_L(1,k),o.TMAG_L(1,k),o.TSIL_L(1,k), ...
                     o.NuLit_L(1,k,:),o.Nreserve_L(1,k),o.Preserve_L(1,k),o.Kreserve_L(1,k), ...
                     o.Careserve_L(1,k),o.Mgreserve_L(1,k),o.Sireserve_L(1,k), ...
                     o.SupN_L(1,k),o.SupP_L(1,k),o.SupK_L(1,k), ...
                     o.SupCa_L(1,k),o.SupMg_L(1,k),o.SupSi_L(1,k), ...
                     o.ISOIL_L(1,k,:)] = Plant_Exports(o.B_L(1,k,:),D1.B_L(1,k,:,c), ...
                        D1.NuLit_L(1,k,:,c),o.Slf_L(1,k),o.Sfr_L(1,k),o.Swm_L(1,k), ...
                        o.Sll_L(1,k),o.Sr_L(1,k), ...
                        o.Rexmy_L(1,k,:),p.Stoich_L(k),p.Mpar_L(k),p.fab_L(k),p.fbe_L(k), ...
                        o.RB_L(1,k,:),o.Nreserve_L(1,k),o.Preserve_L(1,k),o.Kreserve_L(1,k), ...
                        o.Careserve_L(1,k),o.Mgreserve_L(1,k),o.Sireserve_L(1,k), ...
                        o.rNc_L(1,k),o.rPc_L(1,k),o.rKc_L(1,k), ...
                        o.rCac_L(1,k),o.rMgc_L(1,k),o.rSic_L(1,k), ...
                        ManIL_c(k),OPT_SoilBiogeochemistry,id_location);

                    o.SAI_L(1,k) = D1.SAI_L(1,k,c);
                    if p.aSE_L(k) == 2
                        o.hc_L(1,k) = GrassHeight(o.LAI_L(1,k),o.LAIdead_L(1,k));
                    elseif p.aSE_L(k) == 5
                        [o.hc_L(1,k),o.SAI_L(1,k),o.B_L(1,:,:),Ccrown,ZR95_L, ...
                         o.Nreserve_L(1,:),o.Preserve_L(1,:),o.Kreserve_L(1,:), ...
                         o.Careserve_L(1,:),o.Mgreserve_L(1,:),o.Sireserve_L(1,:), ...
                         AgrHarNut_cc,RfL_Zs] = CropHeightType( ...
                            o.LAI_L(1,k),o.LAIdead_L(1,k),k,ZR95_L,o.B_L(1,:,:), ...
                            p.Zs,p.CASE_ROOT,Ccrown,Cbare+p.Cwat+p.Curb+p.Crock, ...
                            o.Nreserve_L(1,:),o.Preserve_L(1,:),o.Kreserve_L(1,:), ...
                            o.Careserve_L(1,:),o.Mgreserve_L(1,:),o.Sireserve_L(1,:), ...
                            ManIL_c,p.Mpar_L,p.VegL_Param_Dyn,OPT_SoilBiogeochemistry);
                        o.AgrHarNut(1,:) = Plant_Harvest_Reserve_Export( ...
                            o.AgrHarNut(1,:),AgrHarNut_cc,Ccrown(k));
                        o.SeedInput(1,:) = Plant_Seed_Input_Accumulator( ...
                            o.SeedInput(1,:),SeedInput_cc,Ccrown(k));
                    else
                        o.hc_L(1,k) = D1.hc_L(1,k,c);
                    end
                else
                    o.LAI_L(1,k)=0; o.B_L(1,k,:)=0; o.NPP_L(1,k)=0; o.ANPP_L(1,k)=0;
                    o.Rg_L(1,k)=0; o.RA_L(1,k)=0; o.Rms_L(1,k)=0; o.Rmr_L(1,k)=0;
                    o.PHE_S_L(1,k)=0; o.dflo_L(1,k)=0; o.AgeL_L(1,k)=0;
                    o.e_rel_L(1,k)=0; o.e_relN_L(1,k)=0; o.SAI_L(1,k)=0; o.hc_L(1,k)=0;
                    o.LAIdead_L(1,k)=D1.LAIdead_L(1,k,c);
                    o.Sr_L(1,k)=0; o.Slf_L(1,k)=0; o.Sfr_L(1,k)=0; o.Swm_L(1,k)=0;
                    o.Sll_L(1,k)=0; o.Rrootl_L(1,k)=0; o.Rexmy_L(1,k,:)=0;
                    o.AgeDL_L(1,k)=0; o.FNC_L(1,k)=1;
                    o.Nreserve_L(1,k)=0; o.Preserve_L(1,k)=0; o.Kreserve_L(1,k)=0;
                    o.Careserve_L(1,k)=0; o.Mgreserve_L(1,k)=0; o.Sireserve_L(1,k)=0;
                    o.rNc_L(1,k)=1; o.rPc_L(1,k)=1; o.rKc_L(1,k)=1;
                    o.rCac_L(1,k)=1; o.rMgc_L(1,k)=1; o.rSic_L(1,k)=1;
                    o.Bfac_dayL(1,k)=0; o.Bfac_weekL(1,k)=0; o.NPPI_L(1,k)=0;
                    o.TdpI_L(1,k)=0; o.NupI_L(1,k,:)=0; o.RB_L(1,k,:)=0;
                    o.TexC_L(1,k)=0; o.TexN_L(1,k)=0; o.TexP_L(1,k)=0; o.TexK_L(1,k)=0;
                    o.TexCa_L(1,k)=0; o.TexMg_L(1,k)=0; o.TexSi_L(1,k)=0;
                    o.TNIT_L(1,k)=0; o.TPHO_L(1,k)=0; o.TPOT_L(1,k)=0;
                    o.TCAL_L(1,k)=0; o.TMAG_L(1,k)=0; o.TSIL_L(1,k)=0;
                    o.ISOIL_L(1,k,:)=0; o.NuLit_L(1,k,:)=0;
                end
            end

            %% Crops may resize Ccrown / rooting depth: persist per column
            Ccrown_c(c,:) = Ccrown;
            ZR95_H_c(c,:) = ZR95_H;  ZR95_L_c(c,:) = ZR95_L;
            RfH_Zs_c{c} = RfH_Zs;    RfL_Zs_c{c} = RfL_Zs;
            if OPT_VCA == 0
                o.Ccrown_t(1,:) = Ccrown;
                Cbare_c(c) = 1 - sum(Ccrown) - p.Cwat - p.Curb - p.Crock;
                Check_Land_Cover_Fractions(p.Crock,p.Curb,p.Cwat,Cbare_c(c),Ccrown);
            end
            ManIH(:,c) = ManIH_c; ManIL(:,c) = ManIL_c;
            dout{c} = o;
            %% This square's carbon-reserve events, replayed on the client
            %% in square order (Plant_Carbon_Deficit).
            carbon_n{c} = Plant_Carbon_Deficit('take');
        end
        clear D1 W1
        dout = [dout{:}];
        LAI_L(j,:,:) = cat(3,dout.LAI_L);
        B_L(j,:,:,:) = cat(4,dout.B_L);
        NPP_L(j,:,:) = cat(3,dout.NPP_L);
        Rg_L(j,:,:) = cat(3,dout.Rg_L);
        RA_L(j,:,:) = cat(3,dout.RA_L);
        Rmc_L(j,:,:) = cat(3,dout.Rmc_L);
        ANPP_L(j,:,:) = cat(3,dout.ANPP_L);
        Rms_L(j,:,:) = cat(3,dout.Rms_L);
        Rmr_L(j,:,:) = cat(3,dout.Rmr_L);
        PHE_S_L(j,:,:) = cat(3,dout.PHE_S_L);
        dflo_L(j,:,:) = cat(3,dout.dflo_L);
        AgeL_L(j,:,:) = cat(3,dout.AgeL_L);
        e_rel_L(j,:,:) = cat(3,dout.e_rel_L);
        SAI_L(j,:,:) = cat(3,dout.SAI_L);
        hc_L(j,:,:) = cat(3,dout.hc_L);
        e_relN_L(j,:,:) = cat(3,dout.e_relN_L);
        BA_L(j,:,:) = cat(3,dout.BA_L);
        LAI_H(j,:,:) = cat(3,dout.LAI_H);
        B_H(j,:,:,:) = cat(4,dout.B_H);
        NPP_H(j,:,:) = cat(3,dout.NPP_H);
        Rg_H(j,:,:) = cat(3,dout.Rg_H);
        RA_H(j,:,:) = cat(3,dout.RA_H);
        Rms_H(j,:,:) = cat(3,dout.Rms_H);
        Rmr_H(j,:,:) = cat(3,dout.Rmr_H);
        ANPP_H(j,:,:) = cat(3,dout.ANPP_H);
        PHE_S_H(j,:,:) = cat(3,dout.PHE_S_H);
        dflo_H(j,:,:) = cat(3,dout.dflo_H);
        Rmc_H(j,:,:) = cat(3,dout.Rmc_H);
        AgeL_H(j,:,:) = cat(3,dout.AgeL_H);
        e_rel_H(j,:,:) = cat(3,dout.e_rel_H);
        SAI_H(j,:,:) = cat(3,dout.SAI_H);
        hc_H(j,:,:) = cat(3,dout.hc_H);
        e_relN_H(j,:,:) = cat(3,dout.e_relN_H);
        BA_H(j,:,:) = cat(3,dout.BA_H);
        Rrootl_H(j,:,:) = cat(3,dout.Rrootl_H);
        Rrootl_L(j,:,:) = cat(3,dout.Rrootl_L);
        Bfac_dayH(j,:,:) = cat(3,dout.Bfac_dayH);
        Bfac_weekH(j,:,:) = cat(3,dout.Bfac_weekH);
        NPPI_H(j,:,:) = cat(3,dout.NPPI_H);
        TdpI_H(j,:,:) = cat(3,dout.TdpI_H);
        Bfac_dayL(j,:,:) = cat(3,dout.Bfac_dayL);
        Bfac_weekL(j,:,:) = cat(3,dout.Bfac_weekL);
        NPPI_L(j,:,:) = cat(3,dout.NPPI_L);
        TdpI_L(j,:,:) = cat(3,dout.TdpI_L);
        Sr_H(j,:,:) = cat(3,dout.Sr_H);
        Slf_H(j,:,:) = cat(3,dout.Slf_H);
        Sfr_H(j,:,:) = cat(3,dout.Sfr_H);
        Swm_H(j,:,:) = cat(3,dout.Swm_H);
        Sll_H(j,:,:) = cat(3,dout.Sll_H);
        Sr_L(j,:,:) = cat(3,dout.Sr_L);
        Slf_L(j,:,:) = cat(3,dout.Slf_L);
        Sfr_L(j,:,:) = cat(3,dout.Sfr_L);
        Swm_L(j,:,:) = cat(3,dout.Swm_L);
        Sll_L(j,:,:) = cat(3,dout.Sll_L);
        LAIdead_H(j,:,:) = cat(3,dout.LAIdead_H);
        LAIdead_L(j,:,:) = cat(3,dout.LAIdead_L);
        Rexmy_H(j,:,:,:) = cat(4,dout.Rexmy_H);
        AgeDL_H(j,:,:) = cat(3,dout.AgeDL_H);
        Nreserve_H(j,:,:) = cat(3,dout.Nreserve_H);
        Preserve_H(j,:,:) = cat(3,dout.Preserve_H);
        Kreserve_H(j,:,:) = cat(3,dout.Kreserve_H);
        Careserve_H(j,:,:) = cat(3,dout.Careserve_H);
        Mgreserve_H(j,:,:) = cat(3,dout.Mgreserve_H);
        Sireserve_H(j,:,:) = cat(3,dout.Sireserve_H);
        rNc_H(j,:,:) = cat(3,dout.rNc_H);
        rPc_H(j,:,:) = cat(3,dout.rPc_H);
        rKc_H(j,:,:) = cat(3,dout.rKc_H);
        rCac_H(j,:,:) = cat(3,dout.rCac_H);
        rMgc_H(j,:,:) = cat(3,dout.rMgc_H);
        rSic_H(j,:,:) = cat(3,dout.rSic_H);
        Rexmy_L(j,:,:,:) = cat(4,dout.Rexmy_L);
        AgeDL_L(j,:,:) = cat(3,dout.AgeDL_L);
        Nreserve_L(j,:,:) = cat(3,dout.Nreserve_L);
        Preserve_L(j,:,:) = cat(3,dout.Preserve_L);
        Kreserve_L(j,:,:) = cat(3,dout.Kreserve_L);
        Careserve_L(j,:,:) = cat(3,dout.Careserve_L);
        Mgreserve_L(j,:,:) = cat(3,dout.Mgreserve_L);
        Sireserve_L(j,:,:) = cat(3,dout.Sireserve_L);
        rNc_L(j,:,:) = cat(3,dout.rNc_L);
        rPc_L(j,:,:) = cat(3,dout.rPc_L);
        rKc_L(j,:,:) = cat(3,dout.rKc_L);
        rCac_L(j,:,:) = cat(3,dout.rCac_L);
        rMgc_L(j,:,:) = cat(3,dout.rMgc_L);
        rSic_L(j,:,:) = cat(3,dout.rSic_L);
        NBLeaf_H(j,:,:) = cat(3,dout.NBLeaf_H);
        PARI_H(j,:,:,:) = cat(4,dout.PARI_H);
        NBLI_H(j,:,:) = cat(3,dout.NBLI_H);
        NBLeaf_L(j,:,:) = cat(3,dout.NBLeaf_L);
        PARI_L(j,:,:,:) = cat(4,dout.PARI_L);
        NBLI_L(j,:,:) = cat(3,dout.NBLI_L);
        NupI_H(j,:,:,:) = cat(4,dout.NupI_H);
        NupI_L(j,:,:,:) = cat(4,dout.NupI_L);
        NavlI(j,:,:) = cat(3,dout.NavlI);
        RB_L(j,:,:,:) = cat(4,dout.RB_L);
        RB_H(j,:,:,:) = cat(4,dout.RB_H);
        NuLit_H(j,:,:,:) = cat(4,dout.NuLit_H);
        NuLit_L(j,:,:,:) = cat(4,dout.NuLit_L);
        BLit(j,:,:) = cat(3,dout.BLit);
        Nuptake_H(j,:,:) = cat(3,dout.Nuptake_H);
        NH4uptake_H(j,:,:) = cat(3,dout.NH4uptake_H);
        NO3uptake_H(j,:,:) = cat(3,dout.NO3uptake_H);
        Puptake_H(j,:,:) = cat(3,dout.Puptake_H);
        Kuptake_H(j,:,:) = cat(3,dout.Kuptake_H);
        Cauptake_H(j,:,:) = cat(3,dout.Cauptake_H);
        Mguptake_H(j,:,:) = cat(3,dout.Mguptake_H);
        Siuptake_H(j,:,:) = cat(3,dout.Siuptake_H);
        FNC_H(j,:,:) = cat(3,dout.FNC_H);
        TexC_H(j,:,:) = cat(3,dout.TexC_H);
        TexN_H(j,:,:) = cat(3,dout.TexN_H);
        TexP_H(j,:,:) = cat(3,dout.TexP_H);
        TexK_H(j,:,:) = cat(3,dout.TexK_H);
        TNIT_H(j,:,:) = cat(3,dout.TNIT_H);
        TPHO_H(j,:,:) = cat(3,dout.TPHO_H);
        TPOT_H(j,:,:) = cat(3,dout.TPOT_H);
        SupN_H(j,:,:) = cat(3,dout.SupN_H);
        SupP_H(j,:,:) = cat(3,dout.SupP_H);
        SupK_H(j,:,:) = cat(3,dout.SupK_H);
        TexCa_H(j,:,:) = cat(3,dout.TexCa_H);
        TexMg_H(j,:,:) = cat(3,dout.TexMg_H);
        TexSi_H(j,:,:) = cat(3,dout.TexSi_H);
        TCAL_H(j,:,:) = cat(3,dout.TCAL_H);
        TMAG_H(j,:,:) = cat(3,dout.TMAG_H);
        TSIL_H(j,:,:) = cat(3,dout.TSIL_H);
        SupCa_H(j,:,:) = cat(3,dout.SupCa_H);
        SupMg_H(j,:,:) = cat(3,dout.SupMg_H);
        SupSi_H(j,:,:) = cat(3,dout.SupSi_H);
        Nuptake_L(j,:,:) = cat(3,dout.Nuptake_L);
        NH4uptake_L(j,:,:) = cat(3,dout.NH4uptake_L);
        NO3uptake_L(j,:,:) = cat(3,dout.NO3uptake_L);
        Puptake_L(j,:,:) = cat(3,dout.Puptake_L);
        Kuptake_L(j,:,:) = cat(3,dout.Kuptake_L);
        Cauptake_L(j,:,:) = cat(3,dout.Cauptake_L);
        Mguptake_L(j,:,:) = cat(3,dout.Mguptake_L);
        Siuptake_L(j,:,:) = cat(3,dout.Siuptake_L);
        FNC_L(j,:,:) = cat(3,dout.FNC_L);
        TexC_L(j,:,:) = cat(3,dout.TexC_L);
        TexN_L(j,:,:) = cat(3,dout.TexN_L);
        TexP_L(j,:,:) = cat(3,dout.TexP_L);
        TexK_L(j,:,:) = cat(3,dout.TexK_L);
        TNIT_L(j,:,:) = cat(3,dout.TNIT_L);
        TPHO_L(j,:,:) = cat(3,dout.TPHO_L);
        TPOT_L(j,:,:) = cat(3,dout.TPOT_L);
        SupN_L(j,:,:) = cat(3,dout.SupN_L);
        SupP_L(j,:,:) = cat(3,dout.SupP_L);
        SupK_L(j,:,:) = cat(3,dout.SupK_L);
        TexCa_L(j,:,:) = cat(3,dout.TexCa_L);
        TexMg_L(j,:,:) = cat(3,dout.TexMg_L);
        TexSi_L(j,:,:) = cat(3,dout.TexSi_L);
        TCAL_L(j,:,:) = cat(3,dout.TCAL_L);
        TMAG_L(j,:,:) = cat(3,dout.TMAG_L);
        TSIL_L(j,:,:) = cat(3,dout.TSIL_L);
        SupCa_L(j,:,:) = cat(3,dout.SupCa_L);
        SupMg_L(j,:,:) = cat(3,dout.SupMg_L);
        SupSi_L(j,:,:) = cat(3,dout.SupSi_L);
        ISOIL_H(j,:,:,:) = cat(4,dout.ISOIL_H);
        ISOIL_L(j,:,:,:) = cat(4,dout.ISOIL_L);
        P(j,:,:) = cat(3,dout.P);
        R_litter(j,:,:) = cat(3,dout.R_litter);
        R_microbe(j,:,:) = cat(3,dout.R_microbe);
        R_litter_sur(j,:,:) = cat(3,dout.R_litter_sur);
        R_ew(j,:,:) = cat(3,dout.R_ew);
        VOL(j,:,:) = cat(3,dout.VOL);
        N2flx(j,:,:) = cat(3,dout.N2flx);
        N2Oflx(j,:,:) = cat(3,dout.N2Oflx);
        NITflx(j,:,:) = cat(3,dout.NITflx);
        DENITflx(j,:,:) = cat(3,dout.DENITflx);
        N2Oflx_nit(j,:,:) = cat(3,dout.N2Oflx_nit);
        FertN_d(j,:,:) = cat(3,dout.FertN_d);
        DepN_d(j,:,:) = cat(3,dout.DepN_d);
        PH_d(j,:,:) = cat(3,dout.PH_d);
        PH_pore_d(j,:,:) = cat(3,dout.PH_pore_d);
        PH_bulk_d(j,:,:) = cat(3,dout.PH_bulk_d);
        Min_N(j,:,:) = cat(3,dout.Min_N);
        Min_P(j,:,:) = cat(3,dout.Min_P);
        R_bacteria(j,:,:) = cat(3,dout.R_bacteria);
        RmycAM(j,:,:) = cat(3,dout.RmycAM);
        RmycEM(j,:,:) = cat(3,dout.RmycEM);
        Prod_B(j,:,:) = cat(3,dout.Prod_B);
        Prod_F(j,:,:) = cat(3,dout.Prod_F);
        BfixN(j,:,:) = cat(3,dout.BfixN);
        LitFirEmi(j,:,:) = cat(3,dout.LitFirEmi);
        LEAK_NH4(j,:,:) = cat(3,dout.LEAK_NH4);
        LEAK_NO3(j,:,:) = cat(3,dout.LEAK_NO3);
        LEAK_P(j,:,:) = cat(3,dout.LEAK_P);
        LEAK_K(j,:,:) = cat(3,dout.LEAK_K);
        LEAK_DOC(j,:,:) = cat(3,dout.LEAK_DOC);
        LEAK_DON(j,:,:) = cat(3,dout.LEAK_DON);
        LEAK_DOP(j,:,:) = cat(3,dout.LEAK_DOP);
        LEAK_Ca(j,:,:) = cat(3,dout.LEAK_Ca);
        LEAK_Mg(j,:,:) = cat(3,dout.LEAK_Mg);
        LEAK_Si(j,:,:) = cat(3,dout.LEAK_Si);
        RexmyI(j,:,:) = cat(3,dout.RexmyI);
        AgrHarNut(j,:,:) = cat(3,dout.AgrHarNut);
        SeedInput(j,:,:) = cat(3,dout.SeedInput);
        Ccrown_t(j,:,:) = cat(3,dout.Ccrown_t);
        AgePl_H(j,:,:) = cat(3,dout.AgePl_H);
        AgePl_L(j,:,:) = cat(3,dout.AgePl_L);
        Tden_H(j,:,:) = cat(3,dout.Tden_H);
        Tden_L(j,:,:) = cat(3,dout.Tden_L);
        TBio_Ht(j,:,:) = cat(3,dout.TBio_Ht);
        TBio_Lt(j,:,:) = cat(3,dout.TBio_Lt);
        clear dout
        for c = 1:N_col
            Plant_Carbon_Deficit('replay',carbon_n{c});
        end
        clear carbon_n
    end

    %% ------------------------------------------------------------------
    %% HOURLY BLOCK  --  vertical physics, one square per worker
    %% ------------------------------------------------------------------
    %% Radiation on every cell at once, once an hour, before the squares:
    %% with topography each cell sees the beam through its own slope and
    %% aspect and the sky through its own sky-view factor (v1.0's one
    %% definition: slope/aspect beam projection, sky-view diffuse and a
    %% terrain-reflected term); without it every cell keeps the station
    %% series and the parameter-file SvF, the identical doubles the call
    %% took before, with no per-step allocation.
    if topography_attached
        [sab1_t,sab2_t,sad1_t,sad2_t,parb_t,pard_t] = ...
            Topographic_Radiation_Projection( ...
            SAB1(i),SAB2(i),SAD1(i),SAD2(i),PARB(i),PARD(i), ...
            Datam(i,:),DeltaGMT,Lon,Lat,topo_t_bef,topo_t_aft, ...
            topo.cell_slope,topo.cell_aspect_deg,topo.cell_svf);
    else
        sab1_t = no_topo_rad; sab2_t = no_topo_rad; sad1_t = no_topo_rad;
        sad2_t = no_topo_rad; parb_t = no_topo_rad; pard_t = no_topo_rad;
    end
    %% Each square reads hour i-1 (row sp) and the 24 hours before i, and
    %% writes hour i into the ...n arrays, copied into row si after the loop.
    %% The hour's forcing goes to the workers as the values themselves.
    Datam_i = Datam(i,:); Ta_pd = Ta(pdind);
    Pr_i = Pr(i); Ta_i = Ta(i); Ds_i = Ds(i); Ws_i = Ws(i); IrD_i = IrD(i);
    ea_i = ea(i); N_i = N(i); Pre_i = Pre(i); Ca_i = Ca(i); Tdew_i = Tdew(i);
    SAB1_i = SAB1(i); SAB2_i = SAB2(i); SAD1_i = SAD1(i); SAD2_i = SAD2(i);
    PARB_i = PARB(i); PARD_i = PARD(i);
    Vn = zeros(1,ms,N_col);
    On = zeros(1,ms,N_col);
    Qi_outn = zeros(1,ms,N_col);
    WTRn = zeros(1,ms,N_col);
    POTn = zeros(1,ms,N_col);
    Tdpn = zeros(1,ms,N_col);
    Vicen = zeros(1,ms,N_col);
    Oicen = zeros(1,ms,N_col);
    OFn = zeros(1,1,N_col);
    OSn = zeros(1,1,N_col);
    ZWTn = zeros(1,1,N_col);
    Pr_snon = zeros(1,1,N_col);
    Pr_liqn = zeros(1,1,N_col);
    ern = zeros(1,1,N_col);
    Tsn = zeros(1,1,N_col);
    ran = zeros(1,1,N_col);
    r_soiln = zeros(1,1,N_col);
    b_soiln = zeros(1,1,N_col);
    Ticen = zeros(1,1,N_col);
    Rnn = zeros(1,1,N_col);
    Hn = zeros(1,1,N_col);
    QEn = zeros(1,1,N_col);
    Qvn = zeros(1,1,N_col);
    Lphon = zeros(1,1,N_col);
    EGn = zeros(1,1,N_col);
    ESNn = zeros(1,1,N_col);
    ESN_Inn = zeros(1,1,N_col);
    EWATn = zeros(1,1,N_col);
    EIn_urbn = zeros(1,1,N_col);
    EIn_rockn = zeros(1,1,N_col);
    dw_SNOn = zeros(1,1,N_col);
    Gn = zeros(1,1,N_col);
    SWEn = zeros(1,1,N_col);
    SNDn = zeros(1,1,N_col);
    rosn = zeros(1,1,N_col);
    In_SWEn = zeros(1,1,N_col);
    SP_wcn = zeros(1,1,N_col);
    WR_SPn = zeros(1,1,N_col);
    U_SWEn = zeros(1,1,N_col);
    NIn_SWEn = zeros(1,1,N_col);
    dQn = zeros(1,1,N_col);
    Qfmn = zeros(1,1,N_col);
    t_slsn = zeros(1,1,N_col);
    DQn = zeros(1,1,N_col);
    DTn = zeros(1,1,N_col);
    In_urbn = zeros(1,1,N_col);
    In_rockn = zeros(1,1,N_col);
    SE_rockn = zeros(1,1,N_col);
    SE_urbn = zeros(1,1,N_col);
    Csnon = zeros(1,1,N_col);
    WISn = zeros(1,1,N_col);
    Lkn = zeros(1,1,N_col);
    fn = zeros(1,1,N_col);
    Rhn = zeros(1,1,N_col);
    Rdn = zeros(1,1,N_col);
    NDVIn = zeros(1,1,N_col);
    alp_soiln = zeros(1,1,N_col);
    tau_snon = zeros(1,1,N_col);
    e_snon = zeros(1,1,N_col);
    ALBn = zeros(1,1,N_col);
    EKn = zeros(1,1,N_col);
    dQVEGn = zeros(1,1,N_col);
    TsVn = zeros(1,1,N_col);
    HVn = zeros(1,1,N_col);
    QEVn = zeros(1,1,N_col);
    Cicen = zeros(1,1,N_col);
    Lk_watn = zeros(1,1,N_col);
    Lk_rockn = zeros(1,1,N_col);
    EICEn = zeros(1,1,N_col);
    WATn = zeros(1,1,N_col);
    ICEn = zeros(1,1,N_col);
    ICE_Dn = zeros(1,1,N_col);
    WR_IPn = zeros(1,1,N_col);
    NIcen = zeros(1,1,N_col);
    Cicewn = zeros(1,1,N_col);
    IP_wcn = zeros(1,1,N_col);
    Csnown = zeros(1,1,N_col);
    FROCKn = zeros(1,1,N_col);
    Imeltn = zeros(1,1,N_col);
    Smeltn = zeros(1,1,N_col);
    Tdampn = zeros(1,1,N_col);
    Gfinn = zeros(1,1,N_col);
    In_Littern = zeros(1,1,N_col);
    ELittern = zeros(1,1,N_col);
    Ws_undern = zeros(1,1,N_col);
    r_littern = zeros(1,cc,N_col);
    In_Hn = zeros(1,cc,N_col);
    In_Ln = zeros(1,cc,N_col);
    OHn = zeros(1,cc,N_col);
    OLn = zeros(1,cc,N_col);
    Psi_s_Hn = zeros(1,cc,N_col);
    Psi_s_Ln = zeros(1,cc,N_col);
    Tdp_Hn = zeros(1,cc,N_col);
    Tdp_Ln = zeros(1,cc,N_col);
    rap_Hn = zeros(1,cc,N_col);
    rap_Ln = zeros(1,cc,N_col);
    Dr_Hn = zeros(1,cc,N_col);
    Dr_Ln = zeros(1,cc,N_col);
    T_Hn = zeros(1,cc,N_col);
    T_Ln = zeros(1,cc,N_col);
    EIn_Hn = zeros(1,cc,N_col);
    EIn_Ln = zeros(1,cc,N_col);
    rb_Hn = zeros(1,cc,N_col);
    rb_Ln = zeros(1,cc,N_col);
    An_Hn = zeros(1,cc,N_col);
    Rdark_Hn = zeros(1,cc,N_col);
    An_Ln = zeros(1,cc,N_col);
    Rdark_Ln = zeros(1,cc,N_col);
    Ci_sunHn = zeros(1,cc,N_col);
    Ci_sunLn = zeros(1,cc,N_col);
    Ci_shdHn = zeros(1,cc,N_col);
    Ci_shdLn = zeros(1,cc,N_col);
    rs_sunHn = zeros(1,cc,N_col);
    rs_sunLn = zeros(1,cc,N_col);
    rs_shdHn = zeros(1,cc,N_col);
    rs_shdLn = zeros(1,cc,N_col);
    Vx_Hn = zeros(1,cc,N_col);
    Vl_Hn = zeros(1,cc,N_col);
    Vx_Ln = zeros(1,cc,N_col);
    Vl_Ln = zeros(1,cc,N_col);
    Psi_x_Hn = zeros(1,cc,N_col);
    Psi_l_Hn = zeros(1,cc,N_col);
    Psi_x_Ln = zeros(1,cc,N_col);
    Psi_l_Ln = zeros(1,cc,N_col);
    gsr_Hn = zeros(1,cc,N_col);
    Jsx_Hn = zeros(1,cc,N_col);
    Jxl_Hn = zeros(1,cc,N_col);
    Kleaf_Hn = zeros(1,cc,N_col);
    Kx_Hn = zeros(1,cc,N_col);
    gsr_Ln = zeros(1,cc,N_col);
    Jsx_Ln = zeros(1,cc,N_col);
    Jxl_Ln = zeros(1,cc,N_col);
    Kleaf_Ln = zeros(1,cc,N_col);
    Kx_Ln = zeros(1,cc,N_col);
    fapar_Hn = zeros(1,cc,N_col);
    fapar_Ln = zeros(1,cc,N_col);
    SIF_Hn = zeros(1,cc,N_col);
    SIF_Ln = zeros(1,cc,N_col);
    Tdebn = zeros(1,max(1,length(PAR{1}.Zs_deb)-1),N_col);
    q_runonn = zeros(1,1,N_col);
    Qi_inn = zeros(1,ms,N_col);
    parfor (c = 1:N_col,parallel_workers)
        p = PAR{c};
        Ccrown = Ccrown_c(c,:);
        if N_col > 1
            Qi_in_c = zeros(1,ms);
            Slo_top_cell = 0;
            Slo_pot_cell = zeros(1,ms);
        else
            Qi_in_c = Qi_out(sp,:,c);
            Slo_top_cell = p.Slo_top;
            Slo_pot_cell = p.Slo_pot;
        end
        %% Run-on is what the router delivered last step; nothing without
        %% topography.
        if topography_attached
            q_runon_c = q_runon_next(c);
            sab1_c = sab1_t(c); sab2_c = sab2_t(c); sad1_c = sad1_t(c);
            sad2_c = sad2_t(c); parb_c = parb_t(c); pard_c = pard_t(c);
            svf_c  = cell_svf(c);
        else
            q_runon_c = 0;
            sab1_c = SAB1_i; sab2_c = SAB2_i; sad1_c = SAD1_i;
            sad2_c = SAD2_i; parb_c = PARB_i; pard_c = PARD_i;
            svf_c  = p.SvF;
        end

        [Vn(1,:,c),Vicen(1,:,c),On(1,:,c),Oicen(1,:,c),ZWTn(1,1,c),OFn(1,1,c),OSn(1,1,c), ...
         OHn(1,:,c),OLn(1,:,c),Psi_s_Hn(1,:,c),Psi_s_Ln(1,:,c),Rdn(1,1,c),Qi_outn(1,:,c),WTRn(1,:,c), ...
         Rhn(1,1,c),Lkn(1,1,c),fn(1,1,c),WISn(1,1,c),Tsn(1,1,c),Pr_snon(1,1,c),Pr_liqn(1,1,c),Csnon(1,1,c), ...
         Cicen(1,1,c),NDVIn(1,1,c),rb_Hn(1,:,c),rb_Ln(1,:,c),rs_sunHn(1,:,c), ...
         rs_sunLn(1,:,c),rs_shdHn(1,:,c),rs_shdLn(1,:,c),r_littern(1,:,c), ...
         An_Ln(1,:,c),An_Hn(1,:,c),Rdark_Ln(1,:,c),Rdark_Hn(1,:,c), ...
         Ci_sunHn(1,:,c),Ci_sunLn(1,:,c),Ci_shdHn(1,:,c),Ci_shdLn(1,:,c), ...
         rap_Hn(1,:,c),rap_Ln(1,:,c),r_soiln(1,1,c),b_soiln(1,1,c),alp_soiln(1,1,c),ran(1,1,c),Rnn(1,1,c), ...
         Hn(1,1,c),QEn(1,1,c),Qvn(1,1,c),Lphon(1,1,c),T_Hn(1,:,c),T_Ln(1,:,c),EIn_Hn(1,:,c),EIn_Ln(1,:,c), ...
         EGn(1,1,c),ESNn(1,1,c),ESN_Inn(1,1,c),ELittern(1,1,c),EWATn(1,1,c),EICEn(1,1,c),EIn_urbn(1,1,c), ...
         EIn_rockn(1,1,c),dw_SNOn(1,1,c), ...
         Gn(1,1,c),Gfinn(1,1,c),Tdpn(1,:,c),Tdebn(1,:,c),Tdampn(1,1,c),Ticen(1,1,c), ...
         Tdp_Hn(1,:,c),Tdp_Ln(1,:,c),SWEn(1,1,c),SNDn(1,1,c),rosn(1,1,c),In_SWEn(1,1,c),SP_wcn(1,1,c), ...
         WR_SPn(1,1,c),U_SWEn(1,1,c),NIn_SWEn(1,1,c),dQn(1,1,c),Qfmn(1,1,c),t_slsn(1,1,c),DQn(1,1,c),DTn(1,1,c), ...
         WATn(1,1,c),ICEn(1,1,c),ICE_Dn(1,1,c),IP_wcn(1,1,c),WR_IPn(1,1,c),NIcen(1,1,c),Cicewn(1,1,c), ...
         Csnown(1,1,c),FROCKn(1,1,c),Imeltn(1,1,c),Smeltn(1,1,c), ...
         In_Hn(1,:,c),In_Ln(1,:,c),In_Littern(1,1,c),In_urbn(1,1,c),In_rockn(1,1,c), ...
         Dr_Hn(1,:,c),Dr_Ln(1,:,c),SE_rockn(1,1,c),SE_urbn(1,1,c),Lk_watn(1,1,c),Lk_rockn(1,1,c),ern(1,1,c), ...
         gsr_Hn(1,:,c),Psi_x_Hn(1,:,c),Psi_l_Hn(1,:,c),Jsx_Hn(1,:,c),Jxl_Hn(1,:,c), ...
         Kleaf_Hn(1,:,c),Kx_Hn(1,:,c),Vx_Hn(1,:,c),Vl_Hn(1,:,c), ...
         gsr_Ln(1,:,c),Psi_x_Ln(1,:,c),Psi_l_Ln(1,:,c),Jsx_Ln(1,:,c),Jxl_Ln(1,:,c), ...
         Kleaf_Ln(1,:,c),Kx_Ln(1,:,c),Vx_Ln(1,:,c),Vl_Ln(1,:,c), ...
         fapar_Hn(1,:,c),fapar_Ln(1,:,c),SIF_Hn(1,:,c),SIF_Ln(1,:,c), ...
         snow_alb_new,tau_snon(1,1,c),e_snon(1,1,c),Ws_undern(1,1,c),dQVEGn(1,1,c),HVn(1,1,c), ...
         QEVn(1,1,c),TsVn(1,1,c),EKn(1,1,c),POTn(1,:,c)] = ...
         HYDROLOGIC_UNIT(V(sp,:,c),Oice(sp,:,c),p.aR,p.Zs, ...
            p.EvL_Zs,p.Inf_Zs,p.Zinf,RfH_Zs_c{c},RfL_Zs_c{c},p.dz,p.Dz,ms,p.Kbot, ...
            Pr_i,Ta_i,Ds_i,Ws_i,p.zatm,Ts(sp,1,c),IrD_i,dt,dth,ea_i,N_i,Pre_i,Tstm0_c(c), ...
            LAI_H(j,:,c),SAI_H(j,:,c),LAI_L(j,:,c),SAI_L(j,:,c),LAIdead_H(j,:,c),LAIdead_L(j,:,c), ...
            Rrootl_H(j,:,c),Rrootl_L(j,:,c),BLit(j,:,c),p.Sllit,p.Kct, ...
            Datam_i,DeltaGMT,Lon,Lat,t_bef,t_aft,h_S, ...
            Ccrown,Cbare_c(c),p.Crock,p.Curb,p.Cwat, ...
            p.Soil_Param,p.Interc_Param,p.SnowIce_Param,p.VegH_Param,p.VegL_Param, ...
            p.Zs_deb,p.Deb_Par, ...
            ZR95_H_c(c,:),ZR95_L_c(c,:), ...
            sab1_c,sab2_c,sad1_c,sad2_c,parb_c,pard_c,svf_c, ...
            SND(sp,1,c),snow_alb_c{c},p.Color_Class,p.OM_H,p.OM_L, ...
            p.PFT_opt_H,p.PFT_opt_L,hc_H(j,:,c),hc_L(j,:,c),p.d_leaf_H,p.d_leaf_L, ...
            Ca_i,Oa,Ci_sunH(sp,:,c),Ci_shdH(sp,:,c),Ci_sunL(sp,:,c),Ci_shdL(sp,:,c), ...
            e_rel_H(j,:,c),e_relN_H(j,:,c),e_rel_L(j,:,c),e_relN_L(j,:,c), ...
            e_sno(sp,1,c),In_H(sp,:,c),In_L(sp,:,c),In_Litter(sp,1,c),In_urb(sp,1,c), ...
            In_rock(sp,1,c),SWE(sp,1,c),In_SWE(sp,1,c), ...
            Tdeb(sp,:,c),Tdp(sp,:,c),Tdamp(sp,1,c),Tice(sp,1,c), ...
            WAT(sp,1,c),ICE(sp,1,c),IP_wc(sp,1,c),ICE_D(sp,1,c),Cicew(sp,1,c), ...
            Vx_H(sp,:,c),Vl_H(sp,:,c),Vx_L(sp,:,c),Vl_L(sp,:,c), ...
            Psi_x_H(sp,:,c),Psi_l_H(sp,:,c),Psi_x_L(sp,:,c),Psi_l_L(sp,:,c), ...
            FROCK(sp,1,c),p.Krock,Ws_under(sp,1,c), ...
            Tdew_i,t_sls(sp,1,c),ros(sp,1,c),SP_wc(sp,1,c),p.fpr,Pr_sno(pd,1,c), ...
            p.In_max_urb,p.In_max_rock,p.K_usle,tau_sno(sp,1,c),Ta_pd, ...
            Slo_top_cell,Slo_pot_cell,p.Asur,p.Ared,p.aTop,EK(sp,1,c),q_runon_c,Qi_in_c, ...
            p.pow_dis,p.a_dis, ...
            p.SPAR,p.SN,OPT_VegSnow,OPT_SoilTemp,OPT_PlantHydr,Opt_CR,Opt_ST, ...
            OPT_SM,OPT_STh,OPT_FR_SOIL,OPT_PH);

        snow_alb_c{c} = snow_alb_new;
        Tstm0_c(c) = 2*Tsn(1,1,c) - Ts(sp,1,c);
        %% Weighted by the radiation this cell received (the station series
        %% without topography, so the same doubles as before).
        STOT = sab1_c+sab2_c+sad1_c+sad2_c;
        ALBn(1,1,c) = sab1_c/STOT*snow_alb_new.dir_vis + sad1_c/STOT*snow_alb_new.dif_vis + ...
                   sab2_c/STOT*snow_alb_new.dir_nir + sad2_c/STOT*snow_alb_new.dif_nir;
        Qi_inn(1,:,c) = Qi_in_c;
        q_runonn(1,1,c) = q_runon_c;
    end
    V(si,:,:) = Vn;
    O(si,:,:) = On;
    Qi_out(si,:,:) = Qi_outn;
    WTR(si,:,:) = WTRn;
    POT(si,:,:) = POTn;
    Tdp(si,:,:) = Tdpn;
    Vice(si,:,:) = Vicen;
    Oice(si,:,:) = Oicen;
    OF(si,:,:) = OFn;
    OS(si,:,:) = OSn;
    ZWT(si,:,:) = ZWTn;
    Pr_sno(si,:,:) = Pr_snon;
    Pr_liq(si,:,:) = Pr_liqn;
    er(si,:,:) = ern;
    Ts(si,:,:) = Tsn;
    ra(si,:,:) = ran;
    r_soil(si,:,:) = r_soiln;
    b_soil(si,:,:) = b_soiln;
    Tice(si,:,:) = Ticen;
    Rn(si,:,:) = Rnn;
    H(si,:,:) = Hn;
    QE(si,:,:) = QEn;
    Qv(si,:,:) = Qvn;
    Lpho(si,:,:) = Lphon;
    EG(si,:,:) = EGn;
    ESN(si,:,:) = ESNn;
    ESN_In(si,:,:) = ESN_Inn;
    EWAT(si,:,:) = EWATn;
    EIn_urb(si,:,:) = EIn_urbn;
    EIn_rock(si,:,:) = EIn_rockn;
    dw_SNO(si,:,:) = dw_SNOn;
    G(si,:,:) = Gn;
    SWE(si,:,:) = SWEn;
    SND(si,:,:) = SNDn;
    ros(si,:,:) = rosn;
    In_SWE(si,:,:) = In_SWEn;
    SP_wc(si,:,:) = SP_wcn;
    WR_SP(si,:,:) = WR_SPn;
    U_SWE(si,:,:) = U_SWEn;
    NIn_SWE(si,:,:) = NIn_SWEn;
    dQ(si,:,:) = dQn;
    Qfm(si,:,:) = Qfmn;
    t_sls(si,:,:) = t_slsn;
    DQ(si,:,:) = DQn;
    DT(si,:,:) = DTn;
    In_urb(si,:,:) = In_urbn;
    In_rock(si,:,:) = In_rockn;
    SE_rock(si,:,:) = SE_rockn;
    SE_urb(si,:,:) = SE_urbn;
    Csno(si,:,:) = Csnon;
    WIS(si,:,:) = WISn;
    Lk(si,:,:) = Lkn;
    f(si,:,:) = fn;
    Rh(si,:,:) = Rhn;
    Rd(si,:,:) = Rdn;
    NDVI(si,:,:) = NDVIn;
    alp_soil(si,:,:) = alp_soiln;
    tau_sno(si,:,:) = tau_snon;
    e_sno(si,:,:) = e_snon;
    ALB(si,:,:) = ALBn;
    EK(si,:,:) = EKn;
    dQVEG(si,:,:) = dQVEGn;
    TsV(si,:,:) = TsVn;
    HV(si,:,:) = HVn;
    QEV(si,:,:) = QEVn;
    Cice(si,:,:) = Cicen;
    Lk_wat(si,:,:) = Lk_watn;
    Lk_rock(si,:,:) = Lk_rockn;
    EICE(si,:,:) = EICEn;
    WAT(si,:,:) = WATn;
    ICE(si,:,:) = ICEn;
    ICE_D(si,:,:) = ICE_Dn;
    WR_IP(si,:,:) = WR_IPn;
    NIce(si,:,:) = NIcen;
    Cicew(si,:,:) = Cicewn;
    IP_wc(si,:,:) = IP_wcn;
    Csnow(si,:,:) = Csnown;
    FROCK(si,:,:) = FROCKn;
    Imelt(si,:,:) = Imeltn;
    Smelt(si,:,:) = Smeltn;
    Tdamp(si,:,:) = Tdampn;
    Gfin(si,:,:) = Gfinn;
    In_Litter(si,:,:) = In_Littern;
    ELitter(si,:,:) = ELittern;
    Ws_under(si,:,:) = Ws_undern;
    r_litter(si,:,:) = r_littern;
    In_H(si,:,:) = In_Hn;
    In_L(si,:,:) = In_Ln;
    OH(si,:,:) = OHn;
    OL(si,:,:) = OLn;
    Psi_s_H(si,:,:) = Psi_s_Hn;
    Psi_s_L(si,:,:) = Psi_s_Ln;
    Tdp_H(si,:,:) = Tdp_Hn;
    Tdp_L(si,:,:) = Tdp_Ln;
    rap_H(si,:,:) = rap_Hn;
    rap_L(si,:,:) = rap_Ln;
    Dr_H(si,:,:) = Dr_Hn;
    Dr_L(si,:,:) = Dr_Ln;
    T_H(si,:,:) = T_Hn;
    T_L(si,:,:) = T_Ln;
    EIn_H(si,:,:) = EIn_Hn;
    EIn_L(si,:,:) = EIn_Ln;
    rb_H(si,:,:) = rb_Hn;
    rb_L(si,:,:) = rb_Ln;
    An_H(si,:,:) = An_Hn;
    Rdark_H(si,:,:) = Rdark_Hn;
    An_L(si,:,:) = An_Ln;
    Rdark_L(si,:,:) = Rdark_Ln;
    Ci_sunH(si,:,:) = Ci_sunHn;
    Ci_sunL(si,:,:) = Ci_sunLn;
    Ci_shdH(si,:,:) = Ci_shdHn;
    Ci_shdL(si,:,:) = Ci_shdLn;
    rs_sunH(si,:,:) = rs_sunHn;
    rs_sunL(si,:,:) = rs_sunLn;
    rs_shdH(si,:,:) = rs_shdHn;
    rs_shdL(si,:,:) = rs_shdLn;
    Vx_H(si,:,:) = Vx_Hn;
    Vl_H(si,:,:) = Vl_Hn;
    Vx_L(si,:,:) = Vx_Ln;
    Vl_L(si,:,:) = Vl_Ln;
    Psi_x_H(si,:,:) = Psi_x_Hn;
    Psi_l_H(si,:,:) = Psi_l_Hn;
    Psi_x_L(si,:,:) = Psi_x_Ln;
    Psi_l_L(si,:,:) = Psi_l_Ln;
    gsr_H(si,:,:) = gsr_Hn;
    Jsx_H(si,:,:) = Jsx_Hn;
    Jxl_H(si,:,:) = Jxl_Hn;
    Kleaf_H(si,:,:) = Kleaf_Hn;
    Kx_H(si,:,:) = Kx_Hn;
    gsr_L(si,:,:) = gsr_Ln;
    Jsx_L(si,:,:) = Jsx_Ln;
    Jxl_L(si,:,:) = Jxl_Ln;
    Kleaf_L(si,:,:) = Kleaf_Ln;
    Kx_L(si,:,:) = Kx_Ln;
    fapar_H(si,:,:) = fapar_Hn;
    fapar_L(si,:,:) = fapar_Ln;
    SIF_H(si,:,:) = SIF_Hn;
    SIF_L(si,:,:) = SIF_Ln;
    Tdeb(si,:,:) = Tdebn;
    q_runon(si,:,:) = q_runonn;
    Qi_in(si,:,:) = Qi_inn;
    clear Vn On Qi_outn WTRn POTn Tdpn Vicen Oicen OFn OSn
    clear ZWTn Pr_snon Pr_liqn ern Tsn ran r_soiln b_soiln Ticen Rnn
    clear Hn QEn Qvn Lphon EGn ESNn ESN_Inn EWATn EIn_urbn EIn_rockn
    clear dw_SNOn Gn SWEn SNDn rosn In_SWEn SP_wcn WR_SPn U_SWEn NIn_SWEn
    clear dQn Qfmn t_slsn DQn DTn In_urbn In_rockn SE_rockn SE_urbn Csnon
    clear WISn Lkn fn Rhn Rdn NDVIn alp_soiln tau_snon e_snon ALBn
    clear EKn dQVEGn TsVn HVn QEVn Cicen Lk_watn Lk_rockn EICEn WATn
    clear ICEn ICE_Dn WR_IPn NIcen Cicewn IP_wcn Csnown FROCKn Imeltn Smeltn
    clear Tdampn Gfinn In_Littern ELittern Ws_undern r_littern In_Hn In_Ln OHn OLn
    clear Psi_s_Hn Psi_s_Ln Tdp_Hn Tdp_Ln rap_Hn rap_Ln Dr_Hn Dr_Ln T_Hn T_Ln
    clear EIn_Hn EIn_Ln rb_Hn rb_Ln An_Hn Rdark_Hn An_Ln Rdark_Ln Ci_sunHn Ci_sunLn
    clear Ci_shdHn Ci_shdLn rs_sunHn rs_sunLn rs_shdHn rs_shdLn Vx_Hn Vl_Hn Vx_Ln Vl_Ln
    clear Psi_x_Hn Psi_l_Hn Psi_x_Ln Psi_l_Ln gsr_Hn Jsx_Hn Jxl_Hn Kleaf_Hn Kx_Hn gsr_Ln
    clear Jsx_Ln Jxl_Ln Kleaf_Ln Kx_Ln fapar_Hn fapar_Ln SIF_Hn SIF_Ln Tdebn q_runonn
    clear Qi_inn
    for c = 1:N_col
        %% Vertical exit at the base of the column [m3]
        exit_water_vert_m3(c) = exit_water_vert_m3(c) + (Lk(si,c)/1000)*A_cell(c)*dth;
        if ~topography_attached
            %% Flat graph: unrouted runoff leaves each cell where it forms.
            exit_water_surface_m3(c) = exit_water_surface_m3(c) + ...
                ((Rh(si,c)+Rd(si,c))/1000)*A_cell(c);
        end
        %% With topography, this step's runoff is routed after the lateral
        %% block and exits only at the outlet cell(s).
    end

    %% ------------------------------------------------------------------
    %% LATERAL COUPLING  --  operator-split, after every column has
    %% taken its vertical step.
    %% ------------------------------------------------------------------
    if N_col > 1
        theta_cell = reshape(O(si,:,:), ms, N_col).';        %% [N_col x ms]
        Oice_cell = reshape(Oice(si,:,:), ms, N_col).';
        T_cell = reshape(Tdp(si,:,:), ms, N_col).';

        %% Biogeo_Leakage defines which T&C pools are dissolved and mobile.
        %% The pools are vertically integrated [gX/m2]; the T&C leakage
        %% formulation treats them as well mixed over the column water.
        P_cells = reshape(P(j,:,:),79,N_col).';
        [~,mobile] = Biogeo_Lateral_Pools(P_cells);
        mobile = mobile .* lateral_fraction;
        C_aq = Biogeo_Lateral_Concentration( ...
            mobile,A_cell,theta_cell,Zlayer_m,bio_active_fraction);

        lambda_cell = zeros(N_col,ms);
        capp_cell = zeros(N_col,ms);
        for c = 1:N_col
            p = PAR{c};
            [lambda_cell(c,:),capp_cell(c,:)] = Soil_Lateral_Thermal_State( ...
                T_cell(c,:),theta_cell(c,:)+Oice_cell(c,:), ...
                Pre(i),p,OPT_FR_SOIL);
        end
        thermal = struct('T_cell',T_cell,'lambda_cell',lambda_cell, ...
            'capp_cell',capp_cell,'Oice_cell',Oice_cell);
        thermal.parameter_cell = PAR;
        thermal.Pre = Pre(i);
        thermal.OPT_FR_SOIL = OPT_FR_SOIL;
        [theta_new, dM_cell, ldiag, T_new, Oice_new] = Lateral_Routing_Step( ...
            topo,theta_cell,C_aq,soilpar_cell,Zlayer_m,dth,thermal, ...
            bio_active_fraction);

        %% Write the redistributed moisture back, keeping V consistent
        %% with O through the column's own residual water content.
        for c = 1:N_col
            p = PAR{c};
            O(si,:,c) = theta_new(c,:);
            Oice(si,:,c) = Oice_new(c,:);
            V(si,:,c) = (theta_new(c,:) - p.Ohy) .* p.dz;
            Vice(si,:,c) = Oice_new(c,:) .* p.dz;
            Tdp(si,:,c) = T_new(c,:);
            [Tdp_H(si,:,c),Tdp_L(si,:,c)] = RootZone_Temp( ...
                T_new(c,:),RfH_Zs_c{c},RfL_Zs_c{c});
        end

        %% Return the advected mass to its T&C pools [gX/m2].
        dM_area = squeeze(sum(dM_cell,2)) ./ A_cell;
        if isfield(lateral_diag,'water_cell_net_mm')
            Q_gain = accumarray(topo.edge_b(:),sum(ldiag.Q_edge,2),[N_col 1]) - ...
                accumarray(topo.edge_a(:),sum(ldiag.Q_edge,2),[N_col 1]);  %% m3/h
            lateral_diag.water_cell_net_mm(j,:) = lateral_diag.water_cell_net_mm(j,:) + ...
                (1000*Q_gain*dth./A_cell).';
            lateral_diag.NO3_cell_net_g_m2(j,:) = lateral_diag.NO3_cell_net_g_m2(j,:) + ...
                dM_area(:,i_lateral_NO3).';
        end
        P_cells = Biogeo_Lateral_Pools(P_cells,dM_area);
        P(j,:,:) = reshape(P_cells.',1,79,N_col);
        transported = unique([lateral_pool_idx{:}]);
        Pj = P(j,transported,:);
        if any(Pj(:) < -1e-10)
            error('MAIN_FRAME_BG_ERW_2D:NegativePool', ...
                'Conservative lateral transport produced a negative pool at i=%d.',i);
        end

        lateral_diag.water_buf_to_main(j) = lateral_diag.water_buf_to_main(j) + ldiag.buf_to_main.water;
        lateral_diag.mass_buf_to_main(j,:) = lateral_diag.mass_buf_to_main(j,:) + ldiag.buf_to_main.net.';
        lateral_diag.NO3_buf_to_main(j)   = lateral_diag.NO3_buf_to_main(j) + ldiag.buf_to_main.net(i_lateral_NO3);
        lateral_diag.sub_steps_max(j)     = max(lateral_diag.sub_steps_max(j), ldiag.sub_steps);
        lateral_diag.heat_buf_to_main_J(j)= lateral_diag.heat_buf_to_main_J(j) + ldiag.buf_to_main.heat_J;
        lateral_diag.heat_residual_J(j)   = lateral_diag.heat_residual_J(j) + ldiag.heat.energy_residual_J;
        lateral_diag.heat_residual_relative(j) = max( ...
            lateral_diag.heat_residual_relative(j), ...
            ldiag.heat.energy_residual_relative);
        lateral_diag.heat_sub_steps_max(j)= max(lateral_diag.heat_sub_steps_max(j),ldiag.heat.sub_steps);
        if OPT_DEM == 1
            lateral_diag.water_edge_abs_m3(j) = lateral_diag.water_edge_abs_m3(j) + ...
                sum(abs(ldiag.Q_edge(:)))*dth;
            Q_face = sum(ldiag.Q_edge,2);                         %% m3/h, a -> b
            Q_cell_net = accumarray(topo.edge_b(:),Q_face,[N_col 1]) - ...
                accumarray(topo.edge_a(:),Q_face,[N_col 1]);      %% m3/h gained
            lateral_diag.water_net_m3(j) = lateral_diag.water_net_m3(j) + ...
                0.5*sum(abs(Q_cell_net))*dth;
        end
    end

    %% ------------------------------------------------------------------
    %% OVERLAND ROUTING (OPT_DEM). Ponded runoff runs down the DEM carrying
    %% its dissolved load, leaves the domain only at the topographic
    %% outlet(s), and is otherwise offered to the next step as run-on.
    %% Ported from v1.0 with its field-origin tracers: on a polygon grid
    %% the perimeter cells are whole squares only partly inside the
    %% polygon, and the field's share of the export is carried, not
    %% inferred from the outlet cell's overlap.
    %%
    %% TWO SOURCES OF SURFACE SOLUTE. Saturation-excess runoff (Rd) is soil
    %% water reaching the surface and carries the active-zone pore-water
    %% concentration, capped at the mobile pool it comes from.
    %% Infiltration-excess runoff (Rh) is ponded water that re-mobilises the
    %% pond store in the share Rh/WIS, WIS being the water Rh is defined
    %% against in HYDROLOGIC_UNIT. What the router delivers goes to the
    %% NEXT step's pond, not into the soil, because its water has not
    %% infiltrated yet.
    %% ------------------------------------------------------------------
    if N_col > 1 && topography_attached
        Rh_c = reshape(Rh(si,1,:),N_col,1);
        Rd_c = reshape(Rd(si,1,:),N_col,1);
        runoff_step_mm = max(Rh_c + Rd_c, 0);
        depth_m = runoff_step_mm/1000;
        P_cells_ov = reshape(P(j,:,:),79,N_col).';
        [~,mobile_now] = Biogeo_Lateral_Pools(P_cells_ov);
        mobile_now = mobile_now .* lateral_fraction;
        theta_ov = reshape(O(si,:,:),ms,N_col).';
        C_soil = Biogeo_Lateral_Concentration( ...
            mobile_now,A_cell,theta_ov,Zlayer_m,bio_active_fraction);
        C_soil = reshape(C_soil(:,1,:),N_col,N_lateral_species);

        WIS_c    = reshape(WIS(si,1,:),N_col,1);
        runon_mm = reshape(q_runon(si,1,:),N_col,1)*dth;
        out_frac = min(1, Rh_c./max(WIS_c,eps));
        out_frac(WIS_c <= 0) = 0;
        pond_out_area = pond_mass_area .* out_frac;
        pond_in_area  = pond_mass_area - pond_out_area;
        pond_in_water_mm = runon_mm.*(1-out_frac);

        soil_out_area = min((Rd_c/1000).*C_soil, mobile_now);
        soil_out_area(~isfinite(soil_out_area) | soil_out_area < 0) = 0;
        C_surf = (soil_out_area + pond_out_area) ./ max(depth_m,eps);
        C_surf(~isfinite(C_surf) | C_surf < 0) = 0;
        removed_area = depth_m .* C_surf;
        soil_share = soil_out_area ./ max(soil_out_area+pond_out_area,eps);
        soil_removed_area = removed_area .* soil_share;

        %% FIELD-ORIGIN provenance as passive tracer columns (v1.0). The Rd
        %% part is this cell's own water, so its field share is the cell's
        %% overlap fraction; the Rh part is surface water whose field share
        %% is the pond's own provenance plus local rain. Solute likewise.
        surf_field_mm = pond_field_water_mm + ...
            max(WIS_c - runon_mm,0).*field_fraction;
        rh_field_frac = surf_field_mm ./ max(WIS_c,eps);
        water_tracer = (Rd_c.*field_fraction + Rh_c.*rh_field_frac) ./ ...
            max(runoff_step_mm,eps);
        pond_field_out_area = pond_field_mass_area .* out_frac;
        mass_tracer = (soil_removed_area.*field_fraction + pond_field_out_area) ./ ...
            max(depth_m,eps);
        water_tracer(~isfinite(water_tracer)) = 0;
        mass_tracer(~isfinite(mass_tracer)) = 0;
        C_route = [C_surf, water_tracer, mass_tracer];
        [q_deliver_mm,C_deliver,~,exit_tracer,~,exit_cell_m3,exit_cell_g] = ...
            Surface_Runoff_Routing(topo,runoff_step_mm,C_route,dth);

        %% Only the pond share that infiltrated joins the soil pools; the
        %% soil-derived load that ran off leaves them. The routed load
        %% stays on the surface with its provenance.
        P_cells_ov = Biogeo_Lateral_Pools(P_cells_ov, pond_in_area - soil_removed_area);
        P(j,:,:) = reshape(P_cells_ov.',1,79,N_col);
        pond_mass_area       = (q_deliver_mm/1000) .* C_deliver(:,1:N_lateral_species);
        pond_field_water_mm  = q_deliver_mm .* C_deliver(:,N_lateral_species+1);
        pond_field_mass_area = (q_deliver_mm/1000) .* ...
            C_deliver(:,N_lateral_species+2:end);

        %% The ERW aqueous stores ride the same water, so they must leave
        %% overland too or the exported solution is not charge balanced.
        if OPT_ERW == 1 && all(cellfun(@(x) isfield(x,'init') && x.init,ERW))
            [ERW,erw_surf] = ERW_Surface_Export(ERW,topo,runoff_step_mm,Rd_c, ...
                out_frac,erw_pond_umol,erw_pond_field_umol,pond_in_water_mm, ...
                field_fraction,dth);
            erw_pond_umol       = erw_surf.delivered_umol;
            erw_pond_field_umol = erw_surf.delivered_field_umol;
            exit_mass_surface_erw_umol = ...
                exit_mass_surface_erw_umol + erw_surf.exit_umol;
            exit_mass_surface_erw_field_umol = ...
                exit_mass_surface_erw_field_umol + erw_surf.exit_field_umol;
        end

        %% q_runon is a RATE [mm/h]; q_deliver_mm is a depth per step.
        q_runon_next          = q_deliver_mm/dth;
        exit_water_surface_m3 = exit_water_surface_m3 + exit_cell_m3;
        exit_mass_surface_g   = exit_mass_surface_g + ...
            sum(exit_cell_g(:,1:N_lateral_species),1).';
        exit_water_field_surface_m3 = exit_water_field_surface_m3 + ...
            exit_tracer(N_lateral_species+1);
        exit_mass_surface_field_g = exit_mass_surface_field_g + ...
            reshape(exit_tracer(N_lateral_species+2:end),[],1);
        lateral_diag.surface_water_out_m3(j) = ...
            lateral_diag.surface_water_out_m3(j) + sum(exit_cell_m3);
        lateral_diag.surface_mass_out_g(j,:) = ...
            lateral_diag.surface_mass_out_g(j,:) + ...
            sum(exit_cell_g(:,1:N_lateral_species),1);
        %% The run-on delivered for the next step's re-infiltration, as in
        %% v1.0: with the outlet export it is every step's runoff exactly
        %% (the router's own closure), so the two ledgers close the
        %% overland budget over a run.
        lateral_diag.surface_reinfilt_m3(j) = ...
            lateral_diag.surface_reinfilt_m3(j) + sum((q_deliver_mm/1000).*A_cell);

        Pj = P(j,transported,:);
        if any(Pj(:) < -1e-9)
            error('MAIN_FRAME_BG_ERW_2D:NegativePool', ...
                'Overland routing produced a negative pool at i=%d.',i);
        end
    end

    if N_col > 1 && OPT_ERW == 1 && ...
            all(cellfun(@(x) isfield(x,'init') && x.init,ERW))
        [ERW,erw_ldiag] = ERW_Lateral_Transport( ...
            ERW,topo,ldiag.Q_edge,theta_cell,Zlayer_m,bio_active_fraction, ...
            dth);
        lateral_diag.erw_mass_buf_to_main_umol(j,:) = ...
            lateral_diag.erw_mass_buf_to_main_umol(j,:) + erw_ldiag.buf_to_main.aqueous_umol.';
        lateral_diag.erw_CO2_buf_to_main_umol(j) = ...
            lateral_diag.erw_CO2_buf_to_main_umol(j) + erw_ldiag.buf_to_main.CO2_umol;
        lateral_diag.erw_mass_residual_umol(j,:) = ...
            lateral_diag.erw_mass_residual_umol(j,:) + erw_ldiag.mass_residual_umol.';
        lateral_diag.erw_mass_abs_residual_umol(j,:) = ...
            lateral_diag.erw_mass_abs_residual_umol(j,:) + abs(erw_ldiag.mass_residual_umol).';
        erw_throughput_umol = sum(abs(erw_ldiag.aqueous_edge_umol),1);
        erw_throughput_umol(4) = erw_throughput_umol(4) + ...
            sum(abs(erw_ldiag.CO2_edge_umol));
        lateral_diag.erw_mass_throughput_umol(j,:) = ...
            lateral_diag.erw_mass_throughput_umol(j,:) + erw_throughput_umol;
    end

    %% ------------------------------------------------------------------
    %% HOURLY FIELD MEANS (the tower-comparison fluxes), hour i complete
    %% ------------------------------------------------------------------
    fh = Field_Hour_Values(si,N_col,cc,n_bio_layers,Ccrown_c, ...
        Rn,H,QE,G,Ts,EG,ELitter,ESN,ESN_In,EWAT,EICE,EIn_urb,EIn_rock, ...
        Lk,Rh,Rd,T_H,T_L,EIn_H,EIn_L,An_H,An_L,Rdark_H,Rdark_L,O,Tdp);
    fhm = Field_Hourly_Means(fh,w_field);
    for k_fh = 1:numel(FIELD_HOURLY_NAMES)
        Field_hourly.(FIELD_HOURLY_NAMES{k_fh})(i,:) = fhm.(FIELD_HOURLY_NAMES{k_fh});
    end

    %% ------------------------------------------------------------------
    %% ERW  --  once per day, after the current lateral transport
    %% ------------------------------------------------------------------
    if OPT_ERW == 1 && ((Datam(i,4) == 1) || i == NN)
        %% Only the rows the day reads go to the workers: the last 50
        %% hours of the hourly records, yesterday's and today's daily rows,
        %% and the CO2 forcing around the day (all of it on the first call,
        %% which seeds each square from the full-length timetable).
        h_lo = max(1,i-HOURLY_WINDOW+1);
        hw = mod((h_lo:i)-1,NH)+1;
        O_w = O(hw,:,:); Lk_w = Lk(hw,:,:); f_w = f(hw,:,:);
        T_H_w = T_H(hw,:,:); T_L_w = T_L(hw,:,:);
        P_w = P(j-1:j,:,:);
        R_microbe_w = R_microbe(j,:,:); R_litter_w = R_litter(j,:,:);
        R_ew_w = R_ew(j,:,:); R_litter_sur_w = R_litter_sur(j,:,:);
        Rmr_H_w = Rmr_H(j,:,:); Rmr_L_w = Rmr_L(j,:,:);
        co2_lo = max(1,i-HOURLY_WINDOW+1); co2_hi = min(NN,i+48);
        CO2_w = struct('first',co2_lo,'ppm',Ca(co2_lo:co2_hi),'full',[]);
        if ~all(cellfun(@(x) isfield(x,'init') && x.init,ERW))
            CO2_w.full = Ca;
        end
        ERW_next = cell(N_col,1);
        parfor (c = 1:N_col,parallel_workers)
            ERW_PARAM_cell = ERW_PARAM_input;
            if isfield(ColCfg,'ERW_PARAM')
                ERW_PARAM_cell = ColCfg(c).ERW_PARAM;
            end
            ERW_next{c} = ERW_Daily_Step_2D(ERW{c}, ColCfg(c), PAR{c}, ...
                i, j, NN, NNd, Datam, Date, CO2_w, ...
                O_w(:,:,c), Lk_w(:,1,c), f_w(:,1,c), T_H_w(:,:,c), T_L_w(:,:,c), ...
                Tdp_bio_c(c), VSUM_c(c), VTSUM_c(c), ...
                R_microbe_w(1,1,c), R_litter_w(1,1,c), R_ew_w(1,1,c), R_litter_sur_w(1,1,c), ...
                Rmr_H_w(1,:,c), Rmr_L_w(1,:,c), Ccrown_c(c,:), ...
                P_w(:,:,c), MOD_PARAM_input, ERW_PARAM_cell, ...
                conv_mol, conv_Al, N_col, h_lo, FULL_RECORDS);
        end
        clear O_w Lk_w f_w T_H_w T_L_w P_w R_microbe_w R_litter_w R_ew_w R_litter_sur_w Rmr_H_w Rmr_L_w CO2_w
        for c = 1:N_col
            ERW{c} = ERW_next{c};
            %% The call before this one, summarised now that its rows hold
            %% their final values (ERW_Daily_Record).
            if isfield(ERW{c},'daily_row') && ~isempty(ERW{c}.daily_row)
                ERW_daily{c} = ERW_Daily_Record('store',ERW_daily{c},ERW{c}.daily_row,NNd);
                ERW{c}.daily_row = [];
            end
            if j == 2 && isfield(ERW{c},'P_initial')
                P(j-1,:,c) = ERW{c}.P_initial;
            end
            %% Cold-start backfill of the pH records (see the seed comment
            %% above the PH_c/PH_bulk_c initialisation). Rows 1:j were
            %% written before this first ERW call from the bulk observation,
            %% so the pore record carried the bulk value. bg row 1 is the
            %% t = 0 ERW initial condition (the first call's block starts at
            %% row 1 of the grid): pH(1) the pore-water pH that
            %% Invert_Bulk_pH recovered from the bulk observation, pH_bulk(1)
            %% the bulk pH Porewater_to_Bulk_pH returns for that state (the
            %% observation itself, to solver tolerance). The daily loop in
            %% ERW_Daily_Step_2D starts at row 2, so row 1 is never
            %% overwritten. Same pairing as SCEPTER: one pore-water H state
            %% solved at t = 0, soil pH derived from it. Records only;
            %% PH_c/PH_bulk_c are updated from the daily step below.
            if j == 2 && isfield(ERW{c},'bg') && isfield(ERW{c}.bg,'pH') && ...
                    isfield(ERW{c}.bg,'pH_bulk') && ...
                    isfinite(ERW{c}.bg.pH(1,1)) && isfinite(ERW{c}.bg.pH_bulk(1,1))
                PH_pore_d(1:j,c) = ERW{c}.bg.pH(1,1);
                PH_bulk_d(1:j,c) = ERW{c}.bg.pH_bulk(1,1);
                PH_d(1:j,c)      = ERW{c}.bg.pH_bulk(1,1);
            end
            if ERW{c}.has_update
                PH_c(c) = ERW{c}.PH;      %% pore-water pH for ERW chemistry
                if isfield(ERW{c},'PH_pore') && ~isempty(ERW{c}.PH_pore)
                    PH_c(c) = ERW{c}.PH_pore;
                end
                if isfield(ERW{c},'PH_bulk') && ~isempty(ERW{c}.PH_bulk)
                    PH_bulk_c(c) = ERW{c}.PH_bulk;
                end
                P(j,60,c) = ERW{c}.P60;  P(j,61,c) = ERW{c}.P61;
                P(j,68,c) = ERW{c}.P68;  P(j,69,c) = ERW{c}.P69;
                P(j,52,c) = ERW{c}.P52;  P(j,53,c) = ERW{c}.P53;
                P(j,76,c) = ERW{c}.P76;
                %% The full-pool path replaces P43; the dissolved-fraction
                %% path returns only new mineral weathering P additively.
                if isfield(ERW{c},'P43')
                    P(j,43,c) = ERW{c}.P43;
                end
            end
        end
        clear ERW_next
    end

    %% Rolling checkpoint after a completed hourly iteration. Chemistry
    %% arrays retain the last completed daily state when saved mid-day.
    if ~isempty(ERW_ckpt_file) && ...
            (toc - ERW_ckpt_last_toc) >= ERW_ckpt_minutes*60
        ERW_ckpt_tmp = [tempname(fileparts(ERW_ckpt_file)) '.mat'];
        try
            if FULL_RECORDS
                save(ERW_ckpt_tmp, '-regexp', '^(?!cleanup$).*', '-v7.3');
            else
                %% The state and the daily stores: every variable but the
                %% daily arrays outside the stores, whose last two rows are
                %% all the run reads again (ckpt_daily_rows).
                ckpt_daily_rows = struct();
                for k_ckpt = 1:numel(DAILY_STATE_ONLY)
                    ckpt_daily_rows.(DAILY_STATE_ONLY{k_ckpt}) = ...
                        eval([DAILY_STATE_ONLY{k_ckpt} '(max(1,j-1):j,:,:,:)']);
                end
                ckpt_daily_rows_first = max(1,j-1);
                ckpt_vars = setdiff(who,[{'cleanup'} DAILY_STATE_ONLY]);
                save(ERW_ckpt_tmp, ckpt_vars{:}, '-v7.3', '-nocompression');
                clear ckpt_daily_rows ckpt_vars
            end
            movefile(ERW_ckpt_tmp, ERW_ckpt_file, 'f');
            fprintf(['  [checkpoint] iter %d / %d saved\n'], i, NN);
        catch ME_ckpt
            if isfile(ERW_ckpt_tmp); delete(ERW_ckpt_tmp); end
            warning('ERW:CheckpointFailed', ...
                'Checkpoint at iter %d failed: %s', i, ME_ckpt.message);
        end
        ERW_ckpt_last_toc = toc;
    end
end

%% The last ERW day is summarised now that its rows hold their final values,
%% and each square's daily chemistry store is trimmed to the days run.
if OPT_ERW == 1
    for c = 1:N_col
        [ERW{c},erw_last_day] = ERW_Daily_Record('close',ERW{c});
        if ~isempty(erw_last_day)
            ERW_daily{c} = ERW_Daily_Record('store',ERW_daily{c},erw_last_day,NNd);
        end
        ERW_daily{c} = ERW_Daily_Record('finish',ERW_daily{c});
    end
    clear erw_last_day
end

%% The particle-class histories were kept in daily blocks during the run
%% (ERW_Daily_Step_2D); the results carry the rows each column kept: those
%% set at initialisation and the last row of every ERW day with full
%% records, of the days still held with daily stores.
if OPT_ERW == 1
    for c = 1:N_col
        ERW{c} = ERW_Particle_History('assemble',ERW{c});
    end
end

%% The hours after the last daily step, which no daily row covers.
i_last_day = find(Datam(1:NN,4) == 1,1,'last');
if isempty(i_last_day); i_last_day = 1; end
tail_rows = mod((i_last_day:NN)-1,NH)+1;
Hourly_tail = struct('hours',[i_last_day NN], ...
    'O',Field_Daily_From_Hourly(O(tail_rows,1:n_bio_layers,:),'mean',dth), ...
    'ZWT',Field_Daily_From_Hourly(ZWT(tail_rows,:,:),'mean',dth), ...
    'Lk',Field_Daily_From_Hourly(Lk(tail_rows,:,:),'rate_total',dth), ...
    'f',Field_Daily_From_Hourly(f(tail_rows,:,:),'rate_total',dth), ...
    'q_runon',Field_Daily_From_Hourly(q_runon(tail_rows,:,:),'rate_total',dth), ...
    'Rh',Field_Daily_From_Hourly(Rh(tail_rows,:,:),'total',dth), ...
    'Rd',Field_Daily_From_Hourly(Rd(tail_rows,:,:),'total',dth));
clear tail_rows

Computational_Time = toc;

%% FINAL STATE of a daily-store run: the hours and days the run would read
%% next, in order, which the x-y comparison checks against the last rows of
%% x's full records. The hourly windows themselves are not saved under the
%% hourly names (a reader would take them for 35 years of hours).
if ~FULL_RECORDS
    Final_state = struct('hour_first',max(1,NN-NH+1),'daily_first',max(1,j-1));
    fs_rows = mod((Final_state.hour_first:NN)-1,NH)+1;
    for k_fs = 1:numel(HOURLY_NAMES)
        Final_state.hourly.(HOURLY_NAMES{k_fs}) = eval([HOURLY_NAMES{k_fs} '(fs_rows,:,:)']);
    end
    for k_fs = 1:numel(DAILY_NAMES)
        Final_state.daily.(DAILY_NAMES{k_fs}) = ...
            eval([DAILY_NAMES{k_fs} '(Final_state.daily_first:j,:,:,:)']);
    end
    clear fs_rows k_fs
end
fprintf('\nCOMPUTATIONAL TIME [h] %.3f   (%.1f ms per column-step)\n', ...
    Computational_Time/3600, 1000*Computational_Time/(NN*N_col));

%% ========================================================================
%% DOMAIN AGGREGATION  --  area-weighted over the sub-plot cells, or over
%% the polygon a polygon grid was built on
%% ========================================================================
if field_grid
    %% Polygon grid (v1.0): each cell weighs by its overlap with the
    %% polygon, and the total area includes any polygon land a majority-
    %% rule grid leaves in dropped squares (carried at the domain mean).
    %% The weights are Field_Aggregate_Areal's, so a mean below is v1.0's.
    [A_field,w_col,~,A_field_total] = Field_Aggregation_Areas(topo);
    Domain.area_m2 = A_field_total;
    Domain.simulated_area_m2 = A_tot;
    Domain.grid_outside_area_m2 = A_tot-sum(A_field);
    if isfield(topo,'field_missed_area_m2')
        Domain.field_missed_area_m2 = A_field_total-sum(A_field);
    end
    Domain.area_main_m2   = sum(A_field(strcmp(topo.cell_region,'main')));
    Domain.area_buffer_m2 = sum(A_field(~strcmp(topo.cell_region,'main')));
else
    w_col = reshape(A_cell/A_tot, 1, N_col);
    Domain.area_m2       = A_tot;
    Domain.area_main_m2  = sum(A_cell(strcmp(topo.cell_region,'main')));
    Domain.area_buffer_m2= sum(A_cell(~strcmp(topo.cell_region,'main')));
end
Domain.N2Oflx        = squeeze(N2Oflx)  * w_col.';
Domain.N2flx         = squeeze(N2flx)   * w_col.';
Domain.VOL            = squeeze(VOL)     * w_col.';
Domain.BfixN          = squeeze(BfixN)   * w_col.';
Domain.FertN          = squeeze(FertN_d) * w_col.';
Domain.DepN           = squeeze(DepN_d)  * w_col.';
Domain.NITflx         = squeeze(NITflx)  * w_col.';
Domain.DENITflx       = squeeze(DENITflx)* w_col.';
Domain.N2Oflx_nit     = squeeze(N2Oflx_nit) * w_col.';
Domain.LEAK_NO3      = squeeze(LEAK_NO3)* w_col.';
Domain.LEAK_NH4      = squeeze(LEAK_NH4)* w_col.';
Domain.PH_bulk_cell_d = squeeze(PH_bulk_d);
Domain.PH_pore_cell_d = squeeze(PH_pore_d);
Domain.PH_bulk_d      = pH_Area_Weighted(Domain.PH_bulk_cell_d,w_col);
Domain.PH_pore_d      = pH_Area_Weighted(Domain.PH_pore_cell_d,w_col);
Domain.PH_bulk_aggregation = 'area_weighted_hydrogen_activity_index_not_pooled_extract';
Domain.PH_d           = Domain.PH_bulk_d;
Domain.NO3_pool      = squeeze(P(:,32,:))* w_col.';
Domain.lateral       = lateral_diag;
Domain.exit_water_vert_m3 = exit_water_vert_m3;
Domain.exit_water_surface_m3 = exit_water_surface_m3;
Domain.exit_water_total_m3 = exit_water_vert_m3+exit_water_surface_m3;
if OPT_DEM == 1
    Domain.topography_attached = topography_attached;
    %% Solute exported at the topographic outlet(s) with overland runoff
    %% [g element], one entry per lateral species, and the ERW aqueous
    %% species [umol] ordered as Domain.lateral.erw_species_names.
    Domain.exit_mass_surface_g = exit_mass_surface_g;
    Domain.exit_mass_surface_species = lateral_species;
    Domain.exit_mass_surface_erw_umol = exit_mass_surface_erw_umol;
    fprintf('Overland: %.4g m3 water left at the outlet(s)\n', ...
        sum(exit_water_surface_m3));
end
if field_grid || topography_attached
    %% Exports that originated inside the polygon (v1.0). Routed runoff
    %% leaves through the OUTLET cell, often a partial perimeter cell, so
    %% with topography its field share is the routed tracer, not that
    %% cell's own overlap fraction. On the plot grid every overlap fraction
    %% is 1 and these equal the totals; they are kept for any DEM run so
    %% its exports read the same on either grid.
    Domain.exit_water_field_vert_m3 = sum(exit_water_vert_m3.*field_fraction);
    if topography_attached
        Domain.exit_water_field_surface_m3 = exit_water_field_surface_m3;
        Domain.exit_mass_surface_field_g = exit_mass_surface_field_g;
        Domain.exit_mass_surface_erw_field_umol = exit_mass_surface_erw_field_umol;
    else
        Domain.exit_water_field_surface_m3 = sum(exit_water_surface_m3.*field_fraction);
    end
    Domain.exit_water_field_total_m3 = Domain.exit_water_field_vert_m3 + ...
        Domain.exit_water_field_surface_m3;
end
if field_grid
    Domain.cell_overlap_area_m2 = A_field;
    Domain.cell_overlap_fraction = field_fraction;
    Domain.field_aggregation = 'polygon_overlap_area';
end
Domain.cell_region   = topo.cell_region;
Domain.cell_area_m2  = A_cell;

fprintf('Domain: %.0f m2 total (%.0f main, %.0f margin)\n', ...
    Domain.area_m2, Domain.area_main_m2, Domain.area_buffer_m2);
if field_grid
    fprintf('Domain: %.0f m2 polygon inside a %.0f m2 simulated grid\n', ...
        Domain.area_m2, Domain.simulated_area_m2);
end
fprintf('Cumulative N2O   : %.4f gN/m2\n', sum(Domain.N2Oflx,'omitnan'));
fprintf('Cumulative NO3   : %.4f gN/m2\n', sum(Domain.LEAK_NO3,'omitnan'));
fprintf('Net lateral water margin->field: %.4f m3\n', ...
    sum(lateral_diag.water_buf_to_main,'omitnan'));
fprintf('Net lateral NO3   margin->field: %.4f g\n', ...
    sum(lateral_diag.NO3_buf_to_main,'omitnan'));

%% ========================================================================
function out = ternary(cond, a, b)
if cond; out = a; else; out = b; end
end
