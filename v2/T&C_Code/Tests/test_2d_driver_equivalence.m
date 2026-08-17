function test_2d_driver_equivalence(n_days)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% test_2d_driver_equivalence
%%
%% The decisive test of the multi-column port.
%%
%% MAIN_FRAME_BG_ERW_2D calls the same unit functions as the
%% single-column MAIN_FRAME_BG_ERW_10mins; the only differences are that
%% state carries a trailing column dimension and that columns exchange
%% water and solute laterally.  On a one-cell topology the lateral block
%% is skipped entirely, so the two frames must produce the SAME numbers.
%% Any divergence is a transcription error in the driver, not physics.
%%
%% Method
%%   1. Run the 1-D frame for a short window via ERW_run and save it.
%%   2. Reload that workspace, rebuild the caller context, and run the
%%      2-D driver on a single-cell topology over the same forcing.
%%   3. Compare state trajectories element by element.
%%
%% Usage
%%   test_2d_driver_equivalence          % 10 simulated days
%%   test_2d_driver_equivalence(20)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 1 || isempty(n_days); n_days = 10; end

here = fileparts(mfilename('fullpath'));
v2   = fileparts(fileparts(here));          % .../v2
code = fullfile(v2,'T&C_Code');
addpath(code, fullfile(code,'Inputs'), fullfile(code,'Parameters'));

fprintf('\n=== test_2d_driver_equivalence (%d days) ===\n', n_days);

%% ---------------------------------------------------------------------
%% 1. Baseline: single-column frame
%% ---------------------------------------------------------------------
scen = 'EQUIV_1D_TMP';
base_file = fullfile(v2,'SIM',scen,'TC_ERW_results_Kantola_crop_control.mat');
if ~isfile(base_file)
    error('test_2d_driver_equivalence:NoBaseline', ...
        ['Baseline not found: %s\nGenerate it first with:\n' ...
         '  TC_CASE=5 TC_SCENARIO=%s TC_NDAYS=%d matlab -batch ERW_run'], ...
        base_file, scen, n_days);
end
R1 = load(base_file);
if R1.NN ~= n_days*24
    error('test_2d_driver_equivalence:WrongBaselineDuration', ...
        'Baseline contains %d h; expected %d h.',R1.NN,n_days*24);
end
fprintf('Baseline loaded (%d h simulated)\n', R1.NN);

%% ---------------------------------------------------------------------
%% 2. Run the 2-D driver at N_col = 1 over identical forcing
%% ---------------------------------------------------------------------
fprintf('Running 2-D driver at N_col = 1...\n');
S2 = run_driver_single(R1, code);

%% ---------------------------------------------------------------------
%% 3. Compare
%% ---------------------------------------------------------------------
%% TWO TOLERANCE CLASSES, for two different reasons.
%%
%%   EXACT (tol = 0).  Water, energy and vegetation state runs through
%%   identical code with identical inputs, so it must agree bit for bit.
%%   Anything else is a transcription error in the driver.
%%
%%   ROUND-OFF (tol = 1e-8 relative).  The ERW chain cannot be held to
%%   bit-exactness once the forcing interpolation is computed a different
%%   way.  interp1 and synchronize both perform linear interpolation but
%%   order the arithmetic differently, so they differ in the last ulp.
%%   ERW_function drives that through fsolve, whose iteration path
%%   depends on the exact input bits, amplifying an ulp to ~1e-9 in pH.
%%   That is nine significant figures of agreement in a quantity measured
%%   to +/-0.1, and it does not grow into a physical difference.
%%   Demanding max|d| = 0 here would be demanding that a nonlinear
%%   iterative solver be bit-reproducible across a change in input
%%   rounding, which it is not.
tol_exact = 0;
tol_chem  = 1e-8;

n_pass = 0; n_fail = 0;

fprintf('\nExact-match state (water, energy, vegetation), tol = 0:\n');
exact_names  = {'O','V','Ts','Tdp','Lk','T_H','T_L','EG','LAI_L','B_L'};
exact_labels = {'O   soil moisture','V   soil water','Ts  surf temp', ...
                'Tdp soil temp','Lk  leakage','T_H transpiration', ...
                'T_L transpiration','EG  evaporation','LAI_L','B_L biomass'};
for k = 1:numel(exact_names)
    ok = compare_field(exact_labels{k}, R1.(exact_names{k}), S2.(exact_names{k}), tol_exact);
    if ok; n_pass = n_pass + 1; else; n_fail = n_fail + 1; end
end

fprintf('\nERW-coupled chemistry, tol = %.0e relative (solver round-off):\n', tol_chem);
chem_names  = {'P','PH_d','N2Oflx','LEAK_NO3','LEAK_NH4','NavlI'};
chem_labels = {'P   BG pools','PH_d pore pH','N2Oflx','LEAK_NO3','LEAK_NH4','NavlI'};
for k = 1:numel(chem_names)
    ok = compare_field(chem_labels{k}, R1.(chem_names{k}), S2.(chem_names{k}), tol_chem);
    if ok; n_pass = n_pass + 1; else; n_fail = n_fail + 1; end
end

fprintf('\n  ---------------------------------------\n');
fprintf('  %d passed, %d failed\n\n', n_pass, n_fail);
if n_fail > 0
    error('test_2d_driver_equivalence: %d comparison(s) failed', n_fail);
end
fprintf('2-D driver reproduces the single-column frame at N_col = 1.\n\n');

end

%% =======================================================================
function ok = compare_field(name, a, b, tol)
a = squeeze(a); b = squeeze(b);
if ~isequal(size(a), size(b))
    fprintf('  FAIL  %-18s shape %s vs %s\n', name, mat2str(size(a)), mat2str(size(b)));
    ok = false; return
end
d  = abs(a(:) - b(:));
d(isnan(a(:)) & isnan(b(:))) = 0;          %% NaN in both is agreement
sc = max(1, max(abs(b(:))));
md = max(d);
if isempty(md) || isnan(md); md = 0; end
rel = md/sc;
if rel <= tol
    fprintf('  PASS  %-18s max|d| = %.3e  (rel %.2e)\n', name, md, rel);
    ok = true;
else
    [~,iw] = max(d);
    fprintf('  FAIL  %-18s max|d| = %.3e  (rel %.2e) first at linear index %d\n', ...
        name, md, rel, iw);
    ok = false;
end
end

%% =======================================================================
function S2 = run_driver_single(R1, code)
%% Rebuild the caller context the 2-D driver expects and run it.  This is
%% a plain local function (no nested functions), so its workspace is not
%% static and the driver script may declare its own variables here.

NN = R1.NN; NNd = R1.NNd; ms = R1.ms; cc = R1.cc; Zs = R1.Zs; %#ok<NASGU>
dt = R1.dt; dth = R1.dth;                                      %#ok<NASGU>
Datam = R1.Datam; Date = R1.Date;                              %#ok<NASGU>
Pr=R1.Pr; Ta=R1.Ta; Ds=R1.Ds; Ws=R1.Ws; ea=R1.ea; N=R1.N; Pre=R1.Pre; %#ok<NASGU>
Tdew=R1.Tdew; SAB1=R1.SAB1; SAB2=R1.SAB2; SAD1=R1.SAD1; SAD2=R1.SAD2;  %#ok<NASGU>
PARB=R1.PARB; PARD=R1.PARD; Ca=R1.Ca; Oa=R1.Oa;                %#ok<NASGU>
Lat=R1.Lat; Lon=R1.Lon; DeltaGMT=R1.DeltaGMT;                  %#ok<NASGU>
t_bef=R1.t_bef; t_aft=R1.t_aft;                                %#ok<NASGU>
id_location = R1.id_location;                                  %#ok<NASGU>
NN_case_study_ERW_MOD = R1.NN_case_study_ERW_MOD;              %#ok<NASGU>
MOD_PARAM_input = R1.MOD_PARAM_input;
ERW_PARAM_input = R1.ERW_PARAM_input;
OPT_ERW = R1.OPT_ERW;                                          %#ok<NASGU>
OPT_SoilBiogeochemistry = R1.OPT_SoilBiogeochemistry;          %#ok<NASGU>
if isfield(R1,'B_out');      B_out = R1.B_out;           end   %#ok<NASGU>
if isfield(R1,'O_out');      O_out = R1.O_out;           end   %#ok<NASGU>
if isfield(R1,'NavlI_out');  NavlI_out = R1.NavlI_out;   end   %#ok<NASGU>
if isfield(R1,'RexmyI_out'); RexmyI_out = R1.RexmyI_out; end   %#ok<NASGU>
if isfield(R1,'NupI_L_out'); NupI_L_out = R1.NupI_L_out; end   %#ok<NASGU>
if isfield(R1,'IrD');        IrD = R1.IrD;               end   %#ok<NASGU>

Parameter_Directory = fullfile(code,'Parameters');
topo = Build_Single_Cell_Topology(3600);                       %#ok<NASGU>
ColCfg = struct( ...
    'PARAM_IC', fullfile(Parameter_Directory, [char(MOD_PARAM_input.MOD_PARAM) '.m']), ...
    'PFT_OVERRIDE', '', ...
    ... %% Mirror the single-column frame's own predicate, which applies rock
    ... %% when diss_f > 0. M_rock_in is the mesocosm's initial rock mass and
    ... %% is 0 for the Kantola field setup, which adds 5000 g/m2 on the
    ... %% time_application schedule instead, so M_rock_in > 0 reported no
    ... %% rock for basalt and calcite. That left cfg.apply_rock false while
    ... %% the chemistry row still carried a feedstock, which is the one
    ... %% combination the 2-D parameter assembly cannot handle.
    'apply_rock', logical(R1.OPT_ERW && ERW_PARAM_input.diss_f > 0)); %#ok<NASGU>

old = cd(code);
cleanup = onCleanup(@() cd(old));
MAIN_FRAME_BG_ERW_2D;

S2 = struct('O',O,'V',V,'Ts',Ts,'Tdp',Tdp,'Lk',Lk,'T_H',T_H,'T_L',T_L, ...
    'EG',EG,'LAI_L',LAI_L,'B_L',B_L,'P',P,'PH_d',PH_d,'N2Oflx',N2Oflx, ...
    'LEAK_NO3',LEAK_NO3,'LEAK_NH4',LEAK_NH4,'NavlI',NavlI);
end
