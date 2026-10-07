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
         '  TC_CASE=5 TC_SCENARIO=%s TC_NDAYS=%d matlab -batch "cd(''v2/run''); ERW_run"'], ...
        base_file, scen, n_days);
end
R1 = load(base_file);
if R1.NN ~= n_days*24
    error('test_2d_driver_equivalence:WrongBaselineDuration', ...
        'Baseline contains %d h; expected %d h.',R1.NN,n_days*24);
end
expected_fingerprint = ERW_Equivalence_Fingerprint(v2,n_days);
if ~isfield(R1,'equivalence_fingerprint') || ...
        ~strcmp(R1.equivalence_fingerprint,expected_fingerprint)
    if isfield(R1,'equivalence_fingerprint')
        report_fingerprint_difference(R1.equivalence_fingerprint, ...
            expected_fingerprint);
    else
        fprintf('The saved baseline carries no equivalence fingerprint.\n');
    end
    error('test_2d_driver_equivalence:StaleBaseline', ...
        ['The saved 1-D baseline does not match the current production code, ' ...
         'parameters, inputs and run entry points. Regenerate and fingerprint it.']);
end
%% The 2-D driver below reuses the baseline's ERW_PARAM_input, so a 1-D
%% frame that rewrites a workbook flag would silently move the comparison
%% off the branch production runs.
wb = readtable(fullfile(code,'Inputs','ERW_MOD_PARAM_parameters.xlsx'),'Sheet','ERW_MOD');
if ~isequal(unique(wb.keyword_add),R1.ERW_PARAM_input.keyword_add)
    error('test_2d_driver_equivalence:BackgroundClosure', ...
        'The 1-D baseline ran keyword_add = %g; production reads %s from the workbook.', ...
        R1.ERW_PARAM_input.keyword_add,mat2str(unique(wb.keyword_add)'));
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
%%
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

fprintf('\nERW-coupled chemistry (relative tol = %.0e):\n',tol_chem);
chem_names  = {'P','PH_d','VOL','NITflx','DENITflx','N2Oflx_nit', ...
    'N2flx','N2Oflx','LEAK_NO3','LEAK_NH4','NavlI'};
chem_labels = {'P   BG pools','PH_d pore pH','VOL','NITflx', ...
    'DENITflx','N2Oflx_nit','N2flx','N2Oflx','LEAK_NO3', ...
    'LEAK_NH4','NavlI'};
for k = 1:numel(chem_names)
    ok = compare_field(chem_labels{k},R1.(chem_names{k}),S2.(chem_names{k}),tol_chem);
    if ok; n_pass = n_pass + 1; else; n_fail = n_fail + 1; end
end

fprintf('\nERW mineral-N handoff at ordinary daily calls (relative tol = 1e-10):\n');
BiogeoPar = Biogeochemistry_Parameter(1);
handoff = { ...
    '1-D NH4 boundary',R1.ERW_biogeochem.NH4_tot,R1.P,R1.Datam, ...
        R1.MOD_PARAM_input.t_res,BiogeoPar.aNH4,31; ...
    '1-D NO3 boundary',R1.ERW_biogeochem.NO3_tot,R1.P,R1.Datam, ...
        R1.MOD_PARAM_input.t_res,BiogeoPar.aNO3,32; ...
    '2-D NH4 boundary',S2.ERW{1}.bg.NH4_tot,S2.P,R1.Datam, ...
        S2.ERW{1}.res,BiogeoPar.aNH4,31; ...
    '2-D NO3 boundary',S2.ERW{1}.bg.NO3_tot,S2.P,R1.Datam, ...
        S2.ERW{1}.res,BiogeoPar.aNO3,32};
for k = 1:size(handoff,1)
    ok = compare_handoff(handoff{k,:});
    if ok; n_pass = n_pass + 1; else; n_fail = n_fail + 1; end
end
ok = compare_field('ERW IC_tot',R1.ERW_biogeochem.IC_tot, ...
    S2.ERW{1}.bg.IC_tot,tol_chem);
if ok; n_pass = n_pass + 1; else; n_fail = n_fail + 1; end

fprintf('\n  ---------------------------------------\n');
fprintf('  %d passed, %d failed\n\n', n_pass, n_fail);
if n_fail > 0
    error('test_2d_driver_equivalence: %d comparison(s) failed', n_fail);
end
fprintf('2-D driver reproduces the single-column frame at N_col = 1.\n\n');

end

%% =======================================================================
function ok = compare_handoff(name,total,P,Datam,res,fraction,pool_column)
daily_i = find(Datam(:,4) == 1);
if numel(daily_i) < 2
    fprintf('  FAIL  %-18s fewer than two ordinary daily calls\n',name);
    ok = false; return
end
anchor_rows = round(res*(daily_i(1:end-1)-1));
pool_rows = (2:numel(daily_i)).';
P = squeeze(P);
if any(anchor_rows < 1) || max(anchor_rows) > size(total,1) || ...
        max(pool_rows) > size(P,1)
    fprintf('  FAIL  %-18s boundary indices are outside the saved state\n',name);
    ok = false; return
end
actual = total(anchor_rows,1);
expected = fraction*P(pool_rows,pool_column)/14*1e6;
if any(~isfinite(actual)) || any(actual < 0) || any(~isfinite(expected))
    fprintf('  FAIL  %-18s nonphysical boundary inventory\n',name);
    ok = false; return
end
rel = max(abs(actual(:)-expected(:))./max(1,abs(expected(:))));
ok = rel < 1e-10;
if ok
    fprintf('  PASS  %-18s max relative error %.3e\n',name,rel);
else
    fprintf('  FAIL  %-18s max relative error %.3e\n',name,rel);
end
end

%% =======================================================================
function ok = compare_field(name, a, b, tol)
a = squeeze(a); b = squeeze(b);
if ~isequal(size(a), size(b))
    fprintf('  FAIL  %-18s shape %s vs %s\n', name, mat2str(size(a)), mat2str(size(b)));
    ok = false; return
end
av = a(:); bv = b(:);
%% One-sided NaN, and +Inf against -Inf, are genuine disagreements.  Reject
%% them before any masking: abs(Inf-Inf) and abs(NaN-x) are both NaN, and
%% forcing a NaN maximum to zero further down would report them as a pass.
if any(xor(isnan(av),isnan(bv)))
    fprintf('  FAIL  %-18s NaN in only one field\n', name);
    ok = false; return
end
if any(isinf(av) & isinf(bv) & (sign(av) ~= sign(bv)))
    fprintf('  FAIL  %-18s opposite infinities\n', name);
    ok = false; return
end
d  = abs(av - bv);
d(isnan(av) & isnan(bv)) = 0;              %% NaN in both is agreement
d(isinf(av) & isinf(bv)) = 0;              %% same-signed Inf is agreement
finite_b = bv(isfinite(bv));
if isempty(finite_b); sc = 1; else; sc = max(1, max(abs(finite_b))); end
md = max(d);
if isempty(md); md = 0; end
rel = md/sc;
if rel <= tol
    fprintf('  PASS  %-18s max|d| = %.3e  (rel %.2e)\n', name, md, rel);
    ok = true;
else
    [~,iw] = max(d);
    first_bad = find(d/sc > tol,1);
    fprintf(['  FAIL  %-18s max|d| = %.3e  (rel %.2e) at linear index %d; ' ...
        'first crossing %d: 1-D %.12g, 2-D %.12g, |d| %.3e\n'], ...
        name, md, rel, iw, first_bad, av(first_bad), bv(first_bad), d(first_bad));
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
    ... %% The single-column frame applies the field record, 2016-2019.
    'last_application_year', 2019, ...
    ... %% No replaced literals and the native saturation: the ensemble's
    ... %% defaults, which must leave the frame bit-identical to the 1-D one.
    'biogeo_overrides', [], ...
    'saturation_exponent', 1, ...
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
    'EG',EG,'LAI_L',LAI_L,'B_L',B_L,'P',P,'PH_d',PH_d,'VOL',VOL, ...
    'NITflx',NITflx,'DENITflx',DENITflx,'N2Oflx_nit',N2Oflx_nit, ...
    'N2flx',N2flx,'N2Oflx',N2Oflx,'LEAK_NO3',LEAK_NO3, ...
    'LEAK_NH4',LEAK_NH4,'NavlI',NavlI,'ERW',{ERW});
end

function report_fingerprint_difference(saved,current)
%% Name what moved between the baseline build and this checkout, so a
%% stale baseline is traced to a file instead of rebuilt by guess.
saved_parts = strsplit(char(saved),'|');
current_parts = strsplit(char(current),'|');
only_saved = setdiff(saved_parts,current_parts);
only_current = setdiff(current_parts,saved_parts);
fprintf('Equivalence fingerprint mismatch: %d baseline entries and %d checkout entries differ.\n', ...
    numel(only_saved),numel(only_current));
for k = 1:numel(only_saved)
    fprintf('  baseline had : %s\n',shorten_hash(only_saved{k}));
end
for k = 1:numel(only_current)
    fprintf('  checkout has : %s\n',shorten_hash(only_current{k}));
end
end

function entry = shorten_hash(entry)
%% Keep the path and the first 12 hex digits of a content hash.
entry = regexprep(entry,'(:[0-9a-f]{12})[0-9a-f]{52}$','$1...');
end
