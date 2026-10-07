function test_cec_pH_dependence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Empirical pH-dependent capacity for TC_GT, with a fixed Appelo capacity.
%%
%% SCEPTER assigns fixed mass-specific CEC to each exchanger solid and varies
%% apparent selectivity with proton loading.  The lumped TC Appelo treatment
%% has no evolving solid inventory, so it uses fixed soil CEC.  The linear
%% site-density relation remains a separately identified TC_GT extension.
%%
%% The biogeochemistry is vertically lumped into a single layer of depth
%% Zbio, so there is one exchanger and one bulk sensitivity.  The parameter
%% is therefore declared mass-specific (cmolc kg-1 per pH unit) and
%% converted onto the areal basis with exactly the factor that built
%% CEC_tot, which keeps it invariant if Zbio changes.
%%
%% Absent, empty or zero must reproduce fixed-CEC behaviour.  A positive
%% Appelo value must fail before initialization or daily exchange begins.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n=== test_cec_pH_dependence ===\n');
n_pass = 0;
here = fileparts(mfilename('fullpath'));
code_dir = fileparts(here);
addpath(code_dir);

%% ---------------------------------------------- 2-D override initialization
frame_source = fileread(fullfile(code_dir,'MAIN_FRAME_BG_ERW_2D.m'));
daily_source = fileread(fullfile(code_dir,'ERW_Daily_Step_2D.m'));
parameter_source = fileread(fullfile(code_dir,'Parameters','ERW_MOD_PARAM_excel.m'));
run_source = fileread(fullfile(fileparts(code_dir),'run','Run_One_Subfield.m'));
override_2d_pos = strfind(frame_source, ...
    'PAR{c}.bulk_pH_params.dCEC_dpH = ColCfg(c).cec_dpH;');
daily_call_pos = strfind(frame_source,'ERW_next{c} = ERW_Daily_Step_2D');
check(~isempty(override_2d_pos) && ~isempty(daily_call_pos) && ...
        override_2d_pos(1) < daily_call_pos(1), ...
    '2-D CEC sensitivity reaches the ERW column before initialization');
daily_copy_pos = strfind(daily_source,'C.bulk_pH_params = p.bulk_pH_params;');
daily_init_pos = strfind(daily_source, ...
    'ERW_MOD_PARAM_excel(t, TT_init, E.TC_const');
check(~isempty(daily_copy_pos) && ~isempty(daily_init_pos) && ...
        daily_copy_pos(1) < daily_init_pos(1), ...
    'ERW initialization receives the overridden bulk-pH parameters');
merge_pos = strfind(parameter_source,'fn_bulk = fieldnames(ERW_TC_const_input.bulk_pH_params);');
invert_pos = strfind(parameter_source,'] = Invert_Bulk_pH(');
check(~isempty(merge_pos) && ~isempty(invert_pos) && merge_pos(1) < invert_pos(1), ...
    '2-D bulk-pH overrides are merged before initial pH inversion');
check(contains(parameter_source, ...
        "(~isfield(bulk_pH_params,'CEC_ref_pH') || isempty(bulk_pH_params.CEC_ref_pH))"), ...
    'an empty CEC reference is initialized on the pore-water pH basis');
check(contains(run_source,'run_provenance.CECpHSensitivity') && ...
        contains(run_source,'run_provenance.CECrefpH') && ...
        contains(run_source,'run_provenance.CEC_cmolc_kg') && ...
        contains(run_source,'cec_params = ERW{c}.const.bulk_pH_params') && ...
        contains(run_source,'ERW_Resolve_CEC_pH_Sensitivity(cec_params,cec_scheme)'), ...
    'saved results retain the effective CEC slope and pore-pH reference');

%% ---------------------------------------------------------------- resolver
[d,pref,dg] = ERW_Resolve_CEC_pH_Sensitivity(struct(),'appelo');
check(d == 0 && ~dg.pH_dependent, 'absent dCEC_dpH is legacy constant CEC');
check(isnan(pref), 'constant CEC does not invent a pH reference');

d = ERW_Resolve_CEC_pH_Sensitivity(struct('dCEC_dpH',[]),'appelo');
check(d == 0, 'empty dCEC_dpH is legacy constant CEC');

[d,~,dg] = ERW_Resolve_CEC_pH_Sensitivity( ...
    struct('dCEC_dpH',2.9,'CEC_ref_pH',5.4),'TC_GT');
check(d == 2.9 && dg.pH_dependent, ...
    'scalar dCEC_dpH is honoured under TC_GT');
check(threw_id(@() ERW_Resolve_CEC_pH_Sensitivity( ...
    struct('dCEC_dpH',2.9,'CEC_ref_pH',5.4),'appelo'), ...
    'ERW:scepterFixedCEC'), ...
    'Appelo rejects a linear change in total site density');

spec = struct('dCEC_dpH',struct('appelo',2.9,'TC_GT',1.4),'CEC_ref_pH',5.4);
check(threw_id(@() ERW_Resolve_CEC_pH_Sensitivity(spec,'appelo'), ...
    'ERW:scepterFixedCEC'), ...
    'scheme-keyed Appelo sensitivity is rejected');
check(ERW_Resolve_CEC_pH_Sensitivity(spec,'TC_GT') == 1.4, ...
    'scheme-keyed dCEC_dpH selects the TC_GT entry under TC_GT');

only_app = struct('dCEC_dpH',struct('appelo',2.9),'CEC_ref_pH',5.4);
check(threw_id(@() ERW_Resolve_CEC_pH_Sensitivity(only_app,'appelo'), ...
    'ERW:scepterFixedCEC'), ...
    'positive Appelo entry cannot enable empirical variable CEC');
check(ERW_Resolve_CEC_pH_Sensitivity(only_app,'TC_GT') == 0, ...
    'scheme absent from the struct falls back to constant CEC');

[~,pref] = ERW_Resolve_CEC_pH_Sensitivity( ...
    struct('dCEC_dpH',2.9,'CEC_ref_pH',5.4),'TC_GT');
check(pref == 5.4, 'CEC_ref_pH override is honoured');
check(threw(@() ERW_Resolve_CEC_pH_Sensitivity( ...
    struct('dCEC_dpH',2.9),'TC_GT')), ...
    'variable CEC without a declared pH basis is rejected');

check(threw(@() ERW_Resolve_CEC_pH_Sensitivity(struct('dCEC_dpH',-1),'appelo')), ...
    'negative dCEC_dpH is rejected');
check(threw(@() ERW_Resolve_CEC_pH_Sensitivity(struct('dCEC_dpH',[1 2]),'appelo')), ...
    'non-scalar dCEC_dpH is rejected');
check(threw(@() ERW_Resolve_CEC_pH_Sensitivity( ...
    struct('dCEC_dpH',struct('appelo',2.9,'gaines',1)),'appelo')), ...
    'unknown scheme field in the struct is rejected');
check(threw(@() ERW_Resolve_CEC_pH_Sensitivity( ...
    struct('dCEC_dpH',2.9,'CEC_ref_pH',12),'TC_GT')), ...
    'out-of-range CEC_ref_pH is rejected');
check(threw(@() ERW_Resolve_CEC_pH_Sensitivity(struct('dCEC_dpH',2.9),'GT')), ...
    'unknown exchange scheme is rejected');

%% ------------------------------------------- areal conversion is Zbio-free
%% ERW_function converts with CEC_tot/CEC_cmolc_kg.  That must equal the
%% factor ERW_MOD_PARAM_excel used to build CEC_tot, so the same declared
%% cmolc kg-1 per pH gives the same molc m-2 per pH at any layer depth.
CEC_cmolc_kg = 13;  rho_bulk = 1295*1e3;  conv_mol = 1e6;  dCEC = 2.9;
for Zrw = [0.30 0.15 0.50]
    CEC_tot_z = CEC_cmolc_kg*1e-5*rho_bulk*Zrw*conv_mol;
    via_ratio = dCEC*(CEC_tot_z/CEC_cmolc_kg);
    via_factor = dCEC*1e-5*rho_bulk*Zrw*conv_mol;
    check(abs(via_ratio-via_factor) <= 1e-6*via_factor, ...
        sprintf('areal conversion matches the CEC_tot factor at Zbio %.0f mm', Zrw*1000));
    check(abs(via_ratio/CEC_tot_z - dCEC/CEC_cmolc_kg) <= 1e-12, ...
        sprintf('fractional sensitivity is depth-invariant at Zbio %.0f mm', Zrw*1000));
end

%% ----------------------------------------------- capacity is not clamped
CEC_ref = 5e7;
dCEC_areal = dCEC*(CEC_ref/CEC_cmolc_kg);
H_ref = 10^(-5.4)*conv_mol;
check(ERW_CEC_At_H(H_ref,conv_mol,CEC_ref,dCEC_areal,5.4) == CEC_ref, ...
    'effective CEC equals the prescribed capacity at its reference pH');
H_low = conv_mol;
check(ERW_CEC_At_H(H_low,conv_mol,CEC_ref,dCEC_areal,5.4) < 0, ...
    'out-of-domain linear CEC remains nonpositive for the physical gate to reject');

%% ------------------------------------------------- the site calibration
%% This pins the conversion arithmetic and a first-order estimate of the
%% sensitivity the Energy Farm needs.  It deliberately does NOT pin a
%% calibrated value: dCEC_dpH is determined by sweeping CECpHSensitivity
%% against the observed basalt response, and this test must not encode a
%% guess as though it were a result.
soil_kg_m2 = 1295*0.30;
check(abs(soil_kg_m2-388.5) < 0.1, 'lumped layer carries 388.5 kg soil m-2');
per_cmolc = 1e-2*soil_kg_m2;               %% molc m-2 per pH per cmolc kg-1
check(abs(per_cmolc-3.885) < 1e-9, ...
    '1 cmolc kg-1 per pH is 3.885 molc m-2 per pH in the lumped layer');

%% Buffer capacity scales inversely with the pH excursion for a fixed
%% alkalinity input, so matching observations needs the implemented
%% capacity multiplied by dpH_modelled/dpH_observed.
buffer_now = 12.7;  dpH_mod = 0.996;  dpH_obs = 0.364;
need = buffer_now*dpH_mod/dpH_obs - buffer_now;
est = need/per_cmolc;
check(est > 1 && est < 10, sprintf( ...
    'implied sensitivity %.1f cmolc kg-1 per pH is within the variable-charge range', est));
fprintf('  first-order target: %.1f molc m-2 per pH extra = %.1f cmolc kg-1 per pH\n', ...
    need, est);

%% ----------------------------------------- numerical: constant vs variable
%% Reuses the mass-conserving Appelo fixture, which is self-contained.
[K_Na_H,K_K_H,K_Ca_H2,K_Mg_H2,K_Al_H3,alpha] = Resolve_Appelo_Constants([]);
CEC_tot = 3.34e7;
soil_water_l = 61.9;
s = 0.4034;
bulk_volume_l = soil_water_l/s;
cation_tot = [1.04525e7,2.42372e6,1.81797e6,2.45141e4];
carbonate = [0.447,4.68e-5,0.01,0.034,1e6];
appelo = [K_Na_H,K_K_H,K_Ca_H2,K_Mg_H2,K_Al_H3,alpha];
gas = [25,15,1/144];
al3_fn = @(H) 10^(-1.2)*H.^3;
k1 = carbonate(1); k2 = carbonate(2); k_w = carbonate(3); k_H = carbonate(4);

H_ref = 0.8;  CO2_ref = 66.3;
pH_ref = -log10(H_ref/carbonate(5));
f_ref = ERW_Appelo_State_At_H(H_ref,cation_tot,al3_fn(H_ref), ...
    CEC_tot,soil_water_l,alpha,K_Na_H,K_K_H,K_Ca_H2,K_Mg_H2,K_Al_H3);
Alk_ref = k1*CO2_ref/H_ref + 2*k1*k2*CO2_ref/H_ref^2 - H_ref + k_w/H_ref;
Alk_tot = sum(f_ref(1:4))*CEC_tot + Alk_ref*soil_water_l;
IC_ref = CO2_ref*bulk_volume_l*(s*(1+k1/H_ref+k1*k2/H_ref^2)+(1-s)/k_H);
IC_pre = IC_ref + gas(1)*(CO2_ref/k_H-gas(2))*gas(3);

%% Referenced at the baseline pH so the unamended state is untouched and
%% only the amended state can differ.
frac = dCEC/CEC_cmolc_kg;
cec_var = @(H) ERW_CEC_At_H(H,carbonate(5),CEC_tot,frac*CEC_tot,pH_ref);

solve = @(A,fn) ERW_Appelo_Mass_Conserving_State(A,IC_pre,cation_tot, ...
    soil_water_l,bulk_volume_l,s,carbonate,CEC_tot,appelo,gas,al3_fn,[],[],fn);

st_legacy = ERW_Appelo_Mass_Conserving_State(Alk_tot,IC_pre,cation_tot, ...
    soil_water_l,bulk_volume_l,s,carbonate,CEC_tot,appelo,gas,al3_fn);
st_const = solve(Alk_tot,@(H) CEC_tot);
check(st_legacy.is_valid && st_const.is_valid, 'baseline states solve');
check(st_legacy.H == st_const.H, ...
    'an explicit constant cec_fn is bit-identical to the legacy call');

st_var = solve(Alk_tot,cec_var);
check(st_var.is_valid, 'variable-charge baseline solves');
check(abs(-log10(st_var.H/carbonate(5)) + log10(st_const.H/carbonate(5))) < 1e-6, ...
    'at CEC_ref_pH the variable exchanger reproduces the constant one');

%% Perturb the conserved alkalinity pool and compare the pH excursion.
%% Variable charge buffers in BOTH directions: sites appear as pH rises and
%% absorb base, and disappear as pH falls, releasing the base they held. So
%% either sign of perturbation exercises the same mechanism.
%%
%% Perturbation size must respect what the fixture can actually hold. Its
%% total base-cation charge is 2Ca+2Mg+K+Na, and the baseline already parks
%% nearly all of it on the exchanger, so there is very little headroom to ADD
%% base and a great deal to remove it. Sizing the perturbation in absolute
%% umol_c rather than relative to that inventory asks for a state with more
%% exchanger charge than there are cations to supply, and every solve
%% correctly returns is_valid == false.
cat_charge = 2*cation_tot(1)+2*cation_tot(2)+cation_tot(3)+cation_tot(4);
mags = cat_charge*[3e-4 1e-3 3e-3 1e-2 3e-2 1e-1];
cands = [-mags, mags];
best = [];  n_tried = 0;  n_solved = 0;
for dAlk = cands
    n_tried = n_tried + 1;
    a = solve(Alk_tot+dAlk,@(H) CEC_tot);
    b = solve(Alk_tot+dAlk,cec_var);
    if ~(a.is_valid && b.is_valid), continue, end
    n_solved = n_solved + 1;
    dc = -log10(a.H/carbonate(5)) - pH_ref;
    dv = -log10(b.H/carbonate(5)) - pH_ref;
    if ~isfinite(dc) || ~isfinite(dv), continue, end
    if abs(dc) >= 0.02 && abs(dc) <= 2.0
        if isempty(best) || abs(dc) > abs(best(1))
            best = [dc dv dAlk];
        end
        if abs(dc) >= 0.2, break, end
    end
end
check(~isempty(best), sprintf( ...
    'found a usable pH excursion (%d perturbations tried, %d solved)', ...
    n_tried, n_solved));
dpH_const = best(1);  dpH_var = best(2);  dAlk_used = best(3);
check(abs(dpH_var) < abs(dpH_const), sprintf( ...
    'variable charge damps the response (%+.3f vs %+.3f pH)', dpH_var, dpH_const));
check(sign(dpH_var) == sign(dpH_const), ...
    'damping does not reverse the direction of the response');
fprintf('  dAlk %+.4g umol_c: dpH constant %+.3f, variable %+.3f (ratio %.2f)\n', ...
    dAlk_used, dpH_const, dpH_var, dpH_var/dpH_const);

%% Sites must still close and cation mass must still be conserved against a
%% site density that moved. Asserted, not skipped: a silent skip here would
%% let the whole conservation argument pass vacuously.
b = solve(Alk_tot+dAlk_used,cec_var);
check(b.is_valid, 'the chosen perturbation re-solves under variable charge');
cec_at = cec_var(b.H);
z = [2 2 1 1];
inv = b.concentration.*soil_water_l + b.fraction(1:4).*cec_at./z;
check(max(abs(inv-cation_tot)./cation_tot) <= 1e-9, ...
    'cation inventory closes against the pH-dependent site density');
check(abs(sum(b.fraction)-1) <= 1e-10, ...
    'exchange sites close against the pH-dependent site density');
if dpH_var > 0
    check(cec_at > CEC_tot, 'site density rose with pH, as intended');
else
    check(cec_at < CEC_tot, 'site density fell with pH, as intended');
end
check(abs(cec_at-CEC_tot)/CEC_tot > 1e-6, ...
    'the site density actually moved, so the comparison is not vacuous');

fprintf('  ---------------------------------------\n');
fprintf('  %d passed, 0 failed\n', n_pass);

    function check(ok, label)
        assert(ok, 'test_cec_pH_dependence:Failed', '%s', label);
        n_pass = n_pass + 1;
        fprintf('  PASS  %s\n', label);
    end
end

function tf = threw(fn)
tf = false;
try
    fn();
catch
    tf = true;
end
end

function tf = threw_id(fn,id)
tf = false;
try
    fn();
catch ME
    tf = strcmp(ME.identifier,id);
end
end
