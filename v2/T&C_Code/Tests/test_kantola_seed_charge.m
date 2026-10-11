function test_kantola_seed_charge
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Charge and ledger consistency of the Kantola ERW seed.
%%
%% The seed must hand the daily solver a pore water that can exist.  Two
%% things made that impossible, and both were live bugs:
%%
%%   1. Handing ERW the whole labile mineral pool instead of its
%%      dissolved fraction.  P(:,43) carries ~183 gP/m2 of labile P; all
%%      of it in solution is 54 mM phosphate, so anions outweigh cations
%%      and no nonnegative counter-anion exists at any pH.  The frame now
%%      applies the model's own solubility coefficients (BiogeoPar.aP,
%%      aNO3, aNH4) -- the same ones Biogeo_Leakage and Biogeo_uptake
%%      already apply to those pools.
%%
%%   2. Assigning exchange fractions by fiat.  The old convention gave
%%      sodium every site the prescribed Ca/Mg/K did not claim (f_Na = 0.76,
%%      an ESP of 76%) and inverted Gaines-Thomas for the concentration,
%%      producing 0.5 M Na whose only job was balancing that phosphate.
%%      The seed now solves the Appelo site sum against the spin-up
%%      solution, with Al set by Al(OH)3 solubility.
%%
%% This runs the REAL excel seed path in the production Kantola
%% configuration and asserts what the daily solver needs on day 1: a
%% nonnegative counter-anion that closes WITHOUT its max(An,0) floor, a
%% positive alkalinity pool, and the ledger identity
%%   Alk_tot(1) = Alk(1)*Vw + R_alk(1)
%% telescoping exactly.  A floored An breaks that identity and the solver
%% then inherits an Alk_tot it can only host as a hyper-acid f_Al = 1
%% state.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n=== test_kantola_seed_charge ===\n');
n_pass = 0;

here = fileparts(mfilename('fullpath'));
v2 = fileparts(fileparts(here));
addpath(fullfile(v2, 'T&C_Code'), fullfile(v2, 'T&C_Code', 'Parameters'));

[cec_probe,cec_diag] = ERW_Resolve_CEC(8.6);
check(cec_probe == 8.6 && ~cec_diag.measured && ~cec_diag.estimated && ...
    strcmp(cec_diag.source,"explicit_prescribed_CEC"), ...
    'Kantola CEC is retained as a prescribed, not measured, boundary');

simfile = fullfile(v2, 'SIM', 'EQUIV_1D_TMP', ...
    'TC_ERW_results_Kantola_crop_control.mat');
if exist(simfile, 'file') ~= 2
    fprintf('  SKIP  baseline state %s not found\n', simfile);
    return
end
S = load(simfile);

t = 1;
TT = S.TT_m(1,:);
ERW_TC_const = struct();
ERW_TC_const.n = S.ERW_const.n;
ERW_TC_const.rsd_bio = S.ERW_const.rho_bulk/1e3;
ERW_TC_const.conv_mol = S.conv_mol;
ERW_TC_const.conv_Al = S.conv_Al;
ERW_TC_const.PH_in = 6.1;
ERW_TC_const.Zbio = S.ERW_const.Zrw*1000;
ERW_TC_const.M_rock_add = 0;
ERW_TC_const.M_rock_add_CaCO3 = 0;
ERW_TC_const.cov_moll = 1;
% The one-exchanger implementation preserves the SCEPTER relative
% selectivities and matches the observed pretreatment pH/base saturation.
fixed_logKH = [4.71627429367 3.61627429367 8.10254858734 ...
    8.41854858734 12.918822881];
ERW_TC_const.bulk_pH_params = struct('scheme', 'appelo', ...
    'Al_solubility_pin', true, 'initial_pH_mode', 'bulk_observed', ...
    'appelo', struct('logKH', fixed_logKH, 'alpha', 1.3));

% The frame hands ERW only the DISSOLVED part of the labile mineral
% pools (MAIN_FRAME_BG_ERW_10mins / ERW_Daily_Step_2D apply the model's
% own solubility coefficients).  Mirror that here, or the harness feeds
% the solver a pore water that no charge balance can close.
bp_handoff = Biogeochemistry_Parameter(1);
ERW_biogeochem = struct();
ERW_biogeochem.Mg(1) = 2268.73; ERW_biogeochem.Si(1) = 1;
ERW_biogeochem.Ca(1) = 530.237; ERW_biogeochem.K(1) = 7675.04;
ERW_biogeochem.PO4(1) = bp_handoff.aP   * 54053.8;
ERW_biogeochem.NO3(1) = bp_handoff.aNO3 * 481.263;
ERW_biogeochem.NH4(1) = bp_handoff.aNH4 * 3088.98;
ERW_biogeochem.excCa = 530.237; ERW_biogeochem.excMg = 2268.73;
ERW_biogeochem.excK = 7675.04;

xl = fullfile(v2, 'T&C_Code', 'Inputs', 'ERW_MOD_PARAM_parameters.xlsx');
ERW_PARAM_input = readtable(xl, 'Sheet', 'ERW_MOD');
ERW_PARAM_input = ERW_PARAM_input(4,:);          % case 5 (Kantola control)
ERW_PARAM_input.f_CEC_in = string("NaN");
for rc = {'rock_filepath','rock_type','mineral','rock_f_in'}
    if ismember(rc{1}, ERW_PARAM_input.Properties.VariableNames)
        ERW_PARAM_input.(rc{1}) = string("NaN");
    end
end
erw_cols = {'Al_w_in','Si_in','Ca_in','Mg_in','K_in','Na_in','PO4_in','NO3_in','NH4_in'};
for c = 1:numel(erw_cols)
    if ~any(strcmp(ERW_PARAM_input.Properties.VariableNames, erw_cols{c}))
        ERW_PARAM_input.(erw_cols{c}) = NaN;
    end
end
ERW_PARAM_input.opt_TC_biochem = 1;
ERW_PARAM_input.num_rock_app = sum(ERW_TC_const.M_rock_add > 0);

ERW_PARAM_bad = ERW_PARAM_input;
ERW_PARAM_bad.Na_in = 1e9;
assert_error(@() ERW_MOD_PARAM_excel( ...
    t, TT, ERW_TC_const, ERW_biogeochem, ERW_PARAM_bad, 1), ...
    'ERW:InfeasibleKantolaNaExchange');
n_pass = n_pass + 1;

[B, ~, ~, EC, ~, bg_names] = ERW_MOD_PARAM_excel( ...
    t, TT, ERW_TC_const, ERW_biogeochem, ERW_PARAM_input, 1);
check(isfield(EC.bulk_pH_params, 'Al_solubility_pin') && ...
    isequal(EC.bulk_pH_params.Al_solubility_pin, true), ...
    'pin flag survives the seed');

% Field-water volume from a porewater-only pool (NH4_tot = NH4*Vw).
Vw = B.NH4_tot(1)/B.NH4(1);
check(isfinite(Vw) && Vw > 0, 'field-water volume is positive');

%% ---- pools the daily solver receives
check(B.Na(1) >= 0, 'porewater Na is non-negative');
check(B.An(1) >= 0, 'counter-anion is non-negative');

% Kantola recompute must close WITHOUT engaging its max(An,0) floor:
% the raw carbonate-frame value equals the stored An.
An_raw = 2*B.Mg(1) + 2*B.Ca(1) + B.Na(1) + B.K(1) + B.NH4(1) ...
    - B.Alk(1) - B.PO4(1) - B.NO3(1);
check(An_raw > -1e-6*max(1, abs(B.An(1))), ...
    'carbonate-frame counter-anion closes without the zero floor');
check(abs(An_raw - B.An(1)) <= 1e-6*max(1, abs(B.An(1))), ...
    'stored An equals the carbonate-frame closure');

%% ---- the alkalinity ledger telescopes to the hostable identity
ledger = B.Alk(1)*Vw + B.R_alk(1);
check(abs(B.Alk_tot(1) - ledger) <= 1e-6*max(1, abs(B.Alk_tot(1))), ...
    'Alk_tot(1) = Alk(1)*Vw + R_alk(1) (ledger telescopes)');
check(B.Alk_tot(1) > 0, 'Alk_tot(1) is positive (hostable by the CEC)');

%% ---- exchange fractions stay a partition of unity
f_sum = B.f_Ca(1) + B.f_Mg(1) + B.f_K(1) + B.f_Na(1) + B.f_Al(1) + B.f_H(1);
check(abs(f_sum - 1) <= 1e-8, 'exchange fractions sum to 1');
check(all([B.f_Ca(1), B.f_Mg(1), B.f_K(1), B.f_Na(1), B.f_Al(1), B.f_H(1)] >= 0), ...
    'exchange fractions are non-negative');
f_base = B.f_Ca(1) + B.f_Mg(1) + B.f_K(1) + B.f_Na(1);
f_base_seed = sum(EC.bulk_pH_params.initial_f_CEC(1:4));
Ca_tot_seed = ERW_biogeochem.Ca(1)*Vw + ...
    EC.bulk_pH_params.initial_f_CEC(1)*B.CEC_tot(1)/2;
Mg_tot_seed = ERW_biogeochem.Mg(1)*Vw + ...
    EC.bulk_pH_params.initial_f_CEC(2)*B.CEC_tot(1)/2;
K_tot_seed = ERW_biogeochem.K(1)*Vw + ...
    EC.bulk_pH_params.initial_f_CEC(3)*B.CEC_tot(1);
fprintf('  f_base_seed=%.12g  f_base_equilibrium=%.12g  I_Ca=%.12g umol/m2/d\n', ...
    f_base_seed, f_base, EC.I_Ca);
fprintf('  pH_bulk=%.12g  pH_pore=%.12g\n', B.pH_bulk(1), B.pH(1));
check(abs(f_base_seed - 0.80) <= 1e-8, ...
    'base-saturation inventory matches the prescribed pretreatment target');
check(abs(B.Ca_tot(1)-Ca_tot_seed) <= 1e-10*max(1,Ca_tot_seed), ...
    'bulk-pH inversion conserves the prescribed dissolved plus exchangeable Ca inventory');
check(abs(B.Mg_tot(1)-Mg_tot_seed) <= 1e-10*max(1,Mg_tot_seed), ...
    'bulk-pH inversion conserves the prescribed dissolved plus exchangeable Mg inventory');
check(abs(B.K_tot(1)-K_tot_seed) <= 1e-10*max(1,K_tot_seed), ...
    'bulk-pH inversion conserves the prescribed dissolved plus exchangeable K inventory');
check(B.Ca(1) <= Ca_tot_seed/Vw && B.Mg(1) <= Mg_tot_seed/Vw && ...
    B.K(1) <= K_tot_seed/Vw, ...
    'equilibrated dissolved cations cannot exceed their conserved field inventories');
check(abs(B.pH_bulk(1)-6.1) <= 1e-8, ...
    'field bulk pH retains the measured pretreatment value');
check(isfinite(B.pH(1)) && B.pH(1) > 0 && B.pH(1) < 14, ...
    'pore-water pH is a finite physical inversion result');
check(abs(B.pH(1)-6.1) > 0.1, ...
    'measured bulk pH is not reused as the pore-water state');
check(isfield(B,'AlOH3_solid') && isfinite(B.AlOH3_solid(1)) && B.AlOH3_solid(1) >= 0, ...
    'finite secondary Al(OH)3 inventory is retained');
diagnostic_fields = {'Al_inventory_residual','charge_residual','gamma_i'};
for k_diag = 1:numel(diagnostic_fields)
    field = diagnostic_fields{k_diag};
    check(isfield(B,field) && any(strcmp(bg_names,field)), ...
        sprintf('%s is retained by the shared 1-D/2-D write path',field));
end

Al_exchange = B.f_Al(1)*B.CEC_tot(1)*EC.Al_CEC_scale/3;
Al_dissolved = B.Al_w(1)*Vw;
Al_residual_expected = Al_exchange+Al_dissolved+B.AlOH3_solid(1)-B.Al_tot(1);
check(abs(B.Al_inventory_residual(1)-Al_residual_expected) <= ...
    1e-12*max(1,abs(B.Al_tot(1))), ...
    'retained Al residual closes aqueous, exchangeable and solid inventory');
charge_residual_expected = B.HCO3(1)+2*B.CO3(1)-B.H(1)+B.k_w(1)/B.H(1)-B.Alk(1);
check(abs(B.charge_residual(1)-charge_residual_expected) <= ...
    1e-12*max(1,abs(B.Alk(1))), ...
    'retained signed charge residual uses the carbonate charge basis');
appelo_test = [];
if isfield(EC.bulk_pH_params,'appelo')
    appelo_test = EC.bulk_pH_params.appelo;
end
[~,~,~,~,~,alpha_test] = Resolve_Appelo_Constants(appelo_test);
gamma_expected = 10.^(alpha_test.*[1 1 2 2 3].*B.f_H(1));
check(all(isfinite(B.gamma_i(1,:))) && ...
    max(abs(B.gamma_i(1,:)-gamma_expected)) <= 1e-12*max(gamma_expected), ...
    'retained Appelo gamma matches explicit alpha, valence and proton fraction');
[K_Na_H,K_K_H,K_Ca_H2,K_Mg_H2,K_Al_H3] = Resolve_Appelo_Constants(appelo_test);
rho_H = B.f_H(1)/B.H(1)*10^(alpha_test*B.f_H(1));
Al3_appelo = B.Al(1)/EC.Al_CEC_scale;
f_appelo_expected = [K_Ca_H2*B.Ca(1)*rho_H^2, ...
    K_Mg_H2*B.Mg(1)*rho_H^2, K_K_H*B.K(1)*rho_H, ...
    K_Na_H*B.Na(1)*rho_H, K_Al_H3*Al3_appelo*rho_H^3];
check(max(abs([B.f_Ca(1),B.f_Mg(1),B.f_K(1),B.f_Na(1),B.f_Al(1)]- ...
    f_appelo_expected)) <= 1e-10, ...
    'initial aqueous and exchange states satisfy Appelo mass action');
check(EC.bulk_pH_params.pH_offset == 0, ...
    'bulk-pH observation operator has no additive baseline offset');
check(strcmp(EC.bulk_pH_params.initial_pH_mode,'bulk_observed'), ...
    'Kantola seed declares the measured bulk-observation boundary');
check(abs(f_base-0.80) <= 1e-2, ...
    'equilibrium base saturation remains within one percentage point of the prescribed target');
check(~isfield(EC.bulk_pH_params.appelo,'initial_affinity_shift'), ...
    'initialization creates no runtime affinity-calibration state');
check(max(abs(EC.bulk_pH_params.appelo.logKH-fixed_logKH)) <= 1e-12, ...
    'bulk-pH inversion leaves fixed Kantola logKH unchanged');
check(EC.I_Ca == 0, ...
    'Kantola seed has no unconstrained boundary Ca supply');

P0 = 1:79;
Pseed = ERW_Initialize_TC_Cation_Pools(P0, B, EC, TT.s(1));
check(abs(sum(Pseed(60:61))-B.Ca_tot(1)*40/1e6) <= 1e-10, ...
    'T&C Ca solution plus exchange equals the ERW seed inventory');
check(abs(sum(Pseed(68:69))-B.Mg_tot(1)*24.3/1e6) <= 1e-10, ...
    'T&C Mg solution plus exchange equals the ERW seed inventory');
check(abs(sum(Pseed(52:53))-B.K_tot(1)*39.1/1e6) <= 1e-10, ...
    'T&C K solution plus exchange equals the ERW seed inventory');
check(isequal(Pseed([1:51,54:59,62:67,70:79]), ...
    P0([1:51,54:59,62:67,70:79])), ...
    'cation synchronization leaves unrelated T&C pools unchanged');

fprintf('  Na(1)=%.6g umol/L  An(1)=%.6g umolc/L  Alk_tot(1)=%.6g umol/m2\n', ...
    B.Na(1), B.An(1), B.Alk_tot(1));
fprintf('  fixed logKH=[%s]\n', ...
    strtrim(sprintf(' %.12g',EC.bulk_pH_params.appelo.logKH)));
fprintf('  ---------------------------------------\n');
fprintf('  %d passed, 0 failed\n', n_pass);

    function check(ok, label)
        assert(ok, 'test_kantola_seed_charge:Failed', '%s', label);
        n_pass = n_pass + 1;
        fprintf('  PASS  %s\n', label);
    end

    function assert_error(f,expected_id)
        try
            f();
        catch ME
            assert(strcmp(ME.identifier,expected_id), ...
                'Expected %s, got %s: %s',expected_id,ME.identifier,ME.message);
            fprintf('  PASS  rejects invalid input with %s\n',expected_id);
            return
        end
        error('Expected error %s was not raised.',expected_id);
    end

end
