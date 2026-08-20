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
%%      sodium every site the measured Ca/Mg/K did not claim (f_Na = 0.76,
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
% Production pin configuration (MOD_PARAM_Kantola_Maize_Soy lines 5-16).
ERW_TC_const.bulk_pH_params = struct('scheme', 'appelo', ...
    'Al_solubility_pin', true);

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

[B, ~, ~, EC] = ERW_MOD_PARAM_excel( ...
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

fprintf('  Na(1)=%.6g umol/L  An(1)=%.6g umolc/L  Alk_tot(1)=%.6g umol/m2\n', ...
    B.Na(1), B.An(1), B.Alk_tot(1));
fprintf('  ---------------------------------------\n');
fprintf('  %d passed, 0 failed\n', n_pass);

    function check(ok, label)
        assert(ok, 'test_kantola_seed_charge:Failed', '%s', label);
        n_pass = n_pass + 1;
        fprintf('  PASS  %s\n', label);
    end

end
