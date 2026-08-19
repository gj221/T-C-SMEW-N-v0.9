function test_pH_conversion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Bulk-extract pH pipeline in v2 (port of the Model TC_Code machinery).
%%
%% The observed air-dried 1:n extract pH is anchored at simulation start
%% (Invert_Bulk_pH in ERW_MOD_PARAM_excel) and the solved pore-water state
%% is converted back to the bulk extract pH every day
%% (Porewater_to_Bulk_pH in ERW_Daily_Step_2D and MAIN_FRAME_BG_ERW_10mins).
%% Daily and extract solvers share ONE exchange scheme selected by
%% bulk_pH_params.scheme: 'appelo' (SCEPTER/Appelo H-referenced, default)
%% or 'TC_GT' (Gaines-Thomas, Model constants).
%% (SCEPTER/Appelo H-referenced selectivity).
%%
%% DELIBERATE DEVIATION FROM THE MODEL: v2 hands the bulk extract pH to the
%% TC reaction block (E.PH / PH), keeping the pore-water pH in E.PH_pore /
%% ERW_biogeochem.pH for diagnostics.  The Model keeps wall chemistry on
%% pore-water pH and exposes bulk pH only to the observation output.
%% v2 chooses bulk pH so the 2-D driver's PH_d (which feeds BIOGEO_UNIT)
%% matches what the field measured, and keeps PH_d comparable between the
%% 1-D frame and the 2-D driver in test_2d_driver_equivalence.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n=== test_pH_conversion ===\n');
n_pass = 0;

v2 = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(v2, fullfile(v2, 'T&C_Code'));

conv_mol = 1e6;
Al_CEC_scale = 1e3;
k1 = 0.447;
k2 = 4.68e-5;
k_w = 1e-2;
k_nh3 = 5.623413e-4;
k_H = 0.83;
K1_Al = 10; K2_Al = 8; K3_Al = 0.2; K4_Al = 0.631;

%% ---- TC_GT forward/inverse round trip on a prescribed exchange state
f_CEC = [0.647 0.129 0.043 0 0 0.181];
K_CEC = [10^0.1, 10^(-2)*conv_mol, 10^0.38*conv_mol, ...
    10^(-0.86)/conv_mol, 10^(-5.4)*conv_mol];
H_expected = 10^(-5.2) * conv_mol;
Ca = (f_CEC(1)/K_CEC(5)) * (H_expected/f_CEC(6))^2;
Mg = (Ca/f_CEC(1)) * K_CEC(1) * f_CEC(2);
K  = sqrt((Ca/f_CEC(1)) * K_CEC(2)) * f_CEC(3);
Na = 0;
CO2_air = 400/22.41;
CO2_w = 10 * k_H * CO2_air;
[An, ~] = ERW_Residual_An_From_Charge( ...
    H_expected, CO2_w, Ca, Mg, Na, K, 0, 0, 0, 0, ...
    k1, k2, k_w, k_nh3, K1_Al, K2_Al, K3_Al, K4_Al, Al_CEC_scale);
check(An >= 0, 'residual anion from charge closure is non-negative');

CEC_tot = 25 * 1e-5 * 1.15e6 * 0.2 * conv_mol;
bpp = struct('n_extract', 2.5, 'K_CEC', K_CEC, ...
    'H_extract_anchor_umolL', H_expected, ...
    'H_pore_anchor_umolL', H_expected, ...
    'H_prev', NaN, 'H_prev_valid', false, 'initial_f_CEC', f_CEC, ...
    'field_water_L_m2', 25);
pH_bulk = Porewater_to_Bulk_pH( ...
    H_expected, Ca, Mg, K, Na, 0, 0, An, 0, 0, 0, CO2_air, ...
    k1, k2, k_w, k_nh3, k_H, K1_Al, K2_Al, K3_Al, K4_Al, ...
    f_CEC(1), f_CEC(2), f_CEC(3), f_CEC(4), f_CEC(5), f_CEC(6), ...
    CEC_tot, 1.15e6, 0.2, conv_mol, bpp, Al_CEC_scale);
check(isfinite(pH_bulk), 'TC_GT forward extraction returns a finite bulk pH');

[H, ~, f_Ca, f_Mg, f_K, f_Na, f_Al, f_H, An_eq, ~, Ca_eq, Mg_eq, K_eq, Na_eq] = ...
    Invert_Bulk_pH(pH_bulk, Ca, Mg, K, Na, 0, 0, 0, CO2_air, CO2_w, ...
    k1, k2, k_w, k_nh3, k_H, K1_Al, K2_Al, K3_Al, K4_Al, ...
    CEC_tot, 1.15e6, 0.2, conv_mol, bpp, 9.66, Al_CEC_scale);
check(abs(-log10(H/conv_mol) - 5.2) < 2e-3, ...
    'TC_GT inverse recovers the pore pH to +/-0.002');
check(norm([f_Ca f_Mg f_K f_Na f_Al f_H] - f_CEC, inf) < 2e-5, ...
    'TC_GT inverse reproduces the prescribed exchange fractions');
check(norm([Ca_eq Mg_eq K_eq Na_eq] - [Ca Mg K Na], inf) < ...
    1e-2*max(norm([Ca Mg K Na], inf), 1) + 1e-9, ...
    'TC_GT inverse leaves the dissolved cation inventory unchanged');
check(An_eq >= 0, 'TC_GT inverse returns a non-negative residual anion');

%% ---- appelo forward/inverse round trip on the same state
[f_Ca_a, f_Mg_a, f_K_a, f_Na_a, f_Al_a, f_H_a] = Exchange_Fractions_GT( ...
    H_expected, Ca, Mg, K, Na, 0, K_CEC, ...
    K1_Al, K2_Al, K3_Al, K4_Al, Al_CEC_scale, 'appelo');
check(abs(sum([f_Ca_a f_Mg_a f_K_a f_Na_a f_Al_a f_H_a]) - 1) < 1e-10, ...
    'appelo exchange fractions sum to unity');
bpp_app = bpp;
bpp_app.scheme = 'appelo';
pH_bulk_app = Porewater_to_Bulk_pH( ...
    H_expected, Ca, Mg, K, Na, 0, 0, An, 0, 0, 0, CO2_air, ...
    k1, k2, k_w, k_nh3, k_H, K1_Al, K2_Al, K3_Al, K4_Al, ...
    f_Ca_a, f_Mg_a, f_K_a, f_Na_a, f_Al_a, f_H_a, ...
    CEC_tot, 1.15e6, 0.2, conv_mol, bpp_app, Al_CEC_scale);
check(isfinite(pH_bulk_app), ...
    'appelo forward extraction returns a finite bulk pH');

%% The inverse repartitions the dissolved and exchangeable inventories at
%% every trial H, so the state returned by the inversion is the fixed
%% point of ITS operator.  The invariant to hold is therefore that the
%% inverse's own outputs forward-extract back to the observed bulk pH
%% (the semantic guarantee the daily driver relies on), plus a loose
%% anchor proximity on the pore pH (H_pore_anchor_umolL = 5.2 here).
[H_app, ~, f_Ca_r, f_Mg_r, f_K_r, f_Na_r, f_Al_r, f_H_r, An_app, ...
    Al_w_app, Ca_eq, Mg_eq, K_eq, Na_eq] = Invert_Bulk_pH( ...
    pH_bulk_app, Ca, Mg, K, Na, 0, 0, 0, CO2_air, CO2_w, ...
    k1, k2, k_w, k_nh3, k_H, K1_Al, K2_Al, K3_Al, K4_Al, ...
    CEC_tot, 1.15e6, 0.2, conv_mol, bpp_app, 9.66, Al_CEC_scale);
check(abs(-log10(H_app/conv_mol) - 5.2) < 0.1, ...
    'appelo inverse stays within 0.1 of the pore pH anchor');
check(An_app >= 0, 'appelo inverse returns a non-negative residual anion');
pH_bulk_rt = Porewater_to_Bulk_pH( ...
    H_app, Ca_eq, Mg_eq, K_eq, Na_eq, Al_w_app, 0, An_app, 0, 0, 0, CO2_air, ...
    k1, k2, k_w, k_nh3, k_H, K1_Al, K2_Al, K3_Al, K4_Al, ...
    f_Ca_r, f_Mg_r, f_K_r, f_Na_r, f_Al_r, f_H_r, ...
    CEC_tot, 1.15e6, 0.2, conv_mol, bpp_app, Al_CEC_scale);
check(abs(pH_bulk_rt - pH_bulk_app) < 2e-3, ...
    'appelo inverse outputs forward-extract back to the observed bulk pH');

%% ---- appelo exchange fractions satisfy the analytic relations
[H_a2, Ca_a2, Mg_a2, K_a2, Na_a2, Al_w_a2] = deal(1, 300, 100, 60, 50, 1000);
[f_Ca, f_Mg, f_K, f_Na, f_Al, f_H] = Exchange_Fractions_GT( ...
    H_a2, Ca_a2, Mg_a2, K_a2, Na_a2, Al_w_a2, ones(1,5), ...
    K1_Al, K2_Al, K3_Al, K4_Al, Al_CEC_scale, 'appelo');
[K_Na_H, K_K_H, K_Ca_H2, K_Mg_H2, K_Al_H3, alpha] = Resolve_Appelo_Constants([]);
den_Al = H_a2^4 + H_a2^3*K1_Al + H_a2^2*K1_Al*K2_Al + ...
    H_a2*K1_Al*K2_Al*K3_Al + K1_Al*K2_Al*K3_Al*K4_Al;
Al3 = H_a2^4/den_Al * Al_w_a2/Al_CEC_scale;
gamma = 10^(alpha*f_H);
rho = f_H/H_a2 * gamma;
expected = [K_Ca_H2*Ca_a2*rho^2, K_Mg_H2*Mg_a2*rho^2, ...
    K_K_H*K_a2*rho, K_Na_H*Na_a2*rho, K_Al_H3*Al3*rho^3, f_H];
check(norm([f_Ca f_Mg f_K f_Na f_Al f_H] - expected, inf) < 1e-10, ...
    'appelo fractions match f = K*[X]*(f_H/[H])^z*gamma^z analytically');

%% ---- invalid scheme is rejected
bpp_bad = bpp;
bpp_bad.scheme = 'bogus';
assert_error(@() Porewater_to_Bulk_pH( ...
    H_expected, Ca, Mg, K, Na, 0, 0, An, 0, 0, 0, CO2_air, ...
    k1, k2, k_w, k_nh3, k_H, K1_Al, K2_Al, K3_Al, K4_Al, ...
    f_CEC(1), f_CEC(2), f_CEC(3), f_CEC(4), f_CEC(5), f_CEC(6), ...
    CEC_tot, 1.15e6, 0.2, conv_mol, bpp_bad, Al_CEC_scale), 'ERW:invalidParam');
n_pass = n_pass + 1;

%% ---- plumbing: drivers hand bulk pH onward, keep pore pH for diagnostics
daily_src = fileread(fullfile(v2, 'T&C_Code', 'ERW_Daily_Step_2D.m'));
check(~isempty(strfind(daily_src, 'E.PH = pH_bulk_now;')), ...
    '2-D daily step exposes E.PH as the bulk extract pH');
check(~isempty(strfind(daily_src, 'E.PH_pore = bg.pH(dm_end,:);')), ...
    '2-D daily step keeps the pore-water pH in E.PH_pore');
frame_src = fileread(fullfile(v2, 'T&C_Code', 'MAIN_FRAME_BG_ERW_10mins.m'));
check(~isempty(strfind(frame_src, 'PH = pH_bulk_now;')), ...
    '1-D frame TC handoff is the bulk extract pH');
check(~isempty(strfind(frame_src, 'PH_pore_now = ERW_biogeochem.pH(dm_end,:);')), ...
    '1-D frame keeps the pore-water pH for diagnostics');
excel_src = fileread(fullfile(v2, 'T&C_Code', 'Parameters', 'ERW_MOD_PARAM_excel.m'));
check(~isempty(strfind(excel_src, 'Invert_Bulk_pH(')), ...
    'seed pH anchoring inverts the observed bulk pH');
erw_src = fileread(fullfile(v2, 'T&C_Code', 'ERW_function.m'));
check(~isempty(strfind(erw_src, 'scheme_is_appelo')), ...
    'daily solver selects the exchange scheme');
check(~isempty(strfind(erw_src, 'exchange_relations(p)')), ...
    'daily solver uses the scheme-selected exchange relations');

%% ---- the ERW_MOD workbook carries the extraction ratio column
wb = fullfile(v2, 'T&C_Code', 'Inputs', 'ERW_MOD_PARAM_parameters.xlsx');
check(exist(wb, 'file') == 2, 'workbook exists beside the parameter code');
opts = detectImportOptions(wb, 'Sheet', 'ERW_MOD');
check(any(strcmp(opts.VariableNames, 'pH_DilutionFactor')), ...
    'ERW_MOD sheet has the pH_DilutionFactor column');

fprintf('  ---------------------------------------\n');
fprintf('  %d passed, 0 failed\n', n_pass);

    function check(ok, label)
        assert(ok, 'test_pH_conversion:Failed', '%s', label);
        n_pass = n_pass + 1;
        fprintf('  PASS  %s\n', label);
    end

end

function assert_error(f, expected_id)
try
    f();
catch ME
    assert(strcmp(ME.identifier, expected_id), ...
        'Expected %s, got %s: %s', expected_id, ME.identifier, ME.message);
    fprintf('  PASS  rejects invalid input with %s\n', expected_id);
    return
end
error('Expected error %s was not raised.', expected_id);
end