% Debug: calibrate the Appelo/SCEPTER cation-exchange constants and alpha
% on the real case-5 (Kantola crop control) data, following the SCEPTER
% refit recipe (spinup_inrt2.py / Kanzaki et al. 2024 GMD, Eq. 3-5):
%
%   logKH_i = z_i * ( -log10( f_i^(1/z_i) * H_pw / (fH_base * c_i^(1/z_i)) )
%                      + alpha * fH_base )
%
% with measured exchangeable fractions f_i (fraction of CEC), porewater
% free concentrations c_i [M], porewater H+ [M], assumed H-site fraction
% fH_base, and surface-charge exponent alpha. Then re-runs the seed
% inversion with (a) each calibrated logK set and (b) aqueous Al pinned by
% amorphous Al(OH)3 solubility (SCEPTER: logK = 10.8, [Al3+] = 10^(10.8-3pH))
% instead of the runaway Al_w charge closure.
clear; clc;
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..')); cd(here);

run_file = fullfile(pwd,'..','..','SIM','EQUIV_1D_TMP','TC_ERW_results_Kantola_crop_control.mat');
S = load(run_file);
bpp0 = S.bulk_pH_params_now;
conv_mol = S.conv_mol; AlCEC = S.conv_Al;

k1 = S.ERW_biogeochem.k1(1); k2 = S.ERW_biogeochem.k2(1);
kw = S.ERW_biogeochem.k_w(1); knh3 = S.ERW_biogeochem.k_nh3(1);
kH = S.ERW_biogeochem.k_H(1);
CO2_w = S.ERW_biogeochem.CO2_w(1); CO2_atm = S.ERW_biogeochem.CO2_air(1);
K1 = S.ERW_const.K1; K2 = S.ERW_const.K2; K3 = S.ERW_const.K3; K4 = S.ERW_const.K4;
rho_bulk = S.ERW_const.rho_bulk; Zrw = S.ERW_const.Zrw;
CEC_tot = S.ERW_biogeochem.CEC_tot(1);
pH_obs = bpp0.observed_pH_initial;
K_AlOH3_sc = 10.8;  % amorphous Al(OH)3, PHREEQC.DAT (SCEPTER)
% K_AlOH3_gb = 8.11;  % gibbsite (alternative)

% ---- 1. measured inputs (case 5) ----
f_meas = bpp0.initial_f_CEC(:).';   % [Ca Mg K Na Al H] fractions of CEC
% porewater concentrations consistent with the exchange law at the seed
% (appelo scheme, from debug_scheme_compare):  [umol/L]
c_pw = [179.651 1319.4 2696.52 541925];   % Ca Mg K Na
c_pw_M = c_pw * 1e-6;                     % [M]
pH_pw = 5.8722;
H_pw = 10^(-pH_pw);                       % [M]
fH_base = 0.22;                           % SCEPTER refit assumption fH ~ 0.25
z = [2 2 1 1];
fprintf('input: f0=[%s] c_pw(M)=[%s] pH_pw=%.4f\n', ...
    sprintf('%.4f ', f_meas), sprintf('%.3e ', c_pw_M), pH_pw);

% ---- 2. SCEPTER inversion recipe over (alpha, fH_base) ----
fprintf('\nlogKH grid (Na K Ca Mg) vs (alpha, fH_base);  SCEPTER defaults 5.9 4.8 10.47 10.786\n');
alpha_list = [0 1.3 2.0 3.4 5.0];
fH_list   = [0.1 0.15 0.22 0.25];
for fh = fH_list
    fprintf(' fH_base=%.2f:', fh);
    for al = alpha_list
        logKH = z .* ( -log10( f_meas(1:4).^(1./z) .* H_pw ./ (fh * c_pw_M.^(1./z)) ) + al * fh );
        fprintf('  a=%.1f[%s]', al, sprintf('%.1f ', logKH));
    end
    fprintf('\n');
end

% ---- 3. seed inversion: SCEPTER-style solubility-pinned Al instead of ----
%        the workbook Al_w (which is empty -> NaN -> runaway closure)    ----
% [Al3+] = 10^(10.8 - 3*pH) at the bulk anchor pH; convert to model Al_w
% units (nmol/L) via the speciation fraction (model umol/L basis).
for fh = fH_list(end)  % fH_base for the chosen calibrated set
    for al = [1.3 3.4]
        logKH = z .* ( -log10( f_meas(1:4).^(1./z) .* H_pw ./ (fh * c_pw_M.^(1./z)) ) + al * fh );
        bpp = struct('K_CEC', bpp0.K_CEC, 'CEC_cmolc_kg', bpp0.CEC_cmolc_kg, ...
            'n_extract', bpp0.n_extract, 'field_water_L_m2', bpp0.field_water_L_m2, ...
            'initial_f_CEC', f_meas, 'scheme', 'appelo', ...
            'appelo', struct('logKH', [logKH(4) logKH(3) logKH(1) logKH(2) 16.47], 'alpha', al));
        H_umol = 10^(-pH_obs) * conv_mol;   % model H basis [umol/L]
        Al3_free_umolL = 10^(K_AlOH3_sc - 3*pH_obs) * conv_mol;
        den = H_umol^4 + H_umol^3*K1 + H_umol^2*K1*K2 + H_umol*K1*K2*K3 + K1*K2*K3*K4;
        fr_al3 = H_umol^4/den;
        Al_w_sol = Al3_free_umolL / fr_al3;      % total Al so free Al3+ = solubility
        bpp.initial_Al_w = Al_w_sol * AlCEC;     % to model nanomol/L units
        fprintf('\n[alpha=%.1f] solubility-pinned Al_w0 = %.4e nmol/L (was 2.16e8)\n', al, bpp.initial_Al_w);
        try
            [Hp, ~, fCa, fMg, fK, fNa, fAl, fH, An, Alw, Ca, Mg, K, Na] = Invert_Bulk_pH( ...
                pH_obs, c_pw(1), c_pw(2), c_pw(3), c_pw(4), ...
                S.ERW_biogeochem.NO3(1), S.ERW_biogeochem.PO4(1), S.ERW_biogeochem.NH4(1), ...
                CO2_atm, CO2_w, k1, k2, kw, knh3, kH, ...
                K1, K2, K3, K4, CEC_tot, rho_bulk, Zrw, conv_mol, ...
                bpp, 9.66, AlCEC);
            fprintf('   seed: pore pH=%.4f f=[%s] An=%g Al_w=%.3e Ca=%g Mg=%g K=%g Na=%g\n', ...
                -log10(Hp/conv_mol), sprintf('%.4f ', [fCa fMg fK fNa fAl fH]), ...
                An, Alw, Ca, Mg, K, Na);
            Al_aq = Alw/AlCEC;
            den2 = Hp^4 + Hp^3*K1 + Hp^2*K1*K2 + Hp*K1*K2*K3 + K1*K2*K3*K4;
            Al3_aq = Hp^4/den2*Al_aq;
            pb = Porewater_to_Bulk_pH(Hp, Ca, Mg, K, Na, Al_aq, Al3_aq, ...
                An, S.ERW_biogeochem.NO3(1), S.ERW_biogeochem.PO4(1), S.ERW_biogeochem.NH4(1), ...
                CO2_atm, k1, k2, kw, knh3, kH, ...
                K1, K2, K3, K4, fCa, fMg, fK, fNa, fAl, fH, ...
                CEC_tot, rho_bulk, Zrw, conv_mol, bpp, AlCEC);
            fprintf('   forward closure: bulk pH = %.4f (target 6.1, d=%.2e)\n', pb, abs(pb-pH_obs));
        catch ME
            fprintf('   Invert FAILED: %s\n', ME.message);
        end
    end
end