% Debug: Gaines-Thomas vs Appelo on the real case-5 (Kantola crop control)
% constants. Compares (1) workbook input sanity, (2) the seed fixed point
% under each scheme via Invert_Bulk_pH, (3) the forward extract map.
clear; clc;
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..')); cd(here);

run_file = fullfile(pwd,'..','..','SIM','EQUIV_1D_TMP','TC_ERW_results_Kantola_crop_control.mat');
S = load(run_file);
bpp0 = S.bulk_pH_params_now;
conv_mol = S.conv_mol; AlCEC = S.conv_Al;
fprintf('ERW_const fields: %s\n', strjoin(fieldnames(S.ERW_const), ', '));
fprintf('TT_m vars: %s\n', strjoin(S.TT_m.Properties.VariableNames, ', '));

% ---- 1. workbook case-5 row (row 4: Kantola crop control) ----
xl = fullfile(pwd,'..','..','T&C_Code','Inputs','ERW_MOD_PARAM_parameters.xlsx');
opts = detectImportOptions(xl, 'Sheet', 'ERW_MOD');
T = readtable(xl, opts);
hdr = opts.VariableNames;
fprintf('\nworkbook row 4 values:\n');
want = {'Na_in','Ca_in','Mg_in','K_in','Si_in','Al_w_in','CEC_tot_in','f_CEC_in','pH_DilutionFactor','soil'};
for w = 1:numel(want)
    if any(strcmp(hdr, want{w}))
        v = T.(want{w})(4);
        if iscell(v); v = v{1}; end
        fprintf('  %-18s = %s\n', want{w}, mat2str(double(v)));
    end
end

% ---- 2. constants from the run workspace ----
k1 = S.ERW_biogeochem.k1(1); k2 = S.ERW_biogeochem.k2(1);
kw = S.ERW_biogeochem.k_w(1); knh3 = S.ERW_biogeochem.k_nh3(1);
kH = S.ERW_biogeochem.k_H(1);
CO2_w = S.ERW_biogeochem.CO2_w(1);
CO2_atm = S.ERW_biogeochem.CO2_air(1);
K1 = S.ERW_const.K1; K2 = S.ERW_const.K2; K3 = S.ERW_const.K3; K4 = S.ERW_const.K4;
K_AlOH3 = 9.66;
rho_bulk = S.ERW_const.rho_bulk;
Zrw = S.ERW_const.Zrw;
CEC_tot = S.ERW_biogeochem.CEC_tot(1);
f0 = bpp0.initial_f_CEC(:).';
Al_w0 = bpp0.initial_Al_w;
pH_obs = bpp0.observed_pH_initial;
fprintf('\nconstants: rho_bulk=%g Zrw=%g CEC_tot=%g f0=[%s] Al_w0=%g pH_obs=%g\n', ...
    rho_bulk, Zrw, CEC_tot, sprintf('%.4f ', f0), Al_w0, pH_obs);

% day-1 dissolved state (from the saved run, where inputs come from the
% spinup state; the workbook *_in columns are informational/empty)
Ca0 = S.ERW_biogeochem.Ca(1); Mg0 = S.ERW_biogeochem.Mg(1);
K0 = S.ERW_biogeochem.K(1); Na0 = S.ERW_biogeochem.Na(1);
NO3 = S.ERW_biogeochem.NO3(1); PO4 = S.ERW_biogeochem.PO4(1); NH4 = S.ERW_biogeochem.NH4(1);
legacy_An = 2*Mg0 + 2*Ca0 + Na0 + K0 + NH4 - (k1*CO2_w/(10^(-pH_obs)*conv_mol)) - PO4 - NO3;
fprintf('day-1 dissolved state: Ca=%g Mg=%g K=%g Na=%g NO3=%g PO4=%g NH4=%g\n', ...
    Ca0, Mg0, K0, Na0, NO3, PO4, NH4);
fprintf('  legacy_An (no Al, C=6.1): %g     An(1) observed: %g\n', ...
    legacy_An, S.ERW_biogeochem.An(1));
fprintf('  Al_w(1)=%g nmol/L, Al(1)=%g, pH(1)=%.4f, pH_bulk(1)=%.4f\n', ...
    S.ERW_biogeochem.Al_w(1), S.ERW_biogeochem.Al(1), ...
    S.ERW_biogeochem.pH(1), S.ERW_biogeochem.pH_bulk(1));

% ---- 3. seed fixed point under each scheme ----
bpp_gt = struct('K_CEC', bpp0.K_CEC, 'CEC_cmolc_kg', bpp0.CEC_cmolc_kg, ...
    'n_extract', bpp0.n_extract, 'field_water_L_m2', bpp0.field_water_L_m2, ...
    'initial_f_CEC', f0, 'initial_Al_w', Al_w0);
bpp_ap = bpp_gt;
bpp_ap.scheme = 'appelo';
bpp_ap.appelo = struct('logKH', [5.9 4.8 10.47 10.786 16.47], 'alpha', 3.4);

for sc = 1:2
    if sc == 1; bpp = bpp_gt; tag = 'TC_GT '; else; bpp = bpp_ap; tag = 'APPELO'; end
    try
        [Hp, ~, fCa, fMg, fK, fNa, fAl, fH, An, Alw, Ca, Mg, K, Na] = Invert_Bulk_pH( ...
            pH_obs, Ca0, Mg0, K0, Na0, NO3, PO4, NH4, ...
            CO2_atm, CO2_w, k1, k2, kw, knh3, kH, ...
            K1, K2, K3, K4, CEC_tot, rho_bulk, Zrw, conv_mol, ...
            bpp, K_AlOH3, AlCEC);
        fprintf('\n[%s] seed fixed point: pore pH=%.4f  f_CEC=[%s]\n', tag, ...
            -log10(Hp/conv_mol), sprintf('%.4f ', [fCa fMg fK fNa fAl fH]));
        fprintf('      An=%g  Al_w=%g  Ca=%g Mg=%g K=%g Na=%g\n', An, Alw, Ca, Mg, K, Na);
    catch ME
        fprintf('\n[%s] Invert FAILED: %s\n', tag, ME.message);
    end
end

% ---- 4. forward extract map over pore pH (fixed workbook state+fractions) ----
Al_aq = Al_w0/AlCEC;
Al3_den = (10^(-pH_obs)*conv_mol)^4 + (10^(-pH_obs)*conv_mol)^3*K1 + ...
    (10^(-pH_obs)*conv_mol)^2*K1*K2 + (10^(-pH_obs)*conv_mol)*K1*K2*K3 + K1*K2*K3*K4;
Al3_aq = (10^(-pH_obs)*conv_mol)^4/Al3_den*Al_aq;
for sc = 1:2
    if sc == 1; bpp = bpp_gt; tag = 'TC_GT '; else; bpp = bpp_ap; tag = 'APPELO'; end
    fprintf('\n[%s] forward map: pore pH -> bulk pH\n', tag);
    for pHpw = 2:0.5:12
        try
            pb = Porewater_to_Bulk_pH(10^(-pHpw)*conv_mol, Ca0, Mg0, K0, Na0, Al_aq, Al3_aq, ...
                legacy_An, NO3, PO4, NH4, CO2_atm, k1, k2, kw, knh3, kH, ...
                K1, K2, K3, K4, f0(1), f0(2), f0(3), f0(4), f0(5), f0(6), ...
                CEC_tot, rho_bulk, Zrw, conv_mol, bpp, AlCEC);
            fprintf('   %5.2f -> %6.3f\n', pHpw, pb);
        catch ME
            fprintf('   %5.2f -> ERR: %s\n', pHpw, ME.identifier);
        end
    end
end