% Debug: replicate the seed Invert trial at the actual seed H with the
% actual case-5 workspace values, and dissect the An closure.
clear; clc;
run_file = fullfile('..','..','SIM','EQUIV_1D_TMP','TC_ERW_results_Kantola_crop_control.mat');
run_file = fullfile(pwd, run_file);
S = load(run_file);

H_pw = S.ERW_biogeochem.H(1);            % umol/L
pH_pw = -log10(H_pw/S.conv_mol);
Ca0 = S.ERW_biogeochem.Ca(1); Mg0 = S.ERW_biogeochem.Mg(1);
K0 = S.ERW_biogeochem.K(1); Na0 = S.ERW_biogeochem.Na(1);
NO3 = S.ERW_biogeochem.NO3(1); PO4 = S.ERW_biogeochem.PO4(1); NH4 = S.ERW_biogeochem.NH4(1);
bpp = S.bulk_pH_params_now;
f0 = bpp.initial_f_CEC;
Al_w0 = bpp.initial_Al_w;
CEC_tot = S.ERW_biogeochem.CEC_tot(1);
field_water = bpp.field_water_L_m2;
conv_mol = S.conv_mol;
AlCEC = S.conv_Al;
fprintf('bpp fields: %s\n', strjoin(fieldnames(bpp), ', '));
if isfield(bpp,'K1_Al'); K1=bpp.K1_Al; K2=bpp.K2_Al; K3=bpp.K3_Al; K4=bpp.K4_Al;
else
    % Excel hardcodes pK1=5,pK2=5.1,pK3=6.7,pK4=6.2 (Weil & Brady), K=10^-pK*conv_mol
    K1 = 10^-(5)*conv_mol; K2 = 10^-(5.1)*conv_mol;
    K3 = 10^-(6.7)*conv_mol; K4 = 10^-(6.2)*conv_mol;
    fprintf('K1..K4 from excel hardcoded pK (%.4g %.4g %.4g %.4g umol/L)\n', K1,K2,K3,K4);
end
k1 = S.ERW_biogeochem.k1(1); k2 = S.ERW_biogeochem.k2(1);
kw = S.ERW_biogeochem.k_w(1); knh3 = S.ERW_biogeochem.k_nh3(1);
CO2_w = S.ERW_biogeochem.CO2_w(1);

fprintf('pH_pw=%6.4f  H=%10.5g umol/L\n', pH_pw, H_pw);
fprintf('Ca0=%g Mg0=%g K0=%g Na0=%g NO3=%g PO4=%g NH4=%g\n', Ca0,Mg0,K0,Na0,NO3,PO4,NH4);
fprintf('CEC_tot=%g umolc/m2  field_water=%g L/m2  cec_per_L=%g\n', CEC_tot, field_water, CEC_tot/field_water);
fprintf('f0=[%s]\n', sprintf('%.4f ', f0));
fprintf('Al_w0=%g  AlCEC=%g\n', Al_w0, AlCEC);
fprintf('K_CEC=%s\n', sprintf('%.4g ', bpp.K_CEC));
fprintf('k1=%.5g k2=%.5g kw=%.5g knh3=%.5g CO2_w=%g\n', k1,k2,kw,knh3,CO2_w);

cec_per_L = CEC_tot/field_water;
Ca_total = Ca0 + cec_per_L*f0(1)/2;
Mg_total = Mg0 + cec_per_L*f0(2)/2;
K_total  = K0  + cec_per_L*f0(3);
Na_total = Na0 + cec_per_L*f0(4);
Al_total = Al_w0 + cec_per_L*f0(5)*AlCEC/3;
Al_den = H_pw^4 + H_pw^3*K1 + H_pw^2*K1*K2 + H_pw*K1*K2*K3 + K1*K2*K3*K4;
Al3_frac = H_pw^4/Al_den;
fprintf('totals: Ca=%g Mg=%g K=%g Na=%g Al=%g  Al3_frac=%g\n', ...
    Ca_total, Mg_total, K_total, Na_total, Al_total, Al3_frac);

% replicate exchange_site_residual bisection (same loop as Invert_Bulk_pH)
site_residual = @(fCa) exch_res(fCa, H_pw, Ca_total, Mg_total, K_total, Na_total, Al_total, ...
    cec_per_L, bpp.K_CEC, Al3_frac, AlCEC);
f_lo = 1e-12;
f_hi = min(1-1e-12, 2*Ca_total/cec_per_L - 1e-12);
g_lo = site_residual(f_lo); g_hi = site_residual(f_hi);
fprintf('f_hi=%g  g_lo=%g  g_hi=%g\n', f_hi, g_lo, g_hi);
for i = 1:100
    fm = 0.5*(f_lo+f_hi);
    gm = site_residual(fm);
    if abs(gm) < 1e-12; break; end
    if sign(gm) == sign(g_lo); f_lo = fm; g_lo = gm; else; f_hi = fm; end
end
f_Ca = 0.5*(f_lo+f_hi);
[R, f_Mg, f_K, f_Na, f_Al, f_H] = site_residual(f_Ca);
Ca = Ca_total - cec_per_L*f_Ca/2;
Mg = Mg_total - cec_per_L*f_Mg/2;
K  = K_total  - cec_per_L*f_K;
Na = Na_total - cec_per_L*f_Na;
Al_w = Al_total - cec_per_L*f_Al*AlCEC/3;
Al = Al3_frac*Al_w;
fprintf('partition: f_Ca=%g f_Mg=%g f_K=%g f_Na=%g f_Al=%g f_H=%g\n', f_Ca,f_Mg,f_K,f_Na,f_Al,f_H);
fprintf('dissolved: Ca=%g Mg=%g K=%g Na=%g Al_w=%g Al=%g\n', Ca,Mg,K,Na,Al_w,Al);

% An closure components (ERW_Residual_An_From_Charge internals)
Al_w_aq = Al_w/AlCEC;
OH = kw/H_pw; HCO3 = k1*CO2_w/H_pw; CO3 = k2*k1*CO2_w/H_pw^2;
Al3 = (H_pw^4/Al_den)*Al_w_aq;
AlOH = (H_pw^3*K1/Al_den)*Al_w_aq;
AlOH2 = (H_pw^2*K1*K2/Al_den)*Al_w_aq;
AlOH3 = (H_pw*K1*K2*K3/Al_den)*Al_w_aq;
AlOH4 = max(Al_w_aq-(Al3+AlOH+AlOH2+AlOH3),0);
Al_net = 3*Al3 + 2*AlOH + AlOH2 - AlOH4;
An = H_pw + 2*Ca + 2*Mg + Na + K + NH4 + Al_net - (OH + HCO3 + 2*CO3 + NO3 + PO4);
fprintf('\nAn closure:\n');
fprintf('  H=%g 2Ca=%g 2Mg=%g Na=%g K=%g NH4=%g\n', H_pw, 2*Ca, 2*Mg, Na, K, NH4);
fprintf('  Al3=%g AlOH=%g AlOH2=%g AlOH3=%g AlOH4=%g Al_net=%g\n', Al3,AlOH,AlOH2,AlOH3,AlOH4,Al_net);
fprintf('  OH=%g HCO3=%g 2CO3=%g NO3=%g PO4=%g\n', OH,HCO3,2*CO3,NO3,PO4);
fprintf('  => An = %g  (expected ~%.5g)\n', An, S.ERW_biogeochem.An(1));

function [residual, f_Mg, f_K, f_Na, f_Al, f_H] = exch_res(f_Ca, H, Ca_total, Mg_total, ...
    K_total, Na_total, Al_total, cec_per_L, K_CEC, Al3_frac, AlCEC)
Ca = Ca_total - cec_per_L*f_Ca/2;
if ~(isfinite(Ca) && Ca > 0)
    residual = NaN; f_Mg = NaN; f_K = NaN; f_Na = NaN; f_Al = NaN; f_H = NaN;
    return;
end
KC = K_CEC;
A_Mg = f_Ca/(KC(1)*Ca);
A_K  = sqrt(f_Ca/(KC(2)*Ca));
A_Na = sqrt(f_Ca/(KC(3)*Ca));
A_Al = (Al3_frac/AlCEC)*sqrt(f_Ca^3/(KC(4)*Ca^3));
f_Mg = A_Mg*Mg_total/(1 + A_Mg*cec_per_L/2);
f_K  = A_K *K_total /(1 + A_K *cec_per_L);
f_Na = A_Na*Na_total/(1 + A_Na*cec_per_L);
f_Al = A_Al*Al_total/(1 + A_Al*cec_per_L*AlCEC/3);
f_H  = H*sqrt(f_Ca/(KC(5)*Ca));
residual = f_Ca + f_Mg + f_K + f_Na + f_Al + f_H - 1;
end