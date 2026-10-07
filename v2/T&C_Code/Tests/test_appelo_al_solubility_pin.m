function test_appelo_al_solubility_pin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Appelo exchange must use free Al3+ directly when amorphous Al(OH)3
%% fixes its aqueous activity.  The solubility pin is not a finite total-Al
%% inventory and therefore must not introduce an inventory denominator.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n=== test_appelo_al_solubility_pin ===\n');

v2 = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(v2, fullfile(v2,'T&C_Code'));

conv_mol = 1e6;
Al_CEC_scale = 1e3;
k1 = 0.447; k2 = 4.68e-5; k_w = 1e-2;
k_nh3 = 5.623413e-4; k_H = 0.83;
K1_Al = 10; K2_Al = 8; K3_Al = 0.2; K4_Al = 0.631;
K_CEC = [10^0.1,10^(-2)*conv_mol,10^0.38*conv_mol, ...
    10^(-0.86)/conv_mol,10^(-5.4)*conv_mol];

H0 = 10^(-5.2)*conv_mol;
Ca = 198; Mg = 50; K = 75; Na = 30;
Al3_pin0 = 10^(-1.2)*H0^3;
den0 = H0^4+H0^3*K1_Al+H0^2*K1_Al*K2_Al+ ...
    H0*K1_Al*K2_Al*K3_Al+K1_Al*K2_Al*K3_Al*K4_Al;
Al_w0 = Al3_pin0/(H0^4/den0)*Al_CEC_scale;
[fCa,fMg,fK,fNa,fAl,fH] = Exchange_Fractions_GT( ...
    H0,Ca,Mg,K,Na,Al_w0,K_CEC,K1_Al,K2_Al,K3_Al,K4_Al, ...
    Al_CEC_scale,'appelo');

CO2_air = 400/22.41;
CO2_w = 10*k_H*CO2_air;
An = ERW_Residual_An_From_Charge(H0,CO2_w,Ca,Mg,Na,K,0,0,0,Al_w0, ...
    k1,k2,k_w,k_nh3,K1_Al,K2_Al,K3_Al,K4_Al,Al_CEC_scale);
CEC_tot = 25*1e-5*1.15e6*0.2*conv_mol;
bpp = struct('n_extract',2.5,'K_CEC',K_CEC,'scheme','appelo', ...
    'Al_solubility_pin',true,'H_extract_anchor_umolL',H0, ...
    'H_pore_anchor_umolL',H0,'H_prev',NaN,'H_prev_valid',false, ...
    'initial_f_CEC',[fCa fMg fK fNa fAl fH], ...
    'initial_Al_w',Al_w0,'field_water_L_m2',25, ...
    'DIC_field_umol_L',CO2_w*(1+k1/H0+k1*k2/H0^2));

pH_bulk = Porewater_to_Bulk_pH(H0,Ca,Mg,K,Na,Al_w0, ...
    Al3_pin0*Al_CEC_scale,An,0,0,0,CO2_air,k1,k2,k_w,k_nh3,k_H, ...
    K1_Al,K2_Al,K3_Al,K4_Al,fCa,fMg,fK,fNa,fAl,fH,CEC_tot, ...
    1.15e6,0.2,conv_mol,bpp,Al_CEC_scale);
[H,~,~,~,~,~,fAl_eq,fH_eq,~,Al_w_eq] = Invert_Bulk_pH( ...
    pH_bulk,Ca,Mg,K,Na,0,0,0,CO2_air,CO2_w,k1,k2,k_w,k_nh3,k_H, ...
    K1_Al,K2_Al,K3_Al,K4_Al,CEC_tot,1.15e6,0.2,conv_mol,bpp, ...
    10.8,Al_CEC_scale);

den = H^4+H^3*K1_Al+H^2*K1_Al*K2_Al+H*K1_Al*K2_Al*K3_Al+ ...
    K1_Al*K2_Al*K3_Al*K4_Al;
Al3 = (H^4/den)*(Al_w_eq/Al_CEC_scale);
Al3_expected = 10^(-1.2)*H^3;
[~,~,~,~,K_Al_H3,alpha] = Resolve_Appelo_Constants([]);
fAl_expected = K_Al_H3*Al3*(fH_eq/H)^3*10^(3*alpha*fH_eq);

assert(abs(Al3-Al3_expected) <= 1e-10*max(1,Al3_expected), ...
    'Amorphous Al(OH)3 solubility does not fix free Al3+.');
assert(abs(fAl_eq-fAl_expected) <= 1e-10*max(1,fAl_expected), ...
    'Pinned free Al3+ does not satisfy Appelo exchange mass action.');
assert(fAl_eq >= 0 && fAl_eq <= 1, ...
    'Pinned Appelo exchangeable Al is outside the physical site fraction.');

fprintf('  PASS  Al(OH)3 solubility and Appelo Al exchange share one free-Al state\n');
end
