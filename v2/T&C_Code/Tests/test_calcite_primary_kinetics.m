function test_calcite_primary_kinetics
% Primary calcite uses the mineral PSD pathway with calcite chemistry.

test_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(test_dir,'..'));

[f_kin,f_direct,rock_f] = ERW_Feedstock_Pathway("calcite",1);
assert(f_kin == 1 && f_direct == 0 && isequal(rock_f,1), ...
    'Pure calcite must remain in the kinetic mineral pool.');

minerals = ["albite","calcite","quartz"];
fractions = [0.60,0.05,0.30];
[f_kin,f_direct,rock_f] = ERW_Feedstock_Pathway(minerals,fractions);
assert(f_kin == 1 && f_direct == 0 && isequal(rock_f,fractions), ...
    'Mixed feedstock mineral fractions must remain unchanged.');

conv_mol = 1e6;
K_CaCO3 = 10^(-8.35)*conv_mol^2;
Ca = 1000;
CO3 = 0.25*K_CaCO3/Ca;
H = 10^(-6)*conv_mol;
CO2_air = 400/22.41;
M_calcite = 130;
MM_calcite = 100.087/conv_mol;
dtw = 1/144;
k_H = 10^(-0.3)*86400*conv_mol;
k_w = 10^(-5.81)*86400*conv_mol;
k_CO3 = 10^(-3.48)*86400*conv_mol;

[omega,Wr,EW] = ERW_Calcite_Kinetics(Ca,CO3,K_CaCO3,H,CO2_air,conv_mol, ...
    0.30,1,k_H,k_w,k_CO3,1,1,20,1,M_calcite,MM_calcite,dtw);
assert(abs(omega-0.25) < 1e-12,'Calcite saturation is dimensionally incorrect.');
assert(Wr > 0 && EW > 0,'Undersaturated calcite must dissolve.');

[~,~,EW_double_area] = ERW_Calcite_Kinetics(Ca,CO3,K_CaCO3,H,CO2_air,conv_mol, ...
    0.30,1,k_H,k_w,k_CO3,1,1,40,1,M_calcite,MM_calcite,dtw);
assert(abs(EW_double_area/EW-2) < 1e-12, ...
    'Calcite dissolution must scale with reactive surface area.');

[omega_sat,Wr_sat,EW_sat] = ERW_Calcite_Kinetics(Ca,K_CaCO3/Ca, ...
    K_CaCO3,H,CO2_air,conv_mol,0.30,1,k_H,k_w,k_CO3,1,1,20,1, ...
    M_calcite,MM_calcite,dtw);
assert(abs(omega_sat-1) < 1e-12 && Wr_sat == 0 && EW_sat == 0, ...
    'Primary calcite dissolution must stop at saturation.');

[~,Wr_cap,EW_cap] = ERW_Calcite_Kinetics(Ca,0,K_CaCO3,H,CO2_air,conv_mol, ...
    1,1,1e20,1e20,1e20,1,1,20,1,M_calcite,MM_calcite,dtw);
EW_max = M_calcite/MM_calcite/dtw;
assert(EW_cap <= EW_max*(1+1e-12) && ...
    abs(Wr_cap*20-EW_cap) <= max(1,EW_cap)*1e-12, ...
    'Calcite flux must not overdraw mass or particle shrinkage.');

[~,Wr_CO2,~] = ERW_Calcite_Kinetics(Ca,CO3,K_CaCO3,H,CO2_air,conv_mol, ...
    0.30,1,0,0,k_CO3,1,1,20,1,M_calcite,MM_calcite,dtw);
[~,Wr_2CO2,~] = ERW_Calcite_Kinetics(Ca,CO3,K_CaCO3,H,2*CO2_air,conv_mol, ...
    0.30,1,0,0,k_CO3,1,1,20,1,M_calcite,MM_calcite,dtw);
assert(abs(Wr_2CO2/Wr_CO2-2) < 1e-12, ...
    'Calcite carbonate-mechanism kinetics must scale with pCO2.');
expected_CO2_rate = 0.30*k_CO3*(400/1e6)*(1-omega);
assert(abs(Wr_CO2-expected_CO2_rate) <= expected_CO2_rate*1e-12, ...
    'Soil-gas concentration must be converted to dimensionless pCO2.');

did_reject = false;
try
    ERW_Calcite_Kinetics(Ca,CO3,K_CaCO3,0,CO2_air,conv_mol, ...
        0.30,1,k_H,k_w,k_CO3,-1,1,20,1,M_calcite,MM_calcite,dtw);
catch ME
    did_reject = strcmp(ME.identifier,'ERW:invalidCalciteKinetics');
end
assert(did_reject, ...
    'Zero H with a negative reaction order must be rejected.');

k1_c = 10^(-6.35)*conv_mol;
k2_c = 10^(-10.33)*conv_mol;
kw_c = 10^(-14)*conv_mol^2;
water_l = 100;
DIC0 = 2000;
H0 = 10^(-7)*conv_mol;
CO2_0 = DIC0/(1+k1_c/H0+k1_c*k2_c/H0^2);
CO3_0 = k1_c*k2_c*CO2_0/H0^2;
Alk0 = k1_c*CO2_0/H0+2*CO3_0-H0+kw_c/H0;
Omega0 = Ca*CO3_0/K_CaCO3;
[EW_limited,W_limited,diag_diss] = ERW_Calcite_Reaction_Limit( ...
    2e5,1e5,1,Omega0,Ca,DIC0,Alk0,water_l,k1_c,k2_c,kw_c,K_CaCO3,conv_mol);
assert(diag_diss.limited && diag_diss.Omega_final <= 1+1e-8 && ...
    abs(EW_limited/W_limited-2) < 1e-10, ...
    'Combined calcite dissolution must stop at saturation without changing pathway shares.');

Ca_fine = Ca;
DIC_fine = DIC0;
Alk_fine = Alk0;
Omega_fine = Omega0;
extent_fine = 0;
for it = 1:10
    [EW_fine,W_fine,diag_fine] = ERW_Calcite_Reaction_Limit( ...
        2e5,1e5,0.1,Omega_fine,Ca_fine,DIC_fine,Alk_fine,water_l, ...
        k1_c,k2_c,kw_c,K_CaCO3,conv_mol);
    extent_fine = extent_fine+(EW_fine+W_fine)*0.1;
    Ca_fine = diag_fine.Ca_final;
    DIC_fine = diag_fine.DIC_final;
    Alk_fine = diag_fine.Alk_final;
    Omega_fine = diag_fine.Omega_final;
end
assert(abs(extent_fine-diag_diss.extent) <= max(1,diag_diss.extent)*1e-8 && ...
    Omega_fine <= 1+1e-8, ...
    'Calcite equilibrium extent must be insensitive to reaction-step splitting.');

H0 = 10^(-9)*conv_mol;
CO2_0 = DIC0/(1+k1_c/H0+k1_c*k2_c/H0^2);
CO3_0 = k1_c*k2_c*CO2_0/H0^2;
Alk0 = k1_c*CO2_0/H0+2*CO3_0-H0+kw_c/H0;
Omega0 = Ca*CO3_0/K_CaCO3;
[EW_precip,W_precip,diag_precip] = ERW_Calcite_Reaction_Limit( ...
    0,-2e5,1,Omega0,Ca,DIC0,Alk0,water_l,k1_c,k2_c,kw_c,K_CaCO3,conv_mol);
assert(EW_precip == 0 && W_precip < 0 && diag_precip.limited && ...
    diag_precip.Omega_final >= 1-1e-8 && diag_precip.Ca_final >= 0 && ...
    diag_precip.DIC_final >= 0, ...
    'Secondary calcite precipitation must stop at saturation and preserve aqueous inventory.');

did_reject = false;
try
    ERW_Calcite_Reaction_Limit(1,-1,1,Omega0,Ca,DIC0,Alk0,water_l, ...
        k1_c,k2_c,kw_c,K_CaCO3,conv_mol);
catch ME
    did_reject = strcmp(ME.identifier,'ERW:calcitePathConflict');
end
assert(did_reject,'Opposing primary and secondary calcite fluxes must be rejected.');

[solid_total,solid_primary] = ERW_Total_Calcite_Solid( ...
    [60,130,30],minerals,[263,100.087,60]/conv_mol,2e5);
assert(abs(solid_primary-130/MM_calcite) < 1e-8 && ...
    abs(solid_total-solid_primary-2e5) < 1e-8, ...
    'Field primary and secondary calcite inventory is incomplete.');

[solid_control,primary_control] = ERW_Total_Calcite_Solid( ...
    0,strings(1,0),0,2e5);
assert(solid_control == 2e5 && primary_control == 0, ...
    'Control placeholder mineral state must retain secondary carbonate.');

fprintf('PASS: primary calcite uses PSD kinetics, saturation and finite mass.\n');
end
