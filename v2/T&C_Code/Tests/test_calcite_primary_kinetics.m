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
    0.30,1,k_H,k_w,k_CO3,1,1,20,1,M_calcite,MM_calcite,dtw,1);
assert(abs(omega-0.25) < 1e-12,'Calcite saturation is dimensionally incorrect.');
assert(Wr > 0 && EW > 0,'Undersaturated calcite must dissolve.');

[~,~,EW_double_area] = ERW_Calcite_Kinetics(Ca,CO3,K_CaCO3,H,CO2_air,conv_mol, ...
    0.30,1,k_H,k_w,k_CO3,1,1,40,1,M_calcite,MM_calcite,dtw,1);
assert(abs(EW_double_area/EW-2) < 1e-12, ...
    'Calcite dissolution must scale with reactive surface area.');

[omega_sat,Wr_sat,EW_sat] = ERW_Calcite_Kinetics(Ca,K_CaCO3/Ca, ...
    K_CaCO3,H,CO2_air,conv_mol,0.30,1,k_H,k_w,k_CO3,1,1,20,1, ...
    M_calcite,MM_calcite,dtw,1);
assert(abs(omega_sat-1) < 1e-12 && Wr_sat == 0 && EW_sat == 0, ...
    'Primary calcite dissolution must stop at saturation.');

% Activity basis: the caller passes gamma_IAP = gamma_div^2,
% the single factor that turns the code's Ca*CO3 product into the ion
% activity product - the stored CO3 is an activity times gamma_mono^2,
% see POST_LEGACY_EQUATION_AUDIT.md section 3b. Omega is LINEAR in that
% factor, not quadratic, and saturation is reached at the
% correspondingly higher concentration product.
gamma_IAP_test = 0.5;
[omega_act,Wr_act,~] = ERW_Calcite_Kinetics(Ca,CO3,K_CaCO3,H,CO2_air,conv_mol, ...
    0.30,1,k_H,k_w,k_CO3,1,1,20,1,M_calcite,MM_calcite,dtw,gamma_IAP_test);
assert(abs(omega_act-gamma_IAP_test*0.25) < 1e-12 && Wr_act > Wr, ...
    'Calcite saturation on the activity basis must be gamma_IAP times the concentration ratio.');
[omega_act_sat,Wr_act_sat,~] = ERW_Calcite_Kinetics(Ca,K_CaCO3/Ca/gamma_IAP_test, ...
    K_CaCO3,H,CO2_air,conv_mol,0.30,1,k_H,k_w,k_CO3,1,1,20,1, ...
    M_calcite,MM_calcite,dtw,gamma_IAP_test);
assert(abs(omega_act_sat-1) < 1e-12 && Wr_act_sat == 0, ...
    'Activity-basis saturation must stop dissolution at Ca*CO3 = K/gamma_IAP.');
did_reject = false;
try
    ERW_Calcite_Kinetics(Ca,CO3,K_CaCO3,H,CO2_air,conv_mol, ...
        0.30,1,k_H,k_w,k_CO3,1,1,20,1,M_calcite,MM_calcite,dtw,1.5);
catch ME
    did_reject = strcmp(ME.identifier,'ERW:invalidCalciteKinetics');
end
assert(did_reject,'An activity coefficient above 1 must be rejected.');

[~,Wr_cap,EW_cap] = ERW_Calcite_Kinetics(Ca,0,K_CaCO3,H,CO2_air,conv_mol, ...
    1,1,1e20,1e20,1e20,1,1,20,1,M_calcite,MM_calcite,dtw,1);
EW_max = M_calcite/MM_calcite/dtw;
assert(EW_cap <= EW_max*(1+1e-12) && ...
    abs(Wr_cap*20-EW_cap) <= max(1,EW_cap)*1e-12, ...
    'Calcite flux must not overdraw mass or particle shrinkage.');

[~,Wr_CO2,~] = ERW_Calcite_Kinetics(Ca,CO3,K_CaCO3,H,CO2_air,conv_mol, ...
    0.30,1,0,0,k_CO3,1,1,20,1,M_calcite,MM_calcite,dtw,1);
[~,Wr_2CO2,~] = ERW_Calcite_Kinetics(Ca,CO3,K_CaCO3,H,2*CO2_air,conv_mol, ...
    0.30,1,0,0,k_CO3,1,1,20,1,M_calcite,MM_calcite,dtw,1);
assert(abs(Wr_2CO2/Wr_CO2-2) < 1e-12, ...
    'Calcite carbonate-mechanism kinetics must scale with pCO2.');
expected_CO2_rate = 0.30*k_CO3*(400/1e6)*(1-omega);
assert(abs(Wr_CO2-expected_CO2_rate) <= expected_CO2_rate*1e-12, ...
    'Soil-gas concentration must be converted to dimensionless pCO2.');

did_reject = false;
try
    ERW_Calcite_Kinetics(Ca,CO3,K_CaCO3,0,CO2_air,conv_mol, ...
        0.30,1,k_H,k_w,k_CO3,-1,1,20,1,M_calcite,MM_calcite,dtw,1);
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
    2e5,1e5,1,Omega0,Ca,DIC0,Alk0,water_l,k1_c,k2_c,kw_c,K_CaCO3,conv_mol,1);
assert(diag_diss.limited && diag_diss.Omega_final <= 1+1e-8 && ...
    abs(EW_limited/W_limited-2) < 1e-10, ...
    'Combined calcite dissolution must stop at saturation without changing pathway shares.');

% The limiter on the activity basis: with gamma_IAP < 1 the same start state
% is further from equilibrium, so the bound admits a larger extent and
% stops where gamma_IAP*Ca*CO3/K = 1 (a concentration product above K).
Omega0_act = gamma_IAP_test*Omega0;
[EW_act,W_act,diag_act] = ERW_Calcite_Reaction_Limit( ...
    2e5,1e5,1,Omega0_act,Ca,DIC0,Alk0,water_l,k1_c,k2_c,kw_c,K_CaCO3,conv_mol,gamma_IAP_test);
assert(diag_act.limited && diag_act.Omega_final <= 1+1e-8 && ...
    diag_act.Omega_final >= 1-1e-6 && diag_act.extent > diag_diss.extent && ...
    abs(EW_act/W_act-2) < 1e-10, ...
    'The activity-basis bound must stop at gamma_IAP*Ca*CO3/K = 1 with a larger extent.');

Ca_fine = Ca;
DIC_fine = DIC0;
Alk_fine = Alk0;
Omega_fine = Omega0;
extent_fine = 0;
for it = 1:10
    [EW_fine,W_fine,diag_fine] = ERW_Calcite_Reaction_Limit( ...
        2e5,1e5,0.1,Omega_fine,Ca_fine,DIC_fine,Alk_fine,water_l, ...
        k1_c,k2_c,kw_c,K_CaCO3,conv_mol,1);
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
    0,-2e5,1,Omega0,Ca,DIC0,Alk0,water_l,k1_c,k2_c,kw_c,K_CaCO3,conv_mol,1);
assert(EW_precip == 0 && W_precip < 0 && diag_precip.limited && ...
    diag_precip.Omega_final >= 1-1e-8 && diag_precip.Ca_final >= 0 && ...
    diag_precip.DIC_final >= 0, ...
    'Secondary calcite precipitation must stop at saturation and preserve aqueous inventory.');

did_reject = false;
try
    ERW_Calcite_Reaction_Limit(1,-1,1,Omega0,Ca,DIC0,Alk0,water_l, ...
        k1_c,k2_c,kw_c,K_CaCO3,conv_mol,1);
catch ME
    did_reject = strcmp(ME.identifier,'ERW:calcitePathConflict');
end
assert(did_reject,'Opposing primary and secondary calcite fluxes must be rejected.');

%% Regression: precipitation that would consume the whole dissolved
%% inventory. The limiter caps the extent at min(Ca,DIC)*water and forms
%% the trial state as DIC + (-DIC*water)/water, which rounds a few ulps
%% below zero for some water volumes. That is an exhausted pool, not an
%% invalid state; the limiter must then bisect to saturation instead of
%% refusing the step (the failure seen in the margin_prairie runs of the
%% calibrated all-12 campaign on 2026-09-13).
%% Two pool sizes: pore-water concentrations, where the rounding residual
%% is far inside the historical absolute tolerance, and a molar pool, where
%% one ulp exceeds that tolerance and only the ulp-based allowance applies.
for exhaustion_case = {[101.3 2000],[1e6 2e6]}
    DIC_low = exhaustion_case{1}(1);
    Ca_high = exhaustion_case{1}(2);
    H0 = 10^(-9.3)*conv_mol;
    CO2_0 = DIC_low/(1+k1_c/H0+k1_c*k2_c/H0^2);
    CO3_0 = k1_c*k2_c*CO2_0/H0^2;
    Alk_low = k1_c*CO2_0/H0+2*CO3_0-H0+kw_c/H0;
    Omega_high = Ca_high*CO3_0/K_CaCO3;
    assert(Omega_high > 1,'The exhaustion regression case must start supersaturated.');
    rounding_volumes = [];
    for water_trial = 0.1:0.1:60
        if DIC_low+(-(DIC_low*water_trial))/water_trial < 0
            rounding_volumes(end+1) = water_trial; %#ok<AGROW>
        end
    end
    assert(~isempty(rounding_volumes), ...
        'At least one water volume must round the exhausted pool below zero.');
    for water_trial = [rounding_volumes(1) rounding_volumes(end) 1 10]
        inventory = DIC_low*water_trial;
        [EW_exh,W_exh,diag_exh] = ERW_Calcite_Reaction_Limit( ...
            0,-1e3*inventory,1,Omega_high,Ca_high,DIC_low,Alk_low,water_trial, ...
            k1_c,k2_c,kw_c,K_CaCO3,conv_mol,1);
        assert(EW_exh == 0 && W_exh < 0 && diag_exh.limited && ...
            -diag_exh.extent < inventory && ...
            abs(diag_exh.Omega_final-1) <= 1e-8 && ...
            diag_exh.DIC_final > 0 && diag_exh.Ca_final > 0, ...
            'Precipitation exceeding the dissolved inventory must stop at saturation.');
    end
end

[solid_total,solid_primary] = ERW_Total_Calcite_Solid( ...
    [60,130,30],minerals,[263,100.087,60]/conv_mol,2e5);
assert(abs(solid_primary-130/MM_calcite) < 1e-8 && ...
    abs(solid_total-solid_primary-2e5) < 1e-8, ...
    'Field primary and secondary calcite inventory is incomplete.');

[solid_control,primary_control] = ERW_Total_Calcite_Solid( ...
    0,strings(1,0),0,2e5);
assert(solid_control == 2e5 && primary_control == 0, ...
    'Control placeholder mineral state must retain secondary carbonate.');

% The production callers must condition the constants at the step's ionic
% strength and form the concentration-basis activity factor at every
% calcite saturation site, main step and initialisation alike, and must
% hand the rate law the proton activity.
code_dir = fullfile(test_dir,'..');
erw_source = fileread(fullfile(code_dir,'ERW_function.m'));
for needed = {'[k_cond,K_AlOH3,gamma_z] = ERW_Conditional_Constants( ...', ...
        'gamma_mono = gamma_z(1);', ...
        'gamma_div = gamma_z(2);', ...
        'gamma_IAP = gamma_div^2;', ...
        'Omega_CaCO3 = gamma_IAP*Ca*CO3/K_CaCO3;', ...
        'Omega_MgCO3 = gamma_IAP*Mg*CO3/K_MgCO3;', ...
        'ERW_Calcite_Kinetics(Ca,CO3,K_CaCO3,gamma_mono*H,CO2_air,conv_mol,s,diss_f, ...', ...
        'MM_min_row(min_j),dtw,gamma_IAP);', ...
        'n*Zrw*s*1000,k1,k2,k_w,K_CaCO3,conv_mol,gamma_IAP);'}
    assert(contains(erw_source,needed{1}), ...
        'ERW_function.m must carry: %s',needed{1});
end
param_source = fileread(fullfile(code_dir,'Parameters','ERW_MOD_PARAM_excel.m'));
for needed = {'[k_seed,K_AlOH3_s,gamma_seed] = ERW_Conditional_Constants( ...', ...
        'gamma_mono_0 = gamma_seed(1);', ...
        'gamma_div_0 = gamma_seed(2);', ...
        'gamma_IAP_0 = gamma_div_0^2;', ...
        'Omega_CaCO3(1) = gamma_IAP_0*Ca(1)*CO3(1)/K_CaCO3;', ...
        'Omega_MgCO3(1) = gamma_IAP_0*Mg(1)*CO3(1)/K_MgCO3;', ...
        'Ca(1),CO3(1),K_CaCO3,gamma_mono_0*H(1),CO2_air(1),conv_mol,s(1),diss_f, ...', ...
        'M_min(1,min_j),MM_min(min_j),dtw_init,gamma_IAP_0);', ...
        'n(1)*Zrw*s(1)*1000,k1_s,k2_s,k_w_s,K_CaCO3,conv_mol,gamma_IAP_0);'}
    assert(contains(param_source,needed{1}), ...
        'ERW_MOD_PARAM_excel.m must carry: %s',needed{1});
end
for stale = {'gamma_div/gamma_mono^2','gamma_div_0/gamma_mono_0^2'}
    assert(~contains(erw_source,stale{1}) && ~contains(param_source,stale{1}), ...
        'The hybrid-CO3 factor %s must not survive the conditional-constant basis.',stale{1});
end

fprintf('PASS: primary calcite uses PSD kinetics, saturation and finite mass.\n');
end
