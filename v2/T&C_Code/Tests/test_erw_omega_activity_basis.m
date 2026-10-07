function test_erw_omega_activity_basis
%%% The carbonate saturation index is formed on the ion activity product.
%%%
%%% ERW_function keeps every species a concentration and conditions each
%%% mass-action constant at the step's ionic strength
%%% (ERW_Conditional_Constants: k1/g1^2, k2/g2, k_w/g1^2, ...).  CO3 formed
%%% from those constants is therefore a concentration like Ca and Mg, and
%%% the calcite ion activity product is
%%%
%%%     a_Ca*a_CO3 = (gamma_div*Ca) * (gamma_div*CO3)
%%%
%%% so Omega = gamma_IAP*Ca*CO3/K_CaCO3 with gamma_IAP = gamma_div^2 and the
%%% thermodynamic K_CaCO3.  The 2026-09-23 revision paired thermodynamic
%%% k1, k2 with concentrations, which made CO3 an activity times gamma_mono^2
%%% and needed gamma_div/gamma_mono^2; that hybrid factor must not survive
%%% the conditional basis, where it would suppress Omega by gamma_mono^6.

test_dir = fileparts(mfilename('fullpath'));
code_dir = fullfile(test_dir,'..');
addpath(code_dir);

%% --- source pins: one conditioning site per solve, one factor per index
erw = fileread(fullfile(code_dir,'ERW_function.m'));
par = fileread(fullfile(code_dir,'Parameters','ERW_MOD_PARAM_excel.m'));
kin = fileread(fullfile(code_dir,'ERW_Calcite_Kinetics.m'));
lim = fileread(fullfile(code_dir,'ERW_Calcite_Reaction_Limit.m'));

assert(contains(erw,'[k_cond,K_AlOH3,gamma_z] = ERW_Conditional_Constants( ...') && ...
    contains(erw,'k1 = k_cond(1); k2 = k_cond(2); k_w = k_cond(3);') && ...
    contains(erw,'gamma_IAP = gamma_div^2;'), ...
    'ERW_function.m must condition its constants and form gamma_IAP as gamma_div^2.');
assert(contains(erw,'Omega_CaCO3 = gamma_IAP*Ca*CO3/K_CaCO3;') && ...
    contains(erw,'Omega_MgCO3 = gamma_IAP*Mg*CO3/K_MgCO3;'), ...
    'Both carbonate saturation indices must carry gamma_IAP.');
assert(contains(erw,'ERW_biogeochem.k1  = k1;') && ...
    contains(erw,'ERW_biogeochem.k_w = k_w;') && ...
    contains(par,'k_w(1) = k_w_s;'), ...
    ['The record must carry the conditioned constants, so saved ' ...
     'concentrations re-speciate exactly.']);
assert(contains(erw,'k1_thermo*CO2_w_rain/H_rain'), ...
    'Rain is dilute and must keep the thermodynamic constants.');
assert(contains(par,'gamma_IAP_0 = gamma_div_0^2;') && ...
    contains(par,'Omega_CaCO3(1) = gamma_IAP_0*Ca(1)*CO3(1)/K_CaCO3;') && ...
    contains(par,'Omega_MgCO3(1) = gamma_IAP_0*Mg(1)*CO3(1)/K_MgCO3;') && ...
    contains(par,'HCO3(1)= k1_s*CO2_w(1)/H(1);'), ...
    'The initial condition must use the same basis as every later step.');
assert(contains(kin,'Omega = gamma_IAP*Ca*CO3/K_CaCO3;'), ...
    'ERW_Calcite_Kinetics must form Omega on the ion activity product.');
assert(contains(lim,'Omega = gamma_IAP*Ca_out*CO3_trial/K_CaCO3;'), ...
    'ERW_Calcite_Reaction_Limit must form its trial Omega on the same basis.');
assert(~contains(kin,'gamma_div/gamma_mono^2') && ...
    ~contains(lim,'gamma_div/gamma_mono^2'), ...
    'Calcite helper guidance must use the conditional-constant activity factor.');
for src = {erw, par}
    assert(~contains(src{1},'gamma_div/gamma_mono^2') && ...
        ~contains(src{1},'gamma_div_0/gamma_mono_0^2'), ...
        'The hybrid-CO3 factor must not survive the conditional basis.');
end

%% --- Davies identity: log10(gamma) scales as z^2 with one shared f(I), so
%%     gamma_div = gamma_mono^4 exactly.  Pinned so that swapping in an
%%     ion-specific activity model (b-dot, Pitzer) fails here rather than
%%     silently changing Omega.
for I = [0 1e-4 1e-3 5e-3 1e-2 2e-2 5e-2 1e-1]
    g1 = ERW_Activity_Coefficient(1,I);
    g2 = ERW_Activity_Coefficient(2,I);
    assert(abs(g2 - g1^4) <= 8*eps(max(g2,1)), ...
        'Davies requires gamma(2) == gamma(1)^4 at I = %g.',I);
    assert(g2^2 > 0 && g2^2 <= 1, 'gamma_IAP must lie in (0,1] at I = %g.',I);
end

%% --- the kinetics routine is linear in gamma_IAP, and gamma_IAP is required
Ca = 1.0e3; CO3 = 5.0; K_CaCO3 = 3.31e3; H = 1.0;
CO2_air = 1.0e4; conv_mol = 1e6; s = 0.7; diss_f = 1;
k_H = 1e-9; k_w = 1e-9; k_CO2 = 1e-9; n_H = 1.0; n_CO2 = 0;
SA = 1.0; rock_fraction = 0.1; M_calcite = 10; MM_calcite = 100.087; dtw = 1;
args = {Ca,CO3,K_CaCO3,H,CO2_air,conv_mol,s,diss_f,k_H,k_w,k_CO2, ...
    n_H,n_CO2,SA,rock_fraction,M_calcite,MM_calcite,dtw};
threw = false;
try
    ERW_Calcite_Kinetics(args{:});
catch err
    threw = strcmp(err.identifier,'ERW:invalidCalciteKinetics');
end
assert(threw, 'Omitting gamma_IAP must raise ERW:invalidCalciteKinetics.');
Omega_unit = ERW_Calcite_Kinetics(args{:},1);
assert(abs(Omega_unit - Ca*CO3/K_CaCO3) <= 8*eps(Omega_unit), ...
    'gamma_IAP = 1, passed explicitly, is the ideal-solution Omega.');
g = ERW_Activity_Coefficient(2,2e-2)^2;
Omega_g = ERW_Calcite_Kinetics(args{:},g);
assert(abs(Omega_g - g*Omega_unit) <= 64*eps(max(Omega_g,1)), ...
    'Omega must be linear in gamma_IAP.');

%% --- the end-to-end statement: speciated with the conditional constants,
%%     gamma_div^2*Ca*CO3 IS the ion activity product, reached here by an
%%     independent route (a_CO3 built from a_H and a_CO2 with the
%%     thermodynamic K1, K2).
c = 1e6;
K1 = 10^-6.35; K2 = 10^-10.33; Kw = 1e-14;
aH = 1e-6; aCO2 = 3.4e-4; Ca_c = 1.0e-3;
for I = [5e-3 2e-2 5e-2 0.126]
    [kc,~,gz] = ERW_Conditional_Constants( ...
        [K1*c,K2*c,Kw*c^2,1e-5*c,1e-5*c,1e-7*c,1e-6*c],[],I);
    H_c = aH/gz(1)*c;                          % solver H, umol/L
    CO3_code = kc(1)*kc(2)*(aCO2*c)/H_c^2;      % solver CO3, umol/L
    IAP_code = gz(2)^2*(Ca_c*c)*CO3_code/c^2;
    IAP_true = (gz(2)*Ca_c)*(K1*K2*aCO2/aH^2);
    assert(abs(IAP_code/IAP_true - 1) <= 1e-12, ...
        'gamma_div^2*Ca*CO3 must equal a_Ca*a_CO3 at I = %g.',I);
    assert(abs(CO3_code/c - K1*K2*aCO2/aH^2/gz(2)) <= 1e-12*CO3_code/c, ...
        'The solver CO3 must be the carbonate concentration a_CO3/gamma_div at I = %g.',I);
end

%% Calcite solubility follows soil temperature like k1, k2 and k_w do.
log_Ksp = @(T) -171.9065-0.077993*T+2839.319./T+71.595*log10(T);
assert(abs(log_Ksp(298.15)+8.48) < 5e-3 && abs(log_Ksp(283.15)+8.41) < 5e-3, ...
    'The calcite polynomial must reproduce Plummer and Busenberg (1982).');
assert(contains(erw,'-171.9065-0.077993*T_K') && contains(par,'-171.9065-0.077993*T_K') && ...
    ~contains(par,'K_CaCO3 = 10^(-8.35)') && contains(par,'''Omega_CaCO3''') && ...
    contains(erw,'ERW_biogeochem.Omega_CaCO3 = Omega_CaCO3;'), ...
    'Calcite solubility must follow soil temperature and each step must record its Omega.');

fprintf('PASS: carbonate saturation is formed on the ion activity product.\n');
end
