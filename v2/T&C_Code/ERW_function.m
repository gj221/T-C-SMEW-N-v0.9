%% ERW module

function [ERW_data, ERW_TC, ERW_biogeochem, ERW_weathering, ERW_plant,ERW_var, ERW_const] =...
    ERW_function(n, s_pre, s, Lw, Tw, Iw, ...
        Zrw, r_het, r_aut, temp_soil, K_CEC, CEC_tot, ...
        M_rock_in, k_diss_H, k_diss_w, k_diss_OH, E_H, E_w, E_OH, R, T_ref, K1, K2, K3, K4, ...
        rho_rock, b, K_CaCO3, K_MgCO3, MM_min, ...
        min_st, dtw, Al_CEC_scale, conv_mol, ...
        ~, ~, Na, ~, An,...
        fsolve_options, Dw_0, D_0, pH_rain, CO2_atm, CO2_air, Z_CO2, ...
        Ca_tot,Mg_tot, K_tot,Na_tot,An_tot,Al_tot,Si_tot, ...
        Al, AlOH4, ~, ...
        R_alk,H, ...
        EW, IC_tot, Wr, psd, d, lamb, SSA,...
        Fs, DIC, ...
        I_Ca, I_Mg, I_Na, I_K, I_An, ~, W_CaCO3, W_MgCO3, j_only, ...
        CaCO3, MgCO3,...
        f_CEC,f_Ca,f_Mg,f_K,f_Na,f_Al,f_H, ...
        r_CaCO3,r_MgCO3, tau_CaCO3, tau_MgCO3,...
        M_min, M_iner, a, psd_0, d_0, lamb_0, time_step, mineral, K_sp, diss_f,n_H,n_OH,...
        ~,...
        dCa_tot,dMg_tot,dK_tot,dSi_tot,...
        ~, M_rock_add, psd_perc_in, rock_f_in, ~, rock_app,...
        f_Lk_biogeo,f_T_biogeo,...
        PO4_tot,~,dPO4_tot,...
        NO3_tot,~,dNO3_tot,...
        NH4_tot,~,dNH4_tot,...
        keyword_add,...
        CaCO3_rock,...
        K_AlOH3, bulk_pH_params)

%%%% INPUT

%     pH: pH value of the solution (unitless)
%     CO2_air: Concentration of CO2 in soil air (??mol/L). The partial
%              pressure is CO2_air*22.41/conv_mol [atm], the conversion
%              ERW_Calcite_Kinetics uses; 44,623 ??mol/L is one atmosphere.
%     Alk: Alkalinity (mol/L)
%     DIC: Dissolved inorganic carbon concentration (??mol/L)
%     Fs: Advection influx of CO2 from the air into the soil (mol-conv/m??)
%     An: Concentration of anions (??mol/L)
%     CO2_w: Concentration of dissolved CO2 in water (??mol/L)
%     Omega_CaCO3: Saturation state of calcium carbonate (unitless)
%     Omega_MgCO3: Saturation state of magnesium carbonate (unitless)
%     CaCO3: Concentration of calcium carbonate (??mol/L)
%     MgCO3: Concentration of magnesium carbonate (??mol/L)
%     AlOH: Concentration of aluminum hydroxide (mol/L)
%     AlOH2: Concentration of aluminum hydroxide (mol/L)
%     AlOH3: Concentration of aluminum hydroxide (mol/L)
%     AlOH4: Concentration of aluminum hydroxide (mol/L)
%     Al_w: Concentration of weathered aluminium (mol/L)
%     Al: Concentration of aluminum (mol/L)
%     Mg: Concentration of magnesium (mol/L)
%     Ca: Concentration of calcium (mol/L)
%     Na: Concentration of sodium (mol/L)
%     K: Concentration of potassium (mol/L)
%     Si: Concentration of silicon (??mol/L)
%     f_Al: Fraction of aluminum
%     f_Mg: Fraction of magnesium
%     f_Na: Fraction of sodium
%     f_K: Fraction of potassium
%     f_H: Fraction of hydrogen
%     f_Ca: Fraction of calcium
%     M_min: Mass of rock (g/m2)
%     CO2_w_rain: CO2 concentration in rainwater (??mol/L)
%     H_rain: Concentration of hydrogen ions in rainwater (??mol/L)
%     DIC_rain: Dissolved inorganic carbon concentration in rainwater (??mol/L)
%     s: Soil moisture (-)
%     Tw: Tanspiration rate (m/d)
%     Lw: Leaching rate (m/d)
%     r_het: Heterotrophic respiration rate (mol/m^2/s)
%     r_aut: Autotrophic respiration rate (mol/m^2/s)
%     temp_soil: Soil temperature (??C)
%     Dw: Water diffusivity (m2/d)
%     f_CEC: Fraction of cation exchange capacity
%     Mg_tot: Total magnesium concentration (mol)
%     Ca_tot: Total calcium concentration (mol)
%     K_tot: Total potassium concentration (mol)
%     Na_tot: Total sodium concentration (mol)
%     An_tot: Total anions (mol)
%     Al_tot: Total aluminum concentration (mol)
%     Si_tot: Total silicon concentration (mol)
%     IC_tot: Total inorganic carbon concentration (mol)
%     R_alk: Residual-alkalinity (??eq/L)
%     H: Concentration of hydrogen ions (??mol/L)
%     Z_CO2: Diffusion path length for the soil-atmosphere CO2 flux (m),
%            set to Zrw/2 in ERW_MOD_PARAM_excel.m. Not a concentration.


%%% OUTPUT

%%% `ERW_data` Structure:
% `pH`: pH value of the solution (unitless)
% `CO2_air`: Concentration of CO2 in air (??mol/L)
% `HCO3`: Concentration of bicarbonate (HCO3-) ions (??mol/L)
% `CO3`: Concentration of carbonate (CO3%-) ions (??mol/L)
% `DIC`: Dissolved inorganic carbon concentration (??mol/L)
% `Fs`: Advection influx of CO2 from the air into the soil (mol-conv/m??)
% `An`: COncentration of anions (??mol/L)
% `Alk`: Alkalinity (??eq/L)
% `CO2_w`: Concentration of dissolved CO2 in water (??mol/L)
% `Omega_CaCO3`: Saturation state of calcium carbonate (unitless)
% `Omega_MgCO3`: Saturation state of magnesium carbonate (unitless)
% `CaCO3`: Concentration of calcium carbonate (??mol/L)
% `MgCO3`: Concentration of magnesium carbonate (??mol/L)
% `AlOH`: Concentration of aluminum hydroxide (??mol/L)
% `AlOH2`: Concentration of aluminum hydroxide (??mol/L)
% `AlOH3`: Concentration of aluminum hydroxide (??mol/L)
% `AlOH4`: Concentration of aluminum hydroxide (??mol/L)
% `Al_w`: Concentration of weathered aluminium (??mol/L)
% `Al`: Concentration of aluminum (??mol/L)
% `Mg`: Concentration of magnesium (??mol/L)
% `Ca`: Concentration of calcium (??mol/L)
% `Na`: Concentration of sodium (??mol/L)
% `K`: Concentration of potassium (??mol/L)
% `Si`: Concentration of silicon (??mol/L)
% `f_Al`: Fraction of aluminum
% `f_Mg`: Fraction of magnesium
% `f_Na`: Fraction of sodium
% `f_K`: Fraction of potassium
% `f_H`: Fraction of hydrogen
% `f_Ca`: Fraction of calcium
% `M_min`: Mineralization rate (mol/m^2/s)
% `CO2_w_rain`: CO2 concentration in rainwater (??mol/L)
% `H_rain`: Concentration of hydrogen ions in rainwater (??mol/L)
% `DIC_rain`: Dissolved inorganic carbon concentration in rainwater (??mol/L)

%%% `ERW_TC` Structure:
% `s`: Soil moisture (-)
% `Tw`: Tanspiration rate (m/d)
% `Lw`: Leaching rate (m/d)
% `r_het`: Heterotrophic respiration rate (mol/m^2/s)
% `r_aut`: Autotrophic respiration rate (mol/m^2/s)
% `temp_soil`: Soil temperature (??C)

%%% `ERW_biogeochem` Structure:
% `Dw`: Water diffusivity (m2/d)
% `Fs`: Advection influx of CO2 from the air into the soil (mol-conv/m??)
% `CO2_air`: Concentration of CO2 in soil air (??mol/L)
% `temp_soil`: Soil temperature (??C)
% `f_CEC`: Fraction of cation exchange capacity
% `Mg_tot`: Total magnesium concentration (??mol/L)
% `Ca_tot`: Total calcium concentration (??mol/L)
% `K_tot`: Total potassium concentration (??mol/L)
% `Na_tot`: Total sodium concentration (??mol/L)
% `An_tot`: Total anions (??eq/L)
% `Al_tot`: Total aluminum concentration (??mol/L)
% `Si_tot`: Total silicon concentration (??mol/L)
% `IC_tot`: Total inorganic carbon concentration (??mol/L)
% `R_alk`: Alkalinity (??eq/L)
% `H`: Concentration of hydrogen ions (??mol/L)
% `Z_CO2`: Diffusion path length for the soil-atmosphere CO2 flux (m),
%          set to Zrw/2. Enters the gas flux as D/Z_CO2; not a concentration.
% `DIC`: Dissolved inorganic carbon concentration (??mol/L)
% `AlOH`, `AlOH2`, `AlOH3`, `AlOH4`: Concentration of aluminum hydroxide species (??mol/L)
% `Alk`: Alkalinity (??eq/L)
% `CO2_w`: Concentration of dissolved CO2 in water (??mol/L)
% `Omega_CaCO3`: Saturation state of calcium carbonate (unitless)
% `Omega_MgCO3`: Saturation state of magnesium carbonate (unitless)
% `CaCO3`: Concentration of calcium carbonate (??mol/L)
% `MgCO3`: Concentration of magnesium carbonate (??mol/L)
% `Al_w`: Concentration of aluminum in water (??mol/L)
% `Al`: Concentration of aluminum (??mol/L)
% `Mg`: Concentration of magnesium (??mol/L)
% `Ca`: Concentration of calcium (??mol/L)
% `Na`: Concentration of sodium (??mol/L)
% `K`: Concentration of potassium (??mol/L)
% `An`: Concentration of anions (??eq/L)
% `Si`: Concentration of silicon (??mol/L)
% `f_Al`, `f_Mg`, `f_Na`, `f_K`, `f_H`, `f_Ca`: Fraction of respective elements
% `M_min`: Mass of mineral (g/mol)

%%% `ERW_weathering` Structure:
% `SSA`: Specific surface area (m^2/g)
% `Omega`: Saturation state (unitless)
% `EW`: Enhanced weathering rate (umol mineral/m^2/day)
% `Wr`: Weathering rate (mol/m^2/s)
% `lamb`: Weathering constant (unitless)
% `W_CaCO3`: Weathering rate of calcium carbonate (mol/m^2/s)
% `W_MgCO3`: Weathering rate of magnesium carbonate (mol/m^2/s)

%%% `ERW_plant` Structure:
% `GF_mn`: Growth factor (unitless)

%%% `ERW_var` Structure:
% `errors`: Errors (unitless)

%%

%%% Ionic-strength closure of the pore-water chemistry.  ERW_function_at_I
%%% solves the step with every mass-action constant conditioned at a given
%%% ionic strength (ERW_Conditional_Constants); NaN takes it from the
%%% incoming state.  When the solved state's own ionic strength differs from
%%% it by more than 1 %, the step is solved again from the SAME inputs at the
%%% solved value, so no flux or amendment is applied twice.  A 1 % error in I
%%% moves the reported pH by less than 0.001.
activity_I = NaN;
for activity_pass = 1:3
    [ERW_data, ERW_TC, ERW_biogeochem, ERW_weathering, ERW_plant, ERW_var, ERW_const] = ...
        ERW_function_at_I(activity_I, n, s_pre, s, Lw, Tw, Iw, ...
        Zrw, r_het, r_aut, temp_soil, K_CEC, CEC_tot, ...
        M_rock_in, k_diss_H, k_diss_w, k_diss_OH, E_H, E_w, E_OH, R, T_ref, K1, K2, K3, K4, ...
        rho_rock, b, K_CaCO3, K_MgCO3, MM_min, ...
        min_st, dtw, Al_CEC_scale, conv_mol, ...
        [], [], Na, [], An,...
        fsolve_options, Dw_0, D_0, pH_rain, CO2_atm, CO2_air, Z_CO2, ...
        Ca_tot,Mg_tot, K_tot,Na_tot,An_tot,Al_tot,Si_tot, ...
        Al, AlOH4, [], ...
        R_alk,H, ...
        EW, IC_tot, Wr, psd, d, lamb, SSA,...
        Fs, DIC, ...
        I_Ca, I_Mg, I_Na, I_K, I_An, [], W_CaCO3, W_MgCO3, j_only, ...
        CaCO3, MgCO3,...
        f_CEC,f_Ca,f_Mg,f_K,f_Na,f_Al,f_H, ...
        r_CaCO3,r_MgCO3, tau_CaCO3, tau_MgCO3,...
        M_min, M_iner, a, psd_0, d_0, lamb_0, time_step, mineral, K_sp, diss_f,n_H,n_OH,...
        [],...
        dCa_tot,dMg_tot,dK_tot,dSi_tot,...
        [], M_rock_add, psd_perc_in, rock_f_in, [], rock_app,...
        f_Lk_biogeo,f_T_biogeo,...
        PO4_tot,[],dPO4_tot,...
        NO3_tot,[],dNO3_tot,...
        NH4_tot,[],dNH4_tot,...
        keyword_add,...
        CaCO3_rock,...
        K_AlOH3, bulk_pH_params);
    I_used = ERW_biogeochem.activity_I;
    I_solved = ERW_biogeochem.ionic_str;
    if abs(I_solved-I_used) <= 1e-2*max(I_solved,1e-4)
        break
    end
    activity_I = I_solved;
end
if abs(I_solved-I_used) > 1e-2*max(I_solved,1e-4)
    error('ERW:activityClosure', ...
        'Step ionic strength did not close in %d passes: used %.6g, solved %.6g mol/L.', ...
        activity_pass,I_used,I_solved);
end
ERW_biogeochem.activity_passes = activity_pass;
end

%%% One step at a fixed activity ionic strength; inputs as ERW_function.
function [ERW_data, ERW_TC, ERW_biogeochem, ERW_weathering, ERW_plant,ERW_var, ERW_const] =...
    ERW_function_at_I(activity_I_in, n, s_pre, s, Lw, Tw, Iw, ...
        Zrw, r_het, r_aut, temp_soil, K_CEC, CEC_tot, ...
        M_rock_in, k_diss_H, k_diss_w, k_diss_OH, E_H, E_w, E_OH, R, T_ref, K1, K2, K3, K4, ...
        rho_rock, b, K_CaCO3, K_MgCO3, MM_min, ...
        min_st, dtw, Al_CEC_scale, conv_mol, ...
        ~, ~, Na, ~, An,...
        fsolve_options, Dw_0, D_0, pH_rain, CO2_atm, CO2_air, Z_CO2, ...
        Ca_tot,Mg_tot, K_tot,Na_tot,An_tot,Al_tot,Si_tot, ...
        Al, AlOH4, ~, ...
        R_alk,H, ...
        EW, IC_tot, Wr, psd, d, lamb, SSA,...
        Fs, DIC, ...
        I_Ca, I_Mg, I_Na, I_K, I_An, ~, W_CaCO3, W_MgCO3, j_only, ...
        CaCO3, MgCO3,...
        f_CEC,f_Ca,f_Mg,f_K,f_Na,f_Al,f_H, ...
        r_CaCO3,r_MgCO3, tau_CaCO3, tau_MgCO3,...
        M_min, M_iner, a, psd_0, d_0, lamb_0, time_step, mineral, K_sp, diss_f,n_H,n_OH,...
        ~,...
        dCa_tot,dMg_tot,dK_tot,dSi_tot,...
        ~, M_rock_add, psd_perc_in, rock_f_in, ~, rock_app,...
        f_Lk_biogeo,f_T_biogeo,...
        PO4_tot,~,dPO4_tot,...
        NO3_tot,~,dNO3_tot,...
        NH4_tot,~,dNH4_tot,...
        keyword_add,...
        CaCO3_rock,...
        K_AlOH3, bulk_pH_params)

% CO2 advection due to moisture variation
ADV = 0;
if s <= s_pre
    ADV = n*Zrw*1000*(s-s_pre)*CO2_atm;
elseif s > s_pre
    ADV = n*Zrw*1000*(s-s_pre)*CO2_air;
end

%% state variables

% soil temperature constant
T_K = temp_soil+273.15;
if ~isfinite(T_K) || T_K < 100 || T_K > 400
    error('ERW:InvalidTemperature', ...
        'Step %d: invalid soil temperature %.6g K.',time_step,T_K);
end

% solute diffusivity in water
Dw = Dw_0*(n*s)^2;% Archie 1942, Grathwohl 1998 (book)

% soil CO2 diffusivity
D = D_0*(1-s)^(10/3)*n^(4/3);%Mill-Quirk (1961)

% species-dependent mineral dissolution rate
if rock_app == 0
    k_diss_H_t(1,:) = 0;
    k_diss_w_t(1,:) = 0;
    k_diss_OH_t(1,:) = 0;
else
    k_diss_H_t(1,:) = k_diss_H.*exp(-E_H.*1000/R*(1/T_K-1/T_ref));
    k_diss_w_t(1,:) = k_diss_w.*exp(-E_w.*1000/R*(1/T_K-1/T_ref));
    k_diss_OH_t(1,:) = k_diss_OH.*exp(-E_OH.*1000/R*(1/T_K-1/T_ref));
end


% carbonate speciation
pk1 = -(-356.309-0.0609*T_K+21834.37/T_K+126.8339*log10(T_K)-1684915/(T_K)^2); % acid dissociation constants
pk2 = -(-107.887-0.032528*(T_K)+5151.79/(T_K)+38.92561*log10(T_K)-563713.9/(T_K)^2);
pk_w = -(-283.971+13323/(T_K)-0.0507*(T_K)+102.24447*log10(T_K)-1119669/(T_K)^2); % acid-water dissociation constant
k1 = 10^(-pk1)*conv_mol; % dissociation constant
k2 = 10^(-pk2)*conv_mol;
pk_CaCO3 = -(-171.9065-0.077993*T_K+2839.319/T_K+71.595*log10(T_K)); % calcite solubility, Plummer & Busenberg (1982)
K_CaCO3 = 10^(-pk_CaCO3)*conv_mol^2;
K_MgCO3 = 10^(-5.27)*exp(17170/R*(1/T_K-1/298.15))*conv_mol^2; % nesquehonite solubility, Harrison et al. (2019); magnesite does not precipitate at soil temperature

if ~isfinite(k1) || k1 <= 0 || ~isfinite(k2) || k2 <= 0
    error('ERW:InvalidEquilibriumConstant', ...
        'Step %d: invalid carbonate constants k1=%.6g, k2=%.6g.',time_step,k1,k2);
end

% water dissociation constant
k_w = 10^(-pk_w)*conv_mol^2;
if ~isfinite(k_w) || k_w <= 0
    error('ERW:InvalidEquilibriumConstant', ...
        'Step %d: invalid water constant k_w=%.6g.',time_step,k_w);
end

% ammonia dissociation constant
pk_nh3 = 0.09018 + 2729.92/T_K;
k_nh3 = 10^(-pk_nh3)*conv_mol;
if ~isfinite(k_nh3) || k_nh3 <= 0
    error('ERW:InvalidEquilibriumConstant', ...
        'Step %d: invalid ammonia constant k_nh3=%.6g.',time_step,k_nh3);
end

% Henry - constant for CO2 solubility
k_H = 0.83*exp(2400*(1/T_K-1/T_ref));


%CEC Gaines-Thomas constants
K_CEC_ = num2cell(K_CEC);
[K_Ca_Mg,K_Ca_K,K_Ca_Na,K_Ca_Al,K_Ca_H]= K_CEC_{:};

%%% Cation-exchange scheme used by the daily and extract solvers.
scheme_str = 'TC_GT';
if isstruct(bulk_pH_params) && isfield(bulk_pH_params, 'scheme') && ...
        ~isempty(bulk_pH_params.scheme)
    scheme_str = bulk_pH_params.scheme;
    if isstring(scheme_str) && isscalar(scheme_str)
        scheme_str = char(scheme_str);
    end
    if ~ischar(scheme_str) || ~any(strcmpi(scheme_str, {'TC_GT','appelo'}))
        error('ERW:invalidParam', 'Exchange scheme must be TC_GT or appelo.');
    end
end
scheme_is_appelo = strcmpi(scheme_str, 'appelo');
appelo_K = [];
if scheme_is_appelo
    appelo_params = [];
    if isfield(bulk_pH_params, 'appelo')
        appelo_params = bulk_pH_params.appelo;
    end
    [K_Na_H, K_K_H, K_Ca_H2, K_Mg_H2, K_Al_H3, appelo_alpha] = ...
        Resolve_Appelo_Constants(appelo_params);
    appelo_K = [K_Al_H3,K_Mg_H2,K_Na_H,K_K_H,K_Ca_H2];
else
    appelo_alpha = [];
end

%%% Al(OH)3 solubility pin (SCEPTER 'amal' semantics).  When enabled under
%%% the Appelo scheme (bulk_pH_params.Al_solubility_pin == true), aqueous
%%% Al is fixed by amorphous Al(OH)3 solubility:
%%%   Al(OH)3 + 3H+ = Al3+ + 3H2O,   *Kso = [Al3+]/[H+]^3 = 10^K_AlOH3
%%%   Al3_model = 10^K_AlOH3 * (H_umol/conv_mol)^3 * conv_mol * Al_CEC_scale
%%%             = 10^(K_AlOH3-12) * H_umol^3 * Al_CEC_scale
%%%             = 10^(-1.2) * H_umol^3 * Al_CEC_scale   (amorphous K_AlOH3 = 10.8)
%%% Al_tot includes dissolved, exchangeable and secondary Al(OH)3.  The
%%% saturation target is inventory-limited, so the phase can precipitate or
%%% redissolve without creating Al.
al_pin = scheme_is_appelo && isstruct(bulk_pH_params) && ...
    isfield(bulk_pH_params, 'Al_solubility_pin') && ...
    isequal(bulk_pH_params.Al_solubility_pin, true);
if al_pin
    % Amorphous Al(OH)3 phase (PHREEQC.DAT / SCEPTER 'amal'): override the
    % kinetic microcrystalline default (9.66) with the equilibrium phase
    % constant verified against the SCEPTER source (scepter.f90 amal case).
    K_AlOH3 = 10.8;
end

%%% Ionic-strength-conditional constants.  The species variables stay
%%% concentrations (the charge and mass balances sum them), so every
%%% mass-action constant applied to them is converted once, here, at the
%%% step's ionic strength (ERW_Conditional_Constants, Davies): k1/g1^2,
%%% k2/g2, k_w/g1^2, the four Al hydrolysis steps and the Al(OH)3 *Kso.
%%% k_H (neutral CO2) and k_nh3 (NH4+ = NH3 + H+) need no factor.
%%% Exchange selectivities stay on the concentration basis on which SCEPTER
%%% states them and on which alpha and logKH were fitted: scepter.f90 calls
%%% get_maqads_all_v4 with maqx and pro, and coefs_v2 applies no
%%% ionic-strength term to keqiex.  Rain is dilute and keeps the
%%% thermodynamic k1 and k2.  ERW_biogeochem reports the conditioned k1, k2
%%% and k_w, so saved concentrations re-speciate exactly; the bulk-extract
%%% operator forms its own lab-temperature constants and conditions them at
%%% the extract's ionic strength.
k1_thermo = k1;
k2_thermo = k2;
if isnan(activity_I_in)
    activity_I = ERW_Incoming_Ionic_Strength(Ca_tot,Mg_tot,K_tot,Na_tot, ...
        An_tot,NO3_tot,NH4_tot,PO4_tot,f_Ca,f_Mg,f_K,f_Na,CEC_tot,H,Al, ...
        AlOH4,n*Zrw*s_pre*1000,j_only,Al_CEC_scale,conv_mol);
else
    activity_I = activity_I_in;
end
[k_cond,K_AlOH3,gamma_z] = ERW_Conditional_Constants( ...
    [k1,k2,k_w,K1,K2,K3,K4],K_AlOH3,activity_I);
k1 = k_cond(1); k2 = k_cond(2); k_w = k_cond(3);
K1 = k_cond(4); K2 = k_cond(5); K3 = k_cond(6); K4 = k_cond(7);
gamma_mono = gamma_z(1);
gamma_div = gamma_z(2);

al3_pin_fn = @(Hp) 10^(K_AlOH3-12) * max(real(Hp), 1e-12)^3 * Al_CEC_scale;
al_w_sat_fn = @(Hp) Initialize_ERW_Equilibrium(Hp, K1, K2, K3, K4, ...
    conv_mol, K_AlOH3, Al_CEC_scale);


%% RAINWATER

CO2_w_rain = k_H*CO2_atm; % [umol/l] Henry's law - dissolved CO2 in rainwater

if ~(isfinite(pH_rain) && pH_rain > 0)
    error('ERW:invalidRainPH', 'Rain pH must be finite and positive.');
end
H_rain = 10^(-pH_rain)*conv_mol; % [umol/l]

DIC_rain = CO2_w_rain+k1_thermo*CO2_w_rain/H_rain+k2_thermo*k1_thermo*CO2_w_rain/(H_rain^2); % dissolved inorganic C in rainwater [umol/l]

%% MASS BALANCE UPDATE FROM TC-BG
% update mass balances for ERW [umol]
Ca_tot_prev = Ca_tot(j_only);
Mg_tot_prev = Mg_tot(j_only);
K_tot_prev = K_tot(j_only);
Na_tot_prev = Na_tot(j_only);
An_tot_prev = An_tot(j_only);
Al_tot_prev = Al_tot(j_only);
Si_tot_prev = Si_tot(j_only);
IC_tot_prev = IC_tot(j_only);
PO4_tot_prev = PO4_tot(j_only);
NO3_tot_prev = NO3_tot(j_only);
NH4_tot_prev = NH4_tot(j_only);
pool0 = [Ca_tot_prev, Mg_tot_prev, K_tot_prev, Na_tot_prev, An_tot_prev, ...
    Al_tot_prev, IC_tot_prev, PO4_tot_prev, NO3_tot_prev, NH4_tot_prev];
inventory_names = {'Ca','Mg','K','Na','An','Al','Si','IC','PO4','NO3','NH4'};
inventory_initial = [Ca_tot_prev, Mg_tot_prev, K_tot_prev, Na_tot_prev, ...
    An_tot_prev, Al_tot_prev, Si_tot_prev, IC_tot_prev, PO4_tot_prev, ...
    NO3_tot_prev, NH4_tot_prev];
initial_bad = ~isfinite(inventory_initial) | inventory_initial < 0;
if any(initial_bad)
    first = find(initial_bad,1);
    error('ERW:PoolInventory', ...
        ['Incoming %s total inventory must be finite and nonnegative; got %.6g. ' ...
         '%d of %d inventories are invalid.'], ...
        inventory_names{first},inventory_initial(first),sum(initial_bad), ...
        numel(inventory_initial));
end
%% I_Ca / I_Mg are the upper-boundary base-cation supply (see
%% ERW_MOD_PARAM_excel).  They are a boundary condition, not weathering,
%% and their alkalinity is booked separately in dAlk_boundary below.
Ca_tot= Ca_tot(j_only)+(I_Ca+sum(min_st(:,1)'.*EW(j_only,:))+W_CaCO3(j_only)+dCa_tot(1))*dtw;
Mg_tot= Mg_tot(j_only)+(I_Mg+sum(min_st(:,2)'.*EW(j_only,:))+W_MgCO3(j_only)+dMg_tot(1))*dtw;
K_tot= K_tot(j_only)+(sum(min_st(:,3)'.*EW(j_only,:))+dK_tot(1))*dtw;
Si_tot= Si_tot(j_only)+(sum(min_st(:,6)'.*EW(j_only,:))+dSi_tot(1))*dtw;

P_weathering_release = sum(min_st(:,7)'.*EW(j_only,:))*dtw;
if ~(isscalar(P_weathering_release) && isfinite(P_weathering_release) && ...
        P_weathering_release >= 0)
    error('ERW:invalidPhosphorusRelease', ...
        'Weathered phosphorus release must be finite and nonnegative.');
end
PO4_tot= PO4_tot(j_only)+P_weathering_release+dPO4_tot(1)*dtw;
NO3_tot= NO3_tot(j_only)+dNO3_tot(1)*dtw;
NH4_tot= NH4_tot(j_only)+dNH4_tot(1)*dtw;
%% NO3 and NH4 move linearly between T&C's daily states, which the incoming
%% check holds nonnegative, so below zero they are sub-step round-off.
NO3_tot(NO3_tot < 0) = 0;
NH4_tot(NH4_tot < 0) = 0;

%% explicit mass balances
% Biogeo Leak and plant uptake
% Al leaches as the mobile species (free Al3+ and Al(OH)4-) in the drained
% water. Al, AlOH4 and Al_tot are all on the solver's Al_CEC_scale basis
% (Al_w0 seed, AlOH4 closure and the Al_tot balance below), so the mobile
% inventory (Al+AlOH4)*water_l times f_Lk_biogeo is already on Al_tot's
% basis, the same construction as Leak_Na and Leak_An. The legacy TC-SMEW
% line carried an extra factor 1e-3 on that term, which put it on the
% umol basis against the other two terms and held Al leaching about a
% thousand times below the drainage the other species see; that factor is
% removed here (POST_LEGACY_EQUATION_AUDIT.md section 3a).
if keyword_add == 0
    Leak_Na = max(0,min([Na_tot(j_only)*f_Lk_biogeo, Na_tot(j_only), Lw(j_only)*1000*Na(j_only)]));
    Leak_Al = max(0,min([(Al(j_only)+AlOH4(j_only))*n*Zrw*s*1000*f_Lk_biogeo, Al_tot(j_only), Lw(j_only)*1000*(Al(j_only)+AlOH4(j_only))]));
    Leak_An = max(0,min([An_tot(j_only)*f_Lk_biogeo, An_tot(j_only), Lw(j_only)*1000*An(j_only)]));
    Leak_DIC = max(0,min([IC_tot(j_only)*f_Lk_biogeo, IC_tot(j_only), Lw(j_only)*1000*DIC(j_only)]));

    Na_uptake = max(0,min([Na_tot(j_only)*f_T_biogeo, Na_tot(j_only), Tw(j_only)*1000*Na(j_only)]));
    An_uptake = max(0,min([An_tot(j_only)*f_T_biogeo, An_tot(j_only), Tw(j_only)*1000*An(j_only)]));
else
    Leak_Na = max(0,min([Na_tot(j_only)*f_Lk_biogeo, Na_tot(j_only)]));
    Leak_Al = max(0,min([(Al(j_only)+AlOH4(j_only))*n*Zrw*s*1000*f_Lk_biogeo, Al_tot(j_only)]));
    Leak_An = max(0,min([An_tot(j_only)*f_Lk_biogeo, An_tot(j_only)]));
    Leak_DIC = max(0,min([IC_tot(j_only)*f_Lk_biogeo, IC_tot(j_only)]));

    Na_uptake = max(0,min([Na_tot(j_only)*f_T_biogeo, Na_tot(j_only)]));
    An_uptake = max(0,min([An_tot(j_only)*f_T_biogeo, An_tot(j_only)]));
end



% explicit mass balances (within Biogeo Active Layer only)
Al_weathering = sum(min_st(:,5)'.*EW(j_only,:))*Al_CEC_scale;
Na_tot= Na_tot(j_only)+(I_Na+sum(min_st(:,4)'.*EW(j_only,:))-Leak_Na-Na_uptake)*dtw;
Al_tot = Al_tot(j_only)+(Al_weathering-Leak_Al)*dtw;
if ~(isfinite(Al_tot) && Al_tot >= 0)
    error('ERW:AlInventory', 'Al processes overdraw the total Al inventory.');
end
An_tot= An_tot(j_only)+(I_An-Leak_An-An_uptake)*dtw;
if isempty(mineral)
    mineral_C = zeros(size(min_st,1),1);
else
    mineral_C = ERW_Mineral_DIC_Stoichiometry(mineral);
end
if numel(mineral_C) ~= size(min_st,1)
    error('ERW:MineralDICSize','Mineral DIC stoichiometry does not match min_st.');
end
dIC_mineral = sum(mineral_C'.*EW(j_only,:))*dtw;
IC_tot= IC_tot(j_only)+dIC_mineral+Iw*1000*DIC_rain*dtw-ADV+ ...
    (W_CaCO3(j_only)+W_MgCO3(j_only)+r_het(j_only)+r_aut(j_only)-Leak_DIC).*dtw;
if IC_tot < 0
    error('ERW:DICReactionOverdraw', ...
        'Step %d: non-gas DIC processes overdraw the pool by %.3e umol m^{-2}.', ...
        time_step,-IC_tot);
end
pool_final = [Ca_tot, Mg_tot, K_tot, Na_tot, An_tot, Al_tot, IC_tot, PO4_tot, NO3_tot, NH4_tot];
%% Negative inventories make the carbonate and charge solve unphysical.
inventory_final = [Ca_tot, Mg_tot, K_tot, Na_tot, An_tot, Al_tot, ...
    Si_tot, IC_tot, PO4_tot, NO3_tot, NH4_tot];
final_bad = ~isfinite(inventory_final) | inventory_final < 0;
if any(final_bad)
    first = find(final_bad,1);
    error('ERW:PoolInventory', ...
        ['Updated %s total inventory must be finite and nonnegative; got %.6g. ' ...
         '%d of %d inventories are invalid.'], ...
        inventory_names{first},inventory_final(first),sum(final_bad), ...
        numel(inventory_final));
end
pool_target0 = ERW_NoSolid_Pool_Target(pool0);
pool_target_final = ERW_NoSolid_Pool_Target(pool_final);
Alk_tot = pool_target_final(1);
soil_water_l = n*Zrw*s*1000;
if ~isfinite(soil_water_l) || soil_water_l <= 0
    error('ERW:invalidState', 'Reactive soil-water volume must be positive.');
end
if ~(isscalar(CEC_tot) && isfinite(CEC_tot) && CEC_tot > 0)
    error('ERW:invalidParam', 'CEC_tot must be a positive finite scalar.');
end

%%% pH-dependent (variable) charge on the exchanger.
%%% CEC_tot is the prescribed effective site density at CEC_ref_pH.  Organic matter and
%%% oxide surfaces gain exchange sites as they deprotonate, so the density
%%% rises with pH; holding it constant leaves a rising-pH simulation with no
%%% sink for the released base cations and overstates the pH response to an
%%% amendment.  The sensitivity is declared per unit soil mass and converted
%%% here with exactly the factor that built CEC_tot from CEC_cmolc_kg, so it
%%% is invariant to the lumped biogeochemical layer depth Zbio.
[dCEC_dpH_cmolc_kg, CEC_ref_pH, cec_pH_diag] = ...
    ERW_Resolve_CEC_pH_Sensitivity(bulk_pH_params, scheme_str);
dCEC_areal = 0;
if dCEC_dpH_cmolc_kg > 0
    if ~(isstruct(bulk_pH_params) && isfield(bulk_pH_params,'CEC_cmolc_kg') && ...
            isscalar(bulk_pH_params.CEC_cmolc_kg) && ...
            isfinite(bulk_pH_params.CEC_cmolc_kg) && bulk_pH_params.CEC_cmolc_kg > 0)
        error('ERW:invalidParam', ...
            ['dCEC_dpH requires bulk_pH_params.CEC_cmolc_kg to convert the ' ...
             'mass-specific sensitivity onto the areal CEC_tot basis.']);
    end
    dCEC_areal = dCEC_dpH_cmolc_kg*(CEC_tot/bulk_pH_params.CEC_cmolc_kg);
end
cec_is_pH_dependent = dCEC_areal > 0;
if cec_is_pH_dependent
    CEC_at_H = @(Hx) ERW_CEC_At_H(Hx,conv_mol,CEC_tot,dCEC_areal,CEC_ref_pH);
else
    CEC_at_H = @(Hx) CEC_tot;
end
if ~(isscalar(Al_CEC_scale) && isfinite(Al_CEC_scale) && Al_CEC_scale > 0)
    error('ERW:invalidParam', 'Al_CEC_scale must be a positive finite scalar.');
end
if cec_is_pH_dependent
    al3_inventory_fn = @(p) ERW_Al_Equilibrium_Target(p(3),p(11),Al_tot, ...
        soil_water_l,CEC_at_H(p(3)),Al_CEC_scale,K_AlOH3,K1,K2,K3,K4);
else
    al3_inventory_fn = @(p) ERW_Al_Equilibrium_Target(p(3),p(11),Al_tot, ...
        soil_water_l,CEC_tot,Al_CEC_scale,K_AlOH3,K1,K2,K3,K4);
end

%% Alk_tot budget diagnostic - source partitioning
% Alk_tot = 2*Ca + 2*Mg + Na + K + NH4 - An - NO3 - PO4
% Track the Alk equivalent of each process group
Alk_tot_prev = pool_target0(1);
min_alk_sil = 2*min_st(:,1) + 2*min_st(:,2) + min_st(:,3) + min_st(:,4);
min_alk_carb = min_alk_sil.*(mineral_C > 0);
min_alk_sil(mineral_C > 0) = 0;
dAlk_sil = sum(min_alk_sil' .* EW(j_only,:)) * dtw;
dAlk_PO4_w = -P_weathering_release;
dAlk_sil = dAlk_sil + dAlk_PO4_w;
dAlk_carb = (sum(min_alk_carb' .* EW(j_only,:)) + ...
    2*W_CaCO3(j_only) + 2*W_MgCO3(j_only)) * dtw;
dAlk_N = (dNH4_tot(1) - dNO3_tot(1)) * dtw;
dAlk_bgc_CaMgK = (2*dCa_tot(1) + 2*dMg_tot(1) + dK_tot(1)) * dtw;
dAlk_bgc_P = -dPO4_tot(1) * dtw;
dAlk_dep = (I_Na - I_An) * dtw;
%% Boundary-supplied base cations carry 2 charges each.  Kept as its own
%% term, never folded into dAlk_sil or dAlk_carb, so that carbon-removal
%% attribution can subtract alkalinity the column was given rather than
%% weathered.
dAlk_boundary = (2*I_Ca + 2*I_Mg) * dtw;
dAlk_leach = (-Leak_Na + Leak_An) * dtw;
dAlk_uptake = (-Na_uptake + An_uptake) * dtw;
dAlk_check = dAlk_sil + dAlk_carb + dAlk_N + dAlk_bgc_CaMgK + dAlk_bgc_P + dAlk_dep + dAlk_boundary + dAlk_leach + dAlk_uptake;
dAlk_total = Alk_tot - Alk_tot_prev;
if abs(dAlk_check - dAlk_total) > max(1e-6, 1e-6*abs(dAlk_total))
    warning('ERW:AlkBudgetClosure', 'Alk_tot budget mismatch: total=%.6g, sum(parts)=%.6g, diff=%.6g', ...
        dAlk_total, dAlk_check, dAlk_check - dAlk_total);
end

% Store in ERW_biogeochem for daily accumulation
ERW_biogeochem.Alk_tot_prev = Alk_tot_prev;
ERW_biogeochem.dAlk_sil = dAlk_sil;
ERW_biogeochem.dAlk_carb = dAlk_carb;
ERW_biogeochem.dAlk_N = dAlk_N;
ERW_biogeochem.dAlk_bgc_CaMgK = dAlk_bgc_CaMgK;
ERW_biogeochem.dAlk_bgc_P = dAlk_bgc_P;
ERW_biogeochem.dAlk_dep = dAlk_dep;
ERW_biogeochem.dAlk_boundary = dAlk_boundary;
ERW_biogeochem.dAlk_leach = dAlk_leach;
ERW_biogeochem.dAlk_uptake = dAlk_uptake;
ERW_biogeochem.dAlk_check = dAlk_check;
ERW_biogeochem.dAlk_total = dAlk_total;

% Define the equations for carbonate alkalinity and ion concentrations
% Al mass balance: dissolved + exchange = target Al pool (legacy no-solid balance)
%%% Cation-exchange relations.
if scheme_is_appelo
    %% SCEPTER/Appelo (1994) H-referenced selectivity, f_H as reference:
    %%   f[X-Na] = K_Na*[Na]*(f_H/[H])*gamma, f[X2-Ca] = K_Ca*[Ca]*(f_H/[H])^2*gamma^2,
    %%   f[X3-Al] = K_Al*[Al3]*(f_H/[H])^3*gamma^3, gamma = 10^(alpha*f_H).
    %% The exchange-product relations are formulated in log space: the
    %% Al-(f_H/H)^3 coupling spans ~50 orders of magnitude between the
    %% seed and the Al-saturated root, which no linear residual plus a
    %% capped equation scale can navigate in double precision.  The log
    %% form has the identical root set (monotone transform) with O(1)
    %% conditioning.  Floors keep iterates finite (f_H ~ 1e-18 at the
    %% Al-saturated root is below every linear threshold).
    %% The floor is applied to both sides of each relation.  The left
    %% side clips the exchanger fraction at fl_app, so when the
    %% selectivity product on the right predicts a fraction below fl_app
    %% (an exchanger holding none of that cation to double precision) an
    %% unclipped right side leaves a residual ln(fl_app/f_pred) that no
    %% admissible fraction can remove, and the step is refused at the
    %% exact physical root.  That happens for Al in an unfertilised
    %% column whose pore water turns strongly alkaline (pore pH 11 at
    %% the refused step) with the Al(OH)3 inventory dissolved as
    %% Al(OH)4-: free Al3+ is inventory-limited, f_Al scales with
    %% H*f_H^3 and falls through 1e-24.  Clipping the right side at
    %% log(fl_app) makes every fraction at or below the floor satisfy
    %% the relation; above the floor max() returns the sum unchanged,
    %% so every step that converged before is bit-identical.
    fl_app = 1e-24;
    log_fl_app = log(fl_app);
    exchange_relations = @(p) [...
        log(max(p(11),fl_app))-max(log(K_Al_H3)+log(max(p(6)/Al_CEC_scale,fl_app)) ...
            +3*log(max(p(15)/p(3),fl_app))+3*appelo_alpha*log(10)*p(15),log_fl_app); ...
        log(max(p(12),fl_app))-max(log(K_Mg_H2)+log(max(p(7),fl_app)) ...
            +2*log(max(p(15)/p(3),fl_app))+2*appelo_alpha*log(10)*p(15),log_fl_app); ...
        log(max(p(13),fl_app))-max(log(K_Na_H)+log(max(p(9),fl_app)) ...
            +log(max(p(15)/p(3),fl_app))+appelo_alpha*log(10)*p(15),log_fl_app); ...
        log(max(p(14),fl_app))-max(log(K_K_H)+log(max(p(10),fl_app)) ...
            +log(max(p(15)/p(3),fl_app))+appelo_alpha*log(10)*p(15),log_fl_app); ...
        log(max(p(16),fl_app))-max(log(K_Ca_H2)+log(max(p(8),fl_app)) ...
            +2*log(max(p(15)/p(3),fl_app))+2*appelo_alpha*log(10)*p(15),log_fl_app); ...
        1-(p(15)+p(11)+p(12)+p(13)+p(14)+p(16))];
else
    exchange_relations = @(p) [...
        p(11)-(p(6)/Al_CEC_scale)*(p(16)^3/(K_Ca_Al*p(8)^3))^(1/2); ...
        p(12)-p(7)*(p(16)/(K_Ca_Mg*p(8))); ...
        p(13)-p(9)*(p(16)/(K_Ca_Na*p(8)))^(1/2); ...
        p(14)-p(10)*(p(16)/(K_Ca_K*p(8)))^(1/2); ...
        p(15)-(p(3))*(p(16)/(K_Ca_H*p(8)))^(1/2); ...
        1-(p(16)+p(11)+p(12)+p(13)+p(14)+p(15))];
end
%%% Two literal forms rather than one handle-threaded form: the constant-CEC
%%% branch is the original expression, so a run without dCEC_dpH is identical
%%% arithmetic and pays no per-residual function-call cost.  CEC enters the
%%% exchanger-charge row and the four cation mass balances, so a pH-dependent
%%% site density both raises the base charge the exchanger holds and removes
%%% those cations from solution; the Jacobian term this introduces is picked
%%% up by the finite-difference Jacobian lsqnonlin already builds.
al_charge_at = @(H,Al_w) (Al_w/Al_CEC_scale)* ...
    (3*H^4+2*K1*H^3+K1*K2*H^2-K1*K2*K3*K4)/ ...
    (H^4+K1*H^3+K1*K2*H^2+K1*K2*K3*H+K1*K2*K3*K4);
if cec_is_pH_dependent
make_equations_no_solid = @(target) @(p)[(target(1)-p(4))+(n*Zrw*s*1000)*al_charge_at(p(3),p(5))-p(1)*(n*Zrw*s*1000);
  target(2)-(D/Z_CO2*(p(2)/k_H-CO2_atm)*1000)*dtw- ...
      (p(2)*(1+k1/(p(3))+k2*k1/((p(3))^2))*s+(p(2)/k_H)*(1-s))*(n*Zrw*1000);
  (k1*p(2)/(p(3))+2*k1*k2*p(2)/((p(3))^2)-(p(3))+k_w/(p(3)))-p(1);
  p(4)-(p(12)+p(16)+p(13)+p(14))*CEC_at_H(p(3));
  (1-al_pin)*(p(5)*n*Zrw*s*1000+(p(11)/3)*CEC_at_H(p(3))*Al_CEC_scale-target(3)) + ...
  al_pin*(p(6)-al3_inventory_fn(p));
  p(6)-((p(3))^4/((p(3))^4+(p(3))^3*K1+(p(3))^2*K1*K2+(p(3))*K1*K2*K3+K1*K2*K3*K4))*p(5);
  p(7)*n*Zrw*s*1000+p(12)/2*CEC_at_H(p(3))-target(4);
  p(8)*n*Zrw*s*1000+p(16)/2*CEC_at_H(p(3))-target(5);
  p(9)*n*Zrw*s*1000+p(13)*CEC_at_H(p(3))-target(6);
  p(10)*n*Zrw*s*1000+p(14)*CEC_at_H(p(3))-target(7);
  exchange_relations(p)];
else
make_equations_no_solid = @(target) @(p)[(target(1)-p(4))+(n*Zrw*s*1000)*al_charge_at(p(3),p(5))-p(1)*(n*Zrw*s*1000);
  target(2)-(D/Z_CO2*(p(2)/k_H-CO2_atm)*1000)*dtw- ...
      (p(2)*(1+k1/(p(3))+k2*k1/((p(3))^2))*s+(p(2)/k_H)*(1-s))*(n*Zrw*1000);
  (k1*p(2)/(p(3))+2*k1*k2*p(2)/((p(3))^2)-(p(3))+k_w/(p(3)))-p(1);
  p(4)-(p(12)+p(16)+p(13)+p(14))*CEC_tot;
  (1-al_pin)*(p(5)*n*Zrw*s*1000+(p(11)/3)*CEC_tot*Al_CEC_scale-target(3)) + ...
  al_pin*(p(6)-al3_inventory_fn(p));
  p(6)-((p(3))^4/((p(3))^4+(p(3))^3*K1+(p(3))^2*K1*K2+(p(3))*K1*K2*K3+K1*K2*K3*K4))*p(5);
  p(7)*n*Zrw*s*1000+p(12)/2*CEC_tot-target(4);
  p(8)*n*Zrw*s*1000+p(16)/2*CEC_tot-target(5);
  p(9)*n*Zrw*s*1000+p(13)*CEC_tot-target(6);
  p(10)*n*Zrw*s*1000+p(14)*CEC_tot-target(7);
  exchange_relations(p)];
end
equations_no_solid = make_equations_no_solid(pool_target_final);
equation_names = {'Alk_tot','IC_tot','Alk_def','R_alk_CEC','Al_tot','Al_spec', ...
    'Mg_tot','Ca_tot','Na_tot','K_tot','f_Al','f_Mg','f_Na','f_K','f_H','CEC_sum'};

% initial guess
j_idx = min(j_only, length(H));
if j_idx < 1, j_idx = 1; end

H0 = H(j_idx);
if H0 < 1e-10 || isnan(H0) || ~isfinite(H0)
    error('ERW:InvalidHydrogenState', ...
        'Step %d: invalid hydrogen state H=%.6g umol L^{-1}.',time_step,H0);
end

%% Exchange-alkalinity cap retry.  R_alk is the alkalinity parked on the
%% CEC.  It can never exceed what the total pool can supply: when the
%% base-cation fractions carried in from yesterday would claim more than
%% Alk_tot, the porewater alkalinity (Alk_tot - R_alk) would go negative
%% and the carbonate system has no physical solution.  The normal path
%% below runs first; only if it fails the mass/physicality check do we
%% desorb the excess by scaling the base-cation fractions down to
%% (1-delta)*Alk_tot and re-solve.  The mass-balance equations then
%% re-partition the released cations into the porewater, so every atom is
%% conserved.  This is the mirror image of the seed An_tot fix on the
%% alkaline side.  Doing it as a fallback (not a pre-solve threshold)
%% keeps the common path bit-identical for the 1D/2D equivalence test.
R_alk_orig = R_alk(j_idx);
f_Ca_orig = f_Ca(j_idx); f_Mg_orig = f_Mg(j_idx);
f_K_orig  = f_K (j_idx); f_Na_orig = f_Na(j_idx);
erw_solve_attempt = 0;
while true
erw_solve_attempt = erw_solve_attempt + 1;

CO2_w0 = IC_tot/(n*Zrw*1000)*1/(s*(1+k1/H0+k2*k1/(H0^2))+(1-s)/k_H);
if isnan(CO2_w0) || ~isfinite(CO2_w0) || s < 1e-6
    CO2_w0 = max(k_H * CO2_atm, 1e-3);
end

R_alk0 = R_alk(j_idx);
CEC_seed = CEC_at_H(H0);
if ~(isfinite(CEC_seed) && CEC_seed > 0)
    error('ERW:nonpositiveCEC', ...
        'Step %d: variable CEC is nonpositive at the incoming pore-water pH.',time_step);
end
if al_pin
    Al_den0 = H0^4+H0^3*K1+H0^2*K1*K2+H0*K1*K2*K3+K1*K2*K3*K4;
    Al3_0 = ERW_Al_Equilibrium_Target(H0,f_Al(j_idx),Al_tot,soil_water_l, ...
        CEC_at_H(H0),Al_CEC_scale,K_AlOH3,K1,K2,K3,K4);
    Al_w0 = Al3_0/max(H0^4/Al_den0,eps);
else
    Al_w0 = max((Al_tot-(f_Al(j_idx)/3)*CEC_seed*Al_CEC_scale)/(n*Zrw*s*1000), 0);
end
Alk0 = (Alk_tot - R_alk0)/(n*Zrw*s*1000)+al_charge_at(H0,Al_w0);
if isnan(Alk0) || ~isfinite(Alk0) || s < 1e-6
    Alk0 = 1e3;
end
Mg0 =(Mg_tot-f_Mg(j_idx)/2*CEC_seed)/(n*Zrw*s*1000);
Na0 =(Na_tot-f_Na(j_idx)*CEC_seed)/(n*Zrw*s*1000);
Ca0 =(Ca_tot-f_Ca(j_idx)/2*CEC_seed)/(n*Zrw*s*1000);
K0 =(K_tot-f_K(j_idx)*CEC_seed)/(n*Zrw*s*1000);
conc_floor = 1e-3;
Al0 =(H0^4/(H0^4+H0^3*K1+H0^2*K1*K2+H0*K1*K2*K3+K1*K2*K3*K4))*Al_w0;
Na0   = max(0, Na0);

% Robust root selection for porewater H: scan charge balance for physical root
charge_res_fn = @(XH) (k1*CO2_w0/XH + 2*k1*k2*CO2_w0/(XH^2) - XH + k_w/XH) - ...
    ((Alk_tot-R_alk0)/soil_water_l+al_charge_at(XH,Al_w0));

% Root selection follows the previous porewater state passed in H.
H_anchor = H0;

% The carbonate charge residual is strictly monotone in positive H, so one
% root inside the physical pH interval is required.
pH_scan = linspace(1e-8,14-1e-8,281);
H_scan = 10.^(-pH_scan) * conv_mol;
%%% The 'fast' and 'mex' kernels (ERW_Chemistry_Kernel_Mode) evaluate
%%% charge_res_fn's expression at the same points in the same order and
%%% return the scan and the bisection's final variables bit for bit, so
%%% everything below, ERW_data included, is unchanged.  The mode is read
%%% where it is used, so the default path adds no variable to ERW_data.
switch ERW_Chemistry_Kernel_Mode(bulk_pH_params)
    case 'matlab'
        R_scan = nan(size(H_scan));
        for ii = 1:numel(H_scan)
            R_scan(ii) = charge_res_fn(H_scan(ii));
        end
    case 'fast'
        %% lo_r, iter, mid_H, R_mid and ii are read only by ERW_data's who().
        [R_scan,lo_H,hi_H,lo_r,iter,mid_H,R_mid] = ERW_Charge_Root_Kernel( ...
            H_scan,k1,k2,k_w,CO2_w0,Alk_tot,R_alk0,soil_water_l,Al_w0, ...
            Al_CEC_scale,K1,K2,K3,K4); %#ok<ASGLU>
        ii = numel(H_scan); %#ok<NASGU>
    case 'mex'
        [R_scan,lo_H,hi_H,lo_r,iter,mid_H,R_mid] = ERW_Charge_Root_Kernel_mex( ...
            H_scan,k1,k2,k_w,CO2_w0,Alk_tot,R_alk0,soil_water_l,Al_w0, ...
            Al_CEC_scale,K1,K2,K3,K4); %#ok<ASGLU>
        ii = numel(H_scan); %#ok<NASGU>
end
valid_scan = isfinite(R_scan);
sign_change = find(valid_scan(1:end-1) & valid_scan(2:end) & ...
    (R_scan(1:end-1) == 0 | R_scan(2:end) == 0 | sign(R_scan(1:end-1)) ~= sign(R_scan(2:end))));
sign_change = sign_change([true,diff(sign_change) > 1]);

if numel(sign_change) == 1
    if strcmp(ERW_Chemistry_Kernel_Mode(bulk_pH_params),'matlab')
    lo_H = H_scan(sign_change);
    hi_H = H_scan(sign_change+1);
    % Refine with bisection (bracketed, guaranteed convergent)
    lo_r = charge_res_fn(lo_H);
    for iter = 1:60
        mid_H = 0.5*(lo_H + hi_H);
        R_mid = charge_res_fn(mid_H);
        if ~isfinite(R_mid)
            hi_H = mid_H;
        elseif abs(R_mid) < max(1e-8, 1e-10) || abs(hi_H-lo_H)/max(mid_H,eps) < 1e-12
            break
        elseif sign(R_mid) == sign(lo_r)
            lo_H = mid_H;
            lo_r = R_mid;
        else
            hi_H = mid_H;
        end
    end
    end
    H_robust = 0.5*(lo_H + hi_H);
    H0_2 = H_robust;
else
    error('ERW:nonUniquePorePH', ...
        'Step %d: expected one pore-water charge root in 0 < pH < 14, found %d.', ...
        time_step,numel(sign_change));
end


% solution 1
x0 = [Alk0, CO2_w0, H0, R_alk0, Al_w0, Al0, Mg0, Ca0, Na0, K0, f_Al(j_only), f_Mg(j_only), f_Na(j_only), f_K(j_only), f_H(j_only), f_Ca(j_only)]';
try
    K_CEC_vec = [K_Ca_Mg, K_Ca_K, K_Ca_Na, K_Ca_Al, K_Ca_H];
    [f_Ca_est, f_Mg_est, f_K_est, f_Na_est, f_Al_est, f_H_est] = Exchange_Fractions_GT( ...
        max(x0(3), 1e-12), max(x0(8), 1e-12), ...
        x0(7), x0(10), x0(9), x0(5), ...
        K_CEC_vec, K1, K2, K3, K4, Al_CEC_scale, scheme_str, appelo_params);
    x0(11:16) = [f_Al_est; f_Mg_est; f_Na_est; f_K_est; f_H_est; f_Ca_est];
catch
    sum_f = sum(x0(11:16));
    if sum_f > 0 && abs(sum_f - 1) > 1e-10
        x0(11:16) = x0(11:16) / sum_f;
    end
end
x0_pool0 = ERW_NoSolid_Target_Seed(pool_target0, R_alk0, H0, ...
    f_Al(j_only), f_Mg(j_only), f_Na(j_only), f_K(j_only), f_H(j_only), f_Ca(j_only), ...
    soil_water_l, n, Zrw, s, CEC_tot, Al_CEC_scale, K1, K2, K3, K4, k1, k2, k_H);
equation_scale_no_solid = ERW_NoSolid_Equation_Scale(pool_target_final, x0, ...
    soil_water_l, CEC_tot, Al_CEC_scale, k_w, K_CEC, conc_floor, H_anchor, ...
    scheme_str, appelo_alpha, appelo_K, al_pin, al3_pin_fn);
equations_scaled_no_solid = @(p) equations_no_solid(p)./equation_scale_no_solid;

valid_initial = all(isfinite(x0)) && all(x0 >= 0);
if valid_initial
    try
        equations_test = equations_scaled_no_solid(x0);
        valid_initial = all(isfinite(equations_test));
    catch ME
        valid_initial = false;
    end
end

if ~valid_initial
    % warning('ERW: Invalid initial guess - x0 contains: Alk0=%g, CO2_w0=%g, H0=%g, R_alk0=%g, Al_w0=%g, Al0=%g, Mg0=%g, Ca0=%g, Na0=%g, K0=%g', ...
    %     Alk0, CO2_w0, H0, R_alk0, Al_w0, Al0, Mg0, Ca0, Na0, K0);
end

res_threshold = 1e-3;  % dimensionless ERW algebraic residual tolerance
if cec_is_pH_dependent
    f_Al_max = 1;
else
    f_Al_max = min(1, 3*Al_tot/(CEC_tot*Al_CEC_scale));
end
if ~(isfinite(f_Al_max) && f_Al_max >= 0)
    f_Al_max = 0;
end
x0_alt = real([Alk0, CO2_w0, H0_2, R_alk0, Al_w0, Al0, Mg0, Ca0, Na0, K0, f_Al(j_only), f_Mg(j_only), f_Na(j_only), f_K(j_only), f_H(j_only), f_Ca(j_only)]');
%% Put the alternate seed ON the branch it is meant to open.  H0_2 is the
%% carbonate charge-balance root, but the line above pairs it with the
%% PREVIOUS step's exchange fractions and the concentrations derived at
%% the previous H.  While the chemistry drifts those agree closely enough
%% to matter little; on an amendment day they disagree by three pH units,
%% so the "alternate" seed is not a high-pH state at all and lsqnonlin
%% simply falls back into the local minimum around x0.  Re-solve the
%% Appelo competitive partition at H0_2 against the current pool totals
%% so the seed satisfies the exchange laws, the site sum and the R_alk
%% definition at its own pH.
if scheme_is_appelo && isfinite(H0_2) && H0_2 > 0 && H0_2 ~= H0
    if al_pin
        Al3_alt = al3_pin_fn(H0_2)/Al_CEC_scale;
    else
        Al3_alt = max(Al0, 0)/Al_CEC_scale;
    end
    [f_alt, conc_alt] = ERW_Appelo_State_At_H(H0_2, ...
        [Ca_tot, Mg_tot, K_tot, Na_tot], Al3_alt, CEC_at_H(H0_2), soil_water_l, ...
        appelo_alpha, K_Na_H, K_K_H, K_Ca_H2, K_Mg_H2, K_Al_H3);
    if all(isfinite(f_alt)) && all(f_alt >= 0) && all(isfinite(conc_alt)) && all(conc_alt >= 0)
        R_alk_alt = (f_alt(1)+f_alt(2)+f_alt(3)+f_alt(4))*CEC_at_H(H0_2);
        Al_w_alt = Al3_alt*Al_CEC_scale * ...
            (H0_2^4+H0_2^3*K1+H0_2^2*K1*K2+H0_2*K1*K2*K3+K1*K2*K3*K4)/H0_2^4;
        Alk0_alt = (Alk_tot - R_alk_alt)/soil_water_l+al_charge_at(H0_2,Al_w_alt);
        x0_alt = real([Alk0_alt, CO2_w0, H0_2, R_alk_alt, Al_w_alt, ...
            Al3_alt*Al_CEC_scale, conc_alt(2), conc_alt(1), conc_alt(4), conc_alt(3), ...
            f_alt(5), f_alt(2), f_alt(4), f_alt(3), f_alt(6), f_alt(1)]');
    end
end
x0_alt = max(x0_alt, 0);

physical_check = @(p) ERW_Is_Physical_Solution(p, pool_target_final, ...
    soil_water_l, CEC_at_H(p(3)), Al_CEC_scale, al_pin);
%% The pH search box has to be able to reach the charge-balance root.
%% H_anchor is the previous step's H, a continuity heuristic that holds
%% while the chemistry drifts, and solve_erw_branch turns it into hard
%% bounds of H_anchor*10^(-+pH_tol).  That breaks on the day an amendment
%% is applied, where liming moves pore-water pH by several units within
%% one step: at the first Kantola application the previous step sits at
%% pH 5.2 while the carbonate charge balance (already bisected above into
%% H0_2) puts the root at pH 8.4, so a +-1 unit box excludes the answer
%% outright.  The solver then returns the best point inside the box,
%% which cannot conserve cation mass, and the column aborts at exactly
%% the application step every year.
%% Widen the box so it spans both the continuity anchor and the located
%% root.  This only enlarges the feasible set; which branch is accepted
%% is still decided by the physical check, the charge-residual gate and
%% the h_ok proximity test.
pH_tol_step = [];
if isfinite(H0) && H0 > 0 && isfinite(H0_2) && H0_2 > 0
    pH_tol_step = max(1.0, abs(log10(H0_2) - log10(H0)) + 0.5);
end
if scheme_is_appelo && al_pin
    carbonate = [k1,k2,k_w,k_H,conv_mol];
    appelo = [K_Na_H,K_K_H,K_Ca_H2,K_Mg_H2,K_Al_H3,appelo_alpha];
    gas_exchange = [D/Z_CO2*1000,CO2_atm,dtw];
    finite_al = struct('Al_tot',Al_tot,'Al_CEC_scale',Al_CEC_scale, ...
        'K_AlOH3',K_AlOH3,'K1',K1,'K2',K2,'K3',K3,'K4',K4);
    %%% Under 'fast' or 'mex' the same solve runs through
    %%% ERW_Appelo_Seed_Kernel, which returns this state bit for bit; it
    %%% takes al3_fn as the constants of al3_pin_fn and needs the constant
    %%% exchanger capacity, the only one Appelo admits.
    if strcmp(ERW_Chemistry_Kernel_Mode(bulk_pH_params),'matlab') || ...
            cec_is_pH_dependent
    reduced = ERW_Appelo_Mass_Conserving_State(Alk_tot,IC_tot, ...
        [Ca_tot,Mg_tot,K_tot,Na_tot],soil_water_l,n*Zrw*1000,s, ...
        carbonate,CEC_tot,appelo,gas_exchange, ...
        @(Hp) al3_pin_fn(Hp)/Al_CEC_scale,H_anchor,finite_al,CEC_at_H);
    else
    reduced = ERW_Appelo_Seed(ERW_Chemistry_Kernel_Mode(bulk_pH_params), ...
        Alk_tot,IC_tot, ...
        [Ca_tot,Mg_tot,K_tot,Na_tot],soil_water_l,n*Zrw*1000,s, ...
        carbonate,CEC_tot,appelo,gas_exchange, ...
        @(Hp) al3_pin_fn(Hp)/Al_CEC_scale,H_anchor,finite_al,CEC_at_H, ...
        K_AlOH3,Al_CEC_scale);
    end
    if reduced.is_valid
        Al_den_reduced = reduced.H^4+reduced.H^3*K1+reduced.H^2*K1*K2+ ...
            reduced.H*K1*K2*K3+K1*K2*K3*K4;
        Al_w_reduced = reduced.Al3*Al_CEC_scale*Al_den_reduced/reduced.H^4;
        x0 = [reduced.Alk,reduced.CO2_w,reduced.H,reduced.R_alk, ...
            Al_w_reduced,reduced.Al3*Al_CEC_scale, ...
            reduced.concentration(2),reduced.concentration(1), ...
            reduced.concentration(4),reduced.concentration(3), ...
            reduced.fraction(5),reduced.fraction(2),reduced.fraction(4), ...
            reduced.fraction(3),reduced.fraction(6),reduced.fraction(1)]';
        x0_alt = x0;
        valid_initial = true;
    end
end
trial = solve_erw_branch(equations_no_solid, equation_scale_no_solid, ...
    x0, x0_alt, valid_initial, fsolve_options, conc_floor, H_anchor, f_Al_max, ...
    pH_tol_step, physical_check);

sol = real(trial.sol);
fsolve_errors = trial.scaled_residual;
max_scaled_res = trial.max_scaled_res;
raw_residual = trial.raw_residual;
max_raw_res = trial.max_raw_res;

% Renormalize CEC fractions to sum to exactly 1, preserving Gaines-Thomas ratios
f_sum = sum(sol(11:16));
if abs(f_sum - 1) > 1e-14 && f_sum > 0
    sol(11:16) = sol(11:16) / f_sum;
end

% Check aqueous charge independently of the carbonate equation residual.
aqueous_charge_res_fn = @(p) p(3)+2*p(8)+2*p(7)+p(9)+p(10)+ ...
    NH4_tot/soil_water_l+al_charge_at(p(3),p(5))- ...
    (k_w/p(3)+k1*p(2)/p(3)+2*k1*k2*p(2)/p(3)^2+ ...
    (NO3_tot+PO4_tot+An_tot)/soil_water_l);
charge_raw = abs(aqueous_charge_res_fn(sol));
charge_tol = max(1e-3, 1e-6 * abs(sol(1)));
if charge_raw > charge_tol && trial.valid && max_scaled_res <= res_threshold
    % Charge balance is bad despite passing scaled check; re-solve with
    % tighter pH bounds to force the solver to the physical branch.
    trial2 = solve_erw_branch(equations_no_solid, equation_scale_no_solid, ...
        x0, x0_alt, valid_initial, fsolve_options, conc_floor, H_anchor, f_Al_max, 0.75, physical_check);
    if trial2.valid && ((trial2.physical && ~trial.physical) || ...
            (trial2.physical == trial.physical && trial2.max_scaled_res < trial.max_scaled_res))
        trial = trial2;
        sol = real(trial.sol);
        fsolve_errors = trial.scaled_residual;
        max_scaled_res = trial.max_scaled_res;
        raw_residual = trial.raw_residual;
        max_raw_res = trial.max_raw_res;
        f_sum = sum(sol(11:16));
        if abs(f_sum - 1) > 1e-14 && f_sum > 0
            sol(11:16) = sol(11:16) / f_sum;
        end
        charge_raw = abs(aqueous_charge_res_fn(sol));
    end
end

physics = ERW_Physical_Solution_Residuals(sol, pool_target_final, ...
    soil_water_l, CEC_at_H(sol(3)), Al_CEC_scale, al_pin);
if physics.is_valid
    break
end
if erw_solve_attempt == 1
    %% First attempt produced no physical solution: release the excess
    %% alkalinity parked on the CEC and re-solve from the capped seeds.
    %% Each fraction is also capped at what its total pool can hold
    %% (the Al analogue, f_Al_max, is already enforced inside the
    %% branch solver); the mass equations then re-partition the released
    %% cations into the porewater, so every atom is conserved.
    alk_cap_frac = 0.995;
    %% The cap is only defined for a positive alkalinity pool: with
    %% Alk_tot <= 0 no nonnegative R_alk can satisfy Alk0 < Alk_tot, and
    %% scaling by a negative factor would flip R_alk's sign and push the
    %% base-cation fractions negative.  Floor the cap at zero (full CEC
    %% release) so the re-solve stays physical and, if that state is
    %% still unhostable, the attempt-2 mass check reports it honestly.
    R_alk(j_idx) = min(R_alk_orig, alk_cap_frac*max(Alk_tot, 0));
    if R_alk_orig > 0
        f_cap = R_alk(j_idx)/R_alk_orig;
    else
        f_cap = 0;
    end
    f_Ca(j_idx) = min(f_Ca_orig*f_cap, 2*Ca_tot/CEC_seed);
    f_Mg(j_idx) = min(f_Mg_orig*f_cap, 2*Mg_tot/CEC_seed);
    f_K (j_idx) = min(f_K_orig*f_cap, K_tot/CEC_seed);
    f_Na(j_idx) = min(f_Na_orig*f_cap, Na_tot/CEC_seed);
    if isfinite(f_Al_max) && f_Al_max >= 0
        f_Al(j_idx) = min(f_Al(j_idx), f_Al_max);
    end
    if j_idx ~= j_only
        R_alk(j_only) = R_alk(j_idx);
        f_Ca(j_only) = f_Ca(j_idx); f_Mg(j_only) = f_Mg(j_idx);
        f_K (j_only) = f_K (j_idx); f_Na(j_only) = f_Na(j_idx);
        f_Al(j_only) = f_Al(j_idx);
    end
    fprintf('[ERW cap] step=%d R_alk %.6g -> %.6g (Alk_tot %.6g, factor %.4f); f_Ca %.4f->%.4f (max %.4f)\n', ...
        time_step, R_alk_orig, R_alk(j_idx), Alk_tot, f_cap, ...
        f_Ca_orig, f_Ca(j_idx), 2*Ca_tot/CEC_seed);
    if R_alk(j_idx) ~= R_alk_orig || f_Ca(j_idx) ~= f_Ca_orig || ...
            f_Mg(j_idx) ~= f_Mg_orig || f_K(j_idx) ~= f_K_orig || f_Na(j_idx) ~= f_Na_orig
        continue
    end
end
if isempty(sign_change)
    root_info = 'none';
else
    root_info = sprintf('pH %.2f', pH_scan(sign_change(1)));
end
fprintf('[ERW fail] step=%d s=%.4f swl=%.3g Zrw=%.3g n=%.3g CEC_tot=%.3g Al_CEC=%.3g\n', ...
    time_step, s, soil_water_l, Zrw, n, CEC_tot, Al_CEC_scale);
fprintf('[ERW fail] seeds H0=%.3g H0_2=%.3g Alk0=%.3g CO2_w0=%.3g R_alk0=%.3g root=%s\n', ...
    H0, H0_2, Alk0, CO2_w0, R_alk0, root_info);
fprintf('[ERW fail] sol=[%s]\n', sprintf('%.6g ', sol));
fprintf('[ERW fail] targets=[%s]\n', sprintf('%.6g ', pool_target_final));
mass_raw = physics.cation_mass_raw_residual;
mass_rel = physics.cation_mass_relative_residual;
error('ERW:NoPhysicalSolution', ...
    ['Step %d: no mass-conserving exchange solution (attempt %d).' ...
    '\n  Site sum residual=%.3e (tol %.0e)' ...
    '\n  Ca:  raw %+.3e umol m^{-2}, rel %.3e (tol %.0e)' ...
    '\n  Mg:  raw %+.3e umol m^{-2}, rel %.3e' ...
    '\n  K:   raw %+.3e umol m^{-2}, rel %.3e' ...
    '\n  Na:  raw %+.3e umol m^{-2}, rel %.3e' ...
    '\n  Al:  raw %+.3e umol m^{-2}, rel %.3e'], ...
    time_step, erw_solve_attempt, physics.site_sum_residual, ...
    physics.site_absolute_tolerance, ...
    mass_raw(1), mass_rel(1), physics.mass_relative_tolerance, ...
    mass_raw(2), mass_rel(2), ...
    mass_raw(3), mass_rel(3), ...
    mass_raw(4), mass_rel(4), ...
    mass_raw(5), mass_rel(5));
end

state_neg_tol = 1e-8;
bad_state = any(~isfinite(sol(:))) || sol(4) < -state_neg_tol || ...
    any(sol([2 3 5:16]) < -state_neg_tol);
[~,max_res_idx] = max(abs(fsolve_errors(:)));
dominant_residual_name = equation_names{max_res_idx};
no_converge = max_scaled_res > res_threshold || charge_raw > charge_tol || bad_state;
candidate_summary = erw_candidate_history_summary(trial.history, equation_names);
suppress_no_converge_warning = false;
if isstruct(bulk_pH_params) && isfield(bulk_pH_params, 'suppress_no_converge_warning')
    suppress_no_converge_warning = logical(bulk_pH_params.suppress_no_converge_warning);
end
if no_converge
    if ~suppress_no_converge_warning
        fprintf('[ERW diag] step=%d  s=%.4f  soil_water_l=%.3g  Al_tot=%.3g  Al_w0=%.3g  res=%.3g  H0=%.3e\n', ...
            time_step, s, soil_water_l, Al_tot, Al_w0, max_scaled_res, H0);
        warning('ERW:NoConverge', ...
            'Step %d: solver did not converge (dominant=%s, max scaled|res|=%.2e, raw dominant=%.2e, max raw|res|=%.2e).', ...
            time_step, dominant_residual_name, max_scaled_res, abs(raw_residual(max_res_idx)), max_raw_res);
        fprintf('[ERW candidate summary] step=%d %s\n', time_step, candidate_summary);
        if erw_candidate_debug_enabled()
            print_erw_candidate_history(time_step, trial.history, equation_names);
        end
    end
    error('ERW:NoConverge', ...
        'Step %d: refusing to write an unconverged ERW state (dominant=%s, max scaled|res|=%.2e).', ...
        time_step,dominant_residual_name,max_scaled_res);
end

Al_den_sol = sol(3)^4+sol(3)^3*K1+sol(3)^2*K1*K2+sol(3)*K1*K2*K3+K1*K2*K3*K4;
sol(6) = (sol(3)^4/Al_den_sol)*sol(5);
Al_w_sat = Initialize_ERW_Equilibrium(sol(3), K1, K2, K3, K4, conv_mol, K_AlOH3, Al_CEC_scale);
if isfinite(Al_w_sat) && Al_w_sat > 0
    Omega_AlOH3 = sol(5) / Al_w_sat;
else
    Omega_AlOH3 = NaN;
end
[~,AlOH3_solid,Al_inventory_residual] = ERW_Al_Equilibrium_Target(sol(3),sol(11),Al_tot, ...
    soil_water_l,CEC_at_H(sol(3)),Al_CEC_scale,K_AlOH3,K1,K2,K3,K4);
% update variables for current time step
CO2_w = sol(2);
H = sol(3);
R_alk = sol(4);
Al_w = sol(5);
Al = sol(6);
Mg = sol(7);
Ca = sol(8);
Na = sol(9);
K = sol(10);
f_Al = sol(11);
f_Mg = sol(12);
f_Na = sol(13);
f_K = sol(14);
f_H = sol(15);
f_Ca = sol(16);
f_CEC(1,:) = [f_Ca,f_Mg,f_K,f_Na,f_Al,f_H];
gamma_i = nan(1,5);
if scheme_is_appelo
    gamma_i = 10.^(appelo_alpha.*[1 1 2 2 3].*f_H);
end

% pH and C.  H is a concentration; an electrode reads the activity pH.
pH_conc = -log10(H/conv_mol);
pH = pH_conc-log10(gamma_mono);
CO2_air = CO2_w/k_H;%[umol/l]
HCO3 = k1*CO2_w/H;
CO3 = k2*k1*CO2_w/(H^2);
DIC = CO2_w+HCO3+CO3;
Alk = HCO3+2*CO3-H+k_w/H;

%Al speciation
AlOH =(H^3*K1/(H^4+H^3*K1+H^2*K1*K2+H*K1*K2*K3+K1*K2*K3*K4))*Al_w;%[umol/l]
AlOH2 =(H^2*K1*K2/(H^4+H^3*K1+H^2*K1*K2+H*K1*K2*K3+K1*K2*K3*K4))*Al_w;
AlOH3 =(H*K1*K2*K3/(H^4+H^3*K1+H^2*K1*K2+H*K1*K2*K3+K1*K2*K3*K4))*Al_w;
AlOH4 = Al_w-(Al+AlOH+AlOH2+AlOH3);

%concentrations
Si = Si_tot/(n*Zrw*s*1000);%[umol-conv/l]
An = An_tot/(n*Zrw*s*1000);%[umol_c-conv/l]

PO4 = PO4_tot/(n*Zrw*s*1000);%[umol-conv/l]
if ~(isfinite(NH4_tot) && NH4_tot >= 0)
    error('ERW:NH4Inventory', 'NH4 total inventory must be finite and nonnegative.');
end
NH4 = NH4_tot/(n*Zrw*s*1000);
NO3 = NO3_tot/(n*Zrw*s*1000);
charge_residual = H+2*Ca+2*Mg+Na+K+NH4+ ...
    (3*Al+2*AlOH+AlOH2-AlOH4)/Al_CEC_scale- ...
    (k_w/H+HCO3+2*CO3+NO3+PO4+An);

% Ionic strength [mol/L] of the solved pore water (ERW_Ionic_Strength, Al
% on its umol/L basis).  The wrapper checks it against activity_I, the
% value the constants were conditioned at.  With conditional constants
% HCO3 and CO3 are concentrations like Ca and Mg, so the carbonate ion
% activity product is gamma_div^2*Ca*CO3, with gamma_div taken at
% activity_I; ERW_Mineral_Saturation applies its single-ion coefficients
% at the same activity_I.
ionic_str = ERW_Ionic_Strength(Ca,Mg,Na,K,Al/Al_CEC_scale,H,NH4,NO3,HCO3,CO3,An,conv_mol, ...
    PO4,AlOH/Al_CEC_scale,AlOH2/Al_CEC_scale,AlOH4/Al_CEC_scale,k_w/H);
gamma_IAP = gamma_div^2;

%CO2 diff flux
Fs = D/(Z_CO2)*(CO2_air-CO2_atm)*1000;%[umol/d]
IC_tot = IC_tot-Fs*dtw;
if ~isfinite(IC_tot) || IC_tot < 0
    error('ERW:DICGasClosure', ...
        'Step %d: implicit gas exchange produced invalid IC_tot %.3e umol m^{-2}.', ...
        time_step,IC_tot);
end

% Carbonate minerals
CaCO3 = CaCO3 - W_CaCO3*dtw + CaCO3_rock;
MgCO3 = MgCO3 - W_MgCO3*dtw;

%Carbonate weathering
Omega_CaCO3 = gamma_IAP*Ca*CO3/K_CaCO3;%[-] activity basis
Omega_MgCO3 = gamma_IAP*Mg*CO3/K_MgCO3;
% update on current time step
[W_CaCO3,W_MgCO3]= carb_W(CaCO3,MgCO3,Omega_CaCO3,Omega_MgCO3,s,Zrw,r_CaCO3,r_MgCO3,tau_CaCO3,tau_MgCO3);
if W_CaCO3 > CaCO3/dtw || W_MgCO3 > MgCO3/dtw
    W_CaCO3 = min(W_CaCO3,CaCO3/dtw);
    W_MgCO3 = min(W_MgCO3,MgCO3/dtw);
end

%Silicate weathering
Omega = zeros(size(EW));
if rock_app > 0
    MM_min_row = ERW_Mineral_Molar_Mass_Row(MM_min,mineral);

    % mineral mass variation
    M_min_available = M_min(j_only,:) + rock_f_in*M_rock_add;
    M_min_loss = EW(j_only,:).*MM_min_row*dtw;
    M_min_tol = 1e-12*max(1,max(M_min_available));
    if any(M_min_loss < -M_min_tol) || any(M_min_loss > M_min_available + M_min_tol)
        error('ERW:mineralInventory', 'Primary-mineral dissolution exceeds available inventory.');
    end
    M_min(j_only,:) = M_min_available - M_min_loss;
    exhausted = M_min_loss >= M_min_available;
    M_min(j_only,exhausted) = 0;
    M_iner = M_iner + M_rock_add*(1-sum(rock_f_in));
    M_rock(j_only,:) = sum(M_min(j_only,:)) + M_iner;%[g]
    if M_rock(j_only,:) > 0
        rock_f(j_only,:) = M_min(j_only,:)./M_rock(j_only,:);
    else
        rock_f(j_only,:) = zeros(size(M_min(j_only,:)));
    end

    if M_rock_add > 0

        %diameter classes
        d(j_only,:,rock_app)= d_0; % particle diameter [m]
        delta_d(j_only,:,rock_app)=[d(j_only,1,rock_app) diff(d(j_only,:,rock_app))]; % change in particle diamter [m]
        psd(j_only,:,rock_app) = psd_perc_in.*M_rock_add(j_only)./delta_d(j_only,:,rock_app); % [g/m]

        %surface area
        lamb(j_only,:,rock_app)= a(rock_app)*d(j_only,:,rock_app).^b;% accounts for topogrpay and porosity in relation to grain diameter [-]
        SSA(j_only,:,rock_app)= 6./(d(j_only,:,rock_app).*rho_rock).*lamb(j_only,:,rock_app);% specific surface area [m2/g]

    end


    for t_appl = 1:size(d,3)

        %diameter class variation
        d_shrink = sum(Wr(j_only,:).*rock_f(j_only,:).*MM_min_row/rho_rock);%[m/d]
        d_loss = 2*d_shrink.*lamb(j_only,:,t_appl)*dtw;
        d_loss = min(d(j_only,:,t_appl),d_loss);
        d(j_only,:,t_appl) = d(j_only,:,t_appl)-d_loss;

        % initial psd and lamb
        psd_0_in = psd_0(find(psd_0(:,end,t_appl)>0,1),:,t_appl);
        lamb_0_in = lamb_0(find(lamb_0(:,end,t_appl)>0,1),:,t_appl);

        if isempty(psd_0_in)
            psd_0_in = psd(j_only,:,t_appl);
            lamb_0_in = lamb(j_only,:,t_appl);
        end

        % for any d values greater than 0
        ind_check = find(d(j_only,:,t_appl)>0);
        if ~isempty(ind_check)
            lamb(j_only,ind_check,t_appl) = a(t_appl)*d(j_only,ind_check,t_appl).^b; %[-]
            SSA(j_only,ind_check,t_appl) = 6./(d(j_only,ind_check,t_appl).*rho_rock).*lamb(j_only,ind_check,t_appl); % [m2/g]
            psd(j_only,ind_check,t_appl)= psd_0_in(ind_check).*(d(j_only,ind_check,t_appl)./d_0(ind_check)).^3.*(lamb_0_in(ind_check)./lamb(j_only,ind_check,t_appl));%[g/m]
            delta_d(j_only,:,t_appl) = [d(j_only,1,t_appl) diff(d(j_only,:,t_appl))];
        else
            delta_d(j_only,:,t_appl) = zeros(1,size(d,2));
        end
    end


    SA(j_only) = sum(SSA(j_only,:,:) .* psd(j_only,:,:) .* delta_d(j_only,:,:),'all'); % [m2]

    %%% Proton-promoted and base terms take the proton ACTIVITY, gamma_mono*H.
    for min_j = 1:length(mineral)
        if M_min(j_only,min_j) > 0 % minimum solvable [g]
            if strcmpi(mineral(min_j),'calcite')
                [Omega(j_only,min_j),Wr(j_only,min_j),EW(j_only,min_j)] = ...
                    ERW_Calcite_Kinetics(Ca,CO3,K_CaCO3,gamma_mono*H,CO2_air,conv_mol,s,diss_f, ...
                    k_diss_H_t(j_only,min_j),k_diss_w_t(j_only,min_j), ...
                    k_diss_OH_t(j_only,min_j),n_H(min_j),n_OH(min_j), ...
                    SA(j_only),rock_f(j_only,min_j),M_min(j_only,min_j), ...
                    MM_min_row(min_j),dtw,gamma_IAP);
            else
                [Omega(j_only,min_j),Omega_defined] = ERW_Mineral_Saturation(mineral(min_j), K, Ca, Mg, Na, Si, Al, PO4, H, K_sp(min_j), conv_mol, Al_CEC_scale, activity_I); %[-]
                if Omega_defined
                    dissolution_affinity = max(0,1-Omega(j_only,min_j));
                else
                    dissolution_affinity = 1; % missing aqueous products: far-from-equilibrium kinetic assumption
                end
                Wr(j_only,min_j) = s*diss_f*(k_diss_H_t(j_only,min_j)*(gamma_mono*H/conv_mol)^n_H(min_j) + k_diss_w_t(j_only,min_j) + k_diss_OH_t(j_only,min_j)*(gamma_mono*H/conv_mol)^n_OH(min_j))*dissolution_affinity; % [mol/m2 d]
                EW(j_only,min_j) = Wr(j_only,min_j)*SA(j_only)*rock_f(j_only,min_j); % [umol/d]
                EW_inventory_max = M_min(j_only,min_j)/MM_min_row(min_j)/dtw;
                if EW(j_only,min_j) > EW_inventory_max
                    EW(j_only,min_j) = EW_inventory_max;
                    area_fraction = SA(j_only)*rock_f(j_only,min_j);
                    if area_fraction > 0
                        Wr(j_only,min_j) = EW(j_only,min_j)/area_fraction;
                    else
                        Wr(j_only,min_j) = 0;
                    end
                end
            end
        else
            M_min(j_only,min_j) = 0;
            Omega(j_only,min_j) = 1;
            Wr(j_only,min_j) = 0;
            EW(j_only,min_j) = 0;
        end
    end
else
    Omega = 0;
    M_min = 0;
    Wr = 0;
    EW = 0;

    d = 0;
    delta_d = 0;
    lamb = 0;
    SSA = 0;
    psd = 0;
    SA = 0;

    M_iner = 0;
    M_rock = 0;

end

calcite_j = find(strcmpi(string(mineral),'calcite'),1);
EW_calcite = 0;
if rock_app > 0 && ~isempty(calcite_j)
    EW_calcite = EW(j_only,calcite_j);
end
[EW_calcite,W_CaCO3,calcite_limit] = ERW_Calcite_Reaction_Limit( ...
    EW_calcite,W_CaCO3,dtw,Omega_CaCO3,Ca,DIC,Alk, ...
    n*Zrw*s*1000,k1,k2,k_w,K_CaCO3,conv_mol,gamma_IAP);
% Nesquehonite has the same 1:1 carbonate stoichiometry, so the same
% closed-water bound applies to it, on the carbonate state calcite leaves.
[~,W_MgCO3] = ERW_Calcite_Reaction_Limit( ...
    0,W_MgCO3,dtw,Omega_MgCO3,Mg,calcite_limit.DIC_final,calcite_limit.Alk_final, ...
    n*Zrw*s*1000,k1,k2,k_w,K_MgCO3,conv_mol,gamma_IAP);
if rock_app > 0 && ~isempty(calcite_j)
    EW(j_only,calcite_j) = EW_calcite;
    area_fraction = SA(j_only)*rock_f(j_only,calcite_j);
    if area_fraction > 0
        Wr(j_only,calcite_j) = EW_calcite/area_fraction;
    else
        Wr(j_only,calcite_j) = 0;
    end
end



%% save variables
% save vairables at current time step
ERW_data = struct();
%%% One eval per variable costs about a third of the step.  'fast' and 'mex'
%%% read every variable in one eval and build the same struct, field for
%%% field, its own partial copy in ERW_data.ERW_data included
%%% (erw_workspace_struct), and build nothing for a caller that keeps no
%%% ERW_data (bulk_pH_params.skip_erw_data, set by both frames).
if strcmp(ERW_Chemistry_Kernel_Mode(bulk_pH_params),'matlab')
vars = who();

for ind_vars =1:length(vars)
    ERW_data.(vars{ind_vars})=eval(vars{ind_vars});
end
elseif ~(isfield(bulk_pH_params,'skip_erw_data') && ...
        isequal(bulk_pH_params.skip_erw_data,true))
    vars = who();
    erw_data_values = eval(['{' sprintf('%s;',vars{:}) '}']);
    ERW_data = erw_workspace_struct(vars,erw_data_values);
end

% could group variables by:
% this may not needed here as they should be already saved in TC
ERW_TC = struct();
ERW_TC.s = s;
ERW_TC.Tw = Tw;
ERW_TC.Lw = Lw;
ERW_TC.Tw = Tw;
ERW_TC.r_het = r_het;
ERW_TC.r_aut = r_aut;
ERW_TC.temp_soil = temp_soil;
ERW_TC.pH = pH;
%ERW_TC.n = n;

% - ERW_biogeochem - f_CEC, f_H, pH, etc.
% ERW_biogeochem = struct(); - created in main file
ERW_biogeochem.pH = pH;
ERW_biogeochem.pH_conc = pH_conc;
ERW_biogeochem.ionic_str = ionic_str;
ERW_biogeochem.activity_I = activity_I;
ERW_biogeochem.gamma_mono = gamma_mono;
%%% Implied soil-gas CO2 partial pressure [atm].  22.41 L/mol is the 0 degC
%%% molar volume, kept for consistency with CO2_atm = Ca/22.41; at soil
%%% temperature it puts the implied Henry constant about 9 % high.
ERW_biogeochem.pCO2_atm = CO2_air*22.41/conv_mol;
ERW_biogeochem.Dw = Dw;
ERW_biogeochem.Fs = Fs;
ERW_biogeochem.CO2_air = CO2_air;
ERW_biogeochem.temp_soil = temp_soil;
ERW_biogeochem.f_CEC = f_CEC;
ERW_biogeochem.Mg_tot = Mg_tot;
ERW_biogeochem.Ca_tot = Ca_tot;
ERW_biogeochem.K_tot = K_tot;
ERW_biogeochem.Na_tot = Na_tot;
ERW_biogeochem.An_tot = An_tot;
ERW_biogeochem.Al_tot = Al_tot;
ERW_biogeochem.Si_tot = Si_tot;
ERW_biogeochem.IC_tot = IC_tot;
ERW_biogeochem.CEC_tot = CEC_tot;
ERW_biogeochem.CEC_eff = CEC_at_H(H);
ERW_biogeochem.CEC_pH_diag = cec_pH_diag;
ERW_biogeochem.K_CEC = K_CEC;
ERW_biogeochem.R_alk = R_alk;
ERW_biogeochem.H = H;
ERW_biogeochem.Z_CO2 = Z_CO2;
ERW_biogeochem.DIC = DIC;
ERW_biogeochem.AlOH = AlOH;
ERW_biogeochem.AlOH2 = AlOH2;
ERW_biogeochem.AlOH3 = AlOH3;
ERW_biogeochem.AlOH4 = AlOH4;
ERW_biogeochem.Alk = Alk;
ERW_biogeochem.CO2_w = CO2_w;
ERW_biogeochem.Omega_CaCO3 = Omega_CaCO3;
ERW_biogeochem.Omega_MgCO3 = Omega_MgCO3;
ERW_biogeochem.CaCO3 = CaCO3;
ERW_biogeochem.MgCO3 = MgCO3;
ERW_biogeochem.Al_w = Al_w;
ERW_biogeochem.Al_w_sat = Al_w_sat;
ERW_biogeochem.AlOH3_solid = AlOH3_solid;
ERW_biogeochem.Al_inventory_residual = Al_inventory_residual;
ERW_biogeochem.charge_residual = charge_residual;
ERW_biogeochem.gamma_i = gamma_i;
ERW_biogeochem.Al = Al;
ERW_biogeochem.Mg = Mg;
ERW_biogeochem.Ca = Ca;
ERW_biogeochem.Na = Na;
ERW_biogeochem.K = K;
ERW_biogeochem.An = An;
ERW_biogeochem.Si = Si;
ERW_biogeochem.f_Al = f_Al;
ERW_biogeochem.f_Mg = f_Mg;
ERW_biogeochem.f_Na = f_Na;
ERW_biogeochem.f_K = f_K;
ERW_biogeochem.f_H = f_H;
ERW_biogeochem.f_Ca = f_Ca;
ERW_biogeochem.M_min = M_min;
ERW_biogeochem.I_K = I_K;
ERW_biogeochem.Alk_tot = Alk_tot;
ERW_biogeochem.CO2_w = CO2_w;
ERW_biogeochem.k_H = k_H;
ERW_biogeochem.k1  = k1;
ERW_biogeochem.k2  = k2;
ERW_biogeochem.k_w = k_w;
ERW_biogeochem.k_nh3 = k_nh3;
ERW_biogeochem.ADV = ADV;
ERW_biogeochem.DIC_rain=DIC_rain;
ERW_biogeochem.D=D;
ERW_biogeochem.CO3=CO3;
ERW_biogeochem.HCO3 = HCO3;
ERW_biogeochem.PO4 = PO4;
ERW_biogeochem.PO4_tot = PO4_tot;
ERW_biogeochem.NO3_tot = NO3_tot;
ERW_biogeochem.NH4_tot = NH4_tot;
ERW_biogeochem.NO3 = NO3;
ERW_biogeochem.NH4 = NH4;
ERW_biogeochem.Leak_Na = Leak_Na;
ERW_biogeochem.Leak_Al = Leak_Al;
ERW_biogeochem.Leak_An = Leak_An;
ERW_biogeochem.Leak_DIC = Leak_DIC;
ERW_biogeochem.An_uptake = An_uptake;
ERW_biogeochem.Na_uptake = Na_uptake;

% - ERW_weathering - e.g., Omega, W_CaCO3
% ERW_weathering = struct();  - created in main file
ERW_weathering.SSA = SSA;
ERW_weathering.Omega = Omega;
ERW_weathering.psd = psd;
ERW_weathering.delta_d = delta_d;
ERW_weathering.d = d;
ERW_weathering.EW = EW;
ERW_weathering.P_weathering_release = P_weathering_release;
ERW_weathering.Wr = Wr;
ERW_weathering.lamb = lamb;
ERW_weathering.W_CaCO3 = W_CaCO3;
% ERW_weathering.W_CaCO3_rock = W_CaCO3_rock;
ERW_weathering.W_MgCO3 = W_MgCO3;
ERW_weathering.Omega_AlOH3 = Omega_AlOH3;
ERW_weathering.a = a;
ERW_weathering.SA = SA;

% - ERW_plant
ERW_plant = struct();
%ERW_plant.conc = conc;
% ERW_plant.UP_Mg = UP_Mg;
% ERW_plant.UP_Ca = UP_Ca;
% ERW_plant.UP_K = UP_K;
% ERW_plant.UP_Si = UP_Si;
% ERW_plant.plant_H = plant_H;

% - ERW_var - e.g., other variables
ERW_var = struct();
ERW_var.errors = fsolve_errors;
ERW_var.max_scaled_res = max_scaled_res;
ERW_var.max_raw_res = max_raw_res;
ERW_var.raw_residual = raw_residual;
ERW_var.no_converge = no_converge;
ERW_var.dominant_residual_name = dominant_residual_name;
ERW_var.candidate_summary = candidate_summary;
ERW_var.res_threshold = res_threshold;
ERW_var.physics = physics;
ERW_var.M_rock = M_rock;

ERW_const = struct();
ERW_const.M_iner = M_iner;
ERW_const.bulk_pH_params = bulk_pH_params;


end


%% Sub-functions

function scale_Ktot = ERW_Ktot_Equation_Scale(K_tot, K0, soil_water_l, CEC_tot)
Ktot_scale_cap = 1e6;
terms = abs([K_tot, K0*soil_water_l, CEC_tot, 1]);
terms = terms(isfinite(terms) & terms > 0);
if isempty(terms)
    scale_Ktot = 1;
else
    scale_Ktot = min(max(terms), Ktot_scale_cap);
end
end

function scale_Altot = ERW_Altot_Equation_Scale(Al_tot, Al_w0, soil_water_l, al_pin, al3_pin_fn)
Altot_scale_cap = 1e6;
if nargin >= 4 && ~isempty(al_pin) && al_pin
    terms = abs([real(al3_pin_fn(1)), 1]);
    scale_Altot = min(max(terms(isfinite(terms) & terms > 0)), Altot_scale_cap);
    if ~(scale_Altot > 0)
        scale_Altot = 1;
    end
    return
end
terms = abs([Al_tot, Al_w0*soil_water_l, 1]);
terms = terms(isfinite(terms) & terms > 0);
if isempty(terms)
    scale_Altot = 1;
else
    scale_Altot = min(max(terms), Altot_scale_cap);
end
end

function scale_GT = ERW_GT_Equation_Scale(seed, K_CEC, conc_floor, H_anchor, ...
    Al_CEC_scale, scheme_str, app_alpha, appelo_K)
if nargin < 3 || isempty(conc_floor), conc_floor = 1e-3; end
if nargin < 4 || isempty(H_anchor), H_anchor = []; end
if nargin < 5 || isempty(Al_CEC_scale), Al_CEC_scale = 1e6; end
if nargin < 6 || isempty(scheme_str), scheme_str = 'TC_GT'; end
if nargin < 7 || isempty(app_alpha), app_alpha = 3.4; end
if nargin < 8, appelo_K = []; end
if ischar(scheme_str) && strcmpi(scheme_str, 'appelo')
    scale_GT = ERW_Appelo_Equation_Scale( ...
        seed, appelo_K, conc_floor, Al_CEC_scale, app_alpha);
    return
end
pH_tol_default = 1.0;
Ca_s = max(abs(seed(8)), 1e-12);
fCa_s = max(abs(seed(16)), 1e-12);
Ca_min = max(conc_floor, 1e-12);
fCa_max = 1.0;
H_max = [];
if ~isempty(H_anchor) && H_anchor > 0
    H_max = H_anchor * 10^pH_tol_default;
end
scale_GT = ones(6, 1);
wr11 = abs(seed(6)/Al_CEC_scale)*sqrt(fCa_s^3/max(K_CEC(4)*Ca_s^3, eps));
wr12 = abs(seed(7))*fCa_s/max(K_CEC(1)*Ca_s, eps);
wr13 = abs(seed(9))*sqrt(fCa_s/max(K_CEC(3)*Ca_s, eps));
wr14 = abs(seed(10))*sqrt(fCa_s/max(K_CEC(2)*Ca_s, eps));
wr15 = abs(seed(3))*sqrt(fCa_s/max(K_CEC(5)*Ca_s, eps));
ww11 = abs(seed(6)/Al_CEC_scale)*sqrt(fCa_max^3/max(K_CEC(4)*Ca_min^3, eps));
ww12 = abs(seed(7))*fCa_max/max(K_CEC(1)*Ca_min, eps);
ww13 = abs(seed(9))*sqrt(fCa_max/max(K_CEC(3)*Ca_min, eps));
ww14 = abs(seed(10))*sqrt(fCa_max/max(K_CEC(2)*Ca_min, eps));
if ~isempty(H_max)
    ww15 = H_max*sqrt(fCa_max/max(K_CEC(5)*Ca_min, eps));
else
    ww15 = wr15;
end
scale_GT(1) = min(max(abs([seed(11), wr11, ww11, 1])), 1e6);
scale_GT(2) = min(max(abs([seed(12), wr12, ww12, 1])), 1e6);
scale_GT(3) = min(max(abs([seed(13), wr13, ww13, 1])), 1e6);
scale_GT(4) = min(max(abs([seed(14), wr14, ww14, 1])), 1e6);
scale_GT(5) = min(max(abs([seed(15), wr15, ww15, 1])), 1e6);
scale_GT(6) = 1;
end

function scale_App = ERW_Appelo_Equation_Scale( ...
    seed, appelo_K, conc_floor, Al_CEC_scale, app_alpha)
if nargin < 3 || isempty(conc_floor), conc_floor = 1e-3; end
if nargin < 4 || isempty(Al_CEC_scale), Al_CEC_scale = 1e6; end
if nargin < 5 || isempty(app_alpha), app_alpha = 3.4; end
if numel(appelo_K) ~= 5 || any(~isfinite(appelo_K)) || any(appelo_K <= 0)
    error('ERW:invalidParam', 'Appelo equation scaling needs five positive constants.');
end
K_Al_H3 = appelo_K(1); K_Mg_H2 = appelo_K(2); K_Na_H = appelo_K(3);
K_K_H = appelo_K(4); K_Ca_H2 = appelo_K(5);
%% Log-space exchange residuals: the per-row scale is the magnitude of
%% the log residual at the seed (or 1), so iterates enter the solver
%% already O(1).  ln(K) terms and ln10*alpha*f_H are O(1..1e2); the cap
%% only guards pathologically malformed seeds.  The right side is
%% clipped at log(fl_app) exactly as in exchange_relations, so the row
%% scale is the residual the solver actually sees at the seed.
fl_app = 1e-24;
log_fl_app = log(fl_app);
logK = [log(K_Al_H3), log(K_Mg_H2), log(K_Na_H), log(K_K_H), log(K_Ca_H2)];
r_s = max(seed(15)/max(abs(seed(3)), 1e-12), fl_app);
scalogs = [
    abs(log(max(seed(11),fl_app)) - max(logK(1) + log(max(seed(6)/Al_CEC_scale,fl_app)) ...
        + 3*log(r_s) + 3*app_alpha*log(10)*max(seed(15),0), log_fl_app));
    abs(log(max(seed(12),fl_app)) - max(logK(2) + log(max(seed(7),fl_app)) ...
        + 2*log(r_s) + 2*app_alpha*log(10)*max(seed(15),0), log_fl_app));
    abs(log(max(seed(13),fl_app)) - max(logK(3) + log(max(seed(9),fl_app)) ...
        + log(r_s) + app_alpha*log(10)*max(seed(15),0), log_fl_app));
    abs(log(max(seed(14),fl_app)) - max(logK(4) + log(max(seed(10),fl_app)) ...
        + log(r_s) + app_alpha*log(10)*max(seed(15),0), log_fl_app));
    abs(log(max(seed(16),fl_app)) - max(logK(5) + log(max(seed(8),fl_app)) ...
        + 2*log(r_s) + 2*app_alpha*log(10)*max(seed(15),0), log_fl_app))];
scalogs(~isfinite(scalogs)) = 1;
scale_App = min(max([scalogs, ones(5,1)], [], 2), 1e3);
scale_App(6) = 1;
end

function target = ERW_NoSolid_Pool_Target(pool)
target = [2*pool(2)+2*pool(1)+pool(4)+pool(3)+pool(10)-pool(5)-pool(9)-pool(8), ...
    pool(7), pool(6), pool(2), pool(1), pool(4), pool(3)];
end

function seed = ERW_NoSolid_Target_Seed(target, R_alk_seed, H_seed, ...
    f_Al_seed, f_Mg_seed, f_Na_seed, f_K_seed, f_H_seed, f_Ca_seed, ...
    soil_water_l, n, Zrw, s, CEC_tot, Al_CEC_scale, K1, K2, K3, K4, k1, k2, k_H, ...
    al_pin, al_w_sat_fn)

if nargin < 23 || isempty(al_pin), al_pin = false; end
if nargin < 24 || isempty(al_w_sat_fn), al_w_sat_fn = @(Hp) NaN; end

H_seed = max(real(H_seed), 1e-10);
CO2_den = s*(1+k1/H_seed+k2*k1/(H_seed^2))+(1-s)/k_H;
if ~(isfinite(CO2_den) && CO2_den > 0)
    CO2_den = 1;
end
CO2_seed = target(2)/(n*Zrw*1000)/CO2_den;
if ~(isfinite(CO2_seed) && CO2_seed > 0)
    CO2_seed = 1e-3;
end

Alk_seed = (target(1)-R_alk_seed)/soil_water_l;
if ~isfinite(Alk_seed)
    Alk_seed = 1e3;
end

Al_w_seed = max((target(3)-(f_Al_seed/3)*CEC_tot*Al_CEC_scale)/soil_water_l, 0);
if al_pin
    Al_w_seed = max(al_w_sat_fn(H_seed), 0);
end
Mg_seed = max((target(4)-f_Mg_seed/2*CEC_tot)/soil_water_l, 0);
Ca_seed = max((target(5)-f_Ca_seed/2*CEC_tot)/soil_water_l, 0);
Na_seed = max((target(6)-f_Na_seed*CEC_tot)/soil_water_l, 0);
K_seed = max((target(7)-f_K_seed*CEC_tot)/soil_water_l, 0);

Al_den = H_seed^4+H_seed^3*K1+H_seed^2*K1*K2+H_seed*K1*K2*K3+K1*K2*K3*K4;
if ~(isfinite(Al_den) && Al_den > 0)
    Al_den = 1;
end
Al_seed = max((H_seed^4/Al_den)*Al_w_seed, 0);

seed = real([Alk_seed, CO2_seed, H_seed, R_alk_seed, Al_w_seed, Al_seed, ...
    Mg_seed, Ca_seed, Na_seed, K_seed, f_Al_seed, f_Mg_seed, f_Na_seed, ...
    f_K_seed, f_H_seed, f_Ca_seed]');
seed(~isfinite(seed)) = 1e-6;
seed = max(seed, 0);
end

function equation_scale = ERW_NoSolid_Equation_Scale(target, seed, soil_water_l, ...
    CEC_tot, Al_CEC_scale, k_w, K_CEC, conc_floor, H_anchor, scheme_str, ...
    app_alpha, appelo_K, al_pin, al3_pin_fn)
if nargin < 8 || isempty(conc_floor), conc_floor = 1e-3; end
if nargin < 9 || isempty(H_anchor), H_anchor = []; end
if nargin < 10 || isempty(scheme_str), scheme_str = 'TC_GT'; end
if nargin < 11 || isempty(app_alpha), app_alpha = []; end
if nargin < 12, appelo_K = []; end
if nargin < 13 || isempty(al_pin), al_pin = false; end
if nargin < 14 || isempty(al3_pin_fn), al3_pin_fn = @(Hp) NaN; end
H_seed = max(abs(seed(3)), eps);
equation_scale = [
    max(abs([target(1), seed(1), seed(1)*soil_water_l, 1]));
    max(abs([target(2), seed(2)*soil_water_l, 1]));
    max(abs([seed(1), seed(2), seed(3), k_w/H_seed, target(1)/soil_water_l, 1]));
    max(abs([seed(4), CEC_tot, 1])) * 100;
    ERW_Altot_Equation_Scale(target(3), seed(5), soil_water_l, al_pin, al3_pin_fn);
    max(abs([seed(6), seed(5), 1]));
    max(abs([target(4), seed(7)*soil_water_l, 1]));
    max(abs([target(5), seed(8)*soil_water_l, 1]));
    max(abs([target(6), seed(9)*soil_water_l, 1]));
    ERW_Ktot_Equation_Scale(target(7), seed(10), soil_water_l, CEC_tot);
    ERW_GT_Equation_Scale(seed, K_CEC, conc_floor, H_anchor, ...
        Al_CEC_scale, scheme_str, app_alpha, appelo_K)];
equation_scale(~isfinite(equation_scale) | equation_scale <= 0) = 1;
end

function trial = solve_erw_pool_continuation(varargin)
    trial = struct('sol', NaN(16,1), 'raw_residual', Inf(16,1), ...
        'scaled_residual', Inf(16,1), 'max_raw_res', Inf, ...
        'max_scaled_res', Inf, 'valid', false, ...
        'label', 'pool_continuation_removed');
end

function trial = solve_erw_branch(equations, equation_scale, x0, x0_alt, ...
    valid_initial, fsolve_options, conc_floor, H_anchor, f_Al_max, pH_tol, physical_check)

if nargin < 10 || isempty(pH_tol)
    pH_tol = 1.0;
end
if nargin < 11 || isempty(physical_check)
    physical_check = @(~) true;
end
if nargin < 9 || isempty(f_Al_max) || ~isfinite(f_Al_max)
    f_Al_max = 1;
end
f_Al_max = max(0, min(1, f_Al_max));
if nargin < 8 || isempty(H_anchor)
    H_anchor = [];
end

h_ok = @(Hv) isempty(H_anchor) || (Hv > 0 && abs(log10(Hv) - log10(H_anchor)) <= pH_tol);
eq_scaled = @(p) equations(p)./equation_scale;

trial = make_erw_trial(x0, equations, equation_scale, f_Al_max, 'initial_x0');
trial.physical = physical_check(trial.sol);
trial_history = trial;
res_threshold = 1e-3;

if trial.valid && trial.max_scaled_res <= res_threshold && ...
        h_ok(trial.sol(3)) && trial.physical
    trial.history = trial_history;
    return
end

lb = [-Inf, 0, 1e-12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]';
ub = [Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,1,1,1,1,1,1]';
ub(11) = f_Al_max;
if ~isempty(H_anchor)
    lb(3) = max(lb(3), H_anchor * 10^(-pH_tol));
    ub(3) = min(ub(3), H_anchor * 10^(pH_tol));
end

%% Solver stopping tolerances.  The acceptance gate for a step is the
%% scaled-residual threshold res_threshold = 1e-3 plus the physicality
%% and mass-conservation checks (raw <= 1e-3, rel <= 1e-4, site <= 1e-5),
%% so a converged scaled residual around 1e-10 is five orders of margin
%% below the gate.  Driving lsqnonlin to 1e-14/1e-12 wastes the long tail
%% of iterations where the residual stalls well below double precision;
%% under the Appelo H-referenced exchange relations that tail is the
%% dominant per-step cost.  The relaxed tolerances are still orders of
%% magnitude below the gate, and any candidate that stalls above the gate
%% falls through the unchanged x0_alt / capped-fallback / error chain.
%% The options never change, so they are built once per MATLAB process
%% (once per thread worker) instead of on every call: optimoptions costs
%% about as much as a short solve, and this runs at every 10-minute step.
persistent opts_trr_cached
if isempty(opts_trr_cached)
    opts_trr_cached = optimoptions('lsqnonlin','Display','off','Algorithm','trust-region-reflective', ...
        'FunctionTolerance', 1e-10, 'StepTolerance', 1e-10, ...
        'OptimalityTolerance', 1e-10, 'MaxIterations', 200, ...
        'MaxFunctionEvaluations', 100000);
end
opts_trr = opts_trr_cached;
% Single TR-R lsqnonlin attempt (standard RTM practice)
try
    sol_try = lsqnonlin(eq_scaled, min(max(real(x0), lb), ub), lb, ub, opts_trr);
    candidate = make_erw_trial(sol_try, equations, equation_scale, f_Al_max, 'scaled_trr');
    candidate.physical = physical_check(candidate.sol);
    trial_history(end+1) = candidate;
    trial = better_erw_trial(trial, candidate, h_ok, res_threshold);
catch
end

% Fallback: if best solution is not physical, retry with alternative initial guess
candidate_physical = isfield(trial, 'physical') && trial.physical;
candidate_valid = isfield(trial, 'valid') && trial.valid;
if ~candidate_physical && exist('x0_alt', 'var') && ~isempty(x0_alt) && all(isfinite(x0_alt))
    x0_alt_clamped = min(max(real(x0_alt), lb), ub);
    try
        sol_try2 = lsqnonlin(eq_scaled, x0_alt_clamped, lb, ub, opts_trr);
        candidate2 = make_erw_trial(sol_try2, equations, equation_scale, f_Al_max, 'scaled_trr_alt');
        candidate2.physical = physical_check(candidate2.sol);
        trial_history(end+1) = candidate2;
        trial = better_erw_trial(trial, candidate2, h_ok, res_threshold);
    catch
    end
end

trial.history = trial_history;
end



function trial = make_erw_trial(sol, equations, equation_scale, f_Al_max, label)

if nargin < 5 || isempty(label)
    label = 'candidate';
end

trial = struct();
trial.label = char(label);
trial.sol = real(sol(:));
trial.raw_residual = inf(16,1);
trial.scaled_residual = inf(16,1);
trial.max_raw_res = inf;
trial.max_scaled_res = inf;
trial.valid = false;

if numel(trial.sol) ~= 16 || any(~isfinite(trial.sol(:))) || any(abs(imag(sol(:))) > 0)
    return
end

try
    raw_residual = equations(trial.sol);
    scaled_residual = raw_residual./equation_scale;
catch
    return
end

trial.raw_residual = raw_residual;
trial.scaled_residual = scaled_residual;
trial.max_raw_res = max(abs(raw_residual(:)));
trial.max_scaled_res = max(abs(scaled_residual(:)));
state_neg_tol = 1e-8;
trial.valid = all(isfinite(raw_residual(:))) && all(isfinite(scaled_residual(:))) && ...
    trial.sol(4) >= -state_neg_tol && all(trial.sol([2 3 5:16]) >= -state_neg_tol) && ...
    trial.sol(11) <= f_Al_max + state_neg_tol;
end

function best = better_erw_trial(best, candidate, h_ok, res_threshold)

if nargin < 3 || isempty(h_ok)
    h_ok = @(~) true;
end
if nargin < 4 || isempty(res_threshold)
    res_threshold = 1e-3;
end

if ~candidate.valid
    if ~best.valid && candidate.max_scaled_res < best.max_scaled_res
        best = candidate;
    end
    return
end
if ~best.valid
    best = candidate;
    return
end

candidate_physical = isfield(candidate, 'physical') && candidate.physical;
best_physical = isfield(best, 'physical') && best.physical;
if candidate_physical && ~best_physical
    best = candidate;
    return
elseif best_physical && ~candidate_physical
    return
end

candidate_pass = candidate.max_scaled_res <= res_threshold;
best_pass = best.max_scaled_res <= res_threshold;
if candidate_pass && ~best_pass
    best = candidate;
    return
elseif best_pass && ~candidate_pass
    return
end

res_tie = max(1e-12, 1e-6*max([1, abs(best.max_scaled_res), abs(candidate.max_scaled_res)]));
res_delta = candidate.max_scaled_res - best.max_scaled_res;
if res_delta < -res_tie
    best = candidate;
elseif abs(res_delta) <= res_tie
    candidate_h_ok = candidate.valid && h_ok(candidate.sol(3));
    best_h_ok = best.valid && h_ok(best.sol(3));
    if candidate_h_ok && ~best_h_ok
        best = candidate;
    elseif candidate_h_ok == best_h_ok && candidate.max_raw_res < best.max_raw_res
        best = candidate;
    end
end
end
function tf = ERW_Is_Physical_Solution(sol, target, soil_water_l, CEC_tot, Al_CEC_scale, al_pin)

if nargin < 6 || isempty(al_pin), al_pin = false; end
physics = ERW_Physical_Solution_Residuals(sol, target, soil_water_l, CEC_tot, Al_CEC_scale, al_pin);
tf = physics.is_valid;
end

function physics = ERW_Physical_Solution_Residuals(sol, target, soil_water_l, CEC_tot, Al_CEC_scale, al_pin)

if nargin < 6 || isempty(al_pin), al_pin = false; end

mass_abs_tol = 1e-3;
mass_rel_tol = 1e-4;
site_abs_tol = 1e-5;

cation_inventory = [
    sol(8)*soil_water_l + sol(16)*CEC_tot/2;
    sol(7)*soil_water_l + sol(12)*CEC_tot/2;
    sol(10)*soil_water_l + sol(14)*CEC_tot;
    sol(9)*soil_water_l + sol(13)*CEC_tot;
    sol(5)*soil_water_l + sol(11)*CEC_tot*Al_CEC_scale/3];
cation_target = [target(5); target(4); target(7); target(6); target(3)];
mass_raw = cation_inventory - cation_target;
mass_relative = abs(mass_raw)./max(abs(cation_target), mass_abs_tol);
if al_pin
    % A negative aqueous-plus-exchange residual is retained Al(OH)3.
    % A positive residual is an inventory overdraw and remains invalid.
    mass_raw(5) = max(mass_raw(5),0);
    mass_relative(5) = abs(mass_raw(5))/max(abs(cation_target(5)),mass_abs_tol);
end
site_residual = 1 - sum(sol(11:16));

physics = struct();
physics.site_sum_residual = site_residual;
physics.cation_mass_raw_residual = mass_raw;
physics.cation_mass_relative_residual = mass_relative;
physics.mass_absolute_tolerance = mass_abs_tol;
physics.mass_relative_tolerance = mass_rel_tol;
physics.site_absolute_tolerance = site_abs_tol;
physics.is_valid = numel(sol) == 16 && all(isfinite(sol)) && ...
    isscalar(CEC_tot) && isfinite(CEC_tot) && CEC_tot > 0 && ...
    sol(3) > 0 && all(sol([2 4:16]) >= 0) && ...
    abs(site_residual) <= site_abs_tol && ...
    all(abs(mass_raw) <= mass_abs_tol | mass_relative <= mass_rel_tol);
end

function tf = erw_candidate_debug_enabled()

flag = strtrim(getenv('ERW_DEBUG_CANDIDATES'));
tf = any(strcmpi(flag, {'1','true','yes','on'}));
end

function summary = erw_candidate_history_summary(history, equation_names)
if isempty(history)
    summary = 'no candidates';
    return
end
parts = cell(1, numel(history));
for ii = 1:numel(history)
    candidate = history(ii);
    label = sprintf('candidate_%d', ii);
    if isfield(candidate, 'label') && ~isempty(candidate.label)
        label = candidate.label;
    end
    dominant = 'unknown';
    if isfield(candidate, 'scaled_residual') && numel(candidate.scaled_residual) == numel(equation_names)
        [~,idx] = max(abs(candidate.scaled_residual(:)));
        dominant = equation_names{idx};
    end
    valid = false;
    if isfield(candidate, 'valid')
        valid = logical(candidate.valid);
    end
    max_scaled = Inf;
    if isfield(candidate, 'max_scaled_res')
        max_scaled = candidate.max_scaled_res;
    end
    parts{ii} = sprintf('%s valid=%d res=%.3g dom=%s', label, valid, max_scaled, dominant);
end
summary = strjoin(parts, ' | ');
end

function print_erw_candidate_history(time_step, history, equation_names)

if isempty(history)
    return
end
fprintf('[ERW candidate diag] step=%d candidates=%d\n', time_step, numel(history));
fprintf('  label              valid H_umol_L max_scaled dominant   Alk_def_raw Ca_raw Mg_raw CEC_sum_raw max_raw\n');
for ii = 1:numel(history)
    candidate = history(ii);
    H_umol_L = NaN;
    if isfield(candidate, 'sol') && numel(candidate.sol) >= 3
        H_umol_L = candidate.sol(3);
    end
    dominant = 'n/a';
    if isfield(candidate, 'scaled_residual') && numel(candidate.scaled_residual) == numel(equation_names)
        [~, idx] = max(abs(candidate.scaled_residual(:)));
        dominant = equation_names{idx};
    end
    label = 'candidate';
    if isfield(candidate, 'label') && ~isempty(candidate.label)
        label = candidate.label;
    end
    fprintf('  %-18s %5d %.3e %.3e %-10s %.3e %.3e %.3e %.3e %.3e\n', ...
        label, candidate.valid, H_umol_L, candidate.max_scaled_res, dominant, ...
        erw_residual_at(candidate.raw_residual, 3), ...
        erw_residual_at(candidate.raw_residual, 8), ...
        erw_residual_at(candidate.raw_residual, 7), ...
        erw_residual_at(candidate.raw_residual, 16), ...
        candidate.max_raw_res);
end
end

function val = erw_residual_at(residual, idx)

val = NaN;
if numel(residual) >= idx
    val = residual(idx);
end
end
function [W_CaCO3,W_MgCO3] = carb_W(CaCO3,MgCO3,Omega_CaCO3,Omega_MgCO3,s,Zrw,r_CaCO3,r_MgCO3,tau_CaCO3,tau_MgCO3)

    % CaCO3
    if Omega_CaCO3<=1
        W_CaCO3 = s.*CaCO3.*(1-Omega_CaCO3)./tau_CaCO3; % dissolution
    else
        W_CaCO3 = r_CaCO3*Zrw.*(1-Omega_CaCO3); % precipitation
    end

    % MgCO3
    if Omega_MgCO3 <= 1
        W_MgCO3 = s.*MgCO3.*(1-Omega_MgCO3)./tau_MgCO3; % dissolution
    else
        W_MgCO3 = r_MgCO3*Zrw.*(1-Omega_MgCO3); % precipitation
    end

end

function S = erw_workspace_struct(names,values)
%%% The struct ERW_function_at_I's field-by-field loop builds from who():
%%% the fields in names' order with values' entries.  The loop meets
%%% ERW_data itself on its way, and stores the struct built so far, so
%%% that entry is the struct of the fields before it.
self = find(strcmp(names,'ERW_data'),1);
if ~isempty(self)
    partial = struct();
    for k = 1:self-1
        partial.(names{k}) = values{k};
    end
    values{self} = partial;
end
S = cell2struct(values,names,1);
end

function I = ERW_Incoming_Ionic_Strength(Ca_tot,Mg_tot,K_tot,Na_tot,An_tot, ...
    NO3_tot,NH4_tot,PO4_tot,f_Ca,f_Mg,f_K,f_Na,CEC_tot,H,Al,AlOH4, ...
    water_l,j_only,Al_CEC_scale,conv_mol)
%%% Ionic strength [mol/L] of the incoming state, from the previous step's
%%% pools, exchange fractions and water volume.  Carbonate is not solved
%%% yet, so the anion charge the explicit anions leave open is counted as
%%% monovalent bicarbonate (CO3-- is minor below pH 8.3).  It is only the
%%% first estimate: ERW_function re-solves at the solved ionic strength
%%% whenever this one is more than 1 % off.
pick = @(x) x(min(j_only,numel(x)));
if ~(isfinite(water_l) && water_l > 0)
    I = 0;
    return
end
CEC = pick(CEC_tot);
Ca = max(pick(Ca_tot)-pick(f_Ca)/2*CEC,0)/water_l;
Mg = max(pick(Mg_tot)-pick(f_Mg)/2*CEC,0)/water_l;
Kc = max(pick(K_tot)-pick(f_K)*CEC,0)/water_l;
Na = max(pick(Na_tot)-pick(f_Na)*CEC,0)/water_l;
An = max(pick(An_tot),0)/water_l;
NO3 = max(pick(NO3_tot),0)/water_l;
NH4 = max(pick(NH4_tot),0)/water_l;
PO4 = max(pick(PO4_tot),0)/water_l;
Hc = max(pick(H),0);
Al3 = max(pick(Al),0)/Al_CEC_scale;
AlOH4c = max(pick(AlOH4),0)/Al_CEC_scale;
cation_charge = Hc+2*Ca+2*Mg+Kc+Na+NH4+3*Al3;
open_anion = max(cation_charge-(An+NO3+PO4+AlOH4c),0);
I = 0.5*(Hc+4*Ca+4*Mg+Kc+Na+NH4+9*Al3+An+NO3+PO4+AlOH4c+open_anion)/conv_mol;
if ~(isfinite(I) && I >= 0)
    I = 0;
end
end
