%% ERW module

function [ERW_data, ERW_TC, ERW_biogeochem, ERW_weathering, ERW_plant,ERW_var, ERW_const] =...
    ERW_function(n, s_pre, s, Lw, Tw, Iw, ...
        Zrw, r_het, r_aut, temp_soil, K_CEC, CEC_tot, ...
        M_rock_in, k_diss_H, k_diss_w, k_diss_OH, E_H, E_w, E_OH, R, T_ref, K1, K2, K3, K4, ...
        rho_rock, b, K_CaCO3, K_MgCO3, MM_min, ...
        min_st, dtw, conv_Al, conv_mol, ...
        Ca, Mg, Na, K, An,...
        fsolve_options, Dw_0, D_0, Alk_rain, CO2_atm, CO2_air, Z_CO2, ...
        Ca_tot,Mg_tot, K_tot,Na_tot,An_tot,Al_tot,Si_tot, ...
        Al, AlOH4, Si, ...
        R_alk,H, ...
        EW, IC_tot, Wr, psd, d, lamb, SSA,...
        Fs, DIC, ...
        I_Ca, I_Mg, I_Na, I_K, I_An, I_Si, W_CaCO3, W_MgCO3, j_only, ...
        CaCO3, MgCO3,...
        f_CEC,f_Ca,f_Mg,f_K,f_Na,f_Al,f_H, ...
        r_CaCO3,r_MgCO3, tau_CaCO3, tau_MgCO3,...
        M_min, M_iner, a, psd_0, d_0, lamb_0, time_step, mineral, K_sp, diss_f,n_H,n_OH,...
        option_TC_biochem,...
        dCa_tot,dMg_tot,dK_tot,dSi_tot,...
        cov_moll, M_rock_add, psd_perc_in, rock_f_in, SSA_in, rock_app,...
        f_Lk_biogeo,f_T_biogeo,...
        PO4_tot,PO4,dPO4_tot,...
        NO3_tot,NO3,dNO3_tot,...
        NH4_tot,NH4,dNH4_tot,...
        keyword_add,...
        CaCO3_rock)

%%%% INPUT

%     pH: pH value of the solution (unitless)
%     CO2_air: Concentration of CO2 in soil air (µmol/L)
%     Alk: Alkalinity (mol/L)
%     DIC: Dissolved inorganic carbon concentration (µmol/L)
%     Fs: Advection influx of CO2 from the air into the soil (mol-conv/m²)
%     An: Concentration of anions (µmol/L)
%     CO2_w: Concentration of dissolved CO2 in water (µmol/L)
%     Omega_CaCO3: Saturation state of calcium carbonate (unitless)
%     Omega_MgCO3: Saturation state of magnesium carbonate (unitless)
%     CaCO3: Concentration of calcium carbonate (µmol/L)
%     MgCO3: Concentration of magnesium carbonate (µmol/L)
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
%     Si: Concentration of silicon (µmol/L)
%     f_Al: Fraction of aluminum
%     f_Mg: Fraction of magnesium
%     f_Na: Fraction of sodium
%     f_K: Fraction of potassium
%     f_H: Fraction of hydrogen
%     f_Ca: Fraction of calcium
%     M_min: Mass of rock (g/m2)
%     CO2_w_rain: CO2 concentration in rainwater (µmol/L)
%     H_rain: Concentration of hydrogen ions in rainwater (µmol/L)
%     DIC_rain: Dissolved inorganic carbon concentration in rainwater (µmol/L)
%     s: Soil moisture (-)
%     Tw: Tanspiration rate (m/d)
%     Lw: Leaching rate (m/d)
%     r_het: Heterotrophic respiration rate (mol/m^2/s)
%     r_aut: Autotrophic respiration rate (mol/m^2/s)
%     temp_soil: Soil temperature (°C)
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
%     R_alk: Residual-alkalinity (µeq/L)
%     H: Concentration of hydrogen ions (µmol/L)
%     Z_CO2: Concentration of dissolved CO2 in water (µmol/L)


%%% OUTPUT

%%% `ERW_data` Structure:
% `pH`: pH value of the solution (unitless)
% `CO2_air`: Concentration of CO2 in air (µmol/L)
% `HCO3`: Concentration of bicarbonate (HCO3-) ions (µmol/L)
% `CO3`: Concentration of carbonate (CO3%-) ions (µmol/L)
% `DIC`: Dissolved inorganic carbon concentration (µmol/L)
% `Fs`: Advection influx of CO2 from the air into the soil (mol-conv/m²)
% `An`: COncentration of anions (µmol/L)
% `Alk`: Alkalinity (µeq/L)
% `CO2_w`: Concentration of dissolved CO2 in water (µmol/L)
% `Omega_CaCO3`: Saturation state of calcium carbonate (unitless)
% `Omega_MgCO3`: Saturation state of magnesium carbonate (unitless)
% `CaCO3`: Concentration of calcium carbonate (µmol/L)
% `MgCO3`: Concentration of magnesium carbonate (µmol/L)
% `AlOH`: Concentration of aluminum hydroxide (µmol/L)
% `AlOH2`: Concentration of aluminum hydroxide (µmol/L)
% `AlOH3`: Concentration of aluminum hydroxide (µmol/L)
% `AlOH4`: Concentration of aluminum hydroxide (µmol/L)
% `Al_w`: Concentration of weathered aluminium (µmol/L)
% `Al`: Concentration of aluminum (µmol/L)
% `Mg`: Concentration of magnesium (µmol/L)
% `Ca`: Concentration of calcium (µmol/L)
% `Na`: Concentration of sodium (µmol/L)
% `K`: Concentration of potassium (µmol/L)
% `Si`: Concentration of silicon (µmol/L)
% `f_Al`: Fraction of aluminum
% `f_Mg`: Fraction of magnesium
% `f_Na`: Fraction of sodium
% `f_K`: Fraction of potassium
% `f_H`: Fraction of hydrogen
% `f_Ca`: Fraction of calcium
% `M_min`: Mineralization rate (mol/m^2/s)
% `CO2_w_rain`: CO2 concentration in rainwater (µmol/L)
% `H_rain`: Concentration of hydrogen ions in rainwater (µmol/L)
% `DIC_rain`: Dissolved inorganic carbon concentration in rainwater (µmol/L)

%%% `ERW_TC` Structure:
% `s`: Soil moisture (-)
% `Tw`: Tanspiration rate (m/d)
% `Lw`: Leaching rate (m/d)
% `r_het`: Heterotrophic respiration rate (mol/m^2/s)
% `r_aut`: Autotrophic respiration rate (mol/m^2/s)
% `temp_soil`: Soil temperature (°C)

%%% `ERW_biogeochem` Structure:
% `Dw`: Water diffusivity (m2/d)
% `Fs`: Advection influx of CO2 from the air into the soil (mol-conv/m²)
% `CO2_air`: Concentration of CO2 in soil air (µmol/L)
% `temp_soil`: Soil temperature (°C)
% `f_CEC`: Fraction of cation exchange capacity
% `Mg_tot`: Total magnesium concentration (µmol/L)
% `Ca_tot`: Total calcium concentration (µmol/L)
% `K_tot`: Total potassium concentration (µmol/L)
% `Na_tot`: Total sodium concentration (µmol/L)
% `An_tot`: Total anions (µeq/L)
% `Al_tot`: Total aluminum concentration (µmol/L)
% `Si_tot`: Total silicon concentration (µmol/L)
% `IC_tot`: Total inorganic carbon concentration (µmol/L)
% `R_alk`: Alkalinity (µeq/L)
% `H`: Concentration of hydrogen ions (µmol/L)
% `Z_CO2`: Concentration of dissolved CO2 in water (µmol/L)
% `DIC`: Dissolved inorganic carbon concentration (µmol/L)
% `AlOH`, `AlOH2`, `AlOH3`, `AlOH4`: Concentration of aluminum hydroxide species (µmol/L)
% `Alk`: Alkalinity (µeq/L)
% `CO2_w`: Concentration of dissolved CO2 in water (µmol/L)
% `Omega_CaCO3`: Saturation state of calcium carbonate (unitless)
% `Omega_MgCO3`: Saturation state of magnesium carbonate (unitless)
% `CaCO3`: Concentration of calcium carbonate (µmol/L)
% `MgCO3`: Concentration of magnesium carbonate (µmol/L)
% `Al_w`: Concentration of aluminum in water (µmol/L)
% `Al`: Concentration of aluminum (µmol/L)
% `Mg`: Concentration of magnesium (µmol/L)
% `Ca`: Concentration of calcium (µmol/L)
% `Na`: Concentration of sodium (µmol/L)
% `K`: Concentration of potassium (µmol/L)
% `An`: Concentration of anions (µeq/L)
% `Si`: Concentration of silicon (µmol/L)
% `f_Al`, `f_Mg`, `f_Na`, `f_K`, `f_H`, `f_Ca`: Fraction of respective elements
% `M_min`: Mass of mineral (g/mol)

%%% `ERW_weathering` Structure:
% `SSA`: Specific surface area (m^2/g)
% `Omega`: Saturation state (unitless)
% `EW`: Enhanced weathering rate (mol/m^2/s)
% `Wr`: Weathering rate (mol/m^2/s)
% `lamb`: Weathering constant (unitless)
% `W_CaCO3`: Weathering rate of calcium carbonate (mol/m^2/s)
% `W_MgCO3`: Weathering rate of magnesium carbonate (mol/m^2/s)

%%% `ERW_plant` Structure:
% `GF_mn`: Growth factor (unitless)

%%% `ERW_var` Structure:
% `errors`: Errors (unitless)

%%

% CO2 advection due to moisture variation
if s <= s_pre
    ADV = n*Zrw*1000*(s-s_pre)*CO2_atm; 
elseif s > s_pre
    ADV = n*Zrw*1000*(s-s_pre)*CO2_air;
end

%% state variables

% soil temperature constant
T_K = temp_soil+273.15;

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

% water dissociation constant
k_w = 10^(-pk_w)*conv_mol^2;

% Henry - constant for CO2 solubility
k_H = 0.83*exp(2400*(1/T_K-1/T_ref));


%CEC Gaines-Thomas constants
K_CEC_ = num2cell(K_CEC);
[K_Ca_Mg,K_Ca_K,K_Ca_Na,K_Ca_Al,K_Ca_H]= K_CEC_{:};


%% RAINWATER

CO2_w_rain = k_H*CO2_atm; % [umol/l] Henry's law - dissolved CO2 in rainwater

% solving for rainwater acidity
equations = @(p)Alk_rain-(k1*CO2_w_rain/p+2*k1*k2*CO2_w_rain/(p^2)-p+k_w/p);

[H_rain, ~, exitflag,~] = fsolve(equations,10^-6*conv_mol, fsolve_options); % [umol/l]
if exitflag<=0
    error('Equation is not solved for H_rain')
end

DIC_rain = CO2_w_rain+k1*CO2_w_rain/H_rain+k2*k1*CO2_w_rain/(H_rain^2); % dissolved inorganic C in rainwater [umol/l]

%% MASS BALANCE UPDATE FROM TC-BG
% update mass balances for ERW [umol]
Ca_tot= Ca_tot(j_only)+(sum(min_st(:,1)'.*EW(j_only,:))+W_CaCO3(j_only)+dCa_tot(1))*dtw;
Mg_tot= Mg_tot(j_only)+(sum(min_st(:,2)'.*EW(j_only,:))+W_MgCO3(j_only)+dMg_tot(1))*dtw;
K_tot= K_tot(j_only)+(sum(min_st(:,3)'.*EW(j_only,:))+dK_tot(1))*dtw;
Si_tot= Si_tot(j_only)+(sum(min_st(:,6)'.*EW(j_only,:))+dSi_tot(1))*dtw;

PO4_tot= PO4_tot(j_only)+(sum(min_st(:,7)'.*EW(j_only,:))+dPO4_tot(1))*dtw;
NO3_tot= NO3_tot(j_only)+dNO3_tot(1)*dtw;
NH4_tot= NH4_tot(j_only)+dNH4_tot(1)*dtw;

%% explicit mass balances
% Biogeo Leak and plant uptake
if keyword_add == 0
    Leak_Na = max(0,min([Na_tot(j_only)*f_Lk_biogeo, Na_tot(j_only), Lw(j_only)*1000*Na(j_only)]));
    Leak_Al = max(0,min([(Al(j_only)+AlOH4(j_only))*n*Zrw*s*1000*f_Lk_biogeo*1e-3, Al_tot(j_only), Lw(j_only)*1000*(Al(j_only)+AlOH4(j_only))]));
    Leak_An = max(0,min([An_tot(j_only)*f_Lk_biogeo, An_tot(j_only), Lw(j_only)*1000*An(j_only)]));
    Leak_DIC = max(0,min([IC_tot(j_only)*f_Lk_biogeo, IC_tot(j_only), Lw(j_only)*1000*DIC(j_only)]));

    Na_uptake = max(0,min([Na_tot(j_only)*f_T_biogeo, Na_tot(j_only), Tw(j_only)*1000*Na(j_only)]));
    An_uptake = max(0,min([An_tot(j_only)*f_T_biogeo, An_tot(j_only), Tw(j_only)*1000*An(j_only)]));
else
    Leak_Na = max(0,min([Na_tot(j_only)*f_Lk_biogeo, Na_tot(j_only)]));
    Leak_Al = max(0,min([(Al(j_only)+AlOH4(j_only))*n*Zrw*s*1000*f_Lk_biogeo*1e-3, Al_tot(j_only)]));
    Leak_An = max(0,min([An_tot(j_only)*f_Lk_biogeo, An_tot(j_only)]));
    Leak_DIC = max(0,min([IC_tot(j_only)*f_Lk_biogeo, IC_tot(j_only)]));

    Na_uptake = max(0,min([Na_tot(j_only)*f_T_biogeo, Na_tot(j_only)]));
    An_uptake = max(0,min([An_tot(j_only)*f_T_biogeo, An_tot(j_only)]));
end



% explicit mass balances (within Biogeo Active Layer only)
Na_tot= Na_tot(j_only)+(I_Na+sum(min_st(:,4)'.*EW(j_only,:))-Leak_Na-Na_uptake)*dtw;
Al_tot= Al_tot(j_only)+(sum(min_st(:,5)'.*EW(j_only,:))*conv_Al-Leak_Al)*dtw;
An_tot= An_tot(j_only)+(I_An-Leak_An-An_uptake)*dtw;
Alk_tot= 2*Mg_tot+2*Ca_tot+Na_tot+K_tot+NH4_tot-An_tot-NO3_tot-PO4_tot;
IC_tot= IC_tot(j_only)+Iw*1000*DIC_rain-ADV+(W_CaCO3(j_only)+W_MgCO3(j_only)+sum(strcmp(mineral,'calcite').*EW(j_only,:))+r_het(j_only)+r_aut(j_only)-Fs(j_only)-Leak_DIC).*dtw; %% grain-route calcite releases 1 mol C per mol dissolved


% Define the equations for alkalnity and ion concentrations directly
equations =@(p)[(Alk_tot-p(4))-p(1)*(n*Zrw*s*1000);
  IC_tot-(p(2)*(1+k1/(p(3))+k2*k1/((p(3))^2))*s+(p(2)/k_H)*(1-s))*(n*Zrw*1000);
  (k1*p(2)/(p(3))+2*k1*k2*p(2)/((p(3))^2)-(p(3))+k_w/(p(3)))-p(1);
  p(4)-(p(12)+p(16)+p(13)+p(14))*CEC_tot;
  p(5)*n*Zrw*s*1000+(p(11)/3)*CEC_tot*conv_Al-Al_tot;
  p(6)-((p(3))^4/((p(3))^4+(p(3))^3*K1+(p(3))^2*K1*K2+(p(3))*K1*K2*K3+K1*K2*K3*K4))*p(5);
  p(7)*n*Zrw*s*1000+p(12)/2*CEC_tot-Mg_tot;
  p(8)*n*Zrw*s*1000+p(16)/2*CEC_tot-Ca_tot;
  p(9)*n*Zrw*s*1000+p(13)*CEC_tot-Na_tot;
  p(10)*n*Zrw*s*1000+p(14)*CEC_tot-K_tot;
  p(11)-(p(6)/conv_Al)*(p(16)^3/(K_Ca_Al*p(8)^3))^(1/2);
  p(12)-p(7)*(p(16)/(K_Ca_Mg*p(8)));
  p(13)-p(9)*(p(16)/(K_Ca_Na*p(8)))^(1/2);
  p(14)-p(10)*(p(16)/(K_Ca_K*p(8)))^(1/2);
  p(15)-(p(3))*(p(16)/(K_Ca_H*p(8)))^(1/2);
  1-(p(16)+p(11)+p(12)+p(13)+p(14)+p(15))];

% initial guess
Alk0 =(Alk_tot-R_alk(j_only))/(n*Zrw*s*1000);
CO2_w0 = IC_tot/(n*Zrw*1000)*1/(s*(1+k1/H(j_only)+k2*k1/(H(j_only)^2))+(1-s)/k_H);
R_alk0 = R_alk(j_only);
Al_w0 =(Al_tot-(f_Al(j_only)/3)*CEC_tot*conv_Al)/(n*Zrw*s*1000);
Al0 =(H(j_only)^4/(H(j_only)^4+H(j_only)^3*K1+H(j_only)^2*K1*K2+H(j_only)*K1*K2*K3+K1*K2*K3*K4))*Al_w0;
Mg0 =(Mg_tot-f_Mg(j_only)/2*CEC_tot)/(n*Zrw*s*1000);
Na0 =(Na_tot-f_Na(j_only)*CEC_tot)/(n*Zrw*s*1000);
Ca0 =(Ca_tot-f_Ca(j_only)/2*CEC_tot)/(n*Zrw*s*1000);
K0 =(K_tot-f_K(j_only)*CEC_tot)/(n*Zrw*s*1000);
H0 = H(j_only);

eqH =@(p)(k1*CO2_w0/p+2*k1*k2*CO2_w0/(p^2)-p+k_w/p)-Alk0;
H0_2 = fsolve(eqH, H(j_only), fsolve_options);


% solution 1
x0 = [Alk0, CO2_w0, H0, R_alk0, Al_w0, Al0, Mg0, Ca0, Na0, K0, f_Al(j_only), f_Mg(j_only), f_Na(j_only), f_K(j_only), f_H(j_only), f_Ca(j_only)]';
options = optimoptions(fsolve_options, 'FunctionTolerance', 1e-12,'MaxIterations',100000);
[sol, fval] = fsolve(equations, x0, options);
fsolve_errors = fval; % residuals

% solution 2
res_threshold = 1e-1;
if any(abs(fsolve_errors(:)) > res_threshold) || any(imag(sol))
    x0 = [Alk0, CO2_w0, H0_2, R_alk0, Al_w0, Al0, Mg0, Ca0, Na0, K0, f_Al(j_only), f_Mg(j_only), f_Na(j_only), f_K(j_only), f_H(j_only), f_Ca(j_only)];
    options = optimoptions(fsolve_options, 'FunctionTolerance', 1e-14,'MaxIterations',100000);
    [sol, fval] = fsolve(equations, x0, options);

    fsolve_errors = fval; % residuals

end

function [c,ceq] = unitdisk(p)
    ceq(1) = (Alk_tot-p(4))-p(1)*(n*Zrw*s*1000);
    ceq(2) = IC_tot-(p(2)*(1+k1/(p(3))+k2*k1/((p(3))^2))*s+(p(2)/k_H)*(1-s))*(n*Zrw*1000);
    ceq(3) = (k1*p(2)/(p(3))+2*k1*k2*p(2)/((p(3))^2)-(p(3))+k_w/(p(3)))-p(1);
    ceq(4) = p(4)-(p(12)+p(16)+p(13)+p(14))*CEC_tot;
    ceq(5) = p(5)*n*Zrw*s*1000+(p(11)/3)*CEC_tot*conv_Al-Al_tot;
    ceq(6) = p(6)-((p(3))^4/((p(3))^4+(p(3))^3*K1+(p(3))^2*K1*K2+(p(3))*K1*K2*K3+K1*K2*K3*K4))*p(5);
    ceq(7) = p(7)*n*Zrw*s*1000+p(12)/2*CEC_tot-Mg_tot;
    ceq(8) = p(8)*n*Zrw*s*1000+p(16)/2*CEC_tot-Ca_tot;
    ceq(9) = p(9)*n*Zrw*s*1000+p(13)*CEC_tot-Na_tot;
    ceq(10) = p(10)*n*Zrw*s*1000+p(14)*CEC_tot-K_tot;
    ceq(11) = p(11)-(p(6)/conv_Al)*(p(16)^3/(K_Ca_Al*p(8)^3))^(1/2);
    ceq(12) = p(12)-p(7)*(p(16)/(K_Ca_Mg*p(8)));
    ceq(13) = p(13)-p(9)*(p(16)/(K_Ca_Na*p(8)))^(1/2);
    ceq(14) = p(14)-p(10)*(p(16)/(K_Ca_K*p(8)))^(1/2);
    ceq(15) = p(15)-(p(3))*(p(16)/(K_Ca_H*p(8)))^(1/2);
    ceq(16) = 1-(p(16)+p(11)+p(12)+p(13)+p(14)+p(15));
    c = [];
end

if any(abs(fsolve_errors(:)) > res_threshold) || any(sol<0) || any(imag(sol))
    
    % solution 3
    x0 = [Alk0, CO2_w0, H0, R_alk0, Al_w0, Al0, Mg0, Ca0, Na0, K0, f_Al(j_only), f_Mg(j_only), f_Na(j_only), f_K(j_only), f_H(j_only), f_Ca(j_only)]';
    lb = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]';
    ub = [Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,1,1,1,1,1,1]';
    
    options2 = optimoptions('lsqnonlin','Display','off','Algorithm','interior-point','FunctionTolerance', 1e-14,'MaxIterations',100000);
    [sol, ~, residual,~,~,lambda,~] = lsqnonlin(equations, x0, lb, ub, [],[],[],[], @unitdisk, options2);

    if any(abs(residual(:)) > res_threshold) || any(imag(sol))
        % solution 4
        x0 = [Alk0, CO2_w0, H0_2, R_alk0, Al_w0, Al0, Mg0, Ca0, Na0, K0, f_Al(j_only), f_Mg(j_only), f_Na(j_only), f_K(j_only), f_H(j_only), f_Ca(j_only)];
        [sol, ~, residual,exitflag,~,lambda,~] = lsqnonlin(equations, x0, lb, ub, [],[],[],[], @unitdisk, options2);

        if any(abs(residual(:)) > res_threshold) || any(imag(sol))
            disp("Solution not converging at step "+time_step); 
        end

    end
    
end

sol = max(1e-6,real(sol));

% update variables for current time step
Alk = sol(1);
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

% pH and C
pH =-log10(H/conv_mol);%[-]
CO2_air = CO2_w/k_H;%[umol/l] 
HCO3 = k1*CO2_w/H;
CO3 = k2*k1*CO2_w/(H^2);
DIC = CO2_w+HCO3+CO3;

%Al speciation
AlOH =(H^3*K1/(H^4+H^3*K1+H^2*K1*K2+H*K1*K2*K3+K1*K2*K3*K4))*Al_w;%[umol/l]
AlOH2 =(H^2*K1*K2/(H^4+H^3*K1+H^2*K1*K2+H*K1*K2*K3+K1*K2*K3*K4))*Al_w;
AlOH3 =(H*K1*K2*K3/(H^4+H^3*K1+H^2*K1*K2+H*K1*K2*K3+K1*K2*K3*K4))*Al_w;
AlOH4 = Al_w-(Al+AlOH+AlOH2+AlOH3);

%concentrations
Si = Si_tot/(n*Zrw*s*1000);%[umol-conv/l]
An = An_tot/(n*Zrw*s*1000);%[umol_c-conv/l]

PO4 = PO4_tot/(n*Zrw*s*1000);%[umol-conv/l]
NH4 = NH4_tot/(n*Zrw*s*1000);
NO3 = NO3_tot/(n*Zrw*s*1000);

%CO2 diff flux
Fs = D/(Z_CO2)*(CO2_air-CO2_atm)*1000;%[umol/d]

% Carbonate minerals
CaCO3 = CaCO3 - W_CaCO3*dtw + CaCO3_rock;
MgCO3 = MgCO3 - W_MgCO3*dtw;

%Carbonate weathering
Omega_CaCO3 = Ca*CO3/K_CaCO3;%[-]
Omega_MgCO3 = Mg*CO3/K_MgCO3;
% update on current time step
[W_CaCO3,W_MgCO3]= carb_W(CaCO3,MgCO3,Omega_CaCO3,Omega_MgCO3,s,Zrw,r_CaCO3,r_MgCO3,tau_CaCO3,tau_MgCO3);

%Silicate weathering
if rock_app > 0

    % mineral mass variation
    M_min(j_only,:)= M_min(j_only,:)-EW(j_only,:).*MM_min(j_only,:)*dtw + rock_f_in*M_rock_add; %[g]
    M_iner = M_iner + M_rock_add*(1-sum(rock_f_in));
    M_rock(j_only,:) = sum(M_min(j_only,:)) + M_iner;%[g]
    rock_f(j_only,:) = M_min(j_only,:)./M_rock(j_only,:);

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
        d_shrink = sum(Wr(j_only,:).*MM_min(j_only,:)/rho_rock);%[m/d]
        d(j_only,:,t_appl) = d(j_only,:,t_appl)-2*d_shrink* lamb(j_only,:,t_appl)* dtw;
        
        % for any d values less than 0
        d(d(j_only,:,t_appl)<0) = 0;

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
            SSA(j_only,ind_check,t_appl) = 6./(d(j_only,ind_check,t_appl)*rho_rock).*lamb(j_only,ind_check,t_appl); % [m2/g]
            psd(j_only,ind_check,t_appl)= psd_0_in(ind_check)*(d(j_only,ind_check,t_appl)/ d_0(ind_check))^3*(lamb_0_in(ind_check)/lamb(j_only,ind_check,t_appl));%[g/m]
            delta_d(j_only,:,t_appl) = [d(j_only,1,t_appl) diff(d(j_only,:,t_appl))];
        else
            delta_d(j_only,:,t_appl) = zeros(1,size(d,2));
        end
    end


    SA(j_only) = sum(SSA(j_only,:,:) .* psd(j_only,:,:) .* delta_d(j_only,:,:),'all'); % [m2]
    
    for min_j = 1:length(mineral)
        if M_min(j_only,min_j) > 0 % minimum solvable [g]  
            Omega(j_only,min_j) = Omega_sil(mineral(min_j), K, Ca, Mg, Na, Si, Al, PO4, CO3, H, K_sp(min_j), conv_mol); %[-]
            Wr(j_only,min_j) = s*diss_f*(k_diss_H_t(j_only,min_j)*(H/conv_mol)^n_H(min_j) + k_diss_w_t(j_only,min_j) + k_diss_OH_t(j_only,min_j)*(H/conv_mol)^n_OH(min_j))*(1-Omega(j_only,min_j)); % [mol/m2 d]
            EW(j_only,min_j) = Wr(j_only,min_j)*SA(j_only)*rock_f(j_only,min_j); % [umol/d]
        else
            M_min(j_only,min_j) = 0;
            Omega(j_only,min_j) = 1;
            Wr(j_only,min_j) = 0;
            EW(j_only,min_j) = 0;
        end
    end
else
    MM_min = 0;
    n_H = 0;
    n_OH = 0;
    K_sp = 0;
    min_st = [0 0 0 0 0 0];
    k_diss_H = 0;
    k_diss_w = 0;
    k_diss_OH = 0;
    E_H = 0;
    E_w = 0;
    E_OH = 0;
    
    k_diss_H_t = 0;
    k_diss_w_t = 0;
    k_diss_OH_t = 0;
    Omega = 0;
    M_min = 0;
    rock_f = 0;
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



%% save variables
% save vairables at current time step
ERW_data = struct();
vars = who();

for ind_vars =1:length(vars)
    ERW_data.(vars{ind_vars})=eval(vars{ind_vars});
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
ERW_weathering.Wr = Wr;
ERW_weathering.lamb = lamb;
ERW_weathering.W_CaCO3 = W_CaCO3;
% ERW_weathering.W_CaCO3_rock = W_CaCO3_rock;
ERW_weathering.W_MgCO3 = W_MgCO3;
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
ERW_var.M_rock = M_rock;

ERW_const = struct();
ERW_const.M_iner = M_iner;


end


%% Sub-functions

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

function [Omega] = Omega_sil(mineral, K, Ca, Mg, Na, Si, Al, PO4, CO3, H, K_sp, conv_mol)
    %%% Silicate saturation index (Omega)
    if strcmp(mineral, 'forsterite')
        Omega = min(1,(Mg/conv_mol).^(1/2).*(Si/conv_mol).^(1/4)/(H/conv_mol)/K_sp);
    elseif strcmp(mineral, 'wollastonite')
        Omega = min(1,(Ca/conv_mol).^(1/2).*(Si/conv_mol).^(1/2)/(H/conv_mol)/K_sp);
    elseif strcmp(mineral, 'diopside')
        % CaMgSi2O6 + (H+) -> 1/4 Ca++ + 1/4Mg++ 1/2 H2SiO3
        Omega = min(1,(Ca/conv_mol).^(1/4)*(Mg/conv_mol).^(1/4)*(Si/conv_mol).^(1/2)/(H/conv_mol)/K_sp);
    elseif strcmp(mineral, 'albite')
        % NaAlSiO3 + (H+) + 5/2 H2O -> 1/2 kaolinite + 2HSiO3 + (Na+)
        Omega = min(1,(Na/conv_mol).*(Si/conv_mol).^2/(H/conv_mol)/K_sp);
    elseif strcmp(mineral, 'anorthite')
        % 1/2 CaAl2Si2O8+ (H+) + 1/2 H2O -> 1/2 kaolinite + 1/2 (Ca++)
        Omega = min(1,(Ca/conv_mol).^(1/2)/(H/conv_mol)/K_sp);
    elseif strcmp(mineral, 'apatite')
        % 1/7 Ca5(PO4)3(OH) + (H+) -> 5/7 (Ca+2) + 3/7 (H2PO4-)
        Omega = min(1,(Ca/conv_mol).^(5/7).*(PO4/conv_mol).^(3/7)/(H/conv_mol)/K_sp);
    elseif strcmp(mineral, 'Fe_forsterite')
        % 1/4 FeMgSiO4 + (H+) -> 1/4 Fe(++) + 1/4 (Mg++) + 1/4 H4SiO4
        Omega = min(1,(Mg/conv_mol).^(1/4).*(Si/conv_mol).^(1/4)/(H/conv_mol)/K_sp);
    elseif strcmp(mineral, 'alkali_feldspar')
        % 1/4.12 xxx + (H+) -> 0.41/4.12 K(+) + 0.56/4.12 (Na+) + 0.03/4.12 Ca(++) + 1.03/4.12 Al(3+) + 2.97/4.12 H4SiO4
        Omega = 0;%min(1,(K/conv_mol).^(0.0995).*(Na/conv_mol).^(0.1359).*(Ca/conv_mol).^(0.0073).*(Al/conv_mol).^(0.2500).*(Si/conv_mol).^(0.7209)/(H/conv_mol)/K_sp);
    elseif strcmp(mineral, 'chlorite')
        % (Mg0.63Fe0.37)5Al(Si3Al1.2)O10 (OH)8
        % 1/18.45 chlorite + (H+) -> 3.15/18.45 (Mg++) + 2.2/18.45 (Al+++) + 3/18.45 H4SiO4
        Omega = min(1,(Mg/conv_mol).^(3.15/18.45).*(Al/conv_mol).^(2.2/18.45).*(Si/conv_mol).^(3/18.45)/(H/conv_mol)/K_sp);
    elseif strcmp(mineral, 'epidote')
        % Ca2(Al2.53Fe1.09)SiO4(SiO7) O(OH)
        % 1/14.86 epidote + (H+) -> 2/14.86 (Ca++) + 2.53/14.86 Al(+++) + 2/14.86 (H4SiO4)
        Omega = min(1,(Ca/conv_mol).^(2/14.86).*(Al/conv_mol).^(2.53/14.86).*(Si/conv_mol).^(2/14.86)/(H/conv_mol)/K_sp);
    elseif strcmp(mineral, 'zoisite')
        % Ca2Al3Si3O12(OH)
        % (1/13)Zoisite + H+ ⇌ 3/13 Al+3 + 2/13 Ca+2 + 3/13 H4SiO4 + 1/13 H2O
        Omega = min(1,(Ca/conv_mol).^(2/13).*(Al/conv_mol).^(3/13).*(Si/conv_mol).^(3/13)/(H/conv_mol)/K_sp);
    elseif strcmp(mineral, 'ferroactinolite')
        % Ca2(Mg0.75Fe0.25)5Si8.5O22(OH)2
        % 1/11.5 ferroactinolite + (H+) -> 2/11.5 (Ca++) + 3.75/11.5 (Mg++)+ 8.5/11.5 H4SiO4
        Omega = min(1,(Ca/conv_mol).^(2/11.5).*(Mg/conv_mol).^(3.75/11.5).*(Si/conv_mol).^(8.5/11.5)/(H/conv_mol)/K_sp);
    elseif strcmp(mineral, 'quartz')
        % 1/4 SiO2 + 1/2 (H2O) -> 1/4 H4SiO4
        Omega = min(1,(Si/conv_mol).^(1/4)/(H/conv_mol)/K_sp);
    elseif strcmp(mineral, 'actinolite')
        % 1/15 actinolite + (H+)-> 0.5/15 (Al+++) + 2/15 (Ca++)
        Omega = min(1,(Ca/conv_mol).^(2/15).*(Al/conv_mol).^(0.5/15)/(H/conv_mol)/K_sp);
    elseif strcmp(mineral, 'calcite')
        % CaCO3 -> (Ca++) + (CO3--) ; Omega = aCa*aCO3/Ksp (Plummer et al 1978; Morse and Arvidson 2003)
        Omega = min(1,(Ca/conv_mol).*(CO3/conv_mol)/K_sp);
    elseif ismember(mineral, {'labradorite', 'augite', 'alkali_feldspar'})
        Omega = 0;
    else
        error('Unknown mineral');
    end


end



