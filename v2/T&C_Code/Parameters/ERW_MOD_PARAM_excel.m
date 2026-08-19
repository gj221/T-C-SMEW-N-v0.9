%%%%%%%%%%%%%%%%%%% PARAMETERS AND INITIAL CONDITION %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ENHANCED ROCK WEATHERING PARAMETER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%

function [ERW_biogeochem,ERW_plant,ERW_weathering,ERW_const,ERW_var, ERW_biogeochem_varnames,ERW_weathering_varnames,ERW_plant_varnames,ERW_const_varnames,ERW_var_varnames,ERW_all] = ...
    ERW_MOD_PARAM_excel...
    (t_input, TT_m_input,ERW_TC_const_input, ERW_biogeochem_input,...
    excel_input,opt_ERW_preprocessing)

t = t_input;
if numel(t) > 1
    dtw_init = t(2) - t(1);
elseif numel(t) == 1
    dtw_init = t(1);
else
    error('ERW_MOD_PARAM_excel:badTime', 'ERW time vector must not be empty.');
end
if ~(isfinite(dtw_init) && dtw_init > 0)
    error('ERW_MOD_PARAM_excel:badTime', 'ERW time step must be finite and positive.');
end

% excel_input.opt_TC_biochem: 
    % 1 - w TC biochem on; 
    % 0 - w/o TC biochem input
    
% excel_input.opt_crop: 
    % 1 - w crop; 
    % 0 - w/o crop
    
% excel_input.keyword_add: background inputs of cations and anions
    % 1 = balance background losses. 
    % 0 = no addition.

%% Preallocating the variables
pH = zeros(length(t), 1);
H = zeros(length(t), 1);
D = zeros(length(t), 1);
Dw = zeros(length(t), 1);

Ca_tot = zeros(length(t), 1);
Ca = zeros(length(t), 1);
UP_Ca = zeros(length(t), 1);

Omega_CaCO3 = zeros(length(t), 1);
CaCO3 = zeros(length(t), 1);
W_CaCO3 = zeros(length(t), 1);
W_CaCO3_rock = zeros(length(t), 1);

Mg_tot = zeros(length(t), 1);
Mg = zeros(length(t), 1);
UP_Mg = zeros(length(t), 1);

Omega_MgCO3 = zeros(length(t), 1);
MgCO3 = zeros(length(t), 1);
W_MgCO3 = zeros(length(t), 1);

K_tot = zeros(length(t), 1);
UP_K = zeros(length(t), 1);

Na_tot = zeros(length(t), 1);
Na = zeros(length(t), 1);

Si = zeros(length(t), 1);
Si_tot = zeros(length(t), 1);
UP_Si = zeros(length(t), 1);

PO4 = zeros(length(t), 1);
PO4_tot = zeros(length(t), 1);

NO3 = zeros(length(t), 1);
NO3_tot = zeros(length(t), 1);

NH4 = zeros(length(t), 1);
NH4_tot = zeros(length(t), 1);

K = zeros(length(t), 1);

An = zeros(length(t), 1);
An_tot = zeros(length(t), 1);
R_alk = zeros(length(t), 1);
Alk_tot = zeros(length(t), 1);
Alk = zeros(length(t), 1);

CO2_air = zeros(length(t), 1);
IC_tot = zeros(length(t), 1);
CO2_w = zeros(length(t), 1);
HCO3 = zeros(length(t), 1);
CO3 = zeros(length(t), 1);
DIC = zeros(length(t), 1);
Fs = zeros(length(t), 1);
ADV = zeros(length(t), 1);
DIC_rain = zeros(length(t), 1);

AlOH = zeros(length(t), 1);
AlOH2 = zeros(length(t), 1);
AlOH3 = zeros(length(t), 1);
AlOH4 = zeros(length(t), 1);

Al_w = zeros(length(t), 1);
Al   = zeros(length(t), 1);
Al_tot = zeros(length(t), 1);
Al_w_sat = zeros(length(t), 1);
Omega_AlOH3 = zeros(length(t), 1);

M_rock = zeros(length(t), 1);
SA = zeros(length(t), 1);
EW = zeros(length(t), 1);

plant_H = zeros(length(t), 1);

f_Ca = zeros(length(t), 1);
f_Mg = zeros(length(t), 1);
f_K = zeros(length(t), 1);
f_Na = zeros(length(t), 1);
f_Al = zeros(length(t), 1);
f_H = zeros(length(t), 1);
f_CEC_size = 6;
f_CEC = zeros(length(t),f_CEC_size);

fsolve_errors = zeros(length(t),16);


pk1 = zeros(length(t),1);
pk2 = zeros(length(t), 1);
pk_w = zeros(length(t), 1);
pk_nh3 = zeros(length(t), 1);
k1 = zeros(length(t), 1);
k2 = zeros(length(t), 1);
k_w = zeros(length(t), 1);
k_nh3 = zeros(length(t), 1);
k_H = zeros(length(t), 1);

T_K = zeros(length(t), 1);

Leak_Na = zeros(length(t), 1);
Leak_Al = zeros(length(t), 1);
Leak_An = zeros(length(t), 1);
Leak_DIC = zeros(length(t), 1);

An_uptake = zeros(length(t), 1);
Na_uptake = zeros(length(t), 1);

fsolve_options = optimoptions('fsolve','Display','off');

%% input from TC

% constants
pH_in_raw = ERW_TC_const_input.PH_in;
if ~isfinite(pH_in_raw) || pH_in_raw <= 0
    error('Initial pH must be finite and positive.');
end
pH_in = pH_in_raw;
Zbio = ERW_TC_const_input.Zbio;
n = ERW_TC_const_input.n;
rsd_bio = ERW_TC_const_input.rsd_bio;
conv_mol = ERW_TC_const_input.conv_mol;
conv_Al = ERW_TC_const_input.conv_Al;
Porg = NaN;
if isfield(ERW_TC_const_input,'Porg'); Porg = ERW_TC_const_input.Porg; end
Pcla = NaN;
if isfield(ERW_TC_const_input,'Pcla'); Pcla = ERW_TC_const_input.Pcla; end
Al_CEC_scale = conv_Al;
M_rock_add = ERW_TC_const_input.M_rock_add;

% varying varaibles
temp_soil = TT_m_input.temp_soil;
CO2_atm = TT_m_input.CO2_atm;
Tw = TT_m_input.Tw;
Lw = TT_m_input.Lw;
s = TT_m_input.s;  

%% input from parameter excel
par_names = excel_input.Properties.VariableNames;
for ind_vars = 3:length(par_names)
    val = excel_input.(par_names{ind_vars});
    var_name = par_names{ind_vars};

    %%% Workspace-eval safety: reject any xlsx column name that is not
    %%% a valid MATLAB identifier before the eval below runs.  This
    %%% closes the "xlsx controls the LHS of an eval" injection vector
    %%% that the original loop had; the RHS (val) is still passed
    %%% through eval, so the user-controlled xlsx can still inject
    %%% code on the RHS, but it can no longer control the LHS or the
    %%% shape of the resulting workspace.
    if ~isvarname(var_name)
        error('ERW_MOD_PARAM_excel:badVarname', ...
            'xlsx column name %s is not a valid MATLAB variable identifier', var_name);
    end

    if iscell(val)
        val = val{1};
    end
    if strcmp(var_name, 'mineral')
        val = parse_mineral_input(val);
    elseif any(strcmp(var_name, {'f_CEC_in','rock_f_in','d_in','psd_in'}))
        val = parse_numeric_input(val);
    end

    eval([var_name,'=val;'])
end

if isfield(ERW_TC_const_input,'rain_pH_vwm') && ~isempty(ERW_TC_const_input.rain_pH_vwm)
    rain_pH_vals = ERW_TC_const_input.rain_pH_vwm(:);
    rain_pH_idx = find(isfinite(rain_pH_vals) & rain_pH_vals > 0, 1, 'first');
    if ~isempty(rain_pH_idx)
        pH_rain = rain_pH_vals(rain_pH_idx);
    end
end
if ~exist('pH_rain','var') || ~(isfinite(pH_rain) && pH_rain > 0)
    % Legacy Kantola workbook carries Alk_rain (umol/l) but no rain-pH column:
    % fall back to the canonical equilibrium-rain pH 5.6 (CO2-saturated rain).
    pH_rain = 5.6;
end

% Derive Na_in from rain steady state when not explicitly provided or zero.
% Na_rain_ss = I_Na_rain [umol/m2/d] / (1000 * mean(Tw+Lw) [m/d -> L/m2/d]).
% Explicit non-negative Na_in from Excel takes precedence.
if ~exist('Na_in','var') || ~(isfinite(Na_in) && Na_in > 0)
    Na_rain_src = zeros(size(Tw));
    if isfield(ERW_TC_const_input, 'I_Na_source') && ~isempty(ERW_TC_const_input.I_Na_source)
        raw = ERW_TC_const_input.I_Na_source(:);
        Na_rain_src(1:min(end,numel(raw))) = raw(1:min(numel(raw),end));
    elseif isfield(ERW_TC_const_input, 'rain_Na_flux_umol_m2_d')
        raw = ERW_TC_const_input.rain_Na_flux_umol_m2_d(:);
        Na_rain_src(1:min(end,numel(raw))) = raw(1:min(numel(raw),end));
    end
    Na_rain_src(~isfinite(Na_rain_src) | Na_rain_src < 0) = 0;
    I_Na_rain_mean = mean(Na_rain_src, 'omitnan');
    mean_water_loss_Lperm2d = ERW_Seed_Water_Loss(ERW_TC_const_input, TT_m_input);
    if I_Na_rain_mean > 0
        Na_rain_ss = max(0, I_Na_rain_mean / mean_water_loss_Lperm2d);
    else
        Na_rain_ss = 0;
    end
    Na_in = Na_rain_ss;
end
Na_in = max(0, Na_in);

if ~exist('mineral','var') || isempty(mineral) || ...
        (isnumeric(mineral) && all(isnan(mineral(:))))
    rock_f_in = NaN;
end
calcite_carbonate_input = false;
if active_rock_input(rock_f_in) && numel(mineral) == 1 && ...
        strcmpi(mineral(1), 'calcite') && M_rock_in > 0
    CaCO3_in = CaCO3_in + M_rock_in/100.087*conv_mol;  %% [umol m-2]
    rock_f_in = NaN;
    calcite_carbonate_input = true;
end


%% soil properties
Zrw = Zbio/1000; % Depth of the active Biogeochemistry zone which weathering takes place in [m]
rho_bulk = rsd_bio*1e3; % soil dry mass bulk density [g/ m3]

%% initial conditions

% CEC 
[CEC_tot_in, cec_diag] = Resolve_CEC(CEC_tot_in, Pcla, Porg);
CEC_tot = CEC_tot_in*1e-5*rho_bulk*Zrw*conv_mol*ones(length(t), 1); % [umol_c] 

% pH
pH(1)= pH_in;

% Carbonate minerals (considered as an additional pool)
CaCO3(1) = CaCO3_in; %[umol-conv]
MgCO3(1) = MgCO3_in;


%%% CEC constants (Gaines-Thomas)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Soil-type selectivity values
% Base values 0-30 cm - https://edepot.wur.nl/31605
% Coefficient estimates can be improved with soil-water coupled measurements
if ismember(soil, {'sand', 'loamy sand', 'sandy loam'})
    K_Ca_Mg = 10^(0.56); % [-]
    K_Ca_K = 10^(-1.16) * conv_mol; % [conc]
    K_Ca_Na = 10^(0.75) * conv_mol; % [conc]
    K_Ca_H = 10^(-5) * conv_mol; % [conc] works with: 0.5 1e-9 *conv_mol, tePas et al
    K_Ca_Al = 10^(-1.7) / conv_mol; % [conc^-1] heterovalent (make it higher to favor Ca adsorbed)
elseif ismember(soil, {'loam', 'silty loam', 'silt'})
    K_Ca_Mg = 10^(0.1); % [-]
    K_Ca_K = 10^(-2) * conv_mol;
    K_Ca_Na = 10^(0.38) * conv_mol;
    K_Ca_H = 10^(-5.4) * conv_mol;
    K_Ca_Al = 10^(-0.86) / conv_mol;
elseif ismember(soil, {'clay', 'clay loam', 'silty clay'})
    K_Ca_Mg = 10^(0.39); % [-]
    K_Ca_K = 10^(-2.42) * conv_mol;
    K_Ca_Na = 10^(0.774) * conv_mol;
    K_Ca_H = 10^(-6.67) * conv_mol;
    K_Ca_Al = 10^(-0.2) / conv_mol;
else
    error('Unknown soil type');
end


%%% Aluminium speciation constants (pag. 398 Weil and Brady)
pK1 = 5; % acid dissociation constant [-]
pK2 = 5.1; 
pK3 = 6.7;
pK4 = 6.2;
K1 = 10^(-pK1)*conv_mol; % speciation constant [micromol]
K2 = 10^(-pK2)*conv_mol;
K3 = 10^(-pK3)*conv_mol;
K4 = 10^(-pK4)*conv_mol;
K_Al = [K1, K2, K3, K4]; % aluminium speciation constants
% log *Kso for Al(OH)3(s) + 3H+ = Al3+ + 3H2O at 25 C (mol/L convention).
% End-members: crystalline gibbsite 8.04 (Lindsay 1979, Chemical Equilibria in
% Soils, p. 51; cf. 8.11 Nordstrom & May 1996, Rev. Mineral. 31),
% amorphous Al(OH)3 10.8 (Lindsay 1979).
% 9.66 is representative of microcrystalline/proto-gibbsite — the kinetically
% active Al-hydroxide phase in agricultural soils receiving silicate amendments.
K_AlOH3 = 9.66;

% hydrogen ion conc
H(1) = 10^(-pH(1))*conv_mol; % [umol/l]

% molar masses [g/mol]
MM_Mg = 24/conv_mol; % magnesium
MM_Ca = 40/conv_mol; % calcium
MM_Na = 23/conv_mol; % sodium
MM_K = 39/conv_mol; % potassium
MM_Si = 28/conv_mol; % silicon 
MM_C = 12/conv_mol; % carbon
MM_Anions = 62/conv_mol; % average of associated anions
MM_Al = 27/conv_mol; % aluminium


if excel_input.opt_TC_biochem == 1
    Mg_in = ERW_biogeochem_input.Mg(1);
    Si_in = ERW_biogeochem_input.Si(1);
    Ca_in = ERW_biogeochem_input.Ca(1);
    K_in = ERW_biogeochem_input.K(1);
    PO4_in = ERW_biogeochem_input.PO4(1);
    NO3_in = ERW_biogeochem_input.NO3(1);
    NH4_in = ERW_biogeochem_input.NH4(1);
    
end

% or read from excel directly
Ca(1) = Ca_in;
Mg(1) = Mg_in;
K(1) = K_in;
Si(1) = Si_in;
Na(1) = Na_in;
Al_w(1) = Al_w_in;

PO4(1) = PO4_in;
NO3(1) = NO3_in;
NH4(1) = NH4_in;

if ischar(f_CEC_in)
    f_CEC_vals = str2num(f_CEC_in); %#ok<ST2NM>
else
    f_CEC_vals = f_CEC_in;
end
if ~isempty(f_CEC_vals) && any(isfinite(f_CEC_vals(:)))
    if numel(f_CEC_vals) ~= 6 || any(~isfinite(f_CEC_vals(:))) || ...
            any(f_CEC_vals(:) < 0) || abs(sum(f_CEC_vals) - 1) > 1e-3
        error('ERW:InvalidCECFractions', ...
            'f_CEC_in must contain six nonnegative fractions that sum to one.');
    end
    f_Ca(1) = f_CEC_vals(1);
    f_Mg(1) = f_CEC_vals(2);
    f_K(1) = f_CEC_vals(3);
    f_Na(1) = f_CEC_vals(4);
    f_Al(1) = f_CEC_vals(5);
    f_H(1) = f_CEC_vals(6);
elseif ~isempty(f_CEC_vals) && ~all(isnan(f_CEC_vals(:)))
    error('ERW:InvalidCECFractions', ...
        'f_CEC_in must be blank/NaN or contain six finite fractions.');
end
% pre-process excel inputs
if opt_ERW_preprocessing == 1
    vars = who();
    ERW_preprocess_data = struct();

    for ind_vars =1:length(vars)
        vn = vars{ind_vars};
        if isvarname(vn)
            ERW_preprocess_data.(vn)=eval(vn);
        end
    end

    [preprocess_vars, ERW_preprocess_output] = ERW_pre_processing(vars,ERW_preprocess_data);

    for ind_vars =1:length(preprocess_vars)
        vn = preprocess_vars{ind_vars};
        if isvarname(vn)
            eval([vn,'=ERW_preprocess_output.',vn,';'])
        end
    end
end

if active_rock_input(rock_f_in)
    if ~isfinite(diss_f)
        diss_f = 1;
    end
    if missing_numeric_vector(d_in) || missing_numeric_vector(psd_in)
        if any(strcmpi(string(mineral), "calcite")) && isfinite(SSA_in) && SSA_in > 0
            d_in = 1e-4;      %% [m] single fine-powder class; SSA_in sets reactive area
            psd_in = 1;
        else
            error('Active rock amendment requires d_in/psd_in or a rock PSD file.');
        end
    end
end

f_CEC(1,:) = [f_Ca(1),f_Mg(1),f_K(1),f_Na(1),f_Al(1),f_H(1)];
if abs(sum(f_CEC(1,:)) - 1) > 1e-3
    error('Sum of fractions must be 1');
end
K_CEC = [K_Ca_Mg, K_Ca_K, K_Ca_Na, K_Ca_Al, K_Ca_H];

%%% Soil-class Gaines-Thomas constants.
K_CEC_bulk = [K_Ca_Mg, K_Ca_K, K_Ca_Na, K_Ca_Al, K_Ca_H];
bulk_pH_params.K_CEC = K_CEC_bulk;
bulk_pH_params.cec_diag = cec_diag;
bulk_pH_params.CEC_cmolc_kg = CEC_tot_in;
bulk_pH_params.H_extract_anchor_umolL = NaN;
bulk_pH_params.H_pore_anchor_umolL = NaN;
bulk_pH_params.observed_pH_initial = pH_in;
% Field water converts the complete dissolved inventory between umol/m2
% and umol/L before that inventory is transferred to the air-dried extract.
bulk_pH_params.field_water_L_m2 = n(1) * s(1) * Zrw * 1000;
bulk_pH_params.H_prev = NaN;
bulk_pH_params.H_prev_valid = false;
if isfield(ERW_TC_const_input,'bulk_pH_params')
    fn_bulk = fieldnames(ERW_TC_const_input.bulk_pH_params);
    for i_bulk = 1:numel(fn_bulk)
        bulk_pH_params.(fn_bulk{i_bulk}) = ERW_TC_const_input.bulk_pH_params.(fn_bulk{i_bulk});
    end
end

% n_extract is set exclusively from pH_DilutionFactor in the Excel parameter
% table (L water / kg dry soil).  Runs after the merge so Excel always wins.
% Methods: 1:1 water:soil (US, n=1), 1:2.5 (ISO 10390, n=2.5), 1:5 (n=5).
if exist('pH_DilutionFactor','var') && isnumeric(pH_DilutionFactor) && ...
        isfinite(pH_DilutionFactor) && pH_DilutionFactor > 0
    bulk_pH_params.n_extract = pH_DilutionFactor;
end


%% Rock Applications
% rock surface fractality (Beerling 2020)
% consolidates measurements taken at different scales  - accounts
% uncertainties in grain topography and porositywhich influence mass
% transfer rates from rock grains to flowing soil water 
% Brantley, S. L. & Mellott, N. P. Surface area and porosity of primary silicate minerals. Am. Mineral. 85, 1767–1783 (2000)
b = 0.35; % [-] surface-area size-scaling exponent
a = (1/(2*1e-10))^b; %[1/m^b]


%% initialise rock properties

if active_rock_input(rock_f_in)

    number_min = length(mineral);
    rock_f_in = rock_f_in(:).';
    if numel(rock_f_in) > number_min
        rock_f_in = rock_f_in(1:number_min);
    elseif numel(rock_f_in) < number_min
        rock_f_in(numel(rock_f_in)+1:number_min) = 0;
    end
    
    MM_min = zeros(number_min,1);
    n_H = zeros(number_min,1);
    n_OH = zeros(number_min,1);
    K_sp = zeros(number_min,1);
    min_st = zeros(number_min,7);
    k_diss_H = zeros(number_min,1);
    k_diss_w = zeros(number_min,1);
    k_diss_OH = zeros(number_min,1);
    E_H = zeros(number_min,1);
    E_w = zeros(number_min,1);
    E_OH = zeros(number_min,1);
    
    k_diss_H_t = zeros(length(t), number_min);
    k_diss_w_t = zeros(length(t), number_min);
    k_diss_OH_t = zeros(length(t), number_min);
    Omega = zeros(length(t), number_min);
    M_min = zeros(length(t), number_min);
    rock_f = zeros(length(t), number_min);
    Wr = zeros(length(t), number_min);
    EW = zeros(length(t), number_min);
    
    n_d_cl = length(d_in);

    if n_d_cl >= 1
        d = zeros(length(t), n_d_cl, num_rock_app);
        delta_d = zeros(length(t), n_d_cl, num_rock_app);
        lamb = zeros(length(t), n_d_cl, num_rock_app);
        SSA = zeros(length(t), n_d_cl, num_rock_app);
        psd = zeros(length(t), n_d_cl, num_rock_app);
        a = zeros(1, num_rock_app);
    end

    for min_j = 1:length(mineral)
        [MM_min(min_j),k_diss_H(min_j),k_diss_w(min_j),k_diss_OH(min_j),E_H(min_j),E_w(min_j),E_OH(min_j),n_H(min_j),n_OH(min_j),min_st(min_j,:),K_sp(min_j)] = min_const(mineral(min_j),conv_mol);
    end

else

    number_min = 1;
    
    MM_min = zeros(number_min,1);
    n_H = zeros(number_min,1);
    n_OH = zeros(number_min,1);
    K_sp = zeros(number_min,1);
    min_st = zeros(number_min,7);
    k_diss_H = zeros(number_min,1);
    k_diss_w = zeros(number_min,1);
    k_diss_OH = zeros(number_min,1);
    E_H = zeros(number_min,1);
    E_w = zeros(number_min,1);
    E_OH = zeros(number_min,1);
    
    k_diss_H_t = zeros(length(t), number_min);
    k_diss_w_t = zeros(length(t), number_min);
    k_diss_OH_t = zeros(length(t), number_min);
    Omega = zeros(length(t), number_min);
    M_min = zeros(length(t), number_min);
    rock_f = zeros(length(t), number_min);
    Wr = zeros(length(t), number_min);
    EW = zeros(length(t), number_min);

    d = zeros(length(t), number_min);
    delta_d = zeros(length(t), number_min);
    lamb = zeros(length(t), number_min);
    SSA = zeros(length(t), number_min);
    psd = zeros(length(t), number_min);

    M_iner = 0;
    psd_perc_in = 0;
    rock_f_in = 0;
    d_in = 0;
    SSA_in = 0; 
    a = 0;

end



%% other parameters
% temperature and gas constants 
R = 8.314/conv_mol; % [J mol-1 K-1]: universal gas constant
T_ref = 25*ones(length(t),1)+273.15; % [K]: temperature standard conditions

% soil CO2 diffusivity
D_0 = 1.6e-5*3600*24; %free-air diffusion [m2/d]
D(1) = D_0*(1-s(1,1))^(10/3)*n(1,1)^(4/3); % [m²/d] - Mill, T., & Quirk, J. (1961). Permeability of porous solids. [Link](https://pubs.rsc.org/en/content/articlelanding/1961/tf/tf9615701200)

% solute diffusivity in soil water
Dw_0 = 1e-9*3600*24; %[m2/d]
Dw(1) = Dw_0*(n(1)*s(1))^2; % Archie 1942, Grathwohl 1998 (book)

% soil temperature [K]
T_K(1) = temp_soil(1)+273.15;

% Carbonate weathering [mol-conv/d]
% carbonate speciation
pk1(1) = -(-356.309-0.0609*T_K(1)+21834.37/T_K(1)+126.8339*log10(T_K(1))-1684915/(T_K(1))^2);
pk2(1) = -(-107.887-0.032528*(T_K(1))+5151.79/(T_K(1))+38.92561*log10(T_K(1))-563713.9/(T_K(1))^2);
pk_w(1) = -(-283.971+13323/(T_K(1))-0.0507*(T_K(1))+102.24447*log10(T_K(1))-1119669/(T_K(1))^2);
k1(1) = 10^(-pk1(1))*conv_mol;
k2(1) = 10^(-pk2(1))*conv_mol;
%water
k_w(1) = 10^(-pk_w(1))*conv_mol^2;
% ammonium acid dissociation NH4+ = NH3(aq) + H+ ; pKa(T) Emerson et al. (1975)
pk_nh3(1) = 0.09018 + 2729.92/T_K(1);
k_nh3(1) = 10^(-pk_nh3(1))*conv_mol;
% Henry
k_H(1) = 0.83*exp(2400*(1/T_K(1)-1/T_ref(1)));

%CHANGE DEPENDING ON PFT OF INTEREST
% default tissue fractions remain within the cited agronomic ranges
v_f_Ca = 0.005; % suggested 2%, range 0.1 - 5 %, Weil and Bredy 2017
v_f_Mg = 0.001; % suggested, 0.5%, range 0.1 - 1 %, Weil and Bredy 2017
v_f_K = 0.01; % suggested 2%, range 1 - 5 %, Weil and Bredy 2017
v_f_Si = 0.01; % suggested 5%, range 1 - 10%, Epstein 1994, PNAS
    
% percentage of element in plant tissue
dry_perc = 0.1;%percent of dry mass

% amount of nutrient per unit dry biomass
xi_Ca = dry_perc*v_f_Ca/MM_Ca;%[mol-conv/g_biomass]
xi_Mg = dry_perc*v_f_Mg/MM_Mg;
xi_Si = dry_perc*v_f_Si/MM_Si;

% root area index - Porporato, A., & Yin, J. (2022), Ecohydrology
%RAI = 9.8; % temperate deciduous [m2/m2]
% RAI = 11; % temperate spruce 
RAI = 79.1; % temperate grassland

% root diameter - Jackson et al., 1997 - https://www.pnas.org/doi/10.1073/pnas.94.14.7362?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%20%200pubmed#sec-2
%root_d = 0.58*1e-3; % trees [m]
root_d = 0.22*1e-3; % grasses

% carrying capacity - doi: 10.1111/j.1365-2486.2009.02146.x
%k_v = 202*1e2; % international average [g/m2]
k_v = 1000*1e2; % temperate forest


%Carb weathering constants
% Carbonate solutibility products
% https://booksite.elsevier.com/9780120885305/appendices/Web_Appendices.pdf
K_CaCO3 = 10^(-8.35)*conv_mol^2; % Calcite [mol^2]
K_MgCO3 = 10^(-7.46)*conv_mol^2; % Magnesite [mol^2]

% precipitation rate
% https://nora.nerc.ac.uk/id/eprint/511084/1/Kirk%20et%20al%202015%20Geochmica%20et%20Cosmochimica%20Acta.pdf
r_CaCO3 = 3*1e7*(1e-9*1e6/(24*3600))*conv_mol; % [mol-conv/d]
r_MgCO3 = 1e7*(1e-9*1e6/(24*3600))*conv_mol; % [mol-conv/d]



% rainwater  
CO2_w_rain(1) = k_H(1)*CO2_atm(1); % [mol/l] Henry's law

%pCO2 - partial pressure of CO2 [Pa]
if Zrw<=0.3
    Z_CO2=Zrw/2;
else
    Z_CO2 = 0.15;
end

CO2_air(1) = max(TT_m_input.r_het+TT_m_input.r_aut)/(D(1)*1000/Z_CO2)+CO2_atm(1);

Fs(1)= D(1)/Z_CO2*(CO2_air(1)-CO2_atm(1))*1000; %[umol-conv/d]
CO2_w(1)= k_H(1)*CO2_air(1); %[umol-conv/l] Henry's law

% The deposition strong-anion budget is always computed: it supplies the
% Phase-3 forcing flux (I_An_source) and the charge-residual diagnostic below.
anion_budget = ERW_Initial_Anion_Budget(ERW_TC_const_input);

% Bulk-extract exchange constants stay GLOBAL (no per-site back-calc).
K_Ca_H_bulk = K_Ca_H;
K_CEC_bulk = [K_Ca_Mg, K_Ca_K, K_Ca_Na, K_Ca_Al, K_Ca_H_bulk];
bulk_pH_params.K_CEC = K_CEC_bulk;
f_CEC_seed = [f_Ca(1), f_Mg(1), f_K(1), f_Na(1), f_Al(1), f_H(1)];
%%% Keep the spin-up exchange seed when supplied.
if ~isfield(bulk_pH_params, 'initial_f_CEC') ...
        && all(isfinite(f_CEC_seed)) && all(f_CEC_seed >= 0) && sum(f_CEC_seed) > 0
    bulk_pH_params.initial_f_CEC = f_CEC_seed / sum(f_CEC_seed);
end
if ~isfield(bulk_pH_params, 'initial_Al_w')
    if isfinite(Al_w(1)) && Al_w(1) >= 0
        bulk_pH_params.initial_Al_w = Al_w(1);
    elseif f_Al(1) == 0
        bulk_pH_params.initial_Al_w = 0;
    end
end

% Seed pH anchoring.  When an observed bulk extract pH and a dilution factor
% are available, invert the forward extraction operator (Porewater_to_Bulk_pH)
% to find the pore-water H whose forward extraction matches the observed bulk pH
% (Invert_Bulk_pH).  Missing extraction data or a nonphysical root is fatal.
if ~isfield(bulk_pH_params,'n_extract') || ...
        ~(isfinite(bulk_pH_params.n_extract) && bulk_pH_params.n_extract > 0)
    error('ERW:missingParam', ...
        'Observed bulk pH initialization requires a positive extraction ratio.');
end
[H(1), bulk_pH_params, f_Ca(1), f_Mg(1), f_K(1), f_Na(1), f_Al(1), f_H(1), ...
    An(1), Al_w(1), Ca(1), Mg(1), K(1), Na(1)] = Invert_Bulk_pH( ...
    pH_in, Ca(1), Mg(1), K(1), Na(1), NO3(1), PO4(1), NH4(1), ...
    CO2_atm(1), CO2_w(1), k1(1), k2(1), k_w(1), k_nh3(1), k_H(1), ...
    K1, K2, K3, K4, CEC_tot(1), rho_bulk, Zrw, conv_mol, ...
    bulk_pH_params, K_AlOH3, Al_CEC_scale);

% Appelo + Al-solubility-pin re-baseline of the legacy Na pool.  The
% legacy 'Kantola' seed convention (ERW_pre_processing) fills the sodium
% slot f_Na = f_CEC_sat - f_Ca - f_Mg - f_K ~ 0.8 and derives
% Na_w = sqrt((Ca/f_Ca)*K_Ca_Na)*f_Na from the pre-correction
% Gaines-Thomas form, seeding a ~0.3-0.7 M porewater Na pool (with the
% matching exchange inventory) that the corrected Appelo selectivity
% cannot host: at the Al-saturated pinned root (f_H ~ 1e-17) the law
% f_Na = K_Na_H*Na_w*(f_H/H)*10^(alpha*f_H) caps exchange Na at O(1e-9)
% while the pool demands f_Na ~ 1e-4 - no joint solution exists.
% Re-seed Na_w to the rain-supported steady state (zero when no rain-Na
% flux is supplied) and rescale f_Na law-consistently (linear in Na_w at
% fixed H and f_H, so the inversion's law ratios and the CEC sum are
% preserved); the pinned run then starts from a state the daily
% conservation + exchange system can actually host.  The frame-side
% companion disables the self-inflating background pump (keyword_add)
% for pinned runs.
if isfield(bulk_pH_params, 'Al_solubility_pin') && ...
        isequal(bulk_pH_params.Al_solubility_pin, true)
    Na_legacy = Na(1);
    if exist('Na_rain_ss', 'var') && isfinite(Na_rain_ss) && Na_rain_ss > 0
        Na_target = Na_rain_ss;
    else
        Na_target = 0;
    end
    % Only salt-paired Na can leave: charge moves as neutral salt pairs
    % (SCEPTER get_soilpH_time.py transfers inventories to the extract
    % run as Na2O/NaCl and hard-stops on Na < Cl - a strong anion with
    % no pairable cation is an inadmissible state).  The carbonate-frame
    % counter-anion recompute below (An = 2Mg+2Ca+Na+K+NH4-Alk-PO4-NO3)
    % needs Na >= Alk+PO4+NO3-2Mg-2Ca-K-NH4 to close with An >= 0;
    % dropping below that floor leaves the modelled anion excess
    % (Kantola porewater PO4 above all) without countercharge, the
    % max(An,0) floor breaks the Alk_tot = R_alk + Alk*V identity, and
    % the daily solver inherits a negative Alk_tot it can only host as a
    % hyper-acid f_Al = 1 state (day-2 hard fail).  The retained Na is
    % the minimal electroneutral remainder of the legacy pool.
    Alk_seed = k1(1)*CO2_w(1)/H(1) + 2*k2(1)*k1(1)*CO2_w(1)/(H(1)^2) ...
        - H(1) + k_w(1)/H(1);
    Na_floor = max(0, Alk_seed + PO4(1) + NO3(1) ...
        - (2*Mg(1) + 2*Ca(1) + K(1) + NH4(1)));
    Na(1) = min(max(Na_target, Na_floor), max(Na_legacy, Na_target));
    if isfinite(Na_legacy) && Na_legacy > 0
        f_Na_new = f_Na(1) * (Na(1)/Na_legacy);
    else
        f_Na_new = f_Na(1);
    end
    f_renorm = 1 - f_Na(1) + f_Na_new;
    if isfinite(f_renorm) && f_renorm > 1e-6
        f_Ca(1) = f_Ca(1)/f_renorm;
        f_Mg(1) = f_Mg(1)/f_renorm;
        f_K(1)  = f_K(1)/f_renorm;
        f_Al(1) = f_Al(1)/f_renorm;
        f_H(1)  = f_H(1)/f_renorm;
        f_Na(1) = f_Na_new/f_renorm;
    else
        % Degenerate legacy seed (exchange Na ~ full CEC): drop the slot.
        f_Na(1) = 0;
    end
    % The removed dissolved Na leaves as neutral salt, taking the same
    % equivalents off the inverter's counter-anion (exact in that frame:
    % only the Na term of its charge closure changed).  The Kantola
    % convention below overwrites An(1) with its own carbonate-frame
    % recompute, which the Na floor above keeps nonnegative.
    if isfinite(Na_legacy) && Na_legacy > Na(1)
        An(1) = max(0, An(1) - (Na_legacy - Na(1)));
    end
end

pH(1) = -log10(H(1)/conv_mol);
Al_den_seed = H(1)^4 + H(1)^3*K1 + H(1)^2*K1*K2 + ...
    H(1)*K1*K2*K3 + K1*K2*K3*K4;
Al(1) = (H(1)^4/Al_den_seed) * Al_w(1);
f_CEC(1,:) = [f_Ca(1), f_Mg(1), f_K(1), f_Na(1), f_Al(1), f_H(1)];
bulk_pH_params.H_init = 10^(-pH_in) * conv_mol;
bulk_pH_params.H_extract_anchor_umolL = 10^(-pH_in) * conv_mol;
bulk_pH_params.H_pore_anchor_umolL = H(1);

%carbonate system
HCO3(1)= k1(1)*CO2_w(1)/H(1); % bicarbonate [umol/l]
CO3(1)= k2(1)*k1(1)*CO2_w(1)/(H(1)^2); % carbonate [umol/l]
DIC(1)= HCO3(1)+CO3(1)+CO2_w(1); % dissolved inorganic C (includes bicarbonate, carbonate and dissolved CO2 (carbonic acid))
IC_tot(1)=(DIC(1)*s(1)+CO2_air(1)*(1-s(1)))*(n(1)*Zrw*1000);% total inorganic C [umol]

%Alkalinity
Alk(1)= HCO3(1)+2*CO3(1)-H(1)+k_w(1)/H(1); % [umol/L]

% Al from Ca-Al GT exchange at initial pore-water pH
% If exchange inputs fail, use Al(OH)3 saturation.
if isnan(Al_w_in) && ~(isfinite(Al_w(1)) && Al_w(1) >= 0)
    Al_w_eq = Initialize_ERW_Equilibrium(H(1), Ca(1), f_Ca(1), f_Al(1), ...
        K_Ca_Al, K1, K2, K3, K4, Al_CEC_scale);
    if isnan(Al_w_eq)
        Al_w_eq = Initialize_ERW_Equilibrium(H(1), K1, K2, K3, K4, ...
            conv_mol, K_AlOH3, Al_CEC_scale);
    end
    if ~isfinite(Al_w_eq) || Al_w_eq <= 0
        error('Cannot initialize Al_w from exchange or Al(OH)3 equilibrium.');
    end
    Al_w(1) = Al_w_eq;
end

% aluminium speciation
Al(1)=(H(1)^4/(H(1)^4+H(1)^3*K1+H(1)^2*K1*K2+H(1)*K1*K2*K3+K1*K2*K3*K4))*Al_w(1);% [mol/l]
AlOH(1)=(H(1)^3*K1/(H(1)^4+H(1)^3*K1+H(1)^2*K1* K2+H(1)*K1*K2*K3+K1*K2*K3*K4))* Al_w(1);
AlOH2(1)=(H(1)^2*K1*K2/(H(1)^4+H(1)^3*K1+H(1)^2*K1*K2+H(1)*K1*K2*K3+K1*K2*K3*K4))*Al_w(1);
AlOH3(1)=(H(1)*K1*K2*K3/(H(1)^4+H(1)^3*K1+H(1)^2*K1*K2+H(1)*K1*K2*K3+K1*K2*K3*K4))*Al_w(1);
AlOH4(1)= Al_w(1)-(Al(1)+AlOH(1)+AlOH2(1)+AlOH3(1));

Al_w_sat(1) = Initialize_ERW_Equilibrium(H(1), K1, K2, K3, K4, conv_mol, K_AlOH3, Al_CEC_scale);
Omega_AlOH3(1) = Al_w(1) / max(Al_w_sat(1), eps);

% Initial An is the unmodelled strong-anion countercharge returned by the
% observed-pH inverse.  The deposition budget below supplies future fluxes only.
%
% Kantola convention: the daily solver closes the charge balance with the
% carbonate-only species (eq3: HCO3+2CO3-H+k_w/H), and the legacy seed An
% was the carbonate-frame countercharge.  The bulk-pH inverse returns its
% own counter-anion from its fuller charge accounting, which on high-Al
% soils (Kantola Al pools ~2e10 umol/m2) makes Alk_tot negative and leaves
% the daily balance with no physical root at near-neutral pH.  Restore the
% legacy carbonate-frame counter-anion so that
%   Alk_tot = R_alk + Alk*porewater_volume
% stays positive and the daily system converges.
if string(excel_input.fX_conc_convention) == "Kantola"
    An(1) = 2*Mg(1)+2*Ca(1)+Na(1)+K(1)+NH4(1)-Alk(1)-PO4(1)-NO3(1);
    An(1) = max(An(1), 0);
end

R_alk(1)=(f_Mg(1)+f_Ca(1)+f_Na(1)+f_K(1))*CEC_tot(1); % [umol_c]
Ca_tot(1) = Ca(1)*n(1)*s(1)*Zrw*1000+f_Ca(1)/2*CEC_tot(1);
Mg_tot(1) = Mg(1)*n(1)*s(1)*Zrw*1000+f_Mg(1)/2 *CEC_tot(1);
K_tot(1) = K(1)*n(1)*s(1)*Zrw*1000+f_K(1)*CEC_tot(1);
Si_tot(1) = Si(1)*n(1)*Zrw*s(1)*1000;
Na_tot(1) = Na(1)*n(1)*s(1)*Zrw*1000+ f_Na(1)* CEC_tot(1);
Al_tot(1) = Al_w(1)*n(1)*Zrw*s(1)*1000+(f_Al(1)/3)*CEC_tot(1)*Al_CEC_scale;
An_tot(1) = An(1)*n(1)*s(1)*Zrw*1000;

% Phosphorus ions (solution)
PO4_tot(1) = PO4(1)*n(1)*Zrw*s(1)*1000;

% Nitrate ions (solution)
NO3_tot(1) = NO3(1)*n(1)*Zrw*s(1)*1000;

% Ammonium ions (solution)
NH4_tot(1) = NH4(1)*n(1)*Zrw*s(1)*1000;

% alkalinity total
Alk_tot(1) = 2*Mg_tot(1)+2*Ca_tot(1)+Na_tot(1)+K_tot(1)+NH4_tot(1)-NO3_tot(1)-PO4_tot(1)-An_tot(1);

% Background inputs (rain, fertilizer, litterfall, background weathering)
I_An_source = anion_budget.I_An_source;
I_An = anion_budget.I_An_initial;
I_Na_source = 0;
if isfield(ERW_TC_const_input, 'I_Na_source') && ~isempty(ERW_TC_const_input.I_Na_source)
    I_Na_source = ERW_TC_const_input.I_Na_source(:);
elseif isfield(ERW_TC_const_input, 'rain_Na_flux_umol_m2_d')
    I_Na_source = ERW_TC_const_input.rain_Na_flux_umol_m2_d;
end
I_Na_source(~isfinite(I_Na_source) | I_Na_source < 0) = 0;
I_Na = first_finite_source(I_Na_source);

I_Ca = 0;
I_Mg = 0;
I_K = 0;
I_Si = 0;

%Carbonate weathering
Omega_CaCO3(1) =Ca(1)*CO3(1)/K_CaCO3;% carbonate saturation index [-]
Omega_MgCO3(1) = Mg(1)*CO3(1)/K_MgCO3;

% dissolution timescale
tau_CaCO3 = 30; % [d]
tau_MgCO3 = 40;

% In soil, precipitates form as discontinuous coatings on the surfaces of soil pores,
% so the precipitation surface area and geometry are indeterminate.
% https://nora.nerc.ac.uk/id/eprint/511084/1/Kirk%20et%20al%202015%20Geochmica%20et%20Cosmochimica%20Acta.pdf
% CaCO3
if Omega_CaCO3(1)<=1
    W_CaCO3(1) = s(1)*CaCO3(1)*(1-Omega_CaCO3(1))/tau_CaCO3; % dissolution [umol/d]
else
    W_CaCO3(1) = r_CaCO3*Zrw*(1-Omega_CaCO3(1)); % precipitation [umol/d]
end

% MgCO3
if Omega_MgCO3(1) <= 1
    W_MgCO3(1) = s(1)*MgCO3(1)*(1-Omega_MgCO3(1))/tau_MgCO3; % dissolution [umol/d]
else
    W_MgCO3(1) = r_MgCO3*Zrw*(1-Omega_MgCO3(1)); % precipitation [umol/d]
end

if active_rock_input(rock_f_in) && num_rock_app > 0

    k_diss_H_t(1,:) = k_diss_H.*exp(-E_H.*1000/R*(1/T_K(1,1)-1/T_ref(1))); % temperature-dependent dissociation constant with respect to hydrogen [umol/m²/d]
    k_diss_w_t(1,:) = k_diss_w.*exp(-E_w.*1000/R*(1/T_K(1,1)-1/T_ref(1))); % temperature-dependent dissociation constant with respect to water [umol/m²/d]
    k_diss_OH_t(1,:) = k_diss_OH.*exp(-E_OH.*1000/R*(1/T_K(1,1)-1/T_ref(1))); % temperature-dependent dissociation constant with respect to hydroxide [umol/m²/d]

    %rock composition
    M_rock(1) = M_rock_in; % rock mass [g/m2]
    rock_f(1,:)= rock_f_in;
    M_min(1,:)= rock_f(1,:)*M_rock(1); % mineral mass [g/m2]
    M_iner(1) = M_rock(1)*(1-sum(rock_f(1,:))); % inert rock mass??[g/m2]
    
    %particle sediment diameter
    delta_d_in = [d_in(1,1) diff(d_in(1,:))];
    psd_perc_in = psd_in; % particle class weight (%)
    
    %refinment of fractal constant based on SSA measured
    if SSA_in > 0
        if M_rock_in > 0
            M_rock_add_amt = M_rock_in;
        else
            M_rock_add_amt = M_rock_add(M_rock_add>0);
        end
        psd_in = psd_perc_in.*M_rock_add_amt./delta_d_in; % [g/m]
        a =(SSA_in*rho_rock*M_rock_add_amt/6)./sum(d_in.^(b-1).*psd_in.*delta_d_in, 2);% [m^-b]
    end

    if M_rock_in >0 || M_rock_add(1)>0
        %diameter classes
        d(1,:,1)= d_in; % particle diameter [m]
        delta_d(1,:,1)=[d(1,1,1) diff(d(1,:,1))]; % change in particle diamter [m]
        psd(1,:,1) = psd_perc_in.*M_rock(1)./delta_d(1,:,1); % [g/m]

        %surface area
        lamb(1,:,1)= a(1)*d(1,:,1).^b;% accounts for topogrpay and porosity in relation to grain diameter [-]
        SSA(1,:,1)= 6./(d(1,:,1).*rho_rock).*lamb(1,:,1);% specific surface area [m2/g]

        SA(1,:)= sum(SSA(1,:,:).*psd(1,:,:).*delta_d(1,:,:),'all');% surface area [m2]
        
        for min_j = 1:length(mineral)
            Omega(1,min_j) = Omega_sil(mineral(min_j), K(1), Ca(1), Mg(1), Na(1), Si(1), Al(1), PO4(1), H(1), K_sp(min_j), conv_mol);
            Wr(1, min_j) = s(1).*diss_f.*(k_diss_H_t(1,min_j).*(H(1)/conv_mol).^n_H(min_j)+k_diss_w_t(1,min_j)+k_diss_OH_t(1,min_j).*(H(1)/conv_mol).^n_OH(min_j)).*(1-Omega(1,min_j)); % background weathering rate [mol-conv/ m2 d]
            EW(1, min_j) = Wr(1, min_j).*SA(1).*rock_f(1, min_j); % enhanced weathering rate [umol/d]
        end

    end

end

%%
% structure containing biogeochemical variables
ERW_biogeochem_varnames = {'pH';'CO2_air';'CO2_atm';'Ca';'Mg';'K';'Na';'An';'Al';'AlOH';'AlOH2';'AlOH3';'AlOH4';'Si';'PO4';'NO3';'NH4';'H';
    'Ca_tot';'Mg_tot';'K_tot';'Na_tot';'An_tot';'Al_tot';'Si_tot';'PO4_tot';'NO3_tot';'NH4_tot';'IC_tot';'CEC_tot';'K_CEC';
    'R_alk';'Fs';'DIC';'CaCO3';'MgCO3';'f_Ca';'f_Mg';'f_K';'f_Na';'f_Al';'f_H';'f_CEC';'M_min';'CO2_w';'Alk';'Alk_tot';'k_H';'ADV';'DIC_rain';'D';'CO3';...
    'Al_w';'Al_w_sat';'HCO3';'Leak_Na';'Leak_Al';'Leak_An';'Leak_DIC';'Na_uptake';'An_uptake';...
    'k1';'k2';'k_w';'k_nh3'};

ERW_biogeochem = struct();
for var_index = 1:size(ERW_biogeochem_varnames,1)
    vn = cell2mat(ERW_biogeochem_varnames(var_index));
    if exist(vn, 'var')
        ERW_biogeochem.(vn) = eval(vn);
    end
end
ERW_biogeochem.pH_bulk = nan(length(t), 1);
ERW_biogeochem.pH_bulk(1) = pH_in;

% structure containing weathering variables
ERW_weathering_varnames = {'EW';'Wr';'psd';'d';'lamb';'SSA';
    'W_CaCO3';'W_MgCO3';'Omega';'Omega_AlOH3';'delta_d';'SA'};

ERW_weathering = struct();
for var_index = 1:size(ERW_weathering_varnames,1)
    vn = cell2mat(ERW_weathering_varnames(var_index));
    if exist(vn, 'var')
        ERW_weathering.(vn) = eval(vn);
    end
end

% structure containing weathering variables
ERW_plant_varnames = {'plant_H'}; %{'UP_Mg';'UP_Ca';'UP_Si';'UP_K';'plant_H'}; 

ERW_plant = struct();
for var_index = 1:size(ERW_plant_varnames,1)
    vn = cell2mat(ERW_plant_varnames(var_index));
    if exist(vn, 'var')
        ERW_plant.(vn) = eval(vn);
    end
end

% structure containing constants
ERW_const_varnames = {'n'; 'Zrw'; 'K_CEC'; 'M_rock_in'; 'k_diss_H'; 'k_diss_w'; 'k_diss_OH'; 'E_H'; 
    'E_w'; 'E_OH'; 'R'; 'T_ref'; 'K1'; 'K2'; 'K3'; 'K4'; 'rho_rock'; 'b'; 'K_CaCO3'; 'K_MgCO3'; 'MM_min'; 
    'min_st'; 'Dw_0'; 'D_0'; 'pH_rain'; 'Z_CO2'; 'I_Ca'; 'I_Mg'; 'I_Na'; 'I_Na_source'; 'I_An'; 'I_An_source'; 'I_Si';'I_K';
    'r_CaCO3'; 'r_MgCO3'; 'tau_CaCO3'; 'tau_MgCO3';
    'conv_mol';
    'M_iner'; 'K_sp'; 'diss_f'; 'n_H'; 'n_OH';'psd_perc_in';'rock_f_in';'d_in';'SSA_in'; 'a'; 'rho_bulk'; 'Al_CEC_scale'; 'K_AlOH3'};

ERW_const = struct();
for var_index = 1:size(ERW_const_varnames,1)
    vn = cell2mat(ERW_const_varnames(var_index));
    if exist(vn, 'var')
        ERW_const.(vn) = eval(vn);
    end
end

ERW_const.bulk_pH_params = bulk_pH_params;

% structure containing other variables not otherwise prescribed to a strucutre 
ERW_var_varnames = {'fsolve_options';'mineral'};

ERW_var = struct();
for var_index = 1:size(ERW_var_varnames,1)
    vn = cell2mat(ERW_var_varnames(var_index));
    if exist(vn, 'var')
        ERW_var.(vn) = eval(vn);
    end
end


% all variables for checking
all_vars = who();
ERW_all = struct();
for var_index = 1:size(all_vars,1)
    vn = cell2mat(all_vars(var_index));
    if isvarname(vn) && exist(vn, 'var')
        ERW_all.(vn) = eval(vn);
    end
end

end

%% sub-functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Q = ERW_Seed_Water_Loss(ERW_TC_const_input, TT_m_input)
%%% Annual-mean water loss [L m-2 d-1] behind the steady-state solute seeds,
%%% supplied by the spin-up over a complete annual cycle.  When the spin-up
%%% handoff is missing (direct transient start), the mean driving flux
%%% (Tw + Lw, m/d -> L/m2/d) over the simulation window is used instead.
if isstruct(ERW_TC_const_input) && isfield(ERW_TC_const_input,'seed_water_loss_l_m2_d')
    Q = ERW_TC_const_input.seed_water_loss_l_m2_d;
else
    Q = NaN;
end
if ~(isscalar(Q) && isfinite(Q) && Q > 0)
    Q = NaN;
    if isstruct(TT_m_input) && isfield(TT_m_input,'Tw') && isfield(TT_m_input,'Lw')
        Tw_v = TT_m_input.Tw(:);
        Lw_v = TT_m_input.Lw(:);
        if ~isempty(Tw_v) && ~isempty(Lw_v) && all(isfinite(Tw_v)) && all(isfinite(Lw_v))
            Q = mean(Tw_v + Lw_v) * 1000;
        end
    end
end
if ~(isscalar(Q) && isfinite(Q) && Q > 0)
    Q = 1;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function budget = ERW_Initial_Anion_Budget(ERW_TC_const_input)
I_An_source = 0;
if isfield(ERW_TC_const_input, 'I_An_source') && ~isempty(ERW_TC_const_input.I_An_source)
    I_An_source = ERW_TC_const_input.I_An_source(:);
elseif isfield(ERW_TC_const_input, 'rain_strong_anion_flux_umolc_m2_d')
    I_An_source = ERW_TC_const_input.rain_strong_anion_flux_umolc_m2_d;
end
I_An_source(~isfinite(I_An_source) | I_An_source < 0) = 0;
budget = struct();
budget.I_An_source = I_An_source;
budget.I_An_initial = first_finite_source(I_An_source);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function I0 = first_finite_source(I_An_source)
idx = find(isfinite(I_An_source) & I_An_source >= 0, 1, 'first');
if isempty(idx)
    I0 = 0;
else
    I0 = I_An_source(idx);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = parse_numeric_input(val)

if isstring(val) && (isempty(val) || ismissing(val(1)) || strlength(strtrim(val(1))) == 0)
    val = NaN;
    return
end
if isstring(val)
    val = char(val);
end
if ischar(val)
    if strcmpi(strtrim(val), 'NaN')
        val = NaN;
    else
        val = str2num(val); %#ok<ST2NM>
        if isempty(val)
            val = NaN;
        end
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mineral = parse_mineral_input(val)

if iscell(val)
    val = val{1};
end
if isstring(val) && (isempty(val) || ismissing(val(1)) || strlength(strtrim(val(1))) == 0)
    mineral = strings(1,0);
    return
end
if isstring(val)
    val = char(val);
end
if isnumeric(val)
    mineral = strings(1,0);
    return
end
txt = strtrim(val);
if strcmpi(txt, 'NaN')
    mineral = strings(1,0);
    return
end
txt = strrep(txt, '[', '');
txt = strrep(txt, ']', '');
txt = strrep(txt, '"', '');
txt = strrep(txt, '''', '');
parts = strtrim(strsplit(txt, ','));
parts = parts(~cellfun(@isempty, parts));
mineral = string(parts);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tf = active_rock_input(rock_f_in)

tf = isnumeric(rock_f_in) && ~isempty(rock_f_in) && ...
    all(isfinite(rock_f_in(:))) && any(rock_f_in(:) > 0);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tf = missing_numeric_vector(x)

tf = ~isnumeric(x) || isempty(x) || all(isnan(x(:)));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [MM_min,k_diss_H,k_diss_w,k_diss_OH,E_H,E_w,...
    E_OH,n_H,n_OH,min_st,K_sp] = min_const(mineral,conv_mol)

    if strcmp(mineral, 'forsterite') % Mg2SiO4
        MM_min = 140/conv_mol; % molar mass of mineral [g/mol-conv]
        k_diss_H = 10^(-6.85)*24*3600*conv_mol; % Dissolution rate of the mineral with respect to hydrogen [mol-conv m-2 d-1]
        k_diss_w = 10^(-10.64)*24*3600*conv_mol; % Dissolution rate of the mineral with respect to water [mol-conv m-2 d-1]
        k_diss_OH = 0; % Dissolution rate of the mineral with respect to hydroxide [mol-conv m-2 d-1]
        E_H = 67.2/conv_mol; % [kJ/mol-conv]: activation energies with respect to hydrogen
        E_w = 79/conv_mol; % [kJ/mol-conv]: activation energies with respect to water
        E_OH = 0; % [kJ/mol-conv]: activation energies with respect to hydroxide 
        n_H = 0.47; % reaction order with respect to hydrogen
        n_OH = 1; % % reaction order with respect to hydroxide
        min_st = [0, 2, 0, 0, 0, 1, 0]; % Stochiometric coefficients [Ca, Mg, K, Na, Al, Si, P]
        K_sp = 10^(7.11); % Solubility product of the mineral [-] https://booksite.elsevier.com/9780120885305/appendices/Web_Appendices.pdf
    elseif strcmp(mineral, 'Fe_forsterite') % FeMgSiO4
        MM_min = 172/conv_mol; % [g/mol-conv]
        k_diss_H = 10^(-6.85)*24*3600*conv_mol; % (-5.37) [mol-conv m-2 d-1]
        k_diss_w = 10^(-10.64)*24*3600*conv_mol; % [mol-conv m-2 d-1]
        k_diss_OH = 0;
        E_H = 67/conv_mol; % [kJ/mol-conv]: activation energies
        E_w = 79/conv_mol; % [kJ/mol-conv]
        E_OH = 0;
        n_H = 0.47; % reaction order
        n_OH = 1;
        min_st = [0, 1, 0, 0, 0, 1, 0]; % Stochiometric coefficients [Ca, Mg, K, Na, Al, Si, P]
        K_sp = NaN;
    elseif strcmp(mineral, 'wollastonite') % CaSiO3
        MM_min = 116/conv_mol; % [g/mol-conv]
        k_diss_H = 10^(-5.37)*24*3600*conv_mol; % [mol-conv m-2 d-1]
        k_diss_w = 10^(-8.88)*24*3600*conv_mol; % [mol-conv m-2 d-1]
        k_diss_OH = 0;
        E_H = 54.7/conv_mol; % [kJ/mol-conv]: activation energies
        E_w = 54.7/conv_mol; % [kJ/mol-conv]
        E_OH = 0;
        n_H = 0.4; % reaction order
        n_OH = 1;
        min_st = [1, 0, 0, 0, 0, 1, 0]; % Stochiometric coefficients [Ca, Mg, K, Na, Al, Si, P]
        K_sp = 10^(6.82);
    elseif strcmp(mineral, 'albite') % NaAlSiO3 - Plagioclase (NaAlSi3O8)
        MM_min = 263/conv_mol; % [g/mol-conv]
        k_diss_H = 10^(-10.16)*24*3600*conv_mol; % [mol-conv m-2 d-1]
        k_diss_w = 10^(-12.56)*24*3600*conv_mol; % [mol-conv m-2 d-1]
        k_diss_OH = 10^(-15.6)*24*3600*conv_mol;
        E_H = 65/conv_mol; % [kJ/mol-conv]: activation energies
        E_w = 70/conv_mol; % [kJ/mol-conv]
        E_OH = 71/conv_mol;
        n_H = 0.457; % reaction order
        n_OH = -0.572;
        min_st = [0, 0, 0, 1, 1, 3, 0]; % Stochiometric coefficients [Ca, Mg, K, Na, Al, Si, P]
        K_sp = 10^(-0.68);
    elseif strcmp(mineral, 'anorthite') % CaAl2Si2O8
        MM_min = 278/conv_mol; % [g/mol-conv]
        k_diss_H = 10^(-3.5)*24*3600*conv_mol; % [mol-conv m-2 d-1]
        k_diss_w = 10^(-9.12)*24*3600*conv_mol; % [mol-conv m-2 d-1]
        k_diss_OH = 0;
        E_H = 16.6/conv_mol; % [kJ/mol-conv]: activation energies
        E_w = 17.8/conv_mol; % [kJ/mol-conv]
        E_OH = 0;
        n_H = 1.4; % reaction order
        n_OH = 1;
        min_st = [1, 0, 0, 0, 2, 2, 0]; % Stochiometric coefficients [Ca, Mg, K, Na, Al, Si, P]
        K_sp = 10^(9.83);
    elseif strcmp(mineral, 'labradorite') % Na0.4Ca0.6Al1.6Si2.4O8 (http://webmineral.com/data/Labradorite.shtml)
        MM_min = 272/conv_mol; % [g/mol-conv]
        k_diss_H = 10^(-7.87)*24*3600*conv_mol; % [mol-conv m-2 d-1]
        k_diss_w = 10^(-10.91)*24*3600*conv_mol; % [mol-conv m-2 d-1]
        k_diss_OH = 0;
        E_H = 42/conv_mol; % [kJ/mol-conv]: activation energies
        E_w = 45/conv_mol; % [kJ/mol-conv]
        E_OH = 0;
        n_H = 0.6; % reaction order
        n_OH = 1;
        min_st = [0.6, 0, 0, 0.4, 1.6, 2.4, 0]; % Stochiometric coefficients [Ca, Mg, K, Na, Al, Si, P]
        K_sp = NaN;
    elseif strcmp(mineral, 'augite') % Ca0.9Na0.1Mg0.9Fe0.2Al0.4Ti0.1Si1.9O6 (http://webmineral.com)
        MM_min = 236/conv_mol; % [g/mol-conv]
        k_diss_H = 10^(-6.82)*24*3600*conv_mol; % [mol-conv m-2 d-1]
        k_diss_w = 10^(-11.97)*24*3600*conv_mol; % [mol-conv m-2 d-1]
        k_diss_OH = 0;
        E_H = 78/conv_mol; % [kJ/mol-conv]: activation energies
        E_w = 78/conv_mol; % [kJ/mol-conv]
        E_OH = 0;
        n_H = 0.7; % reaction order
        n_OH = 1;
        min_st = [0.9, 0.9, 0, 0.1, 0.4, 1.9, 0]; % Stochiometric coefficients [Ca, Mg, K, Na, Al, Si, P]
        K_sp = NaN;
    elseif strcmp(mineral, 'diopside') % MgCaSi2O6
        MM_min = 216/conv_mol; % [g/mol-conv]
        k_diss_H = 10^(-6.36)*24*3600*conv_mol; % [mol-conv m-2 d-1]
        k_diss_w = 10^(-11.11)*24*3600*conv_mol; % [mol-conv m-2 d-1]
        k_diss_OH = 0;
        E_H = 96/conv_mol; % [kJ/mol-conv]: activation energies
        E_w = 40/conv_mol; % [kJ/mol-conv]
        E_OH = 0;
        n_H = 0.71; % reaction order
        n_OH = 1;
        min_st = [1, 1, 0, 0, 0, 2, 0]; % Stochiometric coefficients [Ca, Mg, K, Na, Al, Si, P]
        K_sp = 10^(5.30);
    elseif strcmp(mineral, 'apatite') % Ca5(PO4)3(OH)
        MM_min = 502/conv_mol/5; % [g/mol-conv]
        k_diss_H = 10^(-4.29)*24*3600*conv_mol; % [mol-conv m-2 d-1]
        k_diss_w = 10^(-6)*24*3600*conv_mol;
        k_diss_OH = 0;
        E_H = 250/conv_mol; % [kJ/mol-conv]: activation energies
        E_w = 250/conv_mol;
        E_OH = 0;
        n_H = 0.17; % reaction order
        n_OH = 1;
        min_st = [5, 0, 0, 0, 0, 0, 3]/5; % Stochiometric coefficients [Ca, Mg, K, Na, Al, Si, P]
        K_sp = 10^(14.34); % https://thermoddem.brgm.fr/species/hydroxyapatitenatur
    elseif strcmp(mineral, 'chlorite') % (Mg0.63Fe0.37)5Al(Si3Al1.2)O10 (OH)8 - Lewis et al. (2021)
        MM_min = 620/conv_mol; % [g/mol-conv]
        k_diss_H = 10^(-11.11)*24*3600*conv_mol; % [mol-conv m-2 d-1]
        k_diss_w = 10^(-12.52)*24*3600*conv_mol;
        k_diss_OH = 0;
        E_H = 88/conv_mol; % [kJ/mol-conv]: activation energies
        E_w = 88/conv_mol;
        E_OH = 0;
        n_H = 0.5; % reaction order
        n_OH = 1;
        min_st = [0, 3.15, 0, 0, 2.2, 3, 0]; % Stochiometric coefficients [Ca, Mg, K, Na, Al, Si, P]]
        K_sp = 10^(61.34); %https://thermoddem.brgm.fr/species/chloritecca-2
    elseif strcmp(mineral, 'epidote') % Ca2(Al2.53Fe1.09)SiO4(SiO7) O(OH) - Lewis et al. (2021)
        MM_min = 474/conv_mol; % [g/mol-conv]
        k_diss_H = 10^(-10.6)*24*3600*conv_mol; % [mol-conv m-2 d-1]
        k_diss_w = 10^(-11.99)*24*3600*conv_mol;
        k_diss_OH = 10^(-17.33)*24*3600*conv_mol;
        E_H = 71.1/conv_mol; % [kJ/mol-conv]: activation energies
        E_w = 70.7/conv_mol;
        E_OH = 79.1/conv_mol;
        n_H = 0.338; % reaction order
        n_OH = -0.556;
        min_st = [2, 0, 0, 0, 2.53, 2, 0]; % Stochiometric coefficients [Ca, Mg, K, Na, Al, Si, P]]
        K_sp = 10^(32.23); % https://thermoddem.brgm.fr/species/epidote
    elseif strcmp(mineral, 'zoisite') % zoisite Ca2Al3Si3O12(OH)
        MM_min = 454.357/conv_mol; % [g/mol-conv]
        k_diss_H = 10^(-7.50)*24*3600*conv_mol; % [mol-conv m-2 d-1]
        k_diss_w = 10^(-11.2)*24*3600*conv_mol;
        k_diss_OH = 0;
        E_H = 66.1/conv_mol; % [kJ/mol-conv]: activation energies
        E_w = 66.1/conv_mol;
        E_OH = 0;
        n_H = 0.5; % reaction order
        n_OH = 1;
        min_st = [2, 0, 0, 0, 3, 3, 0]; % Stochiometric coefficients [Ca, Mg, K, Na, Al, Si, P]]
        K_sp = 10^(43.85); % https://thermoddem.brgm.fr/species/epidote
    elseif strcmp(mineral, 'ferroactinolite') % Ca2(Mg0.75Fe0.25)5Si8.5O22(OH)2 - Amphibole  - Lewis et al. (2021)
        MM_min = 865.8/conv_mol; % [g/mol-conv]
        k_diss_H = 10^(-8.4)*24*3600*conv_mol; % [mol-conv m-2 d-1]
        k_diss_w = 10^(-10.6)*24*3600*conv_mol;
        k_diss_OH = 0;
        E_H = 18.9/conv_mol; % [kJ/mol-conv]: activation energies
        E_w = 94.4/conv_mol;
        E_OH = 0;
        n_H = 0.7; % reaction order
        n_OH = 1;
        min_st = [2, 3.75, 0, 0, 0, 8.5, 0]; % Stochiometric coefficients [Ca, Mg, K, Na, Al, Si, P]]
        K_sp = 10^(7.3); % https://thermoddem.brgm.fr/species/actinolite
    elseif strcmp(mineral, 'quartz') % SiO2 - Lewis et al. (2021)
        MM_min = 60/conv_mol; % [g/mol-conv]
        k_diss_H = 10^(-11.36)*24*3600*conv_mol; % [mol-conv m-2 d-1]
        k_diss_w = 10^(-13.4)*24*3600*conv_mol; % [mol-conv m-2 d-1]
        k_diss_OH = 10^(-16.29)*24*3600*conv_mol;
        E_H = 90.9/conv_mol; % [kJ/mol-conv]: activation energies
        E_w = 90.9/conv_mol; % [kJ/mol-conv] % Yoshiki Kanzaki et al.(2022)
        E_OH = 90.9/conv_mol;
        n_H = 0.309; % reaction order
        n_OH = -0.411;
        min_st = [0, 0, 0, 0, 0, 1, 0]; % Stochiometric coefficients [Ca, Mg, K, Na, Al, Si, P]
        K_sp = 10^(-4); 
    elseif strcmp(mineral, 'actinolite') % Ca2(Mg2.25Fe2.5Al0.25)(Si7.75Al0.25)O22(OH)2 - https://thermoddem.brgm.fr/species/actinolite
        MM_min = 892/conv_mol; % [g/mol-conv]
        k_diss_H = 10^(-8.4)*24*3600*conv_mol; % [mol-conv m-2 d-1]
        k_diss_w = 10^(-10.6)*24*3600*conv_mol;
        k_diss_OH = 0;
        E_H = 18.9/conv_mol; % [kJ/mol-conv]: activation energies
        E_w = 94.4/conv_mol;
        E_OH = 0;
        n_H = 0.7; % reaction order
        n_OH = 1;
        min_st = [2, 2.25, 0, 0, 0.5, 7.75, 0]; % Stochiometric coefficients [Ca, Mg, K, Na, Al, Si, P]
        K_sp = 10^(7.13); 
    elseif strcmp(mineral, 'calcite') % CaCO3
        MM_min = 100/conv_mol; % [g/mol-conv]
        k_diss_H = 10^(-0.3)*24*3600*conv_mol; % [mol-conv m-2 d-1]
        k_diss_w = 10^(-5.81)*24*3600*conv_mol;
        k_diss_OH = 10^(-3.48)*24*3600*conv_mol;
        E_H = 14.4/conv_mol; % [kJ/mol-conv]: activation energies
        E_w = 23.5/conv_mol;
        E_OH = 35.4/conv_mol;
        n_H = 1; % reaction order
        n_OH = 1;
        min_st = [1, 0, 0, 0, 0, 0, 0]; % Stochiometric coefficients [Ca, Mg, K, Na, Al, Si, P]
        K_sp = 10^(-8.35); 
    else
        error('No data for this mineral');
    end

end

function [Omega] = Omega_sil(mineral, K, Ca, Mg, Na, Si, Al, PO4, H, K_sp, conv_mol)
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
    elseif ismember(mineral, {'labradorite', 'augite', 'calcite'})
        Omega = 0;
    else
        error('Unknown mineral');
    end

end

function [CEC, diag] = Resolve_CEC(CEC, Pcla, Porg)
%%% Use a positive supplied CEC as measured.  Otherwise use the
%%% texture/organic-C pedotransfer estimate.
diag = struct('source', "explicit_CEC", 'estimated', false, ...
    'measured', true, 'texture_organic_C_ptf_cmolc_kg', NaN);
valid_soil_fractions = isscalar(Pcla) && isfinite(Pcla) && Pcla >= 0 && Pcla <= 1 && ...
    isscalar(Porg) && isfinite(Porg) && Porg >= 0 && Porg <= 1 && Pcla + Porg < 1;
if valid_soil_fractions
    clay_pct = 100*Pcla;
    organic_C_pct = 100*Porg/1.724;
    diag.texture_organic_C_ptf_cmolc_kg = 0.7*clay_pct + 3.5*organic_C_pct;
end
if ~(isscalar(CEC) && isfinite(CEC) && CEC > 0)
    if ~valid_soil_fractions
        error('ERW:InvalidCEC', ...
            'Missing CEC requires valid clay and organic-matter mass fractions.');
    end
    CEC = diag.texture_organic_C_ptf_cmolc_kg;
    if ~(isfinite(CEC) && CEC > 0)
        error('ERW:InvalidCEC', ...
            'Texture and organic matter did not give a positive CEC estimate.');
    end
    diag.source = "texture_organic_C_ptf";
    diag.estimated = true;
    diag.measured = false;
end
diag.CEC_cmolc_kg = CEC;
end
