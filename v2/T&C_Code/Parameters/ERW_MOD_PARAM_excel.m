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
P_weathering_release = zeros(length(t), 1);

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
AlOH3_solid = zeros(length(t), 1);
Omega_AlOH3 = zeros(length(t), 1);
Al_inventory_residual = nan(length(t), 1);
charge_residual = nan(length(t), 1);
dAlk_sil = zeros(length(t), 1);
dAlk_carb = zeros(length(t), 1);
dAlk_check = zeros(length(t), 1);
dAlk_total = zeros(length(t), 1);
gamma_i = nan(length(t), 5);
pH_conc = zeros(length(t), 1);
ionic_str = zeros(length(t), 1);
activity_I = zeros(length(t), 1);
gamma_mono = zeros(length(t), 1);
pCO2_atm = zeros(length(t), 1);
activity_passes = zeros(length(t), 1);

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
M_rock_add_CaCO3 = 0;
if isfield(ERW_TC_const_input,'M_rock_add_CaCO3')
    M_rock_add_CaCO3 = ERW_TC_const_input.M_rock_add_CaCO3;
end

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
    if any(~isfinite(Na_rain_src)) || any(Na_rain_src < 0)
        error('ERW:invalidBoundary', 'Rainwater Na supply must be finite and nonnegative.');
    end
    I_Na_rain_mean = mean(Na_rain_src);
    mean_water_loss_Lperm2d = ERW_Seed_Water_Loss(ERW_TC_const_input, TT_m_input);
    if I_Na_rain_mean > 0
        if ~(isfinite(mean_water_loss_Lperm2d) && mean_water_loss_Lperm2d > 0)
            error('ERW:invalidBoundary', ...
                'A positive water-loss flux is required for the steady rainwater Na seed.');
        end
        Na_rain_ss = I_Na_rain_mean / mean_water_loss_Lperm2d;
    else
        Na_rain_ss = 0;
    end
    Na_in = Na_rain_ss;
end
if ~(isfinite(Na_in) && Na_in >= 0)
    error('ERW:invalidBoundary', 'Seed dissolved Na must be finite and nonnegative.');
end

if ~exist('mineral','var') || isempty(mineral) || ...
        (isnumeric(mineral) && all(isnan(mineral(:))))
    rock_f_in = NaN;
end
%% soil properties
Zrw = Zbio/1000; % Depth of the active Biogeochemistry zone which weathering takes place in [m]
rho_bulk = rsd_bio*1e3; % soil dry mass bulk density [g/ m3]

%% initial conditions

% CEC 
[CEC_tot_in, cec_diag] = ERW_Resolve_CEC(CEC_tot_in);
CEC_tot = CEC_tot_in*1e-5*rho_bulk*Zrw*conv_mol*ones(length(t), 1); % [umol_c] 

% pH
pH(1)= pH_in;

% Carbonate minerals (considered as an additional pool)
CaCO3(1) = CaCO3_in + M_rock_add_CaCO3(1); %[umol-conv]
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
use_amorphous_AlOH3 = isfield(ERW_TC_const_input, 'bulk_pH_params') && ...
    isstruct(ERW_TC_const_input.bulk_pH_params) && ...
    isfield(ERW_TC_const_input.bulk_pH_params, 'Al_solubility_pin') && ...
    isequal(ERW_TC_const_input.bulk_pH_params.Al_solubility_pin, true);
K_AlOH3 = ERW_AlOH3_LogK(use_amorphous_AlOH3);
% When the Al solubility pin is active the daily solver works with the
% amorphous phase (ERW_function overrides K_AlOH3 to 10.8, SCEPTER
% 'amal'), so the seed and its saturation diagnostic must use the same
% phase.  Otherwise Al_w(1) is reported against microcrystalline
% gibbsite and Omega_AlOH3(1) comes out at 10^(10.8-9.66) = 13.8, which
% reads as supersaturation when the solution is exactly at equilibrium.

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
        error('ERW:invalidParam', 'diss_f must be finite for an active amendment.');
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
% Methods: 1:1 water:soil (n=1), 1:2 DI water (n=2), 1:2.5 (n=2.5), 1:5 (n=5).
% Beerling et al. (2023 preprint, arXiv 2307.05343), Soil analyses: soil pH
% measured in a 1:2.5 soil-to-water paste - 10 g air-dried soil sieved <2 mm
% with 25 mL ultrapure water, vortexed 10 s.  The 2024 PNAS
% version cites the method by reference only; its 10 g / 25 mL figure appears
% in the sodium-phosphate Si extraction, which is a different assay.
if exist('pH_DilutionFactor','var') && isnumeric(pH_DilutionFactor) && ...
        isfinite(pH_DilutionFactor) && pH_DilutionFactor > 0
    bulk_pH_params.n_extract = pH_DilutionFactor;
end
if isfield(bulk_pH_params,'invert_initial_bulk_pH')
    error('ERW:deprecatedInitialPHMode', ...
        'Use bulk_pH_params.initial_pH_mode instead of invert_initial_bulk_pH.');
end
if ~isfield(bulk_pH_params,'initial_pH_mode')
    bulk_pH_params.initial_pH_mode = 'state_prescribed';
end
initial_pH_mode = lower(string(bulk_pH_params.initial_pH_mode));
if ~isscalar(initial_pH_mode) || ...
        ~ismember(initial_pH_mode,["state_prescribed","bulk_observed"])
    error('ERW:invalidInitialPHMode', ...
        'initial_pH_mode must be state_prescribed or bulk_observed.');
end
bulk_pH_params.initial_pH_mode = char(initial_pH_mode);
if initial_pH_mode == "state_prescribed"
    bulk_pH_params.prescribed_pore_pH_initial = pH_in;
else
    bulk_pH_params.observed_bulk_pH_initial = pH_in;
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

M_iner = 0; % [g/m2] no inert inventory before an initial or scheduled rock addition
if active_rock_input(rock_f_in)

number_min = length(mineral);
if numel(rock_f_in) ~= number_min
    error('ERW:mineralFractionSize', ...
        'rock_f_in must contain one mass fraction per primary mineral.');
end
    rock_f_in = rock_f_in(:).';
    psd_perc_in = psd_in; % dimensionless particle-class fractions for scheduled additions
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

    psd_perc_in = 0;
    rock_f_in = 0;
    d_in = 0;
    SSA_in = 0; 
    a = 0;

end



%% other parameters
% temperature and gas constants 
R = 8.314; % [J mol-1 K-1]: universal gas constant
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
K_CaCO3 = 10^(-171.9065-0.077993*T_K(1)+2839.319/T_K(1)+71.595*log10(T_K(1)))*conv_mol^2; % Calcite [mol^2], Plummer & Busenberg (1982)
K_MgCO3 = 10^(-5.27)*exp(17170/R*(1/T_K(1)-1/298.15))*conv_mol^2; % Nesquehonite [mol^2], Harrison et al. (2019), van't Hoff with dH_r = -17.17 kJ/mol

% precipitation rate
% https://nora.nerc.ac.uk/id/eprint/511084/1/Kirk%20et%20al%202015%20Geochmica%20et%20Cosmochimica%20Acta.pdf
r_CaCO3 = 5e-10*1e3*(24*3600)*conv_mol; % [umol m^-3 soil d^-1]
r_MgCO3 = r_CaCO3; % nesquehonite; no soil rate is published and it precipitates at equilibrium within hours (Harrison et al., 2013)



% rainwater  
CO2_w_rain(1) = k_H(1)*CO2_atm(1); % [mol/l] Henry's law

%pCO2 - partial pressure of CO2 [Pa]
%%% For uniform production q over a layer of depth L with a no-flux base,
%%% D*C'' = -q gives an excess over the atmosphere of q*L^2/(2D) at the base
%%% and q*L^2/(3D) as the layer mean.  The box flux q*L = D*(C-C_atm)/Z_CO2
%%% with Z_CO2 = Zrw/2 therefore carries the BASE excess, 1.5 times the
%%% layer-mean excess (Z_CO2 = Zrw/3 would carry the mean).
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
if string(excel_input.fX_conc_convention) == "Kantola"
    bulk_pH_params.initial_anion_convention = 'full_charge';
end

% The initialization boundary is explicit: either retain a prescribed field
% state or invert a measured bulk-extract pH.  The two observations are never
% applied to the same hydrogen state.
if ~isfield(bulk_pH_params,'n_extract') || ...
        ~(isfinite(bulk_pH_params.n_extract) && bulk_pH_params.n_extract > 0)
    error('ERW:missingParam', ...
        'Observed bulk pH initialization requires a positive extraction ratio.');
end
%%% Activity closure of the seed.  The seed has to satisfy the conditional
%%% constants (ERW_Conditional_Constants) at its OWN ionic strength, which
%%% depends on the seed, so the inversion and speciation are repeated from
%%% the same inputs until that ionic strength is stationary.  An observed
%%% bulk pH_in is an electrode (activity) pH.  The exported k1, k2, k_w,
%%% K1..K4 and K_AlOH3 stay thermodynamic: the daily solver and the extract
%%% operator condition them at their own ionic strength.
seed_Ca = Ca(1); seed_Mg = Mg(1); seed_K = K(1); seed_Na = Na(1);
seed_An = An(1); seed_Al_w = Al_w(1);
seed_f = [f_Ca(1), f_Mg(1), f_K(1), f_Na(1), f_Al(1), f_H(1)];
seed_bulk_pH_params = bulk_pH_params;
ionic_str_seed = 0;
ionic_str_next = NaN;
seed_closed = false;
for seed_pass = 1:30
    Ca(1) = seed_Ca; Mg(1) = seed_Mg; K(1) = seed_K; Na(1) = seed_Na;
    An(1) = seed_An; Al_w(1) = seed_Al_w;
    f_Ca(1) = seed_f(1); f_Mg(1) = seed_f(2); f_K(1) = seed_f(3);
    f_Na(1) = seed_f(4); f_Al(1) = seed_f(5); f_H(1) = seed_f(6);
    bulk_pH_params = seed_bulk_pH_params;
    [k_seed,K_AlOH3_s,gamma_seed] = ERW_Conditional_Constants( ...
        [k1(1),k2(1),k_w(1),K1,K2,K3,K4],K_AlOH3,ionic_str_seed);
    k1_s = k_seed(1); k2_s = k_seed(2); k_w_s = k_seed(3);
    K1_s = k_seed(4); K2_s = k_seed(5); K3_s = k_seed(6); K4_s = k_seed(7);
    gamma_mono_0 = gamma_seed(1);
    gamma_div_0 = gamma_seed(2);
    if initial_pH_mode == "bulk_observed"
        bulk_pH_params.pore_ionic_strength = ionic_str_seed;
        [H(1), bulk_pH_params, f_Ca(1), f_Mg(1), f_K(1), f_Na(1), f_Al(1), f_H(1), ...
            An(1), Al_w(1), Ca(1), Mg(1), K(1), Na(1)] = Invert_Bulk_pH( ...
            pH_in, Ca(1), Mg(1), K(1), Na(1), NO3(1), PO4(1), NH4(1), ...
            CO2_atm(1), CO2_w(1), k1(1), k2(1), k_w(1), k_nh3(1), k_H(1), ...
            K1, K2, K3, K4, CEC_tot(1), rho_bulk, Zrw, conv_mol, ...
            bulk_pH_params, K_AlOH3, Al_CEC_scale);
    end
    %%% A prescribed pore-water state keeps its hydrogen concentration: the
    %%% pre-processing seed (exchange fractions, Al) was derived at it.  Only
    %%% its reported pH moves to the activity scale.

    % Sodium needs no re-baseline here.  The Kantola seed
    % (ERW_pre_processing) sets pore-water Na from deposition and derives
    % its exchange fraction from the Appelo law, so the inventory handed
    % to the inversion is already the one the exchanger can host.

    pH_conc(1) = -log10(H(1)/conv_mol);
    pH(1) = pH_conc(1)-log10(gamma_mono_0);
    if isfield(bulk_pH_params,'dCEC_dpH') && ...
            (~isfield(bulk_pH_params,'CEC_ref_pH') || isempty(bulk_pH_params.CEC_ref_pH))
        cec_probe = bulk_pH_params;
        cec_probe.CEC_ref_pH = pH_conc(1);
        cec_scheme = 'TC_GT';
        if isfield(cec_probe,'scheme') && ~isempty(cec_probe.scheme)
            cec_scheme = char(string(cec_probe.scheme));
        end
        if ERW_Resolve_CEC_pH_Sensitivity(cec_probe,cec_scheme) > 0
            bulk_pH_params.CEC_ref_pH = pH_conc(1);
        end
        clear cec_probe cec_scheme
    end
    Al_den_seed = H(1)^4 + H(1)^3*K1_s + H(1)^2*K1_s*K2_s + ...
        H(1)*K1_s*K2_s*K3_s + K1_s*K2_s*K3_s*K4_s;
    Al(1) = (H(1)^4/Al_den_seed) * Al_w(1);
    f_CEC(1,:) = [f_Ca(1), f_Mg(1), f_K(1), f_Na(1), f_Al(1), f_H(1)];
    bulk_pH_params.H_init = NaN;
    bulk_pH_params.H_extract_anchor_umolL = NaN;
    bulk_pH_params.H_pore_anchor_umolL = H(1);

    %carbonate system
    HCO3(1)= k1_s*CO2_w(1)/H(1); % bicarbonate [umol/l]
    CO3(1)= k2_s*k1_s*CO2_w(1)/(H(1)^2); % carbonate [umol/l]
    DIC(1)= HCO3(1)+CO3(1)+CO2_w(1); % dissolved inorganic C (includes bicarbonate, carbonate and dissolved CO2 (carbonic acid))
    IC_tot(1)=(DIC(1)*s(1)+CO2_air(1)*(1-s(1)))*(n(1)*Zrw*1000);% total inorganic C [umol]

    %Alkalinity
    Alk(1)= HCO3(1)+2*CO3(1)-H(1)+k_w_s/H(1); % [umol/L]

    % Al from Ca-Al GT exchange at initial pore-water pH
    % If exchange inputs fail, use Al(OH)3 saturation.
    if isnan(Al_w_in) && ~(isfinite(Al_w(1)) && Al_w(1) >= 0)
        Al_w_eq = Initialize_ERW_Equilibrium(H(1), Ca(1), f_Ca(1), f_Al(1), ...
            K_Ca_Al, K1_s, K2_s, K3_s, K4_s, Al_CEC_scale);
        if isnan(Al_w_eq)
            Al_w_eq = Initialize_ERW_Equilibrium(H(1), K1_s, K2_s, K3_s, K4_s, ...
                conv_mol, K_AlOH3_s, Al_CEC_scale);
        end
        if ~isfinite(Al_w_eq) || Al_w_eq <= 0
            error('Cannot initialize Al_w from exchange or Al(OH)3 equilibrium.');
        end
        Al_w(1) = Al_w_eq;
    end

    % aluminium speciation
    Al(1)=(H(1)^4/(H(1)^4+H(1)^3*K1_s+H(1)^2*K1_s*K2_s+H(1)*K1_s*K2_s*K3_s+K1_s*K2_s*K3_s*K4_s))*Al_w(1);% [mol/l]
    AlOH(1)=(H(1)^3*K1_s/(H(1)^4+H(1)^3*K1_s+H(1)^2*K1_s* K2_s+H(1)*K1_s*K2_s*K3_s+K1_s*K2_s*K3_s*K4_s))* Al_w(1);
    AlOH2(1)=(H(1)^2*K1_s*K2_s/(H(1)^4+H(1)^3*K1_s+H(1)^2*K1_s*K2_s+H(1)*K1_s*K2_s*K3_s+K1_s*K2_s*K3_s*K4_s))*Al_w(1);
    AlOH3(1)=(H(1)*K1_s*K2_s*K3_s/(H(1)^4+H(1)^3*K1_s+H(1)^2*K1_s*K2_s+H(1)*K1_s*K2_s*K3_s+K1_s*K2_s*K3_s*K4_s))*Al_w(1);
    AlOH4(1)= Al_w(1)-(Al(1)+AlOH(1)+AlOH2(1)+AlOH3(1));

    Al_w_sat(1) = Initialize_ERW_Equilibrium(H(1), K1_s, K2_s, K3_s, K4_s, conv_mol, K_AlOH3_s, Al_CEC_scale);
    Omega_AlOH3(1) = Al_w(1) / max(Al_w_sat(1), eps);

    % Initial An is the unmodelled strong-anion countercharge returned by the
    % observed-pH inverse.  The deposition budget below supplies future fluxes only.
    %
    if string(excel_input.fX_conc_convention) == "Kantola" && ...
            initial_pH_mode == "state_prescribed"
        An(1) = ERW_Residual_An_From_Charge(H(1),CO2_w(1),Ca(1),Mg(1), ...
            Na(1),K(1),NH4(1),NO3(1),PO4(1),Al_w(1),k1_s,k2_s, ...
            k_w_s,k_nh3(1),K1_s,K2_s,K3_s,K4_s,Al_CEC_scale);
        if ~(isfinite(An(1)) && An(1) >= 0)
            error('ERW:invalidSeed', ...
                'The Kantola seed has a negative residual anion.');
        end
    end

    ionic_str_next = ERW_Ionic_Strength(Ca(1),Mg(1),Na(1),K(1),Al(1)/Al_CEC_scale, ...
        H(1),NH4(1),NO3(1),HCO3(1),CO3(1),An(1),conv_mol,PO4(1), ...
        AlOH(1)/Al_CEC_scale,AlOH2(1)/Al_CEC_scale,AlOH4(1)/Al_CEC_scale,k_w_s/H(1));
    if abs(ionic_str_next-ionic_str_seed) <= 1e-6*max(ionic_str_next,1e-6)
        seed_closed = true;
        break
    end
    ionic_str_seed = ionic_str_next;
end
if ~seed_closed
    error('ERW:seedActivityClosure', ...
        'The seed ionic strength did not close in %d passes (%.6g -> %.6g mol/L).', ...
        seed_pass,ionic_str_seed,ionic_str_next);
end
if isfield(bulk_pH_params,'pore_ionic_strength')
    bulk_pH_params = rmfield(bulk_pH_params,'pore_ionic_strength');
end
ionic_str(1) = ionic_str_next;
activity_I(1) = ionic_str_seed;
gamma_mono(1) = gamma_mono_0;
activity_passes(1) = seed_pass;
%%% The record carries the conditioned constants the seed was speciated
%%% with, as ERW_function records them at every later step.
k1(1) = k1_s;
k2(1) = k2_s;
k_w(1) = k_w_s;

R_alk(1)=(f_Mg(1)+f_Ca(1)+f_Na(1)+f_K(1))*CEC_tot(1); % [umol_c]
Ca_tot(1) = Ca(1)*n(1)*s(1)*Zrw*1000+f_Ca(1)/2*CEC_tot(1);
Mg_tot(1) = Mg(1)*n(1)*s(1)*Zrw*1000+f_Mg(1)/2 *CEC_tot(1);
K_tot(1) = K(1)*n(1)*s(1)*Zrw*1000+f_K(1)*CEC_tot(1);
Si_tot(1) = Si(1)*n(1)*Zrw*s(1)*1000;
Na_tot(1) = Na(1)*n(1)*s(1)*Zrw*1000+ f_Na(1)* CEC_tot(1);
Al_tot(1) = Al_w(1)*n(1)*Zrw*s(1)*1000+(f_Al(1)/3)*CEC_tot(1)*Al_CEC_scale;
An_tot(1) = An(1)*n(1)*s(1)*Zrw*1000;

soil_water_l_0 = n(1)*Zrw*s(1)*1000;
[~,AlOH3_solid(1),Al_inventory_residual(1)] = ERW_Al_Equilibrium_Target( ...
    H(1),f_Al(1),Al_tot(1),soil_water_l_0,CEC_tot(1),Al_CEC_scale, ...
    K_AlOH3_s,K1_s,K2_s,K3_s,K4_s);
charge_residual(1) = H(1)+2*Ca(1)+2*Mg(1)+Na(1)+K(1)+NH4(1)+ ...
    (3*Al(1)+2*AlOH(1)+AlOH2(1)-AlOH4(1))/Al_CEC_scale- ...
    (k_w_s/H(1)+HCO3(1)+2*CO3(1)+NO3(1)+PO4(1)+An(1));
if isfield(bulk_pH_params,'scheme') && strcmpi(bulk_pH_params.scheme,'appelo')
    appelo_seed = [];
    if isfield(bulk_pH_params,'appelo')
        appelo_seed = bulk_pH_params.appelo;
    end
    [~,~,~,~,~,appelo_alpha_seed] = Resolve_Appelo_Constants(appelo_seed);
    gamma_i(1,:) = 10.^(appelo_alpha_seed.*[1 1 2 2 3].*f_H(1));
end

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
if any(~isfinite(I_Na_source(:)) | I_Na_source(:) < 0)
    error('ERW:invalidBoundarySource', ...
        'Sodium boundary forcing must be finite and nonnegative.');
end
I_Na = first_finite_source(I_Na_source);

%%% No site-constrained dissolved Ca or Mg boundary supply is available.
I_Ca = 0;
I_Mg = 0;
I_K = 0;
I_Si = 0;

%Carbonate weathering
% The seed's HCO3 and CO3 come from the conditional constants, so they are
% concentrations like Ca and Mg and the carbonate ion activity product is
% gamma_div^2*Ca*CO3, at the ionic strength the seed closed at -- the
% same basis ERW_function uses at every later step.
pCO2_atm(1) = CO2_air(1)*22.41/conv_mol;
gamma_IAP_0 = gamma_div_0^2;
Omega_CaCO3(1) = gamma_IAP_0*Ca(1)*CO3(1)/K_CaCO3;% carbonate saturation index, activity basis [-]
Omega_MgCO3(1) = gamma_IAP_0*Mg(1)*CO3(1)/K_MgCO3;

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
    M_rock(1) = max(M_rock_in,M_rock_add(1)); % rock mass [g/m2]
    rock_f(1,:)= rock_f_in;
    M_min(1,:)= rock_f(1,:)*M_rock(1); % mineral mass [g/m2]
    M_iner(1) = M_rock(1)*(1-sum(rock_f(1,:))); % inert rock mass??[g/m2]
    
    %particle sediment diameter
    delta_d_in = [d_in(1,1) diff(d_in(1,:))];
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
            if strcmpi(mineral(min_j),'calcite')
                [Omega(1,min_j),Wr(1,min_j),EW(1,min_j)] = ERW_Calcite_Kinetics( ...
                    Ca(1),CO3(1),K_CaCO3,gamma_mono_0*H(1),CO2_air(1),conv_mol,s(1),diss_f, ...
                    k_diss_H_t(1,min_j),k_diss_w_t(1,min_j),k_diss_OH_t(1,min_j), ...
                    n_H(min_j),n_OH(min_j),SA(1),rock_f(1,min_j), ...
                    M_min(1,min_j),MM_min(min_j),dtw_init,gamma_IAP_0);
            else
                [Omega(1,min_j),Omega_defined] = ERW_Mineral_Saturation(mineral(min_j), K(1), Ca(1), Mg(1), Na(1), Si(1), Al(1), PO4(1), H(1), K_sp(min_j), conv_mol, Al_CEC_scale, activity_I(1));
                if Omega_defined
                    dissolution_affinity = max(0,1-Omega(1,min_j));
                else
                    dissolution_affinity = 1; % missing aqueous products: far-from-equilibrium kinetic assumption
                end
                Wr(1, min_j) = s(1).*diss_f.*(k_diss_H_t(1,min_j).*(gamma_mono_0*H(1)/conv_mol).^n_H(min_j)+k_diss_w_t(1,min_j)+k_diss_OH_t(1,min_j).*(gamma_mono_0*H(1)/conv_mol).^n_OH(min_j)).*dissolution_affinity; % background weathering rate [mol-conv/ m2 d]
                EW(1, min_j) = Wr(1, min_j).*SA(1).*rock_f(1, min_j); % enhanced weathering rate [umol/d]
                EW_inventory_max = M_min(1,min_j)/MM_min(min_j)/dtw_init;
                if EW(1,min_j) > EW_inventory_max
                    EW(1,min_j) = EW_inventory_max;
                    area_fraction = SA(1)*rock_f(1,min_j);
                    if area_fraction > 0
                        Wr(1,min_j) = EW(1,min_j)/area_fraction;
                    else
                        Wr(1,min_j) = 0;
                    end
                end
            end
        end

    end

end

calcite_j = find(strcmpi(string(mineral),'calcite'),1);
EW_calcite = 0;
if active_rock_input(rock_f_in) && num_rock_app > 0 && ~isempty(calcite_j)
    EW_calcite = EW(1,calcite_j);
end
[EW_calcite,W_CaCO3(1),calcite_limit] = ERW_Calcite_Reaction_Limit( ...
    EW_calcite,W_CaCO3(1),dtw_init,Omega_CaCO3(1),Ca(1),DIC(1),Alk(1), ...
    n(1)*Zrw*s(1)*1000,k1_s,k2_s,k_w_s,K_CaCO3,conv_mol,gamma_IAP_0);
[~,W_MgCO3(1)] = ERW_Calcite_Reaction_Limit( ...
    0,W_MgCO3(1),dtw_init,Omega_MgCO3(1),Mg(1),calcite_limit.DIC_final,calcite_limit.Alk_final, ...
    n(1)*Zrw*s(1)*1000,k1_s,k2_s,k_w_s,K_MgCO3,conv_mol,gamma_IAP_0);
if active_rock_input(rock_f_in) && num_rock_app > 0 && ~isempty(calcite_j)
    EW(1,calcite_j) = EW_calcite;
    area_fraction = SA(1)*rock_f(1,calcite_j);
    if area_fraction > 0
        Wr(1,calcite_j) = EW_calcite/area_fraction;
    else
        Wr(1,calcite_j) = 0;
    end
end

%%
% structure containing biogeochemical variables
ERW_biogeochem_varnames = {'pH';'CO2_air';'CO2_atm';'Ca';'Mg';'K';'Na';'An';'Al';'AlOH';'AlOH2';'AlOH3';'AlOH4';'Si';'PO4';'NO3';'NH4';'H';
    'Ca_tot';'Mg_tot';'K_tot';'Na_tot';'An_tot';'Al_tot';'Si_tot';'PO4_tot';'NO3_tot';'NH4_tot';'IC_tot';'CEC_tot';'K_CEC';
    'R_alk';'Fs';'DIC';'CaCO3';'MgCO3';'f_Ca';'f_Mg';'f_K';'f_Na';'f_Al';'f_H';'f_CEC';'M_min';'CO2_w';'Alk';'Alk_tot';'k_H';'ADV';'DIC_rain';'D';'CO3';...
    'Al_w';'Al_w_sat';'AlOH3_solid';'Al_inventory_residual';'charge_residual';'dAlk_sil';'dAlk_carb';'dAlk_check';'dAlk_total';'gamma_i';'HCO3';'Leak_Na';'Leak_Al';'Leak_An';'Leak_DIC';'Na_uptake';'An_uptake';...
    'k1';'k2';'k_w';'k_nh3';'pH_conc';'ionic_str';'activity_I';'gamma_mono';'pCO2_atm';'activity_passes';'Omega_CaCO3'};

ERW_biogeochem = struct();
for var_index = 1:size(ERW_biogeochem_varnames,1)
    vn = cell2mat(ERW_biogeochem_varnames(var_index));
    if exist(vn, 'var')
        ERW_biogeochem.(vn) = eval(vn);
    end
end
bulk_pH_params.pH_offset = 0;
bulk_pH_params.DIC_field_umol_L = DIC(1);
[pH_bulk_raw,bulk_pH_diag] = Porewater_to_Bulk_pH( ...
    H(1), Ca(1), Mg(1), K(1), Na(1), Al_w(1), Al(1), An(1), ...
    NO3(1), PO4(1), NH4(1), CO2_atm(1), k1(1), k2(1), k_w(1), ...
    k_nh3(1), k_H(1), K1, K2, K3, K4, ...
    f_Ca(1), f_Mg(1), f_K(1), f_Na(1), f_Al(1), f_H(1), ...
    CEC_tot(1), rho_bulk, Zrw, conv_mol, bulk_pH_params, Al_CEC_scale);
bulk_pH_params.H_prev = bulk_pH_diag.H;
bulk_pH_params.H_prev_valid = true;
ERW_biogeochem.pH_bulk = nan(length(t), 1);
ERW_biogeochem.pH_bulk(1) = pH_bulk_raw;

% structure containing weathering variables
ERW_weathering_varnames = {'EW';'Wr';'psd';'d';'lamb';'SSA';
    'W_CaCO3';'W_MgCO3';'P_weathering_release';'Omega';'Omega_AlOH3';'delta_d';'SA'};

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function budget = ERW_Initial_Anion_Budget(ERW_TC_const_input)
I_An_source = 0;
if isfield(ERW_TC_const_input, 'I_An_source') && ~isempty(ERW_TC_const_input.I_An_source)
    I_An_source = ERW_TC_const_input.I_An_source(:);
elseif isfield(ERW_TC_const_input, 'rain_strong_anion_flux_umolc_m2_d')
    I_An_source = ERW_TC_const_input.rain_strong_anion_flux_umolc_m2_d;
end
if any(~isfinite(I_An_source(:)) | I_An_source(:) < 0)
    error('ERW:invalidBoundarySource', ...
        'Residual-anion boundary forcing must be finite and nonnegative.');
end
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
        E_H = 67.2; % [kJ/mol]: activation energies with respect to hydrogen
        E_w = 79; % [kJ/mol]: activation energies with respect to water
        E_OH = 0; % [kJ/mol]: activation energies with respect to hydroxide
        n_H = 0.47; % reaction order with respect to hydrogen
        n_OH = 1; % % reaction order with respect to hydroxide
        min_st = [0, 2, 0, 0, 0, 1, 0]; % Stochiometric coefficients [Ca, Mg, K, Na, Al, Si, P]
        K_sp = 10^(7.11); % Solubility product of the mineral [-] https://booksite.elsevier.com/9780120885305/appendices/Web_Appendices.pdf
    elseif strcmp(mineral, 'Fe_forsterite') % FeMgSiO4
        MM_min = 172/conv_mol; % [g/mol-conv]
        k_diss_H = 10^(-6.85)*24*3600*conv_mol; % (-5.37) [mol-conv m-2 d-1]
        k_diss_w = 10^(-10.64)*24*3600*conv_mol; % [mol-conv m-2 d-1]
        k_diss_OH = 0;
        E_H = 67; % [kJ/mol]: activation energies
        E_w = 79; % [kJ/mol]
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
        E_H = 54.7; % [kJ/mol]: activation energies
        E_w = 54.7; % [kJ/mol]
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
        E_H = 65; % [kJ/mol]: activation energies
        E_w = 70; % [kJ/mol]
        E_OH = 71;
        n_H = 0.457; % reaction order
        n_OH = -0.572;
        min_st = [0, 0, 0, 1, 1, 3, 0]; % Stochiometric coefficients [Ca, Mg, K, Na, Al, Si, P]
        K_sp = 10^(2.7645/4); % NaAlSi3O8 + 4H+ = Na+ + Al+++ + 3H4SiO4, llnl.dat, per proton
    elseif strcmp(mineral, 'anorthite') % CaAl2Si2O8
        MM_min = 278/conv_mol; % [g/mol-conv]
        k_diss_H = 10^(-3.5)*24*3600*conv_mol; % [mol-conv m-2 d-1]
        k_diss_w = 10^(-9.12)*24*3600*conv_mol; % [mol-conv m-2 d-1]
        k_diss_OH = 0;
        E_H = 16.6; % [kJ/mol]: activation energies
        E_w = 17.8; % [kJ/mol]
        E_OH = 0;
        n_H = 1.4; % reaction order
        n_OH = 1;
        min_st = [1, 0, 0, 0, 2, 2, 0]; % Stochiometric coefficients [Ca, Mg, K, Na, Al, Si, P]
        K_sp = 10^(26.578/8); % CaAl2Si2O8 + 8H+ = Ca++ + 2Al+++ + 2H4SiO4, llnl.dat, per proton
    elseif strcmp(mineral, 'labradorite') % Na0.4Ca0.6Al1.6Si2.4O8 (http://webmineral.com/data/Labradorite.shtml)
        MM_min = 272/conv_mol; % [g/mol-conv]
        k_diss_H = 10^(-7.87)*24*3600*conv_mol; % [mol-conv m-2 d-1]
        k_diss_w = 10^(-10.91)*24*3600*conv_mol; % [mol-conv m-2 d-1]
        k_diss_OH = 0;
        E_H = 42; % [kJ/mol]: activation energies
        E_w = 45; % [kJ/mol]
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
        E_H = 78; % [kJ/mol]: activation energies
        E_w = 78; % [kJ/mol]
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
        E_H = 96; % [kJ/mol]: activation energies
        E_w = 40; % [kJ/mol]
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
        E_H = 250; % [kJ/mol]: activation energies
        E_w = 250;
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
        E_H = 88; % [kJ/mol]: activation energies
        E_w = 88;
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
        E_H = 71.1; % [kJ/mol]: activation energies
        E_w = 70.7;
        E_OH = 79.1;
        n_H = 0.338; % reaction order
        n_OH = -0.556;
        min_st = [2, 0, 0, 0, 2.53, 3, 0]; % Stochiometric coefficients [Ca, Mg, K, Na, Al, Si, P]]
        K_sp = 10^(32.23); % https://thermoddem.brgm.fr/species/epidote
    elseif strcmp(mineral, 'zoisite') % zoisite Ca2Al3Si3O12(OH)
        MM_min = 454.357/conv_mol; % [g/mol-conv]
        k_diss_H = 10^(-7.50)*24*3600*conv_mol; % [mol-conv m-2 d-1]
        k_diss_w = 10^(-11.2)*24*3600*conv_mol;
        k_diss_OH = 0;
        E_H = 66.1; % [kJ/mol]: activation energies
        E_w = 66.1;
        E_OH = 0;
        n_H = 0.5; % reaction order
        n_OH = 1;
        min_st = [2, 0, 0, 0, 3, 3, 0]; % Stochiometric coefficients [Ca, Mg, K, Na, Al, Si, P]]
        K_sp = 10^(43.85); % https://thermoddem.brgm.fr/species/epidote
    elseif strcmp(mineral, 'ferroactinolite') % Ca2(Mg0.75Fe0.25)5Si8.5O22(OH)2 - Amphibole  - Lewis et al. (2021)
        MM_min = 865.8/conv_mol; % [g/mol-conv]
        % Rate law of Hermanska et al. (2022, Chem. Geol. 597, 120807), whose
        % tremolite fit PHREEQC's phreeqc_rates.dat also lists for
        % ferroactinolite. The Palandri and Kharaka (2004) tremolite row used
        % before (10^-8.4, 18.9, 0.7; 10^-10.6, 94.4) is taken from Sverdrup
        % (1990), not fitted to tremolite experiments, and is about 60 times
        % the measured rate at pH 6 (AMPHIBOLE_KINETICS_2026-10-02.md).
        k_diss_H = 10^(-11.3)*24*3600*conv_mol; % [mol-conv m-2 d-1]
        k_diss_w = 10^(-13.1)*24*3600*conv_mol;
        k_diss_OH = 0;
        E_H = 50; % [kJ/mol]: activation energies
        E_w = 48;
        E_OH = 0;
        n_H = 0.2; % reaction order
        n_OH = 1;
        min_st = [2, 3.75, 0, 0, 0, 8.5, 0]; % Stochiometric coefficients [Ca, Mg, K, Na, Al, Si, P]]
        K_sp = 10^(7.3); % https://thermoddem.brgm.fr/species/actinolite
    elseif strcmp(mineral, 'quartz') % SiO2 - Lewis et al. (2021)
        MM_min = 60/conv_mol; % [g/mol-conv]
        k_diss_H = 10^(-11.36)*24*3600*conv_mol; % [mol-conv m-2 d-1]
        k_diss_w = 10^(-13.4)*24*3600*conv_mol; % [mol-conv m-2 d-1]
        k_diss_OH = 10^(-16.29)*24*3600*conv_mol;
        E_H = 90.9; % [kJ/mol]: activation energies
        E_w = 90.9; % [kJ/mol] % Yoshiki Kanzaki et al.(2022)
        E_OH = 90.9;
        n_H = 0.309; % reaction order
        n_OH = -0.411;
        min_st = [0, 0, 0, 0, 0, 1, 0]; % Stochiometric coefficients [Ca, Mg, K, Na, Al, Si, P]
        K_sp = 10^(-4); 
    elseif strcmp(mineral, 'actinolite') % Ca2(Mg2.25Fe2.5Al0.25)(Si7.75Al0.25)O22(OH)2 - https://thermoddem.brgm.fr/species/actinolite
        MM_min = 892/conv_mol; % [g/mol-conv]
        % Same amphibole rate law as ferroactinolite (Hermanska et al. 2022).
        k_diss_H = 10^(-11.3)*24*3600*conv_mol; % [mol-conv m-2 d-1]
        k_diss_w = 10^(-13.1)*24*3600*conv_mol;
        k_diss_OH = 0;
        E_H = 50; % [kJ/mol]: activation energies
        E_w = 48;
        E_OH = 0;
        n_H = 0.2; % reaction order
        n_OH = 1;
        min_st = [2, 2.25, 0, 0, 0.5, 7.75, 0]; % Stochiometric coefficients [Ca, Mg, K, Na, Al, Si, P]
        K_sp = 10^(7.13); 
    elseif strcmp(mineral, 'calcite') % CaCO3
        MM_min = 100.087/conv_mol; % [g/mol-conv]
        k_diss_H = 10^(-0.3)*24*3600*conv_mol; % [mol-conv m-2 d-1]
        k_diss_w = 10^(-5.81)*24*3600*conv_mol;
        k_diss_OH = 10^(-3.48)*24*3600*conv_mol;
        E_H = 14.4; % [kJ/mol]: activation energies
        E_w = 23.5;
        E_OH = 35.4;
        n_H = 1; % reaction order
        n_OH = 1;
        min_st = [1, 0, 0, 0, 0, 0, 0]; % Stochiometric coefficients [Ca, Mg, K, Na, Al, Si, P]
        K_sp = 10^(-8.35); 
    else
        error('No data for this mineral');
    end

end
