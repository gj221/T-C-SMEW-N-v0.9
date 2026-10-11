%#########################################################################
%   Function INVERT_BULK_pH                                             #
%#########################################################################
%%% Invert Porewater_to_Bulk_pH by bisection in log-H space to recover the
%%% pore-water H concentration that maps forward to a measured bulk extract
%%% pH.  This anchors the seed to the observation: the lab 1:n extract pH
%%% (the measured observable) sets the exchange/pH state, and the in-situ
%%% pore-water pH is the value that reproduces it through the forward extract
%%% operator (Kanzaki et al. 2024).
%%%
%%% Called once at initialisation from ERW_MOD_PARAM_excel when n_extract is
%%% set.  All downstream ERW code (carbonate system, Al(OH)3 partition, An
%%% closure) then runs at the returned pore-water H.
%%%
%%% By default, every trial pore-water H re-derives the exchange fractions from
%%% that H via Exchange_Fractions_GT (the same Gaines-Thomas relations the
%%% daily ERW solver enforces).  When initial_f_CEC and field_water_L_m2 are
%%% supplied, the initial dissolved and exchangeable inventories are combined
%%% and repartitioned at every trial H.  The residual anion An is
%%% recomputed to close charge balance.  The returned (H_pw, f_CEC, An, Al_w)
%%% are therefore a self-consistent fixed point of the daily exchange +
%%% charge balance, so the t=0 pre-equilibration is a confirming no-op rather
%%% than a state that must relax away from the observed pH.
%%%
%%% Exchange selectivity is GLOBAL: the forward operator uses the same
%%% empirical soil-class Gaines-Thomas constants as the daily ERW solver.
%%% An optional SCEPTER/Appelo (1994) H-referenced exchange scheme
%%% (bulk_pH_params.scheme = 'appelo' with .appelo.logKH/.alpha)
%%% is shared identically between the exchange solve here,
%%% the forward extract operator, and the daily ERW solver, so the returned
%%% state is a fixed point of the SAME operator that produces pH_bulk_obs.
%%% The default (no .scheme) is the stock TC_GT operator with no fitted
%%% per-site K_Ca_H or alpha.
%%%
%%% INPUTS
%%%   pH_bulk_obs    observed bulk extract pH (e.g. ISO 10390 measurement)
%%%   Ca_pw .. NH4_pw   initial pore-water concentrations [umol/L]
%%%   CO2_air_ref    atmospheric CO2 [umol/L, T&C units] for the extract boundary.
%%%                  Pass CO2_atm (not soil-air CO2): lab pH measurements equilibrate
%%%                  at atmospheric pCO2 regardless of field soil-CO2 conditions.
%%%   CO2_w_ref      soil-dissolved CO2 for pore-water An balance (T&C units)
%%%   k1 k2 kw knh3 kH  carbonate / water / ammonium / Henry constants at reference T
%%%   K1_Al..K4_Al   Al hydrolysis constants [umol/L]
%%%   CEC_tot_ref    total CEC [umolc/m2]
%%%   rho_bulk       soil bulk density [g/m3]
%%%   Zrw            reactive zone depth [m]
%%%   conv_mol       umol-to-mol factor (1e6)
%%%   bulk_pH_params extraction parameter struct (.n_extract,
%%%                  .K_CEC, .field_water_L_m2, ...)
%%%   K_AlOH3        log *Kso for Al(OH)3 dissolution
%%%   Al_CEC_scale   ERW Al storage scaling factor (= conv_Al)
%%%
%%% OUTPUTS
%%%   H_pw           pore-water [H+] that maps forward to pH_bulk_obs [umol/L]
%%%   bulk_pH_params input struct with K_CEC(5) unchanged
%%%   f_Ca_eq..f_H_eq  exchange fractions consistent with H_pw (Gaines-Thomas)
%%%   An_eq          residual strong-acid anion closing charge balance [umolc/L]
%%%   Al_w_eq        total dissolved Al at H_pw [umol/L, ERW model units]

function [H_pw, bulk_pH_params, f_Ca_eq, f_Mg_eq, f_K_eq, f_Na_eq, f_Al_eq, f_H_eq, ...
    An_eq, Al_w_eq, Ca_eq, Mg_eq, K_eq, Na_eq] = Invert_Bulk_pH( ...
    pH_bulk_obs, ...
    Ca_pw, Mg_pw, K_pw, Na_pw, NO3_pw, PO4_pw, NH4_pw, ...
    CO2_air_ref, CO2_w_ref, ...
    k1_ref, k2_ref, k_w_ref, k_nh3_ref, k_H_ref, ...
    K1_Al, K2_Al, K3_Al, K4_Al, ...
    CEC_tot_ref, rho_bulk, Zrw, conv_mol, ...
    bulk_pH_params, K_AlOH3, Al_CEC_scale)

%%% Forward residual is a pure function of trial H and the fixed input
%%% parameters.  Warm-start state must not change the residual during a scan.
warning_state = warning('query', 'ERW:bulkPHLargeOffset');
warning_cleanup = onCleanup(@() warning(warning_state.state, 'ERW:bulkPHLargeOffset'));
warning('off', 'ERW:bulkPHLargeOffset'); % large offsets expected at scan trial values
inverse_bulk_pH_params = bulk_pH_params;
inverse_bulk_pH_params.H_extract_anchor_umolL = ...
    10^(-pH_bulk_obs)*conv_mol;
trial_fn = @(Ht, bpp) Initial_Bulk_pH_Trial( ...
    Ht, Ca_pw, Mg_pw, K_pw, Na_pw, NO3_pw, PO4_pw, NH4_pw, ...
    CO2_air_ref, CO2_w_ref, ...
    k1_ref, k2_ref, k_w_ref, k_nh3_ref, k_H_ref, K1_Al, K2_Al, K3_Al, K4_Al, ...
    CEC_tot_ref, rho_bulk, Zrw, conv_mol, bpp, K_AlOH3, Al_CEC_scale);


%%% Scan log-H grid (pH 0 to 14) for sign changes of (bulk_pH - pH_bulk_obs).
%%% Try a narrow window centred on pH_bulk_obs first (cheap): the model's own
%%% logic already treats a >2.5-unit bulk-vs-pore offset as suspicious
%%% (ERW:bulkPHLargeOffset), so a +-3 window has margin to spare and is a
%%% real prior here, not a possibly-stale carried-over guess (pH_bulk_obs is
%%% the actual measurement, not a warm-started previous state).  Only pay for
%%% the full [0,14] scan when nothing is found nearby.  The narrow window
%%% MUST match the full scan's point density -- a coarser scan can resolve
%%% the spurious Al(OH)3 buffer transition differently and select the wrong
%%% (rejected) branch (see Porewater_to_Bulk_pH.m for the same lesson).
full_logH_lo = log(10^(-14)*conv_mol);
full_logH_hi = log(conv_mol);
full_n_points = 161;
log_step = (full_logH_hi - full_logH_lo) / (full_n_points - 1);
narrow_pH_lo = max(0, pH_bulk_obs - 3);
narrow_pH_hi = min(14, pH_bulk_obs + 3);
narrow_logH_hi = log(10^(-narrow_pH_lo) * conv_mol);
narrow_logH_lo = log(10^(-narrow_pH_hi) * conv_mol);
narrow_n_points = max(2, round((narrow_logH_hi - narrow_logH_lo) / log_step) + 1);
[logH_grid, f_grid] = invert_bulk_pH_residual_scan( ...
    trial_fn, inverse_bulk_pH_params, pH_bulk_obs, ...
    narrow_logH_lo, narrow_logH_hi, narrow_n_points);
idx_root = invert_bulk_pH_sign_changes(f_grid);
if isempty(idx_root)
    [logH_grid, f_grid] = invert_bulk_pH_residual_scan( ...
        trial_fn, inverse_bulk_pH_params, pH_bulk_obs, ...
        full_logH_lo, full_logH_hi, full_n_points);
    idx_root = invert_bulk_pH_sign_changes(f_grid);
end

if isempty(idx_root)
    error('ERW:invertBulkPH', ...
        'Could not invert observed bulk pH %.4g to a pore-water pH.', pH_bulk_obs);
end

%%% Prefer the root continuous with the supplied pore-water anchor.  If this
%%% is the first solve, use the observed extract pH only as the neutral root
%%% selection anchor.  No directional porewater-vs-extract offset is imposed.
mid_logH  = 0.5 * (logH_grid(idx_root) + logH_grid(idx_root+1));
mid_pH_pw = -log10(exp(mid_logH) / conv_mol);
target_pH = pH_bulk_obs;
if isfield(bulk_pH_params, 'H_pore_anchor_umolL') && ...
        isfinite(bulk_pH_params.H_pore_anchor_umolL) && bulk_pH_params.H_pore_anchor_umolL > 0
    target_pH = -log10(bulk_pH_params.H_pore_anchor_umolL / conv_mol);
end
[~, i_best] = min(abs(mid_pH_pw - target_pH));
H_lo = logH_grid(idx_root(i_best));
H_hi = logH_grid(idx_root(i_best)+1);
f_lo = f_grid(idx_root(i_best));
f_hi = f_grid(idx_root(i_best)+1);

%%% Bisect the same fixed-parameter residual in log-H space.
if abs(f_lo) <= abs(f_hi)
    H_best = H_lo;
    f_best = f_lo;
else
    H_best = H_hi;
    f_best = f_hi;
end
for i_inv = 1:80
    H_mid = 0.5*(H_lo + H_hi);
    if H_mid == H_lo || H_mid == H_hi
        break
    end
    pH_bulk_mid = trial_fn(exp(H_mid), inverse_bulk_pH_params);
    f_mid = pH_bulk_mid - pH_bulk_obs;
    if abs(f_mid) < abs(f_best)
        H_best = H_mid;
        f_best = f_mid;
    end
    if abs(f_mid) < 1e-10
        break;
    end
    if sign(f_mid) == sign(f_lo)
        H_lo = H_mid;
        f_lo = f_mid;
    else
        H_hi = H_mid;
        f_hi = f_mid;
    end
end
H_pw = exp(H_best);

%%% Re-evaluate the converged trial to return the self-consistent state.
[pH_bulk_final, trial_state] = trial_fn(H_pw, inverse_bulk_pH_params);
if abs(pH_bulk_final-pH_bulk_obs) > 1e-8
    error('ERW:invertBulkPHNoConvergence', ...
        'Bulk-pH inversion residual %.6g exceeds 1e-8.', pH_bulk_final-pH_bulk_obs);
end
bulk_pH_params.DIC_field_umol_L = trial_state.DIC_field_umol_L;
f_Ca_eq = trial_state.f_Ca;  f_Mg_eq = trial_state.f_Mg;  f_K_eq = trial_state.f_K;
f_Na_eq = trial_state.f_Na;  f_Al_eq = trial_state.f_Al;  f_H_eq = trial_state.f_H;
An_eq   = trial_state.An;    Al_w_eq = trial_state.Al_w;
Ca_eq = trial_state.Ca; Mg_eq = trial_state.Mg;
K_eq = trial_state.K; Na_eq = trial_state.Na;

end % Invert_Bulk_pH

%% =======================================================================
function [logH_grid, f_grid] = invert_bulk_pH_residual_scan( ...
    trial_fn, bulk_pH_params, pH_bulk_obs, logH_lo, logH_hi, n_points)
% Evaluate one fixed-parameter residual over a log-H window.
logH_grid = linspace(logH_lo, logH_hi, n_points);
f_grid = nan(size(logH_grid));
for i_inv = 1:numel(logH_grid)
    try
        pH_bulk = trial_fn(exp(logH_grid(i_inv)), bulk_pH_params);
        f_grid(i_inv) = pH_bulk - pH_bulk_obs;
    catch ME
        if startsWith(ME.identifier, 'ERW:')
            f_grid(i_inv) = NaN;
        else
            rethrow(ME)
        end
    end
end
end

%% =======================================================================
function idx_root = invert_bulk_pH_sign_changes(f_grid)
valid_grid = isfinite(f_grid);
idx_root = find(valid_grid(1:end-1) & valid_grid(2:end) & ...
    (f_grid(1:end-1) == 0 | f_grid(2:end) == 0 | ...
     sign(f_grid(1:end-1)) ~= sign(f_grid(2:end))));
end

%% =======================================================================
function [pH_bulk, trial_state] = Initial_Bulk_pH_Trial( ...
    Ht, Ca_trial, Mg_trial, K_trial, Na_trial, ...
    NO3_trial, PO4_trial, NH4_trial, CO2_air_trial, CO2_w_trial, ...
    k1_trial, k2_trial, k_w_trial, k_nh3_trial, k_H_trial, K1_Al, K2_Al, K3_Al, K4_Al, ...
    CEC_tot_trial, rho_bulk_trial, Zrw_trial, conv_mol_trial, bulk_pH_params_trial, ...
    K_AlOH3_trial, Al_CEC_scale_trial)
%%% Evaluate the forward bulk pH at trial pore-water H = Ht, with exchange
%%% fractions and the residual anion re-derived self-consistently at Ht.
%%% CO2_air_trial  : lab-air CO2 for the extract boundary
%%% CO2_w_trial    : soil-dissolved CO2 (for pore-water alkalinity / An)

%%% Cation-exchange scheme selection, shared with the axial daily solver and
%%% the extract operator so the inversion returns a fixed point of the SAME
%%% operator that produces pH_bulk_obs (see Porewater_to_Bulk_pH.m lines 100+).
exchange_scheme = 'TC_GT';
if isfield(bulk_pH_params_trial, 'scheme') && ~isempty(bulk_pH_params_trial.scheme)
    exchange_scheme = bulk_pH_params_trial.scheme;
    if isstring(exchange_scheme) && isscalar(exchange_scheme)
        exchange_scheme = char(exchange_scheme);
    end
    if ~ischar(exchange_scheme) || ~any(strcmpi(exchange_scheme, {'TC_GT','appelo'}))
        error('ERW:invalidParam', 'Exchange scheme must be TC_GT or appelo.');
    end
end
scheme_is_appelo = strcmpi(exchange_scheme, 'appelo');
appelo_params = [];
if isfield(bulk_pH_params_trial, 'appelo')
    appelo_params = bulk_pH_params_trial.appelo;
end

%%% Effective Ca-H selectivity with the shared Gaines-Thomas constants.
K_CEC_eff = bulk_pH_params_trial.K_CEC;
if scheme_is_appelo
    [K_Na_H_t, K_K_H_t, K_Ca_H2_t, K_Mg_H2_t, K_Al_H3_t, appelo_alpha_t] = ...
        Resolve_Appelo_Constants(appelo_params);
else
    K_Na_H_t = 0; K_K_H_t = 0; K_Ca_H2_t = 0; K_Mg_H2_t = 0; K_Al_H3_t = 0;
    appelo_alpha_t = 0;
end

Al_pin = isfield(bulk_pH_params_trial, 'Al_solubility_pin') && ...
    isequal(bulk_pH_params_trial.Al_solubility_pin, true);
initial_f_CEC = [];
has_conserved_exchange = false;
if isfield(bulk_pH_params_trial, 'initial_f_CEC')
    initial_f_CEC = bulk_pH_params_trial.initial_f_CEC(:).';
    if numel(initial_f_CEC) ~= 6 || any(~isfinite(initial_f_CEC)) || ...
            any(initial_f_CEC < 0) || initial_f_CEC(1) <= 0 || ...
            abs(sum(initial_f_CEC)-1) > 1e-10
        error('ERW:invertBulkPH', ...
            'initial_f_CEC must contain six nonnegative fractions that sum to one with positive Ca.');
    end
    initial_f_CEC = initial_f_CEC ./ sum(initial_f_CEC);
    has_conserved_exchange = true;
end
if has_conserved_exchange
    if ~isfield(bulk_pH_params_trial, 'field_water_L_m2') || ...
            ~(isfinite(bulk_pH_params_trial.field_water_L_m2) && ...
            bulk_pH_params_trial.field_water_L_m2 > 0)
        error('ERW:invertBulkPH', ...
            'Conserved exchange inventories require a positive field_water_L_m2.');
    end
    if any(~isfinite([Ca_trial, Mg_trial, K_trial, Na_trial])) || ...
            any([Ca_trial, Mg_trial, K_trial, Na_trial] < 0)
        error('ERW:invertBulkPH', ...
            'Conserved exchange inventories require finite nonnegative dissolved cations.');
    end
    Al_w_initial = 0;
    if isfield(bulk_pH_params_trial, 'initial_Al_w') && ...
            isfinite(bulk_pH_params_trial.initial_Al_w)
        Al_w_initial = bulk_pH_params_trial.initial_Al_w;
    elseif initial_f_CEC(5) > 0
        error('ERW:invertBulkPH', ...
            'A positive initial Al exchange fraction requires finite initial_Al_w.');
    end
    cec_per_L = CEC_tot_trial / bulk_pH_params_trial.field_water_L_m2;
end
if has_conserved_exchange
    [f_Ca_t, f_Mg_t, f_K_t, f_Na_t, f_Al_t, f_H_t, ...
        Ca_trial, Mg_trial, K_trial, Na_trial, Al_w_trial, Al_trial] = ...
        equilibrate_exchange_totals(Ht, Ca_trial, Mg_trial, K_trial, Na_trial, ...
        Al_w_initial, initial_f_CEC, cec_per_L, K_CEC_eff, ...
        K1_Al, K2_Al, K3_Al, K4_Al, Al_CEC_scale_trial, ...
        scheme_is_appelo, appelo_alpha_t, K_Na_H_t, K_K_H_t, K_Ca_H2_t, K_Mg_H2_t, K_Al_H3_t, ...
        Al_pin);
else
    if any(~isfinite([Ca_trial, Mg_trial, K_trial, Na_trial])) || ...
            any([Ca_trial, Mg_trial, K_trial, Na_trial] < 0) || Ca_trial <= 0
        error('ERW:invertBulkPH', ...
            'Unconstrained inversion requires finite nonnegative cations and positive Ca.');
    end
    Al_w_trial = Initialize_ERW_Equilibrium(Ht, K1_Al, K2_Al, K3_Al, K4_Al, ...
        conv_mol_trial, K_AlOH3_trial, Al_CEC_scale_trial);
    if ~isfinite(Al_w_trial) || Al_w_trial <= 0
        error('ERW:invertBulkPH', ...
            'Cannot initialize trial Al_w from Al(OH)3 equilibrium at H = %.3g umol/L.', Ht);
    end
    Al_den = Ht^4 + Ht^3*K1_Al + Ht^2*K1_Al*K2_Al + ...
        Ht*K1_Al*K2_Al*K3_Al + K1_Al*K2_Al*K3_Al*K4_Al;
    Al_trial = (Ht^4 / Al_den) * Al_w_trial;
    [f_Ca_t, f_Mg_t, f_K_t, f_Na_t, f_Al_t, f_H_t] = Exchange_Fractions_GT( ...
        Ht, Ca_trial, Mg_trial, K_trial, Na_trial, Al_w_trial, ...
        K_CEC_eff, K1_Al, K2_Al, K3_Al, K4_Al, Al_CEC_scale_trial, ...
        exchange_scheme, appelo_params);
end

%%% Use the same strong-anion charge basis as the retained daily state.
anion_convention = 'full_charge';
if isfield(bulk_pH_params_trial, 'initial_anion_convention')
    anion_convention = char(string(bulk_pH_params_trial.initial_anion_convention));
end
if strcmpi(anion_convention, 'carbonate_frame')
    HCO3_trial = k1_trial*CO2_w_trial/Ht;
    CO3_trial = k2_trial*k1_trial*CO2_w_trial/Ht^2;
    Alk_trial = HCO3_trial+2*CO3_trial-Ht+k_w_trial/Ht;
    An_trial = 2*Mg_trial+2*Ca_trial+Na_trial+K_trial+NH4_trial- ...
        Alk_trial-PO4_trial-NO3_trial;
elseif strcmpi(anion_convention, 'full_charge')
    An_trial = ERW_Residual_An_From_Charge( ...
        Ht, CO2_w_trial, Ca_trial, Mg_trial, Na_trial, K_trial, ...
        NH4_trial, NO3_trial, PO4_trial, Al_w_trial, ...
        k1_trial, k2_trial, k_w_trial, k_nh3_trial, K1_Al, K2_Al, K3_Al, K4_Al, ...
        Al_CEC_scale_trial);
else
    error('ERW:invalidParam', ...
        'initial_anion_convention must be full_charge or carbonate_frame.');
end
if ~(isfinite(An_trial) && An_trial >= 0)
    error('ERW:invertBulkPH', ...
        'Trial pH requires a negative residual strong-anion concentration.');
end

%%% Forward bulk pH with the trial field DIC conserved in the closed extract.
bulk_pH_params_forward = bulk_pH_params_trial;
DIC_field_trial = CO2_w_trial * ...
    (1+k1_trial/Ht+k1_trial*k2_trial/Ht^2);
bulk_pH_params_forward.DIC_field_umol_L = DIC_field_trial;
pH_bulk = Porewater_to_Bulk_pH( ...
    Ht, Ca_trial, Mg_trial, K_trial, Na_trial, ...
    Al_w_trial, Al_trial, An_trial, NO3_trial, PO4_trial, NH4_trial, ...
    CO2_air_trial, k1_trial, k2_trial, k_w_trial, k_nh3_trial, k_H_trial, ...
    K1_Al, K2_Al, K3_Al, K4_Al, ...
    f_Ca_t, f_Mg_t, f_K_t, f_Na_t, f_Al_t, f_H_t, ...
    CEC_tot_trial, rho_bulk_trial, Zrw_trial, conv_mol_trial, bulk_pH_params_forward, ...
    Al_CEC_scale_trial);

if nargout > 1
    trial_state = struct('f_Ca', f_Ca_t, 'f_Mg', f_Mg_t, 'f_K', f_K_t, ...
        'f_Na', f_Na_t, 'f_Al', f_Al_t, 'f_H', f_H_t, ...
        'An', An_trial, 'Al_w', Al_w_trial, 'Al', Al_trial, ...
        'Ca', Ca_trial, 'Mg', Mg_trial, 'K', K_trial, 'Na', Na_trial, ...
        'DIC_field_umol_L', DIC_field_trial);
end

end % Initial_Bulk_pH_Trial

function [f_Ca, f_Mg, f_K, f_Na, f_Al, f_H, ...
    Ca, Mg, K, Na, Al_w, Al] = equilibrate_exchange_totals( ...
    H, Ca_initial, Mg_initial, K_initial, Na_initial, Al_w_initial, ...
    f_initial, cec_per_L, K_CEC, K1_Al, K2_Al, K3_Al, K4_Al, Al_CEC_scale, ...
    scheme_is_appelo, appelo_alpha, K_Na_H, K_K_H, K_Ca_H2, K_Mg_H2, K_Al_H3, ...
    Al_solubility_pin)
if scheme_is_appelo
    [f_Ca, f_Mg, f_K, f_Na, f_Al, f_H, ...
        Ca, Mg, K, Na, Al_w, Al] = ...
        equilibrate_exchange_totals_appelo(H, Ca_initial, Mg_initial, K_initial, ...
        Na_initial, Al_w_initial, f_initial, cec_per_L, K1_Al, K2_Al, K3_Al, K4_Al, ...
        Al_CEC_scale, appelo_alpha, K_Na_H, K_K_H, K_Ca_H2, K_Mg_H2, K_Al_H3, ...
        Al_solubility_pin);
    return
end
if numel(K_CEC) ~= 5 || any(~isfinite(K_CEC)) || any(K_CEC <= 0)
    error('ERW:invertBulkPH', 'K_CEC must contain five positive finite constants.');
end
if ~(isfinite(cec_per_L) && cec_per_L > 0)
    error('ERW:invertBulkPH', 'CEC per field-water volume must be positive and finite.');
end
if ~(isfinite(Al_w_initial) && Al_w_initial >= 0) || ...
        ~(isfinite(Al_CEC_scale) && Al_CEC_scale > 0)
    error('ERW:invertBulkPH', 'Initial Al and Al_CEC_scale must be finite and nonnegative.');
end
K_Ca_Mg = K_CEC(1);
K_Ca_K = K_CEC(2);
K_Ca_Na = K_CEC(3);
K_Ca_Al = K_CEC(4);
K_Ca_H = K_CEC(5);
Ca_total = Ca_initial + cec_per_L*f_initial(1)/2;
Mg_total = Mg_initial + cec_per_L*f_initial(2)/2;
K_total = K_initial + cec_per_L*f_initial(3);
Na_total = Na_initial + cec_per_L*f_initial(4);
Al_total = Al_w_initial + cec_per_L*f_initial(5)*Al_CEC_scale/3;
Al_den = H^4 + H^3*K1_Al + H^2*K1_Al*K2_Al + ...
    H*K1_Al*K2_Al*K3_Al + K1_Al*K2_Al*K3_Al*K4_Al;
if ~(isfinite(Al_den) && Al_den > 0 && isfinite(Ca_total) && Ca_total > 0)
    error('ERW:invertBulkPH', 'Invalid conserved inventory or Al hydrolysis denominator.');
end
Al_fraction = H^4/Al_den;
f_Ca_hi = min(1 - 1e-12, 2*Ca_total/cec_per_L - 1e-12);
if f_Ca_hi <= 1e-12
    error('ERW:invertBulkPH', 'The conserved Ca inventory cannot occupy a positive exchange fraction.');
end
site_residual = @(x) exchange_site_residual(x, H, Ca_total, Mg_total, ...
    K_total, Na_total, Al_total, cec_per_L, K_Ca_Mg, K_Ca_K, ...
    K_Ca_Na, K_Ca_Al, K_Ca_H, Al_fraction, Al_CEC_scale);
f_lo = 1e-12;
g_lo = site_residual(f_lo);
g_hi = site_residual(f_Ca_hi);
if ~(isfinite(g_lo) && isfinite(g_hi) && sign(g_lo) ~= sign(g_hi))
    error('ERW:invertBulkPH', 'No physical Gaines-Thomas partition satisfies the conserved inventories.');
end
for i = 1:100
    f_mid = 0.5*(f_lo + f_Ca_hi);
    g_mid = site_residual(f_mid);
    if abs(g_mid) < 1e-12
        break;
    end
    if sign(g_mid) == sign(g_lo)
        f_lo = f_mid;
        g_lo = g_mid;
    else
        f_Ca_hi = f_mid;
    end
end
f_Ca = f_mid;
[~, f_Mg, f_K, f_Na, f_Al, f_H] = exchange_site_residual( ...
    f_Ca, H, Ca_total, Mg_total, K_total, Na_total, Al_total, ...
    cec_per_L, K_Ca_Mg, K_Ca_K, K_Ca_Na, K_Ca_Al, K_Ca_H, ...
    Al_fraction, Al_CEC_scale);
Ca = Ca_total - cec_per_L*f_Ca/2;
Mg = Mg_total - cec_per_L*f_Mg/2;
K = K_total - cec_per_L*f_K;
Na = Na_total - cec_per_L*f_Na;
Al_w = Al_total - cec_per_L*f_Al*Al_CEC_scale/3;
Al = Al_fraction*Al_w;
state = [f_Ca, f_Mg, f_K, f_Na, f_Al, f_H, Ca, Mg, K, Na, Al_w, Al];
if any(~isfinite(state)) || any(state < -1e-9) || abs(sum(state(1:6)) - 1) > 1e-8
    error('ERW:invertBulkPH', 'Conserved exchange equilibration produced an unphysical state.');
end
if any([Ca,Mg,K,Na,Al_w,Al] < 0)
    error('ERW:invertBulkPH', 'Conserved exchange equilibration produced a negative state.');
end
end

function [f_Ca, f_Mg, f_K, f_Na, f_Al, f_H, ...
    Ca, Mg, K, Na, Al_w, Al] = equilibrate_exchange_totals_appelo( ...
    H, Ca_initial, Mg_initial, K_initial, Na_initial, Al_w_initial, ...
    f_initial, cec_per_L, K1_Al, K2_Al, K3_Al, K4_Al, Al_CEC_scale, ...
    appelo_alpha, K_Na_H, K_K_H, K_Ca_H2, K_Mg_H2, K_Al_H3, ...
    Al_solubility_pin)
if nargin < 20
    Al_solubility_pin = false;
end
if ~(isfinite(cec_per_L) && cec_per_L > 0)
    error('ERW:invertBulkPH', 'CEC per field-water volume must be positive and finite.');
end
if any(~isfinite([K_Na_H K_K_H K_Ca_H2 K_Mg_H2 K_Al_H3])) || ...
        any([K_Na_H K_K_H K_Ca_H2 K_Mg_H2 K_Al_H3] <= 0)
    error('ERW:invertBulkPH', 'Appelo selectivity constants must be positive and finite.');
end
if ~(isfinite(appelo_alpha) && appelo_alpha >= 0)
    error('ERW:invertBulkPH', 'Appelo alpha must be a finite nonnegative scalar.');
end
if Al_solubility_pin
    % SCEPTER semantics: aqueous Al is controlled by amorphous Al(OH)3
    % solubility (Al(OH)3 + 3H+ = Al3+ + 3H2O, logK = 10.8) and
    % exchangeable Al derives from solution Al3+; the measured
    % exchange-Al heritage is not conserved.  Units follow the daily
    % solver (ERW_function al3_pin_fn): Al3 [umol-equiv/L] =
    % 10^(K_AlOH3-12)*H^3 = 10^-1.2*H^3 with H in umol/L; the store-unit
    % convention (Al_CEC_scale) enters only when the daily state is
    % packed (p(6) = Al3*Al_CEC_scale), never in the field-water
    % inventory basis used here.
    Al_total = 10^-1.2 * H^3;
else
    Al0_u = Al_w_initial / Al_CEC_scale;
    Al_total = Al0_u + cec_per_L*f_initial(5)/3;
end
Ca_total = Ca_initial + cec_per_L*f_initial(1)/2;
Mg_total = Mg_initial + cec_per_L*f_initial(2)/2;
K_total  = K_initial  + cec_per_L*f_initial(3);
% Sodium is conserved like every other base cation.  The pin is a
% statement about aluminium only (aqueous Al fixed by Al(OH)3
% solubility); tying the Na inventory to it was a workaround for the old
% seed's f_Na ~ 0.8, and with f_Na now derived from the Appelo law the
% exchangeable slot is meaningful and must be carried.  Dropping it
% leaves Na_total too small for the exchanger, the solver drives
% dissolved Na to zero, and the log-space Na exchange residual then has
% no root at all.
Na_total = Na_initial + cec_per_L*f_initial(4);
Al3_frac = H^4 / (H^4 + H^3*K1_Al + H^2*K1_Al*K2_Al + H*K1_Al*K2_Al*K3_Al + ...
    K1_Al*K2_Al*K3_Al*K4_Al);
if ~(isfinite(Al3_frac) && Al3_frac > 0 && isfinite(Ca_total) && Ca_total > 0)
    error('ERW:invertBulkPH', 'Invalid conserved inventory or Al hydrolysis denominator.');
end
site_residual = @(x) appelo_site_residual(x, H, Ca_total, Mg_total, ...
    K_total, Na_total, Al_total, cec_per_L, appelo_alpha, ...
    K_Na_H, K_K_H, K_Ca_H2, K_Mg_H2, K_Al_H3, Al3_frac, ...
    Al_solubility_pin);
f_lo = 1e-24;
f_hi = 1 - 1e-14;
g_lo = site_residual(f_lo);
g_hi = site_residual(f_hi);
if ~(isfinite(g_lo) && isfinite(g_hi) && sign(g_lo) ~= sign(g_hi))
    error('ERW:invertBulkPH', 'No physical Appelo partition satisfies the conserved inventories.');
end
for i = 1:300
    f_mid = min(f_hi, max(f_lo, 0.5*(f_lo + f_hi)));
    g_mid = site_residual(f_mid);
    if abs(g_mid) < 1e-12
        break;
    end
    if sign(g_mid) == sign(g_lo)
        f_lo = f_mid;
        g_lo = g_mid;
    else
        f_hi = f_mid;
    end
end
f_H = 0.5*(f_lo + f_hi);
[~, f_Na, f_K, f_Ca, f_Mg, f_Al] = site_residual(f_H);
Ca = Ca_total - cec_per_L*f_Ca/2;
Mg = Mg_total - cec_per_L*f_Mg/2;
K = K_total - cec_per_L*f_K;
Na = Na_total - cec_per_L*f_Na;
if Al_solubility_pin
    Al_u = Al_total / Al3_frac;
else
    Al_u = Al_total - cec_per_L*f_Al/3;
end
Al_w = Al_u * Al_CEC_scale;
Al = Al3_frac * Al_w;
state = [f_Ca, f_Mg, f_K, f_Na, f_Al, f_H, Ca, Mg, K, Na, Al_w, Al];
if any(~isfinite(state)) || any(state < -1e-9) || abs(sum(state(1:6)) - 1) > 1e-8
    error('ERW:invertBulkPH', 'Conserved Appelo equilibration produced an unphysical state.');
end
if any([Ca,Mg,K,Na,Al_w,Al] < 0)
    error('ERW:invertBulkPH', 'Conserved Appelo equilibration produced a negative state.');
end
end

function [residual, f_Na, f_K, f_Ca, f_Mg, f_Al] = appelo_site_residual( ...
    f_H, H, Ca_total, Mg_total, K_total, Na_total, Al_total, cec_per_L, ...
    appelo_alpha, K_Na_H, K_K_H, K_Ca_H2, K_Mg_H2, K_Al_H3, Al3_frac, ...
    Al_solubility_pin)
if ~(isfinite(f_H) && f_H > 0)
    residual = NaN;
    f_Na = NaN; f_K = NaN; f_Ca = NaN; f_Mg = NaN; f_Al = NaN;
    return
end
gamma = 10^(appelo_alpha*f_H);
rho1 = (f_H/H) * gamma;
rho2 = rho1^2;
rho3 = rho1^3;
A_Na = K_Na_H * rho1;
A_K  = K_K_H  * rho1;
A_Ca = K_Ca_H2* rho2;
A_Mg = K_Mg_H2* rho2;
f_Na = A_Na * Na_total / (1 + A_Na*cec_per_L);
f_K  = A_K  * K_total  / (1 + A_K *cec_per_L);
f_Ca = A_Ca * Ca_total / (1 + A_Ca*cec_per_L/2);
f_Mg = A_Mg * Mg_total / (1 + A_Mg*cec_per_L/2);
if Al_solubility_pin
    f_Al = K_Al_H3 * Al_total * rho3;
else
    A_Al = K_Al_H3 * Al3_frac * rho3;
    f_Al = A_Al * Al_total / (1 + A_Al*cec_per_L/3);
end
residual = f_H + f_Na + f_K + f_Ca + f_Mg + f_Al - 1;
end

function [residual, f_Mg, f_K, f_Na, f_Al, f_H] = exchange_site_residual( ...
    f_Ca, H, Ca_total, Mg_total, K_total, Na_total, Al_total, ...
    cec_per_L, K_Ca_Mg, K_Ca_K, K_Ca_Na, K_Ca_Al, K_Ca_H, ...
    Al_fraction, Al_CEC_scale)
Ca = Ca_total - cec_per_L*f_Ca/2;
if ~(isfinite(Ca) && Ca > 0)
    residual = NaN;
    f_Mg = NaN; f_K = NaN; f_Na = NaN; f_Al = NaN; f_H = NaN;
    return
end
A_Mg = f_Ca/(K_Ca_Mg*Ca);
A_K = sqrt(f_Ca/(K_Ca_K*Ca));
A_Na = sqrt(f_Ca/(K_Ca_Na*Ca));
A_Al = (Al_fraction/Al_CEC_scale)*sqrt(f_Ca^3/(K_Ca_Al*Ca^3));
f_Mg = A_Mg*Mg_total/(1 + A_Mg*cec_per_L/2);
f_K = A_K*K_total/(1 + A_K*cec_per_L);
f_Na = A_Na*Na_total/(1 + A_Na*cec_per_L);
f_Al = A_Al*Al_total/(1 + A_Al*cec_per_L*Al_CEC_scale/3);
f_H = H*sqrt(f_Ca/(K_Ca_H*Ca));
residual = f_Ca + f_Mg + f_K + f_Na + f_Al + f_H - 1;
end
