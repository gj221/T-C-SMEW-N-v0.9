%#########################################################################
%   Subfunction POREWATER_TO_BULK_pH                                    #
%#########################################################################
%%% References  Kanzaki et al. 2024;  ISO 10390
function [pH_bulk, diag] = Porewater_to_Bulk_pH(...
    H_pw, Ca_pw, Mg_pw, Kion_pw, Na_pw, Al_pw, Al3_pw, ...
    An_pw, NO3_pw, PO4_pw, NH4_pw, CO2_air_extract, ...
    k1, k2, k_w, k_nh3, k_H, K1_Al, K2_Al, K3_Al, K4_Al, ...
    f_Ca_0, f_Mg_0, f_K_0, f_Na_0, f_Al_0, f_H_0, ...
    CEC_tot, rho_bulk, Zrw, conv_mol, bulk_pH_params, Al_CEC_scale)
%%%%%% INPUT (concentrations in umol/L unless noted)
%%% H_pw, Ca_pw, Mg_pw, Kion_pw, Na_pw   pore-water cations
%%% Al_pw, Al3_pw                        total dissolved Al, free Al3+
%%% An_pw                                residual bulk anion pool as aqueous concentration [umolc/L]
%%% NO3_pw, PO4_pw, NH4_pw               ERW explicit N/P ions
%%% CO2_air_extract                      atmospheric CO2 [umol/L, T&C units] for the extract
%%%                                       boundary.  Lab pH measurements (1:1, 1:2.5, 1:5
%%%                                       water:soil) equilibrate at atmospheric pCO2, not
%%%                                       soil-air CO2.  Always pass CO2_atm here.
%%% k1, k2, k_w, k_H                     carbonate/water constants
%%% K1_Al..K4_Al                         Al hydrolysis constants [umol/L]
%%% f_Ca_0..f_H_0                        exchange fractions from pore-water solver
%%% CEC_tot                              total CEC [umolc per m2]
%%% rho_bulk                             soil bulk density [g/m3]
%%% Zrw                                  reactive zone depth [m]
%%% conv_mol                             umol-to-mol factor (1e6)
%%% Al_CEC_scale                         ERW solver Al storage scaling (= conv_Al)
%%% bulk_pH_params                       struct with:
%%%   .n_extract        water:soil extraction ratio [required; set via pH_DilutionFactor in Excel]
%%%   .K_CEC            [K_Ca_Mg K_Ca_K K_Ca_Na K_Ca_Al K_Ca_H]
%%%                     Gaines-Thomas constants; units are [-], [umol/L],
%%%                     [umol/L], [L/umol], [umol/L], respectively
%%%   .field_water_L_m2 field dissolved-inventory water volume [L/m2]
%%%   .H_extract_anchor_umolL optional extract-root anchor [umol/L]
%%%   .H_init           optional previous-step H seed [umol/L]
%%%   .pH_init          optional previous-step pH seed [-]
%%%   .scheme           'TC_GT' (default) or H-referenced 'appelo'
%%%   .appelo           struct used when scheme = 'appelo':
%%%                     .logKH  [log10 KH-Na, KH-K, KH-Ca, KH-Mg, KH-Al],
%%%                             SCEPTER defaults [5.9 4.8 10.47 10.786 16.47]
%%%                     .alpha  Appelo exponent, SCEPTER default 3.4
%%%%%% OUTPUT
%%% pH_bulk    bulk soil pH (1:n extract analog) [-]
%%%%%%

%%% Validate the prescribed extraction and exchange parameters.
if ~isfinite(Al_CEC_scale) || Al_CEC_scale <= 0
    error('ERW:invalidParam', 'Al_CEC_scale must be positive.');
end
if ~isfield(bulk_pH_params,'H_prev_valid'); bulk_pH_params.H_prev_valid = false; end
if ~isfield(bulk_pH_params,'n_extract') || ~isfinite(bulk_pH_params.n_extract) || bulk_pH_params.n_extract <= 0
    error('ERW:missingParam', 'bulk_pH_params.n_extract is required; set pH_DilutionFactor in Excel MOD_PARAM for this site.');
end
if ~isfield(bulk_pH_params,'field_water_L_m2') || ...
        ~(isfinite(bulk_pH_params.field_water_L_m2) && bulk_pH_params.field_water_L_m2 > 0)
    error('ERW:missingParam', ...
        'bulk_pH_params.field_water_L_m2 must be positive to conserve the field dissolved inventory.');
end
if ~(isfinite(rho_bulk) && rho_bulk > 0 && isfinite(Zrw) && Zrw > 0)
    error('ERW:invalidParam', 'rho_bulk and Zrw must be positive.');
end

%%% CO2 input check
if any(~isfinite(CO2_air_extract(:))) || any(CO2_air_extract(:) <= 0)
    error('ERW:invalidParam', 'CO2_air_extract must be finite and positive.');
end
if any(CO2_air_extract(:) > 5e4)
    error('ERW:unitCheck', ...
        'CO2_air_extract=%.3g exceeds the physical gas-concentration range.', CO2_air_extract);
end

if ~isfinite(H_pw) || H_pw <= 0
    error('ERW:invalidParam', 'Pore-water H must be positive.');
end
pH_pw = -log10(H_pw / conv_mol);

K_Ca_Mg = bulk_pH_params.K_CEC(1); K_Ca_K  = bulk_pH_params.K_CEC(2); K_Ca_Na = bulk_pH_params.K_CEC(3);
K_Ca_Al = bulk_pH_params.K_CEC(4); K_Ca_H  = bulk_pH_params.K_CEC(5);

%%% SCEPTER uses gamma = 10^(alpha*f_H); its optional beta term is off.
exchange_scheme = 'TC_GT';
if isfield(bulk_pH_params,'scheme') && ~isempty(bulk_pH_params.scheme)
    exchange_scheme = bulk_pH_params.scheme;
    if isstring(exchange_scheme) && isscalar(exchange_scheme)
        exchange_scheme = char(exchange_scheme);
    end
    if ~ischar(exchange_scheme) || ~any(strcmpi(exchange_scheme, {'TC_GT','appelo'}))
        error('ERW:invalidParam', 'Exchange scheme must be TC_GT or appelo.');
    end
end
scheme_is_appelo = strcmpi(exchange_scheme,'appelo');
if scheme_is_appelo
    appelo_params = [];
    if isfield(bulk_pH_params, 'appelo')
        appelo_params = bulk_pH_params.appelo;
    end
    [K_Na_H, K_K_H, K_Ca_H2, K_Mg_H2, K_Al_H3, appelo_alpha] = ...
        Resolve_Appelo_Constants(appelo_params);
else
    appelo_alpha = 0;
    K_Na_H = 0; K_K_H = 0; K_Ca_H2 = 0; K_Mg_H2 = 0; K_Al_H3 = 0;
end

%%% Step 1 - transfer the complete field dissolved inventory to the
%%% air-dried soil slurry (Kanzaki et al. 2024).
rho_b_kg = rho_bulk / 1e3;                                  %% [kg/m3]
soil_mass_kg_m2 = rho_b_kg * Zrw;                           %% [kg/m2]
extract_water_L_m2 = bulk_pH_params.n_extract*soil_mass_kg_m2; %% [L/m2]
dil = bulk_pH_params.field_water_L_m2 / extract_water_L_m2;

H_d   = H_pw    * dil;
Ca_d  = Ca_pw   * dil;
Mg_d  = Mg_pw   * dil;
K_d   = Kion_pw * dil;
Na_d  = Na_pw   * dil;
Al_d  = (Al_pw  / Al_CEC_scale) * dil;  %% convert legacy scaled Al to aqueous Al
Al3_d = (Al3_pw / Al_CEC_scale) * dil;
NO3_d = NO3_pw  * dil;
PO4_d = PO4_pw  * dil;
NH4_d = NH4_pw  * dil;
An_d  = An_pw   * dil;

%%% Extract CO2 from soil air
CO2_ext = k_H * CO2_air_extract;

%%% CEC per litre of extract [umol_c/L]
cec_factor = CEC_tot / extract_water_L_m2;

%%% Charge-balance residual as function of final extract H [umol/L]
R_fn = @(H) charge_residual(H, Ca_d, Mg_d, K_d, Na_d, Al_d, ...
    NO3_d, PO4_d, NH4_d, An_d, CO2_ext, k1, k2, k_w, k_nh3, ...
    K1_Al, K2_Al, K3_Al, K4_Al, ...
    K_Ca_Mg, K_Ca_K, K_Ca_Na, K_Ca_Al, K_Ca_H, ...
    f_Ca_0, f_Mg_0, f_K_0, f_Na_0, f_Al_0, cec_factor, conv_mol, ...
    scheme_is_appelo, appelo_alpha, K_Na_H, K_K_H, K_Ca_H2, K_Mg_H2, K_Al_H3);

%%% Root selection/bracket anchor: H_extract_anchor_umolL when set (parameter input,
%%% keeps the model generalizable across pH ranges); falls back to H_pw
%%% otherwise.  Computed before scanning so a VALIDATED anchor can also
%%% narrow the scan itself (see below), the same warm-start-from-yesterday's
%%% -state approach ERW_function.m uses for its own daily solve
%%% (H_anchor/H_prev), rather than rescanning the full pH 0-14 range from
%%% scratch every call.
has_validated_anchor = isfield(bulk_pH_params,'H_extract_anchor_umolL') && ...
    isfinite(bulk_pH_params.H_extract_anchor_umolL) && bulk_pH_params.H_extract_anchor_umolL > 0;
if has_validated_anchor
    pH_anchor_scan = -log10(bulk_pH_params.H_extract_anchor_umolL / conv_mol);
else
    pH_anchor_scan = pH_pw;
end

%%% Bracket in pH space, then solve with Newton-Raphson.  A VALIDATED anchor
%%% (H_extract_anchor_umolL, e.g. the previous converged day) narrows the scan to a
%%% window around it first -- cheap, and safe because the anchor is known
%%% good.  A cold start keeps the full [2,12] scan so every sign-changing
%%% bracket is available and the one nearest the porewater pH can be selected.
%%% The narrow scan uses a coarser 0.1 pH spacing than the full scan: root
%%% accuracy comes from the Newton solve, so the scan only needs to catch the
%%% sign change in the anchor window and hand Newton a valid bracket.  Using
%%% the full-scan density (0.05 pH) here would roughly double the scan cost
%%% for the common tracked case without any accuracy benefit (cf. SCEPTER,
%%% which converges the charge balance by iteration and never grids pH).
full_pH_lo = 2; full_pH_hi = 12; full_n_points = 201;
narrow_step_pH = 0.1;
found = false;
scan_used = 'none';
anchor_pH_used = pH_anchor_scan;
all_bracket_pH = [];
if has_validated_anchor
    narrow_pH_lo = max(full_pH_lo, pH_anchor_scan - 2);
    narrow_pH_hi = min(full_pH_hi, pH_anchor_scan + 2);
    narrow_n_points = round((narrow_pH_hi - narrow_pH_lo) / narrow_step_pH) + 1;
    [pH_scan, H_scan, R_scan] = bulk_pH_residual_scan(R_fn, conv_mol, ...
        narrow_pH_lo, narrow_pH_hi, narrow_n_points);
    [H_a, H_b, found, all_bracket_pH] = ...
        bulk_pH_select_bracket(pH_scan, H_scan, R_scan, pH_anchor_scan);
    if found
        scan_used = 'narrow';
    end
end
if ~found
    [pH_scan, H_scan, R_scan] = bulk_pH_residual_scan(R_fn, conv_mol, ...
        full_pH_lo, full_pH_hi, full_n_points);
    [H_a, H_b, found, all_bracket_pH] = ...
        bulk_pH_select_bracket(pH_scan, H_scan, R_scan, pH_anchor_scan);
    if found
        scan_used = 'full';
    end
end
if found
    H0 = bulk_pH_initial_guess(bulk_pH_params, H_pw, conv_mol, H_a, H_b);
    [H, solver_diag] = newton_raphson_H(R_fn, H0, H_a, H_b);
else
    error('ERW:noBulkPHRoot', 'Could not solve bulk extract pH from pore-water chemistry.');
end

if ~isfinite(H) || H <= 0
    error('ERW:invalidState', 'Solved bulk extract H must be positive.');
end
[R_final, diag] = R_fn(H);
pH_bulk = diag.pH;
full_bracket_pH = [];
large_offset = abs(pH_bulk - pH_pw) > 2.5;
if large_offset
    %%% Check whether the narrow scan missed a pore-proximal root.
    [pH_scan_full, H_scan_full, R_scan_full] = bulk_pH_residual_scan(R_fn, conv_mol, ...
        full_pH_lo, full_pH_hi, full_n_points);
    [~, ~, ~, full_bracket_pH] = bulk_pH_select_bracket(pH_scan_full, ...
        H_scan_full, R_scan_full, pH_pw);
    full_bracket_pH = full_bracket_pH(:).';
    n_brackets_full = numel(full_bracket_pH);
    n_brackets = numel(all_bracket_pH);
    near_pore = '';
    if ~isempty(full_bracket_pH)
        [d_pore_min, i_pore] = min(abs(full_bracket_pH - pH_pw));
        near_pore = sprintf('; full[2,12] brackets=%d mid=%s nearest-to-pore pH %.2f (%.2f from pore)', ...
            n_brackets_full, full_bracket_list_str(full_bracket_pH), ...
            full_bracket_pH(i_pore), d_pore_min);
    end
    warning('ERW:bulkPHLargeOffset', ...
        ['Bulk pH (%.3f) deviates >2.5 units from pore-water pH (%.2f); ' ...
        'scan=%s anchor_pH=%.2f narrow_brackets=%d narrow_mid=%.2f%s; ' ...
        'check exchange parameters or root selection.'], ...
        pH_bulk, pH_pw, scan_used, anchor_pH_used, numel(all_bracket_pH), ...
        0.5*(-log10(H/conv_mol) + pH_bulk), near_pore);
end
if nargout > 1
    diag.charge_residual = R_final;
    diag.H = H;
    diag.solver = solver_diag.solver;
    diag.iterations = solver_diag.iterations;
    diag.converged = solver_diag.converged;
    diag.initial_H = solver_diag.initial_H;
    diag.newton_fallback_steps = solver_diag.fallback_steps;
    diag.derivative = solver_diag.derivative;
    diag.step_relative = solver_diag.step_relative;
    diag.bracket_H = solver_diag.bracket_H;
    diag.scan_used = scan_used;
    diag.anchor_pH = anchor_pH_used;
    diag.n_brackets = numel(all_bracket_pH);
    diag.bracket_mid_pH = all_bracket_pH;
    diag.full_bracket_mid_pH = full_bracket_pH;
    diag.large_offset = large_offset;
    diag.H_diluted = H_d;
    diag.Al3_diluted_input = Al3_d;
    diag.PO4_diluted_input = PO4_d;
    diag.CO2_air_extract = CO2_air_extract;
    diag.CO2_ext = CO2_ext;
    diag.pH_bulk = pH_bulk;
    diag.pH_pore = pH_pw;

    diag.K_CEC = bulk_pH_params.K_CEC;
    diag.f_CEC_initial = [f_Ca_0, f_Mg_0, f_K_0, f_Na_0, f_Al_0, f_H_0];
    diag.dilution_factor = dil;
    diag.cec_factor = cec_factor;
    diag.field_water_L_m2 = bulk_pH_params.field_water_L_m2;
    diag.extract_water_L_m2 = extract_water_L_m2;
    diag.n_extract = bulk_pH_params.n_extract;
    diag.CEC_tot = CEC_tot;
    diag.rho_bulk = rho_bulk;
    diag.Zrw = Zrw;
    diag.Al_CEC_scale = Al_CEC_scale;
end

return
end

%% =======================================================================
function s = full_bracket_list_str(bracket_pH)
bracket_pH = bracket_pH(:).';
if isempty(bracket_pH)
    s = '';
    return
end
if numel(bracket_pH) <= 4
    s = sprintf('%.2f,', bracket_pH);
    s = s(1:end-1);
else
    s = sprintf('%.2f..%.2f n=%d', bracket_pH(1), bracket_pH(end), numel(bracket_pH));
end
end

%% =======================================================================
function H0 = bulk_pH_initial_guess(bulk_pH_params,H_pw,conv_mol,H_a,H_b)
H0 = H_pw;
if isfield(bulk_pH_params,'H_init') && isfinite(bulk_pH_params.H_init) && bulk_pH_params.H_init > 0
    H0 = bulk_pH_params.H_init;
elseif isfield(bulk_pH_params,'H_prev') && isfinite(bulk_pH_params.H_prev) && bulk_pH_params.H_prev > 0 ...
        && isfield(bulk_pH_params,'H_prev_valid') && bulk_pH_params.H_prev_valid
    H0 = bulk_pH_params.H_prev;
elseif isfield(bulk_pH_params,'pH_init') && isfinite(bulk_pH_params.pH_init)
    H0 = 10^(-bulk_pH_params.pH_init) * conv_mol;
elseif isfield(bulk_pH_params,'pH_prev') && isfinite(bulk_pH_params.pH_prev)
    H0 = 10^(-bulk_pH_params.pH_prev) * conv_mol;
end

H_min = min(H_a,H_b);
H_max = max(H_a,H_b);
if ~isfinite(H0) || H0 <= H_min || H0 >= H_max
    H0 = 0.5*(H_a + H_b);
end
end

%% =======================================================================
function [pH_scan, H_scan, R_scan] = bulk_pH_residual_scan(R_fn, conv_mol, pH_lo, pH_hi, n_points)
%%% Evaluate the charge residual over a pH window.  Called first on a
%%% narrow, anchor-centred window and only on the full [2,12] range as a
%%% fallback, since the narrow scan is ~5x cheaper and covers the common
%%% case where the root sits near the anchor.
pH_scan = linspace(pH_lo, pH_hi, n_points);
H_scan = 10.^(-pH_scan) * conv_mol;
R_scan = nan(size(H_scan));
for ii = 1:numel(H_scan)
    try
        R_scan(ii) = R_fn(H_scan(ii));
    catch ME
        if startsWith(ME.identifier, 'ERW:')
            R_scan(ii) = NaN;
        else
            rethrow(ME)
        end
    end
end
end

%% =======================================================================
function [H_a, H_b, found, bracket_mid_pH] = bulk_pH_select_bracket(pH_scan, H_scan, R_scan, pH_anchor_scan)
% Select the sign-changing bracket nearest the pH anchor.
H_a = NaN; H_b = NaN; found = false;
bracket_mid_pH = [];
valid_scan = isfinite(R_scan);
sign_change = find(valid_scan(1:end-1) & valid_scan(2:end) & ...
    (R_scan(1:end-1) == 0 | R_scan(2:end) == 0 | sign(R_scan(1:end-1)) ~= sign(R_scan(2:end))));
if isempty(sign_change)
    return
end
bracket_mid_pH = 0.5 * (pH_scan(sign_change) + pH_scan(sign_change+1));
[~,i_best] = min(abs(pH_scan(sign_change) - pH_anchor_scan));
H_a = H_scan(sign_change(i_best));
H_b = H_scan(sign_change(i_best)+1);
found = true;
end

%% =======================================================================
function [H, solver_diag] = newton_raphson_H(R_fn,H0,H_a,H_b)
%%% Newton-Raphson on the complete coupled residual; bracket keeps the
%%% same root selected by the original pH scan.

max_iter = 100;
abs_tol  = 1e-8;
rel_tol  = 1e-8;
rel_charge_tol = 1e-7;
fallback_steps = 0;
last_derivative = NaN;
step_rel = Inf;
converged = false;
newton_failed = false;
solver_name = 'newton_raphson';

H_left = H_a;
H_right = H_b;
[R_left, ~] = R_fn(H_left);
R_right = R_fn(H_right);
if ~(isfinite(R_left) && isfinite(R_right)) || sign(R_left) == sign(R_right)
    error('ERW:noBulkPHRoot', 'Could not solve bulk extract pH from pore-water chemistry.');
end

if abs(R_left) < abs_tol
    H = H_left;
    converged = true;
    iter = 0;
elseif abs(R_right) < abs_tol
    H = H_right;
    converged = true;
    iter = 0;
else
    H = H0;
    R = R_fn(H);
    if ~isfinite(R)
        H = 0.5*(H_left + H_right);
        R = R_fn(H);
        fallback_steps = fallback_steps + 1;
    end
    if sign(R) == sign(R_left)
        H_left = H; R_left = R;
    else
        H_right = H;
    end

    for iter = 1:max_iter
        [dR_dH,~] = finite_difference_derivative(R_fn,H,H_left,H_right);
        last_derivative = dR_dH;
        H_min = min(H_left,H_right);
        H_max = max(H_left,H_right);
        H_trial = NaN;
        if isfinite(dR_dH) && dR_dH ~= 0
            % Newton step in log-H (x=ln H, dR/dx=dR/dH*H): H stays positive
            % and a large pH excursion (e.g. calcite addition) is damped to
            % <=1 pH unit per step.
            dlnH = max(min(-R/(dR_dH*H), log(10)), -log(10));
            H_trial = H*exp(dlnH);
        end
        if ~isfinite(H_trial) || H_trial <= H_min || H_trial >= H_max
            H_trial = 0.5*(H_left + H_right);
            fallback_steps = fallback_steps + 1;
        end

        [R_trial, state_trial] = R_fn(H_trial);
        if ~isfinite(R_trial)
            %% Newton wandered into a non-finite region.  The maintained
            %% [H_left,H_right] bracket still has finite, opposite-sign
            %% endpoints, so hand off to the guaranteed bisection fallback.
            newton_failed = true;
            break;
        end
        charge_scale = state_trial.charge_scale;
        step_rel = abs(H_trial - H) / max(abs(H_trial),eps);
        H = H_trial;
        R = R_trial;

        if abs(R) < abs_tol || abs(R) < charge_scale * rel_charge_tol || step_rel < rel_tol
            converged = true;
            break;
        end
        if sign(R) == sign(R_left)
            H_left = H; R_left = R;
        else
            H_right = H;
        end
    end
end

if ~newton_failed && ~converged && abs(H_right-H_left)/max(abs(H),eps) < rel_tol
    converged = true;
end
if ~converged
    %% [H_left,H_right] still brackets the root (opposite-sign residuals),
    %% so bisection converges without aborting the run.
    [H, R, bisect_iter] = bisect_bracket_H(R_fn, H_left, H_right, ...
        abs_tol, rel_tol, rel_charge_tol);
    iter = iter + bisect_iter;
    converged = true;
    solver_name = 'bisection_fallback';
end

solver_diag.solver = solver_name;
solver_diag.iterations = iter;
solver_diag.converged = converged;
solver_diag.initial_H = H0;
solver_diag.fallback_steps = fallback_steps;
solver_diag.derivative = last_derivative;
solver_diag.step_relative = step_rel;
solver_diag.bracket_H = [H_left H_right];
end

%% =======================================================================
function [H, R, iter] = bisect_bracket_H(R_fn, H_left, H_right, abs_tol, rel_tol, rel_charge_tol)
%%% Bisection on a sign-bracketed charge residual; unconditionally
%%% convergent.  Used as the no-abort fallback when Newton-Raphson stalls
%%% or wanders non-finite.  Uses the same tolerances as the Newton solve.
[R_left, ~]  = R_fn(H_left);
[R_right, ~] = R_fn(H_right);
if ~(isfinite(R_left) && isfinite(R_right)) || sign(R_left) == sign(R_right)
    %% No usable sign bracket (should not occur for a pH-scan bracket):
    %% return the finite endpoint with the smaller residual magnitude.
    if isfinite(R_left) && (~isfinite(R_right) || abs(R_left) <= abs(R_right))
        H = H_left; R = R_left;
    else
        H = H_right; R = R_right;
    end
    iter = 0;
    return
end
H = 0.5*(H_left + H_right);
R = R_left;
for iter = 1:200
    H = 0.5*(H_left + H_right);
    [R, state] = R_fn(H);
    if ~isfinite(R)
        %% Discard the half toward the non-finite midpoint, keep H_left side.
        H_right = H;
        continue
    end
    charge_scale = state.charge_scale;
    if abs(R) < abs_tol || abs(R) < charge_scale*rel_charge_tol || ...
            abs(H_right - H_left)/max(abs(H),eps) < rel_tol
        return
    end
    if sign(R) == sign(R_left)
        H_left = H; R_left = R;
    else
        H_right = H;
    end
end
end

%% =======================================================================
function [dR_dH,dH] = finite_difference_derivative(R_fn,H,H_left,H_right)
H_min = min(H_left,H_right);
H_max = max(H_left,H_right);
dH = max(sqrt(eps)*max(abs(H),1), 1e-10*max(abs(H),1));
H_plus = min(H + dH,H_max);
H_minus = max(H - dH,H_min);
if H_plus == H_minus
    dR_dH = NaN;
    return
end
R_plus = R_fn(H_plus);
R_minus = R_fn(H_minus);
if ~(isfinite(R_plus) && isfinite(R_minus))
    dR_dH = NaN;
else
    dR_dH = (R_plus - R_minus)/(H_plus - H_minus);
end
end

%% =======================================================================
function [R, state] = charge_residual(H, Ca_d, Mg_d, K_d, Na_d, Al_d, ...
    NO3_d, PO4_d, NH4_d, An_d, CO2_ext, k1, k2, k_w, ~, ...
    K1_Al, K2_Al, K3_Al, K4_Al, ...
    K_Ca_Mg, K_Ca_K, K_Ca_Na, K_Ca_Al, K_Ca_H, ...
    f_Ca_0, f_Mg_0, f_K_0, f_Na_0, f_Al_0, cec_factor, conv_mol, ...
    scheme_is_appelo, appelo_alpha, K_Na_H, K_K_H, K_Ca_H2, K_Mg_H2, K_Al_H3)

%%% Carbonate at extract pCO2; concentration-based.
HCO3 = k1 * CO2_ext / H;
CO3  = k2 * HCO3 / H;
OH   = k_w / H;

%%% Exchange equilibrium and mass conservation on final extract solution.
state = exchange_state(H, Ca_d, Mg_d, K_d, Na_d, NH4_d, Al_d, ...
    K1_Al, K2_Al, K3_Al, K4_Al, ...
    K_Ca_Mg, K_Ca_K, K_Ca_Na, K_Ca_Al, K_Ca_H, ...
    f_Ca_0, f_Mg_0, f_K_0, f_Na_0, f_Al_0, cec_factor, ...
    scheme_is_appelo, appelo_alpha, K_Na_H, K_K_H, K_Ca_H2, K_Mg_H2, K_Al_H3);

%%% Al speciation
Al3_eq = state.Al3;
AlOH_eq = state.AlOH;
AlOH2_eq = state.AlOH2;
AlOH4_eq = state.AlOH4;

%%% Charge balance [umolc/L].  Exchange partitioning already transfers the aqueous
%%% inventories below; adsorbed site fractions do not enter electroneutrality.
%%% The daily TC-SMEW alkalinity state tracks total ammoniacal N as NH4+.
%%% At the Zhang pH this differs negligibly from explicit NH3 speciation.
NH4_charged = state.NH4;
cations = H + 2*state.Ca + 2*state.Mg + state.K + state.Na + NH4_charged ...
    + 3*Al3_eq + 2*AlOH_eq + AlOH2_eq;
anions  = OH + HCO3 + 2*CO3 + NO3_d + PO4_d + An_d + AlOH4_eq;
R = cations - anions;
state.charge_scale = max(abs(cations) + abs(anions), eps);
state.HCO3 = HCO3;
state.CO3 = CO3;
state.OH = OH;
state.pH = -log10(H / conv_mol);

end

%% =======================================================================
function [Al3,AlOH,AlOH2,AlOH3,AlOH4] = al_hydrolysis_species(...
    Al_eq,H,K1_Al,K2_Al,K3_Al,K4_Al)
%%% AlOH is AlOH2+; AlOH2 is Al(OH)2+; AlOH3 is neutral; AlOH4 is anion.
if ~isfinite(Al_eq) || Al_eq < 0
    Al3 = 0; AlOH = 0; AlOH2 = 0; AlOH3 = 0; AlOH4 = 0;
    return
end
if ~isfinite(H) || H <= 0
    Al3 = Al_eq; AlOH = 0; AlOH2 = 0; AlOH3 = 0; AlOH4 = 0;
    return
end
r1 = K1_Al / H;
r2 = K1_Al*K2_Al / H^2;
r3 = K1_Al*K2_Al*K3_Al / H^3;
r4 = K1_Al*K2_Al*K3_Al*K4_Al / H^4;
den = 1 + r1 + r2 + r3 + r4;
if ~(isfinite(den) && den > 0)
    Al3 = Al_eq; AlOH = 0; AlOH2 = 0; AlOH3 = 0; AlOH4 = 0;
    return
end

Al3   = Al_eq / den;
AlOH  = r1*Al3;
AlOH2 = r2*Al3;
AlOH3 = r3*Al3;
AlOH4 = max(Al_eq - (Al3 + AlOH + AlOH2 + AlOH3), 0);  % floor prevents float-cancellation underflow
end

%% =======================================================================
function state = exchange_state(H, Ca_d, Mg_d, K_d, Na_d, NH4_d, Al_d, ...
    K1_Al, K2_Al, K3_Al, K4_Al, ...
    K_Ca_Mg, K_Ca_K, K_Ca_Na, K_Ca_Al, K_Ca_H, ...
    f_Ca_0, f_Mg_0, f_K_0, f_Na_0, f_Al_0, cec_factor, ...
    scheme_is_appelo, appelo_alpha, K_Na_H, K_K_H, K_Ca_H2, K_Mg_H2, K_Al_H3)

[f_Ca, f_Mg, f_K, f_Na, f_Al, f_H, Ca_eq, Mg_eq, K_eq, Na_eq, Al_eq, mass_residuals] = ...
    partition_exchange_inventory(H, Ca_d, Mg_d, K_d, Na_d, Al_d, ...
    K1_Al, K2_Al, K3_Al, K4_Al, ...
    K_Ca_Mg, K_Ca_K, K_Ca_Na, K_Ca_Al, K_Ca_H, ...
    f_Ca_0, f_Mg_0, f_K_0, f_Na_0, f_Al_0, cec_factor, ...
    scheme_is_appelo, appelo_alpha, K_Na_H, K_K_H, K_Ca_H2, K_Mg_H2, K_Al_H3);

[Al3_eq,AlOH_eq,AlOH2_eq,AlOH3_eq,AlOH4_eq] = al_hydrolysis_species(...
    Al_eq,H,K1_Al,K2_Al,K3_Al,K4_Al);
state.f_Ca = f_Ca;
state.f_Mg = f_Mg;
state.f_K = f_K;
state.f_Na = f_Na;
state.f_Al = f_Al;
state.f_H = f_H;
state.site_sum = f_Ca + f_Mg + f_K + f_Na + f_Al + f_H;
state.Ca = Ca_eq;
state.Mg = Mg_eq;
state.K = K_eq;
state.Na = Na_eq;
state.NH4 = NH4_d;
state.Al = Al_eq;
state.Al3 = Al3_eq;
state.AlOH = AlOH_eq;
state.AlOH2 = AlOH2_eq;
state.AlOH3 = AlOH3_eq;
state.AlOH4 = AlOH4_eq;
if scheme_is_appelo
    gamma_app = 10^(appelo_alpha*f_H);
    state.exchange_residuals = [ ...
        f_Na - K_Na_H * Na_eq * (f_H/H) * gamma_app, ...
        f_K  - K_K_H  * K_eq  * (f_H/H) * gamma_app, ...
        f_Ca - K_Ca_H2* Ca_eq * (f_H/H)^2 * gamma_app^2, ...
        f_Mg - K_Mg_H2* Mg_eq * (f_H/H)^2 * gamma_app^2, ...
        f_Al - K_Al_H3* Al3_eq * (f_H/H)^3 * gamma_app^3];
else
    state.exchange_residuals = [ ...
        f_Mg - (Mg_eq/K_Ca_Mg) * f_Ca / Ca_eq, ...
        f_Na - Na_eq * sqrt(f_Ca / (K_Ca_Na * Ca_eq)), ...
        f_K  - K_eq  * sqrt(f_Ca / (K_Ca_K  * Ca_eq)), ...
        f_Al - (Al3_eq) * sqrt(f_Ca^3 / (K_Ca_Al * Ca_eq^3)), ...
        f_H - H * sqrt(f_Ca / (K_Ca_H * Ca_eq))];
end
state.exchange_mass_residuals_umol_L = mass_residuals;
end

%% =======================================================================
function [f_Ca, f_Mg, f_K, f_Na, f_Al, f_H, Ca, Mg, K, Na, Al, mass_residuals] = ...
    partition_exchange_inventory(H, Ca_d, Mg_d, K_d, Na_d, Al_d, ...
    K1_Al, K2_Al, K3_Al, K4_Al, ...
    K_Ca_Mg, K_Ca_K, K_Ca_Na, K_Ca_Al, K_Ca_H, ...
    f_Ca_0, f_Mg_0, f_K_0, f_Na_0, f_Al_0, cec_factor, ...
    scheme_is_appelo, appelo_alpha, K_Na_H, K_K_H, K_Ca_H2, K_Mg_H2, K_Al_H3)
if scheme_is_appelo
    [f_Ca, f_Mg, f_K, f_Na, f_Al, f_H, Ca, Mg, K, Na, Al, mass_residuals] = ...
        partition_exchange_inventory_appelo(H, Ca_d, Mg_d, K_d, Na_d, Al_d, ...
        K1_Al, K2_Al, K3_Al, K4_Al, ...
        f_Ca_0, f_Mg_0, f_K_0, f_Na_0, f_Al_0, cec_factor, ...
        appelo_alpha, K_Na_H, K_K_H, K_Ca_H2, K_Mg_H2, K_Al_H3);
    return
end
%%% Conserved-inventory Gaines-Thomas partition of the air-dry extract.
%%%
%%% The initial dissolved concentration and the exchangeable pool of each
%%% cation are one conserved inventory.  Exchange fractions are charge
%%% fractions, so divalent Ca/Mg are divided by 2 and trivalent Al by 3
%%% before adding the exchangeable charge cec_factor*f_X_0:
%%%     X_total = X_d + cec_factor * f_X_0 / z_X   [umol/L]
%%% At each trial pore-water H a single free variable f_Ca fixes dissolved
%%% Ca, and the Gaines-Thomas relations plus each cation's own mass balance
%%% give the other fractions; the Ca-H Gaines-Thomas relation
%%% gives f_H.  The scalar closure sum(f_X) = 1 is solved by bisection over
%%% the physical window f_Ca in (0, 2*Ca_total/cec_factor), which keeps
%%% dissolved Ca >= 0 without any concentration floor.  There is no fallback
%%% exchange state: an inventory that cannot fill the sites is a genuine
%%% infeasibility and raises an error rather than returning the seed.  This
%%% replaces the earlier reference-cation-switching routing, whose floors
%%% and seed fallback could manufacture mass at low dissolved concentrations.
if ~(isfinite(cec_factor) && cec_factor > 0)
    error('ERW:partitionExchange', 'cec_factor must be positive and finite.');
end
K_CEC = [K_Ca_Mg, K_Ca_K, K_Ca_Na, K_Ca_Al, K_Ca_H];
if any(~isfinite(K_CEC)) || any(K_CEC <= 0)
    error('ERW:partitionExchange', ...
        'Gaines-Thomas selectivity constants must be positive and finite.');
end
dissolved_in = [Ca_d, Mg_d, K_d, Na_d, Al_d];
if any(~isfinite(dissolved_in)) || any(dissolved_in < 0)
    error('ERW:partitionExchange', ...
        'Dissolved extract concentrations must be finite and nonnegative.');
end
f_initial = [f_Ca_0, f_Mg_0, f_K_0, f_Na_0, f_Al_0];
if any(~isfinite(f_initial)) || any(f_initial < 0)
    error('ERW:partitionExchange', ...
        'Initial exchange fractions must be finite and nonnegative.');
end

Al3_frac = al3_fraction(H, K1_Al, K2_Al, K3_Al, K4_Al);
Ca_total = Ca_d + cec_factor*f_Ca_0/2;
Mg_total = Mg_d + cec_factor*f_Mg_0/2;
K_total  = K_d  + cec_factor*f_K_0;
Na_total = Na_d + cec_factor*f_Na_0;
Al_total = Al_d + cec_factor*f_Al_0/3;
if ~(isfinite(Ca_total) && Ca_total > 0)
    error('ERW:partitionExchange', 'The conserved Ca inventory must be positive.');
end

%%% Dissolved Ca and the proton fraction remain physical only below
%%% their respective f_Ca bounds.
f_lo = 1e-12;
f_inventory_hi = 2*Ca_total/cec_factor;
f_H_hi = K_Ca_H*Ca_total/(H^2 + K_Ca_H*cec_factor/2);
f_hi = min([1 - 1e-12, f_inventory_hi*(1 - 1e-12), ...
    f_H_hi*(1 - 1e-12)]);
if f_hi <= f_lo
    error('ERW:partitionExchange', ...
        'The conserved Ca inventory cannot occupy a positive exchange fraction.');
end

site_residual = @(x) partition_site_residual(x, H, Ca_total, Mg_total, ...
    K_total, Na_total, Al_total, cec_factor, K_Ca_Mg, K_Ca_K, K_Ca_Na, ...
    K_Ca_Al, K_Ca_H, Al3_frac);
g_lo = site_residual(f_lo);
g_hi = site_residual(f_hi);
if ~(isfinite(g_lo) && isfinite(g_hi) && sign(g_lo) ~= sign(g_hi))
    error('ERW:partitionExchange', ...
        'No physical Gaines-Thomas partition satisfies the conserved inventories.');
end
for k = 1:200
    f_mid = 0.5*(f_lo + f_hi);
    g_mid = site_residual(f_mid);
    if abs(g_mid) < 1e-12 || (f_hi - f_lo)/2 < 1e-15
        break;
    end
    if sign(g_mid) == sign(g_lo)
        f_lo = f_mid; g_lo = g_mid;
    else
        f_hi = f_mid;
    end
end
f_Ca = 0.5*(f_lo + f_hi);
[~, f_Mg, f_K, f_Na, f_Al, f_H] = site_residual(f_Ca);

Ca = Ca_total - cec_factor*f_Ca/2;
Mg = Mg_total - cec_factor*f_Mg/2;
K  = K_total  - cec_factor*f_K;
Na = Na_total - cec_factor*f_Na;
Al = Al_total - cec_factor*f_Al/3;

fractions = [f_Ca, f_Mg, f_K, f_Na, f_Al, f_H];
concentrations = [Ca, Mg, K, Na, Al];
conc_tol = 1e-6 * max([1, abs(Ca_total), abs(Mg_total), abs(K_total), ...
    abs(Na_total), abs(Al_total)]);
if any(~isfinite([fractions, concentrations])) || any(fractions < -1e-10) || ...
        abs(sum(fractions) - 1) > 1e-8 || any(concentrations < -conc_tol)
    error('ERW:partitionExchange', ...
        'Conserved exchange partition produced an unphysical state.');
end

%%% Clamp floating-point round-off only; a genuinely negative state errored
%%% above.
fractions = max(fractions, 0);
f_Ca = fractions(1); f_Mg = fractions(2); f_K = fractions(3);
f_Na = fractions(4); f_Al = fractions(5); f_H = fractions(6);
Ca = max(Ca, 0); Mg = max(Mg, 0); K = max(K, 0); Na = max(Na, 0); Al = max(Al, 0);

%%% Per-cation conservation residual (dissolved + exchangeable - initial
%%% total).  Zero to round-off by construction; returned as a solve
%%% diagnostic so a downstream check can confirm no mass was created.
mass_residuals = [ ...
    (Ca + cec_factor*f_Ca/2) - Ca_total, ...
    (Mg + cec_factor*f_Mg/2) - Mg_total, ...
    (K  + cec_factor*f_K)    - K_total, ...
    (Na + cec_factor*f_Na)   - Na_total, ...
    (Al + cec_factor*f_Al/3) - Al_total];
end

%% =======================================================================
function [residual, f_Mg, f_K, f_Na, f_Al, f_H] = partition_site_residual( ...
    f_Ca, H, Ca_total, Mg_total, K_total, Na_total, Al_total, cec_factor, ...
    K_Ca_Mg, K_Ca_K, K_Ca_Na, K_Ca_Al, K_Ca_H, Al3_frac)
%%% Exchange fractions and their site-sum residual at one trial f_Ca, with
%%% each cation's dissolved concentration reduced by its own exchange draw
%%% (mass balance).  Same Gaines-Thomas relations the daily ERW
%%% solver and the extract charge balance enforce.
Ca = Ca_total - cec_factor*f_Ca/2;
if ~(isfinite(Ca) && Ca > 0)
    residual = NaN;
    f_Mg = NaN; f_K = NaN; f_Na = NaN; f_Al = NaN; f_H = NaN;
    return
end
A_Mg = f_Ca / (K_Ca_Mg * Ca);
A_K  = sqrt(f_Ca / (K_Ca_K  * Ca));
A_Na = sqrt(f_Ca / (K_Ca_Na * Ca));
A_Al = Al3_frac * sqrt(f_Ca^3 / (K_Ca_Al * Ca^3));
f_Mg = A_Mg * Mg_total / (1 + A_Mg*cec_factor/2);
f_K  = A_K  * K_total  / (1 + A_K *cec_factor);
f_Na = A_Na * Na_total / (1 + A_Na*cec_factor);
f_Al = A_Al * Al_total / (1 + A_Al*cec_factor/3);
A_H  = H * sqrt(f_Ca / (K_Ca_H * Ca));
f_H  = A_H;
residual = f_Ca + f_Mg + f_K + f_Na + f_Al + f_H - 1;
end

%% =======================================================================
function [f_Ca, f_Mg, f_K, f_Na, f_Al, f_H, Ca, Mg, K, Na, Al, mass_residuals] = ...
    partition_exchange_inventory_appelo(H, Ca_d, Mg_d, K_d, Na_d, Al_d, ...
    K1_Al, K2_Al, K3_Al, K4_Al, ...
    f_Ca_0, f_Mg_0, f_K_0, f_Na_0, f_Al_0, cec_factor, ...
    appelo_alpha, K_Na_H, K_K_H, K_Ca_H2, K_Mg_H2, K_Al_H3)
%%% SCEPTER/Appelo (1994) partition of the air-dry extract.
%%% gamma = 10^(alpha*f_H):
%%%     f[X-Na]  = K_Na *[Na] *(f_H/[H] )  * gamma
%%%     f[X-K]   = K_K  *[K]  *(f_H/[H] )  * gamma
%%%     f[X2-Ca] = K_Ca *[Ca] *(f_H/[H] )^2* gamma^2
%%%     f[X2-Mg] = K_Mg *[Mg] *(f_H/[H] )^2* gamma^2
%%%     f[X3-Al] = K_Al *[Al3]*(f_H/[H] )^3* gamma^3
%%%     site sum  f_H + f[X-Na]+f[X-K]+f[X2-Ca]+f[X2-Mg]+f[X3-Al] = 1
%%% Constants are converted from mol/L to the model's umol/L basis in
%%% Resolve_Appelo_Constants. Each dissolved plus exchange inventory is held.
if ~(isfinite(cec_factor) && cec_factor > 0)
    error('ERW:partitionExchange', 'cec_factor must be positive and finite.');
end
if any(~isfinite([K_Na_H K_K_H K_Ca_H2 K_Mg_H2 K_Al_H3])) || ...
        any([K_Na_H K_K_H K_Ca_H2 K_Mg_H2 K_Al_H3] <= 0)
    error('ERW:partitionExchange', ...
        'Appelo selectivity constants must be positive and finite.');
end
if ~(isfinite(appelo_alpha) && appelo_alpha >= 0)
    error('ERW:partitionExchange', 'Appelo alpha must be a finite nonnegative scalar.');
end
dissolved_in = [Ca_d, Mg_d, K_d, Na_d, Al_d];
if any(~isfinite(dissolved_in)) || any(dissolved_in < 0)
    error('ERW:partitionExchange', ...
        'Dissolved extract concentrations must be finite and nonnegative.');
end
f_initial = [f_Ca_0, f_Mg_0, f_K_0, f_Na_0, f_Al_0];
if any(~isfinite(f_initial)) || any(f_initial < 0)
    error('ERW:partitionExchange', ...
        'Initial exchange fractions must be finite and nonnegative.');
end

Al3_frac = al3_fraction(H, K1_Al, K2_Al, K3_Al, K4_Al);
Ca_total = Ca_d + cec_factor*f_Ca_0/2;
Mg_total = Mg_d + cec_factor*f_Mg_0/2;
K_total  = K_d  + cec_factor*f_K_0;
Na_total = Na_d + cec_factor*f_Na_0;
Al_total = Al_d + cec_factor*f_Al_0/3;
if ~(isfinite(Ca_total) && Ca_total > 0)
    error('ERW:partitionExchange', 'The conserved Ca inventory must be positive.');
end

site_residual = @(x) appelo_site_residual(x, H, Ca_total, Mg_total, ...
    K_total, Na_total, Al_total, cec_factor, appelo_alpha, ...
    K_Na_H, K_K_H, K_Ca_H2, K_Mg_H2, K_Al_H3, Al3_frac);

%%% Solve the site balance over f_H in (0,1).
f_lo = 1e-24;
f_hi = 1 - 1e-14;
g_lo = site_residual(f_lo);
g_hi = site_residual(f_hi);
if ~(isfinite(g_lo) && isfinite(g_hi) && sign(g_lo) ~= sign(g_hi))
    error('ERW:partitionExchange', ...
        'No physical Appelo partition satisfies the conserved inventories.');
end
for k = 1:300
    f_mid = 0.5*(f_lo + f_hi);
    g_mid = site_residual(f_mid);
    if abs(g_mid) < 1e-12
        break;
    end
    if sign(g_mid) == sign(g_lo)
        f_lo = f_mid; g_lo = g_mid;
    else
        f_hi = f_mid;
    end
end
f_H = 0.5*(f_lo + f_hi);
[~, f_Na, f_K, f_Ca, f_Mg, f_Al] = site_residual(f_H);

Ca = Ca_total - cec_factor*f_Ca/2;
Mg = Mg_total - cec_factor*f_Mg/2;
K  = K_total  - cec_factor*f_K;
Na = Na_total - cec_factor*f_Na;
Al = Al_total - cec_factor*f_Al/3;

fractions = [f_Ca, f_Mg, f_K, f_Na, f_Al, f_H];
concentrations = [Ca, Mg, K, Na, Al];
conc_tol = 1e-6 * max([1, abs(Ca_total), abs(Mg_total), abs(K_total), ...
    abs(Na_total), abs(Al_total)]);
if any(~isfinite([fractions, concentrations])) || any(fractions < -1e-10) || ...
        abs(sum(fractions) - 1) > 1e-8 || any(concentrations < -conc_tol)
    error('ERW:partitionExchange', ...
        'Conserved Appelo partition produced an unphysical state.');
end

fractions = max(fractions, 0);
f_Ca = fractions(1); f_Mg = fractions(2); f_K = fractions(3);
f_Na = fractions(4); f_Al = fractions(5); f_H = fractions(6);
Ca = max(Ca, 0); Mg = max(Mg, 0); K = max(K, 0); Na = max(Na, 0); Al = max(Al, 0);

mass_residuals = [ ...
    (Ca + cec_factor*f_Ca/2) - Ca_total, ...
    (Mg + cec_factor*f_Mg/2) - Mg_total, ...
    (K  + cec_factor*f_K)    - K_total, ...
    (Na + cec_factor*f_Na)   - Na_total, ...
    (Al + cec_factor*f_Al/3) - Al_total];
end

%% =======================================================================
function [residual, f_Na, f_K, f_Ca, f_Mg, f_Al] = appelo_site_residual( ...
    f_H, H, Ca_total, Mg_total, K_total, Na_total, Al_total, cec_factor, ...
    appelo_alpha, K_Na_H, K_K_H, K_Ca_H2, K_Mg_H2, K_Al_H3, Al3_frac)
%%% Appelo site balance at one f_H, retaining each cation inventory.
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
A_Al = K_Al_H3 * Al3_frac * rho3;
f_Na = A_Na * Na_total / (1 + A_Na*cec_factor);
f_K  = A_K  * K_total  / (1 + A_K *cec_factor);
f_Ca = A_Ca * Ca_total / (1 + A_Ca*cec_factor/2);
f_Mg = A_Mg * Mg_total / (1 + A_Mg*cec_factor/2);
f_Al = A_Al * Al_total / (1 + A_Al*cec_factor/3);
residual = f_H + f_Na + f_K + f_Ca + f_Mg + f_Al - 1;
end

%% =======================================================================
function a3 = al3_fraction(H,K1_Al,K2_Al,K3_Al,K4_Al)
r1 = K1_Al / H;
r2 = K1_Al*K2_Al / H^2;
r3 = K1_Al*K2_Al*K3_Al / H^3;
r4 = K1_Al*K2_Al*K3_Al*K4_Al / H^4;
den = 1 + r1 + r2 + r3 + r4;
if ~(isfinite(den) && den > 0)
    a3 = 1;
else
    a3 = 1/den;
end
end
