%#########################################################################
%   Subfunction POREWATER_TO_BULK_pH                                    #
%#########################################################################
%%% References  Kanzaki et al. 2024;  ISO 10390
function [pH_bulk, diag] = Porewater_to_Bulk_pH(varargin)
%%% Ionic-strength closure of the 1:n extract.  Porewater_to_Bulk_pH_at_I
%%% solves the slurry, open to the laboratory air, with its constants
%%% conditioned at a given ionic strength; this repeats it from the same
%%% inputs until the extract's own ionic strength is stationary, then
%%% reports the ACTIVITY pH, -log10(g1*[H+]).  That is what the electrode
%%% measured, and what the observed bulk pH and the T&C pH response
%%% functions refer to.  diag.H stays the extract's H+ concentration (the
%%% next day's root anchor); diag.pH_conc is -log10([H+]).
if nargin < 2
    error('ERW:invalidParam','Porewater_to_Bulk_pH requires the full argument list.');
end
bulk_pH_params = varargin{end-1};
if ~isstruct(bulk_pH_params)
    error('ERW:invalidParam','bulk_pH_params must be a struct.');
end
I_extract = 0;
extract_closed = false;
for extract_pass = 1:60
    bulk_pH_params.extract_ionic_strength = I_extract;
    varargin{end-1} = bulk_pH_params;
    [pH_conc, diag] = Porewater_to_Bulk_pH_at_I(varargin{:});
    I_next = diag.ionic_str_extract;
    %%% 1e-6 of I sits above the Newton solve's own noise, which can hold
    %%% the iteration in a two-point cycle about 1e-9 of I apart, and moves
    %%% the reported pH by under 1e-7.
    if abs(I_next-I_extract) <= 1e-14+1e-6*I_next
        extract_closed = true;
        break
    end
    I_extract = I_next;
end
if ~extract_closed
    error('ERW:extractActivityClosure', ...
        'Extract ionic strength did not close in %d passes (%.6g -> %.6g mol/L).', ...
        extract_pass,I_extract,I_next);
end
gamma_mono_extract = ERW_Activity_Coefficient(1,I_extract);
pH_bulk = pH_conc-log10(gamma_mono_extract);
diag.pH_conc = pH_conc;
diag.pH = pH_bulk;
diag.pH_bulk_raw = pH_bulk;
diag.pH_bulk = pH_bulk;
diag.extract_activity_I = I_extract;
diag.gamma_mono_extract = gamma_mono_extract;
diag.extract_activity_passes = extract_pass;
end

%%% Open-slurry solve at one extract ionic strength.
function [pH_bulk, diag] = Porewater_to_Bulk_pH_at_I(...
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
%%% CO2_air_extract                      atmospheric CO2 [umol per litre of air, T&C units:
%%%                                       ppm/22.41].  The paste stands open to this air, so
%%%                                       its dissolved CO2 is k_H(25 C)*CO2_air_extract.
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
%%%   .DIC_field_umol_L field dissolved inorganic carbon [umol/L]; it does
%%%                     not set the reading (the paste is open to the air)
%%%                     and enters only diag.CO2_air_exchange_umol_L
%%%   .CEC_cmolc_kg     prescribed reference CEC [cmolc/kg]
%%%   .dCEC_dpH         local CEC slope [cmolc/kg per pH unit]
%%%   .CEC_ref_pH       pore-water pH where CEC_tot applies [-]
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
f_CEC_0 = [f_Ca_0,f_Mg_0,f_K_0,f_Na_0,f_Al_0,f_H_0];
if any(~isfinite(f_CEC_0)) || any(f_CEC_0 < 0) || ...
        abs(sum(f_CEC_0)-1) > 1e-10
    error('ERW:invalidExchangeFractions', ...
        'Exchange fractions must be finite, nonnegative, and sum to one.');
end

%%% CO2 input check
if ~isscalar(CO2_air_extract) || ~isfinite(CO2_air_extract) || CO2_air_extract <= 0
    error('ERW:invalidParam', 'CO2_air_extract must be a finite positive scalar.');
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
if ~isfield(bulk_pH_params,'DIC_field_umol_L') || ...
        ~(isscalar(bulk_pH_params.DIC_field_umol_L) && ...
        isfinite(bulk_pH_params.DIC_field_umol_L) && ...
        bulk_pH_params.DIC_field_umol_L >= 0)
    error('ERW:missingParam', ...
        'DIC_field_umol_L must be a finite nonnegative scalar.');
end

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
DIC_d = bulk_pH_params.DIC_field_umol_L * dil;

%%% SCEPTER 1.0.2 obtains bulk pH from a separate 25 C laboratory-extraction
%%% simulation with transferred solid and exchanger inputs and the configured
%%% aqueous activity model.  This TC operator is a reduced equilibrium
%%% approximation of the same measurement, and it departs from SCEPTER's
%%% released default in two ways.
%%%
%%% The paste stands open to the laboratory air.  The sample is air-dried
%%% and the suspension is stirred in an open vessel (ISO 10390), so its
%%% dissolved CO2 sits at equilibrium with that air: the CO2 the field pore
%%% water carried has left, and a carbonate-bearing paste takes from the
%%% air what its bicarbonate needs.  SCEPTER closes its laboratory run to
%%% gas and carries the field DIC into it.
%%%
%%% The plot's solid carbonates dissolve.  Calcite and the hydrated Mg
%%% carbonate (nesquehonite) enter the paste with the soil they sit in and
%%% dissolve until the slurry is saturated or the solid is spent; SCEPTER
%%% freezes every solid.  Silicate solids stay inert in the paste, as in
%%% SCEPTER (the measurement is short).
%%%
%%% A soil holding free carbonate therefore reads the carbonate-air
%%% equilibrium of its paste, not the pore water it was sampled with.
T_lab = 298.15;
%%% Henry's constant of ERW_function (dimensionless, water over air) at its
%%% 25 C reference, on the air concentration the callers pass (ppm/22.41):
%%% the dissolved CO2 the model gives any water at equilibrium with that
%%% air at 25 C.
k_H_lab = 0.83;
CO2_lab = k_H_lab*CO2_air_extract;                          %% [umol/L]
pk1_lab = -(-356.309-0.0609*T_lab+21834.37/T_lab+ ...
    126.8339*log10(T_lab)-1684915/T_lab^2);
pk2_lab = -(-107.887-0.032528*T_lab+5151.79/T_lab+ ...
    38.92561*log10(T_lab)-563713.9/T_lab^2);
pk_w_lab = -(-283.971+13323/T_lab-0.0507*T_lab+ ...
    102.24447*log10(T_lab)-1119669/T_lab^2);
k1 = 10^(-pk1_lab)*conv_mol;
k2 = 10^(-pk2_lab)*conv_mol;
k_w = 10^(-pk_w_lab)*conv_mol^2;
k_nh3 = 10^(-(0.09018+2729.92/T_lab))*conv_mol;
%%% The slurry's own ionic strength conditions its constants
%%% (ERW_Conditional_Constants); the Porewater_to_Bulk_pH wrapper closes it.
%%% Exchange stays on the concentration basis, as in ERW_function.
I_extract = 0;
if isfield(bulk_pH_params,'extract_ionic_strength') && ...
        ~isempty(bulk_pH_params.extract_ionic_strength)
    I_extract = bulk_pH_params.extract_ionic_strength;
end
k_extract = ERW_Conditional_Constants( ...
    [k1,k2,k_w,K1_Al,K2_Al,K3_Al,K4_Al],[],I_extract);
k1 = k_extract(1); k2 = k_extract(2); k_w = k_extract(3);
K1_Al = k_extract(4); K2_Al = k_extract(5); K3_Al = k_extract(6); K4_Al = k_extract(7);

%%% The incoming field exchanger fixes each conserved cation inventory.  The
%%% final extract exchanger is evaluated at the trial extract pH, so changing
%%% site density repartitions that fixed inventory rather than changing it.
[dCEC_dpH,CEC_ref_pH] = ...
    ERW_Resolve_CEC_pH_Sensitivity(bulk_pH_params,exchange_scheme);
dCEC_areal = 0;
if dCEC_dpH > 0
    if ~isfield(bulk_pH_params,'CEC_cmolc_kg') || ...
            ~(isscalar(bulk_pH_params.CEC_cmolc_kg) && ...
            isfinite(bulk_pH_params.CEC_cmolc_kg) && ...
            bulk_pH_params.CEC_cmolc_kg > 0)
        error('ERW:invalidParam', ...
            'dCEC_dpH requires a positive finite CEC_cmolc_kg.');
    end
    dCEC_areal = dCEC_dpH*(CEC_tot/bulk_pH_params.CEC_cmolc_kg);
    CEC_at_H = @(Hx) ERW_CEC_At_H(Hx,conv_mol,CEC_tot,dCEC_areal,CEC_ref_pH);
else
    CEC_at_H = @(Hx) CEC_tot + zeros(size(Hx));
end
CEC_field = CEC_at_H(H_pw);
if ~(isfinite(CEC_field) && CEC_field > 0)
    error('ERW:nonpositiveCEC', ...
        'Variable CEC is nonpositive at the incoming pore-water pH.');
end
cec_inventory_factor = CEC_field / extract_water_L_m2;
cec_extract_factor = @(Hx) CEC_at_H(Hx) / extract_water_L_m2;

field_solid_CaCO3_umol_m2 = 0;
if isfield(bulk_pH_params,'solid_CaCO3_umol_m2')
    field_solid_CaCO3_umol_m2 = bulk_pH_params.solid_CaCO3_umol_m2;
end
if ~(isfinite(field_solid_CaCO3_umol_m2) && field_solid_CaCO3_umol_m2 >= 0)
    error('ERW:invalidParam', 'solid_CaCO3_umol_m2 must be finite and nonnegative.');
end
solid_CaCO3_umol_m2 = field_solid_CaCO3_umol_m2;
has_solid_calcite = solid_CaCO3_umol_m2 > 0;
solid_CaCO3_umol_L = solid_CaCO3_umol_m2 / extract_water_L_m2;
%%% Calcite solubility at the 25 C of the paste (Plummer and Busenberg 1982)
%%% as the concentration product the slurry's ionic strength allows:
%%% Ca*CO3 = K/gamma_div^2 at saturation, the activity basis of ERW_function.
K_CaCO3_lab = 10^(-171.9065-0.077993*T_lab+2839.319/T_lab+ ...
    71.595*log10(T_lab))*conv_mol^2;
K_CaCO3_extract = K_CaCO3_lab/ERW_Activity_Coefficient(2,I_extract)^2;
%%% The plot's solid Mg carbonate is the nesquehonite pool of ERW_function
%%% (MgCO3.3H2O).  It is far more soluble than calcite and equilibrates
%%% with water within hours (Harrison et al. 2013, 2019), so it enters the
%%% paste on the same terms: log K = -5.27 at the 25 C of the paste
%%% (Harrison et al. 2019), on the same activity basis.
field_solid_MgCO3_umol_m2 = 0;
if isfield(bulk_pH_params,'solid_MgCO3_umol_m2')
    field_solid_MgCO3_umol_m2 = bulk_pH_params.solid_MgCO3_umol_m2;
end
if ~(isfinite(field_solid_MgCO3_umol_m2) && field_solid_MgCO3_umol_m2 >= 0)
    error('ERW:invalidParam', 'solid_MgCO3_umol_m2 must be finite and nonnegative.');
end
solid_MgCO3_umol_m2 = field_solid_MgCO3_umol_m2;
has_solid_nesquehonite = solid_MgCO3_umol_m2 > 0;
solid_MgCO3_umol_L = solid_MgCO3_umol_m2 / extract_water_L_m2;
K_MgCO3_lab = 10^(-5.27)*conv_mol^2;
K_MgCO3_extract = K_MgCO3_lab/ERW_Activity_Coefficient(2,I_extract)^2;

%%% Charge-balance residual as function of final extract H [umol/L]
R_fn = @(H) charge_residual(H, Ca_d, Mg_d, K_d, Na_d, Al_d, ...
    NO3_d, PO4_d, NH4_d, An_d, CO2_lab, k1, k2, k_w, k_nh3, ...
    K1_Al, K2_Al, K3_Al, K4_Al, ...
    K_Ca_Mg, K_Ca_K, K_Ca_Na, K_Ca_Al, K_Ca_H, ...
    f_Ca_0, f_Mg_0, f_K_0, f_Na_0, f_Al_0, ...
    cec_inventory_factor, cec_extract_factor(H), conv_mol, ...
    scheme_is_appelo, appelo_alpha, K_Na_H, K_K_H, K_Ca_H2, K_Mg_H2, K_Al_H3);
if has_solid_calcite
    %%% The same slurry with xi umol/L of calcite dissolved into it: the Ca
    %%% inventory grows by xi, and the carbonate it brings joins a solution
    %%% whose CO2 the air holds.  The residual then lets xi settle at
    %%% saturation or at the whole solid for every trial H.
    R_at_xi = @(H,xi) charge_residual(H, Ca_d + xi, Mg_d, K_d, Na_d, Al_d, ...
        NO3_d, PO4_d, NH4_d, An_d, CO2_lab, k1, k2, k_w, k_nh3, ...
        K1_Al, K2_Al, K3_Al, K4_Al, ...
        K_Ca_Mg, K_Ca_K, K_Ca_Na, K_Ca_Al, K_Ca_H, ...
        f_Ca_0, f_Mg_0, f_K_0, f_Na_0, f_Al_0, ...
        cec_inventory_factor, cec_extract_factor(H), conv_mol, ...
        scheme_is_appelo, appelo_alpha, K_Na_H, K_K_H, K_Ca_H2, K_Mg_H2, K_Al_H3);
    R_fn = @(H) charge_residual_solid_calcite(H, R_at_xi, ...
        solid_CaCO3_umol_L, K_CaCO3_extract);
end
if has_solid_nesquehonite
    %%% With nesquehonite in the sample the slurry takes x_Ca of calcite and
    %%% x_Mg of nesquehonite: Ca grows by x_Ca and Mg by x_Mg.  The residual
    %%% settles the pair for every trial H.
    R_at_x = @(H,x_Ca,x_Mg) charge_residual(H, Ca_d + x_Ca, Mg_d + x_Mg, K_d, Na_d, Al_d, ...
        NO3_d, PO4_d, NH4_d, An_d, CO2_lab, k1, k2, k_w, k_nh3, ...
        K1_Al, K2_Al, K3_Al, K4_Al, ...
        K_Ca_Mg, K_Ca_K, K_Ca_Na, K_Ca_Al, K_Ca_H, ...
        f_Ca_0, f_Mg_0, f_K_0, f_Na_0, f_Al_0, ...
        cec_inventory_factor, cec_extract_factor(H), conv_mol, ...
        scheme_is_appelo, appelo_alpha, K_Na_H, K_K_H, K_Ca_H2, K_Mg_H2, K_Al_H3);
    R_fn = @(H) charge_residual_solid_carbonates(H, R_at_x, ...
        solid_CaCO3_umol_L, solid_MgCO3_umol_L, ...
        K_CaCO3_extract, K_MgCO3_extract);
end

%%% A previous converged extract root is the only continuity anchor.  A
%%% prescribed pore-water pH or observed bulk pH is not a previous extract
%%% state and cannot silently select among multiple physical roots.
has_previous_root = isfield(bulk_pH_params,'H_prev') && ...
    isscalar(bulk_pH_params.H_prev) && isfinite(bulk_pH_params.H_prev) && ...
    bulk_pH_params.H_prev > 0 && isfield(bulk_pH_params,'H_prev_valid') && ...
    isequal(bulk_pH_params.H_prev_valid,true);
has_inverse_anchor = isfield(bulk_pH_params,'H_extract_anchor_umolL') && ...
    isscalar(bulk_pH_params.H_extract_anchor_umolL) && ...
    isfinite(bulk_pH_params.H_extract_anchor_umolL) && ...
    bulk_pH_params.H_extract_anchor_umolL > 0;
if has_previous_root
    pH_anchor_scan = -log10(bulk_pH_params.H_prev / conv_mol);
elseif has_inverse_anchor
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
if has_previous_root || has_inverse_anchor
    narrow_pH_lo = max(full_pH_lo, pH_anchor_scan - 2);
    narrow_pH_hi = min(full_pH_hi, pH_anchor_scan + 2);
    narrow_n_points = round((narrow_pH_hi - narrow_pH_lo) / narrow_step_pH) + 1;
    if has_solid_nesquehonite
        %%% Each residual of a paste holding Mg carbonate settles two
        %%% dissolution extents, so the same grid is visited outward from
        %%% the anchor and only as far as the first sign change.
        [pH_scan, H_scan, R_scan] = bulk_pH_residual_scan_outward(R_fn, conv_mol, ...
            narrow_pH_lo, narrow_pH_hi, narrow_n_points, pH_anchor_scan);
    else
        [pH_scan, H_scan, R_scan] = bulk_pH_residual_scan(R_fn, conv_mol, ...
            narrow_pH_lo, narrow_pH_hi, narrow_n_points);
    end
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
if ~(has_previous_root || has_inverse_anchor) && numel(all_bracket_pH) ~= 1
    error('ERW:nonUniqueBulkRoot', ...
        'Cold bulk extraction requires one physical root; found %d.', ...
        numel(all_bracket_pH));
end
%%% A solid-bearing paste takes the same scan as a solid-free one.  Its
%%% residual rises monotonically with H like the solid-free residual (the
%%% carbonate dissolved at a trial H, calcite or nesquehonite, adds its
%%% cation's charge, and more of it dissolves at low pH, where the
%%% carbonate ion the air allows is scarce), so the bracket the scan
%%% selects is the root.
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
pH_bulk_raw = diag.pH;
diag.ionic_str_extract = ERW_Ionic_Strength(max(diag.Ca,0),max(diag.Mg,0), ...
    max(diag.Na,0),max(diag.K,0),max(diag.Al3,0),H,max(diag.NH4,0),NO3_d, ...
    diag.HCO3,diag.CO3,An_d,conv_mol,PO4_d,max(diag.AlOH,0), ...
    max(diag.AlOH2,0),max(diag.AlOH4,0),diag.OH);
pH_offset = 0;
if isfield(bulk_pH_params,'pH_offset')
    pH_offset = bulk_pH_params.pH_offset;
    if ~(isscalar(pH_offset) && isfinite(pH_offset))
        error('ERW:invalidParam', 'bulk_pH_params.pH_offset must be a finite scalar.');
    end
    if pH_offset ~= 0
        error('ERW:bulkPHOffsetNotPhysical', ...
            'A fixed pH offset cannot alter the physical bulk-extract operator.');
    end
end
pH_bulk = pH_bulk_raw;
if ~(isfinite(pH_bulk) && pH_bulk > 0 && pH_bulk < 14)
    error('ERW:invalidState', 'Calibrated bulk extract pH must satisfy 0 < pH < 14.');
end
if ~isfield(diag,'calcite_equilibrium_active')
    diag.calcite_equilibrium_active = false;
    diag.Omega_CaCO3 = NaN;
    diag.calcite_remaining_umol_L = 0;
    diag.calcite_Ca_residual_umol_L = 0;
end
diag.calcite_input_umol_m2 = solid_CaCO3_umol_m2;
diag.calcite_remaining_umol_m2 = diag.calcite_remaining_umol_L * extract_water_L_m2;
diag.calcite_Ca_residual_umol_m2 = diag.calcite_Ca_residual_umol_L * extract_water_L_m2;
if ~isfield(diag,'calcite_dissolved_umol_L')
    diag.calcite_dissolved_umol_L = 0;
end
diag.calcite_dissolved_umol_m2 = diag.calcite_dissolved_umol_L * extract_water_L_m2;
if ~isfield(diag,'nesquehonite_equilibrium_active')
    diag.nesquehonite_equilibrium_active = false;
    diag.Omega_MgCO3 = NaN;
    diag.nesquehonite_remaining_umol_L = 0;
    diag.nesquehonite_dissolved_umol_L = 0;
end
diag.nesquehonite_input_umol_m2 = solid_MgCO3_umol_m2;
diag.nesquehonite_remaining_umol_m2 = diag.nesquehonite_remaining_umol_L * extract_water_L_m2;
diag.nesquehonite_dissolved_umol_m2 = diag.nesquehonite_dissolved_umol_L * extract_water_L_m2;
full_bracket_pH = [];
%%% A paste that dissolved carbonate, calcite or Mg carbonate, stands far
%%% above an acid pore water by chemistry, and its residual has one root,
%%% so there is no pore-proximal root to look for.
large_offset = abs(pH_bulk_raw - pH_pw) > 2.5 && ...
    ~(has_solid_nesquehonite && diag.nesquehonite_dissolved_umol_L > 0) && ...
    ~(has_solid_calcite && diag.calcite_dissolved_umol_L > 0);
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
        pH_bulk_raw, pH_pw, scan_used, anchor_pH_used, numel(all_bracket_pH), ...
        0.5*(-log10(H/conv_mol) + pH_bulk_raw), near_pore);
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
    diag.CO2_ext = diag.CO2;
    diag.DIC_ext = diag.DIC;
    %%% Carbon the paste took from the air [umol per litre of extract
    %%% water]: its DIC less the field DIC carried in and the carbonate its
    %%% solids dissolved.  Negative when the pore water's CO2 degasses.
    diag.DIC_transferred_umol_L = DIC_d;
    diag.CO2_air_exchange_umol_L = diag.DIC - DIC_d - ...
        diag.calcite_dissolved_umol_L - diag.nesquehonite_dissolved_umol_L;
    diag.pH_bulk_raw = pH_bulk_raw;
    diag.pH_offset = pH_offset;
    diag.pH_bulk = pH_bulk;
    diag.pH_pore = pH_pw;

    diag.K_CEC = bulk_pH_params.K_CEC;
    diag.f_CEC_initial = [f_Ca_0, f_Mg_0, f_K_0, f_Na_0, f_Al_0, f_H_0];
    diag.dilution_factor = dil;
    diag.cec_factor = cec_extract_factor(H);
    diag.cec_inventory_factor = cec_inventory_factor;
    diag.CEC_field = CEC_field;
    diag.CEC_extract = cec_extract_factor(H)*extract_water_L_m2;
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
function [pH_scan, H_scan, R_scan] = bulk_pH_residual_scan_outward(R_fn, conv_mol, pH_lo, pH_hi, n_points, pH_anchor)
%%% The grid of bulk_pH_residual_scan, evaluated outward from the anchor and
%%% stopped at the first sign change between neighbouring points.  The
%%% residual of a solid-bearing paste rises with H, so that sign change is
%%% the only one and the bracket is the one the whole scan would select.
%%% Points never reached stay NaN, which bulk_pH_select_bracket skips.
pH_scan = linspace(pH_lo, pH_hi, n_points);
H_scan = 10.^(-pH_scan) * conv_mol;
R_scan = nan(size(H_scan));
[~, i_anchor] = min(abs(pH_scan - pH_anchor));
for step = 0:n_points-1
    if step == 0
        visit = i_anchor;
    else
        visit = [i_anchor + step, i_anchor - step];
    end
    for ii = visit(visit >= 1 & visit <= n_points)
        try
            R_scan(ii) = R_fn(H_scan(ii));
        catch ME
            if startsWith(ME.identifier, 'ERW:')
                R_scan(ii) = NaN;
            else
                rethrow(ME)
            end
        end
        for jj = [ii-1, ii+1]
            if jj >= 1 && jj <= n_points && ...
                    isfinite(R_scan(ii)) && isfinite(R_scan(jj)) && ...
                    (R_scan(ii) == 0 || R_scan(jj) == 0 || ...
                    sign(R_scan(ii)) ~= sign(R_scan(jj)))
                return
            end
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
abs_tol  = 1e-10;
rel_tol  = 1e-10;
rel_charge_tol = 1e-9;
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

        if abs(R) < abs_tol || abs(R) < charge_scale * rel_charge_tol
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

if ~converged
    %% [H_left,H_right] still brackets the root (opposite-sign residuals),
    %% so bisection converges without aborting the run.
    [H, R, bisect_iter, bisect_converged] = bisect_bracket_H(R_fn, H_left, H_right, ...
        abs_tol, rel_tol, rel_charge_tol);
    iter = iter + bisect_iter;
    if ~bisect_converged
        error('ERW:bulkPHSolverNoConvergence', ...
            'Bulk extract charge balance did not converge inside its physical bracket.');
    end
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
function [H, R, iter, converged] = bisect_bracket_H(R_fn, H_left, H_right, abs_tol, rel_tol, rel_charge_tol)
%%% Bisection on a sign-bracketed charge residual; unconditionally
%%% convergent.  Used as the no-abort fallback when Newton-Raphson stalls
%%% or wanders non-finite.  Uses the same tolerances as the Newton solve.
[R_left, ~]  = R_fn(H_left);
[R_right, ~] = R_fn(H_right);
converged = false;
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
        converged = abs(R) < abs_tol || abs(R) < charge_scale*rel_charge_tol;
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
    f_Ca_0, f_Mg_0, f_K_0, f_Na_0, f_Al_0, ...
    cec_inventory_factor, cec_factor, conv_mol, ...
    scheme_is_appelo, appelo_alpha, K_Na_H, K_K_H, K_Ca_H2, K_Mg_H2, K_Al_H3)

%%% Open-extract speciation: the laboratory air holds the dissolved CO2 at
%%% CO2_ext, and bicarbonate and carbonate follow it at the trial H.
HCO3 = k1 * CO2_ext / H;
CO3  = k2 * HCO3 / H;
OH   = k_w / H;

%%% Exchange equilibrium and mass conservation on final extract solution.
state = exchange_state(H, Ca_d, Mg_d, K_d, Na_d, NH4_d, Al_d, ...
    K1_Al, K2_Al, K3_Al, K4_Al, ...
    K_Ca_Mg, K_Ca_K, K_Ca_Na, K_Ca_Al, K_Ca_H, ...
    f_Ca_0, f_Mg_0, f_K_0, f_Na_0, f_Al_0, ...
    cec_inventory_factor, cec_factor, ...
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
state.CO2 = CO2_ext;
state.DIC = CO2_ext+HCO3+CO3;
state.OH = OH;
state.pH = -log10(H / conv_mol);

end

%% =======================================================================
function [R, state] = charge_residual_solid_calcite(H, R_at_xi, S, K_sat)
%%% Charge residual of the open slurry at extract H when the soil sample
%%% carries S umol of solid calcite per litre of extract water.  Calcite
%%% dissolves until the slurry is saturated (Ca*CO3 = K_sat, the
%%% concentration product at the slurry's ionic strength) or the solid is
%%% spent.  A slurry already supersaturated before anything dissolves is
%%% left as it is: nothing precipitates within a pH measurement.  At a
%%% given H the air fixes the carbonate ion and the saturation ratio rises
%%% monotonically with the amount dissolved, because the Ca inventory
%%% grows, so the extent is bracketed on [0,S] and bisected.  The extent
%%% returned is the
%%% undersaturated end of the last bracket, so the remaining solid is never
%%% negative.
[R, state] = R_at_xi(H, S);
xi = S;
active = false;
if state.Ca*state.CO3/K_sat > 1
    active = true;
    %%% At a trial H well above the pore water's, the Gaines-Thomas exchanger
    %%% may have no partition until some calcite has dissolved: the cations
    %%% the slurry holds cannot fill its sites.  Such a slurry is short of
    %%% Ca, so it counts as undersaturated and the bracket moves up.  The
    %%% Appelo exchanger has a partition for any inventory and never takes
    %%% this path.
    [R, state, hosted] = paste_state_if_hosted(R_at_xi, H, 0);
    xi = 0;
    if ~hosted || state.Ca*state.CO3/K_sat < 1
        lo = 0;
        hi = S;
        lo_hosted = hosted;
        for it = 1:200
            mid = 0.5*(lo+hi);
            [~, state_mid, mid_hosted] = paste_state_if_hosted(R_at_xi, H, mid);
            if mid_hosted && state_mid.Ca*state_mid.CO3/K_sat > 1
                hi = mid;
            else
                lo = mid;
                lo_hosted = mid_hosted;
            end
            if hi-lo <= 1e-13*S
                break
            end
        end
        xi = lo;
        if ~lo_hosted
            xi = hi;
        end
        [R, state] = R_at_xi(H, xi);
    end
end
state.calcite_equilibrium_active = active;
state.calcite_dissolved_umol_L = xi;
state.calcite_remaining_umol_L = S - xi;
state.Omega_CaCO3 = state.Ca*state.CO3/K_sat;
state.calcite_Ca_residual_umol_L = 0;
end

%% =======================================================================
function [R, state, hosted] = paste_state_if_hosted(R_at_xi, H, xi)
%%% The slurry at extract H with xi of calcite dissolved, or hosted = false
%%% when its exchanger has no partition there (ERW:partitionExchange).
hosted = true;
try
    [R, state] = R_at_xi(H, xi);
catch ME
    if ~strcmp(ME.identifier, 'ERW:partitionExchange')
        rethrow(ME)
    end
    hosted = false;
    R = NaN;
    state = struct('Ca', NaN, 'CO3', NaN);
end
end

%% =======================================================================
function [R, state] = charge_residual_solid_carbonates(H, R_at_x, S_Ca, S_Mg, K_Ca, K_Mg)
%%% Charge residual of the open slurry at extract H when the soil sample
%%% carries S_Ca umol of solid calcite and S_Mg umol of solid nesquehonite
%%% per litre of extract water.  Each solid dissolves until the slurry is
%%% saturated with it or the solid is spent, and nothing precipitates
%%% within a pH measurement.  The two meet on the exchanger, so the
%%% extents are settled together: for each solid
%%%     x = S and Omega <= 1,   0 < x < S and Omega = 1,   or   x = 0 and Omega >= 1.
%%% At a given H the air fixes the carbonate ion, and either saturation
%%% ratio rises with either extent (the cation inventory grows, and the
%%% added cation displaces the other from the exchanger).  The calcite
%%% extent that meets its own condition therefore falls as more
%%% nesquehonite dissolves, and the nesquehonite ratio taken along that
%%% calcite extent still rises, so the pair is one bracketed solve for the
%%% nesquehonite extent with the calcite extent resolved inside it.
[x_Ca, R, state] = paste_calcite_extent(H, R_at_x, S_Ca, S_Mg, K_Ca);
x_Mg = S_Mg;
log_hi = log(state.Mg*state.CO3/K_Mg);
if log_hi > 0
    [x_Ca, R, state] = paste_calcite_extent(H, R_at_x, S_Ca, 0, K_Ca);
    x_Mg = 0;
    log_lo = log(state.Mg*state.CO3/K_Mg);
    if log_lo < 0
        log_sat = @(x) paste_nesquehonite_log_saturation(H, R_at_x, S_Ca, x, K_Ca, K_Mg);
        x_Mg = paste_extent_root(log_sat, S_Mg, log_lo, log_hi);
        [x_Ca, R, state] = paste_calcite_extent(H, R_at_x, S_Ca, x_Mg, K_Ca);
    end
end
state.calcite_equilibrium_active = x_Ca < S_Ca;
state.calcite_dissolved_umol_L = x_Ca;
state.calcite_remaining_umol_L = S_Ca - x_Ca;
state.Omega_CaCO3 = state.Ca*state.CO3/K_Ca;
state.calcite_Ca_residual_umol_L = 0;
state.nesquehonite_equilibrium_active = x_Mg < S_Mg;
state.nesquehonite_dissolved_umol_L = x_Mg;
state.nesquehonite_remaining_umol_L = S_Mg - x_Mg;
state.Omega_MgCO3 = state.Mg*state.CO3/K_Mg;
end

%% =======================================================================
function [x_Ca, R, state] = paste_calcite_extent(H, R_at_x, S_Ca, x_Mg, K_Ca)
%%% Calcite dissolved in the slurry at extract H once x_Mg of nesquehonite
%%% is in it: the whole solid, none of it, or the amount that saturates.
[R, state] = R_at_x(H, S_Ca, x_Mg);
x_Ca = S_Ca;
if S_Ca > 0
    log_hi = log(state.Ca*state.CO3/K_Ca);
    if log_hi > 0
        [R, state] = R_at_x(H, 0, x_Mg);
        x_Ca = 0;
        log_lo = log(state.Ca*state.CO3/K_Ca);
        if log_lo < 0
            log_sat = @(x) paste_calcite_log_saturation(H, R_at_x, x, x_Mg, K_Ca);
            x_Ca = paste_extent_root(log_sat, S_Ca, log_lo, log_hi);
            [R, state] = R_at_x(H, x_Ca, x_Mg);
        end
    end
end
end

%% =======================================================================
function v = paste_calcite_log_saturation(H, R_at_x, x_Ca, x_Mg, K_Ca)
[~, state] = R_at_x(H, x_Ca, x_Mg);
v = log(state.Ca*state.CO3/K_Ca);
end

%% =======================================================================
function v = paste_nesquehonite_log_saturation(H, R_at_x, S_Ca, x_Mg, K_Ca, K_Mg)
[~, ~, state] = paste_calcite_extent(H, R_at_x, S_Ca, x_Mg, K_Ca);
v = log(state.Mg*state.CO3/K_Mg);
end

%% =======================================================================
function x = paste_extent_root(log_sat, S, f_lo, f_hi)
%%% Extent on (0,S) at which a rising log saturation ratio crosses zero,
%%% from its values f_lo < 0 at 0 and f_hi > 0 at S.  Brent's method: the
%%% root stays bracketed throughout, inverse quadratic or secant steps are
%%% taken while they shrink the bracket fast enough and bisection steps
%%% otherwise, and the bracket closes to 1e-13 of S.  The undersaturated
%%% end of the final bracket is returned, so the solid remaining is never
%%% negative and dissolution never leaves the slurry supersaturated.
tol = 1e-13*S;
a = 0;
f_a = f_lo;
b = S;
f_b = f_hi;
c = b;
f_c = f_b;
d = b - a;
e = d;
for it = 1:200
    if (f_b > 0) == (f_c > 0)
        c = a;
        f_c = f_a;
        d = b - a;
        e = d;
    end
    if abs(f_c) < abs(f_b)
        a = b;
        b = c;
        c = a;
        f_a = f_b;
        f_b = f_c;
        f_c = f_a;
    end
    half = 0.5*(c - b);
    if abs(half) <= tol || f_b == 0
        break
    end
    if abs(e) >= tol && abs(f_a) > abs(f_b)
        s = f_b/f_a;
        if a == c
            p = 2*half*s;
            q = 1 - s;
        else
            q = f_a/f_c;
            r = f_b/f_c;
            p = s*(2*half*q*(q - r) - (b - a)*(r - 1));
            q = (q - 1)*(r - 1)*(s - 1);
        end
        if p > 0
            q = -q;
        end
        p = abs(p);
        if 2*p < min(3*half*q - abs(tol*q), abs(e*q))
            e = d;
            d = p/q;
        else
            d = half;
            e = d;
        end
    else
        d = half;
        e = d;
    end
    a = b;
    f_a = f_b;
    if abs(d) > tol
        b = b + d;
    else
        b = b + sign(half)*tol;
    end
    f_b = log_sat(b);
end
if (f_b > 0) == (f_c > 0)
    c = a;
    f_c = f_a;
end
if f_b <= 0
    x = b;
else
    x = c;
end
end


%% =======================================================================
function [R, state] = charge_residual_calcite(H, Ca_d, Mg_d, K_d, Na_d, Al_d, ...
    NO3_d, PO4_d, NH4_d, An_d, CO2_ext, k1, k2, k_w, ...
    K1_Al, K2_Al, K3_Al, K4_Al, ...
    K_Ca_Mg, K_Ca_K, K_Ca_Na, K_Ca_Al, K_Ca_H, ...
    f_Ca_0, f_Mg_0, f_K_0, f_Na_0, f_Al_0, ...
    cec_inventory_factor, cec_factor, conv_mol, ...
    scheme_is_appelo, appelo_alpha, K_Na_H, K_K_H, K_Ca_H2, K_Mg_H2, K_Al_H3, ...
    solid_CaCO3_umol_L, K_CaCO3)

HCO3 = k1 * CO2_ext / H;
CO3  = k2 * HCO3 / H;
OH   = k_w / H;
Ca_total_all = Ca_d + cec_inventory_factor*f_Ca_0/2 + solid_CaCO3_umol_L;
Ca_sat = K_CaCO3 / CO3;

calcite_equilibrium_active = false;
calcite_remaining = 0;
if isfinite(Ca_sat) && Ca_sat > 0 && Ca_sat < Ca_total_all
    state_sat = exchange_state_fixed_ca(H, Ca_sat, Mg_d, K_d, Na_d, NH4_d, Al_d, ...
        K1_Al, K2_Al, K3_Al, K4_Al, ...
        K_Ca_Mg, K_Ca_K, K_Ca_Na, K_Ca_Al, K_Ca_H, ...
        f_Mg_0, f_K_0, f_Na_0, f_Al_0, ...
        cec_inventory_factor, cec_factor, ...
        scheme_is_appelo, appelo_alpha, K_Na_H, K_K_H, K_Ca_H2, K_Mg_H2, K_Al_H3);
    Ca_used = state_sat.Ca + cec_factor*state_sat.f_Ca/2;
    Ca_tol = 1e-10*max(Ca_total_all,1);
    if Ca_used <= Ca_total_all + Ca_tol
        state = state_sat;
        calcite_remaining = max(Ca_total_all - Ca_used,0);
        calcite_equilibrium_active = true;
    end
end
if ~calcite_equilibrium_active
    state = exchange_state(H, Ca_d + solid_CaCO3_umol_L, Mg_d, K_d, Na_d, NH4_d, Al_d, ...
        K1_Al, K2_Al, K3_Al, K4_Al, ...
        K_Ca_Mg, K_Ca_K, K_Ca_Na, K_Ca_Al, K_Ca_H, ...
        f_Ca_0, f_Mg_0, f_K_0, f_Na_0, f_Al_0, ...
        cec_inventory_factor, cec_factor, ...
        scheme_is_appelo, appelo_alpha, K_Na_H, K_K_H, K_Ca_H2, K_Mg_H2, K_Al_H3);
end

cations = H + 2*state.Ca + 2*state.Mg + state.K + state.Na + state.NH4 ...
    + 3*state.Al3 + 2*state.AlOH + state.AlOH2;
anions = OH + HCO3 + 2*CO3 + NO3_d + PO4_d + An_d + state.AlOH4;
R = cations - anions;
state.charge_scale = max(abs(cations) + abs(anions), eps);
state.HCO3 = HCO3;
state.CO3 = CO3;
state.CO2 = CO2_ext;
state.DIC = CO2_ext+HCO3+CO3;
state.OH = OH;
state.pH = -log10(H / conv_mol);
state.calcite_equilibrium_active = calcite_equilibrium_active;
state.calcite_remaining_umol_L = calcite_remaining;
state.Omega_CaCO3 = state.Ca*CO3/K_CaCO3;
state.calcite_Ca_residual_umol_L = ...
    state.Ca + cec_factor*state.f_Ca/2 + calcite_remaining - Ca_total_all;
end

%% =======================================================================
function state = exchange_state_fixed_ca(H, Ca, Mg_d, K_d, Na_d, NH4_d, Al_d, ...
    K1_Al, K2_Al, K3_Al, K4_Al, ...
    K_Ca_Mg, K_Ca_K, K_Ca_Na, K_Ca_Al, K_Ca_H, ...
    f_Mg_0, f_K_0, f_Na_0, f_Al_0, cec_inventory_factor, cec_factor, ...
    scheme_is_appelo, appelo_alpha, K_Na_H, K_K_H, K_Ca_H2, K_Mg_H2, K_Al_H3)

Mg_total = Mg_d + cec_inventory_factor*f_Mg_0/2;
K_total  = K_d  + cec_inventory_factor*f_K_0;
Na_total = Na_d + cec_inventory_factor*f_Na_0;
Al_total = Al_d + cec_inventory_factor*f_Al_0/3;
Al3_frac = al3_fraction(H,K1_Al,K2_Al,K3_Al,K4_Al);

if scheme_is_appelo
    residual_fn = @(f_H) fixed_ca_appelo_residual(f_H,H,Ca,Mg_total,K_total,Na_total,Al_total, ...
        cec_factor,appelo_alpha,K_Na_H,K_K_H,K_Ca_H2,K_Mg_H2,K_Al_H3,Al3_frac);
else
    residual_fn = @(f_Ca) fixed_ca_gt_residual(f_Ca,H,Ca,Mg_total,K_total,Na_total,Al_total, ...
        cec_factor,K_Ca_Mg,K_Ca_K,K_Ca_Na,K_Ca_Al,K_Ca_H,Al3_frac);
end
f_lo = 1e-24;
f_hi = 1 - 1e-14;
g_lo = residual_fn(f_lo);
g_hi = residual_fn(f_hi);
if ~(isfinite(g_lo) && isfinite(g_hi) && sign(g_lo) ~= sign(g_hi))
    error('ERW:partitionExchange', ...
        'No physical exchange partition is compatible with calcite saturation.');
end
for k = 1:300
    f_mid = 0.5*(f_lo + f_hi);
    g_mid = residual_fn(f_mid);
    if abs(g_mid) < 1e-12
        break
    end
    if sign(g_mid) == sign(g_lo)
        f_lo = f_mid;
        g_lo = g_mid;
    else
        f_hi = f_mid;
    end
end
f_root = 0.5*(f_lo + f_hi);
if scheme_is_appelo
    [~,f_Na,f_K,f_Ca,f_Mg,f_Al] = residual_fn(f_root);
    f_H = f_root;
else
    [~,f_Mg,f_K,f_Na,f_Al,f_H] = residual_fn(f_root);
    f_Ca = f_root;
end

Mg = Mg_total - cec_factor*f_Mg/2;
K  = K_total  - cec_factor*f_K;
Na = Na_total - cec_factor*f_Na;
Al = Al_total - cec_factor*f_Al/3;
[Al3,AlOH,AlOH2,AlOH3,AlOH4] = al_hydrolysis_species(Al,H,K1_Al,K2_Al,K3_Al,K4_Al);
fractions = [f_Ca,f_Mg,f_K,f_Na,f_Al,f_H];
concentrations = [Ca,Mg,K,Na,Al];
if any(~isfinite([fractions,concentrations])) || any(fractions < -1e-10) || ...
        abs(sum(fractions)-1) > 1e-8 || any(concentrations < -1e-8)
    error('ERW:partitionExchange', 'Calcite-buffered exchange partition is unphysical.');
end

state.f_Ca=f_Ca; state.f_Mg=f_Mg; state.f_K=f_K; state.f_Na=f_Na;
state.f_Al=f_Al; state.f_H=f_H; state.site_sum=sum(fractions);
state.Ca=Ca; state.Mg=max(Mg,0); state.K=max(K,0); state.Na=max(Na,0);
state.NH4=NH4_d; state.Al=max(Al,0); state.Al3=Al3; state.AlOH=AlOH;
state.AlOH2=AlOH2; state.AlOH3=AlOH3; state.AlOH4=AlOH4;
state.exchange_mass_residuals_umol_L = [0, ...
    state.Mg + cec_factor*f_Mg/2 - Mg_total, ...
    state.K + cec_factor*f_K - K_total, ...
    state.Na + cec_factor*f_Na - Na_total, ...
    state.Al + cec_factor*f_Al/3 - Al_total];
if scheme_is_appelo
    gamma = 10^(appelo_alpha*f_H);
    state.exchange_residuals = [ ...
        f_Na-K_Na_H*state.Na*(f_H/H)*gamma, f_K-K_K_H*state.K*(f_H/H)*gamma, ...
        f_Ca-K_Ca_H2*Ca*(f_H/H)^2*gamma^2, f_Mg-K_Mg_H2*state.Mg*(f_H/H)^2*gamma^2, ...
        f_Al-K_Al_H3*Al3*(f_H/H)^3*gamma^3];
else
    state.exchange_residuals = [ ...
        f_Mg-(state.Mg/K_Ca_Mg)*f_Ca/Ca, ...
        f_Na-state.Na*sqrt(f_Ca/(K_Ca_Na*Ca)), ...
        f_K-state.K*sqrt(f_Ca/(K_Ca_K*Ca)), ...
        f_Al-Al3*sqrt(f_Ca^3/(K_Ca_Al*Ca^3)), ...
        f_H-H*sqrt(f_Ca/(K_Ca_H*Ca))];
end
end

%% =======================================================================
function [residual,f_Mg,f_K,f_Na,f_Al,f_H] = fixed_ca_gt_residual( ...
    f_Ca,H,Ca,Mg_total,K_total,Na_total,Al_total,cec_factor, ...
    K_Ca_Mg,K_Ca_K,K_Ca_Na,K_Ca_Al,K_Ca_H,Al3_frac)
A_Mg = f_Ca/(K_Ca_Mg*Ca);
A_K = sqrt(f_Ca/(K_Ca_K*Ca));
A_Na = sqrt(f_Ca/(K_Ca_Na*Ca));
A_Al = Al3_frac*sqrt(f_Ca^3/(K_Ca_Al*Ca^3));
f_Mg = A_Mg*Mg_total/(1+A_Mg*cec_factor/2);
f_K = A_K*K_total/(1+A_K*cec_factor);
f_Na = A_Na*Na_total/(1+A_Na*cec_factor);
f_Al = A_Al*Al_total/(1+A_Al*cec_factor/3);
f_H = H*sqrt(f_Ca/(K_Ca_H*Ca));
residual = f_Ca+f_Mg+f_K+f_Na+f_Al+f_H-1;
end

%% =======================================================================
function [residual,f_Na,f_K,f_Ca,f_Mg,f_Al] = fixed_ca_appelo_residual( ...
    f_H,H,Ca,Mg_total,K_total,Na_total,Al_total,cec_factor,appelo_alpha, ...
    K_Na_H,K_K_H,K_Ca_H2,K_Mg_H2,K_Al_H3,Al3_frac)
gamma = 10^(appelo_alpha*f_H);
rho1 = (f_H/H)*gamma;
rho2 = rho1^2;
rho3 = rho1^3;
f_Na = K_Na_H*rho1*Na_total/(1+K_Na_H*rho1*cec_factor);
f_K = K_K_H*rho1*K_total/(1+K_K_H*rho1*cec_factor);
f_Ca = K_Ca_H2*Ca*rho2;
f_Mg = K_Mg_H2*rho2*Mg_total/(1+K_Mg_H2*rho2*cec_factor/2);
f_Al = K_Al_H3*Al3_frac*rho3*Al_total/(1+K_Al_H3*Al3_frac*rho3*cec_factor/3);
residual = f_H+f_Na+f_K+f_Ca+f_Mg+f_Al-1;
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
    f_Ca_0, f_Mg_0, f_K_0, f_Na_0, f_Al_0, ...
    cec_inventory_factor, cec_factor, ...
    scheme_is_appelo, appelo_alpha, K_Na_H, K_K_H, K_Ca_H2, K_Mg_H2, K_Al_H3)

[f_Ca, f_Mg, f_K, f_Na, f_Al, f_H, Ca_eq, Mg_eq, K_eq, Na_eq, Al_eq, mass_residuals] = ...
    partition_exchange_inventory(H, Ca_d, Mg_d, K_d, Na_d, Al_d, ...
    K1_Al, K2_Al, K3_Al, K4_Al, ...
    K_Ca_Mg, K_Ca_K, K_Ca_Na, K_Ca_Al, K_Ca_H, ...
    f_Ca_0, f_Mg_0, f_K_0, f_Na_0, f_Al_0, ...
    cec_inventory_factor, cec_factor, ...
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
    f_Ca_0, f_Mg_0, f_K_0, f_Na_0, f_Al_0, ...
    cec_inventory_factor, cec_factor, ...
    scheme_is_appelo, appelo_alpha, K_Na_H, K_K_H, K_Ca_H2, K_Mg_H2, K_Al_H3)
if scheme_is_appelo
    [f_Ca, f_Mg, f_K, f_Na, f_Al, f_H, Ca, Mg, K, Na, Al, mass_residuals] = ...
        partition_exchange_inventory_appelo(H, Ca_d, Mg_d, K_d, Na_d, Al_d, ...
        K1_Al, K2_Al, K3_Al, K4_Al, ...
        f_Ca_0, f_Mg_0, f_K_0, f_Na_0, f_Al_0, ...
        cec_inventory_factor, cec_factor, ...
        appelo_alpha, K_Na_H, K_K_H, K_Ca_H2, K_Mg_H2, K_Al_H3);
    return
end
%%% Conserved-inventory Gaines-Thomas partition of the air-dry extract.
%%%
%%% The initial dissolved concentration and the exchangeable pool of each
%%% cation are one conserved inventory.  Exchange fractions are charge
%%% fractions, so divalent Ca/Mg are divided by 2 and trivalent Al by 3
%%% before adding the incoming exchangeable charge
%%% cec_inventory_factor*f_X_0:
%%%     X_total = X_d + cec_inventory_factor * f_X_0 / z_X   [umol/L]
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
Ca_total = Ca_d + cec_inventory_factor*f_Ca_0/2;
Mg_total = Mg_d + cec_inventory_factor*f_Mg_0/2;
K_total  = K_d  + cec_inventory_factor*f_K_0;
Na_total = Na_d + cec_inventory_factor*f_Na_0;
Al_total = Al_d + cec_inventory_factor*f_Al_0/3;
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
    f_Ca_0, f_Mg_0, f_K_0, f_Na_0, f_Al_0, ...
    cec_inventory_factor, cec_factor, ...
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
Ca_total = Ca_d + cec_inventory_factor*f_Ca_0/2;
Mg_total = Mg_d + cec_inventory_factor*f_Mg_0/2;
K_total  = K_d  + cec_inventory_factor*f_K_0;
Na_total = Na_d + cec_inventory_factor*f_Na_0;
Al_total = Al_d + cec_inventory_factor*f_Al_0/3;
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
