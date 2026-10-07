%#########################################################################
%   Function EXCHANGE_FRACTIONS_GT                                       #
%#########################################################################
%%% Cation-exchange fractions from the Gaines-Thomas site balance, summing
%%% to 1.  These are the SAME relations the daily ERW
%%% solver enforces (ERW_function.m equations p(11)-p(16),
%%% lines 320-324) so any state built with this helper is a fixed point of
%%% the daily exchange equilibrium.  The soil-class selectivity constants are
%%% used as effective Gaines-Thomas constants throughout; no separate
%%% surface-potential coefficient is applied.
%%%
%%% Site balance (Ca taken as the reference cation):
%%%   f_Mg = Mg * f_Ca / (K_Ca_Mg*Ca)
%%%   f_Na = Na * sqrt(f_Ca/(K_Ca_Na*Ca))
%%%   f_K  = K  * sqrt(f_Ca/(K_Ca_K *Ca))
%%%   f_H  = H  * sqrt(f_Ca/(K_Ca_H *Ca))
%%%   f_Al = (Al3/Al_CEC_scale) * sqrt(f_Ca^3/(K_Ca_Al*Ca^3))
%%%   f_Ca + f_Mg + f_K + f_Na + f_Al + f_H = 1
%%% with Al3 = (H^4/den_H)*Al_w the free Al3+ from the hydrolysis chain.
%%%
%%% INPUTS  (concentrations in umol/L; Al_w in ERW model units = Al_CEC_scale x actual)
%%%   H, Ca, Mg, K, Na, Al_w
%%%   K_CEC = [K_Ca_Mg, K_Ca_K, K_Ca_Na, K_Ca_Al, K_Ca_H]
%%%   K1_Al..K4_Al   Al hydrolysis constants [umol/L]
%%%   Al_CEC_scale   ERW solver Al storage scaling (= conv_Al)
%%% OUTPUT
%%%   f_Ca, f_Mg, f_K, f_Na, f_Al, f_H   exchange fractions [-], summing to 1
function [f_Ca, f_Mg, f_K, f_Na, f_Al, f_H] = Exchange_Fractions_GT( ...
    H, Ca, Mg, K, Na, Al_w, K_CEC, K1_Al, K2_Al, K3_Al, K4_Al, Al_CEC_scale, ...
    scheme, appelo_params)


%%% Optional SCEPTER/Appelo (1994) H-referenced exchange.
if nargin < 13 || isempty(scheme)
    scheme = 'TC_GT';
end
if nargin < 14
    appelo_params = [];
end
if isstring(scheme) && isscalar(scheme)
    scheme = char(scheme);
end
if ~ischar(scheme) || ~any(strcmpi(scheme, {'TC_GT','appelo'}))
    error('ERW:invalidParam', 'Exchange scheme must be TC_GT or appelo.');
end
scheme_is_appelo = strcmpi(scheme, 'appelo');
if scheme_is_appelo
    [K_Na_H, K_K_H, K_Ca_H2, K_Mg_H2, K_Al_H3, appelo_alpha] = ...
        Resolve_Appelo_Constants(appelo_params);
    [f_Ca, f_Mg, f_K, f_Na, f_Al, f_H] = ...
        Appelo_Exchange_Fractions(H, Ca, Mg, K, Na, Al_w, Al_CEC_scale, ...
        K1_Al, K2_Al, K3_Al, K4_Al, ...
        K_Na_H, K_K_H, K_Ca_H2, K_Mg_H2, K_Al_H3, appelo_alpha);
    return
end
if ~(isfinite(H) && H > 0) || ~(isfinite(Ca) && Ca > 0)
    error('ERW:invalidState', 'Exchange_Fractions_GT requires positive H and Ca.');
end
if ~(isfinite(Al_CEC_scale) && Al_CEC_scale > 0)
    error('ERW:invalidParam', 'Exchange_Fractions_GT requires positive Al_CEC_scale.');
end

K_Ca_Mg = K_CEC(1); K_Ca_K = K_CEC(2); K_Ca_Na = K_CEC(3);
K_Ca_Al = K_CEC(4); K_Ca_H  = K_CEC(5);

%%% Free Al3+ from the hydrolysis chain (matches ERW_function p(6))
den_H = H^4 + H^3*K1_Al + H^2*K1_Al*K2_Al + H*K1_Al*K2_Al*K3_Al + K1_Al*K2_Al*K3_Al*K4_Al;
if ~(isfinite(den_H) && den_H > 0)
    error('ERW:invalidState', 'Exchange_Fractions_GT: invalid Al hydrolysis denominator.');
end
if ~(isfinite(Al_w) && Al_w >= 0)
    error('ERW:invalidState', 'Exchange_Fractions_GT requires nonnegative dissolved Al.');
end
Al3 = (H^4/den_H) * Al_w;

%%% Gaines-Thomas coefficients (f_X / f_Ca^power) evaluated at the given solution
A_Mg  = Mg  / (K_Ca_Mg  * Ca);
A_Na  = Na  / sqrt(K_Ca_Na * Ca);
A_K   = K   / sqrt(K_Ca_K  * Ca);
A_H   = H   / sqrt(K_Ca_H  * Ca);
A_Al  = (Al3/Al_CEC_scale) / sqrt(K_Ca_Al * Ca^3);

%%% Solve in y=sqrt(f_Ca).  The residual is continuous and strictly increasing
%%% on [0,1], while y remains numerically resolvable when f_Ca tends to zero.
site_sum = @(y) A_Al*y^3 + (1 + A_Mg)*y^2 ...
    + (A_K + A_Na + A_H)*y - 1;
lo = 0;
hi = 1;
for iter = 1:80
    mid = 0.5*(lo + hi);
    if site_sum(mid) > 0
        hi = mid;
    else
        lo = mid;
    end
end
y = 0.5*(lo + hi);
f_Ca = y^2;

f_Mg  = A_Mg * f_Ca;
f_Na  = A_Na * y;
f_K   = A_K  * y;
f_H   = A_H  * y;
f_Al  = A_Al * y^3;

f_tot = f_Ca + f_Mg + f_K + f_Na + f_Al + f_H;
if ~(all(isfinite([f_Ca,f_Mg,f_K,f_Na,f_Al,f_H])) && ...
        all([f_Ca,f_Mg,f_K,f_Na,f_Al,f_H] >= 0) && abs(f_tot - 1) <= 1e-11)
    error('ERW:exchangeClosure', 'Gaines-Thomas exchange fractions do not close.');
end
end

%% =======================================================================
function [f_Ca, f_Mg, f_K, f_Na, f_Al, f_H] = Appelo_Exchange_Fractions( ...
    H, Ca, Mg, K, Na, Al_w, Al_CEC_scale, ...
    K1_Al, K2_Al, K3_Al, K4_Al, ...
    K_Na_H, K_K_H, K_Ca_H2, K_Mg_H2, K_Al_H3, appelo_alpha)
%%% SCEPTER/Appelo (1994) exchange at dissolved concentrations (umol/L).
%%%
%%%     f[X-Na]  = K_Na *[Na] *(f_H/[H] )  * gamma
%%%     f[X-K]   = K_K  *[K]  *(f_H/[H] )  * gamma
%%%     f[X2-Ca] = K_Ca *[Ca] *(f_H/[H] )^2* gamma^2
%%%     f[X2-Mg] = K_Mg *[Mg] *(f_H/[H] )^2* gamma^2
%%%     f[X3-Al] = K_Al *[Al3]*(f_H/[H] )^3* gamma^3
%%%     gamma = 10^(alpha*f_H),  site sum f_H + f[X-Na]+... = 1
%%% Site balance is solved over f_H in (0,1).
if ~(isfinite(H) && H > 0)
    error('ERW:invalidState', 'Appelo exchange requires positive H.');
end
if any(~isfinite([Ca,Mg,K,Na,Al_w])) || any([Ca,Mg,K,Na,Al_w] < 0)
    error('ERW:invalidState', ...
        'Appelo exchange requires finite nonnegative dissolved cations.');
end
if ~(isfinite(Al_CEC_scale) && Al_CEC_scale > 0)
    error('ERW:invalidParam', 'Appelo exchange requires positive Al_CEC_scale.');
end
if ~(isfinite(appelo_alpha) && appelo_alpha >= 0)
    error('ERW:invalidParam', 'Appelo alpha must be a finite nonnegative scalar.');
end
Al3_frac = al3_fraction(H, K1_Al, K2_Al, K3_Al, K4_Al);
Al3 = Al3_frac * max(Al_w, 0) / Al_CEC_scale;

site_residual = @(f_H) appelo_frac_site_residual(f_H, H, Ca, Mg, K, Na, Al3, ...
    K_Na_H, K_K_H, K_Ca_H2, K_Mg_H2, K_Al_H3, appelo_alpha);
f_lo = 1e-24;
f_hi = 1 - 1e-14;
g_lo = site_residual(f_lo);
g_hi = site_residual(f_hi);
if ~(isfinite(g_lo) && isfinite(g_hi) && sign(g_lo) ~= sign(g_hi))
    error('ERW:invalidState', 'No physical Appelo partition at the given solution.');
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
fractions = [f_Ca, f_Mg, f_K, f_Na, f_Al, f_H];
if any(~isfinite(fractions)) || any(fractions < -1e-10) || abs(sum(fractions) - 1) > 1e-8
    error('ERW:invalidState', 'Appelo exchange partition produced an unphysical state.');
end
f_Ca = max(0, f_Ca); f_Mg = max(0, f_Mg); f_K = max(0, f_K);
f_Na = max(0, f_Na); f_Al = max(0, f_Al); f_H = max(0, f_H);
end

function [residual, f_Na, f_K, f_Ca, f_Mg, f_Al] = appelo_frac_site_residual( ...
    f_H, H, Ca, Mg, K, Na, Al3, ...
    K_Na_H, K_K_H, K_Ca_H2, K_Mg_H2, K_Al_H3, appelo_alpha)
gamma = 10^(appelo_alpha*f_H);
rho1 = (f_H/H) * gamma;
f_Na = K_Na_H * Na * rho1;
f_K  = K_K_H  * K  * rho1;
f_Ca = K_Ca_H2* Ca * rho1^2;
f_Mg = K_Mg_H2* Mg * rho1^2;
f_Al = K_Al_H3 * Al3 * rho1^3;
residual = f_H + f_Na + f_K + f_Ca + f_Mg + f_Al - 1;
end

function a3 = al3_fraction(H, K1_Al, K2_Al, K3_Al, K4_Al)
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
