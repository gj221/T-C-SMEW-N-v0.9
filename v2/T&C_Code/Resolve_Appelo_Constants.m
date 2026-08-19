%#########################################################################
%   Function RESOLVE_APPELO_CONSTANTS                                   #
%#########################################################################
%%% SCEPTER/Appelo (1994) H-referenced selectivity constants.
%%%
%%%     .logKH  [log10 KH-Na, KH-K, KH-Ca, KH-Mg, KH-Al],
%%%             SCEPTER defaults [5.9 4.8 10.47 10.786 16.47]
%%%     .alpha  Appelo exponent, SCEPTER default 3.4
%%% logKH is the H-over-cation selectivity of the FULL exchange
%%% reaction, e.g. Ca2+ + 2X-H = X2-Ca + 2H+ has log K = -10.47
%%% (derived in scepter.f90 as log KCa\Na*2 - log KH\Na*2; the Z
%%% multiplicity is already inside the tabulated value).  SCEPTER
%%% therefore uses the mass action (scepter.f90 keqiex_all =
%%% 10^(-logkhaq); get_maqads_all_v3 site-sum solve)
%%%     f[X-ς] = 10^(-logKH) * [ς] * (f_H/[H+])^Z * 10^(alpha*Z*f_H)
%%% where the Appelo (1994) proton-loading term ("log KHX =
%%% 5.9 - 3.4 f(H)") enters ONLY through the gamma factor
%%% 10^(alpha*Z*f_H) applied at use -- nothing of alpha is folded into
%%% the constant, and the sign of logKH is NEGATIVE (H is preferred:
%%% at pH 4 an Appelo organic exchanger is H-loaded, and at Kantola
%%% porewater conditions, pH ~6, the equilibrium is base-saturated
%%% with f_Al ~ 1e-5).  The previous form 10^(+Z*(logKH - alpha))
%%% inverted the selectivity by 10^(2*logKH-Z*alpha) per cation
%%% (10^56 for Al), which Al-saturated the exchanger at pH 6 and made
%%% every base-cation pool thermodynamically unhostable.
%%% OUTPUT (umol/L concentration basis):
%%%   K_Na_H  = 10^(-logKH(1))              [-]
%%%   K_K_H   = 10^(-logKH(2))              [-]
%%%   K_Ca_H2 = 10^(-logKH(3)) * 1e6        [umol/L]
%%%   K_Mg_H2 = 10^(-logKH(4)) * 1e6        [umol/L]
%%%   K_Al_H3 = 10^(-logKH(5)) * 1e12       [(umol/L)^2]
%%%   appelo_alpha                          [-]
function [K_Na_H, K_K_H, K_Ca_H2, K_Mg_H2, K_Al_H3, appelo_alpha] = ...
    Resolve_Appelo_Constants(appelo_params)

appelo_logKH_default = [5.9 4.8 10.47 10.786 16.47];   %% Na K Ca Mg Al
if nargin < 1 || isempty(appelo_params)
    appelo_logKH = appelo_logKH_default;
    appelo_alpha = 3.4;
else
    if ~isstruct(appelo_params) || ~isscalar(appelo_params)
        error('ERW:invalidParam', 'Appelo parameters must be a scalar struct.');
    end
    if isfield(appelo_params, 'logKH')
        appelo_logKH = appelo_params.logKH;
        if ~isnumeric(appelo_logKH) || ~isreal(appelo_logKH) || ...
                numel(appelo_logKH) ~= 5 || any(~isfinite(appelo_logKH))
            error('ERW:invalidParam', ...
                'Appelo logKH must contain five finite real values.');
        end
        appelo_logKH = appelo_logKH(:).';
    else
        appelo_logKH = appelo_logKH_default;
    end
    if isfield(appelo_params, 'alpha')
        appelo_alpha = appelo_params.alpha;
        if ~isnumeric(appelo_alpha) || ~isreal(appelo_alpha) || ...
                ~isscalar(appelo_alpha) || ~isfinite(appelo_alpha) || appelo_alpha < 0
            error('ERW:invalidParam', ...
                'Appelo alpha must be a finite nonnegative scalar.');
        end
    else
        appelo_alpha = 3.4;
    end
end

K_Na_H  = 10^(-appelo_logKH(1));
K_K_H   = 10^(-appelo_logKH(2));
K_Ca_H2 = 10^(-appelo_logKH(3)) * 1e6;
K_Mg_H2 = 10^(-appelo_logKH(4)) * 1e6;
K_Al_H3 = 10^(-appelo_logKH(5)) * 1e12;
constants = [K_Na_H,K_K_H,K_Ca_H2,K_Mg_H2,K_Al_H3];
if any(~isfinite(constants)) || any(constants <= 0)
    error('ERW:invalidParam', ...
        'Appelo logKH values produce non-finite or zero selectivity constants.');
end

end
