%#########################################################################
%   INITIALIZE_ERW_EQUILIBRIUM
%#########################################################################
%%% Initial total dissolved Al
%%% Option 1: Ca-Al Gaines-Thomas exchange
%%% Option 2: Al(OH)3 mineral equilibrium
%%%
%%% The ERW solver stores all Al quantities at Al_CEC_scale times the true
%%% aqueous concentration (conv_Al = 1e3 scaling for numerical conditioning).
%%% Returns Al_w_eq in the same model units.
%%%
%%% Ca-Al Gaines-Thomas form:
%%%   f_Al = (Al3+_model/Al_CEC_scale) * sqrt(f_Ca^3/(K_Ca_Al*Ca^3))
%%%
%%% Dissolution reaction (standard orientation):
%%%   Al(OH)3(s) + 3H+ = Al3+ + 3H2O
%%%   *Kso = [Al3+]_actual / [H+]^3    (mol/L convention)
%%%
%%% Unit conversion to model units (µmol/L stored at Al_CEC_scale × actual):
%%%   Al3+_model = 10^K_AlOH3 * (H_umol/conv_mol)^3 * conv_mol * Al_CEC_scale
%%%
%%% Total dissolved Al from stepwise hydrolysis (Weil & Brady 2017, ch. 9):
%%%   den_H  = H^4 + H^3*K1 + H^2*K1*K2 + H*K1*K2*K3 + K1*K2*K3*K4  [all umol/L]
%%%   Al_w_model = Al3+_model * den_H / H^4
%%%
%%% References:
%%%   Lindsay W.L. (1979) Chemical Equilibria in Soils. Wiley, p. 51.
%%%     Crystalline gibbsite:  log *Kso = 8.04
%%%     Amorphous Al(OH)3:     log *Kso = 10.8
%%%   Nordstrom D.K. & May H.M. (1996) Aqueous equilibrium data for
%%%     mononuclear Al species. Rev. Mineral. 31, 39-80.
%%%     Gibbsite:              log *Kso = 8.11 (+/- 0.09)
%%%   Weil R.R. & Brady N.C. (2017) The Nature and Properties of Soils,
%%%     15th ed., ch. 9 (Al hydrolysis pK chain).
%%%
%%% INPUTS  (all scalar, consistent units)
%%%   Form 1, exchange equilibrium:
%%%   H, Ca, f_Ca, f_Al, K_Ca_Al, K1_Al, K2_Al, K3_Al, K4_Al, Al_CEC_scale
%%%
%%%   Form 2, Al(OH)3 saturation:
%%%   H            hydrogen ion concentration [umol/L]
%%%   K1_Al..K4_Al Al hydrolysis constants [umol/L]  (pK1=5, pK2=5.1, pK3=6.7, pK4=6.2)
%%%   conv_mol     unit conversion mol/L -> umol/L (= 1e6)
%%%   K_AlOH3      log10 *Kso for Al(OH)3 dissolution [mol/L convention].
%%%                Microcrystalline/proto-gibbsite in agricultural soils: 9.66.
%%%                Crystalline gibbsite lower bound: 8.04.  Amorphous upper bound: 10.8.
%%%   Al_CEC_scale ERW solver Al storage scaling (conv_Al = 1e3); result is
%%%                multiplied by this so that Al_w_eq is in model units.
%%% OUTPUT
%%%   Al_w_eq      total dissolved Al [umol/L, stored at Al_CEC_scale x actual];
%%%                NaN if inputs are invalid
function Al_w_eq = Initialize_ERW_Equilibrium(H, varargin)

if nargin == 10
    Ca = varargin{1};
    f_Ca = varargin{2};
    f_Al = varargin{3};
    K_Ca_Al = varargin{4};
    K1_Al = varargin{5};
    K2_Al = varargin{6};
    K3_Al = varargin{7};
    K4_Al = varargin{8};
    Al_CEC_scale = varargin{9};

    if ~(isfinite(H) && H > 0) || ~(isfinite(Ca) && Ca > 0) || ...
            ~(isfinite(f_Ca) && f_Ca > 0) || ~(isfinite(f_Al) && f_Al > 0) || ...
            ~(isfinite(K_Ca_Al) && K_Ca_Al > 0) || ...
            ~(isfinite(Al_CEC_scale) && Al_CEC_scale > 0)
        Al_w_eq = NaN;
        return
    end

    Al3_eq = f_Al * Al_CEC_scale * sqrt(K_Ca_Al * Ca^3 / f_Ca^3);
elseif nargin == 8
    K1_Al = varargin{1};
    K2_Al = varargin{2};
    K3_Al = varargin{3};
    K4_Al = varargin{4};
    conv_mol = varargin{5};
    K_AlOH3 = varargin{6};
    Al_CEC_scale = varargin{7};

    if ~(isfinite(H) && H > 0) || ~(isfinite(conv_mol) && conv_mol > 0) || ...
            ~isfinite(K_AlOH3) || ~(isfinite(Al_CEC_scale) && Al_CEC_scale > 0)
        Al_w_eq = NaN;
        return
    end

    Al3_eq = 10^K_AlOH3 * (H/conv_mol)^3 * conv_mol * Al_CEC_scale;
else
    Al_w_eq = NaN;
    return
end

den_H = H^4 + H^3*K1_Al + H^2*K1_Al*K2_Al + H*K1_Al*K2_Al*K3_Al ...
    + K1_Al*K2_Al*K3_Al*K4_Al;

if ~(isfinite(Al3_eq) && Al3_eq > 0 && isfinite(den_H) && den_H > 0)
    Al_w_eq = NaN;
    return
end

Al_w_eq = Al3_eq * den_H / H^4;
end
