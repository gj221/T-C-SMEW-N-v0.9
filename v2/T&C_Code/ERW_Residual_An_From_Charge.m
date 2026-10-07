function [An, state] = ERW_Residual_An_From_Charge(H, CO2_w, Ca, Mg, ...
    Na, K_aq, NH4, NO3, PO4, Al_w, k1, k2, k_w, k_nh3, K1, K2, K3, K4, ...
    Al_CEC_scale)
%ERW_RESIDUAL_AN_FROM_CHARGE Strong-anion residual from aqueous charge balance.
%
%   All dissolved concentrations are [umol/L].  Al_w follows the ERW
%   model scale, so it is divided by Al_CEC_scale before the aqueous charge
%   balance is evaluated.  The returned An is the non-Al residual anion
%   pool [umolc/L]; Al charge remains explicit in the ERW proton equation.

if nargin < 19 || isempty(Al_CEC_scale)
    Al_CEC_scale = 1;
end
if ~(isfinite(Al_CEC_scale) && Al_CEC_scale > 0)
    error('ERW:invalidParam', 'Al_CEC_scale must be positive.');
end

if ~(isfinite(H) && H > 0)
    error('ERW:invalidState', 'Hydrogen ion concentration must be positive.');
end
conc_inputs = [CO2_w, Ca, Mg, Na, K_aq, NH4, NO3, PO4, Al_w];
eq_inputs = [k1, k2, k_w, k_nh3, K1, K2, K3, K4];
if any(~isfinite([conc_inputs, eq_inputs]))
    error('ERW:invalidState', 'Residual-An inputs must be finite.');
end
if any(conc_inputs < 0)
    error('ERW:invalidState', 'Residual-An concentrations must be non-negative.');
end
if any(eq_inputs <= 0)
    error('ERW:invalidParam', 'Residual-An equilibrium constants must be positive.');
end

Al_w_aq = Al_w / Al_CEC_scale;
Al_den = H^4 + H^3*K1 + H^2*K1*K2 + H*K1*K2*K3 + K1*K2*K3*K4;
if ~(isfinite(Al_den) && Al_den > 0)
    error('ERW:invalidState', 'Invalid Al hydrolysis denominator.');
end

OH = k_w/H;
HCO3 = k1*CO2_w/H;
CO3 = k2*k1*CO2_w/H^2;

Al3 = (H^4/Al_den)*Al_w_aq;
AlOH = (H^3*K1/Al_den)*Al_w_aq;
AlOH2 = (H^2*K1*K2/Al_den)*Al_w_aq;
AlOH3 = (H*K1*K2*K3/Al_den)*Al_w_aq;
AlOH4 = max(Al_w_aq - (Al3 + AlOH + AlOH2 + AlOH3), 0);  % floor prevents float-cancellation underflow

Al_net = 3*Al3 + 2*AlOH + AlOH2 - AlOH4;
% Retain the legacy TC-SMEW strong-ion convention used by Alk_tot.
NH4_charge = NH4;
An = H + 2*Ca + 2*Mg + Na + K_aq + NH4_charge + Al_net - ...
    (OH + HCO3 + 2*CO3 + NO3 + PO4);

if nargout > 1
    state = struct();
    state.OH = OH;
    state.HCO3 = HCO3;
    state.CO3 = CO3;
    state.Al3 = Al3;
    state.AlOH = AlOH;
    state.AlOH2 = AlOH2;
    state.AlOH3 = AlOH3;
    state.AlOH4 = AlOH4;
    state.Al_net = Al_net;
    state.Al_w_aq = Al_w_aq;
    state.NH4_charge = NH4_charge;
end
end
