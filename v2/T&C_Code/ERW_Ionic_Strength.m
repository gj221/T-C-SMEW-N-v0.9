function I = ERW_Ionic_Strength(Ca, Mg, Na, K, Al, H, NH4, NO3, HCO3, CO3, An, conv_mol, ...
    PO4, AlOH, AlOH2, AlOH4, OH)
%ERW_IONIC_STRENGTH  Ionic strength of the ERW pore water, mol/L.
%
%   I = ERW_IONIC_STRENGTH(Ca, Mg, Na, K, Al, H, NH4, NO3, HCO3, CO3, An,
%   conv_mol) computes I = 0.5 * sum(z_i^2 * c_i) in mol/L from the same
%   pore-water ion concentrations [umol-conv/L]. An is the lumped residual anion
%   ERW_function seeds at charge -1 (ERW_pre_processing.m,
%   Invert_Bulk_pH.m); it stands in for chloride, sulfate and other
%   background anions this model does not track individually, so it is
%   included here at z = -1 for consistency with how it enters charge.
%   PO4 is monovalent in the present ERW charge convention. The three charged
%   Al hydroxides and OH are supplied after conv_mol; the legacy 12-input
%   call omits them. Neutral CO2(aq) and Al(OH)3 do not enter ionic strength.
%
%   All concentration inputs are on the model's internal umol-conv/L
%   scale; conv_mol (1e6 in this model) converts to mol/L before the sum.
%   Al is expected on the SAME nmol-vs-umol basis the caller already
%   uses for that species elsewhere in ERW_function (i.e. divided by
%   Al_CEC_scale before this call, exactly as ERW_Residual_An_From_Charge
%   requires of Al_w) -- this function does not know about Al_CEC_scale
%   and will silently compute the wrong I if handed the scaled value.
%
%   Returns a nonnegative scalar or array matching the broadcast shape of
%   the inputs.

if nargin == 12
    [PO4,AlOH,AlOH2,AlOH4,OH] = deal(0);
elseif nargin ~= 17
    error('ERW:invalidIonicStrengthInput','Supply all five minor-ion concentrations.');
end
vals = [Ca(:); Mg(:); Na(:); K(:); Al(:); H(:); NH4(:); NO3(:); HCO3(:); CO3(:); An(:); ...
    PO4(:); AlOH(:); AlOH2(:); AlOH4(:); OH(:)];
if ~(isnumeric(vals) && all(isfinite(vals)))
    error('ERW:invalidIonicStrengthInput', 'All concentration inputs must be finite.');
end
%% A value within rounding distance of zero is a cancellation residue, not
%% a negative amount: Al(OH)4- is the dissolved Al left after its other
%% species, and a pool emptied by a flux of its own size can end a few
%% ulps below zero. A 35-year 2-D run stopped here on -4.5e-16 on
%% 5 October 2026. Inputs above -roundoff_tol (1e-10 umol-conv/L, that is
%% 1e-16 mol/L) therefore count as zero; anything more negative is still
%% an error. Nonnegative inputs take the same path as before, so no
%% result that got past this check changes.
roundoff_tol = 1e-10;
if any(vals < -roundoff_tol)
    error('ERW:negativeConcentration', 'Ionic strength inputs must be nonnegative concentrations; most negative is %.6g.', min(vals));
end
if any(vals < 0)
    [Ca,Mg,Na,K,Al,H,NH4,NO3,HCO3,CO3,An,PO4,AlOH,AlOH2,AlOH4,OH] = clamp_to_zero( ...
        Ca,Mg,Na,K,Al,H,NH4,NO3,HCO3,CO3,An,PO4,AlOH,AlOH2,AlOH4,OH);
end
if ~(isnumeric(conv_mol) && isscalar(conv_mol) && isfinite(conv_mol) && conv_mol > 0)
    error('ERW:invalidConvMol', 'conv_mol must be a positive finite scalar.');
end

z2c = 4*Ca + 4*Mg + Na + K + 9*Al + H + NH4 + NO3 + HCO3 + 4*CO3 + An + ...
    PO4 + 4*AlOH + AlOH2 + AlOH4 + OH;
I = 0.5 * z2c / conv_mol;
end

function varargout = clamp_to_zero(varargin)
varargout = cellfun(@(x) max(x,0),varargin,'UniformOutput',false);
end
