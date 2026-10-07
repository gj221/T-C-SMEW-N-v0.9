function CEC_eff = ERW_CEC_At_H(H,conv_mol,CEC_ref,dCEC_areal,CEC_ref_pH)
%%% Evaluate the local linear variable-charge relation on the areal basis.

if any(~isfinite(H(:))) || any(H(:) <= 0) || ...
        ~(isscalar(conv_mol) && isfinite(conv_mol) && conv_mol > 0) || ...
        ~(isscalar(CEC_ref) && isfinite(CEC_ref) && CEC_ref > 0) || ...
        ~(isscalar(dCEC_areal) && isfinite(dCEC_areal) && dCEC_areal >= 0) || ...
        ~(isscalar(CEC_ref_pH) && isfinite(CEC_ref_pH))
    error('ERW:invalidParam', 'Variable-CEC inputs must be finite with positive H, conv_mol, and CEC_ref.');
end

pH = -log10(H/conv_mol);
CEC_eff = CEC_ref + dCEC_areal*(pH-CEC_ref_pH);
end
