function factor = Biogeo_Denitrification_WFPS_Response(WFPS)
if ~isreal(WFPS) || any(~isfinite(WFPS(:))) || ...
        any(WFPS(:) < 0 | WFPS(:) > 1)
    error('Biogeo_Denitrification_WFPS_Response:InvalidWFPS', ...
        'Water-filled pore space must be finite and satisfy 0 <= WFPS <= 1.');
end

pDEN = 13.^(2.2*WFPS);
factor = 1.4*13.^(-17./pDEN);
end
