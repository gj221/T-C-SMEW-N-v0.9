function factor = Biogeo_Nitrification_PH_Response(PH)
if ~isreal(PH) || any(~isfinite(PH(:))) || any(PH(:) <= 0 | PH(:) >= 14)
    error('Biogeo_Nitrification_PH_Response:InvalidPH', ...
        'Soil pH must be finite and strictly between 0 and 14.');
end

factor = 0.56 + atan(pi*0.45*(PH-5))/pi;
factor = min(1,max(0,factor));
end
