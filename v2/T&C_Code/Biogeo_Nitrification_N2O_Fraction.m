function fraction = Biogeo_Nitrification_N2O_Fraction(PH)
% Fraction of nitrified N emitted as N2O-N after Val Martin et al. (2023).

if ~isreal(PH) || any(~isfinite(PH(:))) || any(PH(:) <= 0 | PH(:) >= 14)
    error('Biogeo_Nitrification_N2O_Fraction:InvalidPH', ...
        'Soil pH must be finite and satisfy 0 < pH < 14.');
end

fraction = min(1,721.86*exp(-2.387*PH));
end
