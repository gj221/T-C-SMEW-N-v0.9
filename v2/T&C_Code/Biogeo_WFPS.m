function WFPS = Biogeo_WFPS(theta,theta_sat)
if ~isreal(theta) || any(~isfinite(theta(:))) || ...
        ~(isscalar(theta_sat) && isreal(theta_sat) && ...
        isfinite(theta_sat) && theta_sat > 0)
    error('Biogeo_WFPS:InvalidWaterContent', ...
        'Water contents must be finite and saturated water content must be positive.');
end

WFPS = theta./theta_sat;
if any(WFPS(:) < 0 | WFPS(:) > 1)
    error('Biogeo_WFPS:InvalidWaterContent', ...
        'Water-filled pore space must satisfy 0 <= WFPS <= 1.');
end
end
