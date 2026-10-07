function VT = Biogeo_Mobile_Water_Volume(Zbio,theta,theta_hyg)
if ~isscalar(Zbio) || ~isfinite(Zbio) || Zbio <= 0
    error('Biogeo_Mobile_Water_Volume:InvalidDepth', ...
        'Biogeochemical depth must be finite and positive.');
end
if any(~isfinite(theta(:))) || any(~isfinite(theta_hyg(:))) || ...
        ~isequal(size(theta),size(theta_hyg))
    error('Biogeo_Mobile_Water_Volume:InvalidWaterContent', ...
        'Water contents must be finite arrays with matching dimensions.');
end

VT = Zbio*max(0,theta-theta_hyg);
end
