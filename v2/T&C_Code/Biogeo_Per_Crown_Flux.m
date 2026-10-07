function per_crown = Biogeo_Per_Crown_Flux(flux,suppression,crown)
% Convert a ground-area flux to PFT area without supplying inactive PFTs.

if ~isequal(size(flux),size(suppression),size(crown))
    error('Biogeo_Per_Crown_Flux:SizeMismatch', ...
        'Flux, suppression and crown arrays must have identical sizes.');
end
if ~isreal(flux) || any(~isfinite(flux(:))) || any(flux(:) < 0)
    error('Biogeo_Per_Crown_Flux:InvalidFlux', ...
        'Ground-area fluxes must be finite and nonnegative.');
end
if ~isreal(suppression) || any(~isfinite(suppression(:))) || ...
        any(suppression(:) < 0 | suppression(:) > 1)
    error('Biogeo_Per_Crown_Flux:InvalidSuppression', ...
        'Suppression fractions must be finite and lie in [0,1].');
end
if ~isreal(crown) || any(~isfinite(crown(:))) || any(crown(:) < 0)
    error('Biogeo_Per_Crown_Flux:InvalidCrown', ...
        'Crown fractions must be finite and nonnegative.');
end

per_crown = zeros(size(flux));
active = crown > 0;
per_crown(active) = flux(active).*(1-suppression(active))./crown(active);
end
