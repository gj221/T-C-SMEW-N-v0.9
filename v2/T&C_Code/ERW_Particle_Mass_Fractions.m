function mass_fraction = ERW_Particle_Mass_Fractions(diameter,cumulative_mass)
%%% Convert a cumulative particle-size curve to discrete class masses.

diameter = double(diameter(:).');
cumulative_mass = double(cumulative_mass(:).');
if isempty(diameter) || numel(diameter) ~= numel(cumulative_mass) || ...
        any(~isfinite(diameter)) || any(~isfinite(cumulative_mass)) || ...
        any(diameter <= 0) || any(diff(diameter) <= 0) || ...
        any(cumulative_mass < 0) || any(diff(cumulative_mass) < 0) || ...
        cumulative_mass(end) <= 0
    error('ERW:invalidParticleDistribution', ...
        'Particle diameters must be positive and ordered; cumulative mass must be finite and nonnegative.');
end

mass_fraction = [cumulative_mass(1),diff(cumulative_mass)];
mass_fraction = mass_fraction/sum(mass_fraction);

end
