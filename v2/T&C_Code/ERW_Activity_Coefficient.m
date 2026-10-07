function gamma = ERW_Activity_Coefficient(z, I)
%ERW_ACTIVITY_COEFFICIENT  Davies-equation single-ion activity coefficient.
%
%   gamma = ERW_ACTIVITY_COEFFICIENT(z, I) returns the activity coefficient
%   for an ion of charge z at ionic strength I [mol/L], using the Davies
%   (1962) extension of Debye-Huckel:
%
%     log10(gamma) = -A * z^2 * ( sqrt(I)/(1+sqrt(I)) - 0.3*I )
%
%   A = 0.509 at 25 degC (the constant is only weakly temperature
%   dependent over the soil range this model covers -- about 0.503 at
%   10 degC, 0.521 at 35 degC -- so a fixed value is a documented
%   simplification, not an omission). z may be negative or positive; only
%   z^2 enters, matching Davies' own convention. Neutral species (z == 0)
%   return gamma == 1 exactly, since the ionic-strength correction is a
%   long-range electrostatic effect that vanishes for uncharged solutes
%   in this formulation.
%
%   Accepts scalar or array z and I (broadcast against each other).
%   Valid for I up to about 0.5 mol/L. Soil pore water at this site's
%   fertilised conditions runs 0.03-0.08 mol/L, and the calcite-amended
%   columns reach about 0.13: inverting the saturation ratios the
%   2026-09-22 run 12 reported against its own recorded Ca, CO2 and pH
%   puts the post-application ionic strength at 0.126 mol/L. Still well
%   inside the Davies range. Returns 1 (no correction) for I <= 0, since a
%   dilute limit has no ionic-strength effect and I is only ever
%   nonnegative in a physical solution -- I < 0 signals an upstream
%   state error, not a chemistry answer, so it is not silently clamped.
%
%   Reference: Davies, C.W. (1962). Ion Association. Butterworths, London.
%   Same activity-coefficient family SCEPTER uses (scheme_act='davies',
%   scepter.f90); no free parameters, nothing calibrated to this site.

if ~(isnumeric(z) && isnumeric(I) && all(isfinite(z(:))) && all(isfinite(I(:))))
    error('ERW:invalidActivityInput', 'z and I must be finite numeric arrays.');
end
if any(I(:) < 0)
    error('ERW:invalidIonicStrength', 'Ionic strength must be nonnegative; got %.6g.', min(I(:)));
end

A = 0.509;
sqrtI = sqrt(I);
log10gamma = -A .* (z.^2) .* (sqrtI./(1+sqrtI) - 0.3*I);
gamma = 10.^log10gamma;
gamma(I == 0) = 1;
end
