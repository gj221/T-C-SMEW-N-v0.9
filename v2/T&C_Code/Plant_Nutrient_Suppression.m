function suppression = Plant_Nutrient_Suppression(reserve,capacity,nonstorage,element_name)
% Preserve the nutrient-stress equation at an empty-plant boundary.

if ~isfinite(nonstorage) || nonstorage < 0
    error('Plant_Nutrient_Suppression:InvalidState', ...
        '%s nonstorage inventory must be finite and nonnegative.',element_name);
end

reserve_fraction = Plant_Reserve_Fraction(reserve,capacity,element_name);
denominator = 0.6*nonstorage + 0.2*capacity;
if denominator == 0
    suppression = reserve_fraction;
    return
end

suppression = min(max((reserve - 0.8*capacity)/denominator,0),1);
end
