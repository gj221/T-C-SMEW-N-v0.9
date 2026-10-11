function [solid_total, solid_primary] = ERW_Total_Calcite_Solid( ...
    M_min, mineral, MM_min, solid_secondary)

names = strtrim(string(mineral));
mass = M_min(:);
MM = MM_min(:);
if ~isfinite(solid_secondary) || solid_secondary < 0
    error('ERW:invalidCalciteInventory', ...
        'Secondary calcite inventory must be finite and nonnegative.');
end
if isempty(names)
    solid_primary = 0;
    solid_total = solid_secondary;
    return
end
if numel(names) ~= numel(mass) || numel(MM) ~= numel(mass) || ...
        any(~isfinite(mass)) || any(mass < 0) || ...
        any(~isfinite(MM)) || any(MM <= 0)
    error('ERW:invalidCalciteInventory', ...
        'Calcite solid inventory requires matching finite nonnegative mineral states.');
end

is_calcite = strcmpi(names,"calcite");
solid_primary = sum(mass(is_calcite)./MM(is_calcite));
solid_total = solid_primary + solid_secondary;

end
