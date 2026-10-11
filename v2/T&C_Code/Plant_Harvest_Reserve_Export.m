function plot_export = Plant_Harvest_Reserve_Export(plot_export,reserve_export,crown_fraction)
% Convert harvested reserve mass from crown area to plot area once.

if ~isequal(size(plot_export),size(reserve_export)) || ...
        any(~isfinite(plot_export)) || any(~isfinite(reserve_export)) || ...
        any(plot_export < 0) || any(reserve_export < 0)
    error('Plant_Harvest_Reserve_Export:InvalidInventory', ...
        'Harvest reserve inventories must be finite, nonnegative and conformable.');
end
if ~isscalar(crown_fraction) || ~isfinite(crown_fraction) || ...
        crown_fraction < 0 || crown_fraction > 1
    error('Plant_Harvest_Reserve_Export:InvalidCrownFraction', ...
        'Crop crown fraction must be finite and within [0,1].');
end

plot_export = plot_export + crown_fraction*reserve_export;
end
