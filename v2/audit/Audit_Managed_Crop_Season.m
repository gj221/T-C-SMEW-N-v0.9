function [growthEnd,harvestRow,harvested,failure] = ...
    Audit_Managed_Crop_Season(dates,lai,plantC,scheduledHarvest)
dates = dates(:);
lai = double(lai(:));
plantC = double(plantC(:));
if ~isdatetime(dates) || ~isscalar(scheduledHarvest) || ...
        ~isdatetime(scheduledHarvest) || numel(lai) ~= numel(dates) || ...
        numel(plantC) ~= numel(dates) || isempty(dates) || ...
        any(~isfinite([lai;plantC])) || any([lai;plantC] < -1e-10)
    error('Audit:ManagedCropInput', ...
        'Managed-crop dates, LAI and plant carbon are invalid.');
end

tol = 1e-10;
transition = find(plantC(1:end-1) > tol & plantC(2:end) <= tol & ...
    lai(2:end) <= tol)+1;
inWindow = dates(transition) >= scheduledHarvest & ...
    dates(transition) <= scheduledHarvest+days(1);
transition = transition(inWindow);

growthEnd = numel(dates);
harvestRow = NaN;
harvested = false;
failure = '';
if ~isempty(transition)
    harvestRow = transition(1);
    growthEnd = harvestRow-1;
    harvested = true;
elseif dates(end) > scheduledHarvest+days(1)
    failure = sprintf('scheduled harvest boundary was not observed by %s', ...
        datestr(scheduledHarvest+days(1),29));
end
end
