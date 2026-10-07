function spread = Audit_Preapplication_Treatment_Spread( ...
    dates,values,firstAmendment)
dates = dates(:);
values = double(values);
if ~isdatetime(dates) || ~isscalar(firstAmendment) || ...
        ~isdatetime(firstAmendment) || size(values,1) ~= numel(dates) || ...
        isempty(dates)
    error('Audit:PreapplicationInput', ...
        'Treatment dates and values are invalid.');
end
mask = dates < firstAmendment;
if ~any(mask)
    error('Audit:PreapplicationWindow', ...
        'No state exists strictly before the first amendment.');
end
selected = values(mask,:);
spread = max(max(selected,[],2)-min(selected,[],2));
end
