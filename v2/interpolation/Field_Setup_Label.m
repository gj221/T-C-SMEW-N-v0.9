function label = Field_Setup_Label(features,indices,setup_index)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Name one connected polygon from whatever the shapefile attributes offer.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for attribute = {'Field','Identifier','Name'}
    if ~isfield(features,attribute{1}); continue; end
    values = string({features(indices).(attribute{1})});
    values = unique(values(strlength(values) > 0),'stable');
    if isempty(values); continue; end
    label = strjoin(values,', ');
    if isscalar(indices) && ~strcmp(attribute{1},'Field')
        label = sprintf('%s %d',label,setup_index);
    end
    return
end
label = sprintf('Polygon %d',setup_index);

end
