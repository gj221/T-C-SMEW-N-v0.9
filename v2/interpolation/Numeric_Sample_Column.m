function values = Numeric_Sample_Column(points,column)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Return one sample table column as double, whatever type it was read as.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~istable(points)
    error('Numeric_Sample_Column:BadTable', ...
        'points must be a table of soil samples.');
end
if ~ismember(char(column),points.Properties.VariableNames)
    error('Numeric_Sample_Column:MissingColumn', ...
        'The sample table has no column named %s.',char(column));
end

raw = points.(char(column));
if isnumeric(raw)
    values = double(raw);
elseif islogical(raw)
    values = double(raw);
elseif iscell(raw) || isstring(raw) || ischar(raw)
    values = str2double(string(raw));
else
    error('Numeric_Sample_Column:UnsupportedType', ...
        'Column %s holds %s, which is not a measured quantity.', ...
        char(column),class(raw));
end
values = reshape(values,[],1);

end
