function [points,metadata] = Read_Sample_Points(sample_file,sheet, ...
        x_column,y_column)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read a soil sample table and keep the rows that carry finite coordinates.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2; sheet = 1; end
if nargin < 3 || isempty(x_column); x_column = 'lng (X)'; end
if nargin < 4 || isempty(y_column); y_column = 'lat (Y)'; end
if isstring(sample_file) && isscalar(sample_file)
    sample_file = char(sample_file);
end
if ~(ischar(sample_file) && isrow(sample_file) && ~isempty(sample_file)) || ...
        ~isfile(sample_file)
    error('Read_Sample_Points:MissingSampleFile', ...
        'sample_file must name an existing workbook or delimited file.');
end

[~,~,extension] = fileparts(sample_file);
if any(strcmpi(extension,{'.xlsx','.xlsm','.xls'}))
    points = readtable(sample_file,'Sheet',sheet, ...
        'VariableNamingRule','preserve');
else
    points = readtable(sample_file,'VariableNamingRule','preserve');
    sheet = '';
end
source_row_count = height(points);

x_column = resolve_column(points,x_column,'Read_Sample_Points');
y_column = resolve_column(points,y_column,'Read_Sample_Points');
points.(x_column) = Numeric_Sample_Column(points,x_column);
points.(y_column) = Numeric_Sample_Column(points,y_column);
located = isfinite(points.(x_column)) & isfinite(points.(y_column));
points = points(located,:);
if isempty(points)
    error('Read_Sample_Points:NoLocatedSamples', ...
        'No row of %s carries a finite coordinate pair.',sample_file);
end

metadata.sample_file = sample_file;
metadata.sheet = sheet;
metadata.x_column = x_column;
metadata.y_column = y_column;
metadata.source_row_count = source_row_count;
metadata.located_row_count = height(points);

end

function name = resolve_column(points,requested,caller)
names = string(points.Properties.VariableNames);
exact = find(names == string(requested),1);
if isempty(exact)
    exact = find(strcmpi(strtrim(names),strtrim(string(requested))),1);
end
if isempty(exact)
    error('%s:MissingColumn',caller);
end
name = char(names(exact));
end
