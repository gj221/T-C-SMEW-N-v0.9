function [field_mean,field_total] = Field_Aggregate_Areal(cell_values,topo)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Aggregate cell areal values to an in-field mean and extensive total
%% over the whole field (Field_Aggregation_Areas).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[area_m2,weights,~,field_area_m2] = Field_Aggregation_Areas(topo);
if ~(isnumeric(cell_values) && ismatrix(cell_values)) || ...
        size(cell_values,2) ~= numel(area_m2)
    error('Field_Aggregate_Areal:SizeMismatch', ...
        'cell_values must have one column per topology cell.');
end

field_mean = cell_values*weights.';
field_total = field_mean*field_area_m2;

end
