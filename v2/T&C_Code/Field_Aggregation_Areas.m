function [area_m2,weights,fraction,field_area_m2] = Field_Aggregation_Areas(topo)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Return the in-field accounting area, weights and full-cell fractions,
%% and the field area that totals are for: the overlaps plus, on a
%% majority-rule grid (Build_Polygon_Grid_Topology), the field land in the
%% dropped squares, which the totals carry at the field mean.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isstruct(topo) || ~isfield(topo,'cell_area_m2')
    error('Field_Aggregation_Areas:BadTopology', ...
        'topo.cell_area_m2 is required.');
end
full_area_m2 = topo.cell_area_m2(:);
if any(~isfinite(full_area_m2) | full_area_m2 <= 0)
    error('Field_Aggregation_Areas:BadCellArea', ...
        'Full simulated cell areas must be finite and positive.');
end

if isfield(topo,'cell_overlap_area_m2')
    area_m2 = topo.cell_overlap_area_m2(:);
else
    area_m2 = full_area_m2;
end
if numel(area_m2) ~= numel(full_area_m2) || any(~isfinite(area_m2))
    error('Field_Aggregation_Areas:BadOverlapArea', ...
        'Cell overlap areas must be finite and match cell_area_m2.');
end

tol = 1e-9*max(full_area_m2);
if any(area_m2 < -tol | area_m2 > full_area_m2+tol)
    error('Field_Aggregation_Areas:BadOverlapArea', ...
        'Each overlap area must lie between zero and its full cell area.');
end
area_m2 = min(full_area_m2,max(0,area_m2));
total_area_m2 = sum(area_m2);
if total_area_m2 <= 0
    error('Field_Aggregation_Areas:EmptyField', ...
        'At least one cell must have positive field overlap.');
end

weights = reshape(area_m2/total_area_m2,1,[]);
fraction = area_m2./full_area_m2;
field_area_m2 = total_area_m2;
if isfield(topo,'field_missed_area_m2')
    missed_area_m2 = topo.field_missed_area_m2;
    if ~(isnumeric(missed_area_m2) && isscalar(missed_area_m2) && ...
            isfinite(missed_area_m2) && missed_area_m2 >= 0)
        error('Field_Aggregation_Areas:BadMissedArea', ...
            'field_missed_area_m2 must be a finite nonnegative scalar.');
    end
    field_area_m2 = total_area_m2+missed_area_m2;
end

end
