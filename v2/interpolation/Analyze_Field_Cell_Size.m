function analysis = Analyze_Field_Cell_Size(field_polygon,cell_sizes_m, ...
        minimum_cell_fractions)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sweep candidate cell widths and score each one against both grid rules.
%%
%% Two different things get called 95% coverage and they do not measure the
%% same error.
%%
%%   aggregate coverage    A/(N*h^2) over every cell with positive overlap.
%%                         Bounds how much simulated area the grid adds in
%%                         total, but says nothing about any single cell: a
%%                         cell 10% inside the polygon still passes.
%%   per-cell fraction     overlap/h^2 for each cell taken on its own. A
%%                         threshold on it bounds how much of any one
%%                         simulated column is physically outside the field,
%%                         and so bounds the scaling correction that column
%%                         needs. It is met by keeping only the cells that
%%                         pass, which leaves a thin band of the polygon
%%                         unsimulated.
%%
%% The band is what the capture column measures, and it is the real cost of
%% the per-cell rule. Overlap fractions are computed once per width and then
%% scored at every requested threshold.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3 || isempty(minimum_cell_fractions)
    minimum_cell_fractions = 0.95;
end
if ~(isa(field_polygon,'polyshape') && isscalar(field_polygon)) || ...
        area(field_polygon) <= 0
    error('Analyze_Field_Cell_Size:EmptyField', ...
        'field_polygon must have finite positive area.');
end
cell_sizes_m = reshape(double(cell_sizes_m),[],1);
if isempty(cell_sizes_m) || ~all(isfinite(cell_sizes_m) & cell_sizes_m > 0)
    error('Analyze_Field_Cell_Size:BadCellSizes', ...
        'cell_sizes_m must hold finite positive widths in metres.');
end
minimum_cell_fractions = reshape(double(minimum_cell_fractions),[],1);
if isempty(minimum_cell_fractions) || ...
        ~all(minimum_cell_fractions > 0 & minimum_cell_fractions <= 1)
    error('Analyze_Field_Cell_Size:BadThreshold', ...
        'minimum_cell_fractions must lie in (0,1].');
end

field_area_m2 = area(field_polygon);
row_count = numel(cell_sizes_m)*numel(minimum_cell_fractions);
cell_size_column = zeros(row_count,1);
threshold = zeros(row_count,1);
grid_columns = zeros(row_count,1);
grid_rows = zeros(row_count,1);
cells_any = zeros(row_count,1);
cells_strict = zeros(row_count,1);
aggregate_coverage = zeros(row_count,1);
captured_fraction = zeros(row_count,1);
strict_outside_fraction = zeros(row_count,1);
minimum_overlap = zeros(row_count,1);
mean_overlap = zeros(row_count,1);
missed_area_m2 = zeros(row_count,1);
missed_band_m = zeros(row_count,1);

boundary_length_m = perimeter(field_polygon);
row = 0;
for k = 1:numel(cell_sizes_m)
    cell_size_m = cell_sizes_m(k);
    oriented_polygon = Orient_Field_Polygon_Grid(field_polygon,cell_size_m);
    [overlap_fraction,nx,ny] = Cell_Overlap_Fractions( ...
        oriented_polygon,cell_size_m);
    retained = overlap_fraction > 0;
    cell_area_m2 = cell_size_m^2;
    for j = 1:numel(minimum_cell_fractions)
        row = row+1;
        strict = overlap_fraction >= minimum_cell_fractions(j);
        captured_area_m2 = sum(overlap_fraction(strict))*cell_area_m2;

        cell_size_column(row) = cell_size_m;
        threshold(row) = minimum_cell_fractions(j);
        grid_columns(row) = nx;
        grid_rows(row) = ny;
        cells_any(row) = sum(retained);
        cells_strict(row) = sum(strict);
        aggregate_coverage(row) = ...
            field_area_m2/(cells_any(row)*cell_area_m2);
        captured_fraction(row) = captured_area_m2/field_area_m2;
        missed_area_m2(row) = field_area_m2-captured_area_m2;
        missed_band_m(row) = missed_area_m2(row)/boundary_length_m;
        if cells_strict(row) > 0
            strict_outside_fraction(row) = ...
                (cells_strict(row)*cell_area_m2-captured_area_m2)/ ...
                (cells_strict(row)*cell_area_m2);
        else
            strict_outside_fraction(row) = NaN;
        end
        minimum_overlap(row) = min(overlap_fraction(retained));
        mean_overlap(row) = mean(overlap_fraction(retained));
    end
end

analysis = table(cell_size_column,threshold,grid_columns,grid_rows, ...
    cells_any,cells_strict,aggregate_coverage,captured_fraction, ...
    strict_outside_fraction,minimum_overlap,mean_overlap,missed_area_m2, ...
    missed_band_m,repmat(field_area_m2,row_count,1), ...
    repmat(boundary_length_m,row_count,1), ...
    'VariableNames',{'CellSizeM','MinimumCellFraction','GridColumns', ...
    'GridRows','CellsAny','CellsStrict','AggregateCoverage', ...
    'CapturedFraction','StrictOutsideFraction','MinOverlapAny', ...
    'MeanOverlapAny','MissedAreaM2','MissedBandM','FieldAreaM2', ...
    'PerimeterM'});

end
