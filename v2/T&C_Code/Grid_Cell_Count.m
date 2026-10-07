function count = Grid_Cell_Count(extent_m,cell_size_m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Number of full equal cells needed to cover a metric extent.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~(isnumeric(extent_m) && isscalar(extent_m) && ...
        isfinite(extent_m) && extent_m >= 0)
    error('Grid_Cell_Count:BadExtent', ...
        'extent_m must be a finite nonnegative scalar in metres.');
end
if ~(isnumeric(cell_size_m) && isscalar(cell_size_m) && ...
        isfinite(cell_size_m) && cell_size_m > 0)
    error('Grid_Cell_Count:BadCellSize', ...
        'cell_size_m must be a finite positive scalar in metres.');
end

ratio = extent_m/cell_size_m;
count = max(1,ceil(ratio));

end
