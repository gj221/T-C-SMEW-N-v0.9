function [cell_size_m,calibration,oriented_polygon] = Calibrate_Field_Grid_Cell_Size(field_polygon)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Minimize retained full cells subject to more than 95% polygon coverage.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~(isa(field_polygon,'polyshape') && isscalar(field_polygon)) || ...
        isempty(field_polygon.Vertices) || area(field_polygon) <= 0
    error('Calibrate_Field_Grid_Cell_Size:EmptyField', ...
        'field_polygon must have finite positive area.');
end

target_fraction = 0.95;
field_area_m2 = area(field_polygon);
vertices = field_polygon.Vertices;
vertices = vertices(all(isfinite(vertices),2),:);
extent_m = max(vertices,[],1)-min(vertices,[],1);
[oriented_polygon,base_calibration] = Orient_Field_Polygon_Grid( ...
    field_polygon,2*hypot(extent_m(1),extent_m(2)));
width_m = base_calibration.bounding_width_m;
height_m = base_calibration.bounding_height_m;
best = struct('cell_size_m',[],'calibration',struct(), ...
    'oriented_polygon',oriented_polygon,'topology',struct());

current_size_m = min(width_m,height_m);
width_divisions = floor(width_m/current_size_m+1e-12)+1;
height_divisions = floor(height_m/current_size_m+1e-12)+1;
while true
    candidate = evaluate_size(current_size_m);
    best = choose_better(best,candidate,target_fraction);

    next_width_size_m = width_m/width_divisions;
    next_height_size_m = height_m/height_divisions;
    next_size_m = max(next_width_size_m,next_height_size_m);
    if ~isempty(best.cell_size_m)
        minimum_future_cells = ceil(field_area_m2/next_size_m^2-1e-12);
        if minimum_future_cells > best.topology.N_cells
            break
        end
    end
    if ~(isfinite(next_size_m) && next_size_m > ...
            eps(max(width_m,height_m)))
        error('Calibrate_Field_Grid_Cell_Size:NoFeasibleGrid', ...
            'A grid with more than 95%% polygon coverage could not be resolved.');
    end
    transition_tol = 1e-12*next_size_m;
    if abs(next_width_size_m-next_size_m) <= transition_tol
        width_divisions = width_divisions+1;
    end
    if abs(next_height_size_m-next_size_m) <= transition_tol
        height_divisions = height_divisions+1;
    end
    current_size_m = next_size_m;
end

best_count = best.topology.N_cells;
cell_size_high = best.cell_size_m;
cell_size_low = 0.5*cell_size_high;
candidate_low = evaluate_size(cell_size_low);
while candidate_low.topology.N_cells <= best_count && ...
        candidate_low.calibration.field_fraction > target_fraction
    best = choose_better(best,candidate_low,target_fraction);
    cell_size_high = cell_size_low;
    cell_size_low = 0.5*cell_size_low;
    candidate_low = evaluate_size(cell_size_low);
end
for iteration = 1:30
    candidate = evaluate_size(0.5*(cell_size_low+cell_size_high));
    if candidate.topology.N_cells <= best_count && ...
            candidate.calibration.field_fraction > target_fraction
        cell_size_high = candidate.cell_size_m;
        best = choose_better(best,candidate,target_fraction);
    else
        cell_size_low = candidate.cell_size_m;
    end
end

try
    exact_topology = Build_Polygon_Grid_Topology( ...
        oriented_polygon,best.cell_size_m);
catch ME
    if ~strcmp(ME.identifier,'Build_Polygon_Grid_Topology:AreaMismatch')
        rethrow(ME)
    end
    safe_candidate = evaluate_size(best.cell_size_m*(1+1e-7));
    if safe_candidate.calibration.field_fraction <= target_fraction || ...
            safe_candidate.topology.N_cells > best.topology.N_cells
        rethrow(ME)
    end
    best = safe_candidate;
    exact_topology = Build_Polygon_Grid_Topology( ...
        oriented_polygon,best.cell_size_m);
end
exact_fraction = exact_topology.field_area_m2/exact_topology.simulated_area_m2;
if exact_topology.N_cells ~= best.topology.N_cells || ...
        exact_fraction <= target_fraction
    error('Calibrate_Field_Grid_Cell_Size:CalibrationMismatch', ...
        'The exact polygon intersections do not satisfy the calibrated grid.');
end
best.topology = exact_topology;
best = update_calibration(best);
cell_size_m = best.cell_size_m;
calibration = best.calibration;
calibration.mode = 'polygon_geometry';
calibration.required_field_fraction = target_fraction;
calibration.orientation_objective = ...
    'minimum_retained_cells_then_outside_area';

    function candidate = evaluate_size(size_m)
        candidate.oriented_polygon = oriented_polygon;
        candidate.calibration = base_calibration;
        candidate.calibration.cell_size_m = size_m;
        candidate.cell_size_m = size_m;
        candidate.topology = Count_Polygon_Grid_Cells( ...
            oriented_polygon,size_m);
        if candidate.topology.simulated_area_m2 < field_area_m2
            candidate.calibration.field_fraction = 0;
            return
        end
        candidate = update_calibration(candidate);
    end

    function candidate = update_calibration(candidate)
        candidate.calibration.grid_nx = candidate.topology.main_grid_nx;
        candidate.calibration.grid_ny = candidate.topology.main_grid_ny;
        candidate.calibration.grid_slot_count = candidate.topology.grid_slot_count;
        candidate.calibration.grid_cell_count = candidate.topology.N_cells;
        candidate.calibration.simulated_area_m2 = ...
            candidate.topology.simulated_area_m2;
        candidate.calibration.field_area_m2 = candidate.topology.field_area_m2;
        candidate.calibration.field_fraction = ...
            candidate.topology.field_area_m2/candidate.topology.simulated_area_m2;
        candidate.calibration.outside_grid_area_m2 = ...
            candidate.topology.simulated_area_m2-candidate.topology.field_area_m2;
        candidate.calibration.outside_grid_fraction = ...
            candidate.calibration.outside_grid_area_m2/ ...
            candidate.topology.simulated_area_m2;
    end

end

function best = choose_better(best,candidate,target_fraction)
if candidate.calibration.field_fraction <= target_fraction
    return
end
if isempty(best.cell_size_m) || ...
        candidate.topology.N_cells < best.topology.N_cells || ...
        (candidate.topology.N_cells == best.topology.N_cells && ...
        candidate.topology.simulated_area_m2 < ...
        best.topology.simulated_area_m2-1e-9)
    best = candidate;
end
end
