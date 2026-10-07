function topo = Build_Projected_Field_Topology( ...
        field_polygon,cell_size_m,metadata,source_feature_indices, ...
        coverage_rule)
%% coverage_rule: which full squares are simulated.
%%   a fraction f   every square that touches the polygon, with more than f
%%                  of their total area inside it (default 0.95 -- the
%%                  production rule that no simulated column is mostly
%%                  outside the field); without cell_size_m the width is
%%                  calibrated to the fewest such squares.  Lowering f is a
%%                  DIAGNOSTIC choice only: coarser cells can then be
%%                  built, but the retained grid covers less of the polygon
%%                  and the field ledgers are for that interior core, not
%%                  for the field.
%%   'majority'     every square more than half inside the polygon, on the
%%                  lattice (rotation and offset) that minimises the
%%                  misrepresented area (Place_Field_Polygon_Grid).
%%                  cell_size_m is required: the width comes from the
%%                  terrain (Run_Topographic_Grid_Convergence), not from
%%                  coverage.  The field land in the dropped squares is
%%                  unsimulated and enters the field totals at the field
%%                  mean (Field_Aggregation_Areas).
if nargin < 5 || isempty(coverage_rule)
    coverage_rule = 0.95;
end
majority = (ischar(coverage_rule) || (isstring(coverage_rule) && ...
    isscalar(coverage_rule))) && strcmp(char(coverage_rule),'majority');
if ~majority && ~(isnumeric(coverage_rule) && isscalar(coverage_rule) && ...
        isfinite(coverage_rule) && coverage_rule > 0 && coverage_rule < 1)
    error('Build_Projected_Field_Topology:BadCoverageRule', ...
        'The coverage rule must be a fraction in (0,1) or ''majority''.');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calibrate and annotate one connected projected field-grid setup.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if majority
    if isempty(cell_size_m)
        error('Build_Projected_Field_Topology:WidthRequired', ...
            ['The majority rule takes a prescribed cell width; set it from ' ...
             'the terrain (Run_Topographic_Grid_Convergence).']);
    end
    [grid_polygon,grid_calibration,grid_offset_m] = ...
        Place_Field_Polygon_Grid(field_polygon,cell_size_m);
    grid_calibration.mode = 'prescribed';
    topo = Build_Polygon_Grid_Topology( ...
        grid_polygon,cell_size_m,'majority',grid_offset_m);
    captured_area_m2 = sum(topo.cell_overlap_area_m2);
    grid_calibration.grid_nx = topo.main_grid_nx;
    grid_calibration.grid_ny = topo.main_grid_ny;
    grid_calibration.grid_slot_count = topo.grid_slot_count;
    grid_calibration.grid_cell_count = topo.N_cells;
    grid_calibration.simulated_area_m2 = topo.simulated_area_m2;
    grid_calibration.field_area_m2 = topo.field_area_m2;
    grid_calibration.field_fraction = captured_area_m2/topo.simulated_area_m2;
    grid_calibration.outside_grid_area_m2 = ...
        topo.simulated_area_m2-captured_area_m2;
    grid_calibration.outside_grid_fraction = ...
        grid_calibration.outside_grid_area_m2/topo.simulated_area_m2;
    grid_calibration.coverage_rule = 'majority';
    grid_calibration.missed_field_area_m2 = topo.field_missed_area_m2;
    grid_calibration.captured_field_fraction = ...
        captured_area_m2/topo.field_area_m2;
    grid_calibration.misrepresented_area_m2 = ...
        grid_calibration.outside_grid_area_m2+topo.field_missed_area_m2;
    grid_calibration.misrepresented_fraction = ...
        grid_calibration.misrepresented_area_m2/topo.field_area_m2;
else
    if isempty(cell_size_m)
        [cell_size_m,grid_calibration,grid_polygon] = ...
            Calibrate_Field_Grid_Cell_Size(field_polygon);
    else
        [grid_polygon,grid_calibration] = ...
            Orient_Field_Polygon_Grid(field_polygon,cell_size_m);
        grid_calibration.mode = 'prescribed';
    end
    topo = Build_Polygon_Grid_Topology(grid_polygon,cell_size_m);
    grid_calibration.grid_nx = topo.main_grid_nx;
    grid_calibration.grid_ny = topo.main_grid_ny;
    grid_calibration.grid_slot_count = topo.grid_slot_count;
    grid_calibration.grid_cell_count = topo.N_cells;
    grid_calibration.simulated_area_m2 = topo.simulated_area_m2;
    grid_calibration.field_area_m2 = topo.field_area_m2;
    grid_calibration.field_fraction = topo.field_area_m2/topo.simulated_area_m2;
    grid_calibration.outside_grid_area_m2 = ...
        topo.simulated_area_m2-topo.field_area_m2;
    grid_calibration.outside_grid_fraction = ...
        grid_calibration.outside_grid_area_m2/topo.simulated_area_m2;
    grid_calibration.required_field_fraction = coverage_rule;
    if grid_calibration.field_fraction <= grid_calibration.required_field_fraction
        error('Build_Projected_Field_Topology:CoverageBelowTarget', ...
            ['%.6g m cells place only %.3f%% of retained full-cell area inside ' ...
             'the polygon; more than %.0f%% is required.'], ...
            cell_size_m,100*grid_calibration.field_fraction, ...
            100*coverage_rule);
    end
end
topo.grid_calibration = grid_calibration;
topo.grid_rotation_degrees = grid_calibration.rotation_degrees;
topo.projected_grid_origin_m = grid_calibration.projected_grid_origin_m;
[topo.cell_xc_projected_m,topo.cell_yc_projected_m] = ...
    Grid_Local_To_Projected(topo.cell_xc,topo.cell_yc,topo);
[topo.field_boundary_projected_x_m,topo.field_boundary_projected_y_m] = ...
    Grid_Local_To_Projected( ...
    topo.field_boundary_x_m,topo.field_boundary_y_m,topo);
topo.source_shapefile = metadata.source_shapefile;
topo.source_crs_name = metadata.source_crs_name;
topo.projected_crs_name = metadata.projected_crs_name;
topo.projected_crs_epsg = metadata.projected_crs_epsg;
topo.shapefile_feature_count = metadata.shapefile_feature_count;
topo.shapefile_included_feature_count = numel(metadata.feature_indices);
topo.source_feature_indices = source_feature_indices;
[is_member,member_index] = ismember( ...
    source_feature_indices,metadata.feature_indices);
assert(all(is_member),'Build_Projected_Field_Topology:FeatureMismatch', ...
    'Region features must occur in the selected shapefile features.');
topo.source_finite_vertex_count = sum( ...
    metadata.source_finite_vertex_count(member_index));
topo.projected_polygon_vertex_count = sum( ...
    metadata.projected_polygon_vertex_count(member_index));
topo.projection_roundtrip_max_error_m = max( ...
    metadata.projection_roundtrip_max_error_m(member_index));
topo.source_to_polygon_max_error_m = max( ...
    metadata.source_to_polygon_max_error_m(member_index));

target_crs = projcrs(metadata.projected_crs_epsg);
[latitude,longitude] = projinv(target_crs, ...
    topo.field_boundary_projected_x_m,topo.field_boundary_projected_y_m);
valid = isfinite(latitude) & isfinite(longitude);
topo.field_latitude_bounds = [min(latitude(valid)) max(latitude(valid))];
topo.field_longitude_bounds = [min(longitude(valid)) max(longitude(valid))];

end
