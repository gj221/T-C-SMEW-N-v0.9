function [oriented_polygon,geometry] = Orient_Field_Polygon_Grid( ...
        field_polygon,cell_size_m,rotation_degrees)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Align a square grid to minimize padded area outside a field polygon.
%% rotation_degrees (optional): rotate the grid by this angle instead
%% (Place_Field_Polygon_Grid).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3; rotation_degrees = []; end
if ~(isa(field_polygon,'polyshape') && isscalar(field_polygon)) || ...
        isempty(field_polygon.Vertices) || area(field_polygon) <= 0
    error('Orient_Field_Polygon_Grid:EmptyField', ...
        'field_polygon must have finite positive area.');
end
if ~isempty(rotation_degrees) && ~(isnumeric(rotation_degrees) && ...
        isscalar(rotation_degrees) && isfinite(rotation_degrees))
    error('Orient_Field_Polygon_Grid:BadRotation', ...
        'rotation_degrees must be a finite scalar.');
end

vertices = field_polygon.Vertices;
finite_vertex = all(isfinite(vertices),2);
separator = all(isnan(vertices),2);
if any(~finite_vertex & ~separator)
    error('Orient_Field_Polygon_Grid:BadField', ...
        'field_polygon vertices must be finite metric coordinates.');
end
finite_vertices = vertices(finite_vertex,:);
if size(finite_vertices,1) < 3
    error('Orient_Field_Polygon_Grid:BadField', ...
        'field_polygon must contain at least three finite vertices.');
end

hull_index = convhull(finite_vertices(:,1),finite_vertices(:,2));
hull_vertices = finite_vertices(hull_index,:);
envelope = polyshape(hull_vertices(:,1),hull_vertices(:,2));
envelope_area_m2 = area(envelope);
envelope_perimeter_m = perimeter(envelope);
if ~(isfinite(envelope_area_m2) && envelope_area_m2 > 0 && ...
        isfinite(envelope_perimeter_m) && envelope_perimeter_m > 0)
    error('Orient_Field_Polygon_Grid:BadEnvelope', ...
        'The polygon convex envelope must have finite positive area and perimeter.');
end

if isempty(cell_size_m)
    cell_size_m = 2*envelope_area_m2/envelope_perimeter_m;
elseif ~(isnumeric(cell_size_m) && isscalar(cell_size_m) && ...
        isfinite(cell_size_m) && cell_size_m > 0)
    error('Orient_Field_Polygon_Grid:BadCellSize', ...
        'cell_size_m must be a finite positive scalar in metres.');
end

origin_m = mean(hull_vertices(1:end-1,:),1);
hull_centered = hull_vertices-origin_m;
edge = diff(hull_vertices,1,1);
sector_angle = mod(atan2(edge(:,2),edge(:,1)),pi/2);
sector_angle = unique(round([0;sector_angle;pi/2],12));
candidate_angle = sector_angle;
for k = 1:numel(sector_angle)-1
    angle_lo = sector_angle(k);
    angle_hi = sector_angle(k+1);
    if angle_hi-angle_lo <= 1e-12; continue; end
    angle_mid = 0.5*(angle_lo+angle_hi);
    c = cos(angle_mid);
    s = sin(angle_mid);
    u = hull_centered(:,1)*c+hull_centered(:,2)*s;
    v = -hull_centered(:,1)*s+hull_centered(:,2)*c;
    [~,u_max] = max(u); [~,u_min] = min(u);
    [~,v_max] = max(v); [~,v_min] = min(v);
    width_vector = hull_centered(u_max,:)-hull_centered(u_min,:);
    height_vector = hull_centered(v_max,:)-hull_centered(v_min,:);
    width_A = width_vector(1);
    width_B = width_vector(2);
    height_A = height_vector(2);
    height_B = -height_vector(1);
    candidate_angle = [candidate_angle; ...
        grid_count_transition_angles( ...
        width_A,width_B,cell_size_m,angle_lo,angle_hi); ...
        grid_count_transition_angles( ...
        height_A,height_B,cell_size_m,angle_lo,angle_hi); ...
        bounding_area_stationary_angles( ...
        width_A,width_B,height_A,height_B,angle_lo,angle_hi)]; %#ok<AGROW>
end
candidate_angle = unique(round(candidate_angle,12));
candidate_angle = unique([candidate_angle; ...
    0.5*(candidate_angle(1:end-1)+candidate_angle(2:end))]);
if ~isempty(rotation_degrees)
    candidate_angle = deg2rad(double(rotation_degrees));
end

best_grid_area = Inf;
best_bounds_area = Inf;
best_angle = [];
best_nx = [];
best_ny = [];
best_width = [];
best_height = [];
for angle = reshape(candidate_angle,1,[])
    c = cos(angle);
    s = sin(angle);
    u = hull_centered(:,1)*c+hull_centered(:,2)*s;
    v = -hull_centered(:,1)*s+hull_centered(:,2)*c;
    width_m = max(u)-min(u);
    height_m = max(v)-min(v);
    nx = Grid_Cell_Count(width_m,cell_size_m);
    ny = Grid_Cell_Count(height_m,cell_size_m);
    grid_area_m2 = nx*ny*cell_size_m^2;
    bounds_area_m2 = width_m*height_m;
    better_grid = grid_area_m2 < best_grid_area-1e-9*cell_size_m^2;
    same_grid = abs(grid_area_m2-best_grid_area) <= 1e-9*cell_size_m^2;
    better_bounds = bounds_area_m2 < best_bounds_area-1e-9*envelope_area_m2;
    if better_grid || (same_grid && better_bounds)
        best_grid_area = grid_area_m2;
        best_bounds_area = bounds_area_m2;
        best_angle = angle;
        best_nx = nx;
        best_ny = ny;
        best_width = width_m;
        best_height = height_m;
    end
end

c = cos(best_angle);
s = sin(best_angle);
centered = vertices-origin_m;
u = centered(:,1)*c+centered(:,2)*s;
v = -centered(:,1)*s+centered(:,2)*c;
oriented_polygon = polyshape(u,v,'Simplify',false);
area_tol = max(1e-6,1e-9*area(field_polygon));
assert(abs(area(oriented_polygon)-area(field_polygon)) < area_tol, ...
    'Orient_Field_Polygon_Grid:AreaMismatch', ...
    'Rigid grid rotation did not preserve polygon area.');

geometry.cell_size_m = cell_size_m;
geometry.rotation_degrees = rad2deg(best_angle);
geometry.projected_grid_origin_m = origin_m;
geometry.grid_nx = best_nx;
geometry.grid_ny = best_ny;
geometry.grid_cell_count = best_nx*best_ny;
geometry.bounding_width_m = best_width;
geometry.bounding_height_m = best_height;
geometry.bounding_area_m2 = best_bounds_area;
geometry.envelope_area_m2 = envelope_area_m2;
geometry.envelope_perimeter_m = envelope_perimeter_m;
geometry.simulated_area_m2 = best_grid_area;
geometry.field_area_m2 = area(field_polygon);
geometry.outside_grid_area_m2 = best_grid_area-geometry.field_area_m2;
geometry.outside_grid_fraction = geometry.outside_grid_area_m2/best_grid_area;
geometry.orientation_objective = 'minimum_full_grid_area';
if ~isempty(rotation_degrees)
    geometry.orientation_objective = 'prescribed_rotation';
end
geometry.orientation_candidate_count = numel(candidate_angle);

end

function angles = grid_count_transition_angles(A,B,cell_size_m,angle_lo,angle_hi)
radius_m = hypot(A,B);
angles = zeros(0,1);
if radius_m <= 0; return; end

maximum_count = floor(radius_m/cell_size_m+1e-10);
phase = atan2(B,A);
for count = 1:maximum_count
    ratio = min(1,count*cell_size_m/radius_m);
    offset = acos(ratio);
    roots = [phase-offset phase+offset];
    for period = -2:2
        shifted = roots+period*2*pi;
        keep = shifted >= angle_lo-1e-12 & shifted <= angle_hi+1e-12;
        angles = [angles;shifted(keep).']; %#ok<AGROW>
    end
end
angles = min(angle_hi,max(angle_lo,angles));
end

function angles = bounding_area_stationary_angles( ...
        width_A,width_B,height_A,height_B,angle_lo,angle_hi)
numerator = width_A*height_B+width_B*height_A;
denominator = width_A*height_A-width_B*height_B;
base = 0.5*atan2(numerator,denominator);
roots = base+(-2:2)*(pi/2);
keep = roots >= angle_lo-1e-12 & roots <= angle_hi+1e-12;
angles = roots(keep).';
angles = min(angle_hi,max(angle_lo,angles));
end
