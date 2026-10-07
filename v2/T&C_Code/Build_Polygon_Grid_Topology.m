function topo = Build_Polygon_Grid_Topology(field_polygon,cell_size_m, ...
        inclusion,grid_offset_m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Build_Polygon_Grid_Topology
%%
%% Construct equal full square cells over a metric field polygon.
%% inclusion 'any' (default): every square with positive polygon overlap is
%% retained and simulated in full.  'majority': only the squares more than
%% half inside the polygon; the field land in the dropped squares is left
%% unsimulated and recorded as field_missed_area_m2.  grid_offset_m [dx dy]
%% in [0,h) shifts the lattice's lower-left corner below and left of the
%% polygon's bounding box (default [0 0]; Place_Field_Polygon_Grid).
%% Polygon-cell intersections are stored separately for field accounting.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3 || isempty(inclusion); inclusion = 'any'; end
if nargin < 4 || isempty(grid_offset_m); grid_offset_m = [0 0]; end
if ~(isa(field_polygon,'polyshape') && isscalar(field_polygon)) || ...
        isempty(field_polygon.Vertices) || area(field_polygon) <= 0
    error('Build_Polygon_Grid_Topology:EmptyField', ...
        'field_polygon must have finite positive area.');
end
if ~(isnumeric(cell_size_m) && isscalar(cell_size_m) && ...
        isfinite(cell_size_m) && cell_size_m > 0)
    error('Build_Polygon_Grid_Topology:BadCellSize', ...
        'cell_size_m must be a finite positive scalar in metres.');
end
if ~((ischar(inclusion) || (isstring(inclusion) && isscalar(inclusion))) && ...
        any(strcmp(char(inclusion),{'any','majority'})))
    error('Build_Polygon_Grid_Topology:BadInclusion', ...
        'inclusion must be ''any'' or ''majority''.');
end
majority = strcmp(char(inclusion),'majority');
if ~(isnumeric(grid_offset_m) && numel(grid_offset_m) == 2 && ...
        all(isfinite(grid_offset_m)) && all(grid_offset_m >= 0) && ...
        all(grid_offset_m < cell_size_m))
    error('Build_Polygon_Grid_Topology:BadOffset', ...
        'grid_offset_m must be two offsets in [0,cell_size_m) metres.');
end
grid_offset_m = reshape(double(grid_offset_m),1,2);

vertices = field_polygon.Vertices;
finite_vertex = all(isfinite(vertices),2);
separator = all(isnan(vertices),2);
if any(~finite_vertex & ~separator)
    error('Build_Polygon_Grid_Topology:BadField', ...
        'field_polygon vertices must be finite metric coordinates.');
end
vertices = vertices(finite_vertex,:);

x_min = min(vertices(:,1))-grid_offset_m(1);
x_max = max(vertices(:,1));
y_min = min(vertices(:,2))-grid_offset_m(2);
y_max = max(vertices(:,2));
nx = Grid_Cell_Count(x_max-x_min,cell_size_m);
ny = Grid_Cell_Count(y_max-y_min,cell_size_m);
N_slots = nx*ny;

cell_row = zeros(N_slots,1);
cell_col = zeros(N_slots,1);
cell_x_lo = zeros(N_slots,1);
cell_x_hi = zeros(N_slots,1);
cell_y_lo = zeros(N_slots,1);
cell_y_hi = zeros(N_slots,1);
cell_overlap_area_m2 = zeros(N_slots,1);

for iy = 1:ny
    for ix = 1:nx
        c = (iy-1)*nx+ix;
        cell_row(c) = iy;
        cell_col(c) = ix;
        cell_x_lo(c) = x_min+(ix-1)*cell_size_m;
        cell_x_hi(c) = cell_x_lo(c)+cell_size_m;
        cell_y_lo(c) = y_min+(iy-1)*cell_size_m;
        cell_y_hi(c) = cell_y_lo(c)+cell_size_m;
        square = polyshape( ...
            [cell_x_lo(c) cell_x_hi(c) cell_x_hi(c) cell_x_lo(c)], ...
            [cell_y_lo(c) cell_y_lo(c) cell_y_hi(c) cell_y_hi(c)]);
        cell_overlap_area_m2(c) = area(intersect(field_polygon,square));
    end
end

area_tol = 1e-12*cell_size_m^2;
cell_overlap_area_m2(cell_overlap_area_m2 < area_tol) = 0;
retain = cell_overlap_area_m2 > 0;
if ~any(retain)
    error('Build_Polygon_Grid_Topology:EmptyGrid', ...
        'No full grid cells intersect the field polygon.');
end
field_area_m2 = area(field_polygon);
overlap_tol = max(1e-6,1e-7*field_area_m2);
if majority
    %% The intersections of every touching square, kept or not, partition
    %% the field; what the dropped squares hold is left unsimulated.
    overlap_area_error_m2 = sum(cell_overlap_area_m2)-field_area_m2;
    if abs(overlap_area_error_m2) >= overlap_tol
        error('Build_Polygon_Grid_Topology:AreaMismatch', ...
            ['Polygon-cell intersections differ from the field area by %.17g m2 ' ...
             'for %.17g m cells.'],overlap_area_error_m2,cell_size_m);
    end
    retain = cell_overlap_area_m2 > 0.5*cell_size_m^2;
    if ~any(retain)
        error('Build_Polygon_Grid_Topology:EmptyGrid', ...
            'No %.6g m grid cell lies more than half inside the field polygon.', ...
            cell_size_m);
    end
end
cell_grid_slot = find(retain);
cell_row = cell_row(retain);
cell_col = cell_col(retain);
cell_x_lo = cell_x_lo(retain);
cell_x_hi = cell_x_hi(retain);
cell_y_lo = cell_y_lo(retain);
cell_y_hi = cell_y_hi(retain);
cell_overlap_area_m2 = cell_overlap_area_m2(retain);
N_cells = numel(cell_grid_slot);
cell_id = (1:N_cells)';
cell_name = arrayfun(@(c) sprintf('Cell_Grid_%04d',c),cell_id, ...
    'UniformOutput',false);
cell_area_m2 = repmat(cell_size_m^2,N_cells,1);
full_cell = abs(cell_overlap_area_m2-cell_area_m2) < area_tol;
cell_overlap_area_m2(full_cell) = cell_area_m2(full_cell);
slot_to_cell = zeros(N_slots,1);
slot_to_cell(cell_grid_slot) = cell_id;
edge_a = zeros(2*N_cells,1);
edge_b = zeros(2*N_cells,1);
e = 0;
for c = 1:N_cells
    slot = cell_grid_slot(c);
    if cell_col(c) < nx && slot_to_cell(slot+1) > 0
        e = e+1;
        edge_a(e) = c;
        edge_b(e) = slot_to_cell(slot+1);
    end
    if cell_row(c) < ny && slot_to_cell(slot+nx) > 0
        e = e+1;
        edge_a(e) = c;
        edge_b(e) = slot_to_cell(slot+nx);
    end
end
N_edges = e;
edge_a = edge_a(1:N_edges);
edge_b = edge_b(1:N_edges);

edge_face_width_m = repmat(cell_size_m,N_edges,1);
edge_face_distance_m = repmat(cell_size_m,N_edges,1);
edge_face_distance_a_m = repmat(cell_size_m/2,N_edges,1);
edge_face_distance_b_m = repmat(cell_size_m/2,N_edges,1);
edge_active = ones(N_edges,1);
adj = sparse([edge_a;edge_b],[edge_b;edge_a],1,N_cells,N_cells) > 0;

outer_face_width_m = zeros(N_cells,1);
for c = 1:N_cells
    slot = cell_grid_slot(c);
    if cell_col(c) == 1 || slot_to_cell(slot-1) == 0
        outer_face_width_m(c) = outer_face_width_m(c)+cell_size_m;
    end
    if cell_col(c) == nx || slot_to_cell(slot+1) == 0
        outer_face_width_m(c) = outer_face_width_m(c)+cell_size_m;
    end
    if cell_row(c) == 1 || slot_to_cell(slot-nx) == 0
        outer_face_width_m(c) = outer_face_width_m(c)+cell_size_m;
    end
    if cell_row(c) == ny || slot_to_cell(slot+nx) == 0
        outer_face_width_m(c) = outer_face_width_m(c)+cell_size_m;
    end
end

simulated_area_m2 = sum(cell_area_m2);
grid_area_tol = max(1e-9,1e-13*field_area_m2);
if (~majority && simulated_area_m2 < field_area_m2) || ...
        any(cell_overlap_area_m2 > cell_area_m2+grid_area_tol)
    error('Build_Polygon_Grid_Topology:AreaMismatch', ...
        'Retained full-cell area does not contain the complete field polygon.');
end
if majority
    overlap_area_error_m2 = 0;
    field_missed_area_m2 = max(0,field_area_m2-sum(cell_overlap_area_m2));
else
    overlap_area_error_m2 = sum(cell_overlap_area_m2)-field_area_m2;
end
if abs(overlap_area_error_m2) >= overlap_tol
    error('Build_Polygon_Grid_Topology:AreaMismatch', ...
        ['Polygon-cell intersections differ from the field area by %.17g m2 ' ...
         'for %.17g m cells.'],overlap_area_error_m2,cell_size_m);
end
%% Polyshape booleans can leave sub-ppm slivers at coincident grid edges.
%% Close only that bounded numerical residual so field accounting is exact.
if overlap_area_error_m2 < 0
    [capacity,correction_cell] = max(cell_area_m2-cell_overlap_area_m2);
    if capacity < -overlap_area_error_m2
        error('Build_Polygon_Grid_Topology:AreaMismatch', ...
            'Numerical overlap residual exceeds the retained-cell capacity.');
    end
    cell_overlap_area_m2(correction_cell) = ...
        cell_overlap_area_m2(correction_cell)-overlap_area_error_m2;
elseif overlap_area_error_m2 > 0
    [available,correction_cell] = max(cell_overlap_area_m2);
    if available < overlap_area_error_m2
        error('Build_Polygon_Grid_Topology:AreaMismatch', ...
            'Numerical overlap residual exceeds the polygon overlap area.');
    end
    cell_overlap_area_m2(correction_cell) = ...
        cell_overlap_area_m2(correction_cell)-overlap_area_error_m2;
end
cell_overlap_fraction = cell_overlap_area_m2./cell_area_m2;
expected_outer_width_m = (4*N_cells-2*N_edges)*cell_size_m;
perimeter_tol = max(1e-9,1e-12*expected_outer_width_m);
assert(abs(sum(outer_face_width_m)-expected_outer_width_m) < perimeter_tol, ...
    'Build_Polygon_Grid_Topology:PerimeterMismatch', ...
    'Grid outer-face widths do not recover the retained-cell perimeter.');

[field_boundary_x_m,field_boundary_y_m] = boundary(field_polygon);
topo.field_side_m = [];
topo.buffer_width_m = 0;
topo.cell_size_m = cell_size_m;
topo.main_grid_nx = nx;
topo.main_grid_ny = ny;
topo.N_cells = N_cells;
topo.grid_slot_count = N_slots;
topo.cell_id = cell_id;
topo.cell_name = cell_name;
topo.cell_grid_slot = cell_grid_slot;
topo.cell_row = cell_row;
topo.cell_col = cell_col;
topo.cell_region = repmat({'main'},N_cells,1);
topo.cell_subtype = repmat({'main'},N_cells,1);
topo.cell_x_lo = cell_x_lo;
topo.cell_x_hi = cell_x_hi;
topo.cell_y_lo = cell_y_lo;
topo.cell_y_hi = cell_y_hi;
topo.cell_xc = 0.5*(cell_x_lo+cell_x_hi);
topo.cell_yc = 0.5*(cell_y_lo+cell_y_hi);
topo.cell_dx = repmat(cell_size_m,N_cells,1);
topo.cell_dy = repmat(cell_size_m,N_cells,1);
topo.cell_area_m2 = cell_area_m2;
topo.cell_overlap_area_m2 = cell_overlap_area_m2;
topo.cell_overlap_fraction = cell_overlap_fraction;
topo.field_area_m2 = field_area_m2;
topo.simulated_area_m2 = simulated_area_m2;
topo.field_boundary_x_m = field_boundary_x_m;
topo.field_boundary_y_m = field_boundary_y_m;
topo.field_region_count = numel(regions(field_polygon));
topo.N_edges = N_edges;
topo.edge_a = edge_a;
topo.edge_b = edge_b;
topo.edge_face_width_m = edge_face_width_m;
topo.edge_face_distance_m = edge_face_distance_m;
topo.edge_face_distance_a_m = edge_face_distance_a_m;
topo.edge_face_distance_b_m = edge_face_distance_b_m;
topo.edge_active = edge_active;
topo.adj = adj;
topo.outer_face_width_m = outer_face_width_m;
if majority
    topo.cell_inclusion = 'majority';
    topo.grid_offset_m = grid_offset_m;
    topo.field_missed_area_m2 = field_missed_area_m2;
end

end
