function topo = Build_Margin_Field_Topology(field_polygon,cell_size_m, ...
        metadata,source_feature_indices,margin_width_m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Build_Margin_Field_Topology
%%
%% A polygon grid for a field and a margin band around it: the polygon-grid
%% counterpart of Build_Cell_Topology's plot with its 5 m margin, so the
%% margin layouts of Run_Matrix_12 (runs 4 to 9) run on the same kind of
%% grid as the whole-field runs. T&C-SMEW-N v1.0 grids whole fields only.
%%
%% The simulated domain is the field polygon widened outward by
%% margin_width_m (polybuffer with mitred joins, so a field corner stays a
%% corner and the band keeps the field's edge directions). v1.0's majority
%% rule grids that domain exactly as it grids a field
%% (Build_Projected_Field_Topology with 'majority': rotation from the edges,
%% offset minimising the misrepresented area). A cell is then 'main' when
%% more than half of it lies inside the field polygon and 'buffer'
%% otherwise: the same majority test, applied to the field instead of the
%% domain.
%%
%% The band is resolved only by cells no wider than it. A wider cell
%% straddles field and margin, and the layouts would place rock and plants
%% by whichever side holds the larger share, so such a grid is refused
%% rather than silently mislabelled.
%%
%% Accounting follows v1.0. cell_overlap_area_m2 is each cell's overlap
%% with the DOMAIN polygon, so the domain results cover the field and its
%% margin, as the 12-column plot's do. cell_field_overlap_area_m2 is the
%% overlap with the field polygon alone, for field-only results.
%%
%%   field_polygon            projected field polygon (polyshape)
%%   cell_size_m              cell width [m] (Run_Topographic_Grid_Convergence)
%%   metadata, source_feature_indices
%%                            as for Build_Projected_Field_Topology
%%   margin_width_m           band width [m], > 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~(isa(field_polygon,'polyshape') && isscalar(field_polygon)) || ...
        isempty(field_polygon.Vertices) || area(field_polygon) <= 0
    error('Build_Margin_Field_Topology:EmptyField', ...
        'field_polygon must have finite positive area.');
end
if ~(isnumeric(margin_width_m) && isscalar(margin_width_m) && ...
        isfinite(margin_width_m) && margin_width_m > 0)
    error('Build_Margin_Field_Topology:BadMargin', ...
        'margin_width_m must be a finite positive width in metres.');
end
if isempty(cell_size_m)
    error('Build_Margin_Field_Topology:WidthRequired', ...
        ['A margin grid takes a prescribed cell width; set it from the ' ...
         'terrain (Run_Topographic_Grid_Convergence).']);
end
if ~(isnumeric(cell_size_m) && isscalar(cell_size_m) && ...
        isfinite(cell_size_m) && cell_size_m > 0)
    error('Build_Margin_Field_Topology:BadCellSize', ...
        'cell_size_m must be a finite positive scalar in metres.');
end
if cell_size_m > margin_width_m*(1+1e-9)
    error('Build_Margin_Field_Topology:MarginNotResolved', ...
        ['%.6g m cells cannot resolve a %.6g m margin: a cell would ' ...
         'straddle the field and its margin. Use cells no wider than ' ...
         'the margin.'],cell_size_m,margin_width_m);
end

domain_polygon = polybuffer(field_polygon,margin_width_m,'JointType','miter');
topo = Build_Projected_Field_Topology(domain_polygon,cell_size_m, ...
    metadata,source_feature_indices,'majority');

%% The field polygon in the grid's own frame (the inverse of
%% Grid_Local_To_Projected), where the squares are axis-aligned and their
%% overlaps are computed as Build_Polygon_Grid_Topology computes them.
angle = deg2rad(topo.grid_rotation_degrees);
origin = topo.projected_grid_origin_m;
field_vertices = field_polygon.Vertices;
x_offset = field_vertices(:,1)-origin(1);
y_offset = field_vertices(:,2)-origin(2);
field_local = polyshape(x_offset*cos(angle)+y_offset*sin(angle), ...
    -x_offset*sin(angle)+y_offset*cos(angle));

cell_field_overlap_area_m2 = zeros(topo.N_cells,1);
for c = 1:topo.N_cells
    square = polyshape( ...
        [topo.cell_x_lo(c) topo.cell_x_hi(c) topo.cell_x_hi(c) topo.cell_x_lo(c)], ...
        [topo.cell_y_lo(c) topo.cell_y_lo(c) topo.cell_y_hi(c) topo.cell_y_hi(c)]);
    cell_field_overlap_area_m2(c) = area(intersect(field_local,square));
end
area_tol = 1e-12*cell_size_m^2;
cell_field_overlap_area_m2(cell_field_overlap_area_m2 < area_tol) = 0;
full_cell = abs(cell_field_overlap_area_m2-topo.cell_area_m2) < area_tol;
cell_field_overlap_area_m2(full_cell) = topo.cell_area_m2(full_cell);

is_main = cell_field_overlap_area_m2 > 0.5*topo.cell_area_m2;
if ~any(is_main)
    error('Build_Margin_Field_Topology:NoFieldCell', ...
        'No %.6g m cell lies more than half inside the field polygon.', ...
        cell_size_m);
end
if all(is_main)
    error('Build_Margin_Field_Topology:NoMarginCell', ...
        'No %.6g m cell lies more than half inside the %.6g m margin.', ...
        cell_size_m,margin_width_m);
end
topo.cell_region(is_main) = {'main'};
topo.cell_region(~is_main) = {'buffer'};
topo.cell_subtype(is_main) = {'main'};
topo.cell_subtype(~is_main) = {'margin'};
topo.buffer_width_m = margin_width_m;
topo.margin_width_m = margin_width_m;
topo.margin_join = 'miter';
topo.cell_field_overlap_area_m2 = cell_field_overlap_area_m2;
topo.field_polygon_area_m2 = area(field_polygon);
topo.domain_polygon_area_m2 = area(domain_polygon);
[topo.field_polygon_boundary_x_m,topo.field_polygon_boundary_y_m] = ...
    boundary(field_local);

end
