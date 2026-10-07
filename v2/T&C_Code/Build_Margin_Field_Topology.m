function topo = Build_Margin_Field_Topology(field_polygon,cell_size_m, ...
        metadata,source_feature_indices,margin_width_m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Build_Margin_Field_Topology
%%
%% A polygon grid whose outer band is the field margin: the polygon-grid
%% counterpart of Build_Cell_Topology, whose 60 m plot carries its margin
%% as the outer 5 m ring (2500 m2 main inside 1100 m2 of margin), so the
%% margin layouts of Run_Matrix_12 (runs 4 to 9) run on the same grid as
%% the whole-field runs. T&C-SMEW-N v1.0 grids whole fields only, every
%% cell 'main'.
%%
%% The grid is v1.0's, unchanged: the field polygon under the majority
%% rule (Build_Projected_Field_Topology with 'majority': rotation from the
%% edges, offset minimising the misrepresented area), so every layout
%% simulates the same squares. Only the labels differ. The interior is the
%% field shrunk inward by margin_width_m (polybuffer with mitred joins, so
%% a corner stays a corner); a square is 'main' when more than half of it
%% lies in the interior and 'buffer' (the margin) otherwise -- the same
%% majority test, applied to the interior instead of the field.
%%
%% The band is resolved only by squares no wider than it. A wider square
%% straddles interior and margin, and the layouts would place rock and
%% plants by whichever side holds the larger share, so such a grid is
%% refused rather than silently mislabelled.
%%
%% Accounting is v1.0's: cell_overlap_area_m2 is each square's overlap
%% with the field, so domain results cover the whole field, margin
%% included, as the plot's do. cell_interior_overlap_area_m2 is the
%% overlap with the interior alone.
%%
%%   field_polygon            projected field polygon (polyshape)
%%   cell_size_m              square width [m] (Run_Topographic_Grid_Convergence)
%%   metadata, source_feature_indices
%%                            as for Build_Projected_Field_Topology
%%   margin_width_m           margin width [m], > 0
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
        ['A margin grid takes a prescribed square width; set it from the ' ...
         'terrain (Run_Topographic_Grid_Convergence).']);
end
if ~(isnumeric(cell_size_m) && isscalar(cell_size_m) && ...
        isfinite(cell_size_m) && cell_size_m > 0)
    error('Build_Margin_Field_Topology:BadCellSize', ...
        'cell_size_m must be a finite positive scalar in metres.');
end
if cell_size_m > margin_width_m*(1+1e-9)
    error('Build_Margin_Field_Topology:MarginNotResolved', ...
        ['%.6g m squares cannot resolve a %.6g m margin: a square would ' ...
         'straddle the margin and the interior. Use squares no wider ' ...
         'than the margin.'],cell_size_m,margin_width_m);
end
interior_polygon = polybuffer(field_polygon,-margin_width_m,'JointType','miter');
if interior_polygon.NumRegions == 0 || area(interior_polygon) <= 0
    error('Build_Margin_Field_Topology:FieldNarrowerThanMargin', ...
        'The field has no interior inside a %.6g m margin.',margin_width_m);
end

topo = Build_Projected_Field_Topology(field_polygon,cell_size_m, ...
    metadata,source_feature_indices,'majority');

%% The interior in the grid's own frame (the inverse of
%% Grid_Local_To_Projected), where the squares are axis-aligned and their
%% overlaps are computed as Build_Polygon_Grid_Topology computes them.
angle = deg2rad(topo.grid_rotation_degrees);
origin = topo.projected_grid_origin_m;
interior_vertices = interior_polygon.Vertices;
x_offset = interior_vertices(:,1)-origin(1);
y_offset = interior_vertices(:,2)-origin(2);
interior_local = polyshape(x_offset*cos(angle)+y_offset*sin(angle), ...
    -x_offset*sin(angle)+y_offset*cos(angle));

cell_interior_overlap_area_m2 = zeros(topo.N_cells,1);
for c = 1:topo.N_cells
    square = polyshape( ...
        [topo.cell_x_lo(c) topo.cell_x_hi(c) topo.cell_x_hi(c) topo.cell_x_lo(c)], ...
        [topo.cell_y_lo(c) topo.cell_y_lo(c) topo.cell_y_hi(c) topo.cell_y_hi(c)]);
    cell_interior_overlap_area_m2(c) = area(intersect(interior_local,square));
end
area_tol = 1e-12*cell_size_m^2;
cell_interior_overlap_area_m2(cell_interior_overlap_area_m2 < area_tol) = 0;
full_cell = abs(cell_interior_overlap_area_m2-topo.cell_area_m2) < area_tol;
cell_interior_overlap_area_m2(full_cell) = topo.cell_area_m2(full_cell);

is_main = cell_interior_overlap_area_m2 > 0.5*topo.cell_area_m2;
if ~any(is_main)
    error('Build_Margin_Field_Topology:NoInteriorCell', ...
        'No %.6g m square lies more than half inside the field''s interior.', ...
        cell_size_m);
end
if all(is_main)
    error('Build_Margin_Field_Topology:NoMarginCell', ...
        'No %.6g m square lies more than half inside the %.6g m margin.', ...
        cell_size_m,margin_width_m);
end
topo.cell_region(is_main) = {'main'};
topo.cell_region(~is_main) = {'buffer'};
topo.cell_subtype(is_main) = {'main'};
topo.cell_subtype(~is_main) = {'margin'};
topo.buffer_width_m = margin_width_m;
topo.margin_width_m = margin_width_m;
topo.margin_join = 'miter';
topo.cell_interior_overlap_area_m2 = cell_interior_overlap_area_m2;
topo.interior_area_m2 = area(interior_polygon);
[topo.interior_boundary_x_m,topo.interior_boundary_y_m] = ...
    boundary(interior_local);

end
