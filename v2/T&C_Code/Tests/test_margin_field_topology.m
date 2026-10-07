function test_margin_field_topology()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The margin of a field grid (Build_Margin_Field_Topology), which puts
%% Run_Matrix_12's margin layouts on v1.0's polygon grid.
%%   1  the 60 m plot on 5 m squares with a 5 m margin: v1.0's 144 squares,
%%      unchanged, 100 interior squares main and the 44 of the outer ring
%%      margin, 2500 and 1100 m2 exactly as the 12-column plot has them
%%   2  the same plot turned 37 degrees in the projected frame
%%   3  an irregular field on 2.5 m squares: v1.0's squares, main squares
%%      more than half inside the interior, margin squares at most half,
%%      the whole interior simulated
%%   4  squares wider than the margin, a zero margin, a missing width, a
%%      field with no interior and a margin under the coverage rule are
%%      refused
%% Building a projected grid needs the Mapping Toolbox (projcrs); without
%% it the test reports a skip. A campaign loads a grid built where the
%% toolbox exists (Build_Field_Grid_Topology, FieldTopologyFile).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(ver('map')) || ~license('test','MAP_Toolbox')
    fprintf('SKIP: test_margin_field_topology needs the Mapping Toolbox.\n');
    return
end
test_directory = fileparts(mfilename('fullpath'));
addpath(fileparts(test_directory));

origin = [397962.6 4435380.4];          %% Energy Farm, EPSG 26916
meta = synthetic_metadata(26916);

%% --- 1 --------------------------------------------------------------------
plot60 = polyshape(origin(1)+[0 60 60 0],origin(2)+[0 0 60 60]);
topo = Build_Margin_Field_Topology(plot60,5,meta,1,5);
v10 = Build_Projected_Field_Topology(plot60,5,meta,1,'majority');
assert_same_squares(topo,v10);
main = strcmp(topo.cell_region,'main');
assert(topo.N_cells == 144 && sum(main) == 100, ...
    'The plot holds 144 squares, 100 of them main (%d, %d).', ...
    topo.N_cells,sum(main));
assert(abs(sum(topo.cell_overlap_area_m2(main))-2500) < 1e-6 && ...
    abs(sum(topo.cell_overlap_area_m2(~main))-1100) < 1e-6, ...
    'Interior and margin are %.9g and %.9g m2, not 2500 and 1100.', ...
    sum(topo.cell_overlap_area_m2(main)),sum(topo.cell_overlap_area_m2(~main)));
plot12 = Build_Cell_Topology(60,5,2,2);
plot12_main = strcmp(plot12.cell_region,'main');
assert(abs(sum(plot12.cell_area_m2(plot12_main))-2500) < 1e-9 && ...
    abs(sum(plot12.cell_area_m2(~plot12_main))-1100) < 1e-9, ...
    'The 12-column plot it mirrors keeps 2500 m2 inside 1100 m2 of margin.');
assert(all(abs(topo.cell_interior_overlap_area_m2(main)-25) < 1e-6) && ...
    all(topo.cell_interior_overlap_area_m2(~main) < 1e-6), ...
    'Main squares lie wholly in the interior and margin squares wholly out.');
assert(abs(topo.interior_area_m2-2500) < 1e-6 && ...
    topo.buffer_width_m == 5 && topo.margin_width_m == 5);
assert(all(strcmp(topo.cell_subtype(~main),'margin')) && ...
    all(strcmp(topo.cell_subtype(main),'main')));
assert(all(topo.outer_face_width_m(~main) > 0) && ...
    all(topo.outer_face_width_m(main) == 0), ...
    'A 5 m margin on 5 m squares is the outer ring of the field.');
[~,weights,~,field_area] = Field_Aggregation_Areas(topo);
assert(abs(sum(weights)-1) < 1e-12 && ...
    abs(sum(weights(main))-2500/3600) < 1e-9 && ...
    abs(field_area-3600) < 1e-6, ...
    'Domain results weigh interior and margin by area, as the plot does.');

%% --- 2 --------------------------------------------------------------------
turn = deg2rad(37);
rotation = [cos(turn) -sin(turn); sin(turn) cos(turn)];
corners = [0 0; 60 0; 60 60; 0 60]*rotation.'+origin;
turned = Build_Margin_Field_Topology( ...
    polyshape(corners(:,1),corners(:,2)),5,meta,1,5);
assert(abs(mod(turned.grid_rotation_degrees,90)-37) < 1e-6 && ...
    turned.N_cells == 144 && sum(strcmp(turned.cell_region,'main')) == 100, ...
    'The turned plot gives %d squares at %.6f degrees.', ...
    turned.N_cells,turned.grid_rotation_degrees);

%% --- 3 --------------------------------------------------------------------
outline = [0 0; 83 0; 95 26; 76 58; 34 62; 8 44; -6 18]+origin;
field = polyshape(outline(:,1),outline(:,2));
irregular = Build_Margin_Field_Topology(field,2.5,meta,1,5);
assert_same_squares(irregular,Build_Projected_Field_Topology(field,2.5,meta,1,'majority'));
main = strcmp(irregular.cell_region,'main');
half = 0.5*irregular.cell_area_m2;
assert(any(main) && any(~main) && ...
    all(irregular.cell_interior_overlap_area_m2(main) > half(main)) && ...
    all(irregular.cell_interior_overlap_area_m2(~main) <= half(~main)+1e-9), ...
    'Squares are labelled by the majority of their area in the interior.');
assert(all(irregular.cell_interior_overlap_area_m2 <= ...
    irregular.cell_overlap_area_m2+1e-6), ...
    'The interior lies inside the field.');
interior = polybuffer(field,-5,'JointType','miter');
assert(abs(sum(irregular.cell_interior_overlap_area_m2)-area(interior)) < ...
    1e-6*area(interior) && abs(irregular.interior_area_m2-area(interior)) < 1e-6, ...
    ['Squares the majority rule drops lie on the field edge, so every ' ...
     'square metre of the interior is simulated.']);

%% --- 4 --------------------------------------------------------------------
expect_error(@() Build_Margin_Field_Topology(plot60,10,meta,1,5), ...
    'Build_Margin_Field_Topology:MarginNotResolved');
expect_error(@() Build_Margin_Field_Topology(plot60,5,meta,1,0), ...
    'Build_Margin_Field_Topology:BadMargin');
expect_error(@() Build_Margin_Field_Topology(plot60,[],meta,1,5), ...
    'Build_Margin_Field_Topology:WidthRequired');
strip = polyshape(origin(1)+[0 60 60 0],origin(2)+[0 0 8 8]);
expect_error(@() Build_Margin_Field_Topology(strip,2,meta,1,5), ...
    'Build_Margin_Field_Topology:FieldNarrowerThanMargin');
expect_error(@() Build_Shapefile_Grid_Setups('none.shp',5,[],0.95,5), ...
    'Build_Shapefile_Grid_Setups:MarginNeedsMajority');

fprintf('test_margin_field_topology: all checks passed\n');

end

function assert_same_squares(topo,v10)
%% The margin only relabels v1.0's grid; the squares, their overlaps with
%% the field and their graph are v1.0's.
for name = {'N_cells','cell_x_lo','cell_x_hi','cell_y_lo','cell_y_hi', ...
        'cell_area_m2','cell_overlap_area_m2','edge_a','edge_b', ...
        'outer_face_width_m','grid_rotation_degrees','projected_grid_origin_m'}
    assert(isequal(topo.(name{1}),v10.(name{1})), ...
        'The margin changed v1.0''s %s.',name{1});
end
end

function meta = synthetic_metadata(epsg)
%% What Read_Shapefile_Polygons returns for one feature already in the
%% projected CRS (Build_Projected_Field_Topology reads these fields).
meta = struct('source_shapefile','synthetic', ...
    'source_crs_name','NAD83 / UTM zone 16N', ...
    'projected_crs_name','NAD83 / UTM zone 16N', ...
    'projected_crs_epsg',epsg,'shapefile_feature_count',1, ...
    'feature_indices',1,'source_finite_vertex_count',4, ...
    'projected_polygon_vertex_count',4, ...
    'projection_roundtrip_max_error_m',0,'source_to_polygon_max_error_m',0);
end

function expect_error(fcn,identifier)
try
    fcn();
catch ME
    assert(strcmp(ME.identifier,identifier), ...
        'Expected %s, got %s.',identifier,ME.identifier);
    return
end
error('Expected error %s was not raised.',identifier);
end
