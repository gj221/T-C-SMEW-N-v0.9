function test_field_grid_majority()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The majority cell rule and its lattice placement.
%%   1  a 10.3 m x 7.2 m rectangle on 1 m cells: the lattice puts two sides
%%      on grid lines, keeps the 70 full squares and drops the 0.3 m column
%%      and the 0.2 m row, so the misrepresented area is exactly the field
%%      land left out, 0.3*7 + 0.2*10 + 0.3*0.2 = 4.16 m2
%%   2  the same rectangle turned 23 degrees in the projected frame: the
%%      grid turns with it and finds the same lattice
%%   3  an irregular field with a hole, turned 31 degrees: the boundary
%%      integral of the placement equals the polygon clipping of the
%%      topology, the touching squares partition the field, every kept
%%      square is more than half inside, and field totals cover the whole
%%      field
%%   4  the default rule is untouched: no offset, every touching square
%%   5  bad rules, widths and offsets are refused
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_directory = fileparts(mfilename('fullpath'));
code_directory = fileparts(test_directory);
addpath(code_directory);

%% --- 1 --------------------------------------------------------------------
rectangle = polyshape([0 10.3 10.3 0],[0 0 7.2 7.2]);
[grid_polygon,geometry,offset] = Place_Field_Polygon_Grid(rectangle,1);
topo = Build_Polygon_Grid_Topology(grid_polygon,1,'majority',offset);
missed = 0.3*7+0.2*10+0.3*0.2;
assert(topo.N_cells == 70 && all(abs(topo.cell_overlap_fraction-1) < 1e-9), ...
    'The rectangle keeps its 70 full squares (%d kept).',topo.N_cells);
assert(abs(topo.field_missed_area_m2-missed) < 1e-9 && ...
    abs(min(geometry.candidate_misrepresented_area_m2)-missed) < 1e-9, ...
    'The rectangle misrepresents %.12g m2, not %.12g.', ...
    topo.field_missed_area_m2,missed);
assert(strcmp(topo.cell_inclusion,'majority') && isequal(topo.grid_offset_m,offset));

%% --- 2 --------------------------------------------------------------------
turn = deg2rad(23);
R = [cos(turn) -sin(turn); sin(turn) cos(turn)];
corners = [0 0; 10.3 0; 10.3 7.2; 0 7.2]*R.'+[512345.6 4598765.4];
turned = polyshape(corners(:,1),corners(:,2));
[turned_polygon,turned_geometry,turned_offset] = Place_Field_Polygon_Grid(turned,1);
turned_topo = Build_Polygon_Grid_Topology(turned_polygon,1,'majority',turned_offset);
assert(abs(mod(turned_geometry.rotation_degrees,90)-23) < 1e-9 && ...
    turned_topo.N_cells == 70 && abs(turned_topo.field_missed_area_m2-missed) < 1e-6, ...
    'The turned rectangle gives %d squares at %.6f degrees.', ...
    turned_topo.N_cells,turned_geometry.rotation_degrees);

%% --- 3 --------------------------------------------------------------------
outer = [0 0; 41 0; 47 13; 38 29; 17 31; 4 22; -3 9];
hole = [15 10; 15 17; 23 17; 23 10];
turn = deg2rad(31);
R = [cos(turn) -sin(turn); sin(turn) cos(turn)];
outer = outer*R.'+[300100 4500200];
hole = hole*R.'+[300100 4500200];
field = polyshape({outer(:,1),hole(:,1)},{outer(:,2),hole(:,2)});
assert(field.NumHoles == 1);
h = 1.7;
[field_polygon,field_geometry,field_offset] = Place_Field_Polygon_Grid(field,h);
[again_polygon,again_geometry,again_offset] = Place_Field_Polygon_Grid(field,h);
assert(isequal(field_offset,again_offset) && ...
    isequal(field_geometry.rotation_degrees,again_geometry.rotation_degrees) && ...
    isequaln(field_polygon.Vertices,again_polygon.Vertices), ...
    'The placement is not reproducible.');
field_topo = Build_Polygon_Grid_Topology(field_polygon,h,'majority',field_offset);
captured = sum(field_topo.cell_overlap_area_m2);
outside = sum(field_topo.cell_area_m2-field_topo.cell_overlap_area_m2);
misrepresented = outside+field_topo.field_missed_area_m2;
[~,best] = min(field_geometry.candidate_misrepresented_area_m2);
assert(abs(field_geometry.candidate_misrepresented_area_m2(best)-misrepresented) < 1e-6, ...
    'Boundary integral %.12g m2 against polygon clipping %.12g m2.', ...
    field_geometry.candidate_misrepresented_area_m2(best),misrepresented);
assert(all(field_topo.cell_overlap_fraction > 0.5) && ...
    abs(captured+field_topo.field_missed_area_m2-area(field)) < 1e-6);
assert(any(abs(mod(field_geometry.candidate_rotation_degrees-31,90)) < 1e-9 | ...
    abs(mod(field_geometry.candidate_rotation_degrees-31,90)-90) < 1e-9), ...
    'The rotations tried must include the directions of the field''s edges.');
everywhere = repmat([2.5;-1],1,field_topo.N_cells);
[field_mean,field_total] = Field_Aggregate_Areal(everywhere,field_topo);
[~,~,~,field_area] = Field_Aggregation_Areas(field_topo);
assert(max(abs(field_mean-[2.5;-1])) < 1e-12 && abs(field_area-area(field)) < 1e-6 && ...
    max(abs(field_total-[2.5;-1]*area(field))) < 1e-5, ...
    'Field totals must cover the whole field, not only the kept squares.');
rough = Build_Polygon_Grid_Topology(field_polygon,h,'majority',[0 0]);
rough_error = sum(rough.cell_area_m2-rough.cell_overlap_area_m2)+rough.field_missed_area_m2;
assert(misrepresented <= rough_error+1e-9, ...
    'The placed lattice misrepresents more than the unshifted one.');

%% --- 4 --------------------------------------------------------------------
any_default = Build_Polygon_Grid_Topology(field_polygon,h);
any_explicit = Build_Polygon_Grid_Topology(field_polygon,h,'any',[0 0]);
assert(isequaln(any_default,any_explicit) && ~isfield(any_default,'field_missed_area_m2'));
assert(abs(sum(any_default.cell_overlap_area_m2)-area(field)) < 1e-6 && ...
    any_default.N_cells > field_topo.N_cells);
[~,fixed_geometry] = Orient_Field_Polygon_Grid(field,h,12.5);
assert(abs(fixed_geometry.rotation_degrees-12.5) < 1e-12 && ...
    strcmp(fixed_geometry.orientation_objective,'prescribed_rotation'));

%% --- 5 --------------------------------------------------------------------
expect_error(@() Build_Polygon_Grid_Topology(field_polygon,h,'half'), ...
    'Build_Polygon_Grid_Topology:BadInclusion');
expect_error(@() Build_Polygon_Grid_Topology(field_polygon,h,'majority',[h 0]), ...
    'Build_Polygon_Grid_Topology:BadOffset');
expect_error(@() Build_Polygon_Grid_Topology(field_polygon,200,'majority'), ...
    'Build_Polygon_Grid_Topology:EmptyGrid');
expect_error(@() Place_Field_Polygon_Grid(field,0), ...
    'Place_Field_Polygon_Grid:BadCellSize');
expect_error(@() Build_Projected_Field_Topology(field,[],struct(),1,'majority'), ...
    'Build_Projected_Field_Topology:WidthRequired');
expect_error(@() Build_Projected_Field_Topology(field,h,struct(),1,1.5), ...
    'Build_Projected_Field_Topology:BadCoverageRule');
expect_error(@() Build_Projected_Field_Topology(field,h,struct(),1,'half'), ...
    'Build_Projected_Field_Topology:BadCoverageRule');
bad = field_topo;
bad.field_missed_area_m2 = -1;
expect_error(@() Field_Aggregation_Areas(bad),'Field_Aggregation_Areas:BadMissedArea');

fprintf('test_field_grid_majority: all checks passed\n');

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
