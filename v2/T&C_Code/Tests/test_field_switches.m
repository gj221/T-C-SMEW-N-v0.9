function test_field_switches
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The field options from v1.0 and their switches (OPT_DEM,
%% OPT_CELL_INPUTS): that topography adds only its own fields, that the
%% overland router is inert without topography and strictly downhill with
%% it, and that Run_One_Subfield refuses a field option it would otherwise
%% ignore, a field grid on a one-column run and a margin layout without a
%% margin, before it reads anything. test_2d_topographic_smoke checks the
%% frame's own refusals; test_total_head_darcy the Darcy half of "off is
%% flat".
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_dir = fileparts(mfilename('fullpath'));
code_dir = fileparts(test_dir);
addpath(code_dir,fullfile(code_dir,'Inputs'),fullfile(code_dir,'Parameters'), ...
    fullfile(fileparts(code_dir),'run'),test_dir);

%% 1. Topography adds its own fields and leaves the rest alone: the
%%    per-cell inputs read only cell ids and coordinates.
topo = Build_Cell_Topology(60,5,2,2);
dem = synthetic_hillslope_dem('plane','Extent',[70 70],'CellSize',5);
t1 = Attach_Topography_To_Topology(topo,dem);
f0 = fieldnames(topo);
for k = 1:numel(f0)
    assert(isequal(topo.(f0{k}),t1.(f0{k})), ...
        'Attaching topography changed %s.',f0{k});
end
assert(isfield(t1,'cell_z_m') && isfield(t1,'edge_overland_weight_ab') && ...
    t1.topography_attached && ~isfield(t1,'cell_xc_projected_m'), ...
    'Topography adds its own fields, marks itself attached and does not georeference.');

%% 2. Without topography the router keeps every cell's runoff where it
%%    formed and exports nothing.
N = topo.N_cells;
runoff = 2*ones(N,1);
C = [ones(N,1) 0.5*ones(N,1)];
[q0,C0,x0,m0] = Surface_Runoff_Routing(topo,runoff,C,1);
assert(isequal(q0,runoff) && isequal(C0,C) && x0 == 0 && all(m0 == 0), ...
    'Without topography every cell keeps its own runoff.');

%% 3. On a sloped plane runoff leaves at the outlet(s), only there, and
%%    never runs uphill.
[q1,~,x1,~,~,x1_cell] = Surface_Runoff_Routing(t1,runoff,C,1);
assert(x1 > 0 && all(x1_cell(~t1.cell_is_outlet) == 0) && sum(q1) < sum(runoff), ...
    'A sloped plane exports runoff at its outlet cells only.');
z = t1.cell_z_m(:);
for e = 1:t1.N_edges
    a = t1.edge_a(e); b = t1.edge_b(e);
    if z(a) > z(b)
        assert(t1.edge_overland_weight_ba(e) == 0,'Edge %d routes uphill.',e);
    elseif z(b) > z(a)
        assert(t1.edge_overland_weight_ab(e) == 0,'Edge %d routes uphill.',e);
    end
end

%% 4. Run_One_Subfield's field options, refused before any input is read.
%%    Any existing file passes the existence check; the refusals come
%%    first.
some_file = [mfilename('fullpath') '.m'];
assert_error(@() Run_One_Subfield(2,'GridCellSizeM',5), ...
    'Run_One_Subfield:IncompleteFieldGrid');
assert_error(@() Run_One_Subfield(2,'FieldDEM',test_dir), ...
    'Run_One_Subfield:IncompleteFieldGrid');
assert_error(@() Run_One_Subfield(2,'MarginWidthM',5), ...
    'Run_One_Subfield:IncompleteFieldGrid');
assert_error(@() Run_One_Subfield(2,'FieldShapefile',some_file, ...
    'FieldTopologyFile',some_file),'Run_One_Subfield:ConflictingFieldGrid');
assert_error(@() Run_One_Subfield(2,'FieldTopologyFile',some_file,'GridNx',3), ...
    'Run_One_Subfield:ConflictingGridOptions');
assert_error(@() Run_One_Subfield(2,'FieldShapefile',some_file, ...
    'FieldMinCoverage',1.5),'Run_One_Subfield:BadOptions');
assert_error(@() Run_One_Subfield(2,'FieldShapefile',fullfile(test_dir,'none.shp')), ...
    'Run_One_Subfield:MissingInputFile');
assert_error(@() Run_One_Subfield(2,'CellAgronomicFile',some_file, ...
    'Pcla',0.1,'Psan',0.2),'Run_One_Subfield:ConflictingTexture');
assert_error(@() Run_One_Subfield(11,'FieldTopologyFile',some_file), ...
    'Run_One_Subfield:BadFieldLayout');
assert_error(@() Run_One_Subfield(5,'FieldShapefile',some_file,'MarginWidthM',0), ...
    'Run_One_Subfield:MarginRequired');

%% 5. The frame's switches come from what the run was given, and a field
%%    grid names its result by its build record.
source = fileread(fullfile(fileparts(code_dir),'run','Run_One_Subfield.m'));
assert(contains(source, ...
    'OPT_DEM = double(isfield(topo,''topography_attached'') && topo.topography_attached);') && ...
    contains(source,'OPT_CELL_INPUTS = double(~isempty(opt.CellAgronomicFile));') && ...
    contains(source,'field_run_tag(topo)'), ...
    'Run_One_Subfield sets the frame''s switches from the grid and the cell table.');

fprintf(['PASS: topography adds only its own fields; the router is inert ' ...
    'without it and strictly downhill with it; Run_One_Subfield refuses ' ...
    'field options it would ignore and grids it cannot use.\n']);
end

function assert_error(fun,identifier)
threw = false;
try
    fun();
catch ME
    threw = strcmp(ME.identifier,identifier);
    if ~threw
        error('Expected %s, got %s: %s',identifier,ME.identifier,ME.message);
    end
end
assert(threw,'Expected error %s, but none was raised.',identifier);
end
