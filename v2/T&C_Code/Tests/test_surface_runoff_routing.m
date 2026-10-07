function test_surface_runoff_routing()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% test_surface_runoff_routing
%%
%%   1  no topography attached -> every cell keeps its own runoff
%%   2  MFD overland weights out of a cell sum to 1 where it has a
%%      downslope neighbour
%%   3  water conservation: delivered + exited == input
%%   4  on a tilted chain generating runoff everywhere, every cell is at
%%      capacity: the cascade delivers everything to the outlet in-step
%%   5  runoff from the ridge onto DRY (absorbing) cells stops at the
%%      first dry cell and is offered for re-infiltration; a wet band
%%      transmits through to the first dry cell below it; a dry OUTLET
%%      absorbs rather than exports
%%   6  solute rides the water, upwind from the contributing cell
%%   7  a convergent hollow: wet flanks drain to a dry valley line
%%   8  an interior pit keeps its water (offered back, nothing exits)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n=== test_surface_runoff_routing ===\n');
here = fileparts(mfilename('fullpath'));
code = fileparts(here);
addpath(code,here);
n_pass = 0; n_fail = 0;
    function check(name,cond,detail)
        if nargin < 3; detail = ''; end
        if cond; n_pass = n_pass+1; fprintf('  PASS  %s\n',name);
        else;     n_fail = n_fail+1; fprintf('  FAIL  %s  %s\n',name,detail); end
    end

%% ---------------------------------------------------------------------
%% 1. No topography -> pass-through
%% ---------------------------------------------------------------------
row = row_topology(5,10);
runoff = [3;1;4;1;5];
[q,~,exit_w] = Surface_Runoff_Routing(row,runoff,[],1);
check('no topography: q_deliver == runoff, nothing exits', ...
    isequal(q,runoff) && exit_w == 0);

%% ---------------------------------------------------------------------
%% 2-4. Tilted chain: 6% slope toward cell n
%% ---------------------------------------------------------------------
n = 8; h = 5;
chain = row_topology(n,h);
z = ((n:-1:1)'-1)*h*0.06;
chain = Attach_Topography_To_Topology(chain,struct('cell_z_m',z));

% weights out of an interior cell sum to 1
c = 4;
w_out = 0;
for e = 1:chain.N_edges
    if chain.edge_a(e) == c; w_out = w_out + chain.edge_overland_weight_ab(e); end
    if chain.edge_b(e) == c; w_out = w_out + chain.edge_overland_weight_ba(e); end
end
check('MFD weights out of an interior cell sum to 1', abs(w_out-1) < 1e-9, ...
    sprintf('sum = %.6f',w_out));

runoff = 2*ones(n,1);
[q,~,exit_w,~,diag] = Surface_Runoff_Routing(chain,runoff,[],1);
A = chain.cell_area_m2(:);
in_w  = sum(runoff/1000.*A);
out_w = sum(q/1000.*A) + exit_w;
check('water conservation (delivered + exited == input)', ...
    abs(in_w-out_w) < 1e-9*in_w, sprintf('in %.6g out %.6g',in_w,out_w));
% every cell generated runoff -> every cell is at capacity -> the whole
% slope is connected and everything reaches the outlet this step
check('a hillslope generating runoff everywhere exports all of it in-step', ...
    abs(exit_w - in_w) < 1e-9*in_w && all(q == 0), ...
    sprintf('exit %.6g of %.6g; max q %.3g',exit_w,in_w,max(q)));
check('the outlet is the downslope end cell', chain.cell_is_outlet(n) && ...
    ~any(chain.cell_is_outlet(1:n-1)));
check('every edge carried water downslope (a->b positive), accumulating', ...
    all(diag.edge_flux_m3 > 0) && all(diff(diag.edge_flux_m3) > 0));
check('every runoff-generating cell is flagged transmitting', ...
    all(diag.transmitting));

%% ---------------------------------------------------------------------
%% 5. Runoff onto dry ground stops at the first absorbing cell
%% ---------------------------------------------------------------------
ridge_only = zeros(n,1); ridge_only(1) = 2;
[q5,~,exit_w5,~,d5] = Surface_Runoff_Routing(chain,ridge_only,[],1);
check('ridge runoff onto a dry slope is delivered to the next cell only', ...
    abs(q5(2) - 2) < 1e-9 && all(q5([1 3:n]) == 0) && exit_w5 == 0, ...
    sprintf('q5 = %s',mat2str(q5.',3)));
check('only the ridge cell is transmitting', ...
    d5.transmitting(1) && ~any(d5.transmitting(2:n)));
wet_band = zeros(n,1); wet_band(1:3) = [1;2;3];
[q5b,~,exit_w5b] = Surface_Runoff_Routing(chain,wet_band,[],1);
check('a wet band transmits through itself to the first dry cell below', ...
    abs(q5b(4) - 6) < 1e-9 && all(q5b([1:3 5:n]) == 0) && exit_w5b == 0, ...
    sprintf('q5b = %s',mat2str(q5b.',3)));
dry_outlet = 2*ones(n,1); dry_outlet(n) = 0;      % whole slope wet, toe dry
[q5c,~,exit_w5c] = Surface_Runoff_Routing(chain,dry_outlet,[],1);
check('a DRY outlet cell absorbs the cascade instead of exporting it', ...
    exit_w5c == 0 && abs(q5c(n) - 2*(n-1)) < 1e-9 && all(q5c(1:n-1) == 0), ...
    sprintf('exit %.3g, q(n) = %.3g',exit_w5c,q5c(n)));

%% ---------------------------------------------------------------------
%% 6. Solute rides the water, upwind from the contributing cell
%% ---------------------------------------------------------------------
C = zeros(n,2);
C(1,1) = 40;                        % only the ridge cell has species 1
[~,C6,~,exit_m6] = Surface_Runoff_Routing(chain,ridge_only,C,1);
check('ridge cell 1 delivers its species-1 load to cell 2 at conc 40', ...
    abs(C6(2,1) - 40) < 1e-9 && all(C6([1 3:n],1) == 0), ...
    sprintf('C6(2,1) = %.6g',C6(2,1)));
check('species-1 does not reach the outlet when the slope below is dry', ...
    exit_m6(1) == 0);
check('species 2 (absent everywhere) moves nothing', abs(exit_m6(2)) < 1e-12);
% confluence mixing: cells 1..3 wet with different concentrations, cell 4 dry
Cb = zeros(n,1); Cb(1:3) = [10;20;30];
[~,C6b,~,~] = Surface_Runoff_Routing(chain,wet_band,Cb,1);
expect = (1*10+2*20+3*30)/6;
check('a confluence delivers the flow-weighted mixed concentration', ...
    abs(C6b(4) - expect) < 1e-9, sprintf('%.4f vs %.4f',C6b(4),expect));

%% ---------------------------------------------------------------------
%% 7. Convergent hollow: wet flanks drain to a dry valley line
%% ---------------------------------------------------------------------
dem = synthetic_hillslope_dem('hollow','Extent',[40 40],'CellSize',2, ...
    'Slope',0.05,'Relief',2);
grid_topo = grid_topology(dem);
grid_topo = Attach_Topography_To_Topology(grid_topo,dem);
xmid = 0.5*(min(grid_topo.cell_xc)+max(grid_topo.cell_xc));
near_valley = abs(grid_topo.cell_xc - xmid) < 3;
far_ridge   = abs(grid_topo.cell_xc - xmid) > 15;
runoff_g = ones(grid_topo.N_cells,1);
runoff_g(near_valley) = 0;                       % the valley floor is dry
[qg,~,exit_wg] = Surface_Runoff_Routing(grid_topo,runoff_g,[],1);
Ag = grid_topo.cell_area_m2(:);
check('hollow: the dry valley line collects the flank runoff', ...
    mean(qg(near_valley)) > 0 && all(qg(far_ridge) == 0), ...
    sprintf('valley %.3f vs ridge %.3f mm',mean(qg(near_valley)),mean(qg(far_ridge))));
check('hollow: water conservation', ...
    abs(sum(runoff_g/1000.*Ag) - (sum(qg/1000.*Ag)+exit_wg)) < 1e-9*sum(runoff_g/1000.*Ag));
[~,~,exit_all] = Surface_Runoff_Routing(grid_topo,ones(grid_topo.N_cells,1),[],1);
check('hollow: a fully wet hollow exports everything through its outlet(s)', ...
    abs(exit_all - sum(Ag)/1000) < 1e-9*sum(Ag)/1000, ...
    sprintf('exit %.4g of %.4g m3',exit_all,sum(Ag)/1000));

%% ---------------------------------------------------------------------
%% 8. Interior pit keeps its water
%% ---------------------------------------------------------------------
pit = row_topology(3,10);
pit = Attach_Topography_To_Topology(pit,struct('cell_z_m',[2;0;2]), ...
    'OutletCells',[]);
% cell 2 is a pit: lower than both neighbours, and (as a strip) it is on
% the perimeter, so derive_outlets will actually call it an outlet. Force
% the interior-pit case with an explicit non-outlet override:
pit.cell_is_outlet = false(3,1);
[qp,~,exit_wp,~,dp] = Surface_Runoff_Routing(pit,[1;1;1],[],1);
check('pit: nothing exits when no cell is an outlet', exit_wp == 0);
check('pit: all water ends up delivered (offered back)', ...
    abs(sum(qp/1000.*pit.cell_area_m2(:)) - sum([1;1;1]/1000.*pit.cell_area_m2(:))) < 1e-9);
check('pit: the low centre cell holds all the water', ...
    abs(qp(2) - 3) < 1e-9 && qp(1) == 0 && qp(3) == 0 && isequal(dp.pit_cells,2));

%% ---------------------------------------------------------------------
%% 9. Sub-field graph with laterally offset buffer strips: the drop must
%%    be a NODAL difference (so a uniform total head is an equilibrium and
%%    the drops sum to zero around every cycle) and the weights must agree
%%    with the sweep order -- a fully wet field routes without stranding
%%    water.  (A face-normal projection of the drop violated both, and the
%%    run aborted on the mass-balance guard as soon as runoff occurred.)
%% ---------------------------------------------------------------------
sub = Build_Cell_Topology(60,5,3,3);
zs = 100 + 0.06*(max(sub.cell_yc)-sub.cell_yc) + 0.5*sin(sub.cell_xc/15);
sub = Attach_Topography_To_Topology(sub,struct('cell_z_m',zs));
check('buffer-strip graph: edge drops are nodal (z_a - z_b exactly)', ...
    isequal(sub.edge_dz_m, zs(sub.edge_a)-zs(sub.edge_b)));
%% every 4-cycle main(i)->buffer->main(j)->main(i) sums to zero by
%% construction of a nodal field; check the whole edge set as a digraph
cyc = 0;
for e = 1:sub.N_edges
    for f = 1:sub.N_edges
        if e == f; continue; end
        %% two edges that close a cycle through a shared pair of cells
        if sub.edge_a(e) == sub.edge_a(f) && sub.edge_b(e) == sub.edge_b(f)
            cyc = max(cyc, abs(sub.edge_dz_m(e)-sub.edge_dz_m(f)));
        end
    end
end
check('buffer-strip graph: parallel faces carry identical drops', cyc < 1e-12);
offset = abs(sub.edge_centre_distance_m - sub.edge_face_distance_m) > 1e-9;
check('buffer-strip graph: offset faces have a longer centre distance, aligned faces do not', ...
    any(offset) && all(sub.edge_centre_distance_m >= sub.edge_face_distance_m - 1e-12), ...
    sprintf('%d of %d faces offset',nnz(offset),sub.N_edges));
flat = Attach_Topography_To_Topology(sub,struct('cell_z_m',zeros(sub.N_cells,1)));
check('buffer-strip graph: a zero-relief DEM leaves every resistance path unchanged', ...
    isequal(flat.edge_face_distance_surface_a_m, sub.edge_face_distance_a_m) && ...
    isequal(flat.edge_face_distance_surface_b_m, sub.edge_face_distance_b_m));
wet_all = 5*ones(sub.N_cells,1);
[qs,~,exit_s,~,ds] = Surface_Runoff_Routing(sub,wet_all,[],1);
As = sub.cell_area_m2(:);
check('buffer-strip graph, fully wet: routes without stranding water', ...
    abs((sum(qs/1000.*As)+exit_s) - sum(wet_all/1000.*As)) < 1e-9*sum(wet_all/1000.*As) && ...
    exit_s > 0 && all(ds.transmitting));
w_out = accumarray([sub.edge_a(:); sub.edge_b(:)], ...
    [sub.edge_overland_weight_ab(:); sub.edge_overland_weight_ba(:)],[sub.N_cells 1]);
check('buffer-strip graph: every non-outlet cell drains (weights sum to 1)', ...
    all(abs(w_out(~sub.cell_is_outlet) - 1) < 1e-9) && ...
    all(w_out(sub.cell_is_outlet) < 1e-9));

%% ---------------------------------------------------------------------
fprintf('  ---------------------------------------\n');
fprintf('  %d passed, %d failed\n\n',n_pass,n_fail);
if n_fail > 0
    error('test_surface_runoff_routing: %d check(s) failed',n_fail);
end
end

%% =======================================================================
function topo = row_topology(n,h)
topo.N_cells = n;
topo.N_edges = n-1;
topo.edge_a = (1:n-1)';
topo.edge_b = (2:n)';
topo.edge_face_width_m      = repmat(h,n-1,1);
topo.edge_face_distance_m   = repmat(h,n-1,1);
topo.edge_face_distance_a_m = repmat(h/2,n-1,1);
topo.edge_face_distance_b_m = repmat(h/2,n-1,1);
topo.edge_active = ones(n-1,1);
topo.cell_area_m2 = repmat(h^2,n,1);
topo.cell_xc = ((1:n)'-0.5)*h;
topo.cell_yc = repmat(h/2,n,1);
topo.outer_face_width_m = repmat(2*h,n,1);
topo.outer_face_width_m([1 n]) = 3*h;
topo.adj = sparse([topo.edge_a;topo.edge_b],[topo.edge_b;topo.edge_a],1,n,n) > 0;
topo.cell_region = repmat({'main'},n,1);
end

%% =======================================================================
function topo = grid_topology(dem)
%% Regular ny x nx cell graph over a synthetic DEM (test helper).
h = dem.options.CellSize;
nx = numel(dem.x); ny = numel(dem.y);
[C,Rr] = meshgrid(1:nx,1:ny);
id = reshape(1:nx*ny,ny,nx);
topo.N_cells = nx*ny;
topo.cell_xc = dem.x(C(:)).';
topo.cell_yc = dem.y(Rr(:)).';
topo.cell_area_m2 = repmat(h^2,nx*ny,1);
ea = []; eb = [];
for j = 1:nx
    for i = 1:ny
        if j < nx; ea(end+1,1) = id(i,j); eb(end+1,1) = id(i,j+1); end %#ok<AGROW>
        if i < ny; ea(end+1,1) = id(i,j); eb(end+1,1) = id(i+1,j); end %#ok<AGROW>
    end
end
topo.N_edges = numel(ea);
topo.edge_a = ea; topo.edge_b = eb;
topo.edge_face_width_m      = repmat(h,topo.N_edges,1);
topo.edge_face_distance_m   = repmat(h,topo.N_edges,1);
topo.edge_face_distance_a_m = repmat(h/2,topo.N_edges,1);
topo.edge_face_distance_b_m = repmat(h/2,topo.N_edges,1);
topo.edge_active = ones(topo.N_edges,1);
topo.outer_face_width_m = zeros(topo.N_cells,1);
for i = 1:ny
    for j = 1:nx
        e = (i==1)+(i==ny)+(j==1)+(j==nx);
        topo.outer_face_width_m(id(i,j)) = e*h;
    end
end
topo.adj = sparse([ea;eb],[eb;ea],1,topo.N_cells,topo.N_cells) > 0;
topo.cell_region = repmat({'main'},topo.N_cells,1);
end
