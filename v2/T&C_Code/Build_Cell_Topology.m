%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Build_Cell_Topology  —  axis-aligned rectangular cell partition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% INPUTS
%%% field_side_m    [m]   field side length L
%%% buffer_width_m  [m]   buffer ring width w
%%% main_grid_nx    [-]   main grid columns (optional, default 4)
%%% main_grid_ny    [-]   main grid rows    (optional, default 4)
%%% OUTPUTS
%%% topo   [-]   struct: cell geometry, edge table, adj matrix
%%%
%%%   Topology (default: 24 cells)
%%%   Main: nx x ny cells (default 4x4 = 16)
%%%   Buffer: 4 side strips + 4 corner squares (8 cells)
%%%
%%%      +-----+-----------------------+-----+
%%%      | Cor |         Buf_N         | Cor |
%%%      +-----+-----------------------+-----+
%%%      |     |    Main   |    Main   |     |
%%%      | Buf |  (4x4=16) |          | Buf |
%%%      |  W  +-----------+-----------+  E  |
%%%      |     |    Main   |    Main   |     |
%%%      +-----+-----------+-----------+-----+
%%%      | Cor |         Buf_S         | Cor |
%%%      +-----+-----------------------+-----+
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function topo = Build_Cell_Topology(field_side_m, buffer_width_m, varargin)

% Parse optional arguments: main_grid_nx, main_grid_ny
if nargin < 1 || isempty(field_side_m)
    error('Build_Cell_Topology: field_side_m must be supplied by caller.');
end
if nargin < 2 || isempty(buffer_width_m)
    error('Build_Cell_Topology: buffer_width_m must be supplied by caller (MOD_PARAM.BufferWidth_m).');
end
main_grid_nx = 4; main_grid_ny = 4;  % defaults
if nargin >= 3 && ~isempty(varargin{1}); main_grid_nx = varargin{1}; end
if nargin >= 4 && ~isempty(varargin{2}); main_grid_ny = varargin{2}; end

L = field_side_m;
w = buffer_width_m;

if ~(isnumeric(L) && isscalar(L) && isfinite(L) && L > 0 && ...
        isnumeric(w) && isscalar(w) && isfinite(w) && w > 0 && L > 2*w)
    error('Build_Cell_Topology:BadDomainSize', ...
        'field_side_m and buffer_width_m must be positive, with field_side_m > 2*buffer_width_m.');
end
if ~(isnumeric(main_grid_nx) && isscalar(main_grid_nx) && isfinite(main_grid_nx) && ...
        main_grid_nx == fix(main_grid_nx) && main_grid_nx >= 1 && ...
        isnumeric(main_grid_ny) && isscalar(main_grid_ny) && isfinite(main_grid_ny) && ...
        main_grid_ny == fix(main_grid_ny) && main_grid_ny >= 1)
    error('Build_Cell_Topology:BadGridSize', ...
        'main_grid_nx and main_grid_ny must be positive integers.');
end

%%%%%  cell coordinates
%%% main area extent (inside buffer ring)
main_x = linspace(w, L-w, main_grid_nx+1);
main_y = linspace(w, L-w, main_grid_ny+1);

% Build main grid cell coordinates
cell_name = {};
xy = [];
n_main = main_grid_nx * main_grid_ny;
idx = 0;

% Main cells (row-major order)
for iy = 1:main_grid_ny
    for ix = 1:main_grid_nx
        idx = idx + 1;
        x_lo = main_x(ix);
        x_hi = main_x(ix+1);
        y_lo = main_y(iy);
        y_hi = main_y(iy+1);
        xy = [xy; x_lo, x_hi, y_lo, y_hi];
        cell_name{idx} = sprintf('Cell_Main_%02d', idx);
    end
end

% Store main cell count for region assignment
N_main = n_main;

% Buffer cells: 4 full-length sides + 4 corners = 8 cells
% North buffer (one full-length strip, width w) spanning x in [w, L-w]
xy = [xy; w, L-w, L-w, L];
cell_name{end+1} = 'Cell_Buf_N';
% South buffer (one full-length strip) spanning x in [w, L-w]
xy = [xy; w, L-w, 0, w];
cell_name{end+1} = 'Cell_Buf_S';
% East buffer (one full-length strip) spanning y in [w, L-w]
xy = [xy; L-w, L, w, L-w];
cell_name{end+1} = 'Cell_Buf_E';
% West buffer (one full-length strip) spanning y in [w, L-w]
xy = [xy; 0, w, w, L-w];
cell_name{end+1} = 'Cell_Buf_W';
% Corner buffers (4)
xy = [xy; 0, w, L-w, L; L-w, L, L-w, L; 0, w, 0, w; L-w, L, 0, w];
cell_name{end+1} = 'Cell_Cor_NW';
cell_name{end+1} = 'Cell_Cor_NE';
cell_name{end+1} = 'Cell_Cor_SW';
cell_name{end+1} = 'Cell_Cor_SE';

N_cells = size(xy,1);

cell_x_lo = xy(:,1); cell_x_hi = xy(:,2);
cell_y_lo = xy(:,3); cell_y_hi = xy(:,4);
cell_dx   = cell_x_hi - cell_x_lo;
cell_dy   = cell_y_hi - cell_y_lo;
cell_xc   = 0.5*(cell_x_lo + cell_x_hi);
cell_yc   = 0.5*(cell_y_lo + cell_y_hi);
cell_area_m2 = cell_dx .* cell_dy;

% Region + subtype labels (first N_main are 'main', rest are 'buffer')
cell_region  = repmat({'buffer'}, N_cells, 1);
cell_region(1:N_main) = {'main'};

% Subtype: main cells, length-wise buffers, corner buffers
cell_subtype = repmat({'buf_mid'}, N_cells, 1);
cell_subtype(1:N_main) = {'main'};
% Last 8 are buffer cells (4 length-wise sides + 4 corners)
cell_subtype(N_main+1:N_main+4) = {'buf_mid'};      % N, S, E, W side strips
cell_subtype(N_main+5:N_cells)  = {'buf_corner'};   % NW, NE, SW, SE corners

% --- Closed-form area sanity check (independent of the table) ---
area_check = (L-2*w)^2 + 4*w*(L-2*w) + 4*w^2;
assert(abs(sum(cell_area_m2) - L^2) < 1e-9*L^2, ...
       'Build_Cell_Topology: cell areas do not sum to L^2');
assert(abs(area_check - L^2) < 1e-9*L^2, ...
       'Build_Cell_Topology: closed-form area identity violated');

%%%%%  edge table
edge_a = []; edge_b = [];
edge_face_width_m   = [];
edge_face_distance_m = [];
edge_face_distance_a_m = [];
edge_face_distance_b_m = [];

tol = 1e-9*L;
for a = 1:N_cells
    for b = a+1:N_cells
        % Vertical shared face: cell_a.x_hi == cell_b.x_lo  (or vice versa)
        if abs(cell_x_hi(a) - cell_x_lo(b)) < tol
            yo_lo = max(cell_y_lo(a), cell_y_lo(b));
            yo_hi = min(cell_y_hi(a), cell_y_hi(b));
            if yo_hi - yo_lo > tol
                edge_a(end+1,1) = a;
                edge_b(end+1,1) = b;
                edge_face_width_m(end+1,1)    = yo_hi - yo_lo;
                edge_face_distance_m(end+1,1) = cell_dx(a)/2 + cell_dx(b)/2;
                edge_face_distance_a_m(end+1,1) = cell_dx(a)/2;
                edge_face_distance_b_m(end+1,1) = cell_dx(b)/2;
            end
        elseif abs(cell_x_hi(b) - cell_x_lo(a)) < tol
            yo_lo = max(cell_y_lo(a), cell_y_lo(b));
            yo_hi = min(cell_y_hi(a), cell_y_hi(b));
            if yo_hi - yo_lo > tol
                edge_a(end+1,1) = a;
                edge_b(end+1,1) = b;
                edge_face_width_m(end+1,1)    = yo_hi - yo_lo;
                edge_face_distance_m(end+1,1) = cell_dx(a)/2 + cell_dx(b)/2;
                edge_face_distance_a_m(end+1,1) = cell_dx(a)/2;
                edge_face_distance_b_m(end+1,1) = cell_dx(b)/2;
            end
        end
        % Horizontal shared face: cell_a.y_hi == cell_b.y_lo (or vv)
        if abs(cell_y_hi(a) - cell_y_lo(b)) < tol
            xo_lo = max(cell_x_lo(a), cell_x_lo(b));
            xo_hi = min(cell_x_hi(a), cell_x_hi(b));
            if xo_hi - xo_lo > tol
                edge_a(end+1,1) = a;
                edge_b(end+1,1) = b;
                edge_face_width_m(end+1,1)    = xo_hi - xo_lo;
                edge_face_distance_m(end+1,1) = cell_dy(a)/2 + cell_dy(b)/2;
                edge_face_distance_a_m(end+1,1) = cell_dy(a)/2;
                edge_face_distance_b_m(end+1,1) = cell_dy(b)/2;
            end
        elseif abs(cell_y_hi(b) - cell_y_lo(a)) < tol
            xo_lo = max(cell_x_lo(a), cell_x_lo(b));
            xo_hi = min(cell_x_hi(a), cell_x_hi(b));
            if xo_hi - xo_lo > tol
                edge_a(end+1,1) = a;
                edge_b(end+1,1) = b;
                edge_face_width_m(end+1,1)    = xo_hi - xo_lo;
                edge_face_distance_m(end+1,1) = cell_dy(a)/2 + cell_dy(b)/2;
                edge_face_distance_a_m(end+1,1) = cell_dy(a)/2;
                edge_face_distance_b_m(end+1,1) = cell_dy(b)/2;
            end
        end
    end
end

N_edges      = numel(edge_a);
edge_active  = ones(N_edges,1);

% For 4x4 (16 main + 8 buffer = 24 cells), there are 48 edges.
% For 2x2 (4 main + 8 buffer = 12 cells), there are 20 edges.
fprintf('Build_Cell_Topology: created %d cells with %d edges\n', N_cells, N_edges);

% Sparse adjacency
adj = sparse([edge_a; edge_b], [edge_b; edge_a], 1, N_cells, N_cells) > 0;

%%%%%  outer-face widths
outer_face_width_m = zeros(N_cells,1);
for c = 1:N_cells
    if abs(cell_x_lo(c) - 0) < tol
        outer_face_width_m(c) = outer_face_width_m(c) + cell_dy(c);
    end
    if abs(cell_x_hi(c) - L) < tol
        outer_face_width_m(c) = outer_face_width_m(c) + cell_dy(c);
    end
    if abs(cell_y_lo(c) - 0) < tol
        outer_face_width_m(c) = outer_face_width_m(c) + cell_dx(c);
    end
    if abs(cell_y_hi(c) - L) < tol
        outer_face_width_m(c) = outer_face_width_m(c) + cell_dx(c);
    end
end
% Sanity: total outer perimeter == 4 L
assert(abs(sum(outer_face_width_m) - 4*L) < tol*N_cells, ...
       'Build_Cell_Topology: outer-face widths do not sum to 4L');
% Sanity: interior main cells have no outer face
assert(all(outer_face_width_m(1:N_main) == 0), ...
       'Build_Cell_Topology: a main cell unexpectedly touches the boundary');

%%%%%  pack output struct
topo.field_side_m       = L;
topo.buffer_width_m     = w;
topo.N_cells            = N_cells;
topo.cell_id            = (1:N_cells)';
topo.cell_name          = cell_name;
topo.cell_region        = cell_region;
topo.cell_subtype       = cell_subtype;
topo.cell_x_lo          = cell_x_lo;
topo.cell_x_hi          = cell_x_hi;
topo.cell_y_lo          = cell_y_lo;
topo.cell_y_hi          = cell_y_hi;
topo.cell_xc            = cell_xc;
topo.cell_yc            = cell_yc;
topo.cell_dx            = cell_dx;
topo.cell_dy            = cell_dy;
topo.cell_area_m2       = cell_area_m2;
topo.N_edges            = N_edges;
topo.edge_a             = edge_a;
topo.edge_b             = edge_b;
topo.edge_face_width_m  = edge_face_width_m;
topo.edge_face_distance_m = edge_face_distance_m;
topo.edge_face_distance_a_m = edge_face_distance_a_m;
topo.edge_face_distance_b_m = edge_face_distance_b_m;
topo.edge_active        = edge_active;
topo.adj                = adj;
topo.outer_face_width_m = outer_face_width_m;

end
