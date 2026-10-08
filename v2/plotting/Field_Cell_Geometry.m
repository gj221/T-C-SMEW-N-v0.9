function G = Field_Cell_Geometry(topo,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Field_Cell_Geometry
%%
%% What the field figures draw of a cell topology (Build_Cell_Topology,
%% Build_Polygon_Grid_Topology, a saved field grid): one polygon per cell
%% (the square, or on a field grid the part of it inside the field), the
%% whole squares, the field outline, the faces between main and margin
%% cells, the cell centres and the north direction.
%%
%%   G = Field_Cell_Geometry(topo)
%%   G = Field_Cell_Geometry(topo,'Coordinates','local','ClipToField',false)
%%
%% Options
%%   'Coordinates'  'auto' (projected when the grid carries its rotation
%%                  and projected origin, else local), 'local' (grid
%%                  metres, axes along the grid) or 'projected' (CRS
%%                  metres less a round offset, north up)
%%   'ClipToField'  draw only the part of each perimeter square inside the
%%                  field polygon (true); the model weighs a square by
%%                  that overlap (Field_Aggregation_Areas)
%%
%% G.faces/G.vertices feed patch; G.face_cell gives each face's cell (a
%% clipped square can leave more than one piece). G.dem_to_plot maps the
%% coordinates of topo.topography_grid to the plot frame.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
p.addParameter('Coordinates','auto',@(x) any(strcmpi(x,{'auto','local','projected'})));
p.addParameter('ClipToField',true,@(x) islogical(x) || isnumeric(x));
p.parse(varargin{:});
o = p.Results;

here = fileparts(mfilename('fullpath'));
if exist('Grid_Local_To_Projected','file') ~= 2
    addpath(fullfile(fileparts(here),'T&C_Code'));
end

N = topo.N_cells;
can_project = isfield(topo,'grid_rotation_degrees') && ...
    isfield(topo,'projected_grid_origin_m') && numel(topo.projected_grid_origin_m) == 2;
coords = lower(o.Coordinates);
if strcmp(coords,'auto')
    if can_project; coords = 'projected'; else; coords = 'local'; end
elseif strcmp(coords,'projected') && ~can_project
    error('Field_Cell_Geometry:NoProjection', ...
        'This grid carries no projected origin and rotation; use local coordinates.');
end

xlo = topo.cell_x_lo(:); xhi = topo.cell_x_hi(:);
ylo = topo.cell_y_lo(:); yhi = topo.cell_y_hi(:);
region = cellstr(topo.cell_region(:));
is_main = strcmp(region,'main');

%% Field outline in the local frame: the polygon of a field grid, the
%% bounding square of the plot grid.
if isfield(topo,'field_boundary_x_m') && ~isempty(topo.field_boundary_x_m)
    fx = topo.field_boundary_x_m(:); fy = topo.field_boundary_y_m(:);
    has_polygon = true;
else
    fx = [min(xlo); max(xhi); max(xhi); min(xlo); min(xlo)];
    fy = [min(ylo); min(ylo); max(yhi); max(yhi); min(ylo)];
    has_polygon = false;
end

if isfield(topo,'cell_overlap_area_m2')
    fraction = topo.cell_overlap_area_m2(:)./topo.cell_area_m2(:);
else
    fraction = ones(N,1);
end

%% Cell faces: whole squares, or the in-field part of partial squares.
square_x = [xlo xhi xhi xlo];
square_y = [ylo ylo yhi yhi];
piece_x = num2cell(square_x,2);
piece_y = num2cell(square_y,2);
piece_cell = num2cell((1:N).');
clip = logical(o.ClipToField) && has_polygon;
if clip
    state = warning('off','all');
    restore = onCleanup(@() warning(state));
    field = polyshape(fx,fy);
    partial = find(fraction < 1-1e-9);
    for k = 1:numel(partial)
        c = partial(k);
        part = intersect(polyshape(square_x(c,:),square_y(c,:)),field);
        part = rmholes(part);
        pieces = regions(part);
        px = cell(numel(pieces),1); py = px;
        for r = 1:numel(pieces)
            [bx,by] = boundary(pieces(r));
            ok = isfinite(bx) & isfinite(by);
            px{r} = bx(ok).'; py{r} = by(ok).';
        end
        piece_x{c} = px; piece_y{c} = py;
        piece_cell{c} = repmat(c,numel(pieces),1);
    end
    clear restore
end
[face_x,face_y,face_cell] = flatten(piece_x,piece_y,piece_cell);

%% Faces between main and margin cells, as segments in the local frame.
zx = []; zy = [];
tol = 1e-9*max(1,max(xhi)-min(xlo));
for e = 1:topo.N_edges
    a = topo.edge_a(e); b = topo.edge_b(e);
    if is_main(a) == is_main(b); continue; end
    if abs(xhi(a)-xlo(b)) < tol || abs(xhi(b)-xlo(a)) < tol
        x = xhi(a)*(abs(xhi(a)-xlo(b)) < tol) + xhi(b)*(abs(xhi(a)-xlo(b)) >= tol);
        y = [max(ylo(a),ylo(b)) min(yhi(a),yhi(b))];
        zx = [zx; x; x; NaN]; zy = [zy; y(:); NaN]; %#ok<AGROW>
    else
        y = yhi(a)*(abs(yhi(a)-ylo(b)) < tol) + yhi(b)*(abs(yhi(a)-ylo(b)) >= tol);
        x = [max(xlo(a),xlo(b)) min(xhi(a),xhi(b))];
        zx = [zx; x(:); NaN]; zy = [zy; y; y; NaN]; %#ok<AGROW>
    end
end

cx = topo.cell_xc(:); cy = topo.cell_yc(:);
offset = [0 0];
if strcmp(coords,'projected')
    to_frame = @(x,y) local_to_projected(x,y,topo);
    [ox,oy] = to_frame(fx,fy);
    offset = [floor(min(ox)/10)*10 floor(min(oy)/10)*10];
    xlabel_text = sprintf('Easting - %d m',offset(1));
    ylabel_text = sprintf('Northing - %d m',offset(2));
    %% The zone, not the EPSG code: v1.0 labels a NAD83 shapefile's grid
    %% WGS84 UTM (326xx) without shifting its coordinates, and the DEM is
    %% carried over the same way, so the numbers are the shapefile's.
    if isfield(topo,'projected_crs_epsg') && ~isempty(topo.projected_crs_epsg)
        epsg = double(topo.projected_crs_epsg);
        if epsg > 32600 && epsg <= 32660
            xlabel_text = sprintf('%s (UTM %dN)',xlabel_text,epsg-32600);
        elseif epsg > 32700 && epsg <= 32760
            xlabel_text = sprintf('%s (UTM %dS)',xlabel_text,epsg-32700);
        elseif epsg > 26900 && epsg <= 26923
            xlabel_text = sprintf('%s (UTM %dN)',xlabel_text,epsg-26900);
        end
    end
    north = [0 1];
else
    to_frame = @(x,y) deal(x,y);
    xlabel_text = 'Distance along the grid x [m]';
    ylabel_text = 'Distance along the grid y [m]';
    if can_project
        t = deg2rad(topo.grid_rotation_degrees);
        north = [sin(t) cos(t)];
    else
        north = [];
    end
end
shift = @(x,y) deal(x-offset(1),y-offset(2));
plot_xy = @(x,y) apply(to_frame,shift,x,y);

[G.outline_x,G.outline_y] = plot_xy(fx,fy);
[G.zone_x,G.zone_y] = plot_xy(zx,zy);
[G.centre_x,G.centre_y] = plot_xy(cx,cy);
[vx,vy] = plot_xy(face_x,face_y);
[G.faces,G.vertices] = to_patch(vx,vy);
G.face_cell = face_cell;
[sx,sy] = plot_xy(square_x,square_y);
G.square_faces = reshape(1:4*N,N,4);
G.square_vertices = [sx(:) sy(:)];
G.coordinates = coords;
G.offset = offset;
G.xlabel = xlabel_text;
G.ylabel = ylabel_text;
G.north = north;
G.clipped = clip;
G.N_cells = N;
G.is_main = is_main;
G.cell_size_m = median(xhi-xlo);
G.extent = [min(G.square_vertices(:,1)) max(G.square_vertices(:,1)) ...
    min(G.square_vertices(:,2)) max(G.square_vertices(:,2))];
G.field_extent = [min(G.outline_x) max(G.outline_x) min(G.outline_y) max(G.outline_y)];

%% The DEM grid of Attach_Topography_To_Topology is in the frame its cell
%% centres were sampled in: projected when the topology has projected
%% centres, local otherwise.
dem_projected = isfield(topo,'cell_xc_projected_m') && isfield(topo,'cell_yc_projected_m');
if dem_projected && strcmp(coords,'local')
    G.dem_to_plot = @(x,y) projected_to_local(x,y,topo);
elseif ~dem_projected && strcmp(coords,'projected')
    G.dem_to_plot = plot_xy;
elseif dem_projected
    G.dem_to_plot = shift;
else
    G.dem_to_plot = @(x,y) deal(x,y);
end
end

%% =======================================================================
function [x,y] = apply(to_frame,shift,x,y)
[x,y] = to_frame(x,y);
[x,y] = shift(x,y);
end

function [x,y] = local_to_projected(x,y,topo)
sz = size(x);
[x,y] = Grid_Local_To_Projected(x(:),y(:),topo);
x = reshape(x,sz); y = reshape(y,sz);
end

function [x,y] = projected_to_local(xp,yp,topo)
t = deg2rad(topo.grid_rotation_degrees);
o = topo.projected_grid_origin_m;
dx = xp-o(1); dy = yp-o(2);
x = cos(t)*dx + sin(t)*dy;
y = -sin(t)*dx + cos(t)*dy;
end

%% =======================================================================
function [fx,fy,fc] = flatten(px,py,pc)
%% Faces as rows padded with NaN: px{c} is a row of vertices, or a cell
%% of rows for a clipped square's pieces.
rows_x = {}; rows_y = {}; fc = [];
for c = 1:numel(px)
    if iscell(px{c})
        rows_x = [rows_x; px{c}(:)]; rows_y = [rows_y; py{c}(:)]; %#ok<AGROW>
    else
        rows_x{end+1,1} = px{c}; rows_y{end+1,1} = py{c}; %#ok<AGROW>
    end
    fc = [fc; pc{c}(:)]; %#ok<AGROW>
end
n = max(cellfun(@numel,rows_x));
fx = nan(numel(rows_x),n); fy = fx;
for k = 1:numel(rows_x)
    fx(k,1:numel(rows_x{k})) = rows_x{k};
    fy(k,1:numel(rows_y{k})) = rows_y{k};
end
end

function [faces,vertices] = to_patch(fx,fy)
%% patch Faces/Vertices from NaN-padded face rows.
[n,m] = size(fx);
vertices = [reshape(fx.',[],1) reshape(fy.',[],1)];
faces = reshape(1:n*m,m,n).';
faces(isnan(fx)) = NaN;
vertices(isnan(vertices(:,1)),:) = 0;
end
