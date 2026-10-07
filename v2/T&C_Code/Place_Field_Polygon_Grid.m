function [oriented_polygon,geometry,offset_m] = Place_Field_Polygon_Grid( ...
        field_polygon,cell_size_m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Place a square grid on a field polygon for the majority cell rule.
%%
%% Under the majority rule a square is simulated when more than half of it
%% lies in the field.  For a given lattice that is the exact minimiser of
%% the misrepresented area E -- field land left out plus outside land
%% taken in, the symmetric difference between the field and the simulated
%% squares -- because a kept square adds its outside part (h^2-a) and a
%% dropped one loses its field part (a):  E = sum over squares of
%% min(a,h^2-a).  E still depends on where the lattice lies, and the
%% lattice is chosen by the same criterion.  Along a boundary edge parallel
%% to the grid the error can vanish, so the rotation is taken from the
%% directions of the field's own edges at least one cell long, modulo 90
%% degrees and to 1e-6 degree (which merges edges parallel to round-off and
%% keeps the rotation the same on every platform).  For each rotation the
%% offset [dx dy] in [0,h)^2, the shift of the lattice's lower-left corner
%% below and left of the polygon's bounding box, comes from a 16 x 16 scan
%% of the offset torus and a compass search from its best point down to
%% 1e-4 h.  Values within 1e-9 of the field area of the best are ties,
%% settled by candidate order, so round-off cannot move the grid between
%% platforms; offsets are rounded to 0.1 mm.
%%
%% The overlap areas come from the boundary integral of the polygon, cell
%% by cell:  a(i,j) = closed integral of clamp(X-x_i,0,h) 1_j(Y) dY, each
%% edge cut at the grid lines.  It is exact for straight edges, so every
%% evaluation is exact without clipping the polygon square by square.
%%
%% Returns the polygon in the chosen grid frame and its geometry
%% (Orient_Field_Polygon_Grid at the chosen rotation, plus the candidate
%% rotations and their misrepresented areas), and the offset for
%% Build_Polygon_Grid_Topology.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~(isa(field_polygon,'polyshape') && isscalar(field_polygon)) || ...
        isempty(field_polygon.Vertices) || area(field_polygon) <= 0
    error('Place_Field_Polygon_Grid:EmptyField', ...
        'field_polygon must have finite positive area.');
end
if ~(isnumeric(cell_size_m) && isscalar(cell_size_m) && ...
        isfinite(cell_size_m) && cell_size_m > 0)
    error('Place_Field_Polygon_Grid:BadCellSize', ...
        'cell_size_m must be a finite positive scalar in metres.');
end

[xa,ya,xb,yb] = polygon_edges(field_polygon.Vertices);
edge_length = hypot(xb-xa,yb-ya);
candidate = edge_length >= cell_size_m;
if ~any(candidate); candidate = edge_length > 0; end
candidate_degrees = unique(mod(round(rad2deg( ...
    atan2(yb(candidate)-ya(candidate),xb(candidate)-xa(candidate))),6),90));
tie_m2 = 1e-9*area(field_polygon);

candidate_error_m2 = zeros(numel(candidate_degrees),1);
best = 0;
for k = 1:numel(candidate_degrees)
    [polygon_k,geometry_k] = Orient_Field_Polygon_Grid( ...
        field_polygon,cell_size_m,candidate_degrees(k));
    [offset_k,candidate_error_m2(k)] = lattice_offset(polygon_k,cell_size_m,tie_m2);
    if best == 0 || candidate_error_m2(k) < candidate_error_m2(best)-tie_m2
        best = k;
        oriented_polygon = polygon_k;
        geometry = geometry_k;
        offset_m = offset_k;
    end
end
geometry.orientation_objective = 'minimum_misrepresented_area';
geometry.orientation_candidate_count = numel(candidate_degrees);
geometry.grid_offset_m = offset_m;
geometry.candidate_rotation_degrees = candidate_degrees;
geometry.candidate_misrepresented_area_m2 = candidate_error_m2;

end

function [offset_m,misrepresented_m2] = lattice_offset(polygon,h,tie_m2)
%% The offset minimising E on one rotation: scan, then compass search.
field_area_m2 = area(polygon);
[xa,ya,xb,yb] = polygon_edges(polygon.Vertices);
moving = ya ~= yb;   %% an edge along a grid row adds nothing to the integral
xa = xa(moving); ya = ya(moving); xb = xb(moving); yb = yb(moving);
finite_vertex = all(isfinite(polygon.Vertices),2);
lower = min(polygon.Vertices(finite_vertex,:),[],1);
upper = max(polygon.Vertices(finite_vertex,:),[],1);

scan_count = 16;
best_offset = [0 0];
best_error = Inf;
for iy = 0:scan_count-1
    for ix = 0:scan_count-1
        offset = [ix iy]*h/scan_count;
        misrepresented = misrepresented_area(offset);
        if misrepresented < best_error-tie_m2
            best_offset = offset;
            best_error = misrepresented;
        end
    end
end
step = h/(2*scan_count);
moves = [1 0; -1 0; 0 1; 0 -1];
while step >= 1e-4*h
    move_error = Inf;
    move_offset = [];
    for k = 1:size(moves,1)
        offset = mod(best_offset+step*moves(k,:),h);
        misrepresented = misrepresented_area(offset);
        if misrepresented < move_error-tie_m2
            move_offset = offset;
            move_error = misrepresented;
        end
    end
    if move_error < best_error-tie_m2
        best_offset = move_offset;
        best_error = move_error;
    else
        step = step/2;
    end
end
offset_m = mod(round(best_offset*1e4)/1e4,h);
misrepresented_m2 = misrepresented_area(offset_m);

    function misrepresented = misrepresented_area(offset)
        cell_area_m2 = grid_overlap_areas(xa,ya,xb,yb,h, ...
            lower(1)-offset(1),lower(2)-offset(2), ...
            Grid_Cell_Count(upper(1)-lower(1)+offset(1),h), ...
            Grid_Cell_Count(upper(2)-lower(2)+offset(2),h));
        total_m2 = sum(cell_area_m2(:));
        if abs(total_m2+field_area_m2) < abs(total_m2-field_area_m2)
            cell_area_m2 = -cell_area_m2;   %% clockwise solid boundaries
            total_m2 = -total_m2;
        end
        if abs(total_m2-field_area_m2) > max(1e-6,1e-7*field_area_m2)
            error('Place_Field_Polygon_Grid:Orientation', ...
                ['The polygon boundaries are not consistently oriented ' ...
                 '(cell areas sum to %.17g m2, the polygon is %.17g m2).'], ...
                total_m2,field_area_m2);
        end
        cell_area_m2 = min(h^2,max(0,cell_area_m2(:)));
        misrepresented = sum(min(cell_area_m2,h^2-cell_area_m2));
    end

end

function area_m2 = grid_overlap_areas(xa,ya,xb,yb,h,x0,y0,nx,ny)
%% Signed overlap of every cell of the nx x ny lattice with corner (x0,y0):
%% each edge is cut at the grid lines into pieces that lie in one cell; a
%% piece adds dY*(mean X - x_i) to its own cell and h*dY to every cell left
%% of it in its row.
own = zeros(ny,nx);
right = zeros(ny,nx);
for e = 1:numel(xa)
    dx = xb(e)-xa(e);
    dy = yb(e)-ya(e);
    t = [0;1];
    if dx ~= 0
        k = ceil((min(xa(e),xb(e))-x0)/h):floor((max(xa(e),xb(e))-x0)/h);
        t = [t;((x0+k(:)*h)-xa(e))/dx]; %#ok<AGROW>
    end
    k = ceil((min(ya(e),yb(e))-y0)/h):floor((max(ya(e),yb(e))-y0)/h);
    t = [t;((y0+k(:)*h)-ya(e))/dy]; %#ok<AGROW>
    t = unique(min(1,max(0,t)));
    t_mid = 0.5*(t(1:end-1)+t(2:end));
    piece_dy = diff(t)*dy;
    x_mid = xa(e)+t_mid*dx;
    y_mid = ya(e)+t_mid*dy;
    i = min(nx,max(1,floor((x_mid-x0)/h)+1));
    j = min(ny,max(1,floor((y_mid-y0)/h)+1));
    slot = j+(i-1)*ny;
    own = own+reshape(accumarray(slot,piece_dy.*(x_mid-(x0+(i-1)*h)), ...
        [nx*ny 1]),ny,nx);
    right = right+reshape(accumarray(slot,h*piece_dy,[nx*ny 1]),ny,nx);
end
area_m2 = own+cumsum(right,2,'reverse')-right;
end

function [xa,ya,xb,yb] = polygon_edges(vertices)
%% Directed edges of every closed boundary (rings separated by NaN rows).
separator = any(~isfinite(vertices),2);
ring = cumsum(separator);
xa = zeros(0,1); ya = xa; xb = xa; yb = xa;
for r = reshape(unique(ring(~separator)),1,[])
    v = vertices(ring == r & ~separator,:);
    if size(v,1) < 3; continue; end
    w = v([2:end 1],:);
    xa = [xa;v(:,1)]; ya = [ya;v(:,2)]; %#ok<AGROW>
    xb = [xb;w(:,1)]; yb = [yb;w(:,2)]; %#ok<AGROW>
end
end
