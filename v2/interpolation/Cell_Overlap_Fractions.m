function [overlap_fraction,nx,ny] = Cell_Overlap_Fractions( ...
        field_polygon,cell_size_m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fraction of each grid square that lies inside the polygon.
%%
%% A square of width h has circumradius h/sqrt(2), so a cell whose centre
%% sits that far inside the polygon is wholly inside it and a cell whose
%% centre sits that far outside is wholly outside. Only the band between the
%% two buffers needs an exact intersection, which keeps a fine sweep cheap
%% without approximating any boundary cell.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~(isa(field_polygon,'polyshape') && isscalar(field_polygon)) || ...
        area(field_polygon) <= 0
    error('Cell_Overlap_Fractions:EmptyField', ...
        'field_polygon must have finite positive area.');
end
if ~(isnumeric(cell_size_m) && isscalar(cell_size_m) && ...
        isfinite(cell_size_m) && cell_size_m > 0)
    error('Cell_Overlap_Fractions:BadCellSize', ...
        'cell_size_m must be a finite positive scalar in metres.');
end

vertices = field_polygon.Vertices;
vertices = vertices(all(isfinite(vertices),2),:);
x_min = min(vertices(:,1));
y_min = min(vertices(:,2));
nx = Grid_Cell_Count(max(vertices(:,1))-x_min,cell_size_m);
ny = Grid_Cell_Count(max(vertices(:,2))-y_min,cell_size_m);
[x_center,y_center] = meshgrid(x_min+((1:nx)-0.5)*cell_size_m, ...
    y_min+((1:ny)-0.5)*cell_size_m);
x_center = x_center(:);
y_center = y_center(:);

circumradius_m = cell_size_m/sqrt(2)*(1+1e-9);
eroded = polybuffer(field_polygon,-circumradius_m);
expanded = polybuffer(field_polygon,circumradius_m);
overlap_fraction = zeros(numel(x_center),1);
fully_inside = Isinterior_In_Batches(eroded,x_center,y_center);
overlap_fraction(fully_inside) = 1;

% Only the cells the first pass left undecided need the second, wider test,
% so the point count checked against the polygon boundary shrinks with the
% fraction of the grid that is unambiguous, not with the grid's own size.
outside_eroded = find(~fully_inside);
in_expanded = Isinterior_In_Batches(expanded, ...
    x_center(outside_eroded),y_center(outside_eroded));
uncertain = outside_eroded(in_expanded);

cell_area_m2 = cell_size_m^2;
for k = reshape(uncertain,1,[])
    square = polyshape( ...
        x_center(k)+0.5*cell_size_m*[-1 1 1 -1], ...
        y_center(k)+0.5*cell_size_m*[-1 -1 1 1]);
    overlap_fraction(k) = ...
        area(intersect(field_polygon,square))/cell_area_m2;
end

end
