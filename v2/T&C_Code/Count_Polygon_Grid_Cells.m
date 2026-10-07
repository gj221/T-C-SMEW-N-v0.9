function summary = Count_Polygon_Grid_Cells(field_polygon,cell_size_m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Count full square cells with positive polygon overlap during calibration.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vertices = field_polygon.Vertices;
vertices = vertices(all(isfinite(vertices),2),:);
x_min = min(vertices(:,1));
x_max = max(vertices(:,1));
y_min = min(vertices(:,2));
y_max = max(vertices(:,2));
nx = Grid_Cell_Count(x_max-x_min,cell_size_m);
ny = Grid_Cell_Count(y_max-y_min,cell_size_m);
[x_center,y_center] = meshgrid( ...
    x_min+((1:nx)-0.5)*cell_size_m, ...
    y_min+((1:ny)-0.5)*cell_size_m);
x_center = x_center(:);
y_center = y_center(:);
active = isinterior(field_polygon,x_center,y_center);

expanded = polybuffer(field_polygon,cell_size_m/sqrt(2)*(1+1e-9));
boundary_candidate = find(~active & ...
    isinterior(expanded,x_center,y_center));
for k = reshape(boundary_candidate,1,[])
    square = polyshape( ...
        x_center(k)+0.5*cell_size_m*[-1 1 1 -1], ...
        y_center(k)+0.5*cell_size_m*[-1 -1 1 1]);
    active(k) = area(intersect(field_polygon,square)) > 0;
end

summary.N_cells = sum(active);
summary.main_grid_nx = nx;
summary.main_grid_ny = ny;
summary.grid_slot_count = nx*ny;
summary.field_area_m2 = area(field_polygon);
summary.simulated_area_m2 = summary.N_cells*cell_size_m^2;
end
