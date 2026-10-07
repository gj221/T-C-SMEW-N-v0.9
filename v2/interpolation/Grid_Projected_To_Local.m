function [x_grid_m,y_grid_m] = Grid_Projected_To_Local(x_projected_m,y_projected_m,topo)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Transform projected CRS coordinates into local rotated-grid coordinates.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isequal(size(x_projected_m),size(y_projected_m))
    error('Grid_Projected_To_Local:SizeMismatch', ...
        'x_projected_m and y_projected_m must have the same size.');
end
if ~isstruct(topo) || ~isfield(topo,'grid_rotation_degrees') || ...
        ~isfield(topo,'projected_grid_origin_m') || ...
        numel(topo.projected_grid_origin_m) ~= 2
    error('Grid_Projected_To_Local:BadTopology', ...
        'The topology must contain its grid rotation and projected origin.');
end

angle = deg2rad(topo.grid_rotation_degrees);
c = cos(angle);
s = sin(angle);
origin = topo.projected_grid_origin_m;
x_offset = x_projected_m-origin(1);
y_offset = y_projected_m-origin(2);
x_grid_m = x_offset*c+y_offset*s;
y_grid_m = -x_offset*s+y_offset*c;

end
