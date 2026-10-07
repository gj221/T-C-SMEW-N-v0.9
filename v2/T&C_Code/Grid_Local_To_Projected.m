function [x_projected_m,y_projected_m] = Grid_Local_To_Projected(x_grid_m,y_grid_m,topo)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Transform local rotated-grid coordinates back to the projected CRS.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isequal(size(x_grid_m),size(y_grid_m))
    error('Grid_Local_To_Projected:SizeMismatch', ...
        'x_grid_m and y_grid_m must have the same size.');
end
if ~isstruct(topo) || ~isfield(topo,'grid_rotation_degrees') || ...
        ~isfield(topo,'projected_grid_origin_m') || ...
        numel(topo.projected_grid_origin_m) ~= 2
    error('Grid_Local_To_Projected:BadTopology', ...
        'The topology must contain its grid rotation and projected origin.');
end

angle = deg2rad(topo.grid_rotation_degrees);
c = cos(angle);
s = sin(angle);
origin = topo.projected_grid_origin_m;
x_projected_m = origin(1)+x_grid_m*c-y_grid_m*s;
y_projected_m = origin(2)+x_grid_m*s+y_grid_m*c;

end
