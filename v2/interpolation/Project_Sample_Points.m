function [samples,report] = Project_Sample_Points(points,field_setup, ...
        x_column,y_column)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Project geographic samples into a field setup and keep those inside it.
%%
%% The samples are carried into the UTM zone the setup already resolved and
%% then into its rotated grid frame, so a sample and a grid cell can be
%% compared without a second projection. Points on the boundary are dropped,
%% matching the containment test used on the Python side.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3 || isempty(x_column); x_column = 'lng (X)'; end
if nargin < 4 || isempty(y_column); y_column = 'lat (Y)'; end
if ~istable(points) || isempty(points)
    error('Project_Sample_Points:BadPoints', ...
        'points must be a non-empty table of soil samples.');
end
if ~isstruct(field_setup) || ~isscalar(field_setup) || ...
        ~isfield(field_setup,'projected_crs_epsg')
    error('Project_Sample_Points:BadFieldSetup', ...
        'field_setup must be one setup from Build_Shapefile_Grid_Setups.');
end

longitude = Numeric_Sample_Column(points,x_column);
latitude = Numeric_Sample_Column(points,y_column);
target_crs = projcrs(field_setup.projected_crs_epsg);
[x_projected_m,y_projected_m] = projfwd(target_crs,latitude,longitude);
[x_grid_m,y_grid_m] = Grid_Projected_To_Local( ...
    x_projected_m,y_projected_m,field_setup);

boundary_polygon = Field_Setup_Polygon(field_setup);
inside = isinterior(boundary_polygon,x_grid_m,y_grid_m);

samples = points;
samples.x_projected_m = x_projected_m;
samples.y_projected_m = y_projected_m;
samples.x_grid_m = x_grid_m;
samples.y_grid_m = y_grid_m;
samples = samples(inside,:);

report.projected_crs_epsg = field_setup.projected_crs_epsg;
report.candidate_count = height(points);
report.inside_count = height(samples);

end
