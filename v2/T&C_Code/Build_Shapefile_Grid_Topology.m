function topo = Build_Shapefile_Grid_Topology( ...
        shape_file,cell_size_m,feature_indices)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read selected polygon records and build one grid over their union.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2; cell_size_m = []; end
if nargin < 3; feature_indices = []; end
[feature_polygons,metadata] = ...
    Read_Shapefile_Polygons(shape_file,feature_indices);

field_polygon = polyshape();
for k = 1:numel(feature_polygons)
    field_polygon = union(field_polygon,feature_polygons{k});
end
if isempty(field_polygon.Vertices) || area(field_polygon) <= 0
    error('Build_Shapefile_Grid_Topology:EmptyField', ...
        'The projected polygon union has no positive area.');
end

topo = Build_Projected_Field_Topology( ...
    field_polygon,cell_size_m,metadata,metadata.feature_indices);

end
