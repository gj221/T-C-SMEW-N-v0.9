function [feature_polygons,metadata] = ...
        Read_Shapefile_Polygons(shape_file,feature_indices)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read selected shapefile polygons into one local metric projection.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2; feature_indices = []; end
if isstring(shape_file) && isscalar(shape_file); shape_file = char(shape_file); end
if ~(ischar(shape_file) && isrow(shape_file) && ~isempty(shape_file)) || ...
        ~isfile(shape_file)
    error('Read_Shapefile_Polygons:MissingShapefile', ...
        'shape_file must name an existing .shp file.');
end
if ~strcmpi(file_extension(shape_file),'.shp')
    error('Read_Shapefile_Polygons:BadFileType', ...
        'shape_file must have a .shp extension.');
end

info = shapeinfo(shape_file);
if ~startsWith(info.ShapeType,'Polygon','IgnoreCase',true)
    error('Read_Shapefile_Polygons:BadGeometry', ...
        'The shapefile must contain polygon geometry.');
end
source_crs = info.CoordinateReferenceSystem;
if isempty(source_crs)
    error('Read_Shapefile_Polygons:MissingCRS', ...
        'The shapefile must define its coordinate reference system.');
end

features = shaperead(shape_file);
feature_count = numel(features);
if feature_count == 0
    error('Read_Shapefile_Polygons:EmptyShapefile', ...
        'The shapefile contains no polygon features.');
end
if isempty(feature_indices)
    feature_indices = 1:feature_count;
elseif ~(isnumeric(feature_indices) && isvector(feature_indices) && ...
        all(isfinite(feature_indices)) && ...
        all(feature_indices == fix(feature_indices)) && ...
        all(feature_indices >= 1 & feature_indices <= feature_count) && ...
        numel(unique(feature_indices)) == numel(feature_indices))
    error('Read_Shapefile_Polygons:BadFeatureIndices', ...
        'feature_indices must contain unique shapefile record indices.');
end
feature_indices = reshape(feature_indices,1,[]);
features = features(feature_indices);

latitude = cell(numel(features),1);
longitude = cell(numel(features),1);
for k = 1:numel(features)
    if isa(source_crs,'geocrs')
        latitude{k} = features(k).Y;
        longitude{k} = features(k).X;
    elseif isa(source_crs,'projcrs')
        [latitude{k},longitude{k}] = projinv( ...
            source_crs,features(k).X,features(k).Y);
    else
        error('Read_Shapefile_Polygons:UnsupportedCRS', ...
            'The shapefile CRS must be geographic or projected.');
    end
end

all_latitude = [latitude{:}];
all_longitude = [longitude{:}];
valid = isfinite(all_latitude) & isfinite(all_longitude);
all_latitude = all_latitude(valid);
all_longitude = all_longitude(valid);
if isempty(all_latitude) || any(all_latitude < -80 | all_latitude > 84) || ...
        range(all_longitude) >= 6
    error('Read_Shapefile_Polygons:UnsupportedExtent', ...
        'The selected polygons must fit inside one local UTM projection.');
end

latitude_center = mean(all_latitude);
longitude_center = mean(all_longitude);
utm_zone = min(60,max(1,floor((longitude_center+180)/6)+1));
if latitude_center >= 0
    projected_crs_epsg = 32600+utm_zone;
else
    projected_crs_epsg = 32700+utm_zone;
end
target_crs = projcrs(projected_crs_epsg);

feature_polygons = cell(numel(features),1);
source_finite_vertex_count = zeros(1,numel(features));
projected_polygon_vertex_count = zeros(1,numel(features));
projection_roundtrip_max_error_m = zeros(1,numel(features));
source_to_polygon_max_error_m = zeros(1,numel(features));
for k = 1:numel(features)
    [x_m,y_m] = projfwd(target_crs,latitude{k},longitude{k});
    valid = isfinite(x_m) & isfinite(y_m);
    feature_polygons{k} = polyshape(x_m,y_m,'Simplify',false, ...
        'KeepCollinearPoints',true);
    if isempty(feature_polygons{k}.Vertices) || area(feature_polygons{k}) <= 0
        error('Read_Shapefile_Polygons:InvalidPolygon', ...
            'Each selected feature must define finite positive-area geometry.');
    end

    [latitude_rt,longitude_rt] = projinv( ...
        target_crs,x_m(valid),y_m(valid));
    [x_rt_m,y_rt_m] = projfwd( ...
        target_crs,latitude_rt,longitude_rt);
    projection_roundtrip_max_error_m(k) = max(hypot( ...
        x_rt_m-x_m(valid),y_rt_m-y_m(valid)));

    polygon_vertices = feature_polygons{k}.Vertices;
    finite_polygon_vertex = all(isfinite(polygon_vertices),2);
    polygon_vertices = polygon_vertices(finite_polygon_vertex,:);
    source_finite_vertex_count(k) = sum(valid);
    projected_polygon_vertex_count(k) = size(polygon_vertices,1);
    source_index = find(valid);
    for j = 1:numel(source_index)
        i = source_index(j);
        source_to_polygon_max_error_m(k) = max( ...
            source_to_polygon_max_error_m(k),min(hypot( ...
            polygon_vertices(:,1)-x_m(i), ...
            polygon_vertices(:,2)-y_m(i))));
    end
    if projection_roundtrip_max_error_m(k) >= 1e-6 || ...
            source_to_polygon_max_error_m(k) >= 1e-9
        error('Read_Shapefile_Polygons:GeometryFidelity', ...
            'Projection or polygon construction changed a source boundary.');
    end
end

metadata.source_shapefile = shape_file;
metadata.source_crs_name = char(source_crs.Name);
metadata.projected_crs_name = char(target_crs.Name);
metadata.projected_crs_epsg = projected_crs_epsg;
metadata.shapefile_feature_count = feature_count;
metadata.feature_indices = feature_indices;
metadata.source_finite_vertex_count = source_finite_vertex_count;
metadata.projected_polygon_vertex_count = projected_polygon_vertex_count;
metadata.projection_roundtrip_max_error_m = ...
    projection_roundtrip_max_error_m;
metadata.source_to_polygon_max_error_m = source_to_polygon_max_error_m;

end

function extension = file_extension(file_path)
[~,~,extension] = fileparts(file_path);
end
