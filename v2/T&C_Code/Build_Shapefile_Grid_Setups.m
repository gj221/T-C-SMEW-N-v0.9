function setups = Build_Shapefile_Grid_Setups( ...
        shape_file,cell_size_m,feature_indices,coverage_rule)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Build one independent full-grid setup per connected polygon region.
%% coverage_rule: a fraction (default 0.95) or 'majority'
%% (Build_Projected_Field_Topology).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2; cell_size_m = []; end
if nargin < 3; feature_indices = []; end
if nargin < 4; coverage_rule = []; end
cache_key = setup_cache_key(shape_file,cell_size_m,feature_indices, ...
    coverage_rule);
persistent cache_keys cache_values
if ~isempty(cache_keys)
    cache_index = find(strcmp(cache_keys,cache_key),1);
    if ~isempty(cache_index)
        setups = cache_values{cache_index};
        return
    end
end
[feature_polygons,metadata] = ...
    Read_Shapefile_Polygons(shape_file,feature_indices);
[field_regions,region_feature_indices] = connected_field_regions( ...
    feature_polygons,metadata.feature_indices);
if isempty(field_regions)
    error('Build_Shapefile_Grid_Setups:EmptyField', ...
        'The selected polygon union has no positive-area regions.');
end

setup_cell = cell(numel(field_regions),1);
sort_key = zeros(numel(field_regions),3);
for k = 1:numel(field_regions)
    region = field_regions{k};
    source_feature_indices = region_feature_indices{k};
    setup_cell{k} = Build_Projected_Field_Topology( ...
        region,cell_size_m,metadata,source_feature_indices, ...
        coverage_rule);
    [center_x,center_y] = centroid(region);
    sort_key(k,:) = [min(source_feature_indices) center_x center_y];
end

[~,order] = sortrows(sort_key,[1 2 3]);
setups = vertcat(setup_cell{order});
for k = 1:numel(setups)
    setups(k).setup_index = k;
    setups(k).setup_count = numel(setups);
end
cache_keys{end+1} = cache_key;
cache_values{end+1} = setups;
if numel(cache_keys) > 8
    cache_keys = cache_keys(end-7:end);
    cache_values = cache_values(end-7:end);
end

end

function key = setup_cache_key(shape_file,cell_size_m,feature_indices, ...
    coverage_rule)
[shape_directory,shape_name] = fileparts(shape_file);
components = cellfun(@(extension) fullfile(shape_directory, ...
    [shape_name extension]),{'.shp','.shx','.dbf','.prj'}, ...
    'UniformOutput',false);
runtime_identity = sprintf('cell_size=%s;features=%s;coverage=%s', ...
    mat2str(cell_size_m,17),mat2str(feature_indices,17), ...
    mat2str(coverage_rule,17));
key = ERW_Checkpoint_Fingerprint(0,0,runtime_identity,components{:});
end

function [field_regions,region_feature_indices] = connected_field_regions( ...
        feature_polygons,feature_indices)
candidate_regions = cell(0,1);
candidate_feature_indices = zeros(0,1);
for k = 1:numel(feature_polygons)
    feature_regions = regions(feature_polygons{k});
    for j = 1:numel(feature_regions)
        candidate_regions{end+1,1} = feature_regions(j); %#ok<AGROW>
        candidate_feature_indices(end+1,1) = feature_indices(k); %#ok<AGROW>
    end
end

region_count = numel(candidate_regions);
adjacency = false(region_count);
bounds = zeros(region_count,4);
for k = 1:region_count
    vertices = candidate_regions{k}.Vertices;
    vertices = vertices(all(isfinite(vertices),2),:);
    bounds(k,:) = [min(vertices(:,1)) max(vertices(:,1)) ...
        min(vertices(:,2)) max(vertices(:,2))];
end
for i = 1:region_count-1
    for j = i+1:region_count
        boxes_overlap = bounds(i,1) <= bounds(j,2) && ...
            bounds(j,1) <= bounds(i,2) && ...
            bounds(i,3) <= bounds(j,4) && ...
            bounds(j,3) <= bounds(i,4);
        if ~boxes_overlap; continue; end
        adjacency(i,j) = isscalar(regions( ...
            union(candidate_regions{i},candidate_regions{j})));
        adjacency(j,i) = adjacency(i,j);
    end
end
component = conncomp(graph(adjacency));

field_regions = cell(max(component),1);
region_feature_indices = cell(max(component),1);
for k = 1:max(component)
    member = find(component == k);
    field_regions{k} = candidate_regions{member(1)};
    for j = reshape(member(2:end),1,[])
        field_regions{k} = union(field_regions{k},candidate_regions{j});
    end
    region_feature_indices{k} = unique( ...
        candidate_feature_indices(member).','stable');
end
end
