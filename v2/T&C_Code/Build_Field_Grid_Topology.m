function topo = Build_Field_Grid_Topology(shape_file,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Build_Field_Grid_Topology
%%
%% The cell grid of one sub-field run on a real field: T&C-SMEW-N v1.0's
%% polygon grid on a field shapefile (Build_Shapefile_Grid_Setups), with
%% this tree's margin band for the margin layouts, and the DEM attached
%% when one is given (Load_Field_DEM, Attach_Topography_To_Topology).
%%
%% Run_One_Subfield calls it for a FieldShapefile run. Saved,
%%     topo = Build_Field_Grid_Topology(...);
%%     save(file,'topo');
%% it is the FieldTopologyFile a run loads instead, which needs no Mapping
%% Toolbox (this function does: shaperead, projcrs, readgeoraster) and
%% gives every run of a campaign the identical grid.
%%
%% Name/value options
%%   'GridCellSizeM'   cell width [m]; the majority rule takes it from the
%%                     terrain (Run_Topographic_Grid_Convergence)
%%   'CoverageRule'    'majority' (default) or a fraction in (0,1), v1.0's
%%                     coverage calibration (whole field only)
%%   'MarginWidthM'    margin band around the field [m]; 0 (default) is
%%                     v1.0's field grid, every cell 'main'
%%   'FeatureIndices'  shapefile features to use (default all)
%%   'PolygonIndex'    connected region to grid (default: the only one)
%%   'FieldDEM'        GeoTIFF, directory of tiles or {paths}; empty = flat
%%   'DEMResolutionM'  DEM spacing [m] (default the cell width, as v1.0:
%%                     attributes are sampled at cell centres, and a finer
%%                     DEM aliases sub-cell microtopography into them)
%%
%% topo.field_grid_build records the inputs, with a SHA-256 of each
%% shapefile component and of each DEM tile.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

opt = struct('GridCellSizeM',[],'CoverageRule','majority', ...
    'MarginWidthM',0,'FeatureIndices',[],'PolygonIndex',[], ...
    'FieldDEM','','DEMResolutionM',[]);
if mod(numel(varargin),2) ~= 0
    error('Build_Field_Grid_Topology:BadOptions', ...
        'Options must be name/value pairs.');
end
for k = 1:2:numel(varargin)
    name = char(varargin{k});
    if ~isfield(opt,name)
        error('Build_Field_Grid_Topology:BadOptions','Unknown option %s.',name);
    end
    opt.(name) = varargin{k+1};
end
if ~(ischar(shape_file) || (isstring(shape_file) && isscalar(shape_file))) || ...
        ~isfile(char(shape_file))
    error('Build_Field_Grid_Topology:MissingShapefile', ...
        'Field shapefile not found: %s',char(string(shape_file)));
end
[~,shape_attributes] = fileattrib(char(shape_file));
shape_file = shape_attributes.Name;
if isempty(opt.CoverageRule); opt.CoverageRule = 'majority'; end
if isstring(opt.CoverageRule); opt.CoverageRule = char(opt.CoverageRule); end
if isempty(opt.MarginWidthM); opt.MarginWidthM = 0; end
if isstring(opt.FieldDEM); opt.FieldDEM = cellstr(opt.FieldDEM); end
if iscell(opt.FieldDEM) && isscalar(opt.FieldDEM); opt.FieldDEM = opt.FieldDEM{1}; end

setups = Build_Shapefile_Grid_Setups(shape_file,opt.GridCellSizeM, ...
    opt.FeatureIndices,opt.CoverageRule,opt.MarginWidthM);
polygon_index = opt.PolygonIndex;
if isempty(polygon_index)
    if numel(setups) ~= 1
        error('Build_Field_Grid_Topology:DisconnectedField', ...
            ['The shapefile holds %d disconnected polygons; choose one ' ...
             'with PolygonIndex.'],numel(setups));
    end
    polygon_index = 1;
elseif ~(isnumeric(polygon_index) && isscalar(polygon_index) && ...
        polygon_index == fix(polygon_index) && polygon_index >= 1 && ...
        polygon_index <= numel(setups))
    error('Build_Field_Grid_Topology:BadPolygonIndex', ...
        'PolygonIndex must be between 1 and %d.',numel(setups));
end
topo = setups(polygon_index);
if topo.N_cells < 2
    error('Build_Field_Grid_Topology:SingleCell', ...
        ['The grid has one cell: a sub-field run needs a lateral graph. ' ...
         'Use a smaller GridCellSizeM.']);
end

dem_files = {};
dem_resolution_m = [];
if ~isempty(opt.FieldDEM)
    dem_resolution_m = opt.DEMResolutionM;
    if isempty(dem_resolution_m); dem_resolution_m = topo.cell_size_m; end
    dem = Load_Field_DEM(opt.FieldDEM,topo,'OutResolutionM',dem_resolution_m);
    topo = Attach_Topography_To_Topology(topo,dem);
    dem_files = reshape(cellstr(dem.source_files),[],1);
end

shape_components = {};
[shape_directory,shape_name] = fileparts(shape_file);
for extension = {'.shp','.shx','.dbf','.prj','.cpg'}
    component = fullfile(shape_directory,[shape_name extension{1}]);
    if isfile(component); shape_components{end+1,1} = component; end %#ok<AGROW>
end
topo.field_grid_build = struct( ...
    'shapefile',shape_file, ...
    'shapefile_sha256',{file_digests(shape_components)}, ...
    'shapefile_components',{shape_components}, ...
    'feature_indices',opt.FeatureIndices, ...
    'polygon_index',polygon_index, ...
    'cell_size_m',topo.cell_size_m, ...
    'coverage_rule',opt.CoverageRule, ...
    'margin_width_m',opt.MarginWidthM, ...
    'dem_files',{dem_files}, ...
    'dem_sha256',{file_digests(dem_files)}, ...
    'dem_resolution_m',dem_resolution_m, ...
    'matlab',version, ...
    'built',char(datetime('now','TimeZone','UTC','Format','yyyy-MM-dd''T''HH:mm:ss''Z''')));

end

function digests = file_digests(files)
%% SHA-256 of each file, read in 64 MB blocks: a 1 m lidar tile is
%% hundreds of MB.
digests = cell(numel(files),1);
for k = 1:numel(files)
    fid = fopen(files{k},'rb');
    if fid < 0
        error('Build_Field_Grid_Topology:FileRead','Cannot read %s',files{k});
    end
    md = javaMethod('getInstance','java.security.MessageDigest','SHA-256');
    while true
        bytes = fread(fid,64*2^20,'*uint8');
        if isempty(bytes); break; end
        md.update(typecast(bytes,'int8'));
    end
    fclose(fid);
    raw = typecast(md.digest(),'uint8');
    digests{k} = lower(reshape(dec2hex(raw,2).',1,[]));
end
end
