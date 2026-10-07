function dem = Load_Field_DEM(source, topo_or_crs, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load_Field_DEM
%%
%% Read one or more projected elevation GeoTIFFs, mosaic them onto a common
%% grid, resample to a chosen spacing, clip to the field, and return a
%% plain gridded elevation model for Field_Topographic_Attributes /
%% Attach_Topography_To_Topology.
%%
%% The DEM must be supplied in a PROJECTED CRS (UTM).  It should match the
%% field topology's CRS; a same-zone datum difference (NAD83 vs WGS84,
%% ~1 m in CONUS) is tolerated with a warning.  Export the tile from the
%% DEM portal in UTM zone 17N for the Bruce Farm / Adam Wilbourne site.
%%
%% INPUTS
%%   source        a GeoTIFF path, a directory of GeoTIFF tiles
%%                 (*.tif / *.tiff), or a {cell} of paths.
%%   topo_or_crs   the field topology struct (uses .projected_crs_epsg and
%%                 the projected cell-centre bounds), an EPSG code, a
%%                 projcrs, or [] to keep the source extent and CRS.
%%
%% Name/value options
%%   'OutResolutionM'  output spacing [m] (default: native, snapped to
%%                     {0.5 1 2 3 5 10 30}).  A coarser output is the
%%                     area mean of the source cells it covers (box
%%                     filter), not a point sample.
%%   'BufferM'         margin kept around the field [m] (default 100).
%%   'BBoxLatLon'      [lonW lonE latS latN] clip window (WGS84 degrees).
%%   'FillMaxCells'    largest NaN void filled by interpolation (default 50).
%%
%% OUTPUT
%%   dem   struct: .Z [ny x nx] elevation [m]; .x [1 x nx], .y [1 x ny]
%%         ascending projected coordinates [m]; .cellsize_m;
%%         .native_resolution_m; .crs_name; .source_files.
%%
%% Requires the Mapping Toolbox (readgeoraster).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

opt = struct('OutResolutionM',[],'BufferM',100,'BBoxLatLon',[], ...
    'FillMaxCells',50);
for k = 1:2:numel(varargin); opt.(varargin{k}) = varargin{k+1}; end
if nargin < 2; topo_or_crs = []; end

files = resolve_source_files(source);
if isempty(files)
    error('Load_Field_DEM:NoFiles','No GeoTIFF found at %s.',string(source));
end

%% --- the crop window, before any tile is read -----------------------
%% A 1 m lidar tile is ~10^8 samples: six of them held as double is
%% ~5 GB, which exhausted the workstation and would exhaust an 8 GB
%% cluster task.  The clip window is known from the target before any
%% data is read (only the first tile's georeference is needed for the
%% source CRS), so every tile is read in its native single precision and
%% cropped to that window plus a margin before it is kept.  Peak memory
%% is then one tile.
info1 = georasterinfo(files{1});
R1 = info1.RasterReference;
if ~(isprop(R1,'ProjectedCRS') && ~isempty(R1.ProjectedCRS))
    error('Load_Field_DEM:GeographicSource', ...
        ['%s is not in a projected CRS. Re-export the DEM from the ' ...
         'portal in a projected (UTM) CRS.'], files{1});
end
src_crs = R1.ProjectedCRS;
[tgt_crs,tgt_epsg,bounds_xy] = resolve_target(topo_or_crs,src_crs,opt);
reproject = ~same_crs(src_crs,tgt_crs);
crop_margin_m = 100;                       % box-filter offsets + reprojection slack
cx = bounds_xy([1 2 2 1]); cy = bounds_xy([3 3 4 4]);
if reproject
    [clat,clon] = projinv(tgt_crs,cx,cy);
    [cx,cy] = projfwd(src_crs,clat,clon);
end
crop_win = [min(cx)-crop_margin_m, max(cx)+crop_margin_m, ...
            min(cy)-crop_margin_m, max(cy)+crop_margin_m];

%% --- read tiles onto ascending-(x,y) grids, cropped -------------------
tiles = cell(numel(files),1);
for k = 1:numel(files)
    [Z,R] = readgeoraster(files{k},'OutputType','single');
    %% readgeoraster does not apply the file's own NoData value, so read
    %% it rather than guessing: -999, -99 and 0 are all in use and would
    %% otherwise enter the interpolant as real elevations.  The sentinel
    %% test stays as a fallback for files that declare nothing.
    have_nd = false;
    try
        nd = georasterinfo(files{k}).MissingDataIndicator;
        if ~isempty(nd); Z = standardizeMissing(Z,nd); have_nd = true; end
    catch
    end
    Z(Z <= -9990 | Z >= 9000) = NaN;              % -9999, -32768, -3.4e38
    if ~have_nd
        %% No declared indicator: also trap the smaller sentinels some
        %% exports use.  0 is left alone (it can be a real elevation), but
        %% is counted so a tile of zeros is not silently believed.
        Z(ismember(Z,[-999 -99])) = NaN;
        if mean(Z(:) == 0) > 0.05
            warning('Load_Field_DEM:ManyZeros', ...
                ['%s declares no NoData value and %.0f%% of its cells are ' ...
                 'exactly 0; if 0 is its NoData sentinel, set them to NaN ' ...
                 'before use.'],files{k},100*mean(Z(:) == 0));
        end
    end
    if ~(isprop(R,'ProjectedCRS') && ~isempty(R.ProjectedCRS))
        error('Load_Field_DEM:GeographicSource', ...
            ['%s is not in a projected CRS. Re-export the DEM from the ' ...
             'portal in a projected (UTM) CRS.'], files{k});
    end
    if ~same_crs(src_crs,R.ProjectedCRS)
        error('Load_Field_DEM:MixedSourceCRS', ...
            '%s is in %s but %s is in %s; the tiles must share one CRS.', ...
            files{k},char(R.ProjectedCRS.Name),files{1},char(src_crs.Name));
    end
    nx = R.RasterSize(2); ny = R.RasterSize(1);
    %% Ascending sample coordinates.  A cell-registered raster stores
    %% pixel CENTRES half a cell inside its world limits; a point
    %% (postings) raster stores the samples ON the limits and reports
    %% SampleSpacingInWorld*, so the half-cell offset must not be applied.
    if isprop(R,'RasterInterpretation') && strcmp(R.RasterInterpretation,'postings')
        dx = R.SampleSpacingInWorldX; dy = R.SampleSpacingInWorldY;
        xs = R.XWorldLimits(1) + (0:nx-1)*dx;
        ys = R.YWorldLimits(1) + (0:ny-1)*dy;
    else
        dx = R.CellExtentInWorldX; dy = R.CellExtentInWorldY;
        xs = R.XWorldLimits(1) + ((1:nx)-0.5)*dx;
        ys = R.YWorldLimits(1) + ((1:ny)-0.5)*dy;
    end
    %% Z rows run north->south when ColumnsStartFrom is 'north', so flip
    %% Z (not ys) to match ascending y.
    if ~isprop(R,'ColumnsStartFrom') || strcmp(R.ColumnsStartFrom,'north')
        Z = flipud(Z);
    end
    if isprop(R,'RowsStartFrom') && strcmp(R.RowsStartFrom,'east')
        Z = fliplr(Z);
    end
    ix = xs >= crop_win(1) & xs <= crop_win(2);
    iy = ys >= crop_win(3) & ys <= crop_win(4);
    %% A tile that keeps fewer than two samples in either direction cannot
    %% be interpolated, and cannot bracket any query point either: the
    %% queries reach at most out_res/2 + native_res beyond the field
    %% window, well inside the crop margin.
    if nnz(ix) < 2 || nnz(iy) < 2
        tiles{k} = [];                        % this tile lies outside the window
        clear Z
        continue
    end
    Z = double(Z(iy,ix)); xs = xs(ix); ys = ys(iy);
    tiles{k} = struct('Z',Z,'x',xs,'y',ys,'res',mean([dx dy]));
    clear Z
end
tiles = tiles(~cellfun(@isempty,tiles));
if isempty(tiles)
    error('Load_Field_DEM:NoOverlap', ...
        ['No DEM tile overlaps the requested window (source x=[%.0f %.0f] ' ...
         'y=[%.0f %.0f]); fetch the tile(s) that cover this field.'], ...
        crop_win(1),crop_win(2),crop_win(3),crop_win(4));
end
native_res = tiles{1}.res;
out_res = opt.OutResolutionM;
if isempty(out_res); out_res = snap_resolution(native_res); end

%% --- target CRS + clip window (in the TARGET CRS): resolved above ------
%% projinv/projfwd carry no datum transformation, so a same-zone datum
%% difference (the USGS tiles are NAD83, the field topology WGS84) passes
%% through as an identity and its ~1-1.5 m offset in CONUS is absorbed as
%% a registration error.  That is tolerable for a field-scale DEM but
%% must not be silent -- the header promises this warning.
if reproject
    da = datum_name(src_crs); db = datum_name(tgt_crs);
    if ~isempty(da) && ~isempty(db) && ~strcmp(da,db)
        warning('Load_Field_DEM:DatumNotTransformed', ...
            ['DEM datum "%s" differs from the target datum "%s"; ' ...
             'projinv/projfwd apply no datum shift, so the DEM keeps a ' ...
             'residual horizontal offset (about 1-1.5 m in CONUS for ' ...
             'NAD83 vs WGS84).'], da, db);
    end
end

xt = (bounds_xy(1)+out_res/2) : out_res : bounds_xy(2);
yt = (bounds_xy(3)+out_res/2) : out_res : bounds_xy(4);
if numel(xt) < 3 || numel(yt) < 3
    error('Load_Field_DEM:TinyExtent', ...
        'Clip window is %g x %g m; check BBoxLatLon / topology bounds.', ...
        diff(bounds_xy(1:2)),diff(bounds_xy(3:4)));
end
[XT,YT] = meshgrid(xt,yt);

%% --- sample tiles, average overlaps ----------------------------
%% When the output is coarser than the source, each output cell is the
%% AREA MEAN of an r x r stencil of source samples inside it (r = the
%% resolution ratio), not a single point sample: point-sampling a 1 m
%% lidar surface onto a 5 m grid aliases the microtopography into the
%% slope / curvature / TWI statistics the flow code depends on.
r = max(1,round(out_res/native_res));
offs = ((1:r)-(r+1)/2)*(out_res/r);
acc = zeros(size(XT)); cnt = zeros(size(XT));
interp = cellfun(@(t) griddedInterpolant({t.y,t.x},t.Z,'linear','none'), ...
    tiles,'UniformOutput',false);
for oi = 1:r
    for oj = 1:r
        Xq = XT+offs(oj); Yq = YT+offs(oi);
        %% target grid points -> source-CRS coordinates for sampling
        if reproject
            [latq,lonq] = projinv(tgt_crs,Xq,Yq);
            [XS,YS]     = projfwd(src_crs,latq,lonq);
        else
            XS = Xq; YS = Yq;
        end
        for k = 1:numel(tiles)
            v = interp{k}(YS,XS);
            ok = isfinite(v);
            acc(ok) = acc(ok) + v(ok);
            cnt(ok) = cnt(ok) + 1;
        end
    end
end
Z = acc./cnt;
Z(cnt == 0) = NaN;
if r == 1 && reproject
    %% keep the source-frame extent diagnostic meaningful below
    [latq,lonq] = projinv(tgt_crs,XT,YT); [XS,YS] = projfwd(src_crs,latq,lonq);
elseif r == 1
    XS = XT; YS = YT;
end
if all(isnan(Z(:)))
    error('Load_Field_DEM:NoOverlap', ...
        ['DEM tiles do not cover the requested area. Tile 1 spans ' ...
         'x=[%.0f %.0f] y=[%.0f %.0f] in its own CRS; the field/box maps ' ...
         'to source x=[%.0f %.0f] y=[%.0f %.0f]. If those disagree the DEM ' ...
         'does not cover this field -- fetch the tile(s) that do.'], ...
        min(tiles{1}.x),max(tiles{1}.x),min(tiles{1}.y),max(tiles{1}.y), ...
        min(XS(:)),max(XS(:)),min(YS(:)),max(YS(:)));
end
Z = fill_voids(Z,opt.FillMaxCells);

dem = struct('Z',Z,'x',xt,'y',yt,'cellsize_m',out_res, ...
    'native_resolution_m',native_res,'epsg',tgt_epsg, ...
    'crs_name',char(tgt_crs.Name),'source_crs_name',char(src_crs.Name), ...
    'reprojected',reproject,'source_files',{files});
fprintf(['Load_Field_DEM: %d tile(s) -> %d x %d at %.2f m, %.1f%% NaN ' ...
    '[%s%s]\n'], numel(files),size(Z,1),size(Z,2),out_res, ...
    100*mean(isnan(Z(:))),char(tgt_crs.Name), ...
    ternary(reproject,sprintf(' <- %s',char(src_crs.Name)),''));
end

%% =======================================================================
function [tgt_crs,tgt_epsg,bounds_xy] = resolve_target(topo_or_crs,src_crs,opt)
tgt_epsg = [];
if isempty(topo_or_crs)
    tgt_crs = src_crs;
elseif isa(topo_or_crs,'projcrs')
    tgt_crs = topo_or_crs;
elseif isnumeric(topo_or_crs) && isscalar(topo_or_crs)
    tgt_epsg = topo_or_crs; tgt_crs = projcrs(tgt_epsg);
elseif isstruct(topo_or_crs) && isfield(topo_or_crs,'projected_crs_epsg')
    tgt_epsg = topo_or_crs.projected_crs_epsg; tgt_crs = projcrs(tgt_epsg);
else
    error('Load_Field_DEM:BadTarget', ...
        'topo_or_crs must be a topology struct, EPSG code, projcrs or [].');
end
if isempty(tgt_epsg); tgt_epsg = epsg_of(tgt_crs); end

if ~isempty(opt.BBoxLatLon)
    bb = opt.BBoxLatLon;
    [x,y] = projfwd(tgt_crs,[bb(3) bb(3) bb(4) bb(4)],[bb(1) bb(2) bb(2) bb(1)]);
    bounds_xy = [min(x) max(x) min(y) max(y)] + opt.BufferM*[-1 1 -1 1];
elseif isstruct(topo_or_crs) && isfield(topo_or_crs,'cell_xc_projected_m')
    x = topo_or_crs.cell_xc_projected_m(:);
    y = topo_or_crs.cell_yc_projected_m(:);
    half = 0;
    if isfield(topo_or_crs,'cell_size_m') && ~isempty(topo_or_crs.cell_size_m)
        half = topo_or_crs.cell_size_m;
    end
    bounds_xy = [min(x) max(x) min(y) max(y)] + (opt.BufferM+half)*[-1 1 -1 1];
else
    error('Load_Field_DEM:NoExtent', ...
        'Pass a topology struct with projected cell centres or a BBoxLatLon.');
end
end

%% =======================================================================
function tf = same_crs(a,b)
try
    tf = strcmp(char(wktstring(a)),char(wktstring(b)));
catch
    tf = strcmp(char(a.Name),char(b.Name));
end
end

%% =======================================================================
function name = datum_name(crs)
name = '';
try
    name = char(crs.GeographicCRS.Datum.Name);
catch
end
end

%% =======================================================================
function e = epsg_of(crs)
e = NaN;
try
    m = regexp(char(wktstring(crs)),'ID\["EPSG",(\d+)\]','tokens');
    if ~isempty(m); e = str2double(m{end}{1}); end
catch
end
end

%% =======================================================================
function files = resolve_source_files(source)
if iscell(source)
    files = source(:).';
elseif ischar(source) || isstring(source)
    source = char(source);
    if isfolder(source)
        d = [dir(fullfile(source,'*.tif')); dir(fullfile(source,'*.tiff'))];
        files = arrayfun(@(f) fullfile(f.folder,f.name),d,'UniformOutput',false).';
    elseif isfile(source)
        files = {source};
    else
        files = {};
    end
else
    files = {};
end
files = files(cellfun(@isfile,files));
end

%% =======================================================================
function r = snap_resolution(res)
ladder = [0.5 1 2 3 5 10 30];
[~,i] = min(abs(ladder-res));
r = ladder(i);
end

%% =======================================================================
function Zf = fill_voids(Z,max_cells)
Zf = Z;
bad = ~isfinite(Zf);
if ~any(bad(:)) || all(bad(:)); return; end
L = label_components(bad);
for lab = 1:max(L(:))
    m = L == lab;
    if nnz(m) > max_cells; continue; end
    [ii,jj] = find(m);
    for k = 1:numel(ii)
        i = ii(k); j = jj(k);
        w = Zf(max(1,i-2):min(end,i+2),max(1,j-2):min(end,j+2));
        w = w(isfinite(w));
        if ~isempty(w); Zf(i,j) = mean(w); end
    end
end
end

%% =======================================================================
function L = label_components(BW)
L = zeros(size(BW)); next = 0;
for j = 1:size(BW,2)
    for i = 1:size(BW,1)
        if BW(i,j) && L(i,j) == 0
            next = next + 1; L(i,j) = next; stack = [i j];
            while ~isempty(stack)
                p = stack(end,:); stack(end,:) = [];
                for d = [-1 0;1 0;0 -1;0 1].'
                    ni = p(1)+d(1); nj = p(2)+d(2);
                    if ni>=1 && ni<=size(BW,1) && nj>=1 && nj<=size(BW,2) ...
                            && BW(ni,nj) && L(ni,nj) == 0
                        L(ni,nj) = next; stack(end+1,:) = [ni nj]; %#ok<AGROW>
                    end
                end
            end
        end
    end
end
end

%% =======================================================================
function out = ternary(cond,a,b)
if cond; out = a; else; out = b; end
end
