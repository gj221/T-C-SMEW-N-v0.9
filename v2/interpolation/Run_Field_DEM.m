function [topo_by_polygon, dem, figure_file] = Run_Field_DEM( ...
    shape_file, dem_source, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run_Field_DEM
%%
%% Ingest a field DEM, attach it to each calibrated polygon grid, and draw
%% the DEM / topographic-attribute preview figure -- elevation, slope,
%% aspect, TWI on the process grid, with the D-infinity drainage direction
%% and the topographic outlet cells marked.  Mirrors the
%% Run_*_Interpolation entry points.
%%
%% The kriged soil-property surfaces are a separate figure: run
%%   Run_Bruce_Interpolation   /  Run_AW_Bruce_Interpolation
%% (continuous surface, no cell mesh).
%%
%% INPUTS
%%   shape_file   field-boundary shapefile (e.g. the control fields)
%%   dem_source   a DEM GeoTIFF, a directory of tiles, or a {cell} of paths
%%                (default: v2/Spatial inputs/Bruce Farm Virginia/DEM)
%%
%% Name/value options
%%   'FeatureIndices'  polygon features to include (default [1 2] for the
%%                     Adam Wilbourne control fields, else all)
%%   'CellSizeM'       process-grid width [m] (default: the calibrated
%%                     coverage width per polygon)
%%   'FigureDirectory' where the PNG/FIG go (default v2/figures)
%%
%% OUTPUTS
%%   topo_by_polygon  cell array of topography-attached topology structs
%%   dem              the ingested gridded elevation model
%%   figure_file      path of the written PNG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

here = fileparts(mfilename('fullpath'));
v2   = fileparts(here);
addpath(here, fullfile(v2,'T&C_Code'));

opt = struct('FeatureIndices',[],'CellSizeM',[], ...
    'FigureDirectory',fullfile(v2,'figures'));
for k = 1:2:numel(varargin); opt.(varargin{k}) = varargin{k+1}; end

if nargin < 1 || isempty(shape_file)
    shape_file = 'AdamWilbourne_Control_Fields.shp';
    if isempty(opt.FeatureIndices); opt.FeatureIndices = [1 2]; end
end
shape_file = char(shape_file);
if ~isfile(shape_file)
    %% resolve a bare name against the tracked Bruce Farm location
    candidates = { ...
        fullfile(v2,'Spatial inputs','Bruce Farm Virginia','Whole farm',shape_file), ...
        fullfile(v2,'Spatial inputs','Bruce Farm Virginia',shape_file)};
    hit = candidates(cellfun(@isfile,candidates));
    if isempty(hit)
        error('Run_Field_DEM:MissingShapefile', ...
            '%s not found (also tried the Bruce Farm folders).', shape_file);
    end
    shape_file = hit{1};
end
if nargin < 2 || isempty(dem_source)
    dem_source = fullfile(v2,'Spatial inputs','Bruce Farm Virginia','DEM');
end
if ~exist(opt.FigureDirectory,'dir'); mkdir(opt.FigureDirectory); end

setups = Build_Shapefile_Grid_Setups(shape_file,opt.CellSizeM,opt.FeatureIndices);
n = numel(setups);
topo_by_polygon = cell(n,1);

fig = figure('Color','w','Position',[80 80 380*4 340*n]);
tl = tiledlayout(fig, n, 4, 'TileSpacing','compact','Padding','compact');
dem = [];

for k = 1:n
    topo = setups(k);
    d = Load_Field_DEM(dem_source, topo, ...
        'OutResolutionM',dem_resolution(opt.CellSizeM,topo));
    if isempty(dem); dem = d; end
    topo = Attach_Topography_To_Topology(topo, d);
    topo_by_polygon{k} = topo;
    g = topo.topography_grid;

    panels = { ...
        g.Z,           sprintf('elevation [m]  (relief %.1f m)', ...
                               max(topo.cell_z_m)-min(topo.cell_z_m)); ...
        100*g.slope,   sprintf('slope [%%]  (mean %.1f)',100*mean(topo.cell_slope)); ...
        g.aspect_deg,  'aspect [deg from N]'; ...
        g.twi,         sprintf('TWI  (%.1f - %.1f)',min(topo.cell_twi),max(topo.cell_twi))};
    for p = 1:4
        ax = nexttile(tl);
        imagesc(ax, g.x, g.y, panels{p,1}); set(ax,'YDir','normal');
        axis(ax,'image'); hold(ax,'on');
        if p == 4; colormap(ax,parula); else; colormap(ax,turbo); end
        colorbar(ax);
        %% field boundary
        plot(ax, topo.field_boundary_projected_x_m, ...
            topo.field_boundary_projected_y_m, 'w-', 'LineWidth',1.5);
        %% outlet cells
        oc = topo.cell_is_outlet;
        scatter(ax, topo.cell_xc_projected_m(oc), topo.cell_yc_projected_m(oc), ...
            12, 'w', 'filled', 'MarkerEdgeColor','k');
        if p == 1
            %% D-infinity drainage direction (subsampled quiver)
            step = max(1, round(numel(g.x)/25));
            [XX,YY] = meshgrid(g.x(1:step:end), g.y(1:step:end));
            U = g.flow_j(1:step:end,1:step:end);
            V = g.flow_i(1:step:end,1:step:end);
            quiver(ax, XX, YY, U, V, 0.6, 'k');
        end
        title(ax, panels{p,2}, 'FontWeight','normal');
        if p == 1
            ylabel(ax, sprintf('%s  (%d cells @ %.2f m)', ...
                Field_Setup_Label_safe(shape_file,setups(k),k), ...
                topo.N_cells, topo.cell_size_m), 'FontWeight','bold');
        end
    end
end
title(tl, sprintf('Field topography from %s', dem_source), ...
    'Interpreter','none');

[~,shp_name] = fileparts(shape_file);
figure_file = fullfile(opt.FigureDirectory, ...
    sprintf('Field_DEM_%s.png', shp_name));
exportgraphics(fig, figure_file, 'Resolution',150);
savefig(fig, strrep(figure_file,'.png','.fig'));
fprintf('Run_Field_DEM: wrote %s\n', figure_file);
end

%% =======================================================================
function s = Field_Setup_Label_safe(shape_file,setup,k)
try
    features = shaperead(shape_file);
    s = Field_Setup_Label(features, setup.source_feature_indices, k);
catch
    s = sprintf('polygon %d', k);
end
end

%% =======================================================================
function r = dem_resolution(requested,topo)
%% Default the DEM spacing to the process-grid cell size rather than the
%% tile's native 1 m.  Terrain attributes are sampled at cell centres, so
%% a DEM far finer than the cell aliases sub-cell microtopography into
%% slope / curvature / TWI (and slope, a gradient magnitude, is biased
%% one-sided by it).  Load_Field_DEM area-averages to whatever spacing it
%% is given.
r = requested;
if isempty(r) && isfield(topo,'cell_size_m') && ~isempty(topo.cell_size_m)
    r = topo.cell_size_m;
end
end
