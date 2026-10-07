function out = Run_Topographic_Grid_Convergence(shape_file, dem_source, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run_Topographic_Grid_Convergence
%%
%% Confirm a process-grid width by convergence rather than by assertion.
%% Over a ladder of widths the DEM is area-averaged to each width
%% (Load_Field_DEM) and the terrain attributes the 2-D flow code depends
%% on are derived at that width (Field_Topographic_Attributes) on the
%% RASTER, masked to the polygon.  That is deliberately grid-free: the
%% question "at what resolution do slope, TWI and curvature stop
%% changing" is a property of the terrain, and the polygon-grid builder's
%% 95 % coverage rule (Build_Projected_Field_Topology) would otherwise
%% forbid every width coarser than the calibrated one.  The graph
%% quantities that DO need a topology -- cell count, outlet count and the
%% share of the polygon that drains to an outlet through the MFD graph --
%% are reported at the widths the builder accepts and NaN elsewhere.
%%
%% Statistics (polygon-masked, DEM pixels weighted equally):
%%   slope       mean and spatial std   -- the Darcy / MFD driver
%%   TWI         mean and spatial std   -- where water converges
%%   |plan curv| mean                   -- tabulated only (see below)
%%   relief      max(z) - min(z)        -- the total head available
%% Slope and curvature are gradient magnitudes, so a coarser grid smooths
%% them low; that bias shrinking as the width decreases IS the
%% convergence being measured.  A statistic has converged at width h when
%% its change to the next finer width, AND every change between finer
%% widths in the ladder, is below 'Tolerance' (monotone criterion): one
%% quiet step between two noisy ones is not convergence.  The verdict is
%% reported per statistic and jointly; the finest width cannot be
%% assessed.
%%
%% Optionally ('RedistributionDays' > 0) a short topographic 2-D run is
%% made at each buildable width through Run_One_Subfield and two
%% redistribution metrics are added: lateral water moved per m2 of the
%% simulated grid per day, and the spatial coefficient of variation of
%% end-of-run water over the whole soil profile.  Expensive; off by default.
%%
%% INPUTS
%%   shape_file   field-boundary shapefile
%%   dem_source   DEM GeoTIFF / tile directory / {cell} of paths
%% Name/value
%%   'FeatureIndices'     polygon features (default all)
%%   'PolygonIndex'       which setup to analyse (default 1)
%%   'Widths'             ladder [m], coarse -> fine (default
%%                        [20 10 5 3.5 2.5 1.5 1]); widths below the DEM's
%%                        native spacing are dropped
%%   'Tolerance'          relative change that counts as converged (0.05)
%%   'CoverageRule'       the grid rule of the graph metrics and the
%%                        redistribution runs (Build_Shapefile_Grid_Setups;
%%                        default the 95 % rule, 'majority' builds every
%%                        width)
%%   'RedistributionDays' 0 (default) = attributes only
%%   'RunId'              Run_Matrix_12 id for the redistribution runs (1)
%%   'OutputDirectory'    default v2/figures
%%
%% OUTPUT struct: .table (one row per width), .converged_width_m (the
%%   coarsest width from which every verdict statistic stays converged),
%%   .converged_width_by_stat (table, one column per statistic; NaN where
%%   the ladder never converges it), .csv_file, .figure_file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

here = fileparts(mfilename('fullpath'));
v2   = fileparts(here);
addpath(here, fullfile(v2,'T&C_Code'), fullfile(v2,'run'));

opt = struct('FeatureIndices',[],'PolygonIndex',1,'Widths',[20 10 5 3.5 2.5 1.5 1], ...
    'Tolerance',0.05,'RedistributionDays',0,'RunId',1, ...
    'OutputDirectory',fullfile(v2,'figures'),'CoverageRule',[]);
for k = 1:2:numel(varargin); opt.(varargin{k}) = varargin{k+1}; end
if ~exist(opt.OutputDirectory,'dir'); mkdir(opt.OutputDirectory); end
widths = sort(opt.Widths(:).','descend');

%% The calibrated setup supplies the polygon (projected) and the DEM's
%% native spacing, which bounds the ladder from below.
probe = Build_Shapefile_Grid_Setups(shape_file,[],opt.FeatureIndices);
ref = probe(opt.PolygonIndex);
dem0 = Load_Field_DEM(dem_source,ref);
native = dem0.native_resolution_m;
drop = widths < native - 1e-9;
if any(drop)
    fprintf(['Run_Topographic_Grid_Convergence: dropping width(s) %s m ' ...
        'below the DEM native spacing %.2f m\n'], mat2str(widths(drop)), native);
    widths = widths(~drop);
end
poly_x = ref.field_boundary_projected_x_m(:);
poly_y = ref.field_boundary_projected_y_m(:);

names = {'width_m','slope_mean','slope_std','twi_mean','twi_std', ...
    'plan_curv_abs_mean','relief_m','n_pixels','n_cells','n_outlets', ...
    'drained_area_fraction','lateral_water_m3_per_m2_day','profile_water_cv','wall_s'};
rows = nan(numel(widths),numel(names));

for i = 1:numel(widths)
    h = widths(i);
    t0 = tic;
    %% --- raster attributes at this width, masked to the polygon ------
    dem = Load_Field_DEM(dem_source,ref,'OutResolutionM',h);
    g = Field_Topographic_Attributes(dem);
    [X,Y] = meshgrid(g.x,g.y);
    in = inpolygon(X,Y,poly_x,poly_y) & isfinite(g.Z);
    if nnz(in) == 0
        rows(i,1) = h; rows(i,8) = 0; rows(i,14) = toc(t0);
        fprintf('  h = %5.2f m: no DEM pixel centre inside the polygon; row left NaN\n',h);
        continue
    end
    slope = g.slope(in); twi = g.twi(in); curv = abs(g.curv_plan(in)); z = g.Z(in);
    twi = twi(isfinite(twi));
    rows(i,1) = h;
    rows(i,2) = mean(slope);        rows(i,3) = std(slope);
    rows(i,4) = mean(twi);          rows(i,5) = std(twi);
    rows(i,6) = mean(curv);
    rows(i,7) = max(z) - min(z);
    rows(i,8) = nnz(in);
    %% --- graph metrics only where the polygon grid can be built -------
    try
        setups = Build_Shapefile_Grid_Setups(shape_file,h,opt.FeatureIndices, ...
            opt.CoverageRule);
        topo = Attach_Topography_To_Topology(setups(opt.PolygonIndex),dem);
        [~,~,ff] = Field_Aggregation_Areas(topo);
        w = ff(:).*topo.cell_area_m2(:); w = w/sum(w);
        rows(i,9)  = topo.N_cells;
        rows(i,10) = nnz(topo.cell_is_outlet);
        rows(i,11) = drained_fraction(topo,w);
        if opt.RedistributionDays > 0
            [rows(i,12),rows(i,13)] = redistribution_metrics(shape_file,dem_source, ...
                opt,h,topo);
        end
        grid_note = sprintf('%d cells, %d outlet(s), %.0f%% drains', ...
            rows(i,9),rows(i,10),100*rows(i,11));
    catch ME
        if strcmp(ME.identifier,'Build_Projected_Field_Topology:CoverageBelowTarget')
            grid_note = 'no grid (coverage rule)';
        else
            rethrow(ME);
        end
    end
    rows(i,14) = toc(t0);
    fprintf(['  h = %5.2f m: %6d px, slope %.4f +/- %.4f, TWI %.2f +/- %.2f, ' ...
        '|curv| %.4f, relief %.2f m; %s; %.0f s\n'], ...
        h, rows(i,8), rows(i,2), rows(i,3), rows(i,4), rows(i,5), rows(i,6), ...
        rows(i,7), grid_note, rows(i,14));
end
T = array2table(rows,'VariableNames',names);

%% Convergence: relative change of each attribute statistic from row i
%% (coarser) to row i+1 (finer).  A width is converged when every
%% statistic changes by less than the tolerance from there to the next
%% finer width; the finest width cannot be assessed.
%% Plan curvature is tabulated but deliberately NOT part of the verdict:
%% a second derivative of a 1 m lidar surface is dominated by sub-metre
%% roughness (furrows, noise) and grows without limit as the width
%% shrinks -- on Control 2 it went 0.012 -> 0.083 1/m from 20 m to 1.5 m
%% while slope, relief and TWI had settled -- so it has no convergence
%% width to find.  The flow physics is driven by slope and TWI.
stat_cols = {'slope_mean','slope_std','twi_mean','twi_std','relief_m'};
rel = nan(height(T),numel(stat_cols));
for i = 1:height(T)-1
    for s = 1:numel(stat_cols)
        a = T.(stat_cols{s})(i); b = T.(stat_cols{s})(i+1);
        rel(i,s) = abs(a-b)/max(abs(b),eps);
    end
end
T.max_rel_change_to_finer = max(rel,[],2);
n = height(T);
within = rel < opt.Tolerance;                 % last row is NaN -> false
from_here = false(n,numel(stat_cols));        % converged from row i down
for s = 1:numel(stat_cols)
    for i = 1:n-1
        from_here(i,s) = all(within(i:n-1,s));
    end
end
T.converged = all(from_here,2);
converged_width = NaN;
idx = find(T.converged,1,'first');
if ~isempty(idx); converged_width = T.width_m(idx); end
by_stat = nan(1,numel(stat_cols));
for s = 1:numel(stat_cols)
    j = find(from_here(:,s),1,'first');
    if ~isempty(j); by_stat(s) = T.width_m(j); end
end
converged_by_stat = array2table(by_stat,'VariableNames',stat_cols);

[~,shp] = fileparts(shape_file);
csv_file = fullfile(opt.OutputDirectory, ...
    sprintf('Grid_Convergence_%s_polygon_%02d.csv',shp,opt.PolygonIndex));
writetable(T,csv_file);

fig = figure('Color','w','Visible','off','Position',[80 80 1100 640]);
try; fig.Theme = 'light'; catch; end       % batch sessions default to dark axes
tl = tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');
panels = {'slope_mean','slope [-]'; 'slope_std','slope spatial std [-]'; ...
    'twi_mean','TWI [-]'; 'twi_std','TWI spatial std [-]'; ...
    'plan_curv_abs_mean','|plan curvature| [1/m]'; 'relief_m','relief [m]'};
for p = 1:size(panels,1)
    ax = nexttile(tl);
    semilogx(ax,T.width_m,T.(panels{p,1}),'o-','LineWidth',1.3); hold(ax,'on');
    %% dashed line: this statistic's own converged width (none for the
    %% tabulated-only curvature and for a statistic that never converges)
    k = find(strcmp(stat_cols,panels{p,1}),1);
    if ~isempty(k) && isfinite(by_stat(k))
        xline(ax,by_stat(k),'--','Color',[0.4 0.4 0.4]);
    end
    set(ax,'XDir','reverse'); grid(ax,'on');
    xlabel(ax,'cell width [m]  (finer \rightarrow)'); ylabel(ax,panels{p,2});
end
title(tl,sprintf(['%s polygon %d: terrain-attribute convergence ' ...
    '(all statistics within %.0f%%: %s; dashed = per-statistic width)'], ...
    shp,opt.PolygonIndex,100*opt.Tolerance,ternary(isfinite(converged_width), ...
    sprintf('from %.2f m',converged_width),'none in ladder')), ...
    'Interpreter','none','FontSize',10);
figure_file = strrep(csv_file,'.csv','.png');
exportgraphics(fig,figure_file,'Resolution',150); close(fig);

out = struct('table',T,'converged_width_m',converged_width, ...
    'converged_width_by_stat',converged_by_stat, ...
    'tolerance',opt.Tolerance,'csv_file',csv_file,'figure_file',figure_file);
fprintf('Run_Topographic_Grid_Convergence: %s\n', ternary(isfinite(converged_width), ...
    sprintf('every verdict statistic stays within %.0f%% of the next finer width from %.2f m down', ...
    100*opt.Tolerance,converged_width), ...
    sprintf('no width in %s m keeps every statistic within %.0f%% through the ladder', ...
    mat2str(widths),100*opt.Tolerance)));
fprintf('  per statistic, coarsest width from which every finer step stays within %.0f%%:\n', ...
    100*opt.Tolerance);
for s = 1:numel(stat_cols)
    fprintf('    %-12s %s\n',stat_cols{s},ternary(isfinite(by_stat(s)), ...
        sprintf('%.2f m',by_stat(s)),'not converged within the ladder'));
end
disp(T);
end

%% =======================================================================
function frac = drained_fraction(topo,w)
%% Share of polygon area whose MFD flow path reaches an outlet: walk the
%% weight digraph downslope from every cell; a cell drains if any path
%% ends at an outlet (a pit does not).
N = topo.N_cells;
src = [topo.edge_a(topo.edge_overland_weight_ab > 0); ...
       topo.edge_b(topo.edge_overland_weight_ba > 0)];
dst = [topo.edge_b(topo.edge_overland_weight_ab > 0); ...
       topo.edge_a(topo.edge_overland_weight_ba > 0)];
reaches = logical(topo.cell_is_outlet(:));
if ~isempty(src)
    G = digraph(src(:),dst(:),ones(numel(src),1),N);
    D = distances(G);                              % N x N hop counts
    outlets = find(reaches);
    reaches = reaches | any(isfinite(D(:,outlets)),2);
end
frac = sum(w(reaches));
end

%% =======================================================================
function [lat_per_area_day,profile_cv] = redistribution_metrics(shape_file,dem_source,opt,h,topo)
%% A short topographic 2-D run at this width; the caller's Run_Matrix
%% id selects the parameter case.  Reads the lateral ledger the driver
%% keeps (Domain.lateral) and the end-of-run profile water per cell.
result_file = Run_One_Subfield(opt.RunId,'FieldShapefile',shape_file, ...
    'FieldFeatureIndices',opt.FeatureIndices,'FieldPolygonIndex',opt.PolygonIndex, ...
    'GridCellSizeM',h,'FieldDEM',dem_source,'DEMResolutionM',h, ...
    'FieldMinCoverage',opt.CoverageRule, ...
    'NDays',opt.RedistributionDays,'SaveResults',true);
S = load(result_file,'Domain','V');
%% net subsurface lateral redistribution: the water that changed cells,
%% counted once (the |Q_edge| throughput would also count exchange back
%% and forth; the buffer/main ledger is zero on a polygon grid)
lat = sum(S.Domain.lateral.water_net_m3,'omitnan');               % m3 moved
%% the throughput is over every simulated cell, so it is per m2 of the
%% simulated grid, not of the field inside it
lat_per_area_day = lat/max(S.Domain.simulated_area_m2,eps)/opt.RedistributionDays;
Vend = squeeze(sum(S.V(end,:,:),2));       % mm per cell, all layers (whole profile)
[~,~,ff] = Field_Aggregation_Areas(topo);
w = ff(:).*topo.cell_area_m2(:); w = w/sum(w);
mu = sum(w.*Vend(:)); sd = sqrt(max(0,sum(w.*Vend(:).^2)-mu^2));
profile_cv = sd/max(mu,eps);
end

%% =======================================================================
function o = ternary(c,a,b); if c; o = a; else; o = b; end; end
