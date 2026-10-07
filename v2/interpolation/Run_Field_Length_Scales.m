function out = Run_Field_Length_Scales(shape_file, dem_source, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run_Field_Length_Scales
%%
%% Compute the two measured length scales that set the topographic 2-D
%% process-grid cell size for a field, and the recommended width:
%%
%%   L_soil   shortest fitted variogram range across the measured soil
%%            attributes (bulk pH, sand, clay, OM, CEC) at the field --
%%            the tightest heterogeneity the grid must resolve.
%%   L_terr   the lag at which the detrended-elevation semivariogram of the
%%            DEM over the field reaches ~95 % sill -- the microtopography
%%            autocorrelation length.
%%   h*       target cell width = max( DEM native spacing,
%%                                     min(L_soil, L_terr) / k ),  k = 3
%%            (>= 3 cells per correlation length), capped by the 95 %
%%            per-cell coverage rule and by any compute cap.
%%
%% This replaces the fixed 80 m "notebook recipe" variogram range for the
%% purpose of grid sizing -- that default is not a Bruce Farm measurement.
%%
%% INPUTS
%%   shape_file   the field-boundary shapefile (e.g. the control fields)
%%   dem_source   a DEM GeoTIFF, a directory of tiles, or a {cell} of paths
%%
%% Name/value options
%%   'FeatureIndices'  polygon features to include (default all)
%%   'SampleFile'      soil-sample workbook (default: the Bruce Farm
%%                     agronomic workbook resolved by Bruce_Farm_Input_File)
%%   'SampleSheet'     workbook sheet (default 2, the whole-grower sheet
%%                     that reaches the control polygons)
%%   'Attributes'      {cellstr} of sample columns to fit (default the
%%                     bulk pH / texture / OM / CEC set, auto-detected)
%%   'NyquistK'        cells per correlation length (default 3)
%%   'ComputeCapCells' largest total cell count the ERW run can afford
%%                     ([] = no cap)
%%   'Seed'            seed of the DEM-point subsample used for the
%%                     terrain variogram (default 0, reproducible)
%%
%% OUTPUT struct with .L_soil_m .L_terr_m .dem_native_m .h_target_m
%%   .h_coverage_min_m .h_recommended_m .per_attribute (table) and the
%%   fitted variograms, plus a printed summary.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

here = fileparts(mfilename('fullpath'));
v2   = fileparts(here);
addpath(here, fullfile(v2,'T&C_Code'));

opt = struct('FeatureIndices',[],'SampleFile','','SampleSheet',2, ...
    'Attributes',{{}},'NyquistK',3,'ComputeCapCells',[],'Seed',0);
for k = 1:2:numel(varargin); opt.(varargin{k}) = varargin{k+1}; end

%% --- field setups + calibrated coverage cell size ------------------
setups = Build_Shapefile_Grid_Setups(shape_file,[],opt.FeatureIndices);
h_coverage = min(arrayfun(@(s) s.cell_size_m, setups));   % smallest that passes 95%

%% --- L_soil: fit a variogram to each measured attribute -----------
sample_file = opt.SampleFile;
if isempty(sample_file)
    sample_file = Bruce_Farm_Input_File(fullfile('Whole farm','Data', ...
        'RJ_msft_yale_geochem - agronomic.xlsx'));
end
points = Read_Sample_Points(sample_file,opt.SampleSheet,'lng (X)','lat (Y)');
attrs = opt.Attributes;
if isempty(attrs)
    %% The same names Krige_Cell_Agronomic_Values' default map knows, so
    %% the two auto-detectors in this module agree on what the site's
    %% attributes are called (the workbook says 'Cation Exchange
    %% Capacity' and 'Organic Matter'; the catalogue export 'CEC Value',
    %% 'Organic Matter (%)', 'Sand Mass (%)', 'Clay Mass (%)').
    candidates = {'pH','pH Value','pH_BL','Bulk pH', ...
        'sand','clay','sand_pct','clay_pct','Sand Mass (%)','Clay Mass (%)', ...
        'Percent Sand','Percent Clay', ...
        'OM','OM_pct','organic matter','Organic Matter (%)', ...
        'CEC','CEC_meq','CEC Value','Cation Exchange Capacity'};
    attrs = candidates(ismember(lower(candidates), lower(points.Properties.VariableNames)));
    %% case-insensitive resolve back to the real column names
    for a = 1:numel(attrs)
        hit = find(strcmpi(points.Properties.VariableNames,attrs{a}),1);
        attrs{a} = points.Properties.VariableNames{hit};
    end
end
if isempty(attrs)
    error('Run_Field_Length_Scales:NoAttributes', ...
        ['Could not auto-detect soil attribute columns in %s. Pass ' ...
         '''Attributes'' explicitly.'], sample_file);
end

%% One variogram PER POLYGON and attribute.  Project_Sample_Points returns
%% x_grid_m/y_grid_m in each setup's OWN local frame (its own origin and
%% grid rotation), so coordinates from different polygons cannot be
%% pooled into one cloud -- the pairwise distances across polygons would
%% be meaningless.  L_soil is the minimum over the resolved fits.
rows = cell(0,5);
L_soil = inf; soil_variograms = struct();
for a = 1:numel(attrs)
    col = attrs{a};
    for k = 1:numel(setups)
        s = Project_Sample_Points(points,setups(k),'lng (X)','lat (Y)');
        v = Numeric_Sample_Column(s,col);
        good = isfinite(v);
        xk = s.x_grid_m(good); yk = s.y_grid_m(good); vk = v(good);
        if numel(vk) < 10
            rows(end+1,:) = {col, k, numel(vk), NaN, 'too few samples'}; %#ok<AGROW>
            continue
        end
        try
            vg = Fit_Empirical_Variogram_Range(xk,yk,vk);
            soil_variograms.(sprintf('%s_polygon_%d', ...
                matlab.lang.makeValidName(col),k)) = vg;
            if vg.resolved
                note = sprintf('rmse %.3g',vg.rmse);
                L_soil = min(L_soil, vg.range_m);
            else
                note = sprintf('unresolved (nugget %.0f%%, lag>=%.0f m)', ...
                    100*vg.nugget_fraction, vg.shortest_lag_m);
            end
            rows(end+1,:) = {col, k, vg.n_samples, vg.range_m, note}; %#ok<AGROW>
        catch e
            rows(end+1,:) = {col, k, numel(vk), NaN, e.message}; %#ok<AGROW>
        end
    end
end
per_attribute = cell2table(rows, ...
    'VariableNames',{'attribute','polygon','n_samples','range_m','note'});
if ~isfinite(L_soil)
    warning('Run_Field_Length_Scales:NoSoilRange', ...
        'No attribute yielded a variogram range; L_soil left as NaN.');
    L_soil = NaN;
end

%% --- L_terr: detrended-elevation variogram of the DEM ------------
dem = Load_Field_DEM(dem_source,setups(1));
[XX,YY] = meshgrid(dem.x,dem.y);
z = dem.Z; ok = isfinite(z);
xr = XX(ok); yr = YY(ok); zr = z(ok);
G = [xr-mean(xr), yr-mean(yr), ones(nnz(ok),1)];
resid = zr - G*(G\zr);                           % remove the regional plane
%% Reproducible subsample for tractability (private stream: the caller's
%% global generator is left alone).
n = numel(zr);
idx = randperm(RandStream('mt19937ar','Seed',opt.Seed), n, min(n,4000));
vg_terr = Fit_Empirical_Variogram_Range(xr(idx), yr(idx), resid(idx), ...
    'NBins',20,'MinPairs',50);
L_terr = vg_terr.range_m;
dem_native = dem.native_resolution_m;

%% --- recommendation -------------------------------------------
%% A length scale only constrains the grid when it is both resolved (the
%% sampling can see it) and shorter than the field itself -- a range
%% longer than the field means "no internal structure here", which
%% relaxes the grid, it does not coarsen it.  h_resolve is the finest
%% width any measured heterogeneity actually demands; Inf if neither
%% L_soil nor L_terr binds, in which case the grid is set by the
%% existing 95%-coverage calibration alone.
field_extent = sqrt(sum(arrayfun(@(s) s.field_area_m2, setups)));
h_resolve = inf;
if isfinite(L_soil) && L_soil < field_extent
    h_resolve = min(h_resolve, L_soil/opt.NyquistK);
end
if isfinite(L_terr) && L_terr < field_extent
    h_resolve = min(h_resolve, L_terr/opt.NyquistK);
end
h_target = max(dem_native, h_coverage);
if isfinite(h_resolve)
    h_target = max(dem_native, min(h_target, h_resolve));
end
h_rec = h_target;
compute_note = '';
if ~isempty(opt.ComputeCapCells)
    A_field = sum(arrayfun(@(s) s.field_area_m2, setups));
    h_cpu = sqrt(A_field/opt.ComputeCapCells);
    if h_cpu > h_rec
        compute_note = sprintf(' compute cap forces %.2f m', h_cpu);
        h_rec = h_cpu;
    end
end

out = struct('L_soil_m',L_soil,'L_terr_m',L_terr,'dem_native_m',dem_native, ...
    'h_target_m',h_target,'h_coverage_min_m',h_coverage, ...
    'h_recommended_m',h_rec,'nyquist_k',opt.NyquistK, ...
    'per_attribute',per_attribute,'soil_variograms',soil_variograms, ...
    'terrain_variogram',vg_terr);

fprintf('\n=== Field length scales: %s ===\n', shape_file);
disp(per_attribute);
fprintf('  field extent (sqrt area)            %6.1f m\n', field_extent);
fprintf('  L_soil (shortest RESOLVED range)    %s\n', range_note(L_soil,field_extent));
fprintf('  L_terr (DEM microtopo range)        %s\n', range_note(L_terr,field_extent));
fprintf('  DEM native spacing                  %6.2f m\n', dem_native);
fprintf('  95%% per-cell coverage calibration    %6.2f m\n', h_coverage);
if isfinite(h_resolve)
    fprintf(['  measured heterogeneity demands  <= %6.2f m (min resolved ' ...
        'range / %d)\n'], h_resolve, opt.NyquistK);
else
    fprintf(['  measured heterogeneity demands     none (no attribute or ' ...
        'terrain range is resolved and shorter than the field)\n']);
end
fprintf('  RECOMMENDED process-grid width      %6.2f m%s\n', h_rec, compute_note);
end

%% =======================================================================
function s = range_note(L,field_extent)
if ~isfinite(L)
    s = '   n/a (unresolved)';
elseif L >= field_extent
    s = sprintf('%6.1f m  (>= field extent -> no internal structure)',L);
else
    s = sprintf('%6.1f m',L);
end
end
