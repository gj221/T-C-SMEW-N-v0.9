function cfg = Delta_Forcing_Config()
% Shared settings for the delta-forcing pipeline.
%
% The pipeline runs across two machines: delta extraction needs internet and
% runs locally, while forcing assembly and the model runs happen on the
% cluster. Both sides read this file so the two cannot drift apart.

% ---- ERA5-Land cell ------------------------------------------------------
% These identify the grid cell the meteorology was downloaded from, and must
% equal the cell centre stored in that file; download_era5land_point.py prints
% the values. Apply_CMIP6_Deltas_To_ERA5_Land validates against them.
cfg.Lat      =  40.1;
cfg.Lon      = -88.2;

% ---- ERA5-Land de-accumulation ------------------------------------------
% Data_Kantola.mat, the forcing Data_Kantola_BG.mat was spun up on and the
% parameters were calibrated with, is this ERA5-Land record with one
% de-accumulation defect: across the 00 UTC reset, the hour ending 01 UTC
% was taken as that hour's accumulation less the previous day's total,
% clipped at zero. That drops about 4 % of the rain and 0.4 % of the
% shortwave. With no new spin-up the projections stay on that basis:
% 'legacy' reproduces it on the whole record, and Run_Build_Forcing_Cluster
% checks the observed years against Data_Kantola.mat; 'timeseries' keeps
% the product's hourly increments.
cfg.Deaccumulation = 'legacy';

% ---- model site ----------------------------------------------------------
% Where the model actually runs. These reach T&C through the forcing .mat and
% feed the solar geometry, so they must match the validated Data_Kantola.mat
% (40.062778, -88.196111, 221 m) or the trajectory runs are not comparable
% with the validated ones. They are deliberately NOT the cell centre.
cfg.SiteLat  =  40.062778;
cfg.SiteLon  = -88.196111;
cfg.Zbas     =  221;      % Energy Farm elevation, m a.s.l.
cfg.DeltaGMT =  -6;       % CST, standard time year-round

% ---- scenarios and ensemble --------------------------------------------
cfg.scenarios = {'ssp126','ssp245','ssp585'};
cfg.models = { ...
    'ACCESS-CM2','CanESM5','EC-Earth3','GFDL-ESM4','INM-CM5-0', ...
    'IPSL-CM6A-LR','MIROC6','MPI-ESM1-2-HR','MRI-ESM2-0','NorESM2-MM'};

cfg.BaselineYears = [1995 2014];   % matches the ERA5-Land resampling pool
cfg.FutureYears   = [2036 2065];   % centred on the ramp target
cfg.PrecipQuantiles = 0.05:0.05:0.95;
cfg.PrecipMethod = 'qdm';          % 'delta' for the plainer variant

% Quantile resolution per variable. Quantile delta mapping is the general
% method and a single mean factor is its degenerate case at one quantile, so
% these are recorded choices rather than omissions. Precipitation gets the
% full grid because heavy events intensify at roughly Clausius-Clapeyron
% while the mean rises far more slowly, so one factor misallocates the change.
% The others get 1 because at this site their projected signal is comparable
% to the inter-model spread -- shortwave especially, a 0.045 seasonal range
% against a 0.047 p10-p90 -- so quantile-resolved factors would fit sampling
% noise. Only pr is wired through; raising another entry needs matching
% changes in Extract_CMIP6_Monthly_Deltas and Apply_CMIP6_Deltas_To_ERA5_Land.
cfg.QuantileResolution = struct( ...
    'pr',      numel(cfg.PrecipQuantiles), ...
    'tas',     1, ...
    'hurs',    1, ...
    'rsds',    1, ...
    'sfcWind', 1);

% Official CMIP6 concentration pathways, one CSV per scenario with columns
% year,ppm, extracted from the Meinshausen et al. 2020 supplement data tables
% (GMD 13:3571, tables T4/T5/T11, World CO2 in ppm) and cross-checked against
% the input4MIPs files on greenhousegases.science.unimelb.edu.au. Resolved
% against the working directory, so keep them beside the pipeline .m files.
% Build_SSP_CO2_Series requires these; there is no approximate fallback.
cfg.ConcentrationFiles = struct( ...
    'ssp126','co2_ssp126_meinshausen2020.csv', ...
    'ssp245','co2_ssp245_meinshausen2020.csv', ...
    'ssp585','co2_ssp585_meinshausen2020.csv');

cfg.nexgddp_root = '/nfs/roberts/project/pi_np363/gj326/incoming/nexgddp_cell';

% ---- run period ---------------------------------------------------------
cfg.observed  = 2016:2024;   % real ERA5-Land, ramp forced to zero
cfg.projected = 2025:2050;   % resampled from the pool below
cfg.pool      = 1995:2014;

cfg.year_map = [cfg.observed(:) cfg.observed(:); ...
                leap_matched_map(cfg.projected,cfg.pool)];
end

function map = leap_matched_map(targets,pool)
% Pair each target year with a source year of the same leap status, cycling
% through the pool so the hourly axis stays continuous across the join.
is_leap = @(y) (mod(y,4)==0 & mod(y,100)~=0) | mod(y,400)==0;
leap_pool = pool(is_leap(pool));
norm_pool = pool(~is_leap(pool));
if isempty(leap_pool) || isempty(norm_pool)
    error('Delta_Forcing_Config:PoolTooNarrow', ...
        'The source pool must contain both leap and non-leap years.');
end
map = zeros(numel(targets),2);
li = 0; ni = 0;
for k = 1:numel(targets)
    t = targets(k);
    if is_leap(t)
        li = li + 1;
        src = leap_pool(mod(li-1,numel(leap_pool))+1);
    else
        ni = ni + 1;
        src = norm_pool(mod(ni-1,numel(norm_pool))+1);
    end
    map(k,:) = [t src];
end
end
