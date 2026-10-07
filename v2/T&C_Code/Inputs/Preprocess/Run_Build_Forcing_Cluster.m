% CLUSTER side of the delta-forcing pipeline. Runs entirely offline.
%
% Inputs to stage here first:
%   forcing/deltas_<scenario>.mat      from Run_Extract_Deltas_Local
%   era5land/energyfarm_1995_2024.nc   from download_era5land_point.py
%   Data_Kantola.mat                   the validation forcing, needed when
%                                      cfg.Deaccumulation is 'legacy'
%
% Produces one T&C forcing .mat per scenario, each with Ca appended from the
% matching concentration pathway. On the legacy basis each is then checked
% against Data_Kantola.mat where they overlap.

addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))),'T&C_Code','Inputs','Preprocess'));

cfg = Delta_Forcing_Config();
era5_source = fullfile(pwd,'era5land','energyfarm_1995_2024.nc');
out_root = fullfile(pwd,'forcing');

if ~isfile(era5_source)
    error('Run_Build_Forcing_Cluster:MissingERA5', ...
        'ERA5-Land source not staged: %s',era5_source);
end

% The builder compares stored coordinates to cfg.Lat/Lon at 1e-8, so catch a
% mismatch here with a clear message rather than deep inside the pipeline.
file_lat = double(ncread(era5_source,'latitude'));
file_lon = double(ncread(era5_source,'longitude'));
if abs(file_lat-cfg.Lat)>1e-8 || abs(file_lon-cfg.Lon)>1e-8
    error('Run_Build_Forcing_Cluster:CoordinateMismatch', ...
        ['ERA5-Land file holds %.4f, %.4f but Delta_Forcing_Config has ' ...
        '%.4f, %.4f. Update the config to the cell centre in the file.'], ...
        file_lat,file_lon,cfg.Lat,cfg.Lon);
end

legacy_file = fullfile(pwd,'Data_Kantola.mat');
if strcmp(cfg.Deaccumulation,'legacy') && ~isfile(legacy_file)
    error('Run_Build_Forcing_Cluster:MissingLegacyForcing', ...
        'Validation forcing not staged: %s',legacy_file);
end

for s = 1:numel(cfg.scenarios)
    sc = cfg.scenarios{s};
    fprintf('\n=== %s ===\n',sc);

    delta_file = fullfile(out_root,sprintf('deltas_%s.mat',sc));
    if ~isfile(delta_file)
        error('Run_Build_Forcing_Cluster:MissingDeltas', ...
            'Deltas not staged: %s',delta_file);
    end

    perturbed_nc = fullfile(out_root,sprintf('era5land_%s_2016_2050.nc',sc));
    report = Apply_CMIP6_Deltas_To_ERA5_Land( ...
        era5_source,perturbed_nc,delta_file, ...
        cfg.Lat,cfg.Lon,cfg.Zbas,cfg.DeltaGMT, ...
        'YearMap',cfg.year_map, ...
        'PassthroughYears',cfg.observed, ...
        'PrecipMethod',cfg.PrecipMethod, ...
        'Deaccumulation',cfg.Deaccumulation, ...
        'SiteLat',cfg.SiteLat,'SiteLon',cfg.SiteLon);
    fprintf('  ramp mean %.3f\n',report.mean_ramp);

    forcing_mat = fullfile(out_root,sprintf('TC_forcing_%s_2016_2050.mat',sc));
    % The perturbed file now carries the site coordinates, so the builder
    % validates against those and writes them into the T&C forcing.
    Build_TC_Forcing_From_ERA5_Land_Timeseries( ...
        perturbed_nc,forcing_mat,cfg.SiteLat,cfg.SiteLon,cfg.Zbas, ...
        cfg.DeltaGMT);

    % Ca must come from the same pathway as the climate deltas.
    F = load(forcing_mat);
    Ca = Build_SSP_CO2_Series(F.Date,sc, ...
        'ConcentrationFile',cfg.ConcentrationFiles.(sc));
    save(forcing_mat,'Ca','-append');
    fprintf('  appended Ca to %s\n',forcing_mat);

    % On the legacy basis the observed years must reproduce the forcing the
    % spin-up and the calibration used.
    if strcmp(cfg.Deaccumulation,'legacy')
        Check_Forcing_Against_Legacy(forcing_mat,legacy_file);
    end
end

fprintf('\nDone. Three forcing sets written to %s\n',out_root);
