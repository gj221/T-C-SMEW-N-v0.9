% LOCAL side of the delta-forcing pipeline. Needs internet; run on a machine
% with outbound access to NCCS, not on a cluster compute node.
%
% Produces three small .mat files (a few kB each) to transfer to the cluster,
% where Run_Build_Forcing_Cluster consumes them offline.
%
% Prerequisite: ERA5-Land downloaded separately with
% Inputs/Preprocess/download_era5land_point.py, which also needs internet.

addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))),'T&C_Code','Inputs','Preprocess'));

cfg = Delta_Forcing_Config();
out_root = fullfile(pwd,'forcing');
if ~isfolder(out_root); mkdir(out_root); end

fprintf('Extracting deltas for %d scenarios, %d models.\n', ...
    numel(cfg.scenarios),numel(cfg.models));
fprintf('This is network-bound: expect a long run and restartable failures.\n\n');

for s = 1:numel(cfg.scenarios)
    sc = cfg.scenarios{s};
    delta_file = fullfile(out_root,sprintf('deltas_%s.mat',sc));
    if isfile(delta_file)
        fprintf('=== %s: already extracted, skipping ===\n',sc);
        continue
    end
    fprintf('=== %s ===\n',sc);
    Extract_CMIP6_Monthly_Deltas(cfg.nexgddp_root,delta_file, ...
        cfg.Lat,cfg.Lon,sc,cfg.models, ...
        'BaselineYears',cfg.BaselineYears, ...
        'FutureYears',cfg.FutureYears, ...
        'PrecipQuantiles',cfg.PrecipQuantiles);
end

fprintf('\nDone. Transfer these to the cluster alongside the ERA5-Land NetCDF:\n');
for s = 1:numel(cfg.scenarios)
    f = fullfile(out_root,sprintf('deltas_%s.mat',cfg.scenarios{s}));
    d = dir(f);
    if ~isempty(d)
        fprintf('  %s  (%.1f kB)\n',f,d.bytes/1024);
    end
end
