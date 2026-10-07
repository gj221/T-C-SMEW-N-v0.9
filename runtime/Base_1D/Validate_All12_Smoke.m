function Validate_All12_Smoke
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Three-day smoke of the calibrated set on one column and on twelve.
%%
%% Run 10 (one column) and run 1 (whole field, twelve columns with the
%% 5 m margin) are integrated for three days with the calibrated CEC 13
%% parameters, and the saved provenance of both is put through the same
%% assertion production applies, so an option that does not reach the
%% column configuration, or a workbook that no longer says CEC 13, fails
%% here in minutes rather than after a multi-day production run. The
%% twelve-column result is also checked for lateral conservation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

campaign_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(campaign_dir));
v2_dir = fullfile(repo_root,'v2');
two_d_dir = fullfile(repo_root,'runtime', ...
    'Base_2D');
addpath(fullfile(v2_dir,'run'),fullfile(v2_dir,'T&C_Code'), ...
    campaign_dir,two_d_dir);
P = Campaign_Parameters;
campaign_tag = strtrim(getenv('TC_CAMPAIGN_TAG'));
if isempty(campaign_tag)
    error('Bouchet_Kantola_CEC13:SmokeCampaignTag', ...
        'TC_CAMPAIGN_TAG must be set for the smoke validation.');
end

root = tempname(tempdir);
mkdir(root);
old = struct('ckpt',getenv('TC_CKPT_MINUTES'), ...
    'force',getenv('TC_CKPT_FORCE_RESUME'), ...
    'workers',getenv('TC_PARALLEL_WORKERS'), ...
    'stamp',getenv('TC_VALIDATION_FINGERPRINT'));
cleanup = onCleanup(@() restore_environment(root,old));
setenv('TC_CKPT_MINUTES','0');
setenv('TC_CKPT_FORCE_RESUME','');
setenv('TC_PARALLEL_WORKERS','0');
smoke_digest = 'cec13_all12_smoke';
setenv('TC_VALIDATION_FINGERPRINT',smoke_digest);
n_days = 3;

%% One column: run 10 (control, vertical_1d).
dir_1d = fullfile(root,'run10');
mkdir(dir_1d);
file_1d = Run_One_Subfield(10, ...
    'StartDate',P.StartDate,'NDays',n_days, ...
    'AppeloAlpha',P.AppeloAlpha,'AppeloLogKH',P.AppeloLogKH, ...
    'KNmaxd',P.KNmaxd,'Kbot',P.Kbot,'ParallelWorkers',0, ...
    'SaveResults',true,'ResultDirectory',dir_1d);
saved = load(file_1d,'run_provenance');
Assert_Base_Provenance(saved.run_provenance,10,1,0, ...
    smoke_digest,campaign_tag,dir_1d,n_days);
fprintf('Smoke: run 10 (one column, %d days) carries the calibrated identity.\n', ...
    n_days);

%% Twelve columns: run 1 (control, whole_field, 2 x 2 main grid + margin).
dir_2d = fullfile(root,'run01');
mkdir(dir_2d);
file_2d = Run_One_Subfield(1, ...
    'StartDate',P.StartDate,'NDays',n_days, ...
    'GridNx',2,'GridNy',2, ...
    'AppeloAlpha',P.AppeloAlpha,'AppeloLogKH',P.AppeloLogKH, ...
    'KNmaxd',P.KNmaxd,'Kbot',P.Kbot,'ParallelWorkers',0, ...
    'SaveResults',true,'ResultDirectory',dir_2d);
saved = load(file_2d,'run_provenance');
Assert_Base_Provenance(saved.run_provenance,1,12,0, ...
    smoke_digest,campaign_tag,dir_2d,n_days);
fprintf('Smoke: run 1 (twelve columns, %d days) carries the calibrated identity.\n', ...
    n_days);
Check_Lateral_Conservation(file_2d,'Bouchet_Kantola_CEC13');
fprintf('PASS: calibrated CEC 13 smoke on one and twelve columns\n');
end

function restore_environment(root,old)
setenv('TC_CKPT_MINUTES',old.ckpt);
setenv('TC_CKPT_FORCE_RESUME',old.force);
setenv('TC_PARALLEL_WORKERS',old.workers);
setenv('TC_VALIDATION_FINGERPRINT',old.stamp);
if isfolder(root) && startsWith(root,tempdir)
    rmdir(root,'s');
end
end
