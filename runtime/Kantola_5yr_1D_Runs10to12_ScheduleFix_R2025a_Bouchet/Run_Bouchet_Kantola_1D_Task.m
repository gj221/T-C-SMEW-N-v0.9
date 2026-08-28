function Run_Bouchet_Kantola_1D_Task
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run one corrected five-year Kantola 1-D treatment on Bouchet.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

campaign_tag = 'Kantola_5yr_1D_Runs10to12_ScheduleFix_R2025a_Bouchet';
campaign_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(campaign_dir));
run_dir = fullfile(repo_root,'v2','run');
code_dir = fullfile(repo_root,'v2','T&C_Code');
addpath(run_dir,code_dir,fullfile(code_dir,'Inputs'), ...
    fullfile(code_dir,'Parameters'));

assert_single_thread_runtime;
run_id = read_positive_integer('TC_RUN_ID');
if run_id < 10 || run_id > 12
    error('Bouchet_Kantola_1D:BadRunId', ...
        'TC_RUN_ID must select run 10, 11 or 12.');
end
if read_positive_integer('SLURM_CPUS_PER_TASK') ~= 1
    error('Bouchet_Kantola_1D:BadCpuCount', ...
        'Each production treatment must reserve exactly one CPU.');
end
if ~strcmp(strtrim(getenv('TC_PARALLEL_WORKERS')),'0')
    error('Bouchet_Kantola_1D:BadWorkerCount', ...
        'TC_PARALLEL_WORKERS must be 0 for a one-column run.');
end
if ~strcmp(strtrim(getenv('TC_CKPT_MINUTES')),'30')
    error('Bouchet_Kantola_1D:BadCheckpointInterval', ...
        'TC_CKPT_MINUTES must be 30.');
end
if ~isempty(strtrim(getenv('TC_CKPT_FORCE_RESUME')))
    error('Bouchet_Kantola_1D:ForcedResumeForbidden', ...
        'TC_CKPT_FORCE_RESUME must remain unset for production.');
end

required_files = { ...
    fullfile(code_dir,'Inputs','ERW_MOD_PARAM_parameters.xlsx'), ...
    fullfile(code_dir,'Inputs','Data_Kantola.mat'), ...
    fullfile(code_dir,'Inputs','Data_Kantola_BG.mat'), ...
    fullfile(code_dir,'Inputs','Ca_Data.mat'), ...
    fullfile(code_dir,'Inputs','All_deposition_data.mat'), ...
    fullfile(code_dir,'Inputs','blueridge_CDF_data.json')};
for k = 1:numel(required_files)
    if ~isfile(required_files{k})
        error('Bouchet_Kantola_1D:MissingDeploymentFile', ...
            'Missing required deployment file: %s',required_files{k});
    end
end

results_dir = fullfile(repo_root,'v2','SIM','Subfield_12',campaign_tag);
logs_dir = fullfile(repo_root,'v2','logs','Subfield_12',campaign_tag);
stamp_file = fullfile(logs_dir,'validation_stamp.txt');
if ~isfile(stamp_file)
    error('Bouchet_Kantola_1D:MissingValidation', ...
        'Run the R2025a validation job before production.');
end
expected_stamp = Bouchet_Campaign_Validation_Fingerprint( ...
    fullfile(repo_root,'v2'),campaign_dir,campaign_tag,1,0);
if ~strcmp(strtrim(fileread(stamp_file)),expected_stamp)
    error('Bouchet_Kantola_1D:StaleValidation', ...
        'Code, inputs, tests or numerical runtime changed after validation.');
end
%% The stamp this process verified travels into the saved result so a later
%% job can tell which validation produced each MAT.  It travels as a digest
%% because the full fingerprint is far past setenv's 32766-character limit.
stamp_digest = Bouchet_Fingerprint_Digest(expected_stamp);
setenv('TC_VALIDATION_FINGERPRINT',stamp_digest);
if ~exist(results_dir,'dir'); mkdir(results_dir); end
if ~exist(logs_dir,'dir'); mkdir(logs_dir); end

%% A missing final MAT only means that no run has finished.  Take the lock
%% before anything else touches the status file, checkpoint or result path,
%% so a duplicate submission fails here rather than racing.
lock_cleanup = Bouchet_Acquire_Run_Lock( ...
    fullfile(logs_dir,'locks',sprintf('run_%02d.lock',run_id))); %#ok<NASGU>

if strcmp(Bouchet_Check_Existing_Result( ...
        results_dir,run_id,stamp_digest,1,1826),'complete')
    fprintf('Run %d needs no work under the current validation.\n',run_id);
    return
end

worker_log = fullfile(logs_dir,sprintf('worker_run%02d.log',run_id));
diary(worker_log);
diary_cleanup = onCleanup(@() diary('off'));
fprintf('Campaign: %s\nRun: %d\nNode: %s\n', ...
    campaign_tag,run_id,getenv('SLURMD_NODENAME'));

status_file = fullfile(logs_dir,sprintf('run_%02d.status',run_id));
Run_Subfield_Worker(run_id,status_file, ...
    'NDays',1826,'StartDate',datenum(2016,1,1), ...
    'AppeloAlpha',1.3,'ParallelWorkers',0,'SaveResults',true, ...
    'ResultDirectory',results_dir);

status = fileread(status_file);
if ~contains(status,'state=complete')
    error('Bouchet_Kantola_1D:IncompleteStatus', ...
        'Run %d did not write a complete status.',run_id);
end
remaining = dir(fullfile(results_dir,'checkpoints', ...
    sprintf('ckpt_run%02d_*.mat',run_id)));
if ~isempty(remaining)
    error('Bouchet_Kantola_1D:CheckpointRemains', ...
        'Run %d completed but retained a checkpoint.',run_id);
end

%% Jobs execute from the shared source tree, which stays writable while
%% they run.  MATLAB loads functions lazily, so an upload part-way through
%% a multi-day run can change the chemistry mid-trajectory.  Recompute the
%% gate and reject the result if anything moved.
if ~strcmp(Bouchet_Campaign_Validation_Fingerprint( ...
        fullfile(repo_root,'v2'),campaign_dir,campaign_tag,1,0),expected_stamp)
    error('Bouchet_Kantola_1D:SourceChangedDuringRun', ...
        ['Code, inputs, tests or numerical runtime changed while run %d ' ...
         'was executing. Discard its result.'],run_id);
end
if ~strcmp(Bouchet_Check_Existing_Result( ...
        results_dir,run_id,stamp_digest,1,1826),'complete')
    error('Bouchet_Kantola_1D:MissingResult', ...
        'Run %d reported completion but wrote no matching final MAT.',run_id);
end

end

function assert_single_thread_runtime
matlab_release = version('-release');
if ~strcmp(matlab_release,'2025a')
    error('Bouchet_Kantola_1D:WrongMatlabRelease', ...
        'This campaign requires MATLAB R2025a; loaded R%s.',matlab_release);
end
fprintf('MATLAB: R%s\n',matlab_release);
Bouchet_Single_Thread_Runtime('Bouchet_Kantola_1D');
end

function value = read_positive_integer(name)
raw = strtrim(getenv(name));
value = str2double(raw);
if isempty(raw) || ~(isfinite(value) && value == fix(value) && value >= 1)
    error('Bouchet_Kantola_1D:BadEnvironment', ...
        '%s must be a positive integer.',name);
end
end
