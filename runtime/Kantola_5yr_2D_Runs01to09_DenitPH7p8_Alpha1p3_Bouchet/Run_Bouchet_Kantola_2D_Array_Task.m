function Run_Bouchet_Kantola_2D_Array_Task
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run one Slurm array element for the five-year Kantola 2-D campaign.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

campaign_tag = 'Kantola_5yr_2D_Runs01to09_DenitPH7p8_Alpha1p3_Bouchet';
campaign_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(campaign_dir));
run_dir = fullfile(repo_root,'v2','run');
code_dir = fullfile(repo_root,'v2','T&C_Code');
addpath(run_dir,code_dir,fullfile(code_dir,'Inputs'),fullfile(code_dir,'Parameters'));

matlab_release = version('-release');
if ~strcmp(matlab_release,'2025a')
    error('Bouchet_Kantola_2D:WrongMatlabRelease', ...
        'This campaign requires MATLAB R2025a; loaded R%s.',matlab_release);
end
Bouchet_Single_Thread_Runtime('Bouchet_Kantola_2D');

required_files = { ...
    fullfile(code_dir,'Inputs','ERW_MOD_PARAM_parameters.xlsx'), ...
    fullfile(code_dir,'Inputs','Data_Kantola.mat'), ...
    fullfile(code_dir,'Inputs','Data_Kantola_BG.mat'), ...
    fullfile(code_dir,'Inputs','Ca_Data.mat'), ...
    fullfile(code_dir,'Inputs','All_deposition_data.mat'), ...
    fullfile(code_dir,'Inputs','blueridge_CDF_data.json'), ...
    fullfile(repo_root,'v2','SIM','EQUIV_1D_TMP', ...
        'TC_ERW_results_Kantola_crop_control.mat')};
for k = 1:numel(required_files)
    if ~isfile(required_files{k})
        error('Bouchet_Kantola_2D:MissingDeploymentFile', ...
            'Missing required deployment file: %s',required_files{k});
    end
end

run_id = read_positive_integer('SLURM_ARRAY_TASK_ID');
if run_id > 9
    error('Bouchet_Kantola_2D:BadRunId', ...
        'SLURM_ARRAY_TASK_ID must select run 1 through 9.');
end
allocated_cpus = read_positive_integer('SLURM_CPUS_PER_TASK');
workers = Resolve_2D_Parfor_Workers(12);
if workers ~= 6
    error('Bouchet_Kantola_2D:BadWorkerCount', ...
        'TC_PARALLEL_WORKERS must be 6 for this campaign.');
end
if allocated_cpus < workers
    error('Bouchet_Kantola_2D:InsufficientCpus', ...
        'Reserve at least %d CPUs.',workers);
end
if ~strcmp(getenv('TC_CKPT_MINUTES'),'10')
    error('Bouchet_Kantola_2D:BadCheckpointInterval', ...
        'TC_CKPT_MINUTES must be 10.');
end
if ~isempty(getenv('TC_CKPT_FORCE_RESUME'))
    error('Bouchet_Kantola_2D:ForcedResumeForbidden', ...
        'TC_CKPT_FORCE_RESUME must remain unset for production.');
end
if ~license('test','Distrib_Computing_Toolbox')
    error('Bouchet_Kantola_2D:ParallelToolboxUnavailable', ...
        'Parallel Computing Toolbox is required for this campaign.');
end
validate_only = strtrim(getenv('TC_VALIDATE_ONLY'));
if ~isempty(validate_only) && ~strcmp(validate_only,'1')
    error('Bouchet_Kantola_2D:BadValidationMode', ...
        'TC_VALIDATE_ONLY must be unset or 1.');
end
if strcmp(validate_only,'1')
    run_validation(code_dir,repo_root,campaign_dir,campaign_tag);
    return
end

results_dir = fullfile(repo_root,'v2','SIM','Subfield_12',campaign_tag);
logs_dir = fullfile(repo_root,'v2','logs','Subfield_12',campaign_tag);
stamp_file = fullfile(logs_dir,'validation_stamp.txt');
if ~isfile(stamp_file)
    error('Bouchet_Kantola_2D:MissingValidation', ...
        'Run the R2025a validation job before production.');
end
expected_stamp = Bouchet_Campaign_Validation_Fingerprint( ...
    fullfile(repo_root,'v2'),campaign_dir,campaign_tag,12,6, ...
    fullfile(repo_root,'runtime', ...
        'Kantola_5yr_1D_Runs10to12_ScheduleFix_R2025a_Bouchet', ...
        'Validate_Bouchet_1D_Preapplication_Identity.m'));
if ~strcmp(strtrim(fileread(stamp_file)),expected_stamp)
    error('Bouchet_Kantola_2D:StaleValidation', ...
        'Code, inputs, tests or numerical runtime changed after validation.');
end
%% The stamp this process verified travels into the saved result so a later
%% job can tell which validation produced each MAT.
setenv('TC_VALIDATION_FINGERPRINT',expected_stamp);
if ~exist(results_dir,'dir'); mkdir(results_dir); end
if ~exist(logs_dir,'dir'); mkdir(logs_dir); end

%% A missing final MAT only means that no run has finished.  Take the lock
%% before anything else touches the status file, checkpoint or result path,
%% so a requeued or duplicate array element fails here rather than racing.
lock_cleanup = Bouchet_Acquire_Run_Lock( ...
    fullfile(logs_dir,'locks',sprintf('run_%02d.lock',run_id))); %#ok<NASGU>

if strcmp(Bouchet_Check_Existing_Result( ...
        results_dir,run_id,expected_stamp,12,1826),'complete')
    fprintf('Run %d needs no work under the current validation.\n',run_id);
    return
end

worker_log = fullfile(logs_dir,sprintf('worker_run%02d.log',run_id));
diary(worker_log);
diary_cleanup = onCleanup(@() diary('off'));
fprintf('Campaign: %s\nSlurm run: %d\nWorkers: %d\n', ...
    campaign_tag,run_id,workers);

status_file = fullfile(logs_dir,sprintf('run_%02d.status',run_id));
Run_Subfield_Worker(run_id,status_file, ...
    'NDays',1826,'StartDate',datenum(2016,1,1), ...
    'GridNx',2,'GridNy',2,'AppeloAlpha',1.3, ...
    'ParallelWorkers',workers,'SaveResults',true, ...
    'ResultDirectory',results_dir);

status = fileread(status_file);
if ~contains(status,'state=complete')
    error('Bouchet_Kantola_2D:IncompleteStatus', ...
        'Run %d did not write a complete status.',run_id);
end

%% Jobs execute from the shared source tree, which stays writable while
%% they run.  MATLAB loads functions lazily, so an upload part-way through
%% a multi-day run can change the chemistry mid-trajectory.  Recompute the
%% gate and reject the result if anything moved.
if ~strcmp(Bouchet_Campaign_Validation_Fingerprint( ...
        fullfile(repo_root,'v2'),campaign_dir,campaign_tag,12,6, ...
        fullfile(repo_root,'runtime', ...
            'Kantola_5yr_1D_Runs10to12_ScheduleFix_R2025a_Bouchet', ...
            'Validate_Bouchet_1D_Preapplication_Identity.m')),expected_stamp)
    error('Bouchet_Kantola_2D:SourceChangedDuringRun', ...
        ['Code, inputs, tests or numerical runtime changed while run %d ' ...
         'was executing. Discard its result.'],run_id);
end
if ~strcmp(Bouchet_Check_Existing_Result( ...
        results_dir,run_id,expected_stamp,12,1826),'complete')
    error('Bouchet_Kantola_2D:MissingResult', ...
        'Run %d reported completion but wrote no matching final MAT.',run_id);
end

end

function run_validation(code_dir,repo_root,campaign_dir,campaign_tag)
test_dir = fullfile(code_dir,'Tests');
one_d_runtime = fullfile(repo_root,'runtime', ...
    'Kantola_5yr_1D_Runs10to12_ScheduleFix_R2025a_Bouchet');
addpath(test_dir,one_d_runtime);
old_dir = pwd;
cleanup = onCleanup(@() cd(old_dir));
cd(test_dir);
%% Hashing the test files proves only that they exist.  The 2-D gate runs
%% the same authoritative suite as the 1-D gate, so chemistry, nutrient,
%% conservation and error-path coverage stands behind both campaigns.
run_post_legacy_equation_audit_tests(1);
test_biogeochemistry_io_path_portability;
test_bouchet_linux_portability(campaign_dir);
Validate_Bouchet_1D_Preapplication_Identity;
logs_dir = fullfile(repo_root,'v2','logs','Subfield_12',campaign_tag);
fingerprint = Bouchet_Campaign_Validation_Fingerprint( ...
    fullfile(repo_root,'v2'),campaign_dir,campaign_tag,12,6, ...
    fullfile(one_d_runtime,'Validate_Bouchet_1D_Preapplication_Identity.m'));
Bouchet_Write_Validation_Stamp( ...
    fullfile(logs_dir,'validation_stamp.txt'),fingerprint);
fprintf('PASS: Bouchet R2025a parallel 2-D validation\n');
end

function value = read_positive_integer(name)
raw = strtrim(getenv(name));
value = str2double(raw);
if isempty(raw) || ~(isfinite(value) && value == fix(value) && value >= 1)
    error('Bouchet_Kantola_2D:BadEnvironment', ...
        '%s must be a positive integer.',name);
end
end
