function Run_Sweep_Task
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run one five-year joint-sweep trajectory on Bouchet.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

campaign_tag = ...
    'Calibration_Sweep_1D';
campaign_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(campaign_dir));
v2_dir = fullfile(repo_root,'v2');
run_dir = fullfile(v2_dir,'run');
code_dir = fullfile(v2_dir,'T&C_Code');
addpath(run_dir,code_dir,fullfile(code_dir,'Inputs'), ...
    fullfile(code_dir,'Parameters'),campaign_dir);

assert_single_thread_runtime;
array_task_id = read_positive_integer('SLURM_ARRAY_TASK_ID');
sweep_layout = strtrim(getenv('TC_SWEEP_LAYOUT'));
if isempty(sweep_layout); sweep_layout = 'independent'; end
sweep_scope = strtrim(getenv('TC_SWEEP_SCOPE'));
if isempty(sweep_scope); sweep_scope = 'full'; end
if ~ismember(sweep_scope,{'full','pre_control'})
    error('Bouchet_Joint_Sweep:BadScope', ...
        'TC_SWEEP_SCOPE must be full or pre_control.');
end
case_id_raw = strtrim(getenv('TC_SWEEP_CASE_ID'));
if isempty(case_id_raw)
    task_id = array_task_id;
else
    task_id = read_positive_integer('TC_SWEEP_CASE_ID');
end
cases = Sweep_Cases;
if numel(cases) ~= 336
    error('Bouchet_Joint_Sweep:CaseCount', ...
        'The deployed case table must contain exactly 336 trajectories.');
end
if task_id > numel(cases)
    error('Bouchet_Joint_Sweep:BadTaskId', ...
        'SLURM_ARRAY_TASK_ID must be 1..%d.',numel(cases));
end
case_cfg = cases(task_id);
if strcmp(sweep_layout,'independent')
    valid_layout = task_id == array_task_id;
elseif strcmp(sweep_layout,'packed')
    valid_layout = case_cfg.combination_id == array_task_id;
else
    valid_layout = false;
end
if ~valid_layout
    error('Bouchet_Joint_Sweep:BadLayout', ...
        'The array task, logical case and sweep layout are inconsistent.');
end
start_date = datenum(2016,1,1); %#ok<DATNM>
if strcmp(sweep_scope,'pre_control')
    if case_cfg.id ~= 1 || ~strcmp(sweep_layout,'independent') || ...
            ~isempty(case_id_raw)
        error('Bouchet_Joint_Sweep:BadPreControlCase', ...
            'pre_control requires independent case 1 only.');
    end
    n_days = 366;
else
    n_days = 1826;
end
if read_positive_integer('SLURM_CPUS_PER_TASK') ~= 1
    error('Bouchet_Joint_Sweep:BadCpuCount', ...
        'Each production trajectory must reserve exactly one CPU.');
end
if ~strcmp(strtrim(getenv('TC_PARALLEL_WORKERS')),'0')
    error('Bouchet_Joint_Sweep:BadWorkerCount', ...
        'TC_PARALLEL_WORKERS must be 0 for a one-column run.');
end
if ~strcmp(strtrim(getenv('TC_CKPT_MINUTES')),'30')
    error('Bouchet_Joint_Sweep:BadCheckpointInterval', ...
        'TC_CKPT_MINUTES must be 30.');
end
if ~isempty(strtrim(getenv('TC_CKPT_FORCE_RESUME')))
    error('Bouchet_Joint_Sweep:ForcedResumeForbidden', ...
        'TC_CKPT_FORCE_RESUME must remain unset for production.');
end
if ~strcmp(strtrim(getenv('TC_CAMPAIGN_TAG')),campaign_tag)
    error('Bouchet_Joint_Sweep:CampaignIdentity', ...
        'TC_CAMPAIGN_TAG does not identify this sweep.');
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
        error('Bouchet_Joint_Sweep:MissingDeploymentFile', ...
            'Missing required deployment file: %s',required_files{k});
    end
end

campaign_results_root = fullfile(v2_dir,'SIM','Subfield_12',campaign_tag);
campaign_logs_root = fullfile(v2_dir,'logs','Subfield_12',campaign_tag);
stamp_file = fullfile(campaign_logs_root,'validation_stamp.txt');
if ~isfile(stamp_file)
    error('Bouchet_Joint_Sweep:MissingValidation', ...
        'Run the R2025a validation job before production.');
end
expected_stamp = Sweep_Fingerprint( ...
    v2_dir,campaign_dir,campaign_tag);
if ~strcmp(strtrim(fileread(stamp_file)),expected_stamp)
    error('Bouchet_Joint_Sweep:StaleValidation', ...
        'Code, inputs, tests or numerical runtime changed after validation.');
end
stamp_digest = Fingerprint_Digest(expected_stamp);
setenv('TC_VALIDATION_FINGERPRINT',stamp_digest);

if strcmp(sweep_scope,'pre_control')
    results_root = fullfile(campaign_results_root,'pre_control_gate');
    logs_root = fullfile(campaign_logs_root,'pre_control_gate');
else
    results_root = campaign_results_root;
    logs_root = campaign_logs_root;
end

results_dir = fullfile(results_root,case_cfg.tag);
logs_dir = fullfile(logs_root,case_cfg.tag);
if ~exist(results_dir,'dir'); mkdir(results_dir); end
if ~exist(logs_dir,'dir'); mkdir(logs_dir); end
status_file = fullfile(logs_dir,sprintf( ...
    'run_%02d.status',case_cfg.run_id));
lock_cleanup = Bouchet_Acquire_Run_Lock(fullfile(logs_dir,'locks', ...
    sprintf('run_%02d.lock',case_cfg.run_id))); %#ok<NASGU>

if strcmp(Bouchet_Check_Existing_Result(results_dir,case_cfg.run_id, ...
        stamp_digest,1,n_days),'complete')
    result_file = assert_case_provenance(results_dir,case_cfg,stamp_digest, ...
        campaign_tag,n_days,start_date);
    recovered = dir(fullfile(results_dir,'checkpoints', ...
        sprintf('ckpt_run%02d_*.mat',case_cfg.run_id)));
    for k = 1:numel(recovered)
        delete(fullfile(recovered(k).folder,recovered(k).name));
    end
    write_recovered_status(status_file,result_file);
    fprintf('Sweep task %d (%s) needs no work.\n',task_id,case_cfg.tag);
    return
end

worker_log = fullfile(logs_dir,sprintf( ...
    'worker_run%02d.log',case_cfg.run_id));
diary(worker_log);
diary_cleanup = onCleanup(@() diary('off'));
if isnan(case_cfg.Kbot); kbot_text = 'free drainage';
else;                    kbot_text = sprintf('%.17g mm h^-1',case_cfg.Kbot);
end
fprintf(['Campaign: %s\nScope: %s\nTask: %d\nCase: %s\nRun: %d\n' ...
    'Layout: %s\nCEC: %.17g cmolc kg^-1\nKNmaxd: %.17g d^-1\n' ...
    'Appelo alpha: %.17g\nKbot: %s\nDissolution factor: %.17g\nNode: %s\n'], ...
    campaign_tag,sweep_scope,task_id,case_cfg.tag,case_cfg.run_id, ...
    sweep_layout, ...
    case_cfg.CEC_cmolc_kg, ...
    case_cfg.KNmaxd, ...
    case_cfg.AppeloAlpha,kbot_text,case_cfg.DissolutionFactor, ...
    getenv('SLURMD_NODENAME'));

Run_Subfield_Worker(case_cfg.run_id,status_file, ...
    'NDays',n_days,'StartDate',start_date, ...
    'AppeloAlpha',case_cfg.AppeloAlpha, ...
    'AppeloLogKH',case_cfg.AppeloLogKH, ...
    'KNmaxd',case_cfg.KNmaxd,'Kbot',case_cfg.Kbot, ...
    'DissolutionFactor',case_cfg.DissolutionFactor, ...
    'ParallelWorkers',0,'SaveResults',true, ...
    'ResultDirectory',results_dir);

status = regexp(fileread(status_file),'\r?\n','split');
if isempty(status) || ~strcmp(strtrim(status{1}),'state=complete')
    error('Bouchet_Joint_Sweep:IncompleteStatus', ...
        'Case %s did not write a complete status.',case_cfg.tag);
end
remaining = dir(fullfile(results_dir,'checkpoints', ...
    sprintf('ckpt_run%02d_*.mat',case_cfg.run_id)));
if ~isempty(remaining)
    error('Bouchet_Joint_Sweep:CheckpointRemains', ...
        'Case %s completed but retained a checkpoint.',case_cfg.tag);
end
if ~strcmp(Sweep_Fingerprint( ...
        v2_dir,campaign_dir,campaign_tag),expected_stamp)
    error('Bouchet_Joint_Sweep:SourceChangedDuringRun', ...
        'Code, inputs, tests or numerical runtime changed during case %s.', ...
        case_cfg.tag);
end
if ~strcmp(Bouchet_Check_Existing_Result(results_dir,case_cfg.run_id, ...
        stamp_digest,1,n_days),'complete')
    error('Bouchet_Joint_Sweep:MissingResult', ...
        'Case %s reported completion but wrote no matching MAT.',case_cfg.tag);
end
assert_case_provenance(results_dir,case_cfg,stamp_digest, ...
    campaign_tag,n_days,start_date);

end

function result_file = assert_case_provenance( ...
        results_dir,case_cfg,stamp_digest,campaign_tag,n_days,start_date)
files = dir(fullfile(results_dir,sprintf( ...
    'TC_2D_run%02d_*.mat',case_cfg.run_id)));
files = files(~[files.isdir]);
if numel(files) ~= 1
    error('Bouchet_Joint_Sweep:AmbiguousCaseResult', ...
        'Case %s must contain exactly one result.',case_cfg.tag);
end
result_file = fullfile(files.folder,files.name);
S = load(result_file,'run_provenance');
if ~isfield(S,'run_provenance')
    error('Bouchet_Joint_Sweep:MissingCaseProvenance', ...
        'Result %s has no run_provenance.',result_file);
end
Assert_Sweep_Provenance(S.run_provenance,case_cfg, ...
    stamp_digest,campaign_tag,n_days,start_date,results_dir);
end

function write_recovered_status(status_file,result_file)
tmp_file = [status_file '.tmp'];
fid = fopen(tmp_file,'w');
if fid < 0
    error('Bouchet_Joint_Sweep:StatusWrite', ...
        'Cannot write recovered status %s.',tmp_file);
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid,'state=complete\n');
fprintf(fid,'detail=%s\n',strrep(result_file,newline,' '));
fprintf(fid,'message=recovered from validated result\n');
fprintf(fid,'finished=%s\n',datestr(now,30));
clear cleanup
[ok,message] = movefile(tmp_file,status_file,'f');
if ~ok
    if isfile(tmp_file); delete(tmp_file); end
    error('Bouchet_Joint_Sweep:StatusWrite','%s',message);
end
end

function assert_single_thread_runtime
matlab_release = version('-release');
if ~strcmp(matlab_release,'2025a')
    error('Bouchet_Joint_Sweep:WrongMatlabRelease', ...
        'This campaign requires MATLAB R2025a; loaded R%s.',matlab_release);
end
fprintf('MATLAB: R%s\n',matlab_release);
Bouchet_Single_Thread_Runtime('Bouchet_Joint_Sweep');
end

function value = read_positive_integer(name)
raw = strtrim(getenv(name));
value = str2double(raw);
if isempty(raw) || ~(isfinite(value) && value == fix(value) && value >= 1)
    error('Bouchet_Joint_Sweep:BadEnvironment', ...
        '%s must be a positive integer.',name);
end
end
