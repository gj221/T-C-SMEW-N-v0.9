function Run_Base_1D_Task
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run one five-year Kantola 1-D treatment with the calibrated CEC 13 set.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

campaign_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(campaign_dir));
run_dir = fullfile(repo_root,'v2','run');
code_dir = fullfile(repo_root,'v2','T&C_Code');
addpath(run_dir,code_dir,fullfile(code_dir,'Inputs'), ...
    fullfile(code_dir,'Parameters'),campaign_dir);
[calcite_case,case_suffix] = Calcite_Case;
campaign_tag = ['Base_1D' case_suffix];

assert_single_thread_runtime;
P = Campaign_Parameters;
run_id = read_positive_integer('TC_RUN_ID');
if run_id < 10 || run_id > 12
    error('Bouchet_Kantola_CEC13:BadRunId', ...
        'TC_RUN_ID must select run 10, 11 or 12.');
end
if ~isempty(calcite_case) && run_id ~= 12
    error('Bouchet_Kantola_CEC13:BadRunId', ...
        'A calcite-case campaign runs the calcite treatment, run 12, only.');
end
if read_positive_integer('SLURM_CPUS_PER_TASK') ~= 1
    error('Bouchet_Kantola_CEC13:BadCpuCount', ...
        'Each production treatment must reserve exactly one CPU.');
end
if ~strcmp(strtrim(getenv('TC_PARALLEL_WORKERS')),'0')
    error('Bouchet_Kantola_CEC13:BadWorkerCount', ...
        'TC_PARALLEL_WORKERS must be 0 for a one-column run.');
end
if ~strcmp(strtrim(getenv('TC_CKPT_MINUTES')),'30')
    error('Bouchet_Kantola_CEC13:BadCheckpointInterval', ...
        'TC_CKPT_MINUTES must be 30.');
end
if ~isempty(strtrim(getenv('TC_CKPT_FORCE_RESUME')))
    error('Bouchet_Kantola_CEC13:ForcedResumeForbidden', ...
        'TC_CKPT_FORCE_RESUME must remain unset for production.');
end
if ~strcmp(strtrim(getenv('TC_CAMPAIGN_TAG')),campaign_tag)
    error('Bouchet_Kantola_CEC13:BadCampaignTag', ...
        'TC_CAMPAIGN_TAG must be %s.',campaign_tag);
end

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
        error('Bouchet_Kantola_CEC13:MissingDeploymentFile', ...
            'Missing required deployment file: %s',required_files{k});
    end
end

results_dir = fullfile(repo_root,'v2','SIM','Subfield_12',campaign_tag);
logs_dir = fullfile(repo_root,'v2','logs','Subfield_12',campaign_tag);
stamp_file = fullfile(logs_dir,'validation_stamp.txt');
if ~isfile(stamp_file)
    error('Bouchet_Kantola_CEC13:MissingValidation', ...
        'Run the R2025a validation job before production.');
end
expected_stamp = Bouchet_Campaign_Validation_Fingerprint( ...
    fullfile(repo_root,'v2'),campaign_dir,campaign_tag,1,0);
if ~strcmp(strtrim(fileread(stamp_file)),expected_stamp)
    error('Bouchet_Kantola_CEC13:StaleValidation', ...
        'Code, inputs, tests or numerical runtime changed after validation.');
end
%% The stamp this process verified travels into the saved result so a later
%% job can tell which validation produced each MAT.  It travels as a digest
%% because the full fingerprint is far past setenv's 32766-character limit.
stamp_digest = Fingerprint_Digest(expected_stamp);
setenv('TC_VALIDATION_FINGERPRINT',stamp_digest);
if ~exist(results_dir,'dir'); mkdir(results_dir); end
if ~exist(logs_dir,'dir'); mkdir(logs_dir); end

%% A missing final MAT only means that no run has finished.  Take the lock
%% before anything else touches the status file, checkpoint or result path,
%% so a duplicate submission fails here rather than racing.
lock_cleanup = Bouchet_Acquire_Run_Lock( ...
    fullfile(logs_dir,'locks',sprintf('run_%02d.lock',run_id))); %#ok<NASGU>

if strcmp(Bouchet_Check_Existing_Result( ...
        results_dir,run_id,stamp_digest,1,P.NDays),'complete')
    fprintf('Run %d needs no work under the current validation.\n',run_id);
    return
end

worker_log = fullfile(logs_dir,sprintf('worker_run%02d.log',run_id));
diary(worker_log);
diary_cleanup = onCleanup(@() diary('off'));
fprintf(['Campaign: %s\nRun: %d\nAppelo alpha: %.17g\n' ...
    'KNmaxd: %.17g d^-1\nKbot: %.17g mm h^-1\nCEC: %g cmolc kg^-1\n' ...
    'Node: %s\n'],campaign_tag,run_id,P.AppeloAlpha,P.KNmaxd,P.Kbot, ...
    P.CEC_cmolc_kg,getenv('SLURMD_NODENAME'));

status_file = fullfile(logs_dir,sprintf('run_%02d.status',run_id));
run_options = {};
if ~isempty(calcite_case)
    run_options = {'CalciteCase',calcite_case};
end
Run_Subfield_Worker(run_id,status_file, ...
    'NDays',P.NDays,'StartDate',P.StartDate, ...
    'AppeloAlpha',P.AppeloAlpha,'AppeloLogKH',P.AppeloLogKH, ...
    'KNmaxd',P.KNmaxd,'Kbot',P.Kbot, ...
    'ParallelWorkers',0,'SaveResults',true, ...
    'ResultDirectory',results_dir,run_options{:});

status = fileread(status_file);
if ~contains(status,'state=complete')
    error('Bouchet_Kantola_CEC13:IncompleteStatus', ...
        'Run %d did not write a complete status.',run_id);
end
remaining = dir(fullfile(results_dir,'checkpoints', ...
    sprintf('ckpt_run%02d_*.mat',run_id)));
if ~isempty(remaining)
    error('Bouchet_Kantola_CEC13:CheckpointRemains', ...
        'Run %d completed but retained a checkpoint.',run_id);
end

%% Jobs execute from the shared source tree, which stays writable while
%% they run.  MATLAB loads functions lazily, so an upload part-way through
%% a multi-day run can change the chemistry mid-trajectory.  Recompute the
%% gate and reject the result if anything moved.
if ~strcmp(Bouchet_Campaign_Validation_Fingerprint( ...
        fullfile(repo_root,'v2'),campaign_dir,campaign_tag,1,0),expected_stamp)
    error('Bouchet_Kantola_CEC13:SourceChangedDuringRun', ...
        ['Code, inputs, tests or numerical runtime changed while run %d ' ...
         'was executing. Discard its result.'],run_id);
end
if ~strcmp(Bouchet_Check_Existing_Result( ...
        results_dir,run_id,stamp_digest,1,P.NDays),'complete')
    error('Bouchet_Kantola_CEC13:MissingResult', ...
        'Run %d reported completion but wrote no matching final MAT.',run_id);
end

%% The saved provenance must carry the calibrated identity exactly.
result = dir(fullfile(results_dir,sprintf('TC_2D_run%02d_*.mat',run_id)));
result = result(~[result.isdir]);
result_file = fullfile(result(1).folder,result(1).name);
saved = load(result_file,'run_provenance');
Assert_Base_Provenance(saved.run_provenance,run_id,1,0, ...
    stamp_digest,campaign_tag,results_dir,P.NDays);
%% The provenance does not record the workbook row; the result name does.
if ~isempty(calcite_case) && ~contains(result(1).name,['_' case_suffix '__'])
    error('Bouchet_Kantola_CEC13:WrongCalciteCase', ...
        'Run %d was not saved as calcite case %d: %s',run_id,calcite_case,result(1).name);
end
fprintf('Provenance: run %d carries the calibrated CEC 13 identity.\n',run_id);

%% Runs 10 to 12 recompute the sweep's cases P.sweep_case_ids under the same
%% numerical identity, so their daily states should match the sweep result
%% bit for bit.  Report the comparison; a difference is recorded, not
%% treated as a failure of this run.  The sweep ran no other calcite case.
if isempty(calcite_case)
    report_sweep_identity(repo_root,P,run_id,result_file,logs_dir);
end
end

function report_sweep_identity(repo_root,P,run_id,result_file,logs_dir)
sweep_dir = fullfile(repo_root,'v2','SIM','Subfield_12',P.sweep_campaign, ...
    sprintf('%s_run%d',P.sweep_combination_tag,run_id));
[~,name,ext] = fileparts(result_file);
sweep_file = fullfile(sweep_dir,[name ext]);
note_file = fullfile(logs_dir,sprintf('sweep_identity_run%02d.txt',run_id));
if ~isfile(sweep_file)
    text = sprintf('sweep result not present: %s\n',sweep_file);
else
    fields = {'PH_bulk_d','PH_pore_d','N2Oflx','LEAK_NO3','P','Lk'};
    a = load(result_file,fields{:});
    b = load(sweep_file,fields{:});
    same = cellfun(@(f) isequaln(a.(f),b.(f)),fields);
    if all(same)
        text = sprintf('identical to sweep case %d: %s\n', ...
            P.sweep_case_ids(run_id-9),sweep_file);
    else
        text = sprintf('DIFFERS from sweep case %d in %s: %s\n', ...
            P.sweep_case_ids(run_id-9),strjoin(fields(~same),', '), ...
            sweep_file);
    end
end
fid = fopen(note_file,'w');
if fid > 0
    fprintf(fid,'%s',text);
    fclose(fid);
end
fprintf('Sweep identity: %s',text);
end

function assert_single_thread_runtime
matlab_release = version('-release');
if ~strcmp(matlab_release,'2025a')
    error('Bouchet_Kantola_CEC13:WrongMatlabRelease', ...
        'This campaign requires MATLAB R2025a; loaded R%s.',matlab_release);
end
fprintf('MATLAB: R%s\n',matlab_release);
Bouchet_Single_Thread_Runtime('Bouchet_Kantola_CEC13');
end

function value = read_positive_integer(name)
raw = strtrim(getenv(name));
value = str2double(raw);
if isempty(raw) || ~(isfinite(value) && value == fix(value) && value >= 1)
    error('Bouchet_Kantola_CEC13:BadEnvironment', ...
        '%s must be a positive integer.',name);
end
end
