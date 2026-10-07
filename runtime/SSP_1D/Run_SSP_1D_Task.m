function Run_SSP_1D_Task
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run one Kantola 1-D treatment of an SSP campaign with the calibrated
%% CEC 13 parameter set.
%%
%% TC_SSP selects the forcing file, TC_HORIZON the run length, and both
%% enter the campaign tag, so one directory serves every scenario and both
%% horizons and each campaign gets its own results directory.  See the 2-D
%% array task for the full reasoning.
%%
%% Unlike the five-year campaign this has no sweep-identity report: runs 10
%% to 12 there recompute the sweep's cases under the same forcing and
%% horizon, and neither holds once the forcing is a scenario series.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[scenario,horizon,n_days,last_application_year] = read_scenario_and_horizon;
calcite_case = read_calcite_case;
campaign_tag = sprintf( ...
    'SSP_1D_%s_%s%s', ...
    horizon,scenario,case_suffix(calcite_case));
one_d_tag = ...
    'Base_1D';
campaign_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(campaign_dir));
one_d_dir = fullfile(repo_root,'runtime',one_d_tag);
run_dir = fullfile(repo_root,'v2','run');
code_dir = fullfile(repo_root,'v2','T&C_Code');
addpath(run_dir,code_dir,fullfile(code_dir,'Inputs'), ...
    fullfile(code_dir,'Parameters'),campaign_dir,one_d_dir);

assert_single_thread_runtime;
P = Campaign_Parameters;
%% The calibrated set was adopted on a five-year horizon.  Only the horizon
%% changes here; StartDate and every chemical parameter stay as the sweep
%% adopted them.
P.NDays = n_days;
run_id = read_positive_integer('TC_RUN_ID');
if run_id < 10 || run_id > 12
    error('Bouchet_Kantola_SSP:BadRunId', ...
        'TC_RUN_ID must select run 10, 11 or 12.');
end
if ~isempty(calcite_case) && run_id ~= 12
    error('Bouchet_Kantola_SSP:BadRunId', ...
        'A calcite-case campaign runs the calcite treatment, run 12, only.');
end
if read_positive_integer('SLURM_CPUS_PER_TASK') ~= 1
    error('Bouchet_Kantola_SSP:BadCpuCount', ...
        'Each production treatment must reserve exactly one CPU.');
end
if ~strcmp(strtrim(getenv('TC_PARALLEL_WORKERS')),'0')
    error('Bouchet_Kantola_SSP:BadWorkerCount', ...
        'TC_PARALLEL_WORKERS must be 0 for a one-column run.');
end
if ~strcmp(strtrim(getenv('TC_CKPT_MINUTES')),'30')
    error('Bouchet_Kantola_SSP:BadCheckpointInterval', ...
        'TC_CKPT_MINUTES must be 30.');
end
if ~isempty(strtrim(getenv('TC_CKPT_FORCE_RESUME')))
    error('Bouchet_Kantola_SSP:ForcedResumeForbidden', ...
        'TC_CKPT_FORCE_RESUME must remain unset for production.');
end
if ~strcmp(strtrim(getenv('TC_CAMPAIGN_TAG')),campaign_tag)
    error('Bouchet_Kantola_SSP:BadCampaignTag', ...
        'TC_CAMPAIGN_TAG must be %s.',campaign_tag);
end

required_files = { ...
    fullfile(code_dir,'Inputs','ERW_MOD_PARAM_parameters.xlsx'), ...
    fullfile(code_dir,'Inputs',sprintf('Data_Kantola_%s.mat',scenario)), ...
    fullfile(code_dir,'Inputs','Data_Kantola_BG.mat'), ...
    fullfile(code_dir,'Inputs','Ca_Data.mat'), ...
    fullfile(code_dir,'Inputs','All_deposition_data.mat'), ...
    fullfile(code_dir,'Inputs','blueridge_CDF_data.json'), ...
    fullfile(repo_root,'v2','SIM','EQUIV_1D_TMP', ...
        'TC_ERW_results_Kantola_crop_control.mat'), ...
    fullfile(one_d_dir,'Campaign_Parameters.m'), ...
    fullfile(one_d_dir,'Assert_Base_Provenance.m'), ...
    fullfile(one_d_dir,'Fingerprint_Digest.m')};
for k = 1:numel(required_files)
    if ~isfile(required_files{k})
        error('Bouchet_Kantola_SSP:MissingDeploymentFile', ...
            'Missing required deployment file: %s',required_files{k});
    end
end

%% Kept identical to the 2-D task so both campaigns of a scenario share one
%% definition of what the gate covers.
stamp_dependencies = { ...
    fullfile(one_d_dir,'Campaign_Parameters.m'), ...
    fullfile(one_d_dir,'Assert_Base_Provenance.m'), ...
    fullfile(one_d_dir,'Fingerprint_Digest.m'), ...
    fullfile(one_d_dir,'test_all_12_configuration.m'), ...
    fullfile(campaign_dir,'submit_ssp.sh')};

results_dir = fullfile(repo_root,'v2','SIM','Subfield_12',campaign_tag);
logs_dir = fullfile(repo_root,'v2','logs','Subfield_12',campaign_tag);
stamp_file = fullfile(logs_dir,'validation_stamp.txt');
if ~isfile(stamp_file)
    error('Bouchet_Kantola_SSP:MissingValidation', ...
        'Run the R2025a validation job before production.');
end
expected_stamp = Bouchet_Campaign_Validation_Fingerprint( ...
    fullfile(repo_root,'v2'),campaign_dir,campaign_tag,1,0, ...
    stamp_dependencies{:});
if ~strcmp(strtrim(fileread(stamp_file)),expected_stamp)
    error('Bouchet_Kantola_SSP:StaleValidation', ...
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
%% so a duplicate submission fails here rather than racing.  The cleanup
%% object must stay in scope for the whole run; it releases the lock when
%% this function returns or errors.
lock_cleanup = Bouchet_Acquire_Run_Lock( ...
    fullfile(logs_dir,'locks',sprintf('run_%02d.lock',run_id)));

if strcmp(Bouchet_Check_Existing_Result( ...
        results_dir,run_id,stamp_digest,1,P.NDays),'complete')
    fprintf('Run %d needs no work under the current validation.\n',run_id);
    return
end

worker_log = fullfile(logs_dir,sprintf('worker_run%02d.log',run_id));
diary(worker_log);
diary_cleanup = onCleanup(@() diary('off'));
fprintf(['Campaign: %s\nScenario: %s\nHorizon: %s (%d days)\n' ...
    'Run: %d\nAppelo alpha: %.17g\n' ...
    'KNmaxd: %.17g d^-1\nKbot: %.17g mm h^-1\nCEC: %g cmolc kg^-1\n' ...
    'Node: %s\n'],campaign_tag,scenario,horizon,P.NDays,run_id, ...
    P.AppeloAlpha,P.KNmaxd,P.Kbot,P.CEC_cmolc_kg,getenv('SLURMD_NODENAME'));

status_file = fullfile(logs_dir,sprintf('run_%02d.status',run_id));
run_options = {};
if ~isempty(calcite_case)
    run_options = [run_options {'CalciteCase',calcite_case}];
end
if ~isempty(last_application_year)
    run_options = [run_options {'LastApplicationYear',last_application_year}];
end
Run_Subfield_Worker(run_id,status_file, ...
    'NDays',P.NDays,'StartDate',P.StartDate, ...
    'AppeloAlpha',P.AppeloAlpha,'AppeloLogKH',P.AppeloLogKH, ...
    'KNmaxd',P.KNmaxd,'Kbot',P.Kbot, ...
    'ParallelWorkers',0,'SaveResults',true, ...
    'ResultDirectory',results_dir,'Scenario',scenario,run_options{:});

status = fileread(status_file);
if ~contains(status,'state=complete')
    error('Bouchet_Kantola_SSP:IncompleteStatus', ...
        'Run %d did not write a complete status.',run_id);
end
remaining = dir(fullfile(results_dir,'checkpoints', ...
    sprintf('ckpt_run%02d_*.mat',run_id)));
if ~isempty(remaining)
    error('Bouchet_Kantola_SSP:CheckpointRemains', ...
        'Run %d completed but retained a checkpoint.',run_id);
end

%% Jobs execute from the shared source tree, which stays writable while
%% they run.  MATLAB loads functions lazily, so an upload part-way through
%% a multi-day run can change the chemistry mid-trajectory.  Recompute the
%% gate and reject the result if anything moved.
if ~strcmp(Bouchet_Campaign_Validation_Fingerprint( ...
        fullfile(repo_root,'v2'),campaign_dir,campaign_tag,1,0, ...
        stamp_dependencies{:}),expected_stamp)
    error('Bouchet_Kantola_SSP:SourceChangedDuringRun', ...
        ['Code, inputs, tests or numerical runtime changed while run %d ' ...
         'was executing. Discard its result.'],run_id);
end
if ~strcmp(Bouchet_Check_Existing_Result( ...
        results_dir,run_id,stamp_digest,1,P.NDays),'complete')
    error('Bouchet_Kantola_SSP:MissingResult', ...
        'Run %d reported completion but wrote no matching final MAT.',run_id);
end

%% The saved provenance must carry the calibrated identity exactly.
result = dir(fullfile(results_dir,sprintf('TC_2D_run%02d_*.mat',run_id)));
result = result(~[result.isdir]);
result_file = fullfile(result(1).folder,result(1).name);
saved = load(result_file,'run_provenance');
Assert_Base_Provenance(saved.run_provenance,run_id,1,0, ...
    stamp_digest,campaign_tag,results_dir,P.NDays);
fprintf('Provenance: run %d carries the calibrated CEC 13 identity.\n',run_id);
fprintf('Run %d complete under %s.\n',run_id,campaign_tag);
end

function [scenario,horizon,n_days,last_application_year] = read_scenario_and_horizon
scenario = strtrim(getenv('TC_SSP'));
if ~ismember(scenario,{'ssp126','ssp245','ssp585'})
    error('Bouchet_Kantola_SSP:BadScenario', ...
        ['TC_SSP must be ssp126, ssp245 or ssp585, not ''%s''. It ' ...
        'selects the forcing file Data_Kantola_<scenario>.mat, and is not ' ...
        'defaulted because a default would silently run the historical ' ...
        'forcing under a scenario name.'],scenario);
end
horizon = strtrim(getenv('TC_HORIZON'));
%% 2016-01-01 is the start date the calibrated set fixes.  5yr reaches
%% 2020-12-31 (1826 days, one leap day) and is the span with observations to
%% validate against; 35yr reaches 2050-12-31 (12784 days, nine leap days),
%% which is the full length of the delta forcing.
%% The 35-year projections keep applying rock every 10 November to the
%% end of the forcing (user decision 2026-09-24); the five-year
%% validation keeps the field record, 2016-2019, which is the default.
switch horizon
    case '5yr'
        n_days = 1826;
        last_application_year = [];
    case '35yr'
        n_days = 12784;
        last_application_year = 2050;
    otherwise
        error('Bouchet_Kantola_SSP:BadHorizon', ...
            'TC_HORIZON must be 5yr or 35yr, not ''%s''.',horizon);
end
end

function calcite_case = read_calcite_case
%% Optional campaign axis: an ERW workbook row for the calcite runs only,
%% such as 8, the agronomic 5 t/ha rate. Unset is the standard campaign.
raw = strtrim(getenv('TC_CALCITE_CASE'));
calcite_case = [];
if ~isempty(raw)
    calcite_case = str2double(raw);
    if ~(isfinite(calcite_case) && calcite_case == fix(calcite_case) && calcite_case >= 1)
        error('Bouchet_Kantola_SSP:BadCalciteCase', ...
            'TC_CALCITE_CASE must be unset or a positive integer.');
    end
end
end

function suffix = case_suffix(calcite_case)
suffix = '';
if ~isempty(calcite_case)
    suffix = sprintf('_case%d',calcite_case);
end
end

function assert_single_thread_runtime
matlab_release = version('-release');
if ~strcmp(matlab_release,'2025a')
    error('Bouchet_Kantola_SSP:WrongMatlabRelease', ...
        'This campaign requires MATLAB R2025a; loaded R%s.',matlab_release);
end
fprintf('MATLAB: R%s\n',matlab_release);
Bouchet_Single_Thread_Runtime('Bouchet_Kantola_SSP');
end

function value = read_positive_integer(name)
raw = strtrim(getenv(name));
value = str2double(raw);
if isempty(raw) || ~(isfinite(value) && value == fix(value) && value >= 1)
    error('Bouchet_Kantola_SSP:BadEnvironment', ...
        '%s must be a positive integer.',name);
end
end
