function Run_SSP_2D_Task
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run one Slurm array element of a Kantola 2-D SSP campaign with the
%% calibrated CEC 13 parameter set.
%%
%% One directory serves every scenario and both horizons.  TC_SSP
%% selects the forcing file, TC_HORIZON the run length, and both enter the
%% campaign tag.  That is what gives each campaign its own results
%% directory, which Run_One_Subfield requires of a scenario run because
%% Bouchet_Check_Existing_Result globs one pattern per run and errors on
%% more than one match.  TC_CAMPAIGN_TAG is still checked against the tag
%% derived here, so a submitter that disagrees with this script fails at the
%% gate rather than writing into the wrong tree.
%%
%% The chemistry is untouched.  Every calibrated number still comes from
%% Campaign_Parameters in the five-year 1-D campaign directory;
%% only the horizon is overridden, and only here.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[scenario,horizon,n_days,last_application_year] = read_scenario_and_horizon;
calcite_case = read_calcite_case;
campaign_tag = sprintf( ...
    'SSP_2D_%s_%s%s', ...
    horizon,scenario,case_suffix(calcite_case));
%% one_d_dir is the finished five-year campaign that supplies the calibrated
%% set; ssp_1d_dir is this campaign's own 1-D companion, which holds the
%% shared submitter.
one_d_tag = ...
    'Base_1D';
ssp_1d_tag = ...
    'SSP_1D';
campaign_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(campaign_dir));
one_d_dir = fullfile(repo_root,'runtime',one_d_tag);
ssp_1d_dir = fullfile(repo_root,'runtime',ssp_1d_tag);
run_dir = fullfile(repo_root,'v2','run');
code_dir = fullfile(repo_root,'v2','T&C_Code');
addpath(run_dir,code_dir,fullfile(code_dir,'Inputs'), ...
    fullfile(code_dir,'Parameters'),campaign_dir,one_d_dir);

matlab_release = version('-release');
if ~strcmp(matlab_release,'2025a')
    error('Bouchet_Kantola_SSP:WrongMatlabRelease', ...
        'This campaign requires MATLAB R2025a; loaded R%s.',matlab_release);
end
Bouchet_Single_Thread_Runtime('Bouchet_Kantola_SSP');
P = Campaign_Parameters;
%% The calibrated set was adopted on a five-year horizon.  This campaign runs
%% the same chemistry over a different span, so the horizon is the one field
%% that changes; StartDate and every chemical parameter stay as the sweep
%% adopted them.
P.NDays = n_days;

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

run_id = read_positive_integer('SLURM_ARRAY_TASK_ID');
if run_id > 9
    error('Bouchet_Kantola_SSP:BadRunId', ...
        'SLURM_ARRAY_TASK_ID must select run 1 through 9.');
end
allocated_cpus = read_positive_integer('SLURM_CPUS_PER_TASK');
%% One worker per column. The five-year campaign used six, so every worker
%% carried two of the twelve columns and the daily ERW parfor ran as two
%% sequential batches. Worker count does not change the numbers --
%% test_2d_parfor_equivalence asserts isequaln between serial and parfor
%% output across P, O, ERW and both pH records -- so this is a scheduling
%% choice, not a numerical one. It remains part of the declared identity,
%% which is why the fingerprint and the provenance assertion below receive
%% it rather than a constant.
workers = Resolve_2D_Parfor_Workers(12);
if workers ~= 12
    error('Bouchet_Kantola_SSP:BadWorkerCount', ...
        'TC_PARALLEL_WORKERS must be 12 for this campaign.');
end
if allocated_cpus < workers
    error('Bouchet_Kantola_SSP:InsufficientCpus', ...
        'Reserve at least %d CPUs.',workers);
end
%% Two hours. Checkpoints are a restart mechanism taking no part in the
%% computation, so the interval cannot move a result;
%% test_2d_parfor_equivalence runs with them disabled entirely. Each one
%% writes the whole frame state. In package u's 35-year runs, at thirty
%% minutes, a write took 2-8 minutes on average and grew with the run, to
%% about 15 minutes by half way through the whole-field basalt run, which
%% then spent a fifth of its time writing checkpoints. The cost of two
%% hours is that a job stopped part-way loses up to two hours.
if ~strcmp(getenv('TC_CKPT_MINUTES'),'120')
    error('Bouchet_Kantola_SSP:BadCheckpointInterval', ...
        'TC_CKPT_MINUTES must be 120.');
end
if ~isempty(getenv('TC_CKPT_FORCE_RESUME'))
    error('Bouchet_Kantola_SSP:ForcedResumeForbidden', ...
        'TC_CKPT_FORCE_RESUME must remain unset for production.');
end
if ~strcmp(strtrim(getenv('TC_CAMPAIGN_TAG')),campaign_tag)
    error('Bouchet_Kantola_SSP:BadCampaignTag', ...
        'TC_CAMPAIGN_TAG must be %s.',campaign_tag);
end
if ~license('test','Distrib_Computing_Toolbox')
    error('Bouchet_Kantola_SSP:ParallelToolboxUnavailable', ...
        'Parallel Computing Toolbox is required for this campaign.');
end
validate_only = strtrim(getenv('TC_VALIDATE_ONLY'));
if ~isempty(validate_only) && ~strcmp(validate_only,'1')
    error('Bouchet_Kantola_SSP:BadValidationMode', ...
        'TC_VALIDATE_ONLY must be unset or 1.');
end
%% The five-year 1-D campaign directory supplies the parameter set, the
%% provenance assertion and the digest helper, and its configuration test
%% gates the shared calibrated identity; they are part of this campaign's
%% gate.  This campaign's own submitter is included so that editing it
%% invalidates the stamp.
stamp_dependencies = { ...
    fullfile(one_d_dir,'Campaign_Parameters.m'), ...
    fullfile(one_d_dir,'Assert_Base_Provenance.m'), ...
    fullfile(one_d_dir,'Fingerprint_Digest.m'), ...
    fullfile(one_d_dir,'test_all_12_configuration.m'), ...
    fullfile(ssp_1d_dir,'submit_ssp.sh')};
if strcmp(validate_only,'1')
    run_validation(code_dir,repo_root,campaign_dir,campaign_tag, ...
        one_d_dir,stamp_dependencies,workers);
    return
end
%% Checked after validation: the validation element is array element 1.
if ~isempty(calcite_case) && ~ismember(run_id,[3 6 9])
    error('Bouchet_Kantola_SSP:BadRunId', ...
        'A calcite-case campaign runs the calcite treatments, runs 3, 6 and 9, only.');
end

results_dir = fullfile(repo_root,'v2','SIM','Subfield_12',campaign_tag);
logs_dir = fullfile(repo_root,'v2','logs','Subfield_12',campaign_tag);
stamp_file = fullfile(logs_dir,'validation_stamp.txt');
if ~isfile(stamp_file)
    error('Bouchet_Kantola_SSP:MissingValidation', ...
        'Run the R2025a validation job before production.');
end
expected_stamp = Bouchet_Campaign_Validation_Fingerprint( ...
    fullfile(repo_root,'v2'),campaign_dir,campaign_tag,12,workers, ...
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
%% so a requeued or duplicate array element fails here rather than racing.
%% The cleanup object must stay in scope for the whole run; it releases the
%% lock when this function returns or errors.
lock_cleanup = Bouchet_Acquire_Run_Lock( ...
    fullfile(logs_dir,'locks',sprintf('run_%02d.lock',run_id)));

if strcmp(Bouchet_Check_Existing_Result( ...
        results_dir,run_id,stamp_digest,12,P.NDays),'complete')
    fprintf('Run %d needs no work under the current validation.\n',run_id);
    return
end

worker_log = fullfile(logs_dir,sprintf('worker_run%02d.log',run_id));
diary(worker_log);
diary_cleanup = onCleanup(@() diary('off'));
fprintf(['Campaign: %s\nScenario: %s\nHorizon: %s (%d days)\n' ...
    'Slurm run: %d\nWorkers: %d\nAppelo alpha: %.17g\n' ...
    'KNmaxd: %.17g d^-1\nKbot: %.17g mm h^-1\nCEC: %g cmolc kg^-1\n'], ...
    campaign_tag,scenario,horizon,P.NDays,run_id,workers,P.AppeloAlpha, ...
    P.KNmaxd,P.Kbot,P.CEC_cmolc_kg);

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
    'GridNx',2,'GridNy',2, ...
    'AppeloAlpha',P.AppeloAlpha,'AppeloLogKH',P.AppeloLogKH, ...
    'KNmaxd',P.KNmaxd,'Kbot',P.Kbot, ...
    'ParallelWorkers',workers,'SaveResults',true, ...
    'ResultDirectory',results_dir,'Scenario',scenario,run_options{:});

status = fileread(status_file);
if ~contains(status,'state=complete')
    error('Bouchet_Kantola_SSP:IncompleteStatus', ...
        'Run %d did not write a complete status.',run_id);
end

%% Jobs execute from the shared source tree, which stays writable while
%% they run.  MATLAB loads functions lazily, so an upload part-way through
%% a multi-day run can change the chemistry mid-trajectory.  Recompute the
%% gate and reject the result if anything moved.
if ~strcmp(Bouchet_Campaign_Validation_Fingerprint( ...
        fullfile(repo_root,'v2'),campaign_dir,campaign_tag,12,workers, ...
        stamp_dependencies{:}),expected_stamp)
    error('Bouchet_Kantola_SSP:SourceChangedDuringRun', ...
        ['Code, inputs, tests or numerical runtime changed while run %d ' ...
         'was executing. Discard its result.'],run_id);
end
if ~strcmp(Bouchet_Check_Existing_Result( ...
        results_dir,run_id,stamp_digest,12,P.NDays),'complete')
    error('Bouchet_Kantola_SSP:MissingResult', ...
        'Run %d reported completion but wrote no matching final MAT.',run_id);
end

existing = dir(fullfile(results_dir,sprintf('TC_2D_run%02d_*.mat',run_id)));
existing = existing(~[existing.isdir]);
result_file = fullfile(existing(1).folder,existing(1).name);
saved = load(result_file,'run_provenance');
Assert_Base_Provenance(saved.run_provenance,run_id,12, ...
    workers,stamp_digest,campaign_tag,results_dir,P.NDays);
fprintf('Provenance: run %d carries the calibrated CEC 13 identity.\n',run_id);

%% The lateral scheme is proven conservative for a single step by
%% test_lateral_routing, but a campaign run takes many thousands of steps
%% across twelve columns and nothing checked that the residuals stayed small
%% once accumulated.  The driver records them; this reads them back.  It runs
%% after the result is saved, so a failure flags the run without losing it.
fprintf('Lateral conservation over the whole run:\n');
Check_Lateral_Conservation(result_file,'Bouchet_Kantola_SSP');
fprintf('Run %d complete under %s.\n',run_id,campaign_tag);
end

function run_validation(code_dir,repo_root,campaign_dir,campaign_tag, ...
    one_d_dir,stamp_dependencies,workers)
test_dir = fullfile(code_dir,'Tests');
addpath(test_dir,one_d_dir);
old_dir = pwd;
cleanup = onCleanup(@() cd(old_dir));
cd(test_dir);
%% Production also depends on the complete 1-D audit, the pre-application
%% identity check and the one- and twelve-column smoke, which the 1-D
%% validation runs.  This element checks the 2-D wrapper and the Linux
%% deployment without repeating those on a second node.
test_bouchet_2d_parallel_configuration;
test_biogeochemistry_io_path_portability;
test_bouchet_linux_portability({campaign_dir,one_d_dir});
test_all_12_configuration;
logs_dir = fullfile(repo_root,'v2','logs','Subfield_12',campaign_tag);
fingerprint = Bouchet_Campaign_Validation_Fingerprint( ...
    fullfile(repo_root,'v2'),campaign_dir,campaign_tag,12,workers, ...
    stamp_dependencies{:});
Bouchet_Write_Validation_Stamp( ...
    fullfile(logs_dir,'validation_stamp.txt'),fingerprint);
fprintf('PASS: Bouchet R2025a CEC 13 SSP parallel 2-D validation (%s)\n', ...
    campaign_tag);
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

function value = read_positive_integer(name)
raw = strtrim(getenv(name));
value = str2double(raw);
if isempty(raw) || ~(isfinite(value) && value == fix(value) && value >= 1)
    error('Bouchet_Kantola_SSP:BadEnvironment', ...
        '%s must be a positive integer.',name);
end
end
