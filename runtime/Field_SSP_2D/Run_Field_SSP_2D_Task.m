function Run_Field_SSP_2D_Task
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run one Slurm array element of a 35-year Kantola SSP campaign on a
%% field grid: run 1 to 9 of Run_Matrix_12 on T&C-SMEW-N v1.0's polygon
%% grid of the real field (runtime/Field_2D) under the scenario forcing,
%% rock every 10 November to 2050, as runtime/SSP_2D runs them on the
%% 12-column plot.
%%
%% The calibrated set is Campaign_Parameters with the horizon moved to
%% 35 years, and the run options are SSP_2D's with the field grid in place
%% of the 2x2 plot grid, so the result is the original member of the
%% field-grid uncertainty ensemble (runtime/Field_Uncertainty_2D).
%%
%% The gates are Field_2D's: the grid file beside the deployed root, its
%% content in the validation stamp, TC_FIELD_WORKERS workers for the
%% squares' daily, hourly and ERW work (package y), two-hourly
%% checkpoints, a run lock, the stamp recomputed after the run, the
%% calibrated provenance, the grid's provenance and the lateral
%% conservation read back. The job runs on the day partition and requeues
%% itself before the 24-hour limit (requeue_near_limit.sh); a requeued
%% element continues from its checkpoint under the same lock.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

prefix = 'Bouchet_Kantola_FieldSSP';
one_d_tag = ...
    'Base_1D';
campaign_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(campaign_dir));
one_d_dir = fullfile(repo_root,'runtime',one_d_tag);
run_dir = fullfile(repo_root,'v2','run');
code_dir = fullfile(repo_root,'v2','T&C_Code');
addpath(run_dir,code_dir,fullfile(code_dir,'Inputs'), ...
    fullfile(code_dir,'Parameters'),campaign_dir,one_d_dir);

matlab_release = version('-release');
if ~strcmp(matlab_release,'2025a')
    error([prefix ':WrongMatlabRelease'], ...
        'This campaign requires MATLAB R2025a; loaded R%s.',matlab_release);
end
Bouchet_Single_Thread_Runtime(prefix);
P = Campaign_Parameters;
%% The 35-year horizon of the SSP campaigns: 2016-01-01 to 2050-12-31
%% (12784 days), rock every 10 November to the end of the forcing.
n_days = 12784;
last_application_year = 2050;
P.NDays = n_days;

scenario = read_scenario(prefix);
[calcite_case,case_suffix] = Calcite_Case;
[topology_file,topology_name] = Field_Campaign_Grid('read',repo_root,prefix);
campaign_tag = sprintf('Field_SSP_2D_35yr_%s%s__%s',scenario,case_suffix,topology_name);
topo = Field_Campaign_Grid('load',topology_file,prefix);
n_columns = topo.N_cells;

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
    fullfile(one_d_dir,'Calcite_Case.m'), ...
    fullfile(one_d_dir,'Assert_Base_Provenance.m'), ...
    fullfile(one_d_dir,'Fingerprint_Digest.m'), ...
    topology_file};
for k = 1:numel(required_files)
    if ~isfile(required_files{k})
        error([prefix ':MissingDeploymentFile'], ...
            'Missing required deployment file: %s',required_files{k});
    end
end

run_id = read_positive_integer('SLURM_ARRAY_TASK_ID',prefix);
if run_id > 9
    error([prefix ':BadRunId'],'SLURM_ARRAY_TASK_ID must select run 1 through 9.');
end
allocated_cpus = read_positive_integer('SLURM_CPUS_PER_TASK',prefix);
%% TC_FIELD_WORKERS workers share out the squares' daily vegetation and
%% biogeochemistry, hourly physics and daily ERW chemistry (package y). The
%% count moves no number (test_2d_parfor_equivalence) but is part of the
%% declared identity, the stamp and the provenance.
field_workers = read_positive_integer('TC_FIELD_WORKERS',prefix);
workers = Resolve_2D_Parfor_Workers(n_columns);
if workers ~= min(field_workers,n_columns)
    error([prefix ':BadWorkerCount'], ...
        'TC_PARALLEL_WORKERS must equal TC_FIELD_WORKERS (%d) for this campaign.', ...
        field_workers);
end
if allocated_cpus < workers
    error([prefix ':InsufficientCpus'],'Reserve at least %d CPUs.',workers);
end
%% Package z's ERW chemistry kernels: matlab (the original code) or fast
%% (the same numbers, bit for bit); each result records the one it ran.
kernel = Field_Campaign_Grid('kernel',prefix);
%% Two hours, as the SSP and Field_2D campaigns: each write holds the
%% daily stores and state, so a shorter interval costs more than the work
%% a requeue repeats.
if ~strcmp(getenv('TC_CKPT_MINUTES'),'120')
    error([prefix ':BadCheckpointInterval'],'TC_CKPT_MINUTES must be 120.');
end
if ~isempty(getenv('TC_CKPT_FORCE_RESUME'))
    error([prefix ':ForcedResumeForbidden'], ...
        'TC_CKPT_FORCE_RESUME must remain unset for production.');
end
if ~strcmp(strtrim(getenv('TC_CAMPAIGN_TAG')),campaign_tag)
    error([prefix ':BadCampaignTag'],'TC_CAMPAIGN_TAG must be %s.',campaign_tag);
end
if ~license('test','Distrib_Computing_Toolbox')
    error([prefix ':ParallelToolboxUnavailable'], ...
        'Parallel Computing Toolbox is required for this campaign.');
end
validate_only = strtrim(getenv('TC_VALIDATE_ONLY'));
if ~isempty(validate_only) && ~strcmp(validate_only,'1')
    error([prefix ':BadValidationMode'],'TC_VALIDATE_ONLY must be unset or 1.');
end
%% The 1-D campaign directory supplies the parameter set, the calcite case,
%% the provenance assertion and the digest helper; the scenario forcing
%% and the grid file are this campaign's own inputs. All are in the gate.
stamp_dependencies = { ...
    fullfile(one_d_dir,'Campaign_Parameters.m'), ...
    fullfile(one_d_dir,'Calcite_Case.m'), ...
    fullfile(one_d_dir,'Assert_Base_Provenance.m'), ...
    fullfile(one_d_dir,'Fingerprint_Digest.m'), ...
    fullfile(repo_root,'runtime','Uncertainty_1D','requeue_near_limit.sh'), ...
    fullfile(code_dir,'Inputs',sprintf('Data_Kantola_%s.mat',scenario)), ...
    topology_file};
run_options = {'NDays',n_days,'StartDate',P.StartDate, ...
    'FieldTopologyFile',topology_file, ...
    'AppeloAlpha',P.AppeloAlpha,'AppeloLogKH',P.AppeloLogKH, ...
    'KNmaxd',P.KNmaxd,'Kbot',P.Kbot, ...
    'Scenario',scenario,'LastApplicationYear',last_application_year};
if ~isempty(calcite_case)
    run_options = [run_options {'CalciteCase',calcite_case}];
end
if strcmp(validate_only,'1')
    run_validation(code_dir,repo_root,campaign_dir,campaign_tag,one_d_dir, ...
        stamp_dependencies,workers,topology_file,topo,run_options,calcite_case, ...
        n_days,kernel,prefix);
    return
end
%% Checked after validation: the validation element is array element 1.
if ~isempty(calcite_case) && ~ismember(run_id,[3 6 9])
    error([prefix ':BadRunId'], ...
        'A calcite-case campaign runs the calcite treatments, runs 3, 6 and 9, only.');
end
Field_Campaign_Grid('check_memory',n_columns,n_days,prefix);

results_dir = fullfile(repo_root,'v2','SIM','Subfield_12',campaign_tag);
logs_dir = fullfile(repo_root,'v2','logs','Subfield_12',campaign_tag);
stamp_file = fullfile(logs_dir,'validation_stamp.txt');
if ~isfile(stamp_file)
    error([prefix ':MissingValidation'],'Run the R2025a validation job before production.');
end
expected_stamp = Bouchet_Campaign_Validation_Fingerprint( ...
    fullfile(repo_root,'v2'),campaign_dir,campaign_tag,n_columns,workers, ...
    stamp_dependencies{:});
if ~strcmp(strtrim(fileread(stamp_file)),expected_stamp)
    error([prefix ':StaleValidation'], ...
        'Code, inputs, grid, tests or numerical runtime changed after validation.');
end
stamp_digest = Fingerprint_Digest(expected_stamp);
setenv('TC_VALIDATION_FINGERPRINT',stamp_digest);
if ~exist(results_dir,'dir'); mkdir(results_dir); end
if ~exist(logs_dir,'dir'); mkdir(logs_dir); end

%% A requeued element keeps its job id and reclaims its own lock; any other
%% holder is a live duplicate or a killed job that needs a human.
lock_cleanup = Bouchet_Acquire_Run_Lock( ...
    fullfile(logs_dir,'locks',sprintf('run_%02d.lock',run_id))); %#ok<NASGU>

if strcmp(Bouchet_Check_Existing_Result( ...
        results_dir,run_id,stamp_digest,n_columns,n_days),'complete')
    fprintf('Run %d needs no work under the current validation.\n',run_id);
    return
end

worker_log = fullfile(logs_dir,sprintf('worker_run%02d.log',run_id));
diary(worker_log);
diary_cleanup = onCleanup(@() diary('off')); %#ok<NASGU>
fprintf(['Campaign: %s\nScenario: %s (%d days, rock to %d)\n' ...
    'Grid: %s (%d squares of %g m, topography %d)\nSlurm run: %d (restart %s)\n' ...
    'Workers: %d\nERW kernel: %s\nAppelo alpha: %.17g\nKNmaxd: %.17g d^-1\n' ...
    'Kbot: %.17g mm h^-1\nCEC: %g cmolc kg^-1\nStarted: %s\n'], ...
    campaign_tag,scenario,n_days,last_application_year,topology_file,n_columns, ...
    topo.cell_size_m,isfield(topo,'topography_attached') && topo.topography_attached, ...
    run_id,getenv('SLURM_RESTART_COUNT'),workers,kernel,P.AppeloAlpha,P.KNmaxd,P.Kbot, ...
    P.CEC_cmolc_kg,datestr(now,'yyyy-mm-dd HH:MM:SS'));

status_file = fullfile(logs_dir,sprintf('run_%02d.status',run_id));
Run_Subfield_Worker(run_id,status_file,run_options{:}, ...
    'ParallelWorkers',workers,'SaveResults',true,'ResultDirectory',results_dir);

status = fileread(status_file);
if ~contains(status,'state=complete')
    error([prefix ':IncompleteStatus'],'Run %d did not write a complete status.',run_id);
end
%% Jobs execute from the shared source tree, which stays writable while
%% they run. Recompute the gate and reject the result if anything moved.
if ~strcmp(Bouchet_Campaign_Validation_Fingerprint( ...
        fullfile(repo_root,'v2'),campaign_dir,campaign_tag,n_columns,workers, ...
        stamp_dependencies{:}),expected_stamp)
    error([prefix ':SourceChangedDuringRun'], ...
        ['Code, inputs, grid, tests or numerical runtime changed while run %d ' ...
         'was executing. Discard its result.'],run_id);
end
if ~strcmp(Bouchet_Check_Existing_Result( ...
        results_dir,run_id,stamp_digest,n_columns,n_days),'complete')
    error([prefix ':MissingResult'], ...
        'Run %d reported completion but wrote no matching final MAT.',run_id);
end

existing = dir(fullfile(results_dir,sprintf('TC_2D_run%02d_*.mat',run_id)));
existing = existing(~[existing.isdir]);
result_file = fullfile(existing(1).folder,existing(1).name);
saved = load(result_file,'run_provenance');
Assert_Base_Provenance(saved.run_provenance,run_id,n_columns,workers, ...
    stamp_digest,campaign_tag,results_dir,n_days);
Field_Campaign_Grid('provenance',saved.run_provenance,topology_file,topo,prefix);
Field_Campaign_Grid('kernel_provenance',saved.run_provenance,kernel,prefix);
assert_projection_provenance(saved.run_provenance,existing(1).name,scenario, ...
    calcite_case,last_application_year,prefix);
fprintf('Provenance: run %d carries the calibrated CEC 13 identity on %s under %s.\n', ...
    run_id,topology_name,scenario);

fprintf('Lateral conservation over the whole run:\n');
Check_Lateral_Conservation(result_file,prefix);
fprintf('Run %d complete under %s.\n',run_id,campaign_tag);
end

function run_validation(code_dir,repo_root,campaign_dir,campaign_tag,one_d_dir, ...
    stamp_dependencies,workers,topology_file,topo,run_options,calcite_case, ...
    n_days,kernel,prefix)
test_dir = fullfile(code_dir,'Tests');
addpath(test_dir,one_d_dir);
old_dir = pwd;
cleanup = onCleanup(@() cd(old_dir));
cd(test_dir);
%% The complete audit (the post-legacy suite, with the field-grid and
%% field-campaign tests), the Linux deployment, then the grid itself for
%% one day across the first application under this campaign's scenario,
%% case and calibrated options, saved and read back as production does.
run_post_legacy_equation_audit_tests(1);
cd(test_dir);
test_biogeochemistry_io_path_portability;
test_bouchet_linux_portability({campaign_dir,one_d_dir});
cd(old_dir);
smoke_dir = tempname;
mkdir(smoke_dir);
smoke_cleanup = onCleanup(@() rmdir(smoke_dir,'s'));
smoke_run = 5;
if ~isempty(calcite_case); smoke_run = 6; end
smoke_options = run_options;
smoke_options(1:4) = {'NDays',1,'StartDate','2016-11-10'};
ckpt_minutes = getenv('TC_CKPT_MINUTES');
setenv('TC_CKPT_MINUTES','0');
smoke_file = Run_One_Subfield(smoke_run,smoke_options{:}, ...
    'ParallelWorkers',workers,'SaveResults',true,'ResultDirectory',smoke_dir);
setenv('TC_CKPT_MINUTES',ckpt_minutes);
[~,smoke_name,smoke_extension] = fileparts(smoke_file);
saved = load(smoke_file,'run_provenance');
Field_Campaign_Grid('provenance',saved.run_provenance,topology_file,topo,prefix);
Field_Campaign_Grid('kernel_provenance',saved.run_provenance,kernel,prefix);
fprintf('Smoke result (%d characters, ERW kernel %s): %s%s\n', ...
    numel([smoke_name smoke_extension]),kernel,smoke_name,smoke_extension);
Check_Lateral_Conservation(smoke_file,prefix);
clear smoke_cleanup
fprintf('Memory floor for a production run: %.1f GiB (%d squares, %d days).\n', ...
    Field_Campaign_Grid('memory',topo.N_cells,n_days),topo.N_cells,n_days);
logs_dir = fullfile(repo_root,'v2','logs','Subfield_12',campaign_tag);
fingerprint = Bouchet_Campaign_Validation_Fingerprint( ...
    fullfile(repo_root,'v2'),campaign_dir,campaign_tag,topo.N_cells,workers, ...
    stamp_dependencies{:});
Bouchet_Write_Validation_Stamp(fullfile(logs_dir,'validation_stamp.txt'),fingerprint);
fprintf('PASS: Bouchet R2025a CEC 13 field-grid 35-year validation (%s)\n',campaign_tag);
end

function assert_projection_provenance(p,file_name,scenario,calcite_case, ...
    last_application_year,prefix)
%% The scenario and calcite case are in the file name alone; the rock
%% schedule's last year is in the provenance.
case_tag = '';
if ~isempty(calcite_case); case_tag = sprintf('__case%d',calcite_case); end
ok = isfield(p,'last_application_year') && ...
    p.last_application_year == last_application_year && ...
    strcmp(regexp(file_name,'__ssp[0-9]{3}(?=__)','match','once'),['__' scenario]) && ...
    strcmp(regexp(file_name,'__case[0-9]+(?=__)','match','once'),case_tag);
if ~ok
    error([prefix ':WrongProjection'], ...
        '%s is not a %s%s run with rock to %d.',file_name,scenario,case_tag, ...
        last_application_year);
end
end

function scenario = read_scenario(prefix)
scenario = strtrim(getenv('TC_SSP'));
if ~ismember(scenario,{'ssp126','ssp245','ssp585'})
    error([prefix ':BadScenario'], ...
        ['TC_SSP must be ssp126, ssp245 or ssp585, not ''%s''. It selects the ' ...
        'forcing file Data_Kantola_<scenario>.mat.'],scenario);
end
end

function value = read_positive_integer(name,prefix)
raw = strtrim(getenv(name));
value = str2double(raw);
if isempty(raw) || ~(isfinite(value) && value == fix(value) && value >= 1)
    error([prefix ':BadEnvironment'],'%s must be a positive integer.',name);
end
end
