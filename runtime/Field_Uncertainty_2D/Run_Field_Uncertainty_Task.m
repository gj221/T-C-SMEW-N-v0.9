function Run_Field_Uncertainty_Task(array_index)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% One array element of the Kantola uncertainty ensemble on a field grid:
%% the lower or upper parameter set (Uncertainty_Members, Uncertainty_Knobs)
%% on one run of the central campaign on that grid, over its horizon (the
%% 35-year projection, or the five-year field record), for the campaign
%% the controller's environment names (Field_Uncertainty_Campaign). The
%% original set is the central field campaign's result (runtime/Field_SSP_2D).
%%
%% The gates are Field_2D's, as runtime/Field_SSP_2D applies them: a
%% validation element (TC_VALIDATE_ONLY=1, array index 1) runs the complete
%% audit, the Linux checks and one day of the grid under each set and
%% writes the stamp; production checks the stamp, takes the run lock, runs
%% with TC_FIELD_WORKERS workers and two-hourly checkpoints, recomputes the
%% gate, reads the set's provenance and the grid's back, and checks the
%% lateral conservation. The set's knobs reach every square.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

prefix = 'Kantola_FieldUncertainty';
here = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(here));
code_dir = fullfile(repo_root,'v2','T&C_Code');
addpath(here,fullfile(repo_root,'runtime','Uncertainty_1D'), ...
    fullfile(repo_root,'runtime','Base_1D'),fullfile(repo_root,'runtime','Field_SSP_2D'), ...
    fullfile(repo_root,'v2','run'),code_dir,fullfile(code_dir,'Inputs'), ...
    fullfile(code_dir,'Parameters'));

matlab_release = version('-release');
if ~strcmp(matlab_release,'2025a')
    error([prefix ':WrongMatlabRelease'], ...
        'This campaign requires MATLAB R2025a; loaded R%s.',matlab_release);
end
Bouchet_Single_Thread_Runtime(prefix);
C = Field_Uncertainty_Campaign;
P = Campaign_Parameters;
if ~(isnumeric(array_index) && isscalar(array_index) && ...
        array_index == fix(array_index) && array_index >= 1)
    error([prefix ':BadIndex'],'The array index must be a positive integer.');
end
workers = Resolve_2D_Parfor_Workers(C.n_columns);
if workers ~= C.workers
    error([prefix ':BadWorkerCount'], ...
        'TC_PARALLEL_WORKERS must equal TC_FIELD_WORKERS (%d) for this campaign.',C.workers);
end
if str2double(getenv('SLURM_CPUS_PER_TASK')) < workers
    error([prefix ':InsufficientCpus'],'Reserve at least %d CPUs.',workers);
end
%% Package z's fast ERW chemistry kernels (the original code's numbers, bit
%% for bit), the only ones a field run runs; each result records them.
kernel = Field_Campaign_Grid('kernel',prefix);
if ~strcmp(getenv('TC_CKPT_MINUTES'),'120')
    error([prefix ':BadCheckpointInterval'],'TC_CKPT_MINUTES must be 120.');
end
if ~isempty(getenv('TC_CKPT_FORCE_RESUME'))
    error([prefix ':ForcedResumeForbidden'], ...
        'TC_CKPT_FORCE_RESUME must remain unset for production.');
end
if ~strcmp(strtrim(getenv('TC_CAMPAIGN_TAG')),C.tag)
    error([prefix ':BadCampaignTag'],'TC_CAMPAIGN_TAG must be %s.',C.tag);
end
if ~license('test','Distrib_Computing_Toolbox')
    error([prefix ':ParallelToolboxUnavailable'], ...
        'Parallel Computing Toolbox is required for the 2-D runs.');
end
required_files = { ...
    fullfile(code_dir,'Inputs',C.forcing), ...
    fullfile(repo_root,'v2','SIM','EQUIV_1D_TMP','TC_ERW_results_Kantola_crop_control.mat')};
for k = 1:numel(required_files)
    if ~isfile(required_files{k})
        error([prefix ':MissingDeploymentFile'], ...
            'Missing required deployment file: %s',required_files{k});
    end
end
validate_only = strtrim(getenv('TC_VALIDATE_ONLY'));
if ~isempty(validate_only) && ~strcmp(validate_only,'1')
    error([prefix ':BadValidationMode'],'TC_VALIDATE_ONLY must be unset or 1.');
end
logs_root = fullfile(repo_root,'v2','logs','Subfield_12',C.tag);
stamp_file = fullfile(logs_root,'validation_stamp.txt');
if strcmp(validate_only,'1')
    run_validation(code_dir,repo_root,here,C,P,workers,stamp_file,kernel,prefix);
    return
end
if array_index > numel(C.elements)
    fprintf('Element %d is beyond the %d elements of %s; nothing to do.\n', ...
        array_index,numel(C.elements),C.tag);
    return
end
member = C.elements(array_index).member;
run_id = C.elements(array_index).run_id;
K = Uncertainty_Knobs(member,P);
Field_Campaign_Grid('check_memory',C.n_columns,C.n_days,prefix);

if ~isfile(stamp_file)
    error([prefix ':MissingValidation'],'Run the R2025a validation job before production.');
end
fingerprint = Bouchet_Campaign_Validation_Fingerprint( ...
    fullfile(repo_root,'v2'),here,C.tag,C.n_columns,workers,C.dependencies{:});
if ~strcmp(strtrim(fileread(stamp_file)),fingerprint)
    error([prefix ':StaleValidation'], ...
        'Code, inputs, grid, tests or numerical runtime changed after validation.');
end
digest = Fingerprint_Digest(fingerprint);
setenv('TC_VALIDATION_FINGERPRINT',digest);
results_dir = fullfile(repo_root,'v2','SIM','Subfield_12',C.tag,member.tag);
logs_dir = fullfile(logs_root,member.tag);

%% A duplicate element fails here rather than sharing the status file,
%% checkpoint and result path; a requeued one keeps its job id and
%% reclaims its own lock.
lock_cleanup = Bouchet_Acquire_Run_Lock( ...
    fullfile(logs_dir,'locks',sprintf('run_%02d.lock',run_id))); %#ok<NASGU>

%% A result already there must be this set's, or it is in the way.
[state,result_file] = Uncertainty_Result_State(results_dir,run_id,C,K,C.tag,digest);
switch state
    case 'complete'
        Field_Campaign_Grid('provenance',load_provenance(result_file), ...
            C.topology_file,C.topo,prefix);
        Check_Lateral_Conservation(result_file,prefix);
        fprintf('Set %s run %d (%s) is already complete with matching provenance; nothing to do.\n', ...
            member.tag,run_id,C.tag);
        return
    case 'stale'
        error([prefix ':StaleResult'], ...
            '%s does not carry the provenance of set %s; remove it before resubmitting.', ...
            result_file,member.tag);
end

if ~exist(results_dir,'dir'); mkdir(results_dir); end
worker_log = fullfile(logs_dir,sprintf('worker_run%02d.log',run_id));
diary(worker_log);
diary_cleanup = onCleanup(@() diary('off')); %#ok<NASGU>
fprintf(['Campaign: %s (original: %s)\nSet: %s\nRun: %d\nGrid: %s (%d squares)\n' ...
    'Horizon: %s (%d days, rock to %d)\nWorkers: %d\nERW kernel: %s\n' ...
    'Slurm restart: %s\nStarted: %s\n'], ...
    C.tag,C.original_tag,member.tag,run_id,C.topology_file,C.n_columns,C.horizon, ...
    C.n_days,C.last_application_year,workers,kernel,getenv('SLURM_RESTART_COUNT'), ...
    datestr(now,'yyyy-mm-dd HH:MM:SS'));

status_file = fullfile(logs_dir,sprintf('run_%02d.status',run_id));
options = [{'NDays',C.n_days,'StartDate',P.StartDate}, K.options, C.options, ...
    {'ParallelWorkers',workers,'SaveResults',true,'ResultDirectory',results_dir}];
Run_Subfield_Worker(run_id,status_file,options{:});
if ~contains(fileread(status_file),'state=complete')
    error([prefix ':IncompleteStatus'],'Run %d did not write a complete status.',run_id);
end
%% MATLAB loads functions lazily, so an upload part-way through a
%% multi-day run could change the chemistry mid-trajectory.
if ~strcmp(Bouchet_Campaign_Validation_Fingerprint( ...
        fullfile(repo_root,'v2'),here,C.tag,C.n_columns,workers,C.dependencies{:}),fingerprint)
    error([prefix ':SourceChangedDuringRun'], ...
        'Code, inputs or numerical runtime changed while run %d was executing. Discard its result.',run_id);
end

%% Read the provenance back: the set and the grid must be what was run.
[state,result_file] = Uncertainty_Result_State(results_dir,run_id,C,K,C.tag,digest);
if ~strcmp(state,'complete')
    error([prefix ':Provenance'], ...
        'The saved result of run %d in %s does not carry the provenance of set %s.', ...
        run_id,results_dir,member.tag);
end
p = load_provenance(result_file);
Field_Campaign_Grid('provenance',p,C.topology_file,C.topo,prefix);
Field_Campaign_Grid('kernel_provenance',p,kernel,prefix);
%% After the save, so a leak flags the run without losing it.
fprintf('Lateral conservation over the whole run:\n');
Check_Lateral_Conservation(result_file,prefix);
fprintf('PASS: set %s run %d (%s) saved with matching provenance: %s\n', ...
    member.tag,run_id,C.tag,result_file);
end

function run_validation(code_dir,repo_root,here,C,P,workers,stamp_file,kernel,prefix)
test_dir = fullfile(code_dir,'Tests');
addpath(test_dir);
old_dir = pwd;
cleanup = onCleanup(@() cd(old_dir));
cd(test_dir);
%% The complete audit, the Linux checks of this campaign and the code it
%% borrows, then one day of the grid across the first application under
%% each set, saved and read back as production does: the long knob names
%% of a field grid go through Bounded_File_Tag here.
run_post_legacy_equation_audit_tests(1);
cd(test_dir);
test_biogeochemistry_io_path_portability;
test_bouchet_linux_portability({here,fullfile(repo_root,'runtime','Uncertainty_1D'), ...
    fullfile(repo_root,'runtime','Field_SSP_2D')});
cd(old_dir);
smoke_run = 5;
if ~isempty(C.calcite_case); smoke_run = 6; end
members = Uncertainty_Members;
members = members(~[members.reused]);
C_smoke = C;
C_smoke.n_days = 1;
for m = 1:numel(members)
    K = Uncertainty_Knobs(members(m),P);
    smoke_dir = tempname;
    mkdir(smoke_dir);
    smoke_cleanup = onCleanup(@() rmdir(smoke_dir,'s'));
    ckpt_minutes = getenv('TC_CKPT_MINUTES');
    setenv('TC_CKPT_MINUTES','0');
    options = [{'NDays',1,'StartDate','2016-11-10'}, K.options, C.options, ...
        {'ParallelWorkers',workers,'SaveResults',true,'ResultDirectory',smoke_dir}];
    smoke_file = Run_One_Subfield(smoke_run,options{:});
    setenv('TC_CKPT_MINUTES',ckpt_minutes);
    [~,name,extension] = fileparts(smoke_file);
    p = load_provenance(smoke_file);
    if ~Uncertainty_Result_Matches(p,[name extension],C_smoke,K,C.tag,'')
        error([prefix ':SmokeProvenance'], ...
            'The one-day %s run does not carry the set''s provenance: %s',members(m).tag,smoke_file);
    end
    Field_Campaign_Grid('provenance',p,C.topology_file,C.topo,prefix);
    Field_Campaign_Grid('kernel_provenance',p,kernel,prefix);
    fprintf('Smoke %s (%d characters, ERW kernel %s): %s%s\n',members(m).tag, ...
        numel([name extension]),kernel,name,extension);
    Check_Lateral_Conservation(smoke_file,prefix);
    clear smoke_cleanup
end
fprintf('Memory floor for a production run: %.1f GiB (%d squares, %d days).\n', ...
    Field_Campaign_Grid('memory',C.n_columns,C.n_days),C.n_columns,C.n_days);
fingerprint = Bouchet_Campaign_Validation_Fingerprint( ...
    fullfile(repo_root,'v2'),here,C.tag,C.n_columns,workers,C.dependencies{:});
Bouchet_Write_Validation_Stamp(stamp_file,fingerprint);
fprintf('PASS: Bouchet R2025a field-grid uncertainty validation (%s)\n',C.tag);
end

function p = load_provenance(result_file)
S = load(result_file,'run_provenance');
p = S.run_provenance;
end
