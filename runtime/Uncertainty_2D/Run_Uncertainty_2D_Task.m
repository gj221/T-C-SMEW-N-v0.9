function Run_Uncertainty_2D_Task(array_index)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% One array element of the Kantola 2-D uncertainty ensemble: the lower or
%% upper parameter set on one of runs 1-9 (3, 6 and 9 for a calcite case)
%% on the 2x2 grid, for the campaign the controller's environment names
%% (Uncertainty_Campaign). The original set is not run here; it is the
%% Base_2D or SSP_2D campaign's own result.
%%
%% The path is the SSP 2-D campaign's: the worker count of the original's
%% campaign, a run lock taken before anything touches the result, the
%% gate recomputed after the run, and the accumulated lateral residuals
%% read back. The set's knobs reach every column, buffer columns included,
%% and the saved provenance is checked against the set per column.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

here = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(here));
addpath(here,fullfile(repo_root,'runtime','Uncertainty_1D'), ...
    fullfile(repo_root,'runtime','Base_1D'),fullfile(repo_root,'runtime','Base_2D'), ...
    fullfile(repo_root,'v2','run'),fullfile(repo_root,'v2','T&C_Code'));

C = Uncertainty_Campaign('2D');
if ~(isnumeric(array_index) && isscalar(array_index) && ...
        array_index == fix(array_index) && array_index >= 1)
    error('Kantola_Uncertainty:BadIndex','The array index must be a positive integer.');
end
if array_index > numel(C.elements)
    fprintf('Element %d is beyond the %d elements of %s; nothing to do.\n', ...
        array_index,numel(C.elements),C.tag);
    return
end
member = C.elements(array_index).member;
run_id = C.elements(array_index).run_id;
P = Campaign_Parameters;
K = Uncertainty_Knobs(member,P);

Bouchet_Single_Thread_Runtime('Kantola_Uncertainty');
workers = Resolve_2D_Parfor_Workers(C.n_columns);
if workers ~= C.workers
    error('Kantola_Uncertainty:BadWorkerCount', ...
        'TC_PARALLEL_WORKERS must give %d workers, as the original''s campaign ran.',C.workers);
end
if str2double(getenv('SLURM_CPUS_PER_TASK')) < workers
    error('Kantola_Uncertainty:InsufficientCpus','Reserve at least %d CPUs.',workers);
end
if ~isempty(getenv('TC_CKPT_FORCE_RESUME'))
    error('Kantola_Uncertainty:ForcedResumeForbidden', ...
        'TC_CKPT_FORCE_RESUME must remain unset for production.');
end
if ~license('test','Distrib_Computing_Toolbox')
    error('Kantola_Uncertainty:ParallelToolboxUnavailable', ...
        'Parallel Computing Toolbox is required for the 2-D runs.');
end
results_dir = fullfile(repo_root,'v2','SIM','Subfield_12',C.tag,member.tag);
logs_dir = fullfile(repo_root,'v2','logs','Subfield_12',C.tag,member.tag);
fingerprint = Bouchet_Campaign_Validation_Fingerprint( ...
    fullfile(repo_root,'v2'),here,C.tag,C.n_columns,workers,C.dependencies{:});
digest = Fingerprint_Digest(fingerprint);
setenv('TC_VALIDATION_FINGERPRINT',digest);
setenv('TC_CAMPAIGN_TAG',C.tag);

%% A duplicate or requeued element fails here rather than sharing the
%% status file, checkpoint and result path. The cleanup object releases
%% the lock when this function returns or errors.
lock_cleanup = Bouchet_Acquire_Run_Lock( ...
    fullfile(logs_dir,'locks',sprintf('run_%02d.lock',run_id))); %#ok<NASGU>

%% A result already there must be this set's, or it is in the way.
[state,result_file] = Uncertainty_Result_State(results_dir,run_id,C,K,C.tag,digest);
switch state
    case 'complete'
        Check_Lateral_Conservation(result_file,'Kantola_Uncertainty');
        fprintf('Set %s run %d (%s) is already complete with matching provenance; nothing to do.\n', ...
            member.tag,run_id,C.tag);
        return
    case 'stale'
        error('Kantola_Uncertainty:StaleResult', ...
            '%s does not carry the provenance of set %s; remove it before resubmitting.', ...
            result_file,member.tag);
end

if ~exist(results_dir,'dir'); mkdir(results_dir); end
worker_log = fullfile(logs_dir,sprintf('worker_run%02d.log',run_id));
diary(worker_log);
diary_cleanup = onCleanup(@() diary('off')); %#ok<NASGU>
fprintf('Campaign: %s (original: %s)\nSet: %s\nRun: %d\nHorizon: %s (%d days)\nWorkers: %d\n', ...
    C.tag,C.original_tag,member.tag,run_id,C.horizon,C.n_days,workers);

status_file = fullfile(logs_dir,sprintf('run_%02d.status',run_id));
options = [{'NDays',C.n_days,'StartDate',P.StartDate}, K.options, C.options, ...
    {'ParallelWorkers',workers,'SaveResults',true,'ResultDirectory',results_dir}];
Run_Subfield_Worker(run_id,status_file,options{:});
if ~contains(fileread(status_file),'state=complete')
    error('Kantola_Uncertainty:IncompleteStatus', ...
        'Run %d did not write a complete status.',run_id);
end
%% MATLAB loads functions lazily, so an upload part-way through a
%% multi-day run could change the chemistry mid-trajectory.
if ~strcmp(Bouchet_Campaign_Validation_Fingerprint( ...
        fullfile(repo_root,'v2'),here,C.tag,C.n_columns,workers,C.dependencies{:}),fingerprint)
    error('Kantola_Uncertainty:SourceChangedDuringRun', ...
        'Code, inputs or numerical runtime changed while run %d was executing. Discard its result.',run_id);
end

%% Read the provenance back: the set must be what was run.
[state,result_file] = Uncertainty_Result_State(results_dir,run_id,C,K,C.tag,digest);
if ~strcmp(state,'complete')
    error('Kantola_Uncertainty:Provenance', ...
        'The saved result of run %d in %s does not carry the provenance of set %s.', ...
        run_id,results_dir,member.tag);
end
%% After the save, so a leak flags the run without losing it.
fprintf('Lateral conservation over the whole run:\n');
Check_Lateral_Conservation(result_file,'Kantola_Uncertainty');
fprintf('PASS: set %s run %d (%s) saved with matching provenance: %s\n', ...
    member.tag,run_id,C.tag,result_file);
end
