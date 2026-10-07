function Run_Uncertainty_Task(array_index)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% One array element of the Kantola 1-D uncertainty ensemble: the lower or
%% upper parameter set on one of runs 10-12 (run 12 alone for a calcite
%% case), for the campaign the controller's environment names
%% (Uncertainty_Campaign). The original set is not run here; it is the
%% Base or SSP campaign's own result.
%%
%% Everything the set does not move is the calibrated campaign's, and the
%% horizon, forcing and application years are those of the campaign the
%% set is compared with. Numerics are pinned to one thread as in every
%% Bouchet campaign, and the saved provenance is read back and checked
%% against the set. An element whose result already carries the set's
%% provenance exits at once, so a requeued or resubmitted array only runs
%% what is missing.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

here = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(here));
calibrated = fullfile(repo_root,'runtime','Base_1D');
addpath(here,calibrated,fullfile(repo_root,'v2','run'),fullfile(repo_root,'v2','T&C_Code'));

C = Uncertainty_Campaign('1D');
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
results_dir = fullfile(repo_root,'v2','SIM','Subfield_12',C.tag,member.tag);

Bouchet_Single_Thread_Runtime('Kantola_Uncertainty');
fingerprint = Bouchet_Campaign_Validation_Fingerprint( ...
    fullfile(repo_root,'v2'),here,C.tag,1,0,C.dependencies{:});
digest = Fingerprint_Digest(fingerprint);
setenv('TC_VALIDATION_FINGERPRINT',digest);
setenv('TC_CAMPAIGN_TAG',C.tag);
fprintf('Uncertainty set %s, run %d, %s -> %s\n',member.tag,run_id,C.tag,results_dir);

%% A result already there must be this set's, or it is in the way.
switch Uncertainty_Result_State(results_dir,run_id,C,K,C.tag,digest)
    case 'complete'
        fprintf('Set %s run %d (%s) is already complete with matching provenance; nothing to do.\n', ...
            member.tag,run_id,C.tag);
        return
    case 'stale'
        error('Kantola_Uncertainty:StaleResult', ...
            'Run %d in %s does not carry the provenance of set %s; remove it before resubmitting.', ...
            run_id,results_dir,member.tag);
end

if ~exist(results_dir,'dir'); mkdir(results_dir); end
options = [{'NDays',C.n_days,'StartDate',P.StartDate}, K.options, C.options, ...
    {'ParallelWorkers',0,'SaveResults',true,'ResultDirectory',results_dir}];
Run_One_Subfield(run_id,options{:});
if ~strcmp(Bouchet_Campaign_Validation_Fingerprint( ...
        fullfile(repo_root,'v2'),here,C.tag,1,0,C.dependencies{:}),fingerprint)
    error('Kantola_Uncertainty:SourceChangedDuringRun', ...
        'Code, inputs or numerical runtime changed while the set ran.');
end

%% Read the provenance back: the set must be what was run.
[state,result_file] = Uncertainty_Result_State(results_dir,run_id,C,K,C.tag,digest);
if ~strcmp(state,'complete')
    error('Kantola_Uncertainty:Provenance', ...
        'The saved result of run %d in %s does not carry the provenance of set %s.', ...
        run_id,results_dir,member.tag);
end
fprintf('PASS: set %s run %d (%s) saved with matching provenance: %s\n', ...
    member.tag,run_id,C.tag,result_file);
end
