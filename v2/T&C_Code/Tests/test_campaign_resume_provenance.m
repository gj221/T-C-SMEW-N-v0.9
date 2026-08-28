function test_campaign_resume_provenance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% What a campaign run actually writes to disk, including after a resume.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The static campaign tests read source text. This one runs the solver for
%% a short horizon and inspects the saved MAT, because the resume path
%% restores the whole workspace from the checkpoint and had silently
%% overwritten both the run provenance and the caller's options.

test_dir = fileparts(mfilename('fullpath'));
code_dir = fileparts(test_dir);
v2_dir = fileparts(code_dir);
addpath(fullfile(v2_dir,'run'),code_dir, ...
    fullfile(code_dir,'Inputs'),fullfile(code_dir,'Parameters'));

scratch = tempname;
mkdir(scratch);
cleanup = onCleanup(@() rmdir(scratch,'s')); %#ok<NASGU>

env_names = {'TC_VALIDATION_FINGERPRINT','TC_NUMERICAL_IDENTITY', ...
    'TC_CAMPAIGN_TAG','TC_CKPT_MINUTES','TC_CKPT_FORCE_RESUME', ...
    'SLURM_JOB_ID'};
saved_env = cellfun(@getenv,env_names,'UniformOutput',false);
restore_env = onCleanup(@() restore(env_names,saved_env)); %#ok<NASGU>

stamp = 'TEST-CAMPAIGN-FINGERPRINT';
setenv('TC_VALIDATION_FINGERPRINT',stamp);
setenv('TC_NUMERICAL_IDENTITY','test_identity;parallel_workers=0');
setenv('TC_CAMPAIGN_TAG','TEST_RESUME_PROVENANCE');
setenv('TC_CKPT_FORCE_RESUME','');
n_days = 2;

%% ---------------------------------------------------------------------
%% The gate that production relies on must see a real result
%% ---------------------------------------------------------------------
setenv('TC_CKPT_MINUTES','0');
setenv('SLURM_JOB_ID','9001');
result_file = Run_One_Subfield(10,'NDays',n_days, ...
    'StartDate',datenum(2016,1,1),'AppeloAlpha',1.3, ...
    'ParallelWorkers',0,'SaveResults',true,'ResultDirectory',scratch);
p = getfield(load(result_file,'run_provenance'),'run_provenance');
assert(strcmp(p.validation_fingerprint,stamp), ...
    'The validation stamp must travel into the saved result.');
assert(p.n_days == n_days && p.n_timesteps == n_days*24 && ...
    p.n_columns == 1,'The saved horizon and column count are wrong.');
assert(strcmp(Bouchet_Check_Existing_Result( ...
    scratch,10,stamp,1,n_days),'complete'), ...
    'The campaign gate must accept a matching real result.');
assert_error(@() Bouchet_Check_Existing_Result(scratch,10,'OTHER',1,n_days), ...
    'Bouchet:StaleResult');
assert_error(@() Bouchet_Check_Existing_Result(scratch,10,stamp,1,1826), ...
    'Bouchet:WrongResultHorizon');

%% ---------------------------------------------------------------------
%% A resumed segment must not inherit the checkpoint's identity
%% ---------------------------------------------------------------------
%% The rolling checkpoint is only deleted once results are saved, so a
%% SaveResults = false pass leaves one for the next pass to resume from.
resume_dir = fullfile(scratch,'resume');
mkdir(resume_dir);
setenv('TC_CKPT_MINUTES','0.001');
setenv('SLURM_JOB_ID','9002');
Run_One_Subfield(10,'NDays',n_days,'StartDate',datenum(2016,1,1), ...
    'AppeloAlpha',1.3,'ParallelWorkers',0,'SaveResults',false, ...
    'ResultDirectory',resume_dir);
ckpt = dir(fullfile(resume_dir,'checkpoints','ckpt_run10_*.mat'));
assert(~isempty(ckpt),'The first pass wrote no rolling checkpoint.');

%% The checkpoint may only be replayed under the same numerical runtime.
fp = getfield(load(fullfile(ckpt(1).folder,ckpt(1).name), ...
    'ERW_ckpt_fingerprint'),'ERW_ckpt_fingerprint');
for token = {'matlab=','computer=','blas=','max_threads=', ...
        'parallel_workers='}
    assert(contains(fp,token{1}), ...
        'The checkpoint identity omits %s',token{1});
end

setenv('SLURM_JOB_ID','9003');
resumed_file = Run_One_Subfield(10,'NDays',n_days, ...
    'StartDate',datenum(2016,1,1),'AppeloAlpha',1.3, ...
    'ParallelWorkers',0,'SaveResults',true,'ResultDirectory',resume_dir);
%% SaveResults and ResultDirectory are not in the checkpoint fingerprint,
%% so a restored opt would send this result to the wrong place, or nowhere.
assert(~isempty(resumed_file) && isfile(resumed_file), ...
    'A resumed run must honour the options it was called with.');
assert(strcmp(fileparts(resumed_file),resume_dir), ...
    'A resumed run must write to the directory it was given.');

q = getfield(load(resumed_file,'run_provenance'),'run_provenance');
assert(isfield(q,'execution_segments') && ...
    numel(q.execution_segments) >= 1, ...
    'A resumed run must record its earlier execution segments.');
assert(strcmp(q.slurm_job_id,'9003'), ...
    'The final MAT must report the job that finished the run.');
assert(strcmp(q.execution_segments{1}.slurm_job_id,'9002'), ...
    'The earlier segment must keep its own job id.');
assert(q.resumed_from_iteration == n_days*24, ...
    'The resume point must be recorded.');
assert(strcmp(q.validation_fingerprint,stamp), ...
    'A resumed run must keep its validation stamp.');

%% A saved result must not leave its checkpoint behind.
assert(isempty(dir(fullfile(resume_dir,'checkpoints','ckpt_run10_*.mat'))), ...
    'A completed run must clear its rolling checkpoint.');

fprintf(['PASS: campaign results carry provenance and survive a resume.' ...
    '\n']);
end

function restore(names,values)
for k = 1:numel(names)
    setenv(names{k},values{k});
end
end

function assert_error(action,identifier)
try
    action();
catch ME
    assert(strcmp(ME.identifier,identifier), ...
        'Expected %s but received %s.',identifier,ME.identifier);
    return
end
error('Expected error %s was not raised.',identifier);
end
