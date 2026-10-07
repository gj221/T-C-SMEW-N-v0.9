function test_bouchet_validation_stamp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Campaign gate mechanics: stamps, run locks and existing-result vetting.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

root = tempname;
mkdir(root);
cleanup = onCleanup(@() rmdir(root,'s')); %#ok<NASGU>
stamp = fullfile(root,'nested','validation_stamp.txt');
Bouchet_Write_Validation_Stamp(stamp,'first');
assert(strcmp(strtrim(fileread(stamp)),'first'));
Bouchet_Write_Validation_Stamp(stamp,'second');
assert(strcmp(strtrim(fileread(stamp)),'second'));
assert_error(@() Bouchet_Write_Validation_Stamp(stamp,''), ...
    'Bouchet:ValidationStampInput');

%% ---------------------------------------------------------------------
%% Run locks
%% ---------------------------------------------------------------------
%% The lock, not the absence of a final MAT, is what stops two live jobs
%% from sharing one status file, checkpoint and result path.
old_job = getenv('SLURM_JOB_ID');
restore_job = onCleanup(@() setenv('SLURM_JOB_ID',old_job)); %#ok<NASGU>
setenv('SLURM_JOB_ID','');

lock = fullfile(root,'locks','run_10.lock');
held = Bouchet_Acquire_Run_Lock(lock);
assert(exist(lock,'dir') == 7,'The lock directory must exist while held.');
assert_error(@() Bouchet_Acquire_Run_Lock(lock),'Bouchet:RunLocked');
clear held
assert(exist(lock,'dir') ~= 7,'The lock must be released on cleanup.');
held = Bouchet_Acquire_Run_Lock(lock);
clear held
assert_error(@() Bouchet_Acquire_Run_Lock(''),'Bouchet:RunLockInput');

%% Both campaigns are submitted with --requeue. A requeued job keeps its
%% job id and finds the lock its killed predecessor left behind, so it must
%% reclaim that one and only that one; otherwise the rolling checkpoint is
%% unusable and every requeue needs a human.
orphan = fullfile(root,'locks','run_11.lock');
setenv('SLURM_JOB_ID','12345');
write_owner(orphan,'12345');
held = Bouchet_Acquire_Run_Lock(orphan);
assert(exist(orphan,'dir') == 7,'A requeued job must reclaim its own lock.');
clear held
write_owner(orphan,'99999');
assert_error(@() Bouchet_Acquire_Run_Lock(orphan),'Bouchet:RunLocked');
rmdir(orphan,'s');
setenv('SLURM_JOB_ID','');

%% ---------------------------------------------------------------------
%% Existing results
%% ---------------------------------------------------------------------
%% A filename match says nothing about which code wrote the file, so the
%% saved provenance decides whether a result may stand.
results = fullfile(root,'results');
mkdir(results);
assert(strcmp(Bouchet_Check_Existing_Result(results,10,'stampA',1,1826), ...
    'none'),'An empty campaign directory needs work.');

good = fullfile(results,'TC_2D_run10_control.mat');
write_result(good,'stampA',1,1826);
assert(strcmp(Bouchet_Check_Existing_Result(results,10,'stampA',1,1826), ...
    'complete'),'A matching result must be accepted.');
assert_error(@() Bouchet_Check_Existing_Result(results,10,'stampB',1,1826), ...
    'Bouchet:StaleResult');
assert_error(@() Bouchet_Check_Existing_Result(results,10,'stampA',12,1826), ...
    'Bouchet:WrongResultColumns');
assert_error(@() Bouchet_Check_Existing_Result(results,10,'stampA',1,315), ...
    'Bouchet:WrongResultHorizon');
assert_error(@() Bouchet_Check_Existing_Result(results,10,'',1,1826), ...
    'Bouchet:ExistingResultFingerprint');

duplicate = fullfile(results,'TC_2D_run10_control__alpha1p3.mat');
write_result(duplicate,'stampA',1,1826);
assert_error(@() Bouchet_Check_Existing_Result(results,10,'stampA',1,1826), ...
    'Bouchet:AmbiguousResult');
delete(duplicate);

legacy = fullfile(results,'TC_2D_run11_basalt.mat');
placeholder = 1; %#ok<NASGU>
save(legacy,'placeholder');
assert_error(@() Bouchet_Check_Existing_Result(results,11,'stampA',1,1826), ...
    'Bouchet:UnprovenancedResult');

fprintf(['PASS: Bouchet validation stamps, run locks and result vetting ' ...
    'behave.\n']);
end

function write_owner(lock_dir,job_id)
if ~exist(lock_dir,'dir'); mkdir(lock_dir); end
fid = fopen(fullfile(lock_dir,'owner.txt'),'wt');
fprintf(fid,'job=%s;array_task=;restart=;node=test;pid=1;acquired=now\n', ...
    job_id);
fclose(fid);
end

function write_result(file,fingerprint,n_columns,n_days)
run_provenance = struct( ...
    'validation_fingerprint',fingerprint, ...
    'n_columns',n_columns, ...
    'n_days',n_days, ...
    'n_timesteps',n_days*24); %#ok<NASGU>
save(file,'run_provenance');
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
