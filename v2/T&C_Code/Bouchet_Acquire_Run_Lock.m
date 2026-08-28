function cleanup = Bouchet_Acquire_Run_Lock(lock_dir)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Claim exclusive ownership of one campaign treatment.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The absence of a final MAT only proves that no run has finished. It does
%% not stop a second submission of the same treatment, and two live jobs
%% would then share one status file, one rolling checkpoint and one result
%% path. Directory creation is atomic on both POSIX and NTFS, so exactly
%% one caller can create a given lock and the loser fails here instead of
%% corrupting the trajectory.
%%
%% Both campaigns are submitted with --requeue, and a requeued job keeps
%% its job id while its predecessor is killed without releasing anything.
%% A lock left behind by the same job id is therefore reclaimed, which is
%% what makes the rolling checkpoint usable. A lock from any other job is
%% not: that is either a live duplicate or a crash that needs a human.
%%
%% INPUT
%%   lock_dir  full path of the lock directory for this run
%% OUTPUT
%%   cleanup   onCleanup object; the lock is released when it is destroyed

if ~(ischar(lock_dir) || (isstring(lock_dir) && isscalar(lock_dir))) || ...
        strlength(string(lock_dir)) == 0
    error('Bouchet:RunLockInput','The lock path must be nonempty text.');
end
lock_dir = char(lock_dir);
parent = fileparts(lock_dir);
if ~isempty(parent) && ~exist(parent,'dir'); mkdir(parent); end

owner_file = fullfile(lock_dir,'owner.txt');
[made,msg,msgid] = mkdir(lock_dir);
if ~made
    error('Bouchet:RunLockUnavailable', ...
        'Cannot create the run lock %s: %s',lock_dir,msg);
end
if strcmp(msgid,'MATLAB:MKDIR:DirectoryExists')
    if isfile(owner_file); owner = strtrim(fileread(owner_file));
    else;                  owner = '';
    end
    this_job = strtrim(getenv('SLURM_JOB_ID'));
    held_by = regexp(owner,'job=([^;]*)','tokens','once');
    same_job = ~isempty(this_job) && ~isempty(held_by) && ...
        strcmp(strtrim(held_by{1}),this_job);
    if ~same_job
        if isempty(owner); owner = 'unknown'; end
        error('Bouchet:RunLocked', ...
            ['The run lock %s is already held by {%s}. Another job is ' ...
             'running this treatment, or a previous job was killed. ' ...
             'Remove the lock only after confirming with squeue that no ' ...
             'job is active.'],lock_dir,owner);
    end
    fprintf('Reclaiming the run lock %s left by this job (restart %s).\n', ...
        lock_dir,getenv('SLURM_RESTART_COUNT'));
end

fid = fopen(owner_file,'wt');
if fid >= 0
    fprintf(fid,['job=%s;array_task=%s;restart=%s;node=%s;pid=%d;' ...
        'acquired=%s\n'], ...
        getenv('SLURM_JOB_ID'),getenv('SLURM_ARRAY_TASK_ID'), ...
        getenv('SLURM_RESTART_COUNT'),getenv('SLURMD_NODENAME'), ...
        feature('getpid'),datestr(now,'yyyy-mm-ddTHH:MM:SS'));
    fclose(fid);
end
cleanup = onCleanup(@() release_lock(lock_dir));

end

function release_lock(lock_dir)
if exist(lock_dir,'dir')
    [ok,msg] = rmdir(lock_dir,'s');
    if ~ok
        warning('Bouchet:RunLockRelease', ...
            'Cannot remove the run lock %s: %s',lock_dir,msg);
    end
end
end
