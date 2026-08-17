function results = Sensitivity_Run(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run the vertical_1d trio under each nitrogen parameter set.
%%
%% Nothing in the model knows this exists. Each child is launched with
%% v2/Sensitivity on its path and the set's values in TC_BGP_* environment
%% variables; Biogeochemistry_Parameter picks them up through its one
%% optional hook. Children that are launched without this folder on the
%% path run the model's own parameters.
%%
%% Options
%%   Sets            indices into Sensitivity_Parameters   [all]
%%   Runs            treatment group                       [10:12]
%%   NDays           shortened run for a smoke test        []
%%   MaxConcurrent   simultaneous MATLAB processes         [4]
%%   OutputDirectory root for per-set result folders       [SIM/Sensitivity]
%%
%% Concurrency: each child holds about 2 GB. Four is the safe default on a
%% 16-thread box with roughly 9 GB free; raise it only after checking free
%% memory, because a swapping run is slower than a serial one.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sens_dir = fileparts(mfilename('fullpath'));
v2 = fileparts(sens_dir);
addpath(v2,sens_dir);

opt = struct('Sets',[],'Runs',10:12,'NDays',[],'MaxConcurrent',4, ...
    'OutputDirectory',fullfile(v2,'SIM','Sensitivity'), ...
    'Resume',true,'MinFreeGB',3);
if mod(numel(varargin),2) ~= 0
    error('Sensitivity_Run:BadOptions','Options must be name/value pairs.');
end
for k = 1:2:numel(varargin)
    name = char(varargin{k});
    if ~isfield(opt,name)
        error('Sensitivity_Run:BadOptions','Unknown option %s.',name);
    end
    opt.(name) = varargin{k+1};
end

[allowed,sets] = Sensitivity_Parameters();
if isempty(opt.Sets); opt.Sets = 1:numel(sets); end
sets = sets(opt.Sets);
RM = Run_Matrix_12();
run_ids = opt.Runs(:).';

%% Validate every set before launching anything. A four-hour batch that
%% dies on the last set because of a typo is a bad way to find out.
for s = 1:numel(sets)
    f = setdiff(fieldnames(sets{s}),{'tag','note'});
    bad = setdiff(f,allowed);
    if ~isempty(bad)
        error('Sensitivity_Run:UnknownParameter', ...
            'Set %s changes %s, which is not sweepable.', ...
            sets{s}.tag,strjoin(bad,', '));
    end
    for j = 1:numel(f)
        v = sets{s}.(f{j});
        if ~(isnumeric(v) && isscalar(v) && isfinite(v))
            error('Sensitivity_Run:BadValue', ...
                'Set %s gives a non-finite value for %s.',sets{s}.tag,f{j});
        end
    end
end

if ~isfolder(opt.OutputDirectory); mkdir(opt.OutputDirectory); end
matlab_exe = fullfile(matlabroot,'bin','matlab.exe');

%% Build the full job list: every set crossed with every run in the group.
job = struct('set',{},'tag',{},'run_id',{},'name',{},'dir',{}, ...
    'env',{},'result',{},'log',{});
for s = 1:numel(sets)
    sdir = fullfile(opt.OutputDirectory,sets{s}.tag);
    if ~isfolder(sdir); mkdir(sdir); end
    f = setdiff(fieldnames(sets{s}),{'tag','note'});
    env = cell(1,numel(f));
    for j = 1:numel(f)
        env{j} = sprintf('TC_BGP_%s=%.17g',f{j},sets{s}.(f{j}));
    end
    for r = run_ids
        suffix = '';
        if ~isempty(opt.NDays); suffix = sprintf('__%dd',opt.NDays); end
        job(end+1) = struct('set',s,'tag',sets{s}.tag,'run_id',r, ...
            'name',RM(r).name,'dir',sdir,'env',{env}, ...
            'result',fullfile(sdir,sprintf('TC_2D_run%02d_%s%s.mat', ...
                r,RM(r).name,suffix)), ...
            'log',fullfile(sdir,sprintf('run_%02d_%s%s.log', ...
                r,RM(r).name,suffix))); %#ok<AGROW>
    end
end

fprintf('\n=== Nitrogen sensitivity ===\n');
fprintf('sets: %d   runs per set: %d   jobs: %d   max concurrent: %d\n', ...
    numel(sets),numel(run_ids),numel(job),opt.MaxConcurrent);
for s = 1:numel(sets)
    fprintf('  %-20s %s\n',sets{s}.tag,sets{s}.note);
end
fprintf('results: %s\n\n',opt.OutputDirectory);

%% Checkpointing. Every completed run is already an atomic file on disk:
%% Run_One_Subfield saves to a tempname and then movefiles it, so a killed
%% process leaves either a whole result or none, never a truncated one.
%% What was missing was the ability to pick up again, so a crash six hours
%% in cost the whole batch. Resume skips jobs whose result is present and
%% loadable, which caps the loss at the runs actually in flight.
queue = 1:numel(job);
if opt.Resume
    done_already = false(1,numel(job));
    for q = 1:numel(job)
        done_already(q) = result_is_complete(job(q).result);
    end
    if any(done_already)
        fprintf('  resuming: %d of %d jobs already complete\n', ...
            sum(done_already),numel(job));
        queue = queue(~done_already);
    end
else
    for q = 1:numel(job)
        if isfile(job(q).result); delete(job(q).result); end
    end
end
%% Sweep up temp files orphaned by a previous kill mid-save.
for s = 1:numel(sets)
    sweep_orphans(fullfile(opt.OutputDirectory,sets{s}.tag));
end

running = zeros(0,1);
proc = cell(numel(job),1);
started = zeros(numel(job),1);
failed = zeros(0,1);
manifest = fullfile(opt.OutputDirectory,'progress.csv');
write_manifest(manifest,job,queue,failed);
while ~isempty(queue) || ~isempty(running)
    while ~isempty(queue) && numel(running) < opt.MaxConcurrent
        %% Hold back if memory is short. Launching into swap makes every
        %% running job slower, so waiting is strictly better than starting.
        if ~isempty(running) && free_gb() < opt.MinFreeGB
            fprintf('  [%s] holding at %d concurrent, %.1f GB free\n', ...
                datestr(now,'HH:MM:SS'),numel(running),free_gb());
            break
        end
        q = queue(1); queue(1) = [];
        proc{q} = launch(matlab_exe,job(q),v2,sens_dir,opt);
        started(q) = now;
        running(numel(running)+1,1) = q; %#ok<AGROW>
        fprintf('  [%s] launch %-20s run %d  (%.1f GB free)\n', ...
            datestr(now,'HH:MM:SS'),job(q).tag,job(q).run_id,free_gb());
    end
    pause(15);
    done = false(size(running));
    for i = 1:numel(running)
        q = running(i);
        if process_alive(proc{q}); continue; end
        done(i) = true;
        if isfile(job(q).result)
            fprintf('  [%s] done   %-20s run %d  (%.2f h)\n', ...
                datestr(now,'HH:MM:SS'),job(q).tag,job(q).run_id, ...
                (now-started(q))*24);
        else
            fprintf(2,'  [%s] FAILED %-20s run %d; see %s\n', ...
                datestr(now,'HH:MM:SS'),job(q).tag,job(q).run_id,job(q).log);
            failed(end+1,1) = q; %#ok<AGROW>
        end
    end
    running(done) = [];
    running = running(:);
    write_manifest(manifest,job,queue,failed);
end

results = struct('sets',{sets},'job',job,'failed',failed, ...
    'OutputDirectory',opt.OutputDirectory);
if ~isempty(failed)
    fprintf(2,'\n%d of %d jobs failed.\n',numel(failed),numel(job));
else
    fprintf('\nall %d jobs completed.\n',numel(job));
end
fprintf('analyse with: Analyse_Sensitivity\n');

end

%% =======================================================================
function alive = process_alive(p)
%% .NET on Windows, Java elsewhere; both report liveness differently.
if ispc
    alive = ~p.HasExited;
else
    alive = p.isAlive();
end
end

%% =======================================================================
function p = launch(matlab_exe,j,v2,sens_dir,opt)
%% One detached MATLAB per job, launched the same way Run_All_12_Parallel
%% does. The parameter values are applied with setenv inside the child, so
%% they belong to that process alone; a sibling running a different set
%% cannot see them. v2/Sensitivity goes on the child's path here and only
%% here, which is what arms the hook in Biogeochemistry_Parameter.
setenv_cmd = '';
for k = 1:numel(j.env)
    kv = strsplit(j.env{k},'=');
    setenv_cmd = [setenv_cmd sprintf('setenv(''%s'',''%s'');',kv{1},kv{2})]; %#ok<AGROW>
end
nd = '';
if ~isempty(opt.NDays); nd = sprintf(',''NDays'',%d',opt.NDays); end
batch = sprintf(['cd(''%s''); addpath(''%s'',''%s''); %s' ...
    'Run_One_Subfield(%d,''ResultDirectory'',''%s''%s);'], ...
    v2,v2,sens_dir,setenv_cmd,j.run_id,j.dir,nd);
if ispc
    cmd_exe = fullfile(getenv('SystemRoot'),'System32','cmd.exe');
    args = sprintf('/d /s /c ""%s" -batch "%s" > "%s" 2>&1"', ...
        matlab_exe,batch,j.log);
    si = System.Diagnostics.ProcessStartInfo(cmd_exe,args);
    si.UseShellExecute = false;
    si.CreateNoWindow = true;
    p = System.Diagnostics.Process.Start(si);
else
    cmdlist = java.util.ArrayList;
    cmdlist.add(matlab_exe); cmdlist.add('-batch'); cmdlist.add(batch);
    pb = java.lang.ProcessBuilder(cmdlist);
    pb.redirectErrorStream(true);
    pb.redirectOutput(java.io.File(j.log));
    p = pb.start();
end
end

%% =======================================================================
function ok = result_is_complete(f)
%% A completed run is an atomic file: Run_One_Subfield saves to a tempname
%% and movefiles it into place, so a present file is a whole file. whos
%% reads only the variable table, which catches a corrupt or truncated
%% file without loading the gigabytes of data it contains.
ok = false;
if ~isfile(f); return; end
try
    ok = ~isempty(whos('-file',f));
catch
    ok = false;
end
end

%% =======================================================================
function sweep_orphans(d)
%% Remove tempname leftovers from runs killed mid-save. Every result is
%% named TC_2D_run*, so any other .mat in the set directory is an orphan.
if ~isfolder(d); return; end
names = dir(fullfile(d,'*.mat'));
for k = 1:numel(names)
    if ~names(k).isdir && ~startsWith(names(k).name,'TC_2D_run')
        try
            delete(fullfile(d,names(k).name));
        catch
        end
    end
end
end

%% =======================================================================
function write_manifest(manifest,job,queue,failed)
%% One row per job, refreshed every poll so a crash mid-batch leaves a
%% readable picture of what finished and what did not.
fid = fopen(manifest,'w');
if fid < 0; return; end
fprintf(fid,'tag,run_id,name,result,status\n');
for q = 1:numel(job)
    if any(failed == q)
        st = 'failed';
    elseif isfile(job(q).result)
        st = 'done';
    elseif any(queue == q)
        st = 'pending';
    else
        st = 'running';
    end
    fprintf(fid,'%s,%d,%s,%s,%s\n',job(q).tag,job(q).run_id, ...
        job(q).name,job(q).result,st);
end
fclose(fid);
end

%% =======================================================================
function gb = free_gb()
%% Physical memory available for new arrays, in GB. MATLAB's memory()
%% reports this directly on Windows; elsewhere fall back to /proc/meminfo.
gb = inf;
if ispc
    try
        m = memory;
        gb = m.MemAvailableAllArrays/1e9;
        return
    catch
    end
end
try
    [~,txt] = system('free -b');
    tok = regexp(txt,'Mem:\s+\d+\s+\d+\s+(\d+)','tokens','once');
    if ~isempty(tok); gb = str2double(tok{1})/1e9; end
catch
end
end
