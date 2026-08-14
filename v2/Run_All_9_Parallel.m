function plan = Run_All_9_Parallel(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Launch Run_Matrix_9 as isolated MATLAB processes.
%%
%% Options
%%   MaxConcurrent   maximum live child processes                 [3]
%%   Runs            matrix ids                                   [1:9]
%%   NDays           optional shortened run                       []
%%   GridNx, GridNy  optional main-grid overrides                 []
%%   PollSeconds     status polling interval                      [20]
%%   MaxRuntimeHours stop waiting after this time per child       [72]
%%   DryRun          return and print commands without launching  [false]
%%   Overwrite       rerun when the result already exists         [false]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

opt = struct('MaxConcurrent',3,'Runs',1:9,'NDays',[], ...
    'GridNx',[],'GridNy',[],'PollSeconds',20, ...
    'MaxRuntimeHours',72,'DryRun',false,'Overwrite',false);
if mod(numel(varargin),2) ~= 0
    error('Run_All_9_Parallel:BadOptions','Options must be name/value pairs.');
end
for k = 1:2:numel(varargin)
    name = char(varargin{k});
    if ~isfield(opt,name)
        error('Run_All_9_Parallel:BadOptions','Unknown option %s.',name);
    end
    opt.(name) = varargin{k+1};
end
if ~(isnumeric(opt.MaxConcurrent) && isscalar(opt.MaxConcurrent) && ...
        isfinite(opt.MaxConcurrent) && opt.MaxConcurrent == fix(opt.MaxConcurrent) && ...
        opt.MaxConcurrent >= 1)
    error('Run_All_9_Parallel:BadMaxConcurrent','MaxConcurrent must be positive.');
end
if ~(isnumeric(opt.Runs) && isvector(opt.Runs) && ...
        all(isfinite(opt.Runs)) && all(opt.Runs == fix(opt.Runs)) && ...
        all(ismember(opt.Runs,1:9)) && numel(unique(opt.Runs)) == numel(opt.Runs))
    error('Run_All_9_Parallel:BadRuns','Runs must contain unique integer ids from 1 to 9.');
end
validateattributes(opt.PollSeconds,{'numeric'},{'scalar','positive','finite'});
validateattributes(opt.MaxRuntimeHours,{'numeric'},{'scalar','positive','finite'});
validate_optional_integer(opt.NDays,'NDays');
validate_optional_integer(opt.GridNx,'GridNx');
validate_optional_integer(opt.GridNy,'GridNy');

here = fileparts(mfilename('fullpath'));
RM = Run_Matrix_9();
runs = opt.Runs(:)';
log_dir = fullfile(here,'logs','Subfield_9');
out_dir = fullfile(here,'SIM','Subfield_9');
status_dir = fullfile(log_dir,'status');

extra = {};
if ~isempty(opt.NDays);  extra = [extra {'NDays',opt.NDays}]; end
if ~isempty(opt.GridNx); extra = [extra {'GridNx',opt.GridNx}]; end
if ~isempty(opt.GridNy); extra = [extra {'GridNy',opt.GridNy}]; end

matlab_exe = fullfile(matlabroot,'bin','matlab.exe');
if ~isfile(matlab_exe); matlab_exe = fullfile(matlabroot,'bin','matlab'); end
plan = repmat(struct('run_id',0,'name','','result_file','', ...
    'status_file','','log_file','','batch','','command',''),numel(runs),1);
for q = 1:numel(runs)
    id = runs(q);
    tag = RM(id).name;
    if ~isempty(opt.GridNx) || ~isempty(opt.GridNy)
        nx = RM(id).geom.main_grid_nx;
        ny = RM(id).geom.main_grid_ny;
        if ~isempty(opt.GridNx); nx = opt.GridNx; end
        if ~isempty(opt.GridNy); ny = opt.GridNy; end
        tag = sprintf('%s__grid%dx%d',tag,nx,ny);
    end
    if ~isempty(opt.NDays); tag = sprintf('%s__%dd',tag,opt.NDays); end
    plan(q).run_id = id;
    plan(q).name = RM(id).name;
    plan(q).result_file = fullfile(out_dir, ...
        sprintf('TC_2D_run%02d_%s.mat',id,tag));
    plan(q).status_file = fullfile(status_dir, ...
        sprintf('run_%02d_%s.status',id,tag));
    plan(q).log_file = fullfile(log_dir, ...
        sprintf('run_%02d_%s.log',id,tag));
    args = worker_args(id,plan(q).status_file,extra);
    plan(q).batch = sprintf("cd('%s'); Run_Subfield_Worker(%s)", ...
        strrep(here,"'","''"),args);
    plan(q).command = sprintf('"%s" -batch "%s" > "%s" 2>&1', ...
        matlab_exe,plan(q).batch,plan(q).log_file);
end

fprintf('\n=== Sub-field simulations ===\n');
fprintf('runs: %s   max concurrent: %d\n',mat2str(runs),opt.MaxConcurrent);
fprintf('results: %s\nlogs: %s\n',out_dir,log_dir);
if opt.DryRun
    for q = 1:numel(plan)
        fprintf('  run %d  %s\n    %s\n',plan(q).run_id,plan(q).name,plan(q).command);
    end
    return
end

if ~exist(log_dir,'dir'); mkdir(log_dir); end
if ~exist(status_dir,'dir'); mkdir(status_dir); end
if ~exist(out_dir,'dir'); mkdir(out_dir); end

pending = true(size(plan));
if ~opt.Overwrite
    for q = 1:numel(plan)
        if isfile(plan(q).result_file)
            pending(q) = false;
            fprintf('  existing result, skip run %d  %s\n',plan(q).run_id,plan(q).name);
        end
    end
end
queue = find(pending);
running = zeros(0,1);
started = nan(size(plan));
failed = zeros(0,1);
process = cell(size(plan));
process_guard = cell(size(plan));

while ~isempty(queue) || ~isempty(running)
    while ~isempty(queue) && numel(running) < opt.MaxConcurrent
        q = queue(1); queue(1) = [];
        if isfile(plan(q).status_file); delete(plan(q).status_file); end
        if opt.Overwrite && isfile(plan(q).result_file); delete(plan(q).result_file); end
        process{q} = launch_process(matlab_exe,plan(q).batch,plan(q).log_file);
        process_guard{q} = onCleanup(@() terminate_process(process{q}));
        started(q) = now;
        running(end+1,1) = q; %#ok<AGROW>
        fprintf('  [%s] launch run %d  %s\n',datestr(now,'HH:MM:SS'), ...
            plan(q).run_id,plan(q).name);
    end
    pause(opt.PollSeconds);
    finished = false(size(running));
    for r = 1:numel(running)
        q = running(r);
        state = read_state(plan(q).status_file);
        if strcmp(state,'complete') && isfile(plan(q).result_file)
            fprintf('  [%s] done run %d  %s\n',datestr(now,'HH:MM:SS'), ...
                plan(q).run_id,plan(q).name);
            finished(r) = true;
        elseif strcmp(state,'failed')
            fprintf(2,'  [%s] FAILED run %d  %s; see %s\n', ...
                datestr(now,'HH:MM:SS'),plan(q).run_id,plan(q).name,plan(q).log_file);
            failed(end+1,1) = q; %#ok<AGROW>
            finished(r) = true;
        elseif ~process{q}.isAlive()
            fprintf(2,'  [%s] FAILED run %d  %s (child exit %d); see %s\n', ...
                datestr(now,'HH:MM:SS'),plan(q).run_id,plan(q).name, ...
                process{q}.exitValue(),plan(q).log_file);
            failed(end+1,1) = q; %#ok<AGROW>
            finished(r) = true;
        elseif (now-started(q))*24 > opt.MaxRuntimeHours
            process{q}.destroyForcibly();
            fprintf(2,'  [%s] TIMEOUT run %d  %s after %.1f h\n', ...
                datestr(now,'HH:MM:SS'),plan(q).run_id,plan(q).name,opt.MaxRuntimeHours);
            failed(end+1,1) = q; %#ok<AGROW>
            finished(r) = true;
        end
    end
    running(finished) = [];
end

if ~isempty(failed)
    error('Run_All_9_Parallel:ChildFailed','Run(s) %s failed or timed out.', ...
        mat2str([plan(failed).run_id]));
end
fprintf('All requested sub-field runs completed.\n');

end

function args = worker_args(run_id,status_file,extra)
parts = {sprintf('%d',run_id),sprintf('''%s''',strrep(status_file,'''',''''''))};
for k = 1:2:numel(extra)
    parts{end+1} = sprintf('''%s''',extra{k}); %#ok<AGROW>
    parts{end+1} = mat2str(extra{k+1}); %#ok<AGROW>
end
args = strjoin(parts,',');
end

function process = launch_process(matlab_exe,batch,log_file)
command = java.util.ArrayList;
command.add(matlab_exe);
command.add('-batch');
command.add(batch);
builder = java.lang.ProcessBuilder(command);
builder.redirectErrorStream(true);
builder.redirectOutput(java.io.File(log_file));
process = builder.start();
end

function state = read_state(status_file)
state = '';
if ~isfile(status_file); return; end
txt = fileread(status_file);
tok = regexp(txt,'(?m)^state=([^\r\n]+)','tokens','once');
if ~isempty(tok); state = tok{1}; end
end

function terminate_process(process)
try
    if ~isempty(process) && process.isAlive(); process.destroyForcibly(); end
catch
end
end

function validate_optional_integer(value,name)
if isempty(value); return; end
if ~(isnumeric(value) && isscalar(value) && isfinite(value) && ...
        value == fix(value) && value >= 1)
    error('Run_All_9_Parallel:BadOptions','%s must be a positive integer.',name);
end
end
