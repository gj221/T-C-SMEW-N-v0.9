function plan = Run_All_12_Parallel(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Launch Run_Matrix_12 as isolated MATLAB processes.
%%
%% Options
%%   MaxConcurrent   maximum live child processes                 [2]
%%   Runs            matrix ids                         [10:12 1:9]
%%   NDays           optional shortened run                       []
%%   StartDate      optional forcing start date                    []
%%   GridNx, GridNy  optional main-grid overrides                 []
%%   AppeloAlpha     optional run-local Appelo alpha              []
%%   AppeloLogKH     optional [Na K Ca Mg Al] logKH vector        []
%%   KNmaxd          optional denitrification rate [d^-1]         []
%%   Kbot            optional bottom boundary [mm h^-1]           []
%%   ParallelWorkers daily ERW workers inside each run             []
%%   ResultDirectory optional result and checkpoint directory      ''
%%   Scenario        forcing suffix, Data_<Location>_<Scenario>    ''
%%   CalciteCase     ERW workbook row for the calcite runs only    []
%%   LastApplicationYear  last 10 November rock application     [] (2019)
%%   PreferenceDirectoryRoot optional root for isolated MATLAB prefs ''
%%   PollSeconds     status polling interval                      [20]
%%   MaxRuntimeHours stop waiting after this time per child      [168]
%%   DryRun          return and print commands without launching  [false]
%%   Overwrite       rerun when the result already exists         [false]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

opt = struct('MaxConcurrent',2,'Runs',[10:12 1:9],'NDays',[],'StartDate',[], ...
    'GridNx',[],'GridNy',[],'AppeloAlpha',[],'AppeloLogKH',[], ...
    'KNmaxd',[],'Kbot',[], ...
    'ParallelWorkers',[],'PollSeconds',20, ...
    'MaxRuntimeHours',168,'DryRun',false,'Overwrite',false, ...
    'ResultDirectory','','LogDirectory','','PreferenceDirectoryRoot','', ...
    'Scenario','','CalciteCase',[],'LastApplicationYear',[]);
if mod(numel(varargin),2) ~= 0
    error('Run_All_12_Parallel:BadOptions','Options must be name/value pairs.');
end
for k = 1:2:numel(varargin)
    name = char(varargin{k});
    if ~isfield(opt,name)
        error('Run_All_12_Parallel:BadOptions','Unknown option %s.',name);
    end
    opt.(name) = varargin{k+1};
end
if ~(islogical(opt.DryRun) && isscalar(opt.DryRun) && ...
        islogical(opt.Overwrite) && isscalar(opt.Overwrite))
    error('Run_All_12_Parallel:BadOptions', ...
        'DryRun and Overwrite must be scalar logical values.');
end
if ~(isnumeric(opt.MaxConcurrent) && isscalar(opt.MaxConcurrent) && ...
        isfinite(opt.MaxConcurrent) && opt.MaxConcurrent == fix(opt.MaxConcurrent) && ...
        opt.MaxConcurrent >= 1)
    error('Run_All_12_Parallel:BadMaxConcurrent','MaxConcurrent must be positive.');
end
if ~(isnumeric(opt.Runs) && isvector(opt.Runs) && ...
        all(isfinite(opt.Runs)) && all(opt.Runs == fix(opt.Runs)) && ...
        all(ismember(opt.Runs,1:12)) && numel(unique(opt.Runs)) == numel(opt.Runs))
    error('Run_All_12_Parallel:BadRuns', ...
        'Runs must contain unique integer ids from 1 to 12.');
end
validateattributes(opt.PollSeconds,{'numeric'},{'scalar','positive','finite'});
validateattributes(opt.MaxRuntimeHours,{'numeric'},{'scalar','positive','finite'});
validate_optional_integer(opt.NDays,'NDays');
validate_optional_integer(opt.GridNx,'GridNx');
validate_optional_integer(opt.GridNy,'GridNy');
validate_optional_nonnegative_scalar(opt.AppeloAlpha,'AppeloAlpha');
validate_optional_nonnegative_scalar(opt.KNmaxd,'KNmaxd');
validate_optional_kbot(opt.Kbot);
validate_optional_nonnegative_integer(opt.ParallelWorkers,'ParallelWorkers');
if ~isempty(opt.AppeloLogKH) && ~(isnumeric(opt.AppeloLogKH) && ...
        isvector(opt.AppeloLogKH) && numel(opt.AppeloLogKH) == 5 && ...
        all(isfinite(opt.AppeloLogKH)))
    error('Run_All_12_Parallel:BadOptions', ...
        'AppeloLogKH must be five finite values ordered [Na K Ca Mg Al].');
end
opt.AppeloLogKH = reshape(opt.AppeloLogKH,1,[]);
if ~isempty(opt.StartDate)
    try
        opt.StartDate = datenum(opt.StartDate);
    catch
        error('Run_All_12_Parallel:BadStartDate', ...
            'StartDate must be a scalar date accepted by datenum.');
    end
    if ~(isnumeric(opt.StartDate) && isscalar(opt.StartDate) && isfinite(opt.StartDate))
        error('Run_All_12_Parallel:BadStartDate', ...
            'StartDate must be a scalar date accepted by datenum.');
    end
end
if ~(ischar(opt.ResultDirectory) || ...
        (isstring(opt.ResultDirectory) && isscalar(opt.ResultDirectory))) || ...
        ~(ischar(opt.LogDirectory) || ...
        (isstring(opt.LogDirectory) && isscalar(opt.LogDirectory))) || ...
        ~(ischar(opt.PreferenceDirectoryRoot) || ...
        (isstring(opt.PreferenceDirectoryRoot) && isscalar(opt.PreferenceDirectoryRoot)))
    error('Run_All_12_Parallel:BadOptions', ...
        ['ResultDirectory, LogDirectory and PreferenceDirectoryRoot must be ' ...
        'character vectors or scalar strings.']);
end

%% This launcher and its worker live in v2/run; results, logs and status
%% files belong to the repository root one level up.
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
RM = Run_Matrix_12();
runs = opt.Runs(:)';
requires_explicit_alpha = any(ismember(runs,10:12)) && ...
    (isempty(opt.NDays) || opt.NDays >= 679);
if requires_explicit_alpha && isempty(opt.AppeloAlpha)
    error('Run_All_12_Parallel:MissingAppeloAlpha', ...
        'Runs 10-12 reaching day 679 require an explicit prescribed AppeloAlpha.');
end
log_dir = char(opt.LogDirectory);
if isempty(log_dir); log_dir = fullfile(root,'logs','Subfield_12'); end
out_dir = char(opt.ResultDirectory);
if isempty(out_dir); out_dir = fullfile(root,'SIM','Subfield_12'); end
status_dir = fullfile(log_dir,'status');
pref_root = char(opt.PreferenceDirectoryRoot);

extra = {};
if ~isempty(opt.NDays);  extra = [extra {'NDays',opt.NDays}]; end
if ~isempty(opt.StartDate); extra = [extra {'StartDate',opt.StartDate}]; end
if ~isempty(opt.GridNx); extra = [extra {'GridNx',opt.GridNx}]; end
if ~isempty(opt.GridNy); extra = [extra {'GridNy',opt.GridNy}]; end
if ~isempty(opt.AppeloAlpha)
    extra = [extra {'AppeloAlpha',opt.AppeloAlpha}];
end
if ~isempty(opt.AppeloLogKH)
    extra = [extra {'AppeloLogKH',opt.AppeloLogKH}];
end
if ~isempty(opt.KNmaxd)
    extra = [extra {'KNmaxd',opt.KNmaxd}];
end
if ~isempty(opt.Kbot)
    extra = [extra {'Kbot',opt.Kbot}];
end
if ~isempty(opt.ParallelWorkers)
    extra = [extra {'ParallelWorkers',opt.ParallelWorkers}];
end
if ~isempty(opt.Scenario)
    extra = [extra {'Scenario',opt.Scenario}];
end
if ~isempty(opt.CalciteCase)
    extra = [extra {'CalciteCase',opt.CalciteCase}];
end
if ~isempty(opt.LastApplicationYear)
    extra = [extra {'LastApplicationYear',opt.LastApplicationYear}];
end
if ~isempty(opt.ResultDirectory)
    extra = [extra {'ResultDirectory',out_dir}];
end
matlab_exe = fullfile(matlabroot,'bin','matlab.exe');
if ~isfile(matlab_exe); matlab_exe = fullfile(matlabroot,'bin','matlab'); end
plan = repmat(struct('run_id',0,'name','','result_file','', ...
    'status_file','','log_file','','runtime_directory','', ...
    'preference_directory','','temporary_directory','', ...
    'batch','','command',''),numel(runs),1);
for q = 1:numel(runs)
    id = runs(q);
    tag = result_tag(RM(id),opt);
    plan(q).run_id = id;
    plan(q).name = RM(id).name;
    plan(q).result_file = fullfile(out_dir, ...
        sprintf('TC_2D_run%02d_%s.mat',id,tag));
    plan(q).status_file = fullfile(status_dir, ...
        sprintf('run_%02d_%s.status',id,tag));
    plan(q).log_file = fullfile(log_dir, ...
        sprintf('run_%02d_%s.log',id,tag));
    if ~isempty(pref_root)
        plan(q).runtime_directory = fullfile(pref_root, ...
            sprintf('run_%02d',id));
        plan(q).preference_directory = fullfile(plan(q).runtime_directory,'pref');
        plan(q).temporary_directory = fullfile(plan(q).runtime_directory,'temp');
    end
    args = worker_args(id,plan(q).status_file,extra);
    plan(q).batch = sprintf("cd('%s'); Run_Subfield_Worker(%s)", ...
        strrep(here,"'","''"),args);
    if isempty(plan(q).preference_directory)
        plan(q).command = sprintf('"%s" -batch "%s" > "%s" 2>&1', ...
            matlab_exe,plan(q).batch,plan(q).log_file);
    else
        plan(q).command = sprintf( ...
            ['set "MATLAB_PREFDIR=%s" && set "TEMP=%s" && set "TMP=%s" && ' ...
             '"%s" -batch "%s" > "%s" 2>&1'], ...
            plan(q).preference_directory,plan(q).temporary_directory, ...
            plan(q).temporary_directory,matlab_exe,plan(q).batch,plan(q).log_file);
    end
end

fprintf('\n=== Twelve simulations ===\n');
fprintf('runs: %s   max concurrent: %d\n',mat2str(runs),opt.MaxConcurrent);
fprintf('max runtime per run: %.0f h\n',opt.MaxRuntimeHours);
fprintf('results: %s\nlogs: %s\n',out_dir,log_dir);
if opt.DryRun
    for q = 1:numel(plan)
        fprintf('  run %d  %s\n    %s\n', ...
            plan(q).run_id,plan(q).name,plan(q).command);
    end
    return
end

if ~exist(log_dir,'dir'); mkdir(log_dir); end
if ~exist(status_dir,'dir'); mkdir(status_dir); end
if ~exist(out_dir,'dir'); mkdir(out_dir); end
if ~isempty(pref_root) && ~exist(pref_root,'dir'); mkdir(pref_root); end

pending = true(size(plan));
if ~opt.Overwrite
    for q = 1:numel(plan)
        if isfile(plan(q).result_file)
            pending(q) = false;
            fprintf('  existing result, skip run %d  %s\n', ...
                plan(q).run_id,plan(q).name);
        end
    end
end
queue = find(pending);
running = zeros(0,1);
started = nan(size(plan));
failed = zeros(0,1);
process = cell(size(plan));
process_guard = cell(size(plan));

%% A failure is almost always systemic rather than run-specific, so the
%% remaining scenarios would hit the same wall hours later. Stop the
%% whole campaign on the first failure instead of burning the queue.
abort = false;
while (~isempty(queue) || ~isempty(running)) && ~abort
    while ~isempty(queue) && numel(running) < opt.MaxConcurrent && ~abort
        q = queue(1); queue(1) = [];
        if isfile(plan(q).status_file); delete(plan(q).status_file); end
        if opt.Overwrite && isfile(plan(q).result_file); delete(plan(q).result_file); end
        if ~isempty(plan(q).preference_directory) && ...
                ~exist(plan(q).preference_directory,'dir')
            mkdir(plan(q).preference_directory);
        end
        if ~isempty(plan(q).temporary_directory) && ...
                ~exist(plan(q).temporary_directory,'dir')
            mkdir(plan(q).temporary_directory);
        end
        process{q} = launch_process(matlab_exe,plan(q).batch, ...
            plan(q).log_file,plan(q).preference_directory, ...
            plan(q).temporary_directory);
        process_guard{q} = onCleanup(@() Terminate_Process_Tree(process{q}));
        started(q) = now;
        %% Index from numel rather than end: deleting the last entry from a
        %% 1x1 leaves a 1x0, so end+1 would resolve to row 2 and leave a
        %% zero in slot 1. Only reachable at MaxConcurrent = 1.
        running(numel(running)+1,1) = q; %#ok<AGROW>
        fprintf('  [%s] launch run %d  %s\n',datestr(now,'HH:MM:SS'), ...
            plan(q).run_id,plan(q).name);
    end
    pause(opt.PollSeconds);
    finished = false(size(running));
    for r = 1:numel(running)
        q = running(r);
        state = read_state(plan(q).status_file);
        alive = process_is_alive(process{q});
        if strcmp(state,'complete') && isfile(plan(q).result_file) && ~alive
            fprintf('  [%s] done run %d  %s\n',datestr(now,'HH:MM:SS'), ...
                plan(q).run_id,plan(q).name);
            cleanup_preference_directory(plan(q).runtime_directory);
            finished(r) = true;
        elseif strcmp(state,'failed')
            fprintf(2,'  [%s] FAILED run %d  %s; see %s\n', ...
                datestr(now,'HH:MM:SS'),plan(q).run_id, ...
                plan(q).name,plan(q).log_file);
            failed(end+1,1) = q; %#ok<AGROW>
            finished(r) = true;
            abort = true;
        elseif ~alive
            fprintf(2,'  [%s] FAILED run %d  %s (child exit %d); see %s\n', ...
                datestr(now,'HH:MM:SS'),plan(q).run_id,plan(q).name, ...
                process_exit_value(process{q}),plan(q).log_file);
            failed(end+1,1) = q; %#ok<AGROW>
            finished(r) = true;
        elseif (now-started(q))*24 > opt.MaxRuntimeHours
            Terminate_Process_Tree(process{q});
            fprintf(2,'  [%s] TIMEOUT run %d  %s after %.1f h\n', ...
                datestr(now,'HH:MM:SS'),plan(q).run_id, ...
                plan(q).name,opt.MaxRuntimeHours);
            failed(end+1,1) = q; %#ok<AGROW>
            finished(r) = true;
        end
    end
    running(finished) = [];
    running = running(:);   %% keep column orientation when it empties
end

if abort
    if ~isempty(running)
        fprintf(2,'  [%s] aborting: terminating %d run(s) still active\n', ...
            datestr(now,'HH:MM:SS'),numel(running));
        for r = 1:numel(running)
            Terminate_Process_Tree(process{running(r)});
        end
    end
    if ~isempty(queue)
        fprintf(2,'  not launched: %s\n',mat2str([plan(queue).run_id]));
    end
end

if ~isempty(failed)
    error('Run_All_12_Parallel:ChildFailed', ...
        'Run(s) %s failed or timed out.',mat2str([plan(failed).run_id]));
end
fprintf('All requested simulations completed.\n');
if ~isempty(pref_root) && isfolder(pref_root) && isempty(dir(fullfile(pref_root,'*')))
    rmdir(pref_root);
end

end

function tag = result_tag(R,opt)
tag = R.name;
if ~isempty(opt.Scenario)
    tag = sprintf('%s__%s',tag,char(opt.Scenario));
end
if ~isempty(opt.CalciteCase)
    tag = sprintf('%s__case%d',tag,opt.CalciteCase);
end
if R.use_subdomains && (~isempty(opt.GridNx) || ~isempty(opt.GridNy))
    nx = R.geom.main_grid_nx;
    ny = R.geom.main_grid_ny;
    if ~isempty(opt.GridNx); nx = opt.GridNx; end
    if ~isempty(opt.GridNy); ny = opt.GridNy; end
    tag = sprintf('%s__grid%dx%d',tag,nx,ny);
end
if ~isempty(opt.AppeloAlpha)
    alpha_tag = number_tag(opt.AppeloAlpha);
    tag = sprintf('%s__alpha%s',tag,alpha_tag);
end
if ~isempty(opt.AppeloLogKH)
    logkh_tag = arrayfun(@number_tag,opt.AppeloLogKH,'UniformOutput',false);
    tag = sprintf('%s__logkh%s',tag,strjoin(logkh_tag,'_'));
end
if ~isempty(opt.KNmaxd)
    tag = sprintf('%s__knmaxd%s',tag,number_tag(opt.KNmaxd));
end
if ~isempty(opt.Kbot)
    tag = sprintf('%s__kbot%s',tag,kbot_tag(opt.Kbot));
end
if ~isempty(opt.NDays); tag = sprintf('%s__%dd',tag,opt.NDays); end
if ~isempty(opt.StartDate)
    tag = sprintf('%s__from%s',tag,datestr(opt.StartDate,'yyyymmdd'));
end
end

function args = worker_args(run_id,status_file,extra)
parts = {sprintf('%d',run_id),sprintf('''%s''',strrep(status_file,'''',''''''))};
for k = 1:2:numel(extra)
    parts{end+1} = sprintf('''%s''',extra{k}); %#ok<AGROW>
    value = extra{k+1};
    if ischar(value) || (isstring(value) && isscalar(value))
        parts{end+1} = sprintf('''%s''',strrep(char(value),'''','''''')); %#ok<AGROW>
    else
        parts{end+1} = mat2str(value); %#ok<AGROW>
    end
end
args = strjoin(parts,',');
end

function process = launch_process(matlab_exe,batch,log_file,preference_directory,temporary_directory)
if ispc
    cmd_exe = fullfile(getenv('SystemRoot'),'System32','cmd.exe');
    args = sprintf('/d /s /c ""%s" -batch "%s" > "%s" 2>&1"', ...
        matlab_exe,batch,log_file);
    start_info = System.Diagnostics.ProcessStartInfo(cmd_exe,args);
    start_info.UseShellExecute = false;
    start_info.CreateNoWindow = true;
    if ~isempty(preference_directory)
        if start_info.EnvironmentVariables.ContainsKey('MATLAB_PREFDIR')
            start_info.EnvironmentVariables.Remove('MATLAB_PREFDIR');
        end
        start_info.EnvironmentVariables.Add( ...
            'MATLAB_PREFDIR',preference_directory);
    end
    if ~isempty(temporary_directory)
        for key = {'TEMP','TMP'}
            if start_info.EnvironmentVariables.ContainsKey(key{1})
                start_info.EnvironmentVariables.Remove(key{1});
            end
            start_info.EnvironmentVariables.Add(key{1},temporary_directory);
        end
    end
    process = System.Diagnostics.Process.Start(start_info);
    return
end
command = java.util.ArrayList;
command.add(matlab_exe);
command.add('-batch');
command.add(batch);
builder = java.lang.ProcessBuilder(command);
if ~isempty(preference_directory)
    builder.environment().put('MATLAB_PREFDIR',preference_directory);
end
if ~isempty(temporary_directory)
    builder.environment().put('TEMP',temporary_directory);
    builder.environment().put('TMP',temporary_directory);
end
builder.redirectErrorStream(true);
builder.redirectOutput(java.io.File(log_file));
process = builder.start();
end

function cleanup_preference_directory(preference_directory)
if ~isempty(preference_directory) && isfolder(preference_directory)
    rmdir(preference_directory,'s');
end
end

function tag = number_tag(value)
text = sprintf('%.12g',value);
if str2double(text) ~= value
    text = sprintf('%s_h%s',text,num2hex(value));
end
tag = strrep(text,'.','p');
tag = strrep(tag,'-','m');
tag = strrep(tag,'+','');
end

function tag = kbot_tag(value)
if isnan(value)
    tag = 'free';
else
    tag = number_tag(value);
end
end

function state = read_state(status_file)
state = '';
if ~isfile(status_file); return; end
txt = fileread(status_file);
tok = regexp(txt,'(?m)^state=([^\r\n]+)','tokens','once');
if ~isempty(tok); state = tok{1}; end
end

function alive = process_is_alive(process)
try
    if ispc && isa(process,'System.Diagnostics.Process')
        alive = ~process.HasExited;
    else
        alive = process.isAlive();
    end
catch
    alive = false;
end
end

function value = process_exit_value(process)
try
    if ispc && isa(process,'System.Diagnostics.Process')
        value = process.ExitCode;
    else
        value = process.exitValue();
    end
catch
    value = -1;
end
end

function validate_optional_integer(value,name)
if isempty(value); return; end
if ~(isnumeric(value) && isscalar(value) && isfinite(value) && ...
        value == fix(value) && value >= 1)
    error('Run_All_12_Parallel:BadOptions', ...
        '%s must be a positive integer.',name);
end
end

function validate_optional_nonnegative_scalar(value,name)
if isempty(value); return; end
if ~(isnumeric(value) && isscalar(value) && isfinite(value) && value >= 0)
    error('Run_All_12_Parallel:BadOptions', ...
        '%s must be a finite nonnegative scalar.',name);
end
end

function validate_optional_kbot(value)
if isempty(value); return; end
if ~(isnumeric(value) && isreal(value) && isscalar(value) && ...
        (isnan(value) || (isfinite(value) && value >= 0)))
    error('Run_All_12_Parallel:BadOptions', ...
        'Kbot must be NaN for free drainage or a finite nonnegative scalar.');
end
end

function validate_optional_nonnegative_integer(value,name)
if isempty(value); return; end
if ~(isnumeric(value) && isscalar(value) && isfinite(value) && ...
        value == fix(value) && value >= 0)
    error('Run_All_12_Parallel:BadOptions', ...
        '%s must be a nonnegative integer.',name);
end
end
