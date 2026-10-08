function test_2d_parfor_equivalence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Exact serial/thread-parfor comparison on the 12-column grid, and the
%% daily stores of a multi-square run against its full records.
%%
%% Three days of run 2 (basalt on every square) from 9 November 2016,
%% across the rock application of 10 November, so the ERW chemistry runs
%% in daily blocks on amended columns and the 50-hour hourly windows wrap:
%%   serial    daily stores, ParallelWorkers 0
%%   parallel  daily stores, ParallelWorkers 4: the daily vegetation and
%%             biogeochemistry, the hourly physics and the daily ERW
%%             chemistry each run square by square on the pool
%%   full      full records (TC_FULL_RECORDS=1), serial
%% parallel must equal serial in every saved variable but the run's timing
%% and worker count, and serial must reproduce full's records exactly
%% (Daily_Store_Check): every daily store, hourly field mean, ERW day
%% summary and the final state.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_dir = fileparts(mfilename('fullpath'));
code_dir = fileparts(test_dir);
v2_dir = fileparts(code_dir);
run_dir = fullfile(v2_dir,'run');
addpath(code_dir,run_dir);

if ~license('test','Distrib_Computing_Toolbox')
    error('test_2d_parfor_equivalence:ParallelToolboxUnavailable', ...
        'Parallel Computing Toolbox is required for this test.');
end
root = tempname(tempdir);
serial_dir = fullfile(root,'serial');
parallel_dir = fullfile(root,'parallel');
full_dir = fullfile(root,'full');
mkdir(serial_dir); mkdir(parallel_dir); mkdir(full_dir);
old_workers = getenv('TC_PARALLEL_WORKERS');
old_ckpt = getenv('TC_CKPT_MINUTES');
old_force = getenv('TC_CKPT_FORCE_RESUME');
old_records = getenv('TC_FULL_RECORDS');
cleanup = onCleanup(@() cleanup_test(root,old_workers,old_ckpt,old_force,old_records)); %#ok<NASGU>
setenv('TC_CKPT_MINUTES','0');
setenv('TC_CKPT_FORCE_RESUME','');
setenv('TC_FULL_RECORDS','');

args = {'NDays',3,'StartDate',datenum(2016,11,9), ...
    'GridNx',2,'GridNy',2,'AppeloAlpha',1.3,'SaveResults',true};
t_serial = tic;
serial_file = Run_One_Subfield(2,args{:},'ParallelWorkers',0, ...
    'ResultDirectory',serial_dir);
serial_seconds = toc(t_serial);
t_parallel = tic;
parallel_file = Run_One_Subfield(2,args{:},'ParallelWorkers',4, ...
    'ResultDirectory',parallel_dir);
parallel_seconds = toc(t_parallel);
setenv('TC_FULL_RECORDS','1');
full_file = Run_One_Subfield(2,args{:},'ParallelWorkers',0, ...
    'ResultDirectory',full_dir);
setenv('TC_FULL_RECORDS','');

serial = load(serial_file);
parallel = load(parallel_file);
full = load(full_file);
assert(~serial.FULL_RECORDS && ~parallel.FULL_RECORDS && full.FULL_RECORDS, ...
    'The runs did not keep the records asked for.');
assert(serial.NN > serial.HOURLY_WINDOW, ...
    'The run is too short for the hourly windows to wrap.');
amended = cellfun(@(E) any(E.TC_const.M_rock_add(2:end) > 0),full.ERW);
assert(any(amended),'No square received rock inside the run.');

%% What may differ: the run's timing, its worker count and the identities
%% and paths that carry them.
timing = {'Computational_Time','run_provenance','ERW_ckpt_last_toc', ...
    'ERW_ckpt_runtime_identity','ERW_execution_identity','ERW_ckpt_fingerprint', ...
    'ERW_ckpt_file','ERW_ckpt_tmp','ckpt_dir','tmp_result','result_file', ...
    'parallel_workers','opt'};
names = union(fieldnames(serial),fieldnames(parallel));
for k = 1:numel(names)
    name = names{k};
    assert(isfield(serial,name) && isfield(parallel,name), ...
        'Only one of the serial and parfor results holds %s.',name);
    if ~any(strcmp(name,timing))
        assert(same_value(serial.(name),parallel.(name)), ...
            'Serial/parfor mismatch in %s.',name);
    end
end

%% The daily-store run against the same run with full records.
records = [timing {'FULL_RECORDS','NH'}];
[bad,report] = Daily_Store_Check(full,serial,records);
assert(isempty(bad), ...
    'The daily stores do not reproduce the full records:\n  %s',strjoin(bad,'\n  '));
assert(report.n_erw_squares == 12 && report.n_erw_days == 12*(serial.j) && ...
    report.n_daily_steps == serial.j - 1, ...
    'The check covered %d ERW squares, %d ERW days and %d daily steps.', ...
    report.n_erw_squares,report.n_erw_days,report.n_daily_steps);

fprintf(['PASS: 12-column serial/parfor outputs are exactly equal and the ' ...
    'daily stores reproduce the full records (%d variables, %d ERW days; ' ...
    'serial %.3f s, parfor %.3f s, speedup %.2fx)\n'], ...
    report.n_same,report.n_erw_days,serial_seconds,parallel_seconds, ...
    serial_seconds/parallel_seconds);

end

function tf = same_value(a,b)
%% isequaln, except that function handles compare by their text and the
%% values they captured (two loaded copies of one anonymous function are
%% never isequaln).
if isequaln(a,b); tf = true; return; end
tf = false;
if isa(a,'function_handle') && isa(b,'function_handle')
    if ~strcmp(func2str(a),func2str(b)); return; end
    fa = functions(a); fb = functions(b);
    if isfield(fa,'workspace') ~= isfield(fb,'workspace'); return; end
    tf = ~isfield(fa,'workspace') || same_value(fa.workspace,fb.workspace);
    return
end
if ~strcmp(class(a),class(b)) || ~isequal(size(a),size(b)); return; end
if isstruct(a)
    names = fieldnames(a);
    if ~isempty(setxor(names,fieldnames(b))); return; end
    for k = 1:numel(a)
        for q = 1:numel(names)
            if ~same_value(a(k).(names{q}),b(k).(names{q})); return; end
        end
    end
    tf = true;
elseif iscell(a)
    for k = 1:numel(a)
        if ~same_value(a{k},b{k}); return; end
    end
    tf = true;
end
end

function cleanup_test(root,workers,ckpt,force_resume,records)
setenv('TC_PARALLEL_WORKERS',workers);
setenv('TC_CKPT_MINUTES',ckpt);
setenv('TC_CKPT_FORCE_RESUME',force_resume);
setenv('TC_FULL_RECORDS',records);
pool = gcp('nocreate');
if ~isempty(pool); delete(pool); end
root = char(java.io.File(root).getCanonicalPath());
tmp = char(java.io.File(tempdir).getCanonicalPath());
if isfolder(root) && startsWith(root,[tmp filesep])
    rmdir(root,'s');
end
end
