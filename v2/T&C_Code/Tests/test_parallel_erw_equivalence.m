function test_parallel_erw_equivalence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Serial and parallel ERW paths must produce identical treated states.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(ver('parallel')) || ~license('test','Distrib_Computing_Toolbox')
    error('test_parallel_erw_equivalence:MissingParallelToolbox', ...
        'Parallel Computing Toolbox is required for the Bouchet ERW path.');
end

test_dir = fileparts(mfilename('fullpath'));
v2 = fileparts(fileparts(test_dir));
addpath(fullfile(v2,'run'),fullfile(v2,'T&C_Code'));

parallel_dir = tempname;
serial_dir = tempname;
mkdir(parallel_dir);
mkdir(serial_dir);
old_checkpoint_minutes = getenv('TC_CKPT_MINUTES');
old_parallel_workers = getenv('TC_PARALLEL_WORKERS');
cleanup = onCleanup(@() restore_test_environment(parallel_dir,serial_dir, ...
    old_checkpoint_minutes,old_parallel_workers)); %#ok<NASGU>
setenv('TC_CKPT_MINUTES','0');
setenv('TC_PARALLEL_WORKERS','');

parallel_file = Run_One_Subfield(2,'StartDate','2016-11-10','NDays',1, ...
    'ParallelWorkers',6,'ResultDirectory',parallel_dir);
serial_file = Run_One_Subfield(2,'StartDate','2016-11-10','NDays',1, ...
    'ParallelWorkers',0,'ResultDirectory',serial_dir);

parallel = load(parallel_file);
serial = load(serial_file);
state_names = {'P','O','V','Tdp','PH_d','PH_pore_d','PH_bulk_d', ...
    'N2Oflx','LEAK_NO3','LEAK_NH4','NavlI','ERW','Domain'};
for k = 1:numel(state_names)
    name = state_names{k};
    assert(isfield(parallel,name) && isfield(serial,name), ...
        'Missing comparison state %s.',name);
    if strcmp(name,'ERW')
        %% ERW_var carries the fsolve options object built once in
        %% ERW_MOD_PARAM_excel.  It is solver configuration, not model
        %% state, and an optim.options.Fsolve copied through a thread pool
        %% does not satisfy isequaln even when every setting matches, so
        %% comparing it tests object identity rather than the science.
        %% The settings are compared explicitly instead, and the rest of
        %% the ERW state still has to agree exactly.
        for c = 1:numel(parallel.ERW)
            assert_same_solver_options(parallel.ERW{c},serial.ERW{c},c);
        end
        assert(isequaln(strip_solver_options(parallel.ERW), ...
                        strip_solver_options(serial.ERW)), ...
            'Parallel ERW changed %s.',name);
    else
        assert(isequaln(parallel.(name),serial.(name)), ...
            'Parallel ERW changed %s.',name);
    end
end

for c = 1:numel(parallel.ERW)
    assert(parallel.ERW{c}.TC_const.M_rock_add(1) == 5000 && ...
        any(parallel.ERW{c}.bg.M_min(1,:) > 0), ...
        'The exact-start basalt application was not present in column %d.',c);
end

assert(isempty(gcp('nocreate')), ...
    'A worker pool remained after the parallel run completed.');
fprintf('PASS: six-worker ERW matches the serial treated 2-D state exactly.\n');
end

function E = strip_solver_options(E)
for c = 1:numel(E)
    if isfield(E{c},'var') && isfield(E{c}.var,'fsolve_options')
        E{c}.var = rmfield(E{c}.var,'fsolve_options');
    end
end
end

function assert_same_solver_options(a,b,column)
if ~isfield(a,'var') || ~isfield(a.var,'fsolve_options')
    assert(~isfield(b,'var') || ~isfield(b.var,'fsolve_options'), ...
        'Column %d has solver options in only one path.',column);
    return
end
assert(isfield(b,'var') && isfield(b.var,'fsolve_options'), ...
    'Column %d has solver options in only one path.',column);
pa = a.var.fsolve_options; pb = b.var.fsolve_options;
assert(strcmp(class(pa),class(pb)), ...
    'Column %d solver options class differs: %s vs %s.', ...
    column,class(pa),class(pb));
for prop = {'Algorithm','Display','FunctionTolerance','StepTolerance', ...
        'MaxIterations','MaxFunctionEvaluations'}
    if isprop(pa,prop{1}) && isprop(pb,prop{1})
        assert(isequaln(pa.(prop{1}),pb.(prop{1})), ...
            'Column %d solver option %s differs.',column,prop{1});
    end
end
end

function restore_test_environment(parallel_dir,serial_dir, ...
        checkpoint_minutes,parallel_workers)
setenv('TC_CKPT_MINUTES',checkpoint_minutes);
setenv('TC_PARALLEL_WORKERS',parallel_workers);
pool = gcp('nocreate');
if ~isempty(pool); delete(pool); end
if isfolder(parallel_dir); rmdir(parallel_dir,'s'); end
if isfolder(serial_dir); rmdir(serial_dir,'s'); end
end
