function test_2d_parfor_equivalence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% One-day exact serial/thread-parfor comparison on the 12-column grid.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_dir = fileparts(mfilename('fullpath'));
code_dir = fileparts(test_dir);
v2_dir = fileparts(code_dir);
run_dir = fullfile(v2_dir,'run');
addpath(code_dir,run_dir);

root = tempname(tempdir);
serial_dir = fullfile(root,'serial');
parallel_dir = fullfile(root,'parallel');
mkdir(serial_dir); mkdir(parallel_dir);
old_workers = getenv('TC_PARALLEL_WORKERS');
old_ckpt = getenv('TC_CKPT_MINUTES');
old_force = getenv('TC_CKPT_FORCE_RESUME');
cleanup = onCleanup(@() cleanup_test(root,old_workers,old_ckpt,old_force)); %#ok<NASGU>
setenv('TC_CKPT_MINUTES','0');
setenv('TC_CKPT_FORCE_RESUME','');

args = {'NDays',1,'StartDate',datenum(2016,1,1), ...
    'GridNx',2,'GridNy',2,'AppeloAlpha',1.3,'SaveResults',true};
t_serial = tic;
serial_file = Run_One_Subfield(1,args{:},'ParallelWorkers',0, ...
    'ResultDirectory',serial_dir);
serial_seconds = toc(t_serial);

if ~license('test','Distrib_Computing_Toolbox')
    error('test_2d_parfor_equivalence:ParallelToolboxUnavailable', ...
        'Parallel Computing Toolbox is required for this test.');
end
t_parallel = tic;
parallel_file = Run_One_Subfield(1,args{:},'ParallelWorkers',4, ...
    'ResultDirectory',parallel_dir);
parallel_seconds = toc(t_parallel);

variables = {'Domain','P','O','Tdp','ERW','VOL','NITflx','DENITflx', ...
    'N2Oflx_nit','N2flx','N2Oflx','PH_bulk_d','PH_pore_d'};
serial = load(serial_file,variables{:});
parallel = load(parallel_file,variables{:});
for k = 1:numel(variables)
    name = variables{k};
    assert(isequaln(serial.(name),parallel.(name)), ...
        'Serial/parfor mismatch in %s.',name);
end

fprintf(['PASS: 12-column serial/parfor outputs are exactly equal ' ...
    '(serial %.3f s, parfor %.3f s, speedup %.2fx)\n'], ...
    serial_seconds,parallel_seconds,serial_seconds/parallel_seconds);

end

function cleanup_test(root,workers,ckpt,force_resume)
setenv('TC_PARALLEL_WORKERS',workers);
setenv('TC_CKPT_MINUTES',ckpt);
setenv('TC_CKPT_FORCE_RESUME',force_resume);
pool = gcp('nocreate');
if ~isempty(pool); delete(pool); end
root = char(java.io.File(root).getCanonicalPath());
tmp = char(java.io.File(tempdir).getCanonicalPath());
if isfolder(root) && startsWith(root,[tmp filesep])
    rmdir(root,'s');
end
end
