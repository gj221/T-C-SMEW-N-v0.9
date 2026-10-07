function test_current_12_run_mineralogy_zero_p
%%% Every treatment selected by Run_Matrix_12 has zero mineral P today.

test_dir = fileparts(mfilename('fullpath'));
code_dir = fileparts(test_dir);
v2_dir = fileparts(code_dir);
addpath(fullfile(v2_dir,'run'),code_dir, ...
    fullfile(code_dir,'Inputs'),fullfile(code_dir,'Parameters'));

run_matrix = Run_Matrix_12;
expected_cases = repmat([5 4 7],1,4);
assert(isequal([run_matrix.case_number],expected_cases), ...
    'The twelve-run treatment-to-workbook mapping changed.');

result_dir = tempname;
mkdir(result_dir);
old_ckpt = getenv('TC_CKPT_MINUTES');
old_workers = getenv('TC_PARALLEL_WORKERS');
cleanup = onCleanup(@() restore(result_dir,old_ckpt,old_workers)); %#ok<NASGU>
setenv('TC_CKPT_MINUTES','0');
setenv('TC_PARALLEL_WORKERS','0');

expected_rows = [1 6 1];
for run_id = 10:12
    result_file = Run_One_Subfield(run_id, ...
        'StartDate',datenum(2016,11,8),'NDays',4, ...
        'AppeloAlpha',1.3,'ParallelWorkers',0, ...
        'SaveResults',true,'ResultDirectory',result_dir);
    S = load(result_file,'ERW');
    E = S.ERW{1};
    k = run_id-9;
    assert(size(E.const.min_st,1) == expected_rows(k), ...
        'Run %d did not initialize the expected feedstock minerals.',run_id);
    assert(all(E.const.min_st(:,7) == 0), ...
        'Run %d currently selects a P-bearing mineral.',run_id);
    assert(isfield(E.weath,'P_weathering_release') && ...
        all(E.weath.P_weathering_release(:) == 0), ...
        'Run %d released phosphorus from its current feedstock.',run_id);
end

fprintf('PASS: all twelve mapped runs currently use P-free mineralogy.\n');
end

function restore(result_dir,old_ckpt,old_workers)
setenv('TC_CKPT_MINUTES',old_ckpt);
setenv('TC_PARALLEL_WORKERS',old_workers);
if isfolder(result_dir); rmdir(result_dir,'s'); end
end
