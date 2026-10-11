function test_parallel_erw_configuration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parallel ERW execution contract.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_dir = fileparts(mfilename('fullpath'));
v2 = fileparts(fileparts(test_dir));
addpath(fullfile(v2,'run'),fullfile(v2,'T&C_Code'));

plan = Run_All_12_Parallel('DryRun',true,'Runs',1,'NDays',1, ...
    'ParallelWorkers',4);
assert(contains(plan.command,'''ParallelWorkers'',4'), ...
    'The launcher did not forward the ERW worker count.');

assert_error(@() Run_All_12_Parallel('DryRun',true, ...
    'ParallelWorkers',-1),'Run_All_12_Parallel:BadOptions');
assert_error(@() Run_One_Subfield(1,'SaveResults',false, ...
    'ParallelWorkers',1.5),'Run_One_Subfield:BadOptions');

driver_file = fullfile(v2,'T&C_Code','MAIN_FRAME_BG_ERW_2D.m');
driver_source = fileread(driver_file);
parallel_loops = regexp(driver_source,'parfor\s*\(c\s*=\s*1:N_col', ...
    'match');
assert(numel(parallel_loops) == 1, ...
    'Only the independent daily ERW column loop may use parfor.');
assert(contains(driver_source, ...
    'parfor (c = 1:N_col,parallel_workers)'), ...
    'The ERW parfor loop must be bounded by the requested worker count.');

erw_source = fileread(fullfile(v2,'T&C_Code','ERW_Daily_Step_2D.m'));
assert(~contains(erw_source,'persistent BiogeoPar_erw_handoff'), ...
    'The ERW worker path must not share mutable persistent state.');

fprintf('PASS: parallel ERW execution is isolated from lateral routing.\n');
end

function assert_error(f,expected_id)
try
    f();
catch ME
    assert(strcmp(ME.identifier,expected_id), ...
        'Expected %s, got %s: %s',expected_id,ME.identifier,ME.message);
    return
end
error('Expected error %s was not raised.',expected_id);
end
