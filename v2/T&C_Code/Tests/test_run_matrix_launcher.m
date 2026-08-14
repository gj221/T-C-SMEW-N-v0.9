function test_run_matrix_launcher
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Matrix definition and process-launcher contract.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n=== test_run_matrix_launcher ===\n');
n_pass = 0;

RM = Run_Matrix_9();
check(numel(RM) == 9, 'matrix contains nine entries');
check(isequal({RM.treatment}, repmat({'control','basalt','calcite'},1,3)), ...
    'treatment order is control, basalt, calcite within each layout');
check(isequal([RM.case_number], repmat([5 3 7],1,3)), ...
    'workbook cases are control 5, basalt 3, calcite 7');
check(strcmp(RM(8).margin_pft,'PFT_Prairie_Margin'), ...
    'run 8 assigns the prairie margin PFT');
check(~RM(1).apply_main && ~RM(1).apply_margin && ...
      ~RM(4).apply_main && ~RM(4).apply_margin, ...
    'control placement flags are off');
check(strcmp(RM(1).name,'control__whole_field') && ...
      strcmp(RM(4).name,'control__margin_crop'), ...
    'duplicate controls retain distinct matrix names');

plan = Run_All_9_Parallel('DryRun',true,'Runs',[1 8], ...
    'NDays',3,'GridNx',1,'GridNy',1,'PollSeconds',1);
check(isequal([plan.run_id],[1 8]), 'dry run returns the selected runs');
check(all(contains({plan.command},'Run_Subfield_Worker')), ...
    'every child command uses the status-writing worker');
check(all(contains({plan.command},'''NDays'',3')), ...
    'NDays is forwarded to every child');
check(all(contains({plan.command},'''GridNx'',1')) && ...
      all(contains({plan.command},'''GridNy'',1')), ...
    'grid overrides are forwarded to every child');
check(numel(unique({plan.result_file})) == 2 && ...
      numel(unique({plan.status_file})) == 2, ...
    'each child owns distinct result and status files');
check(all(contains({plan.result_file},'__3d.mat')) && ...
      all(contains({plan.status_file},'__3d.status')), ...
    'shortened runs have smoke-specific result and status names');
full_plan = Run_All_9_Parallel('DryRun',true,'Runs',1);
check(~strcmp(plan(1).result_file,full_plan.result_file), ...
    'a smoke result cannot cause a full run to be skipped');

assert_error(@() Run_All_9_Parallel('DryRun',true,'MaxConcurrent',0), ...
    'Run_All_9_Parallel:BadMaxConcurrent');
n_pass = n_pass + 1;
assert_error(@() Run_All_9_Parallel('DryRun',true,'Runs',[0 1]), ...
    'Run_All_9_Parallel:BadRuns');
n_pass = n_pass + 1;
assert_error(@() Run_One_Subfield(13,'SaveResults',false), ...
    'Run_One_Subfield:BadId');
n_pass = n_pass + 1;
assert_error(@() Run_One_Subfield(1.5,'SaveResults',false), ...
    'Run_One_Subfield:BadId');
n_pass = n_pass + 1;

fprintf('  ---------------------------------------\n');
fprintf('  %d passed, 0 failed\n', n_pass);

    function check(ok, label)
        assert(ok, 'test_run_matrix_launcher:Failed', '%s', label);
        n_pass = n_pass + 1;
        fprintf('  PASS  %s\n', label);
    end

end

function assert_error(f, expected_id)
try
    f();
catch ME
    assert(strcmp(ME.identifier,expected_id), ...
        'Expected %s, got %s: %s', expected_id, ME.identifier, ME.message);
    fprintf('  PASS  rejects invalid input with %s\n', expected_id);
    return
end
error('Expected error %s was not raised.', expected_id);
end
