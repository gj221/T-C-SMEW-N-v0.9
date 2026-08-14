function test_run_matrix_12_launcher
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Twelve-scenario matrix, launcher and plotting contract.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n=== test_run_matrix_12_launcher ===\n');
n_pass = 0;

RM = Run_Matrix_12();
check(numel(RM) == 12, 'matrix contains twelve entries');
check(isequal({RM.treatment},repmat({'control','basalt','calcite'},1,4)), ...
    'each scenario group contains control, basalt and calcite');
check(isequal({RM.layout},repelem( ...
    {'whole_field','margin_crop','margin_prairie','vertical_1d'},3)), ...
    'matrix contains four ordered scenario groups');
check(isequal([RM.case_number],repmat([5 3 7],1,4)), ...
    'all groups use the same three workbook treatment classes');
check(all([RM(1:9).use_subdomains]) && ...
      ~any([RM(10:12).use_subdomains]), ...
    'only runs 10 to 12 omit field subdomains');
check(~RM(10).apply_main && RM(11).apply_main && RM(12).apply_main && ...
      ~any([RM(10:12).apply_margin]), ...
    'one-dimensional treatments apply across the single field column');
check(all(arrayfun(@(r) r.geom.field_side_m^2 == 3600,RM(10:12))), ...
    'one-dimensional scenarios retain the 3600 m2 plot area');
check(all([RM.biogeo_depth_mm] == 2500), ...
    'all twelve scenarios use the common 2500 mm active depth');

RM9 = Run_Matrix_9();
check(numel(RM9) == 9 && isequal({RM9.name},{RM(1:9).name}), ...
    'the nine-case matrix remains a compatible subset');

default_output = evalc("Run_All_12_Parallel('DryRun',true,'Runs',1)");
check(contains(default_output,'max concurrent: 2'), ...
    'launcher defaults to two workers for the 32 GB target machine');

plan = Run_All_12_Parallel('DryRun',true,'Runs',[1 10 12], ...
    'NDays',2,'PollSeconds',1);
check(isequal([plan.run_id],[1 10 12]), ...
    'twelve-case dry run returns the selected scenarios');
check(all(contains({plan.result_file},fullfile('SIM','Subfield_12'))), ...
    'twelve-case results use the canonical output directory');
check(all(contains({plan.command},'Run_Subfield_Worker')), ...
    'every twelve-case command uses the status-writing worker');
check(all(contains({plan.command},'''NDays'',2')), ...
    'shortened duration is forwarded to every child');

plan9 = Run_All_9_Parallel('DryRun',true,'Runs',[1 9]);
check(isequal([plan9.run_id],[1 9]) && ...
      all(contains({plan9.result_file},fullfile('SIM','Subfield_12'))), ...
    'nine-case launcher remains a compatible subset launcher');

assert_error(@() Run_All_12_Parallel('DryRun',true,'Runs',[12 13]), ...
    'Run_All_12_Parallel:BadRuns');
n_pass = n_pass + 1;
assert_error(@() Run_All_9_Parallel('DryRun',true,'Runs',10), ...
    'Run_All_9_Parallel:BadRuns');
n_pass = n_pass + 1;
assert_error(@() Run_All_9_Parallel('Runs'), ...
    'Run_All_9_Parallel:BadOptions');
n_pass = n_pass + 1;
assert_error(@() Run_All_12_Parallel('DryRun',1), ...
    'Run_All_12_Parallel:BadOptions');
n_pass = n_pass + 1;
assert_error(@() Run_All_12_Parallel('DryRun',true,'Overwrite',1), ...
    'Run_All_12_Parallel:BadOptions');
n_pass = n_pass + 1;
assert_error(@() Run_One_Subfield(13,'SaveResults',false), ...
    'Run_One_Subfield:BadId');
n_pass = n_pass + 1;

fprintf('  ---------------------------------------\n');
fprintf('  %d passed, 0 failed\n',n_pass);

    function check(ok,label)
        assert(ok,'test_run_matrix_12_launcher:Failed','%s',label);
        n_pass = n_pass + 1;
        fprintf('  PASS  %s\n',label);
    end

end

function assert_error(f,expected_id)
try
    f();
catch ME
    assert(strcmp(ME.identifier,expected_id), ...
        'Expected %s, got %s: %s',expected_id,ME.identifier,ME.message);
    fprintf('  PASS  rejects invalid input with %s\n',expected_id);
    return
end
error('Expected error %s was not raised.',expected_id);
end
