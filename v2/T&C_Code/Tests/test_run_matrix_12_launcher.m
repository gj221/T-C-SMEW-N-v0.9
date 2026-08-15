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

[default_output,default_plan] = evalc("Run_All_12_Parallel('DryRun',true)");
check(contains(default_output,'max concurrent: 2'), ...
    'launcher defaults to two workers for the 32 GB target machine');
check(contains(default_output,'max runtime per run: 168 h'), ...
    'launcher gives each worker a seven-day timeout buffer');
check(isequal([default_plan.run_id],[10:12 1:9]), ...
    'launcher queues all vertical one-dimensional runs first');

plan = Run_All_12_Parallel('DryRun',true,'Runs',[1 10 12], ...
    'NDays',2,'GridNx',1,'GridNy',1,'PollSeconds',1);
check(isequal([plan.run_id],[1 10 12]), ...
    'twelve-case dry run returns the selected scenarios');
check(all(contains({plan.result_file},fullfile('SIM','Subfield_12'))), ...
    'twelve-case results use the canonical output directory');
check(all(contains({plan.command},'Run_Subfield_Worker')), ...
    'every twelve-case command uses the status-writing worker');
check(all(contains({plan.command},'''NDays'',2')), ...
    'shortened duration is forwarded to every child');
check(all(contains({plan.command},'''GridNx'',1')) && ...
      all(contains({plan.command},'''GridNy'',1')), ...
    'grid overrides are forwarded to every child');
check(numel(unique({plan.result_file})) == 3 && ...
      numel(unique({plan.status_file})) == 3, ...
    'each child owns distinct result and status files');
full_plan = Run_All_12_Parallel('DryRun',true,'Runs',1);
check(~strcmp(plan(1).result_file,full_plan.result_file), ...
    'a smoke result cannot cause a full run to be skipped');
check(exist('Terminate_Process_Tree','file') == 2, ...
    'launcher provides process-tree cleanup for interrupted workers');
if ispc
    [test_process,child_pid] = start_test_process_tree();
    process_cleanup = onCleanup(@() Terminate_Process_Tree(test_process));
    check(~test_process.HasExited && child_pid > 0, ...
        'process-tree test starts a wrapper and child');
    Terminate_Process_Tree(test_process);
    pause(0.2);
    check(test_process.HasExited && ~windows_pid_exists(child_pid), ...
        'interrupt cleanup terminates the wrapper and child process');
    clear process_cleanup
end

assert_error(@() Run_All_12_Parallel('DryRun',true,'Runs',[12 13]), ...
    'Run_All_12_Parallel:BadRuns');
n_pass = n_pass + 1;
assert_error(@() Run_All_12_Parallel('DryRun',true,'MaxConcurrent',0), ...
    'Run_All_12_Parallel:BadMaxConcurrent');
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
assert_error(@() Run_One_Subfield(1.5,'SaveResults',false), ...
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

function [process,child_pid] = start_test_process_tree
system_root = getenv('SystemRoot');
cmd_exe = fullfile(system_root,'System32','cmd.exe');
ping_exe = fullfile(system_root,'System32','PING.EXE');
args = sprintf('/d /s /c ""%s" -n 30 127.0.0.1 >NUL"',ping_exe);
start_info = System.Diagnostics.ProcessStartInfo(cmd_exe,args);
start_info.UseShellExecute = false;
start_info.CreateNoWindow = true;
process = System.Diagnostics.Process.Start(start_info);
pause(0.3);
ps = sprintf(['powershell -NoProfile -Command "$p=(Get-CimInstance ' ...
    'Win32_Process -Filter ''ParentProcessId=%d'').ProcessId; $p"'],process.Id);
[status,out] = system(ps);
assert(status == 0,'Could not inspect the launcher test process tree.');
child_pid = sscanf(out,'%d',1);
if isempty(child_pid); child_pid = 0; end
end

function tf = windows_pid_exists(process_id)
[status,out] = system(sprintf( ...
    'tasklist /FI "PID eq %d" /FO CSV /NH',process_id));
tf = status == 0 && contains(out,sprintf('"%d"',process_id));
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
