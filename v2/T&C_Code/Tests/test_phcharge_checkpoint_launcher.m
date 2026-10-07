function test_phcharge_checkpoint_launcher
test_dir = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(fileparts(test_dir)));
launcher = fullfile(root,'run_phcharge.m');
assert(isfile(launcher),'The chemistry rerun launcher is missing.');
source = fileread(launcher);
assert(contains(source,"setenv('TC_CKPT_MINUTES','30')"), ...
    'Chemistry reruns must write rolling checkpoints.');
assert(contains(source,"fullfile(root,'results')"), ...
    'Chemistry reruns must write to their own result directory.');
assert(~contains(source,'C:/t/omfix'), ...
    'Chemistry reruns must not overwrite active Claude results.');
fprintf('PASS: chemistry reruns enable isolated rolling checkpoints.\n');
end
