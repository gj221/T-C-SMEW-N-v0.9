function test_checkpoint_omits_cleanup
% Checkpoint files must not serialize the caller's directory guard.

test_dir = fileparts(mfilename('fullpath'));
code_file = fullfile(test_dir,'..','MAIN_FRAME_BG_ERW_2D.m');
code = fileread(code_file);
assert(contains(code,"save(ERW_ckpt_tmp, '-regexp', '^(?!cleanup$).*', '-v7.3');"), ...
    'Checkpoint save does not exclude cleanup.');

tmp = [tempname '.mat'];
state_value = 17;
cleanup = onCleanup(@() []); %#ok<NASGU>
lastwarn('');
try
    save(tmp, '-regexp', '^(?!cleanup$).*', '-v7.3');
    saved = who('-file',tmp);
    [warn_msg,~] = lastwarn;
    assert(ismember('state_value',saved),'Numerical state was not saved.');
    assert(~ismember('cleanup',saved),'Cleanup guard entered checkpoint.');
    assert(isempty(warn_msg),'Checkpoint save emitted a warning: %s',warn_msg);
catch ME
    clear cleanup
    if isfile(tmp); delete(tmp); end
    rethrow(ME)
end
clear cleanup
delete(tmp)

fprintf('PASS: checkpoint retains state and omits the cleanup guard.\n');
end
