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

%% A daily-store run (package y) saves a variable list: every variable but
%% the cleanup guard and the daily arrays outside the stores, whose last two
%% rows go into ckpt_daily_rows.
assert(contains(code,"ckpt_vars = setdiff(who,[{'cleanup'} DAILY_STATE_ONLY]);") && ...
    contains(code,"save(ERW_ckpt_tmp, ckpt_vars{:}, '-v7.3', '-nocompression');"), ...
    'The daily-store checkpoint does not exclude cleanup.');
tmp = [tempname '.mat'];
state_value = 17; %#ok<NASGU>
dropped_daily = zeros(5,2); %#ok<NASGU>
DAILY_STATE_ONLY = {'dropped_daily'};
ckpt_daily_rows = struct('dropped_daily',[1 2; 3 4]); %#ok<NASGU>
cleanup = onCleanup(@() []); %#ok<NASGU>
lastwarn('');
try
    ckpt_vars = setdiff(who,[{'cleanup'} DAILY_STATE_ONLY]);
    save(tmp, ckpt_vars{:}, '-v7.3', '-nocompression');
    saved = who('-file',tmp);
    [warn_msg,~] = lastwarn;
    assert(all(ismember({'state_value','ckpt_daily_rows'},saved)), ...
        'Numerical state was not saved.');
    assert(~any(ismember({'cleanup','dropped_daily'},saved)), ...
        'Cleanup guard or an excluded daily array entered the checkpoint.');
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
