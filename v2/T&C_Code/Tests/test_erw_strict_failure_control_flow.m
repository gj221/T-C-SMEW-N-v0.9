function test_erw_strict_failure_control_flow
%%% A hard ERW solve failure must not fall through to dead recovery code.

test_dir = fileparts(mfilename('fullpath'));
code_dir = fullfile(test_dir,'..');
frame_source = fileread(fullfile(code_dir,'MAIN_FRAME_BG_ERW_10mins.m'));
erw_source = fileread(fullfile(code_dir,'ERW_function.m'));

assert(contains(erw_source,"error('ERW:NoConverge'"), ...
    'ERW_function must stop before returning an unconverged state.');
assert(~contains(frame_source,'ERW_Should_Substep') && ...
    ~contains(frame_source,'ERW_Substep_State'), ...
    'The 1-D driver must not retain recovery code that is unreachable after a hard error.');

fprintf('PASS: strict ERW failure control flow has no unreachable recovery branch.\n');
end
