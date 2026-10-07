function test_erw_application_event_alignment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Rock applications align with the timestamped ERW state grid.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dtw = 1/144;
n_steps = 3*144;
idx = ERW_Application_Indices([0 1 2 3],dtw,n_steps);
assert(isequal(idx,[1 145 289]), ...
    'Applications must map onto rows timestamped 0, 1 and 2 days.');
assert(ERW_Application_Indices(314,dtw,315*144) == 45217, ...
    'The first Kantola event must occupy row 45217 from 1 January.');

assert_error(@() ERW_Application_Indices(1+dtw/2,dtw,n_steps), ...
    'ERW:applicationGrid');

code_dir = fileparts(fileparts(mfilename('fullpath')));
daily_source = fileread(fullfile(code_dir,'ERW_Daily_Step_2D.m'));
frame_source = fileread(fullfile(code_dir,'MAIN_FRAME_BG_ERW_10mins.m'));
seed_source = fileread(fullfile(code_dir,'Parameters','ERW_MOD_PARAM_excel.m'));
assert(contains(daily_source,'ERW_Application_Indices') && ...
    contains(frame_source,'ERW_Application_Indices'), ...
    'The shared 1-D and 2-D paths must use the same event-index rule.');
assert(contains(seed_source, ...
    'M_rock(1) = max(M_rock_in,M_rock_add(1))') && ...
    contains(seed_source, ...
    'CaCO3(1) = CaCO3_in + M_rock_add_CaCO3(1)'), ...
    'Initial feedstock and a run-start application must enter the initial state once.');

fprintf('PASS: ERW applications align with the timestamped state grid.\n');
end

function assert_error(action,identifier)
try
    action();
catch ME
    assert(strcmp(ME.identifier,identifier), ...
        'Expected %s but received %s.',identifier,ME.identifier);
    return
end
error('Expected error %s was not raised.',identifier);
end
