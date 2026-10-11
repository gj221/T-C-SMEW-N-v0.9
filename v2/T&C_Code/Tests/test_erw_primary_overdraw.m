function test_erw_primary_overdraw
% Primary-mineral overdraw uses code-native mass units and inventory boundaries.

test_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(test_dir,'..'));

mineral = char('albite','calcite');
MM = [263;100.087]/1e6;
EW = [2e5 3e5;4e5 5e5];
dtw = 1/144;
reacted_expected = EW.*MM.'*dtw;
available = reacted_expected + [0.2 -0.1;0.4 0.3];

[max_overdraw,overdraw,reacted] = ERW_Primary_Overdraw( ...
    EW,dtw,MM,mineral,available);
assert(max(abs(reacted-reacted_expected),[],'all') < 1e-15, ...
    'EW columns did not retain the MM_min mineral ordering.');
assert(max(abs(overdraw-(reacted_expected-available)),[],'all') < 1e-15, ...
    'Overdraw did not use the stored post-reaction primary inventory.');
assert(abs(max_overdraw-0.1) < 1e-15, ...
    'The maximum positive primary overdraw is incorrect.');

secondary_CaCO3 = 1e6;
assert(max_overdraw > max([reacted(:)-(available(:)+secondary_CaCO3);0]), ...
    'Secondary carbonate inventory must not credit primary-mineral supply.');

[inactive_max,inactive_overdraw,inactive_reacted] = ERW_Primary_Overdraw( ...
    zeros(3,1),dtw,0,strings(1,0),zeros(3,1));
assert(inactive_max == 0 && all(inactive_overdraw == 0,'all') && ...
    all(inactive_reacted == 0,'all'), ...
    'The no-primary-mineral placeholder must remain inactive.');

assert_error(@() ERW_Primary_Overdraw(EW,dtw,MM,mineral,available(:,1)), ...
    'ERW:PrimaryOverdrawSize');
assert_error(@() ERW_Primary_Overdraw(EW,dtw,MM(1),mineral,available), ...
    'ERW:MineralMolarMassSize');
assert_error(@() ERW_Primary_Overdraw(-EW,dtw,MM,mineral,available), ...
    'ERW:PrimaryOverdrawFlux');
assert_error(@() ERW_Primary_Overdraw(EW,dtw,MM,mineral,-available), ...
    'ERW:PrimaryOverdrawInventory');

fprintf('PASS: primary overdraw preserves units, ordering and inventory boundaries.\n');
end

function assert_error(fn,id)
try
    fn();
catch ME
    assert(strcmp(ME.identifier,id), ...
        'Expected %s, received %s.',id,ME.identifier);
    return
end
error('Expected error %s was not raised.',id);
end
