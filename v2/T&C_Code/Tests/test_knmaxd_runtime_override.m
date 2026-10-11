function test_knmaxd_runtime_override
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run-local denitrification-rate override without editing shared source.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_dir = fileparts(mfilename('fullpath'));
code_dir = fileparts(test_dir);
addpath(code_dir);

baseline = Biogeochemistry_Parameter(1);
override = Biogeochemistry_Parameter(1,0.0064);
assert(abs(baseline.KNmaxd-0.0271) < eps, ...
    'The default denitrification rate changed.');
assert(abs(override.KNmaxd-0.0064) < eps, ...
    'The run-local denitrification rate was not applied.');

assert_error(@() Biogeochemistry_Parameter(1,-eps), ...
    'Biogeochemistry_Parameter:InvalidKNmaxd');
assert_error(@() Biogeochemistry_Parameter(1,[0.0064 0.0032]), ...
    'Biogeochemistry_Parameter:InvalidKNmaxd');

unit_source = fileread(fullfile(code_dir,'BIOGEO_UNIT.m'));
assert(contains(unit_source, ...
    'Biogeochemistry_Parameter(opt_cons_CUE,KNmaxd_override)'), ...
    'BIOGEO_UNIT does not apply the run-local denitrification rate.');

fprintf('PASS: run-local KNmaxd override is validated and isolated.\n');
end

function assert_error(fun,identifier)
try
    fun();
catch ME
    assert(strcmp(ME.identifier,identifier), ...
        'Expected %s, got %s.',identifier,ME.identifier);
    return
end
error('Expected error %s was not raised.',identifier);
end
