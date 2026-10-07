function test_knmaxd_runtime_override
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run-local denitrification-rate override without editing shared source.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_dir = fileparts(mfilename('fullpath'));
code_dir = fileparts(test_dir);
addpath(code_dir);

baseline = Biogeochemistry_Parameter(1);
override = Biogeochemistry_Parameter(1,0.0064);
assert(abs(baseline.KNmaxd-0.0048) < eps, ...
    'The default denitrification rate changed.');
assert(abs(override.KNmaxd-0.0064) < eps, ...
    'The run-local denitrification rate was not applied.');

assert_error(@() Biogeochemistry_Parameter(1,-eps), ...
    'Biogeochemistry_Parameter:InvalidKNmaxd');
assert_error(@() Biogeochemistry_Parameter(1,[0.0064 0.0032]), ...
    'Biogeochemistry_Parameter:InvalidKNmaxd');

%% The ensemble's decomposition literals ride the same call.
replaced = Biogeochemistry_Parameter(1,[],struct('Vm',528,'mrb',0.0012));
assert(replaced.Vm == 528 && replaced.mrb == 0.0012 && ...
    replaced.Vpc == baseline.Vpc && replaced.KNmaxd == baseline.KNmaxd, ...
    'The run-local decomposition literals were not applied in isolation.');
assert_error(@() Biogeochemistry_Parameter(1,[],struct('KNmaxd',0.0064)), ...
    'Biogeochemistry_Parameter:InvalidOverride');
assert_error(@() Biogeochemistry_Parameter(1,[],struct('k_met_sur',1)), ...
    'Biogeochemistry_Parameter:UnknownOverride');
assert_error(@() Biogeochemistry_Parameter(1,[],struct('Vm',-1)), ...
    'Biogeochemistry_Parameter:InvalidOverride');
assert_error(@() Biogeochemistry_Parameter(1,[],struct('Vm',[1 2])), ...
    'Biogeochemistry_Parameter:InvalidOverride');

unit_source = fileread(fullfile(code_dir,'BIOGEO_UNIT.m'));
assert(contains(unit_source, ...
    'Biogeochemistry_Parameter(opt_cons_CUE,KNmaxd_override,biogeo_overrides)'), ...
    'BIOGEO_UNIT does not apply the run-local denitrification rate and literals.');

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
