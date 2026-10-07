function test_erw_initial_boundary_defaults
% Initialization must use explicit or physically derived boundaries.

test_dir = fileparts(mfilename('fullpath'));
addpath(fileparts(test_dir));

S = struct('seed_water_loss_l_m2_d',2.5);
assert(ERW_Seed_Water_Loss(S,struct()) == 2.5, ...
    'Explicit annual water loss must take precedence.');

forcing = struct('Tw',[0.001 0.002],'Lw',[0.002 0.004]);
assert(abs(ERW_Seed_Water_Loss(struct(),forcing)-4.5) < 1e-12, ...
    'Forcing-derived water loss must convert m/d to L m-2 d-1.');
assert(isnan(ERW_Seed_Water_Loss(struct(),struct())), ...
    'Missing water loss must remain unavailable, not use an arbitrary fallback.');

assert(ERW_AlOH3_LogK(false) == 8.04, ...
    'Crystalline gibbsite must use the Lindsay thermodynamic constant.');
assert(ERW_AlOH3_LogK(true) == 10.8, ...
    'The amorphous phase must use its explicit solubility constant.');

[cec,diag] = ERW_Resolve_CEC(18.5);
assert(cec == 18.5 && ~diag.measured && ~diag.estimated && ...
    strcmp(diag.source,"explicit_prescribed_CEC"), ...
    'A prescribed CEC must pass through unchanged without being called measured.');
assert_throws(@() ERW_Resolve_CEC(NaN),'ERW:InvalidCEC');

fprintf('PASS: ERW initialization has explicit water, Al-phase, and CEC boundaries.\n');
end

function assert_throws(fun,identifier)
threw = false;
try
    fun();
catch ME
    threw = strcmp(ME.identifier,identifier);
end
assert(threw,'Expected error %s.',identifier);
end
