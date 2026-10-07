function test_erw_activity_coefficient
%%% Davies activity coefficients and the pore-water ionic strength that
%%% enter the ERW saturation indices: reference values, boundary
%%% behaviour, and the recorded campaign state that motivated the change.

test_dir = fileparts(mfilename('fullpath'));
code_dir = fullfile(test_dir,'..');
addpath(code_dir);

% Davies (1962), A = 0.509 at 25 degC: hand-computed reference values.
ref = [0.001 1 0.9651; 0.001 2 0.8674; 0.001 3 0.7260; ...
       0.010 1 0.9021; 0.010 2 0.6622; 0.010 3 0.3956; ...
       0.100 1 0.7816; 0.100 2 0.3732; 0.100 3 0.1089];
for r = 1:size(ref,1)
    g = ERW_Activity_Coefficient(ref(r,2),ref(r,1));
    assert(abs(g-ref(r,3)) < 5e-4, ...
        'Davies coefficient for z = %d at I = %.3f must be %.4f (got %.4f).', ...
        ref(r,2),ref(r,1),ref(r,3),g);
end
assert(ERW_Activity_Coefficient(2,0) == 1, ...
    'Zero ionic strength must give gamma = 1 exactly.');
assert(ERW_Activity_Coefficient(0,0.05) == 1, ...
    'A neutral species must have gamma = 1 at any ionic strength.');
assert(ERW_Activity_Coefficient(-2,0.05) == ERW_Activity_Coefficient(2,0.05), ...
    'Only z^2 enters the Davies equation.');
g_path = arrayfun(@(I) ERW_Activity_Coefficient(2,I),linspace(1e-4,0.3,50));
assert(all(diff(g_path) < 0), ...
    'gamma must decrease with ionic strength over the range this model reaches.');
assert(ERW_Activity_Coefficient(3,0.05) < ERW_Activity_Coefficient(2,0.05) && ...
    ERW_Activity_Coefficient(2,0.05) < ERW_Activity_Coefficient(1,0.05), ...
    'Higher charge must be more suppressed at fixed ionic strength.');
did_reject = false;
try
    ERW_Activity_Coefficient(2,-0.01);
catch ME
    did_reject = strcmp(ME.identifier,'ERW:invalidIonicStrength');
end
assert(did_reject,'A negative ionic strength must be rejected, not clamped.');

% Ionic strength I = 0.5*sum(z^2 c) from the model's umol-conv/L inputs.
I_synth = ERW_Ionic_Strength(1000,200,100,50,0,1,0,200,1500,10,50,1e6);
assert(abs(I_synth-0.0033705) < 1e-7, ...
    'Ionic strength of the synthetic state must be 3.3705e-3 mol/L.');
I_minor = ERW_Ionic_Strength(1000,200,100,50,0,1,0,200,1500,10,50,1e6, ...
    30,4,3,2,1);
assert(abs(I_minor-I_synth-0.5*(30+4*4+3+2+1)/1e6) < 1e-12, ...
    'Phosphate, charged Al hydroxides and OH must enter ionic strength.');
did_reject = false;
try
    ERW_Ionic_Strength(-1,0,0,0,0,0,0,0,0,0,0,1e6);
catch ME
    did_reject = strcmp(ME.identifier,'ERW:negativeConcentration');
end
assert(did_reject,'A negative concentration must be rejected.');
% A value a rounding error below zero counts as zero: a 35-year run
% stopped on -4.5e-16 (5 October 2026). Beyond 1e-10 below zero it is
% still rejected.
I_zero = ERW_Ionic_Strength(1000,200,100,50,0,1,0,200,1500,10,50,1e6, ...
    30,4,3,0,1);
I_roundoff = ERW_Ionic_Strength(1000,200,100,50,0,1,0,200,1500,10,50,1e6, ...
    30,4,3,-4.54747e-16,1);
assert(I_roundoff == I_zero, ...
    'A concentration a rounding error below zero must count as zero.');
did_reject = false;
try
    ERW_Ionic_Strength(1000,200,100,50,0,1,0,200,1500,10,50,1e6,30,4,3,-2e-10,1);
catch ME
    did_reject = strcmp(ME.identifier,'ERW:negativeConcentration');
end
assert(did_reject,'A concentration beyond rounding distance below zero must be rejected.');
did_reject = false;
try
    ERW_Ionic_Strength(1,0,0,0,0,0,0,0,0,0,0,0);
catch ME
    did_reject = strcmp(ME.identifier,'ERW:invalidConvMol');
end
assert(did_reject,'A nonpositive conv_mol must be rejected.');

fprintf('PASS: Davies activity coefficients and ionic strength are dimensionally correct.\n');
end
