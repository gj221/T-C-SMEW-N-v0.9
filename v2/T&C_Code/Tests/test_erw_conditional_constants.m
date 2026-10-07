function test_erw_conditional_constants
%%% The conditional constants close every mass-action law on an independent
%%% activity route, from I = 0 to the calcite-amended I = 0.126.
%%%
%%% Concentrations are speciated with the conditional constants, converted
%%% to activities with the Davies coefficients, and the THERMODYNAMIC
%%% constants are then recovered from the activity quotients.  A factor on
%%% the wrong side, a squared divalent coefficient or a lost conv_mol power
%%% all break the recovery by far more than the tolerance.

test_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(test_dir,'..'));

c = 1e6;                                   % conv_mol, umol/L per mol/L
k = [10^-6.35*c, 10^-10.33*c, 1e-14*c^2, 10^-5*c, 10^-5.1*c, 10^-6.7*c, 10^-6.2*c];
logKso = 10.8;

for I = [0 1e-4 0.03 0.126]
    [kc,logKc,g] = ERW_Conditional_Constants(k,logKso,I);
    assert(isequal(g,ERW_Activity_Coefficient(1:3,I)), ...
        'The returned coefficients must be the Davies values at I = %g.',I);
    assert(abs(g(2)-g(1)^4) <= 8*eps && abs(g(3)-g(1)^9) <= 8*eps, ...
        'Davies requires g2 = g1^4 and g3 = g1^9 at I = %g.',I);

    %% carbonate and water, concentrations in umol/L
    H = 2.0; CO2 = 300;
    HCO3 = kc(1)*CO2/H;
    CO3 = kc(2)*HCO3/H;
    OH = kc(3)/H;
    aH = g(1)*H/c; aHCO3 = g(1)*HCO3/c; aCO3 = g(2)*CO3/c;
    close_rel(aH*aHCO3/(CO2/c), k(1)/c, 'K1', I);
    close_rel(aH*aCO3/aHCO3, k(2)/c, 'K2', I);
    close_rel(aH*(g(1)*OH/c), k(3)/c^2, 'Kw', I);

    %% Al hydrolysis chain, charges 3, 2, 1, 0, -1
    Al3 = 5e-3;
    AlOH = kc(4)*Al3/H;
    AlOH2 = kc(5)*AlOH/H;
    AlOH3 = kc(6)*AlOH2/H;
    AlOH4 = kc(7)*AlOH3/H;
    close_rel((g(2)*AlOH/c)*aH/(g(3)*Al3/c), k(4)/c, 'K1_Al', I);
    close_rel((g(1)*AlOH2/c)*aH/(g(2)*AlOH/c), k(5)/c, 'K2_Al', I);
    close_rel((AlOH3/c)*aH/(g(1)*AlOH2/c), k(6)/c, 'K3_Al', I);
    close_rel((g(1)*AlOH4/c)*aH/(AlOH3/c), k(7)/c, 'K4_Al', I);

    %% Al(OH)3 solubility: the concentration the conditional *Kso admits
    Al3_sat = 10^logKc*(H/c)^3;                      % mol/L
    close_rel(g(3)*Al3_sat/aH^3, 10^logKso, '*Kso', I);

    if I == 0
        assert(isequal(kc,k) && logKc == logKso, ...
            'At I = 0 the conditional constants must equal the thermodynamic ones.');
    else
        assert(kc(1) > k(1) && kc(2) > k(2) && kc(3) > k(3) && logKc > logKso, ...
            'Nonideality must raise the concentration quotients at I = %g.',I);
    end
end

expect_error(@() ERW_Conditional_Constants(k,logKso,-1),'ERW:invalidIonicStrength');
expect_error(@() ERW_Conditional_Constants(k(1:6),logKso,0.01),'ERW:invalidParam');
expect_error(@() ERW_Conditional_Constants(k,NaN,0.01),'ERW:invalidParam');
[kc_empty,logKc_empty] = ERW_Conditional_Constants(k,[],0.03);
assert(isempty(logKc_empty) && isequal(kc_empty,ERW_Conditional_Constants(k,logKso,0.03)), ...
    'An empty logKso must leave the other constants unchanged and return empty.');

fprintf('PASS: conditional constants close every mass-action law on the activity basis.\n');
end

function close_rel(actual,expected,label,I)
assert(abs(actual/expected-1) <= 1e-12, ...
    '%s is not recovered at I = %g: %.15g against %.15g.',label,I,actual,expected);
end

function expect_error(fn,identifier)
threw = false;
try
    fn();
catch err
    threw = strcmp(err.identifier,identifier);
end
assert(threw,'Expected error %s.',identifier);
end
