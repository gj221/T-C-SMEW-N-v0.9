function test_denitrification_ph_response
% Gross denitrification and N2:N2O partition use the Energy Farm pH functions.

test_dir = fileparts(mfilename('fullpath'));
addpath(fileparts(test_dir));

pH = [3.5 4 5.5 7 7.8 8];
[gross_factor,n2_to_n2o] = Biogeo_Denitrification_PH_Response(pH);
pH_opt = 7.8;
expected_gross = min(1,max(0.001,(pH-4)/(pH_opt-4)));
expected_ratio = 0.0016*exp(1.006*pH);

assert(max(abs(gross_factor-expected_gross)) < 1e-12, ...
    'Gross denitrification pH response does not use the pH 7.8 breakpoint.');
assert(max(abs(n2_to_n2o-expected_ratio)) < 1e-12, ...
    'N2:N2O pH response does not match Blanc-Betes et al. (2021).');
assert(all(diff(gross_factor) >= 0) && all(diff(n2_to_n2o) >= 0) && ...
    n2_to_n2o(end) > n2_to_n2o(1), ...
    'Higher pH must not reduce gross denitrification or the N2:N2O ratio.');
assert(all(gross_factor >= 0.001 & gross_factor <= 1), ...
    'The gross denitrification pH factor must remain within its configured bounds.');
assert(all(gross_factor(pH <= 4) == 0.001), ...
    'Gross denitrification must use its lower floor at and below pH 4.');
assert(all(abs(gross_factor(pH >= pH_opt)-1) < 1e-12), ...
    'Gross denitrification must plateau at and above pH 7.8.');
assert(all(n2_to_n2o > 0), ...
    'The N2:N2O pH factor must remain positive.');
assert(all(n2_to_n2o(pH <= 4) < 0.1), ...
    'Blanc-Betes Eq. 9 carries no floor; the Eq. 6 floor is on the product.');
assert(Biogeo_Denitrification_Product_Ratio( ...
    n2_to_n2o(1),0.05,50,20) == 0.1, ...
    'The complete N2:N2O product must retain the published 0.1 floor.');

assert_throws(@() Biogeo_Denitrification_PH_Response(0), ...
    'Biogeo_Denitrification_PH_Response:InvalidPH');
assert_throws(@() Biogeo_Denitrification_PH_Response(NaN), ...
    'Biogeo_Denitrification_PH_Response:InvalidPH');

fprintf('PASS: denitrification pH functions use the literature-constrained pH 7.8 breakpoint.\n');
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
