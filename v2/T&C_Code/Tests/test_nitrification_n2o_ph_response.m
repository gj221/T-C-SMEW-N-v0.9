function test_nitrification_n2o_ph_response
% Nitrification N2O yield follows the published pH response.

test_dir = fileparts(mfilename('fullpath'));
addpath(fileparts(test_dir));

pH = [5.8 6.1 6.2];
actual = Biogeo_Nitrification_N2O_Fraction(pH);
expected = 721.86*exp(-2.387*pH);
assert(max(abs(actual-expected)) < 1e-12, ...
    'Nitrification N2O yield does not match Val Martin et al. (2023).');
assert(all(diff(actual) < 0), ...
    'Nitrification N2O yield must decline as pH rises.');
assert(all(actual >= 2.5e-4 & actual <= 7.5e-4), ...
    'Cropland-pH yields must remain in the published 3e-4 to 7e-4 range.');
assert(Biogeo_Nitrification_N2O_Fraction(1) == 1, ...
    'The empirical yield must not exceed one mole N per mole N nitrified.');

assert_throws(@() Biogeo_Nitrification_N2O_Fraction(NaN), ...
    'Biogeo_Nitrification_N2O_Fraction:InvalidPH');
assert_throws(@() Biogeo_Nitrification_N2O_Fraction(14), ...
    'Biogeo_Nitrification_N2O_Fraction:InvalidPH');

fprintf('PASS: nitrification N2O yield follows its published pH response.\n');
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
