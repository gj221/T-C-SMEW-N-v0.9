function test_nitrification_rate_ph_response
% Gross nitrification uses the published Parton pH response.

test_dir = fileparts(mfilename('fullpath'));
addpath(fileparts(test_dir));

pH = [3.5 5 6.1 8 10];
factor = Biogeo_Nitrification_PH_Response(pH);
expected = min(1,max(0,0.56+atan(pi*0.45*(pH-5))/pi));

assert(max(abs(factor-expected)) < 1e-12, ...
    'Gross nitrification pH response does not match Parton et al. (1996).');
assert(all(diff(factor) >= 0) && all(factor >= 0 & factor <= 1), ...
    'Gross nitrification pH response must be monotonic and bounded.');

assert_throws(@() Biogeo_Nitrification_PH_Response(0), ...
    'Biogeo_Nitrification_PH_Response:InvalidPH');
assert_throws(@() Biogeo_Nitrification_PH_Response(NaN), ...
    'Biogeo_Nitrification_PH_Response:InvalidPH');
assert_throws(@() Biogeo_Nitrification_PH_Response(6+1i), ...
    'Biogeo_Nitrification_PH_Response:InvalidPH');

fprintf('PASS: gross nitrification follows the bounded Parton pH response.\n');
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
