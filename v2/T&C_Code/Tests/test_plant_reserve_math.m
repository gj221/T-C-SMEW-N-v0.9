function test_plant_reserve_math
% Plant nutrient reserves must remain finite and tied to physical storage.

code_dir = fileparts(fileparts(mfilename('fullpath')));
addpath(code_dir);

assert(Plant_Reserve_Fraction(0,0,'N') == 0, ...
    'A zero reserve and zero capacity must give a zero allocation fraction.');
assert(abs(Plant_Reserve_Fraction(2,4,'N') - 0.5) < eps, ...
    'Reserve allocation must equal reserve divided by storage capacity.');
assert(Plant_Reserve_Fraction(6,4,'N') == 1, ...
    'Reserve allocation must remain bounded without changing stored mass.');

assert(Plant_Nutrient_Suppression(0,0,0,'N') == 0, ...
    'An empty plant must not produce a NaN nutrient suppression factor.');
expected = max(0,min(1,(2 - 0.8*4)/(0.6*10 + 0.2*4)));
assert(abs(Plant_Nutrient_Suppression(2,4,10,'N') - expected) < eps, ...
    'Nutrient suppression must retain the documented dimensionless equation.');

plant_exports = fileread(fullfile(code_dir,'Plant_Exports.m'));
assert(contains(plant_exports,'if sum(Btm1)==0 || Btm1(4)==0'), ...
    'Previous reserve capacity must retain the TC-SMEW zero-reserve-C boundary.');
assert(contains(plant_exports,'if sum(B)==0 || B(4)==0'), ...
    'Current reserve capacity must retain the TC-SMEW zero-reserve-C boundary.');

assert(Plant_Reserve_Fraction(1,0,'P') == 1, ...
    'Positive reserve at zero capacity must retain the saturated TC-SMEW limit.');
assert(Plant_Reserve_Fraction(1e-6,1e-18,'P') == 1, ...
    'Positive reserve at negligible capacity must retain the saturated TC-SMEW limit.');

expected = min(max(1/(0.6*2),0),1);
assert(abs(Plant_Nutrient_Suppression(1,0,2,'P') - expected) < eps, ...
    'Zero-capacity suppression must retain the TC-SMEW uptake equation.');
assert(Plant_Nutrient_Suppression(1,0,0,'P') == 1, ...
    'Positive reserve with no plant capacity must fully suppress further uptake.');

fprintf('PASS: plant reserve fractions are finite and capacity constrained.\n');
end
