function test_pH_area_weighted
% Domain pH is calculated from hydrogen activity, never arithmetic pH.

here = fileparts(mfilename('fullpath'));
addpath(fileparts(here));

pH = pH_Area_Weighted([4 8],[1 1]);
assert(abs(pH-(-log10(0.5*(10^-4+10^-8)))) < 1e-12, ...
    'Area-weighted pH must equal the hydrogen-activity mean.');
assert(abs(pH-mean([4 8])) > 1, ...
    'Area-weighted pH must not collapse to an arithmetic pH mean.');

series = pH_Area_Weighted([6.1 6.1; 5.9 5.9],[3 1]);
assert(max(abs(series-[6.1;5.9])) < 1e-12, ...
    'Identical cell pH must be invariant to area weights.');

single_cell = pH_Area_Weighted([6.1;5.9],1);
assert(max(abs(single_cell-[6.1;5.9])) < 1e-12, ...
    'A one-cell 2-D domain must retain the 1-D pH series.');

assert_throws(@() pH_Area_Weighted([6 7],[-1 2]),'ERW:InvalidAreaWeights');
assert_throws(@() pH_Area_Weighted([6 NaN],[1 1]),'ERW:InvalidPH');
assert_throws(@() pH_Area_Weighted([6 7 8],[1 1]),'ERW:InvalidAreaWeights');

fprintf('PASS: 2-D domain pH uses area-weighted hydrogen activity.\n');
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
