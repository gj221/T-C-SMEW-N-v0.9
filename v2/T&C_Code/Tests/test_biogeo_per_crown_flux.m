function test_biogeo_per_crown_flux
% Ground-area uptake must not enter an inactive PFT through crown scaling.

test_dir = fileparts(mfilename('fullpath'));
addpath(fileparts(test_dir));

flux = [2 3 4];
suppression = [0.25 0.5 0];
crown = [0.5 0 0.25];
actual = Biogeo_Per_Crown_Flux(flux,suppression,crown);
expected = [3 0 16];
assert(max(abs(actual-expected)) < 1e-12, ...
    'Per-crown conversion is incorrect or supplies an inactive PFT.');

assert_throws(@() Biogeo_Per_Crown_Flux([1 2],0,[1 1]), ...
    'Biogeo_Per_Crown_Flux:SizeMismatch');
assert_throws(@() Biogeo_Per_Crown_Flux([1 -1],[0 0],[1 1]), ...
    'Biogeo_Per_Crown_Flux:InvalidFlux');
assert_throws(@() Biogeo_Per_Crown_Flux([1 1],[0 NaN],[1 1]), ...
    'Biogeo_Per_Crown_Flux:InvalidSuppression');
assert_throws(@() Biogeo_Per_Crown_Flux([1 1],[0 0],[1 -0.1]), ...
    'Biogeo_Per_Crown_Flux:InvalidCrown');

fprintf('PASS: per-crown uptake conversion preserves inactive-PFT boundaries.\n');
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
