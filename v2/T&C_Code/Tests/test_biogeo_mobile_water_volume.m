function test_biogeo_mobile_water_volume
% Mobile water storage must not depend on vertical mesh resolution.

test_dir = fileparts(mfilename('fullpath'));
addpath(fileparts(test_dir));

theta = 0.31;
theta_hyg = 0.12;

v_300 = Biogeo_Mobile_Water_Volume(300,theta,theta_hyg);
v_600 = Biogeo_Mobile_Water_Volume(600,theta,theta_hyg);

assert(abs(v_300-57) < 1e-12, ...
    'Mobile water must equal depth times mobile volumetric water content.');
assert(abs(v_600-2*v_300) < 1e-12, ...
    'Mobile water must scale with physical depth, not layer count.');
assert(Biogeo_Mobile_Water_Volume(300,0.10,0.12) == 0, ...
    'Water at or below the hygroscopic content must not be mobile.');

pool = zeros(1,76);
pool(32) = 10;
par = Biogeochemistry_Parameter(1);
v_profile = 230;
[~,leak_no3] = Biogeo_Leakage(pool,6,v_profile,par);
assert(abs(leak_no3-pool(32)*6/v_profile) < 1e-12, ...
    'Bottom NO3 leakage must use whole-profile mobile water.');

unit_source = fileread(fullfile(fileparts(test_dir),'BIOGEO_UNIT.m'));
profile_volume_call = ['Biogeo_Leakage\s*\(\s*Ptm1\s*,\s*Lk\s*,' ...
    '\s*V\s*,\s*BiogeoPar\s*\)'];
assert(~isempty(regexp(unit_source,profile_volume_call,'once')), ...
    'BIOGEO_UNIT must retain the legacy bottom-leakage water boundary.');

assert_throws(@() Biogeo_Mobile_Water_Volume(0,theta,theta_hyg), ...
    'Biogeo_Mobile_Water_Volume:InvalidDepth');
assert_throws(@() Biogeo_Mobile_Water_Volume(300,NaN,theta_hyg), ...
    'Biogeo_Mobile_Water_Volume:InvalidWaterContent');

fprintf('PASS: biogeochemical mobile water is physical and mesh independent.\n');
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
