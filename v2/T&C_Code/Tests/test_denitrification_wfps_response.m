function test_denitrification_wfps_response
% Denitrification partition uses the DayCent WFPS response once.

here = fileparts(mfilename('fullpath'));
code = fileparts(here);
addpath(code,fullfile(code,'Inputs'),fullfile(code,'Parameters'));

theta = [0 0.18 0.45];
theta_sat = 0.45;
WFPS = Biogeo_WFPS(theta,theta_sat);
assert(max(abs(WFPS-[0 0.4 1])) < 1e-12, ...
    'WFPS must be volumetric water content divided by saturated content.');

expected = 1.4./13.^(17./13.^(2.2*WFPS));
actual = Biogeo_Denitrification_WFPS_Response(WFPS);
assert(max(abs(actual-expected)) < 1e-12, ...
    'Denitrification WFPS response does not match the DayCent equation.');
assert(all(diff(actual) > 0), ...
    'The DayCent N2:N2O moisture factor must increase with WFPS.');

assert_throws(@() Biogeo_WFPS(-eps,theta_sat), ...
    'Biogeo_WFPS:InvalidWaterContent');
assert_throws(@() Biogeo_WFPS(theta_sat+eps,theta_sat), ...
    'Biogeo_WFPS:InvalidWaterContent');
assert_throws(@() Biogeo_Denitrification_WFPS_Response(NaN), ...
    'Biogeo_Denitrification_WFPS_Response:InvalidWFPS');

F = load(fullfile(here,'Fixtures','no3_immobilization_20190906.mat'),'fixture');
args = F.fixture.args;
args{50} = Biogeochemistry_Parameter(1);
for k = 12:47
    args{k} = 0;
end
args{5}(:) = 0;
args{52}(:) = 0;
args{6} = 20;
args{7} = 20;
args{9} = 5.5;
args{10} = 0.8;
args{11} = 0.5;
args{55} = 1;
args{2}(1:34) = max(args{2}(1:34),1e-3);
args{2}(12) = 10*args{50}.Kc*args{3};
args{2}(13) = 0;
args{2}(31) = 0;
args{2}(32) = args{50}.Kn*args{3};

[n2_dry,n2o_dry,denit_dry,nit_n2o_dry,resp_dry] = partition_flux(args,0.4);
[n2_wet,n2o_wet,denit_wet,nit_n2o_wet,resp_wet] = partition_flux(args,0.9);
soil_mass = args{3}*args{4}*1000; %% [g soil/m2]
NO3_ppm = args{2}(32)/soil_mass*1e6; %% [microgram N/g soil]
[~,pH_ratio] = Biogeo_Denitrification_PH_Response(args{9});
expected_dry_ratio = Biogeo_Denitrification_Product_Ratio( ...
    pH_ratio,0.4,NO3_ppm,resp_dry);
expected_wet_ratio = Biogeo_Denitrification_Product_Ratio( ...
    pH_ratio,0.9,NO3_ppm,resp_wet);

assert(denit_dry > 0 && abs(denit_dry-denit_wet) < 1e-12, ...
    'WFPS partition factor must not change gross denitrification.');
assert(abs(n2_dry+n2o_dry-denit_dry) < 1e-12 && ...
       abs(n2_wet+n2o_wet-denit_wet) < 1e-12, ...
    'N2 and N2O partitioning must conserve total denitrified N.');
assert(abs(n2_dry/n2o_dry-expected_dry_ratio) < 1e-12 && ...
       abs(n2_wet/n2o_wet-expected_wet_ratio) < 1e-12, ...
    'Production N2:N2O does not equal the DayCent product ratio.');
assert(n2_wet > n2_dry && n2o_wet < n2o_dry, ...
    'Higher WFPS must shift denitrification from N2O toward N2.');
assert(nit_n2o_dry == 0 && nit_n2o_wet == 0, ...
    'The partition fixture must isolate denitrification N2O.');

fprintf('PASS: DayCent WFPS modifies denitrification partition once and conserves N.\n');
end

function [n2,n2o,denit,nit_n2o,respiration] = partition_flux(args,WFPS)
call_args = args;
call_args{end+1} = WFPS;
[~,R_litter,R_microbe,~,R_ew,~,n2,~,~,~,RmycAM,RmycEM,~,~,n2o,~,denit,nit_n2o] = ...
    BIOGEOCHEMISTRY_DYNAMIC3(call_args{:});
respiration = (R_litter + R_microbe + R_ew - RmycAM - RmycEM)*10; %% [kg C/ha d]
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
