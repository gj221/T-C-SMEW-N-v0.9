function test_denitrification_nitrate_limitation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The TC-SMEW denitrification sink must retain its legacy first-order
%% nitrate dependence. A Monod fraction cannot be added without the
%% denitrifier biomass and maximum-capacity terms used by DNDC.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

here = fileparts(mfilename('fullpath'));
code = fileparts(here);
addpath(code,fullfile(code,'Inputs'),fullfile(code,'Parameters'));
fprintf('\n=== test_denitrification_nitrate_limitation ===\n');

F = load(fullfile(here,'Fixtures','no3_immobilization_20190906.mat'),'fixture');
args = F.fixture.args;
args{end+1} = 0.8;
args{50} = Biogeochemistry_Parameter(1);
for k = 12:47
    args{k} = 0;
end
args{5}(:) = 0;
args{52}(:) = 0;
args{6} = 20;
args{7} = 20;
args{9} = 5;
args{10} = 1;
args{11} = 0.5;
args{55} = 1;

NO3_lo = 1;
NO3_hi = 10;
DOC_concentration = 10*args{50}.Kc;
DOC = DOC_concentration*args{3};

rate_lo = denitrification_rate(args,DOC,NO3_lo);
rate_hi = denitrification_rate(args,DOC,NO3_hi);
observed_ratio = rate_hi/rate_lo;
expected_ratio = NO3_hi/NO3_lo;
fprintf('rates %.6g %.6g, ratio %.6g\n', ...
    rate_lo,rate_hi,observed_ratio);

assert(rate_lo > 0 && rate_hi > 0, ...
    'The nitrate-response fixture must produce finite positive denitrification.');
assert(abs(observed_ratio-expected_ratio) < 1e-8, ...
    ['Denitrification does not retain the legacy first-order nitrate ' ...
    'response: observed ratio %.6g, expected %.6g.'], ...
    observed_ratio,expected_ratio);
fprintf('PASS: first-order nitrate rate ratio %.6g matches %.6g\n', ...
    observed_ratio,expected_ratio);

depths = [0.15 0.30];
fractional_rate = zeros(size(depths));
NO3_concentration = 30;
for k = 1:numel(depths)
    depth_args = args;
    depth_args{3} = depths(k);
    DOC_areal = DOC_concentration*depths(k);
    NO3_areal = NO3_concentration*depths(k);
    fractional_rate(k) = ...
        denitrification_rate(depth_args,DOC_areal,NO3_areal)/NO3_areal;
end
assert(abs(diff(fractional_rate)) < 1e-12, ...
    ['The DOC half-saturation changes with layer depth. Kc must retain ' ...
    'its source volumetric basis.']);
fprintf('PASS: volumetric DOC limitation is independent of layer depth\n');

cap_args = args;
cap_args{50}.KNmaxd = 100;
NO3_available = 1;
cap_rate = denitrification_rate(cap_args,DOC,NO3_available);
assert(abs(cap_rate-NO3_available) < 1e-10, ...
    'The availability cap must prevent denitrification from overdrawing NO3.');
fprintf('PASS: the denitrification availability cap preserves NO3 mass\n');

bio_source = fileread(fullfile(code,'BIOGEOCHEMISTRY_DYNAMIC3.m'));
assert(contains(bio_source,'B(31)= B(31)- NH4_Uptake') && ...
       contains(bio_source,'B(32)=B(32)-NO3_Uptake'), ...
    'Working mineral-N pools must reserve plant uptake once.');
assert(contains(bio_source,'(1-f_org_lea)*Ip_lea_nit - LEAK_NH4'), ...
    'NH4 reaction availability must use the post-uptake working pool.');
assert(contains(bio_source,'N2Oflx_nit - LEAK_NO3'), ...
    'NO3 reaction availability must use the post-uptake working pool.');
fix_source = fileread(fullfile(code,'Biogeo_Bio_fixation.m'));
assert(contains(fix_source,'REXMY(3)/Costfix'), ...
    'Biological N fixation must use only carbon allocated to fixation.');
fprintf('PASS: mineral-N reactions reserve uptake and fixation uses its C allocation\n');
end

function rate = denitrification_rate(args,DOC,NO3)
args{2}(1:34) = max(args{2}(1:34),1e-3);
args{2}(31:32) = 0;
args{2}(12) = DOC;
args{2}(13) = 0;
args{2}(32) = NO3;
out = cell(1,18);
[out{:}] = BIOGEOCHEMISTRY_DYNAMIC3(args{:});
rate = out{17};
end
