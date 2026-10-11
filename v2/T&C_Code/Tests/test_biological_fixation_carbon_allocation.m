function test_biological_fixation_carbon_allocation
% Biological fixation must use the carbon paid to fixation on the same day.

test_dir = fileparts(mfilename('fullpath'));
code_dir = fileparts(test_dir);
addpath(code_dir);

carbon = [2 3 4];
temperature = 25;
cost = -6.25*(exp(-3.62+0.27*temperature*(1-0.5*temperature/25.15))-2);
fixed_n = Biogeo_Bio_fixation(0,0,zeros(1,32),carbon,temperature,1);

assert(abs(fixed_n-carbon(3)/cost) < 1e-12, ...
    'Only same-day carbon allocated to fixation may generate fixed N.');
assert(cost >= 7.5 && cost <= 12.5, ...
    'Fixation cost must remain within the Fisher et al. (2010) bounds.');
assert(Biogeo_Bio_fixation(0,0,zeros(1,32),carbon,temperature,0) == 0, ...
    'Fixation must be zero when the pathway is disabled.');

unit_source = fileread(fullfile(code_dir,'BIOGEO_UNIT.m'));
assert(contains(unit_source, ...
    'Biogeo_Bio_fixation(AAET,LAI,Ptm1,Rexmy,Ts,opt_BfixN)'), ...
    'BIOGEO_UNIT must not supply a moving-average carbon payment to fixation.');

fprintf('PASS: biological fixation uses its same-day carbon allocation.\n');
end
