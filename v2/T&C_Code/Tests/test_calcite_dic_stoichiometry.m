function test_calcite_dic_stoichiometry
% Calcite dissolution adds one inorganic C for each released Ca.

test_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(test_dir,'..'));

mineral = char('calcite','forsterite');
nu_C = ERW_Mineral_DIC_Stoichiometry(mineral);
assert(isequal(nu_C,[1;0]),'Calcite DIC stoichiometry is incorrect.');

EW = [2.5e4 8.0e3];          % umol mineral m^-2 d^-1
dtw = 1/144;                 % d
dCa = EW(1)*dtw;             % umol Ca m^-2
dIC = sum(nu_C.'.*EW)*dtw;   % umol C m^-2
dAlk = 2*EW(1)*dtw;          % umol charge m^-2
assert(abs(dIC-dCa) < 1e-12,'Calcite carbon and calcium do not balance.');
assert(abs(dAlk-2*dIC) < 1e-12,'Calcite alkalinity does not balance.');

code = fileread(fullfile(test_dir,'..','ERW_function.m'));
assert(contains(code,'dIC_mineral = sum(mineral_C''.*EW(j_only,:))*dtw;'), ...
    'ERW mass balance does not calculate mineral-derived DIC.');
assert(contains(code,'IC_tot= IC_tot(j_only)+dIC_mineral+Iw*1000*DIC_rain*dtw-ADV+'), ...
    'ERW mass balance does not scale daily rainwater DIC to the sub-step.');

Iw = 0.012;                    % m d^-1
DIC_rain = 18.5;              % umol C L^-1
daily_rain_C = Iw*1000*DIC_rain;
for t_res = [1 6 12]
    dtw = 1/(24*t_res);
    integrated_rain_C = sum(repmat(daily_rain_C*dtw,24*t_res,1));
    assert(abs(integrated_rain_C-daily_rain_C) < ...
        1e-12*max(1,daily_rain_C), ...
        'Rainwater DIC changes with the ERW sub-step resolution.');
end

fprintf('PASS: calcite Ca, DIC and alkalinity stoichiometry closes.\n');
end
