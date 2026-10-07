function report = Compare_Forcing_Sources(legacy_file, era5_file, DeltaGMT)
% Is the legacy forcing the same ERA5-Land series as the new build, and if
% not, where do they differ?  Data_Kantola.mat matches the new
% energyfarm_1995_2024.nc build on temperature and wind to two decimals
% but carries 3-4 % less precipitation.  A de-accumulation that drops one
% hour of each UTC day loses 1/24 = 4.2 % and would leave that hour of day
% empty, so precipitation is also compared by local hour of day.
%
%   legacy_file  T&C forcing .mat with Date (local datenum), Pr [mm/h], Ta,
%                Ws, SAD1+SAD2+SAB1+SAB2 [W/m2]
%   era5_file    ERA5-Land point NetCDF: valid_time (seconds since 1970,
%                UTC), tp [m/h], t2m [K], u10, v10, ssrd [J/m2 per hour]
%   DeltaGMT     local standard time offset, hours (-6 for the Energy Farm)
%
% Returns the per-year table and the hour-of-day ratios.
if nargin < 1 || isempty(legacy_file)
    legacy_file = fullfile(fileparts(fileparts(mfilename('fullpath'))),'Data_Kantola.mat');
end
if nargin < 2 || isempty(era5_file)
    era5_file = fullfile(fileparts(fileparts(mfilename('fullpath'))),'energyfarm_1995_2024.nc');
end
if nargin < 3 || isempty(DeltaGMT), DeltaGMT = -6; end

L = load(legacy_file, 'Date', 'Pr', 'Ta', 'Ws', 'SAD1', 'SAD2', 'SAB1', 'SAB2');
tL = round(L.Date(:)*24);                       % integer local hours
PL = L.Pr(:); TL = L.Ta(:); WL = L.Ws(:);
RL = L.SAD1(:) + L.SAD2(:) + L.SAB1(:) + L.SAB2(:);

tv = double(ncread(era5_file, 'valid_time'));
utc = datetime(1970,1,1) + seconds(tv(:));
tE = round(datenum(utc + hours(DeltaGMT))*24);
PE = 1000*max(double(ncread(era5_file,'tp')),0); PE = PE(:);
TE = double(ncread(era5_file,'t2m')) - 273.15; TE = TE(:);
WE = hypot(double(ncread(era5_file,'u10')), double(ncread(era5_file,'v10'))); WE = WE(:);
RE = max(double(ncread(era5_file,'ssrd')),0)/3600; RE = RE(:);

[~, iL, iE] = intersect(tL, tE);
if isempty(iL)
    error('Compare_Forcing_Sources:NoOverlap','The two series share no hour.');
end
yr = year(datetime(tL(iL)/24, 'ConvertFrom', 'datenum'));
hod = mod(tL(iL), 24);                              % local hour of day

fprintf('\nlegacy %s vs ERA5-Land %s, %d shared hours\n', legacy_file, era5_file, numel(iL));
fprintf('  year  P_legacy  P_era5  ratio |  Ta_leg  Ta_era5 | Rsw_leg Rsw_era5 ratio | Ws_leg Ws_era5\n');
T = [];
for y = unique(yr).'
    m = yr == y;
    if sum(m) < 24*360, continue, end               % whole years only
    row = [y sum(PL(iL(m))) sum(PE(iE(m))) mean(TL(iL(m))) mean(TE(iE(m))) ...
        mean(RL(iL(m))) mean(RE(iE(m))) mean(WL(iL(m))) mean(WE(iE(m)))];
    T = [T; row]; %#ok<AGROW>
    fprintf('  %4d %9.0f %7.0f %6.3f | %7.2f %8.2f | %7.1f %8.1f %5.3f | %6.2f %7.2f\n', ...
        row(1), row(2), row(3), row(2)/row(3), row(4), row(5), row(6), row(7), row(6)/row(7), row(8), row(9));
end

fprintf('\n  precipitation by local hour of day, whole shared record\n');
fprintf('  hour  legacy_mm  era5_mm  ratio  era5_share\n');
ratio_h = nan(1,24);
for h = 0:23
    m = hod == h;
    a = sum(PL(iL(m))); b = sum(PE(iE(m)));
    ratio_h(h+1) = a/max(b, eps);
    fprintf('  %4d %10.0f %8.0f %6.3f %10.3f\n', h, a, b, ratio_h(h+1), b/sum(PE(iE)));
end
[~, worst] = min(ratio_h);
missing = sum(PE(iE(hod == worst-1))) - sum(PL(iL(hod == worst-1)));
total_gap = sum(PE(iE)) - sum(PL(iL));
fprintf(['  --> hour %d (local, = %02d UTC) carries the largest deficit: %.0f mm of a ' ...
    '%.0f mm total gap (%.0f %%)\n'], worst-1, mod(worst-1-DeltaGMT,24), missing, total_gap, ...
    100*missing/max(total_gap, eps));

report = struct('annual', T, 'annual_columns', ...
    {{'year','P_legacy','P_era5','Ta_legacy','Ta_era5','Rsw_legacy','Rsw_era5','Ws_legacy','Ws_era5'}}, ...
    'hour_ratio', ratio_h, 'worst_local_hour', worst-1, 'worst_hour_gap_mm', missing, ...
    'total_gap_mm', total_gap);
end
