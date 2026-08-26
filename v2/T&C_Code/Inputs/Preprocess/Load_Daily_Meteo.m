%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [met, info] = Load_Daily_Meteo(meteo_file, start_dn, end_dn, Lat, Lon, DeltaGMT)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isdatetime(start_dn); start_dn = datenum(start_dn); end %#ok<DATNM>
if isdatetime(end_dn);   end_dn   = datenum(end_dn);   end %#ok<DATNM>
if ~isfile(meteo_file)
    error('Load_Daily_Meteo:MissingFile', ...
        'Daily meteo file not found: %s', meteo_file);
end

C = readcell(meteo_file,'FileType','text','Delimiter','\t','NumHeaderLines',2);
yr = cell_numeric(C(:,1));
mo = cell_numeric(C(:,2));
da = cell_numeric(C(:,3));

% When start_dn or end_dn is NaN, auto-detect from the first/last complete
% calendar years present in the file (years having both Jan 1 and Dec 31).
if isnan(start_dn) || isnan(end_dn)
    valid_rows = isfinite(yr) & isfinite(mo) & isfinite(da);
    dn_all    = floor(datenum(yr(valid_rows), mo(valid_rows), da(valid_rows)));
    all_yr    = unique(yr(valid_rows));
    dn_jan1   = datenum(all_yr, ones(size(all_yr)),        ones(size(all_yr)));
    dn_dec31  = datenum(all_yr, 12*ones(size(all_yr)),     31*ones(size(all_yr)));
    complete  = ismember(dn_jan1, dn_all) & ismember(dn_dec31, dn_all);
    complete_yr = all_yr(complete);
    if isempty(complete_yr)
        error('Load_Daily_Meteo:NoCompleteYear', ...
            'No complete calendar years (Jan 1–Dec 31) found in %s', meteo_file);
    end
    if isnan(start_dn); start_dn = datenum(complete_yr(1),   1,  1); end
    if isnan(end_dn);   end_dn   = datenum(complete_yr(end), 12, 31); end
end

target_days = (floor(start_dn):floor(end_dn))';

ws_mph      = cell_numeric(C(:,6));
rs_MJ_m2_d  = cell_numeric(C(:,10));
ta_F        = cell_numeric(C(:,16));
rh_pct      = cell_numeric(C(:,22));
tdew_F      = cell_numeric(C(:,24));
pr_in       = cell_numeric(C(:,26));

date_dn = floor(datenum(yr,mo,da));
keep = ismember(date_dn,target_days);
date_dn = date_dn(keep);
ws_mph = ws_mph(keep);
rs_MJ_m2_d = rs_MJ_m2_d(keep);
ta_F = ta_F(keep);
rh_pct = rh_pct(keep);
tdew_F = tdew_F(keep);
pr_in = pr_in(keep);

[date_dn, ia] = unique(date_dn,'stable');
ws_mph = ws_mph(ia);
rs_MJ_m2_d = rs_MJ_m2_d(ia);
ta_F = ta_F(ia);
rh_pct = rh_pct(ia);
tdew_F = tdew_F(ia);
pr_in = pr_in(ia);

[has_day, loc] = ismember(target_days,date_dn);
if any(~has_day)
    error('Load_Daily_Meteo:MissingDays', ...
        'Meteo file is missing %d day(s): %s', ...
        sum(~has_day), format_date_list(target_days(~has_day)));
end

ws_mph = ws_mph(loc);
rs_MJ_m2_d = rs_MJ_m2_d(loc);
ta_F = ta_F(loc);
rh_pct = rh_pct(loc);
tdew_F = tdew_F(loc);
pr_in = pr_in(loc);

raw_missing = ~isfinite(ws_mph) | ~isfinite(rs_MJ_m2_d) | ~isfinite(ta_F) | ...
    ~isfinite(rh_pct) | ~isfinite(tdew_F) | ~isfinite(pr_in);
ws_mph = fill_meteo_series(target_days, ws_mph);
rs_MJ_m2_d = max(fill_meteo_series(target_days, rs_MJ_m2_d),0);
ta_F = fill_meteo_series(target_days, ta_F);
rh_pct = min(max(fill_meteo_series(target_days, rh_pct),0),100);
tdew_F = fill_meteo_series(target_days, tdew_F);
pr_in = max(fill_meteo_series(target_days, pr_in),0);

missing = ~isfinite(ws_mph) | ~isfinite(rs_MJ_m2_d) | ~isfinite(ta_F) | ...
    ~isfinite(rh_pct) | ~isfinite(tdew_F) | ~isfinite(pr_in);
if any(missing)
    error('Load_Daily_Meteo:MissingValues', ...
        'Meteo file has missing required values for %d day(s): %s', ...
        sum(missing), format_date_list(target_days(missing)));
end

n_days = numel(target_days);
n_hours = n_days * 24;
Date = zeros(n_hours,1);
Pr = zeros(n_hours,1);
Ta = zeros(n_hours,1);
Ws = zeros(n_hours,1);
ea = zeros(n_hours,1);
Pre = 1013.25 * ones(n_hours,1);
Tdew = zeros(n_hours,1);
esat = zeros(n_hours,1);
N = zeros(n_hours,1);
Rsw = zeros(n_hours,1);

for dd = 1:n_days
    ii = (dd-1)*24 + (1:24);
    Date(ii) = target_days(dd) + (0:23)'/24;
    Pr(ii) = max(pr_in(dd),0) * 25.4 / 24;
    Ta(ii) = (ta_F(dd) - 32) * 5/9;
    Ws(ii) = max(ws_mph(dd) * 0.44704, 0.01);
    Tdew(ii) = (tdew_F(dd) - 32) * 5/9;
    ea(ii) = 610.8 * exp(17.27*Tdew(ii)./(Tdew(ii)+237.3));
    esat(ii) = 610.8 * exp(17.27*Ta(ii)./(Ta(ii)+237.3));
    N(ii) = min(max(1 - rh_pct(dd)/100,0),1);
    Rsw(ii) = distribute_daily_radiation(target_days(dd), rs_MJ_m2_d(dd), Lat, Lon, DeltaGMT);
end

direct_frac = min(max(1 - N,0.15),0.85);
diffuse_frac = 1 - direct_frac;
SAB1 = 0.5 * Rsw .* direct_frac;
SAB2 = 0.5 * Rsw .* direct_frac;
SAD1 = 0.5 * Rsw .* diffuse_frac;
SAD2 = 0.5 * Rsw .* diffuse_frac;
PARB = 0.45 * (SAB1 + SAB2);
PARD = 0.45 * (SAD1 + SAD2);

met = struct('Date',Date,'Pr',Pr,'Ta',Ta,'Ws',Ws,'ea',ea,'Pre',Pre, ...
    'SAD1',SAD1,'SAD2',SAD2,'SAB1',SAB1,'SAB2',SAB2, ...
    'N',N,'Tdew',Tdew,'esat',esat,'PARB',PARB,'PARD',PARD);

info = struct();
info.source = meteo_file;
info.n_days = n_days;
info.n_hours = n_hours;
info.start_date = target_days(1);
info.end_date = target_days(end);
info.n_gapfilled_days = sum(raw_missing);
end


function vals = distribute_daily_radiation(day_dn, rs_MJ_m2_d, Lat, Lon, DeltaGMT)
vals = zeros(24,1);
if rs_MJ_m2_d <= 0
    return
end
shape = zeros(24,1);
[yy,mo,da] = datevec(day_dn);
for hh = 0:23
    [h_S,~,~,~,~,~] = SetSunVariables([yy,mo,da,hh],DeltaGMT,Lon,Lat,-0.5,0.5);
    shape(hh+1) = max(sin(h_S),0);
end
if sum(shape) <= 0
    return
end
daily_J_m2 = rs_MJ_m2_d * 1e6;
vals = daily_J_m2 * shape / (sum(shape) * 3600);
end


function v = cell_numeric(x)
v = NaN(size(x));
for ii = 1:numel(x)
    xi = x{ii};
    if isnumeric(xi)
        if ~isempty(xi); v(ii) = xi; end
    elseif isstring(xi) || ischar(xi)
        txt = strtrim(char(xi));
        if isempty(txt) || strcmp(txt,'----') || strcmpi(txt,'M')
            v(ii) = NaN;
        else
            v(ii) = str2double(txt);
        end
    end
end
v = v(:);
end


function v = fill_meteo_series(days, v)
missing = ~isfinite(v);
if ~any(missing)
    return
end
good = isfinite(v);
if sum(good) < 2
    error('Load_Daily_Meteo:MissingValues', ...
        'Meteo series cannot be gap-filled because fewer than two values are valid.');
end
v(missing) = interp1(days(good), v(good), days(missing), 'linear', 'extrap');
end


function txt = format_date_list(dn)
n = min(numel(dn), 8);
txt = strjoin(cellstr(datestr(dn(1:n),'yyyy-mm-dd')), ', ');
if numel(dn) > n
    txt = [txt, ', ...'];
end
end
