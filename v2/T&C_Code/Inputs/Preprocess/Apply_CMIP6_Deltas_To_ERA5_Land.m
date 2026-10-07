function report = Apply_CMIP6_Deltas_To_ERA5_Land( ...
    era5_source,output_file,deltas,Lat,Lon,Zbas,DeltaGMT,varargin)
% Apply CMIP6 monthly change factors to an ERA5-Land hourly point series and
% write a perturbed NetCDF in the same format, so that
% Build_TC_Forcing_From_ERA5_Land_Timeseries can consume it unmodified.
%
% deltas is the struct (or .mat path) produced by Extract_CMIP6_Monthly_Deltas.
% Output years are assembled from source years through YearMap, so the run can
% extend past the observed record by resampling baseline years.
%
% Perturbations, all scaled by a linear ramp keyed to the target year:
%   t2m   additive dT
%   d2m   recomputed at relative humidity perturbed by dRH, against warmed t2m
%   tp    multiplicative fP, or quantile delta mapping when PrecipMethod='qdm'
%   ssrd  multiplicative fSW
%   u10,v10  multiplicative fWs, preserving direction
%   sp    unchanged
%
% Name-value options:
%   YearMap          [target source; ...]  required unless the record is used as-is
%   PassthroughYears []      target years written from observations with ramp 0
%   RampReferenceYear  mean(deltas.BaselineYears)   year at which ramp = 0
%   RampTargetYear     mean(deltas.FutureYears)     year at which ramp = 1
%   PrecipMethod     'delta' | 'qdm'
%   Deaccumulation   (required) 'legacy' | 'timeseries'. 'legacy' first
%                    re-applies the de-accumulation Data_Kantola.mat was
%                    built with to tp and ssrd: each hour ending 01 UTC
%                    becomes that hour less the 24 before it, clipped at
%                    zero. 'timeseries' keeps the increments as downloaded.
%                    See Delta_Forcing_Config.
%   SiteLat          (required) latitude written to the perturbed file
%   SiteLon          (required) longitude written to the perturbed file
%                    Lat and Lon above identify the ERA5-Land cell and are
%                    validated against the source. SiteLat and SiteLon are the
%                    location the model runs at, which downstream feeds the
%                    solar geometry. They differ: the cell centre is not the
%                    field. Not defaulted, because defaulting to the cell
%                    would silently reintroduce the mismatch.
%   Verbose          true

narginchk(7,inf);
opt = parse_options(varargin);

era5_source = char(era5_source);
output_file = char(output_file);
if ischar(deltas)
    deltas = load(deltas);
end
validate_deltas(deltas);
validateattributes(Lat,{'numeric'},{'scalar','real','finite','>=',-90,'<=',90});
validateattributes(Lon,{'numeric'},{'scalar','real','finite','>=',-180,'<=',180});
validateattributes(Zbas,{'numeric'},{'scalar','real','finite'});
validateattributes(DeltaGMT,{'numeric'},{'scalar','real','finite','>=',-12,'<=',14});
[output_dir,~,output_ext] = fileparts(output_file);
if ~isempty(output_dir) && ~isfolder(output_dir)
    error('Apply_CMIP6_Deltas_To_ERA5_Land:MissingOutputDirectory', ...
        'Output directory not found: %s',output_dir);
end
if ~strcmpi(output_ext,'.nc')
    error('Apply_CMIP6_Deltas_To_ERA5_Land:BadOutputExtension', ...
        'Output file must have a .nc extension.');
end
if isempty(opt.RampReferenceYear)
    opt.RampReferenceYear = mean(deltas.BaselineYears);
end
if isempty(opt.RampTargetYear)
    opt.RampTargetYear = mean(deltas.FutureYears);
end
if opt.RampTargetYear<=opt.RampReferenceYear
    error('Apply_CMIP6_Deltas_To_ERA5_Land:BadRamp', ...
        'RampTargetYear must be later than RampReferenceYear.');
end

src = read_source(era5_source,DeltaGMT);
check_source_coordinates(src,Lat,Lon);
if strcmp(opt.Deaccumulation,'legacy')
    src.tp = legacy_deaccumulation(src.tp,src.utc);
    src.ssrd = legacy_deaccumulation(src.ssrd,src.utc);
end

if isempty(opt.YearMap)
    years = unique(src.year(:))';
    opt.YearMap = [years(:) years(:)];
end
validate_year_map(opt.YearMap,src,opt.PassthroughYears);

% ---- assemble the output series from source years -----------------------
nmap = size(opt.YearMap,1);
blocks = cell(nmap,1);
for k = 1:nmap
    target = opt.YearMap(k,1);
    source = opt.YearMap(k,2);
    sel = find(src.year==source);
    if isempty(sel)
        error('Apply_CMIP6_Deltas_To_ERA5_Land:MissingSourceYear', ...
            'Source year %d is not present in the ERA5-Land record.',source);
    end
    b = struct();
    b.target = target;
    b.idx = sel;
    b.utc = src.utc(sel) + calyears(target-source);
    b.ramp = ramp_factor(target,opt);
    if ismember(target,opt.PassthroughYears)
        b.ramp = 0;
    end
    blocks{k} = b;
end

idx = vertcat(blocks_field(blocks,'idx'));
utc = vertcat(blocks_field(blocks,'utc'));
ramp = cell2mat(arrayfun(@(k) repmat(blocks{k}.ramp,numel(blocks{k}.idx),1), ...
    (1:nmap)','UniformOutput',false));

[utc,order] = sort(utc);
idx = idx(order);
ramp = ramp(order);

local = utc + hours(DeltaGMT);
local.TimeZone = '';
Date = datenum(local);
if any(abs(diff(Date)-1/24)>1e-8)
    error('Apply_CMIP6_Deltas_To_ERA5_Land:NonHourlyOutput', ...
        ['Assembled output is not a strictly increasing hourly series. ' ...
        'Check YearMap for duplicate or overlapping target years.']);
end
mo = month(local);

t2m = src.t2m(idx);
d2m = src.d2m(idx);
sp = src.sp(idx);
u10 = src.u10(idx);
v10 = src.v10(idx);
tp = src.tp(idx);
ssrd = src.ssrd(idx);

Rsw = max(ssrd,0)/3600;

% ---- temperature and humidity ------------------------------------------
Ta = t2m - 273.15;
Td = d2m - 273.15;
RH = 100*magnus(Td)./magnus(Ta);
RH = min(max(RH,1),100);

Ta_new = Ta + ramp.*column(deltas.dT(mo));
RH_new = RH + ramp.*column(deltas.dRH(mo));
RH_new = min(max(RH_new,1),100);
ea_new = (RH_new/100).*magnus(Ta_new);
Td_new = inverse_magnus(ea_new);
Td_new = min(Td_new,Ta_new);

t2m_new = Ta_new + 273.15;
d2m_new = Td_new + 273.15;

% ---- precipitation ------------------------------------------------------
switch opt.PrecipMethod
    case 'delta'
        tp_new = tp .* (1 + ramp.*(column(deltas.fP(mo))-1));
    case 'qdm'
        if ~isfield(deltas,'fP_quantile')
            error('Apply_CMIP6_Deltas_To_ERA5_Land:MissingQuantiles', ...
                ['PrecipMethod ''qdm'' needs fP_quantile. Re-run ' ...
                'Extract_CMIP6_Monthly_Deltas with PrecipQuantiles set.']);
        end
        tp_new = apply_qdm(tp,mo,ramp,src.tp,month(src.local), ...
            deltas.PrecipQuantiles,deltas.fP_quantile);
end
tp_new = max(tp_new,0);

% ---- shortwave ----------------------------------------------------------
% No clear-sky ceiling. Automatic_Radiation_Partition_I already reconciles
% observed shortwave against clear sky in Estimate_CloudCover, which caps the
% clearness index at 1 and treats such hours as cloudless rather than reducing
% the radiation. Its tenth output is the cloud-adjusted all-sky reconstruction,
% not a clear-sky bound, so a ceiling built from it clipped most daylight hours
% and rewrote the passthrough years, where ramp is 0 and the observation must
% survive untouched. Neither legacy T&C-SMEW nor TeC_BG_2D clips.
Rsw_new = Rsw .* (1 + ramp.*(column(deltas.fSW(mo))-1));
Rsw_new = max(Rsw_new,0);
ssrd_new = Rsw_new*3600;

% ---- wind ---------------------------------------------------------------
wind_scale = 1 + ramp.*(column(deltas.fWs(mo))-1);
u10_new = u10 .* wind_scale;
v10_new = v10 .* wind_scale;

% ---- guards mirroring the downstream builder ---------------------------
if any(~isfinite([t2m_new;d2m_new;sp;u10_new;v10_new;tp_new;ssrd_new]))
    error('Apply_CMIP6_Deltas_To_ERA5_Land:NonFiniteOutput', ...
        'Perturbation produced NaN or infinite values.');
end
if any(t2m_new<=0) || any(d2m_new<=0)
    error('Apply_CMIP6_Deltas_To_ERA5_Land:BadTemperature', ...
        'Perturbed temperatures are not positive in kelvin.');
end
if any(d2m_new>t2m_new+1e-9)
    error('Apply_CMIP6_Deltas_To_ERA5_Land:SupersaturatedOutput', ...
        'Perturbed dew point exceeds air temperature.');
end

% The perturbed file carries the site coordinates, so the downstream
% builder validates against them and writes them into the T&C forcing.
write_era5_like(output_file,utc,opt.SiteLat,opt.SiteLon, ...
    t2m_new,d2m_new,sp,u10_new,v10_new,tp_new,ssrd_new,deltas,opt);

report = struct();
report.n = numel(Date);
report.years = [min(year(local)) max(year(local))];
report.mean_ramp = mean(ramp);
report.scenario = deltas.scenario;

if opt.Verbose
    fprintf('Wrote %s: %d hourly records, %d-%d, scenario %s.\n', ...
        output_file,report.n,report.years(1),report.years(2),deltas.scenario);
end
end

% ---------------------------------------------------------------------

function e = magnus(T_C)
e = 610.8*exp(17.27*T_C./(237.3+T_C));
end

function T_C = inverse_magnus(e)
L = log(max(e,1e-6)/610.8);
T_C = 237.3*L./(17.27-L);
end

function v = column(x)
v = reshape(x,[],1);
end

function out = blocks_field(blocks,name)
out = cellfun(@(b) b.(name)(:),blocks,'UniformOutput',false);
out = vertcat(out{:});
end

function r = ramp_factor(target,opt)
r = (target-opt.RampReferenceYear)/(opt.RampTargetYear-opt.RampReferenceYear);
r = min(max(r,0),1);
end

function x = legacy_deaccumulation(x,utc)
% ERA5-Land accumulates from 00 UTC and its 00 UTC value holds the whole
% previous day, so differencing the accumulations across the reset, as
% Data_Kantola.mat was built, gives the hour ending 01 UTC its own increment
% less the previous day's total, clipped at zero. The previous day is the
% 24 increments before the hour, fewer at the start of the record.
if any(abs(diff(utc)-hours(1))>seconds(1))
    error('Apply_CMIP6_Deltas_To_ERA5_Land:NonHourlySource', ...
        'The legacy de-accumulation needs a gap-free hourly record.');
end
increment = x;
for k = find(hour(utc)==1)'
    x(k) = max(increment(k)-sum(increment(max(k-24,1):k-1)),0);
end
end

function tp_new = apply_qdm(tp,mo,ramp,tp_all,mo_all,probs,ratios)
% Quantile delta mapping on wet hours, per calendar month. The ratios come
% from daily model quantiles and are applied at the hourly quantile position,
% which assumes the daily intensity change carries to sub-daily intensity.
wet = 1e-5;   % 0.01 mm expressed in metres
tp_new = tp;
probs = probs(:)';
for m = 1:12
    ref = tp_all(mo_all==m);
    ref = sort(ref(ref>wet));
    if numel(ref)<50
        continue
    end
    sel = find(mo==m & tp>wet);
    if isempty(sel)
        continue
    end
    p = position_in(ref,tp(sel));
    r = interp1(probs,ratios(m,:),p,'linear','extrap');
    r = min(max(r,0.1),10);
    tp_new(sel) = tp(sel) .* (1 + ramp(sel).*(r(:)-1));
end
end

function p = position_in(sorted_ref,values)
% Empirical CDF position of each value within a sorted reference sample.
% Counted by interpolation on the unique support so a 35-year hourly record
% costs O((n+m)log n) rather than a pairwise comparison.
n = numel(sorted_ref);
u = unique(sorted_ref(:));
cnt = cumsum(histcounts(sorted_ref(:),[u; inf])');
rank = interp1(u,cnt,values(:),'previous','extrap');
rank(values(:)<u(1)) = 0;
rank(values(:)>=u(end)) = n;
p = (rank-0.5)/n;
p = min(max(p,1/(2*n)),1-1/(2*n));
end

function src = read_source(era5_source,DeltaGMT)
if isfile(era5_source)
    files = {era5_source};
elseif isfolder(era5_source)
    listing = dir(fullfile(era5_source,'*.nc'));
    if isempty(listing)
        error('Apply_CMIP6_Deltas_To_ERA5_Land:NoNetCDF', ...
            'No NetCDF files were found in: %s',era5_source);
    end
    [~,ord] = sort({listing.name});
    listing = listing(ord);
    files = arrayfun(@(x) fullfile(x.folder,x.name),listing,'UniformOutput',false);
else
    error('Apply_CMIP6_Deltas_To_ERA5_Land:MissingInput', ...
        'ERA5-Land NetCDF source not found: %s',era5_source);
end

wanted = { ...
    't2m',{'t2m','2m_temperature'}; ...
    'd2m',{'d2m','2m_dewpoint_temperature'}; ...
    'sp',{'sp','surface_pressure'}; ...
    'u10',{'u10','10m_u_component_of_wind'}; ...
    'v10',{'v10','10m_v_component_of_wind'}; ...
    'tp',{'tp','total_precipitation'}; ...
    'ssrd',{'ssrd','surface_solar_radiation_downwards'}};

src = struct();
src.time_units = '';
for k = 1:size(wanted,1)
    [f,name] = find_variable(files,wanted{k,2});
    src.(wanted{k,1}) = double(squeeze(ncread(f,name)));
    src.(wanted{k,1}) = src.(wanted{k,1})(:);
    if isempty(src.time_units)
        tname = find_variable_name(f,{'valid_time','time'});
        src.utc = read_utc(f,tname);
        src.time_units = ncreadatt(f,tname,'units');
        src.lat = double(squeeze(ncread(f,find_variable_name(f,{'latitude','lat'}))));
        src.lon = double(squeeze(ncread(f,find_variable_name(f,{'longitude','lon'}))));
    end
end

n = numel(src.utc);
names = wanted(:,1)';
for k = 1:numel(names)
    if numel(src.(names{k}))~=n
        error('Apply_CMIP6_Deltas_To_ERA5_Land:LengthMismatch', ...
            'ERA5 variables do not share one hourly time axis.');
    end
end
src.local = src.utc + hours(DeltaGMT);
src.local.TimeZone = '';
src.year = year(src.local);
end

function [f,name] = find_variable(files,candidates)
for k = 1:numel(candidates)
    for j = 1:numel(files)
        available = {ncinfo(files{j}).Variables.Name};
        idx = find(strcmpi(available,candidates{k}),1);
        if ~isempty(idx)
            f = files{j};
            name = available{idx};
            return
        end
    end
end
error('Apply_CMIP6_Deltas_To_ERA5_Land:MissingVariable', ...
    'Required ERA5 variable not found. Expected one of: %s', ...
    strjoin(candidates,', '));
end

function name = find_variable_name(f,candidates)
available = {ncinfo(f).Variables.Name};
for k = 1:numel(candidates)
    idx = find(strcmpi(available,candidates{k}),1);
    if ~isempty(idx)
        name = available{idx};
        return
    end
end
error('Apply_CMIP6_Deltas_To_ERA5_Land:MissingVariable', ...
    'Required ERA5 variable not found. Expected one of: %s', ...
    strjoin(candidates,', '));
end

function utc = read_utc(f,tname)
values = double(squeeze(ncread(f,tname)));
units = char(ncreadatt(f,tname,'units'));
parts = regexp(strtrim(units),'^(seconds|hours|days)\s+since\s+(.+)$', ...
    'tokens','once','ignorecase');
if isempty(parts)
    error('Apply_CMIP6_Deltas_To_ERA5_Land:BadTimeUnits', ...
        'Unsupported NetCDF time units: %s',units);
end
origin_text = strrep(strtrim(parts{2}),'T',' ');
origin_text = regexprep(origin_text,'(\.\d+)?\s*(UTC|Z|[+]00:?00)?$','');
formats = {'yyyy-MM-dd HH:mm:ss','yyyy-MM-dd HH:mm','yyyy-MM-dd'};
origin = NaT;
for k = 1:numel(formats)
    try
        origin = datetime(origin_text,'InputFormat',formats{k},'TimeZone','UTC');
        break
    catch
    end
end
if isnat(origin)
    error('Apply_CMIP6_Deltas_To_ERA5_Land:BadTimeOrigin', ...
        'Could not parse NetCDF time origin: %s',origin_text);
end
switch lower(parts{1})
    case 'seconds', utc = origin + seconds(values(:));
    case 'hours',   utc = origin + hours(values(:));
    case 'days',    utc = origin + days(values(:));
end
end

function check_source_coordinates(src,Lat,Lon)
if numel(src.lat)~=1 || numel(src.lon)~=1
    error('Apply_CMIP6_Deltas_To_ERA5_Land:NotPointData', ...
        'ERA5 source must be a one-point time series.');
end
if abs(src.lat-Lat)>1e-8 || ...
        abs(mod(src.lon+180,360)-180 - (mod(Lon+180,360)-180))>1e-8
    error('Apply_CMIP6_Deltas_To_ERA5_Land:CoordinateMismatch', ...
        'ERA5 source coordinates do not match Lat and Lon.');
end
end

function validate_deltas(d)
required = {'dT','dRH','fP','fSW','fWs','BaselineYears','FutureYears','scenario'};
missing = required(~isfield(d,required));
if ~isempty(missing)
    error('Apply_CMIP6_Deltas_To_ERA5_Land:BadDeltas', ...
        'Delta struct is missing: %s',strjoin(missing,', '));
end
for k = 1:5
    if numel(d.(required{k}))~=12
        error('Apply_CMIP6_Deltas_To_ERA5_Land:BadDeltaLength', ...
            'Delta field %s must have 12 monthly values.',required{k});
    end
end
end

function validate_year_map(map,src,passthrough)
validateattributes(map,{'numeric'},{'2d','ncols',2,'integer','finite'});
if numel(unique(map(:,1)))~=size(map,1)
    error('Apply_CMIP6_Deltas_To_ERA5_Land:DuplicateTargetYear', ...
        'YearMap contains a repeated target year.');
end
for k = 1:size(map,1)
    if is_leap(map(k,1))~=is_leap(map(k,2))
        error('Apply_CMIP6_Deltas_To_ERA5_Land:LeapMismatch', ...
            ['Target year %d and source year %d differ in leap status. ' ...
            'Pair leap years with leap years so the hourly axis stays continuous.'], ...
            map(k,1),map(k,2));
    end
end
for y = passthrough(:)'
    row = map(map(:,1)==y,:);
    if isempty(row)
        error('Apply_CMIP6_Deltas_To_ERA5_Land:BadPassthrough', ...
            'Passthrough year %d is not in YearMap.',y);
    end
    if row(1)~=row(2)
        error('Apply_CMIP6_Deltas_To_ERA5_Land:BadPassthrough', ...
            'Passthrough year %d must map to itself.',y);
    end
end
counts = accumarray(findgroups(src.year(:)),1);
if any(counts~=8760 & counts~=8784)
    error('Apply_CMIP6_Deltas_To_ERA5_Land:IncompleteSourceYear', ...
        'The ERA5-Land record contains a partial year; trim it to whole years.');
end
end

function tf = is_leap(y)
tf = (mod(y,4)==0 & mod(y,100)~=0) | mod(y,400)==0;
end

function write_era5_like(f,utc,Lat,Lon,t2m,d2m,sp,u10,v10,tp,ssrd,deltas,opt)
if isfile(f)
    delete(f);
end
n = numel(utc);
origin = datetime(1970,1,1,'TimeZone','UTC');
tvals = seconds(utc-origin);

nccreate(f,'valid_time','Dimensions',{'valid_time',n},'Datatype','double');
ncwrite(f,'valid_time',tvals);
ncwriteatt(f,'valid_time','units','seconds since 1970-01-01 00:00:00');
ncwriteatt(f,'valid_time','calendar','proleptic_gregorian');

nccreate(f,'latitude','Datatype','double');
ncwrite(f,'latitude',Lat);
ncwriteatt(f,'latitude','units','degrees_north');
nccreate(f,'longitude','Datatype','double');
ncwrite(f,'longitude',Lon);
ncwriteatt(f,'longitude','units','degrees_east');

% Units strings must survive the downstream normaliser, which strips
% whitespace and ^*{} before a case-insensitive compare.
spec = { ...
    't2m',t2m,'K'; ...
    'd2m',d2m,'K'; ...
    'sp',sp,'Pa'; ...
    'u10',u10,'m s-1'; ...
    'v10',v10,'m s-1'; ...
    'tp',tp,'m'; ...
    'ssrd',ssrd,'J m-2'};
for k = 1:size(spec,1)
    nccreate(f,spec{k,1},'Dimensions',{'valid_time',n},'Datatype','double');
    ncwrite(f,spec{k,1},spec{k,2});
    ncwriteatt(f,spec{k,1},'units',spec{k,3});
end

ncwriteatt(f,'/','title','ERA5-Land hourly point series with CMIP6 delta perturbation');
ncwriteatt(f,'/','scenario',deltas.scenario);
ncwriteatt(f,'/','delta_models',strjoin(reshape(deltas.models,1,[]),', '));
ncwriteatt(f,'/','delta_baseline',sprintf('%d-%d',deltas.BaselineYears));
ncwriteatt(f,'/','delta_future',sprintf('%d-%d',deltas.FutureYears));
ncwriteatt(f,'/','ramp_reference_year',opt.RampReferenceYear);
ncwriteatt(f,'/','ramp_target_year',opt.RampTargetYear);
ncwriteatt(f,'/','precip_method',opt.PrecipMethod);
ncwriteatt(f,'/','deaccumulation',opt.Deaccumulation);
ncwriteatt(f,'/','created',datestr(now,'yyyy-mm-ddTHH:MM:SS'));
end

function opt = parse_options(args)
opt.YearMap = [];
opt.PassthroughYears = [];
opt.RampReferenceYear = [];
opt.RampTargetYear = [];
opt.PrecipMethod = 'delta';
opt.Deaccumulation = '';
opt.SiteLat = [];
opt.SiteLon = [];
opt.Verbose = true;
if mod(numel(args),2)~=0
    error('Apply_CMIP6_Deltas_To_ERA5_Land:BadOptions', ...
        'Options must be supplied as name-value pairs.');
end
for k = 1:2:numel(args)
    name = validatestring(args{k},{'YearMap','PassthroughYears', ...
        'RampReferenceYear','RampTargetYear','PrecipMethod', ...
        'Deaccumulation','SiteLat','SiteLon','Verbose'});
    opt.(name) = args{k+1};
end
opt.PrecipMethod = validatestring(opt.PrecipMethod,{'delta','qdm'});
if isempty(opt.Deaccumulation)
    error('Apply_CMIP6_Deltas_To_ERA5_Land:MissingDeaccumulation', ...
        ['Deaccumulation is required: ''legacy'' for the basis of ' ...
        'Data_Kantola.mat, ''timeseries'' for the increments as downloaded.']);
end
opt.Deaccumulation = validatestring(opt.Deaccumulation,{'legacy','timeseries'});
if isempty(opt.SiteLat) || isempty(opt.SiteLon)
    error('Apply_CMIP6_Deltas_To_ERA5_Land:MissingSiteCoordinates', ...
        ['SiteLat and SiteLon are required. They are the coordinates the ' ...
        'model runs at, which need not equal the ERA5-Land cell centre that ' ...
        'Lat and Lon identify.']);
end
validateattributes(opt.SiteLat,{'numeric'}, ...
    {'scalar','real','finite','>=',-90,'<=',90});
validateattributes(opt.SiteLon,{'numeric'}, ...
    {'scalar','real','finite','>=',-180,'<=',180});
end
