function Build_TC_Forcing_From_ERA5_Land_Timeseries( ...
    era5_file,output_file,Lat,Lon,Zbas,DeltaGMT)
% Convert one-point ERA5-Land hourly time-series data to T&C forcing.
% Use the CDS time-series product because tp and ssrd are hourly increments.
% Lat and Lon are degrees, Zbas is m a.s.l., and DeltaGMT is a fixed hour offset.

narginchk(6,6);
era5_file = char(era5_file);
output_file = char(output_file);
if ~isfile(era5_file)
    error('Build_TC_Forcing_From_ERA5_Land_Timeseries:MissingInput', ...
        'ERA5-Land NetCDF file not found: %s',era5_file);
end
[output_dir,~,output_ext] = fileparts(output_file);
if ~isempty(output_dir) && ~isfolder(output_dir)
    error('Build_TC_Forcing_From_ERA5_Land_Timeseries:MissingOutputDirectory', ...
        'Output directory not found: %s',output_dir);
end
if ~strcmpi(output_ext,'.mat')
    error('Build_TC_Forcing_From_ERA5_Land_Timeseries:BadOutputExtension', ...
        'Output file must have a .mat extension.');
end
validateattributes(Lat,{'numeric'},{'scalar','real','finite','>=',-90,'<=',90});
validateattributes(Lon,{'numeric'},{'scalar','real','finite','>=',-180,'<=',180});
validateattributes(Zbas,{'numeric'},{'scalar','real','finite'});
validateattributes(DeltaGMT,{'numeric'},{'scalar','real','finite','>=',-12,'<=',14});

info = ncinfo(era5_file);
time_name = find_variable(info,{'valid_time','time'});
t2m_name = find_variable(info,{'t2m','2m_temperature'});
d2m_name = find_variable(info,{'d2m','2m_dewpoint_temperature'});
sp_name = find_variable(info,{'sp','surface_pressure'});
u10_name = find_variable(info,{'u10','10m_u_component_of_wind'});
v10_name = find_variable(info,{'v10','10m_v_component_of_wind'});
tp_name = find_variable(info,{'tp','total_precipitation'});
ssrd_name = find_variable(info,{'ssrd','surface_solar_radiation_downwards'});

require_units(era5_file,t2m_name,{'k','kelvin'});
require_units(era5_file,d2m_name,{'k','kelvin'});
require_units(era5_file,sp_name,{'pa','pascal','pascals'});
require_units(era5_file,u10_name,{'ms-1','m/s'});
require_units(era5_file,v10_name,{'ms-1','m/s'});
require_units(era5_file,tp_name,{'m','metre','meter'});
require_units(era5_file,ssrd_name,{'jm-2','j/m2'});

Date = read_local_date(era5_file,time_name,DeltaGMT);
t2m = read_vector(era5_file,t2m_name);
d2m = read_vector(era5_file,d2m_name);
sp = read_vector(era5_file,sp_name);
u10 = read_vector(era5_file,u10_name);
v10 = read_vector(era5_file,v10_name);
tp = read_vector(era5_file,tp_name);
ssrd = read_vector(era5_file,ssrd_name);

n = numel(Date);
series = {t2m,d2m,sp,u10,v10,tp,ssrd};
if any(cellfun(@numel,series)~=n)
    error('Build_TC_Forcing_From_ERA5_Land_Timeseries:LengthMismatch', ...
        'ERA5 variables do not share one hourly time axis.');
end
if any(~isfinite([t2m;d2m;sp;u10;v10;tp;ssrd]))
    error('Build_TC_Forcing_From_ERA5_Land_Timeseries:MissingValues', ...
        'ERA5 forcing contains NaN or infinite values.');
end
if any(abs(diff(Date)-1/24)>1e-8)
    error('Build_TC_Forcing_From_ERA5_Land_Timeseries:NonHourlyTime', ...
        'ERA5 time must be strictly increasing with one-hour spacing.');
end
if any(t2m<=0) || any(d2m<=0)
    error('Build_TC_Forcing_From_ERA5_Land_Timeseries:BadTemperature', ...
        'ERA5 temperatures must be supplied in kelvin.');
end
if any(sp<=0)
    error('Build_TC_Forcing_From_ERA5_Land_Timeseries:BadPressure', ...
        'ERA5 surface pressure must be positive.');
end
if any(tp < -1e-9)
    error('Build_TC_Forcing_From_ERA5_Land_Timeseries:NegativePrecipitation', ...
        'ERA5 hourly precipitation contains a negative value.');
end
if any(ssrd < -1e-3)
    error('Build_TC_Forcing_From_ERA5_Land_Timeseries:NegativeRadiation', ...
        'ERA5 hourly surface solar radiation contains a negative value.');
end

Ta = t2m - 273.15;
Tdew = d2m - 273.15;
Pre = sp/100;
Ws = sqrt(u10.^2 + v10.^2);
Ws(Ws==0) = 0.01;
Pr = 1000*max(tp,0);
Rsw = max(ssrd,0)/3600;

esat = 610.8*exp(17.27*Ta./(237.3+Ta));
ea = 610.8*exp(17.27*Tdew./(237.3+Tdew));
ea = min(ea,esat);

[~,~,SAD1,SAD2,SAB1,SAB2,PARB,PARD,N] = ...
    Automatic_Radiation_Partition_I( ...
    Date,Lat,Lon,Zbas,DeltaGMT,Pr,Tdew,Rsw,0);
SAD1 = SAD1(:);
SAD2 = SAD2(:);
SAB1 = SAB1(:);
SAB2 = SAB2(:);
PARB = PARB(:);
PARD = PARD(:);
N = N(:);

partition_total = SAD1 + SAD2 + SAB1 + SAB2;
residual = Rsw - partition_total;
nonclosure = abs(residual)>1e-8;
SAD1(nonclosure) = SAD1(nonclosure) + residual(nonclosure)/2;
SAD2(nonclosure) = SAD2(nonclosure) + residual(nonclosure)/2;
if any(~isfinite([SAD1;SAD2;SAB1;SAB2;PARB;PARD;N]))
    error('Build_TC_Forcing_From_ERA5_Land_Timeseries:PartitionFailure', ...
        'Radiation partitioning produced NaN or infinite values.');
end
if any([SAD1;SAD2;SAB1;SAB2;PARB;PARD] < -1e-8)
    error('Build_TC_Forcing_From_ERA5_Land_Timeseries:NegativePartition', ...
        'Radiation partitioning produced a negative flux.');
end
if max(abs(SAD1+SAD2+SAB1+SAB2-Rsw))>1e-8
    error('Build_TC_Forcing_From_ERA5_Land_Timeseries:RadiationNonclosure', ...
        'Partitioned shortwave radiation does not close to ERA5 radiation.');
end

save(output_file,'Date','Pr','Ta','Ws','ea','Pre', ...
    'SAD1','SAD2','SAB1','SAB2','N','Tdew','esat','PARB','PARD', ...
    'Lat','Lon','Zbas','DeltaGMT','-v7');
fprintf('Wrote %s with %d hourly forcing records.\n',output_file,n);
end

function name = find_variable(info,candidates)
available = {info.Variables.Name};
for k = 1:numel(candidates)
    idx = find(strcmpi(available,candidates{k}),1);
    if ~isempty(idx)
        name = available{idx};
        return
    end
end
error('Build_TC_Forcing_From_ERA5_Land_Timeseries:MissingVariable', ...
    'Required ERA5 variable not found. Expected one of: %s', ...
    strjoin(candidates,', '));
end

function require_units(nc_file,name,allowed)
try
    units = ncreadatt(nc_file,name,'units');
catch
    error('Build_TC_Forcing_From_ERA5_Land_Timeseries:MissingUnits', ...
        'NetCDF variable %s has no units attribute.',name);
end
units_normalized = regexprep(char(units),'[\s\^\*\{\}]','');
if ~any(strcmpi(units_normalized,allowed))
    error('Build_TC_Forcing_From_ERA5_Land_Timeseries:UnexpectedUnits', ...
        'Unexpected units for %s: %s',name,char(units));
end
end

function values = read_vector(nc_file,name)
values = squeeze(double(ncread(nc_file,name)));
if ~isvector(values)
    error('Build_TC_Forcing_From_ERA5_Land_Timeseries:NotPointData', ...
        ['Variable %s is not a one-point time series. Download the ERA5-Land ' ...
         'hourly time-series product for one location.'],name);
end
values = values(:);
end

function Date = read_local_date(nc_file,time_name,DeltaGMT)
time_values = read_vector(nc_file,time_name);
try
    units = char(ncreadatt(nc_file,time_name,'units'));
catch
    error('Build_TC_Forcing_From_ERA5_Land_Timeseries:MissingTimeUnits', ...
        'NetCDF time variable has no units attribute.');
end
parts = regexp(strtrim(units), ...
    '^(seconds|hours|days)\s+since\s+(.+)$','tokens','once','ignorecase');
if isempty(parts)
    error('Build_TC_Forcing_From_ERA5_Land_Timeseries:BadTimeUnits', ...
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
    error('Build_TC_Forcing_From_ERA5_Land_Timeseries:BadTimeOrigin', ...
        'Could not parse NetCDF time origin: %s',origin_text);
end
switch lower(parts{1})
    case 'seconds'
        utc_time = origin + seconds(time_values);
    case 'hours'
        utc_time = origin + hours(time_values);
    case 'days'
        utc_time = origin + days(time_values);
end
local_time = utc_time + hours(DeltaGMT);
local_time.TimeZone = '';
Date = datenum(local_time); %#ok<DATNM>
Date = Date(:);
end
