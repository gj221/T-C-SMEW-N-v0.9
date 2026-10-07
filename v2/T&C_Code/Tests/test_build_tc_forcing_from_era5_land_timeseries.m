function test_build_tc_forcing_from_era5_land_timeseries
% ERA5-Land hourly point data must be converted to the T&C forcing contract.

test_dir = fileparts(mfilename('fullpath'));
preprocess_dir = fullfile(fileparts(test_dir),'Inputs','Preprocess');
addpath(preprocess_dir);
path_cleanup = onCleanup(@() rmpath(preprocess_dir));

work_dir = tempname;
mkdir(work_dir);
file_cleanup = onCleanup(@() rmdir(work_dir,'s'));
out_file = fullfile(work_dir,'Data_Test.mat');
temperature_file = fullfile(work_dir,'era5_temperature.nc');
pressure_file = fullfile(work_dir,'era5_pressure_precipitation.nc');
radiation_file = fullfile(work_dir,'era5_radiation.nc');
wind_file = fullfile(work_dir,'era5_wind.nc');

n = 48;
first_hour = (datenum(2020,6,1)-datenum(1970,1,1))*24; %#ok<DATNM>
hours_since_epoch = int64(first_hour+(0:n-1)');
hour_of_day = mod((0:n-1)',24);
local_hour = mod(hour_of_day-6,24);
Rsw = max(600*sin(pi*(local_hour-6)/12),0);
u10 = 3*ones(n,1);
v10 = 4*ones(n,1);
u10(1) = 0;
v10(1) = 0;
tp = 0.001*ones(n,1);
tp(1) = -5e-8;
ssrd = Rsw*3600;
ssrd(1) = -2;
Rsw(1) = 0;

initialize_file(temperature_file,hours_since_epoch,40,-88);
initialize_file(pressure_file,hours_since_epoch,40,-88);
initialize_file(radiation_file,hours_since_epoch,40,-88);
initialize_file(wind_file,hours_since_epoch,40,-88);
write_variable(temperature_file,'t2m',293.15*ones(n,1),'K');
write_variable(temperature_file,'d2m',283.15*ones(n,1),'K');
write_variable(pressure_file,'sp',101325*ones(n,1),'Pa');
write_variable(pressure_file,'tp',tp,'m');
write_variable(radiation_file,'ssrd',ssrd,'J m**-2');
write_variable(wind_file,'u10',u10,'m s**-1');
write_variable(wind_file,'v10',v10,'m s**-1');

Build_TC_Forcing_From_ERA5_Land_Timeseries( ...
    work_dir,out_file,40,-88,220,-6);

forcing = load(out_file);
required = {'Date','Pr','Ta','Ws','ea','Pre','SAD1','SAD2','SAB1','SAB2', ...
    'N','Tdew','esat','PARB','PARD','Lat','Lon','Zbas','DeltaGMT'};
assert(all(isfield(forcing,required)),'Output forcing fields are incomplete.');
assert(numel(forcing.Date)==n,'Output length does not match the ERA5 time axis.');
assert(abs(forcing.Ws(1)-0.01)<1e-12,'Zero wind was not replaced with 0.01 m s^-1.');
assert(abs(forcing.Ws(2)-5)<1e-12,'Wind components were not combined correctly.');
assert(forcing.Pr(1)==0,'Negative precipitation packing noise was not clipped.');
assert(max(abs(forcing.Pr(2:end)-1))<1e-12, ...
    'Precipitation was not converted from m to mm.');
assert(max(abs(forcing.Ta-20))<1e-12,'Air temperature was not converted from K to degC.');
assert(max(abs(forcing.Tdew-10))<1e-12,'Dew point was not converted from K to degC.');
assert(max(abs(forcing.Pre-1013.25))<1e-12,'Pressure was not converted from Pa to hPa.');
assert(all(forcing.ea<=forcing.esat),'Actual vapour pressure exceeds saturation.');
assert(max(abs(forcing.SAD1+forcing.SAD2+forcing.SAB1+forcing.SAB2-Rsw))<1e-8, ...
    'Partitioned shortwave radiation does not close to ERA5 shortwave radiation.');
assert(isequal(datevec(forcing.Date(1)),[2020 5 31 18 0 0]), ...
    'ERA5 UTC time was not converted to local standard time.');

clear file_cleanup path_cleanup
fprintf('PASS: ERA5-Land hourly forcing conversion and radiation partitioning.\n');
end

function initialize_file(nc_file,time_values,Lat,Lon)
write_variable(nc_file,'valid_time',time_values,'hours since 1970-01-01');
nccreate(nc_file,'latitude','Dimensions',{'point',1},'Datatype','double');
ncwrite(nc_file,'latitude',Lat);
ncwriteatt(nc_file,'latitude','units','degrees_north');
nccreate(nc_file,'longitude','Dimensions',{'point',1},'Datatype','double');
ncwrite(nc_file,'longitude',Lon);
ncwriteatt(nc_file,'longitude','units','degrees_east');
end

function write_variable(nc_file,name,values,units)
if isfile(nc_file)
    nccreate(nc_file,name,'Dimensions',{'time',numel(values)},'Datatype',class(values));
else
    nccreate(nc_file,name,'Dimensions',{'time',numel(values)}, ...
        'Datatype',class(values),'Format','netcdf4');
end
ncwrite(nc_file,name,values);
ncwriteatt(nc_file,name,'units',units);
end
