function test_build_tc_forcing_from_era5_land_timeseries
% ERA5-Land hourly point data must be converted to the T&C forcing contract.

test_dir = fileparts(mfilename('fullpath'));
preprocess_dir = fullfile(fileparts(test_dir),'Inputs','Preprocess');
addpath(preprocess_dir);
path_cleanup = onCleanup(@() rmpath(preprocess_dir));

work_dir = tempname;
mkdir(work_dir);
file_cleanup = onCleanup(@() rmdir(work_dir,'s'));
nc_file = fullfile(work_dir,'era5_land_hourly.nc');
out_file = fullfile(work_dir,'Data_Test.mat');

n = 48;
hours_since_epoch = int64((0:n-1)'*3600);
hour_of_day = mod((0:n-1)',24);
local_hour = mod(hour_of_day-6,24);
Rsw = max(600*sin(pi*(local_hour-6)/12),0);
u10 = 3*ones(n,1);
v10 = 4*ones(n,1);
u10(1) = 0;
v10(1) = 0;

write_variable(nc_file,'valid_time',hours_since_epoch,'seconds since 2020-06-01 00:00:00');
write_variable(nc_file,'t2m',293.15*ones(n,1),'K');
write_variable(nc_file,'d2m',283.15*ones(n,1),'K');
write_variable(nc_file,'sp',101325*ones(n,1),'Pa');
write_variable(nc_file,'u10',u10,'m s**-1');
write_variable(nc_file,'v10',v10,'m s**-1');
write_variable(nc_file,'tp',0.001*ones(n,1),'m');
write_variable(nc_file,'ssrd',Rsw*3600,'J m**-2');

Build_TC_Forcing_From_ERA5_Land_Timeseries( ...
    nc_file,out_file,40,-88,220,-6);

forcing = load(out_file);
required = {'Date','Pr','Ta','Ws','ea','Pre','SAD1','SAD2','SAB1','SAB2', ...
    'N','Tdew','esat','PARB','PARD','Lat','Lon','Zbas','DeltaGMT'};
assert(all(isfield(forcing,required)),'Output forcing fields are incomplete.');
assert(numel(forcing.Date)==n,'Output length does not match the ERA5 time axis.');
assert(abs(forcing.Ws(1)-0.01)<1e-12,'Zero wind was not replaced with 0.01 m s^-1.');
assert(abs(forcing.Ws(2)-5)<1e-12,'Wind components were not combined correctly.');
assert(max(abs(forcing.Pr-1))<1e-12,'Precipitation was not converted from m to mm.');
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
