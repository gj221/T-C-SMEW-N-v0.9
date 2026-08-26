function [spinup_file,transient_file] = Preprocess_Bruce_Field_ERA5( ...
    era5_dir,output_dir)
% Build Bruce Field T&C forcing for the spin-up and transient periods.

if nargin<1 || isempty(era5_dir)
    candidates = { ...
        ['C:\Users\Gregory\OneDrive - Yale University\Documents\' ...
         'Centrifugation\Bruce Farm_Virginia\Whole farm\Data\' ...
         'ERA5_bruce_field_VA_2014-2025'], ...
        ['C:\Users\Gregory\OneDrive - Yale University\Documents\' ...
         'Bruce Farm_Virginia\Whole farm\Data\' ...
         'ERA5_bruce_field_VA_2014-2025']};
    found = find(cellfun(@isfolder,candidates),1);
    if isempty(found)
        error('Preprocess_Bruce_Field_ERA5:MissingInputDirectory', ...
            'The configured Bruce Field ERA5 directory was not found.');
    end
    era5_dir = candidates{found};
end
if nargin<2 || isempty(output_dir)
    output_dir = fileparts(fileparts(mfilename('fullpath')));
end
era5_dir = char(era5_dir);
output_dir = char(output_dir);
if ~isfolder(era5_dir)
    error('Preprocess_Bruce_Field_ERA5:MissingInputDirectory', ...
        'Bruce Field ERA5 directory not found: %s',era5_dir);
end
if ~isfolder(output_dir)
    error('Preprocess_Bruce_Field_ERA5:MissingOutputDirectory', ...
        'Output directory not found: %s',output_dir);
end

Lat = 36.7;
Lon = -78.5;
Zbas = 127.4689;  % USGS 3DEP elevation at the ERA5 grid point.
DeltaGMT = -5;

temporary_file = [tempname,'.mat'];
temporary_cleanup = onCleanup(@() delete_if_present(temporary_file));
Build_TC_Forcing_From_ERA5_Land_Timeseries( ...
    era5_dir,temporary_file,Lat,Lon,Zbas,DeltaGMT);
forcing = load(temporary_file);

utc_date = forcing.Date-DeltaGMT/24;
expected_start = datenum(2014,1,1); %#ok<DATNM>
expected_end = datenum(2026,1,1)-1/24; %#ok<DATNM>
if abs(utc_date(1)-expected_start)>1e-10 || ...
        abs(utc_date(end)-expected_end)>1e-10
    error('Preprocess_Bruce_Field_ERA5:UnexpectedCoverage', ...
        'ERA5 source must cover every UTC hour from 2014 through 2025.');
end

spinup_file = fullfile(output_dir, ...
    'Data_Bruce_Field_Spinup_2014_2023.mat');
transient_file = fullfile(output_dir, ...
    'Data_Bruce_Field_Transient_2024_2025.mat');
write_period(forcing,spinup_file,2014,2024,'spinup');
write_period(forcing,transient_file,2024,2026,'transient');

clear temporary_cleanup
fprintf('Spin-up forcing: %s\n',spinup_file);
fprintf('Transient forcing: %s\n',transient_file);
end

function write_period(forcing,output_file,start_year,end_year,period_name)
utc_date = forcing.Date-forcing.DeltaGMT/24;
period_start = datenum(start_year,1,1); %#ok<DATNM>
period_end = datenum(end_year,1,1); %#ok<DATNM>
selected = utc_date>=period_start & utc_date<period_end;
expected_hours = round((period_end-period_start)*24);
if sum(selected)~=expected_hours
    error('Preprocess_Bruce_Field_ERA5:IncompletePeriod', ...
        '%s forcing does not contain every expected UTC hour.',period_name);
end

period = forcing;
names = fieldnames(period);
for k = 1:numel(names)
    value = period.(names{k});
    if isvector(value) && numel(value)==numel(utc_date)
        period.(names{k}) = value(selected);
    end
end
period.ForcingPeriod = period_name;
period.ForcingUTCStart = period_start;
period.ForcingUTCEnd = period_end-1/24;
period.SourceProduct = 'ERA5-Land hourly time-series';
save(output_file,'-struct','period','-v7');
end

function delete_if_present(file)
if isfile(file)
    delete(file);
end
end
