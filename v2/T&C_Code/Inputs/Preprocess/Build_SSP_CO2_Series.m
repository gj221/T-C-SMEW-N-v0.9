function Ca = Build_SSP_CO2_Series(Date,scenario,varargin)
% Atmospheric CO2 concentration [ppm] on a T&C forcing Date vector.
%
% Ca must match the scenario used for the climate deltas. In this model it
% feeds both photosynthesis and, via ERW_TC.CO2_atm = Ca/22.41, the soil
% carbonate system, so a mismatched pathway biases weathering directly.
%
% Name-value options:
%   ConcentrationFile  (required) CSV with columns year,ppm from the CMIP6
%                            concentration dataset (Meinshausen et al. 2020,
%                            GMD 13:3571). There is no built-in fallback: Ca
%                            drives photosynthesis and, via
%                            ERW_TC.CO2_atm = Ca/22.41, the soil carbonate
%                            system, so the pathway must be explicit rather
%                            than approximated.
%   Verbose            true

narginchk(2,inf);
opt = parse_options(varargin);
validateattributes(Date,{'numeric'},{'vector','real','finite','increasing'});
scenario = lower(strrep(char(scenario),'-',''));

if isempty(opt.ConcentrationFile)
    error('Build_SSP_CO2_Series:MissingConcentrationFile', ...
        ['ConcentrationFile is required for %s. Supply the official CMIP6 ' ...
        'series (Meinshausen et al. 2020, GMD 13:3571) as a CSV with columns ' ...
        'year and ppm. The approximate decadal anchors this function used to ' ...
        'fall back on read 8 to 11 ppm low at 2050 for ssp126 and ssp245, so ' ...
        'they were removed rather than left as a silent default.'],scenario);
end
if ~isfile(opt.ConcentrationFile)
    error('Build_SSP_CO2_Series:MissingFile', ...
        'Concentration file not found: %s',opt.ConcentrationFile);
end
t = readtable(opt.ConcentrationFile);
if ~all(ismember({'year','ppm'},lower(t.Properties.VariableNames)))
    error('Build_SSP_CO2_Series:BadFile', ...
        'Concentration file must have columns named year and ppm.');
end
t.Properties.VariableNames = lower(t.Properties.VariableNames);
% The CMIP6 concentrations are annual MEANS, so each value represents the
% middle of its year rather than 1 January. Anchoring them at the year
% boundary shifts the series half a year early and, on a rising pathway,
% reads high by half the annual increment: +1.3 ppm at 2016 and, because
% ssp585 accelerates to +7.2 ppm yr^-1, +3.6 ppm at 2050.
anchors_year = t.year(:)' + 0.5;
anchors_ppm = t.ppm(:)';

yr = decimal_year(Date(:));
if min(yr) < min(anchors_year)-1e-9 || max(yr) > max(anchors_year)+1e-9
    error('Build_SSP_CO2_Series:OutOfRange', ...
        'Date range %.1f-%.1f falls outside the concentration series %d-%d.', ...
        min(yr),max(yr),min(anchors_year),max(anchors_year));
end

Ca = interp1(anchors_year,anchors_ppm,yr,'pchip');
Ca = Ca(:);
if any(~isfinite(Ca)) || any(Ca<=0)
    error('Build_SSP_CO2_Series:BadInterpolation', ...
        'CO2 interpolation produced non-positive or non-finite values.');
end
if opt.Verbose
    fprintf('CO2 series for %s: %.1f ppm at %.1f to %.1f ppm at %.1f.\n', ...
        scenario,Ca(1),yr(1),Ca(end),yr(end));
end
end

% ---------------------------------------------------------------------

function yr = decimal_year(Date)
dv = datevec(Date);
y = dv(:,1);
start_of_year = datenum([y ones(numel(y),2) zeros(numel(y),3)]);
next_year = datenum([y+1 ones(numel(y),2) zeros(numel(y),3)]);
yr = y + (Date(:)-start_of_year)./(next_year-start_of_year);
end

function opt = parse_options(args)
opt.ConcentrationFile = '';
opt.Verbose = true;
if mod(numel(args),2)~=0
    error('Build_SSP_CO2_Series:BadOptions', ...
        'Options must be supplied as name-value pairs.');
end
for k = 1:2:numel(args)
    name = validatestring(args{k},{'ConcentrationFile','Verbose'});
    opt.(name) = args{k+1};
end
end
