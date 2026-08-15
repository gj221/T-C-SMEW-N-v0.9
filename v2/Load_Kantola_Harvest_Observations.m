function T = Load_Kantola_Harvest_Observations(file)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Summarise block-level Kantola harvest carbon observations.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 1 || isempty(file)
    here = fileparts(mfilename('fullpath'));
    file = fullfile(here,'Experimental observations','Kantola carbon', ...
        'maize_biomass_ERW_20.csv');
end
if ~(ischar(file) || (isstring(file) && isscalar(file)))
    error('Load_Kantola_Harvest_Observations:BadFile', ...
        'Observation file must be a character vector or string scalar.');
end
file = char(file);
if ~isfile(file)
    error('Load_Kantola_Harvest_Observations:MissingFile', ...
        'Observation file not found: %s',file);
end

raw = readtable(file,'TextType','string','VariableNamingRule','preserve');
required = {'Year','Crop','Block','Treatment','Aboveground Peak Biomass', ...
    'Belowground Peak Biomass','Harvest Biomass (Grain)'};
if ~all(ismember(required,raw.Properties.VariableNames))
    error('Load_Kantola_Harvest_Observations:BadSchema', ...
        'Harvest observation table does not contain the required columns.');
end
harvest = raw.('Harvest Biomass (Grain)');
if ~isnumeric(raw.Year) || ~isnumeric(raw.Block) || ~isnumeric(harvest) || ...
        any(~isfinite(raw.Year) | raw.Year ~= fix(raw.Year)) || ...
        any(~isfinite(raw.Block) | raw.Block ~= fix(raw.Block)) || ...
        any(~isfinite(harvest) | harvest < 0)
    error('Load_Kantola_Harvest_Observations:BadValues', ...
        'Harvest years, blocks or carbon values are invalid.');
end
crop = lower(strtrim(string(raw.Crop)));
treatment = lower(strtrim(string(raw.Treatment)));
if any(ismissing(crop) | ~ismember(crop,["maize" "soybean"])) || ...
        any(ismissing(treatment) | ~ismember(treatment,["control" "basalt"]))
    error('Load_Kantola_Harvest_Observations:BadValues', ...
        'Harvest crop or treatment labels are invalid.');
end

[group,Year,Crop,Treatment] = findgroups(raw.Year,crop,treatment);
SampleSize = splitapply(@numel,harvest,group);
if any(SampleSize < 2)
    error('Load_Kantola_Harvest_Observations:BadValues', ...
        'At least two independent blocks are required per harvest group.');
end
duplicate_block = splitapply(@(x) numel(unique(x)) ~= numel(x),raw.Block,group);
if any(duplicate_block)
    error('Load_Kantola_Harvest_Observations:DuplicateObservation', ...
        'Harvest blocks must be unique within each treatment-year group.');
end
Mean_gC_m2 = splitapply(@mean,harvest,group);
SEM_gC_m2 = splitapply(@(x) std(x,0)/sqrt(numel(x)),harvest,group);
MetricID = repmat("grain_c",numel(Year),1);
Units = repmat("g C m-2",numel(Year),1);
[~,source_name,source_ext] = fileparts(file);
SourceFile = repmat(string([source_name source_ext]),numel(Year),1);
SourceKind = repmat("raw_block_data",numel(Year),1);
T = table(MetricID,Year,Crop,Treatment,Mean_gC_m2,SEM_gC_m2,Units, ...
    SampleSize,SourceFile,SourceKind);

end
