function observations = Load_Sweep_Observations(v2_dir)
obs_dir = fullfile(v2_dir,'Experimental observations', ...
    'Web plot digitizer');
pH_file = fullfile(obs_dir,'Data_Beerling.mat');
n2o_file = fullfile(obs_dir,'figure3a_n2o_digitised.csv');
maize_file = fullfile(obs_dir,'Kantola_maize_extracted_data.mat');
soybean_file = fullfile(obs_dir,'Kantola_soybean_extracted_data.mat');
files = {pH_file,n2o_file,maize_file,soybean_file};
if ~all(cellfun(@isfile,files))
    error('Bouchet_Joint_Sweep:MissingObservations', ...
        'One or more local pH, N2O or carbon observation files are missing.');
end
observations.files = files;
observations.pH = load_pH_observations(pH_file);
observations.base_sat = load_base_sat_observations(pH_file);
observations.n2o = load_n2o_observations(n2o_file);
observations.carbon.maize = load_carbon_crop(maize_file);
observations.carbon.soybean = load_carbon_crop(soybean_file);
end

function obs = load_pH_observations(file)
S = load(file,'pH_Exp');
required = {'X','control_0','control_10','basalt_0','basalt_10'};
if ~isfield(S,'pH_Exp') || ~isstruct(S.pH_Exp) || ...
        ~all(isfield(S.pH_Exp,required))
    error('Bouchet_Joint_Sweep:BadPHObservations', ...
        'The local bulk-pH observation structure is incomplete.');
end
obs.labels = strip(string(S.pH_Exp.X(:)));
if ~isequal(lower(obs.labels),["pre";"2017";"2018";"2019";"2020"])
    error('Bouchet_Joint_Sweep:BadPHObservations', ...
        'The bulk-pH observation periods are not pre and 2017 to 2020.');
end
for treatment = {'control','basalt'}
    shallow = S.pH_Exp.([treatment{1} '_0']);
    deep = S.pH_Exp.([treatment{1} '_10']);
    shallow = double(shallow(:));
    deep = double(deep(:));
    if ~isreal(shallow) || ~isreal(deep) || ...
            numel(shallow) ~= 5 || numel(deep) ~= 5 || ...
            any(~isfinite([shallow;deep])) || ...
            any([shallow;deep] <= 0 | [shallow;deep] >= 14)
        error('Bouchet_Joint_Sweep:BadPHObservations', ...
            'The %s bulk-pH observations are invalid.',treatment{1});
    end
    obs.(treatment{1}).shallow_0_10 = shallow;
    obs.(treatment{1}).deep_10_30 = deep;
    mixed_H = (10*10.^(-shallow)+20*10.^(-deep))/30;
    profile = -log10(mixed_H);
    aggregation_residual = abs(10.^(-profile)-mixed_H);
    if any(~isfinite(profile)) || ...
            any(profile < min([shallow deep],[],2) | ...
                profile > max([shallow deep],[],2)) || ...
            any(aggregation_residual > 32*eps(mixed_H))
        error('Bouchet_Joint_Sweep:PHProfileAggregation', ...
            'The %s 0 to 30 cm hydrogen-activity average is invalid.', ...
            treatment{1});
    end
    obs.(treatment{1}).profile_0_30 = profile;
end
end

function obs = load_n2o_observations(file)
T = readtable(file,'TextType','string');
required = {'series','source','year','day_of_year','n2o_flux_gN_ha_day'};
if ~all(ismember(required,T.Properties.VariableNames))
    error('Bouchet_Joint_Sweep:BadN2OObservations', ...
        'The digitised N2O table is missing a required column.');
end
T = T(T.source == "Energy Farm measurement" & ...
    ismember(T.year,[2017 2018]),:);
treatment = strings(height(T),1);
treatment(T.series == "EF Cont") = "control";
treatment(T.series == "EF ERW") = "basalt";
keep = treatment ~= "";
date_value = datenum(T.year(keep),1,1)+ ...
    round(T.day_of_year(keep))-1; %#ok<DATNM>
flux = double(T.n2o_flux_gN_ha_day(keep));
treatment = treatment(keep);
if isempty(flux) || any(~isfinite(date_value)) || ...
        any(~isfinite(flux)) || any(flux < 0) || ...
        ~all(ismember(["control" "basalt"],unique(treatment)))
    error('Bouchet_Joint_Sweep:BadN2OObservations', ...
        'The matched 2017 to 2018 N2O observations are invalid.');
end
key = strcat(treatment,"_",string(date_value));
[~,~,group] = unique(key,'stable');
obs.date_datenum = accumarray(group,date_value,[],@mean);
obs.flux_gN_ha_d = accumarray(group,flux,[],@mean);
names = strings(max(group),1);
for k = 1:max(group); names(k) = treatment(find(group == k,1)); end
obs.treatment = names;
[obs.date_datenum,order] = sort(obs.date_datenum);
obs.flux_gN_ha_d = obs.flux_gN_ha_d(order);
obs.treatment = obs.treatment(order);
end

function obs = load_carbon_crop(file)
S = load(file);
for treatment = {'control','basalt'}
    for flux = {'GPP','NEE','Re'}
        name = [treatment{1} '_' flux{1}];
        if ~isfield(S,name) || ~isstruct(S.(name)) || ...
                ~isfield(S.(name),'y_daily')
            error('Bouchet_Joint_Sweep:BadCarbonObservations', ...
                'The carbon observation structure %s is missing.',name);
        end
        value = double(S.(name).y_daily(:));
        if ~isreal(value) || numel(value) ~= 366 || any(~isfinite(value))
            error('Bouchet_Joint_Sweep:BadCarbonObservations', ...
                'The carbon observation curve %s is invalid.',name);
        end
        obs.(treatment{1}).(flux{1}) = value;
        obs.(treatment{1}).([flux{1} '_annual_gC_m2']) = sum(value);
    end
end
end

function obs = load_base_sat_observations(file)
% Fig. 1F base saturation, per cent of CEC, digitised beside the pH. The
% 0 to 30 cm value is the soil-mass average of the two intervals (10 and
% 20 cm at one CEC per kg), a linear quantity unlike pH.
S = load(file,'CEC_Exp');
required = {'X','control_0','control_10','basalt_0','basalt_10'};
if ~isfield(S,'CEC_Exp') || ~isstruct(S.CEC_Exp) || ...
        ~all(isfield(S.CEC_Exp,required))
    error('Bouchet_Joint_Sweep:BadBaseSatObservations', ...
        'The local base-saturation observation structure is incomplete.');
end
obs.labels = strip(string(S.CEC_Exp.X(:)));
if ~isequal(lower(obs.labels),["pre";"2017";"2018";"2019";"2020"])
    error('Bouchet_Joint_Sweep:BadBaseSatObservations', ...
        'The base-saturation observation periods are not pre and 2017 to 2020.');
end
for treatment = {'control','basalt'}
    shallow = double(S.CEC_Exp.([treatment{1} '_0'])(:));
    deep = double(S.CEC_Exp.([treatment{1} '_10'])(:));
    if ~isreal(shallow) || ~isreal(deep) || ...
            numel(shallow) ~= 5 || numel(deep) ~= 5 || ...
            any(~isfinite([shallow;deep])) || ...
            any([shallow;deep] <= 0 | [shallow;deep] > 100)
        error('Bouchet_Joint_Sweep:BadBaseSatObservations', ...
            'The %s base-saturation observations are invalid.',treatment{1});
    end
    obs.(treatment{1}).shallow_0_10 = shallow;
    obs.(treatment{1}).deep_10_30 = deep;
    obs.(treatment{1}).profile_0_30 = (10*shallow+20*deep)/30;
end
end
