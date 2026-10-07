function summary = Kantola_Bulk_pH_Model_Summary(result_file,sample_date,observation_boundary)
S = load(result_file,'Date','ERW');
if ~isfield(S,'ERW') || ~iscell(S.ERW) || isempty(S.ERW)
    error('Calibrate_Kantola_Bulk_pH:BadResult', ...
        'Result lacks the one-column ERW state: %s',result_file);
end
E = S.ERW{1};
if ~isfield(E,'bg') || ~isfield(E.bg,'pH_bulk')
    error('Calibrate_Kantola_Bulk_pH:BadResult', ...
        'Result lacks retained bulk pH: %s',result_file);
end
required = {'f_Na','f_K','f_Ca','f_Mg','f_H','f_Al'};
for k = 1:numel(required)
    if ~isfield(E.bg,required{k})
        error('Calibrate_Kantola_Bulk_pH:BadResult', ...
            'Result lacks exchange fraction %s.',required{k});
    end
end

n = size(E.bg.pH_bulk,1);
if isfield(E,'dm_end') && ~isempty(E.dm_end)
    n = min(n,E.dm_end);
end
pH_all = E.bg.pH_bulk(1:n,1);
fractions_all = zeros(n,numel(required));
for k = 1:numel(required)
    value = E.bg.(required{k});
    fractions_all(:,k) = value(1:n,1);
end
site_sum = sum(fractions_all,2);
if any(isinf(pH_all)) || any(~isfinite(fractions_all(:))) || ...
        any(fractions_all(:) < -1e-12) || ...
        max(abs(site_sum-1)) > 1e-8
    error('Calibrate_Kantola_Bulk_pH:UnphysicalResult', ...
        'Bulk pH or exchange fractions are unphysical in %s.',result_file);
end

retained = find(isfinite(pH_all));
if numel(retained) < 2 || ~isfield(E,'tm') || numel(E.tm) < n
    error('Calibrate_Kantola_Bulk_pH:BadResult', ...
        'Result lacks a usable retained bulk-pH time grid: %s',result_file);
end
date0 = floor(S.Date(1));
model_day = round(double(E.tm(retained))/24);
if any(~isfinite(model_day)) || any(diff(model_day) < 0)
    error('Calibrate_Kantola_Bulk_pH:BadResult', ...
        'Retained bulk-pH rows are not ordered on the ERW time grid: %s',result_file);
end
keep = [diff(model_day) > 0; true];
retained = retained(keep);
model_day = model_day(keep);
daily_date = date0+model_day;
pH_bulk = pH_all(retained);
fractions = fractions_all(retained,:);
if any(pH_bulk <= 0) || any(pH_bulk >= 9)
    error('Calibrate_Kantola_Bulk_pH:UnphysicalResult', ...
        'Bulk pH or exchange fractions are unphysical in %s.',result_file);
end
pre = daily_date < datenum(2016,11,10);
if ~any(pre) || sample_date > daily_date(end)
    error('Calibrate_Kantola_Bulk_pH:ShortResult', ...
        'Result does not span the pre-application and target periods.');
end
sample = find(daily_date == sample_date,1,'last');
if isempty(sample) || sample <= 1
    error('Calibrate_Kantola_Bulk_pH:ShortResult', ...
        'Result lacks the row preceding the second-application boundary.');
end
base_sat_pct = 100*sum(fractions(:,1:4),2);

summary.pre_bulk_pH = -log10(mean(10.^(-pH_bulk(pre))));
summary.pre_base_sat_pct = mean(base_sat_pct(pre));
summary.pre_application_bulk_pH = pH_bulk(sample-1);
summary.post_application_bulk_pH = pH_bulk(sample);
summary.pre_application_base_sat_pct = base_sat_pct(sample-1);
summary.post_application_base_sat_pct = base_sat_pct(sample);
if strcmp(observation_boundary,'pre_application')
    summary.sample_bulk_pH = summary.pre_application_bulk_pH;
    summary.sample_base_sat_pct = summary.pre_application_base_sat_pct;
else
    summary.sample_bulk_pH = summary.post_application_bulk_pH;
    summary.sample_base_sat_pct = summary.post_application_base_sat_pct;
end
summary.exchange_closure_max = max(abs(site_sum-1));
end
