function summary = Kantola_Bulk_pH_Model_Summary(result_file,sample_date,observation_boundary,square)
%% Bulk pH and base saturation of one ERW column (square, 1 by default)
%% at the rows that carry the bulk pH: the initial row and the last row of
%% each ERW call. A result with the ten-minute records reads them there; a
%% multi-square result of package y keeps the same rows as the column's
%% ERW day summaries (ERW_daily, ERW_Daily_Record) and its t = 0 row, and
%% is read from those. Without the records, the exchange-closure check
%% uses each day's largest residual and the sign check the day-end values.
if nargin < 4 || isempty(square); square = 1; end
saved = who('-file',result_file);
S = load(result_file,'Date','ERW');
if ~isfield(S,'ERW') || ~iscell(S.ERW) || numel(S.ERW) < square
    error('Calibrate_Kantola_Bulk_pH:BadResult', ...
        'Result lacks the one-column ERW state: %s',result_file);
end
E = S.ERW{square};
required = {'f_Na','f_K','f_Ca','f_Mg','f_H','f_Al'};
if isfield(E,'tm') || ~ismember('ERW_daily',saved)
    [pH_all,fractions_all,model_time,closure] = from_records(E,required,result_file);
else
    D = load(result_file,'ERW_daily');
    [pH_all,fractions_all,model_time,closure] = ...
        from_daily(E,D.ERW_daily{square},required,result_file);
end
if any(isinf(pH_all)) || any(~isfinite(fractions_all(:))) || ...
        any(fractions_all(:) < -1e-12) || closure > 1e-8
    error('Calibrate_Kantola_Bulk_pH:UnphysicalResult', ...
        'Bulk pH or exchange fractions are unphysical in %s.',result_file);
end

retained = find(isfinite(pH_all));
if numel(retained) < 2 || numel(model_time) < numel(pH_all)
    error('Calibrate_Kantola_Bulk_pH:BadResult', ...
        'Result lacks a usable retained bulk-pH time grid: %s',result_file);
end
date0 = floor(S.Date(1));
model_day = round(double(model_time(retained))/24);
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
summary.exchange_closure_max = closure;
end

function [pH_all,fractions_all,model_time,closure] = from_records(E,required,result_file)
%% Every ten-minute row; the bulk pH is NaN except where it is set.
if ~isfield(E,'bg') || ~isfield(E.bg,'pH_bulk')
    error('Calibrate_Kantola_Bulk_pH:BadResult', ...
        'Result lacks retained bulk pH: %s',result_file);
end
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
closure = max(abs(sum(fractions_all,2)-1));
if ~isfield(E,'tm') || numel(E.tm) < n
    model_time = [];
else
    model_time = E.tm(1:n);
end
end

function [pH_all,fractions_all,model_time,closure] = from_daily(E,D,required,result_file)
%% The t = 0 row, then each ERW call's last row, at the grid time of the row.
if ~isfield(E,'daily_initial') || ~isfield(E,'res') || isempty(D) || ...
        ~isfield(D,'end') || ~isfield(D.end,'pH_bulk') || ~isfield(D,'max')
    error('Calibrate_Kantola_Bulk_pH:BadResult', ...
        'Result lacks retained bulk pH: %s',result_file);
end
B0 = E.daily_initial.bg;
pH_all = [B0.pH_bulk(1,1); D.end.pH_bulk(:,1)];
fractions_all = zeros(numel(pH_all),numel(required));
for k = 1:numel(required)
    if ~isfield(B0,required{k}) || ~isfield(D.end,required{k})
        error('Calibrate_Kantola_Bulk_pH:BadResult', ...
            'Result lacks exchange fraction %s.',required{k});
    end
    fractions_all(:,k) = [B0.(required{k})(1,1); D.end.(required{k})(:,1)];
end
closure = max([abs(sum(fractions_all(1,:))-1); D.max.f_sum(:,1)]);
model_time = ([1; D.rows(:,2)]-1)/E.res;
end
