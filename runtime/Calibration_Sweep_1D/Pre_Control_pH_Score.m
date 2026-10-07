function score = Pre_Control_pH_Score( ...
        dates,bulk_pH,pH_observations)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Score the 2016 control trajectory against the manuscript pre pH.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dates = double(dates(:));
bulk_pH = double(bulk_pH(:));
%% Row 1 is the lone first hour, dated 31 December 2015, and the day of
%% 31 December 2016 ends at 00:00 on 1 January 2017, after a 366-day run.
complete_2016 = (datenum(2016,1,1):datenum(2016,12,30)).'; %#ok<DATNM>
if ~isreal(dates) || ~isreal(bulk_pH) || ...
        numel(bulk_pH) ~= numel(dates) || ...
        ~all(ismember(complete_2016,dates)) || ...
        any(diff(dates) ~= 1)
    error('Bouchet_Joint_Sweep:PreControlWindow', ...
        'The score requires every complete day of 2016 on a daily time axis.');
end
if any(~isfinite(bulk_pH)) || any(bulk_pH <= 0 | bulk_pH >= 14)
    error('Bouchet_Joint_Sweep:InvalidPH', ...
        'The exploratory control trajectory contains invalid bulk pH.');
end
if ~isstruct(pH_observations) || ...
        ~isfield(pH_observations,'control') || ...
        ~all(isfield(pH_observations.control, ...
            {'shallow_0_10','deep_10_30','profile_0_30'}))
    error('Bouchet_Joint_Sweep:BadPHObservations', ...
        'The control bulk-pH observations are incomplete.');
end

shallow = double(pH_observations.control.shallow_0_10(:));
deep = double(pH_observations.control.deep_10_30(:));
profile = double(pH_observations.control.profile_0_30(:));
if ~isreal(shallow) || ~isreal(deep) || ~isreal(profile) || ...
        numel(shallow) ~= 5 || numel(deep) ~= 5 || numel(profile) ~= 5 || ...
        any(~isfinite([shallow;deep;profile]))
    error('Bouchet_Joint_Sweep:BadPHObservations', ...
        'The control bulk-pH observations are invalid.');
end
expected_observed = -log10((10*10.^(-shallow(1)) + ...
    20*10.^(-deep(1)))/30);
if abs(profile(1)-expected_observed) > 64*eps(expected_observed)
    error('Bouchet_Joint_Sweep:PHProfileAggregation', ...
        'The pre 0 to 30 cm observation is not depth aggregated in H space.');
end

pre_mask = dates >= datenum(2016,1,1) & dates < datenum(2016,11,10); %#ok<DATNM>
n_pre_days = nnz(pre_mask);
if n_pre_days ~= 314
    error('Bouchet_Joint_Sweep:PreControlWindow', ...
        'The manuscript pre period must contain exactly 314 complete days.');
end

score.n_days = numel(dates);
score.n_pre_days = n_pre_days;
score.date_start = dates(1);
score.date_end = dates(end);
score.pre_date_end = dates(find(pre_mask,1,'last'));
score.model_bulk_pH_pre = -log10(mean(10.^(-bulk_pH(pre_mask))));
score.observed_bulk_pH_pre_0_10 = shallow(1);
score.observed_bulk_pH_pre_10_30 = deep(1);
score.observed_bulk_pH_pre = profile(1);
score.residual = score.model_bulk_pH_pre-score.observed_bulk_pH_pre;
score.bias = score.residual;
score.rmse = abs(score.residual);
score.residual_0_10 = score.model_bulk_pH_pre-shallow(1);
score.residual_10_30 = score.model_bulk_pH_pre-deep(1);

end
