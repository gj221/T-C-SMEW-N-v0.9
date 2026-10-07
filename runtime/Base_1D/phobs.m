function phobs(res_dir,label)
% August bulk pH of runs 10 (control) and 11 (basalt) against the
% Kantola/Beerling observations. The soil cores were taken inside the
% peak-biomass quadrats, collected in mid August (Beerling et al. 2024,
% Methods; SI Fig. S7). The model is one lumped 0-30 cm box
% (Zbio = 300 mm), so its August mean of daily PH_bulk_d (2016 and
% 2017-2020) is set against BOTH observed depth intervals, 0-10 cm and
% 10-30 cm, separately -- the field has a gradient the box cannot hold, so
% the two intervals bracket where a lumped value can fall. The 0-30 cm
% hydrogen-activity average of the two (the calibration's earlier single
% target) is kept as the third row. Signal = basalt minus control, taken
% day by day. SEM rows are the block-mean standard errors read from the
% vector geometry of Fig. 1E (2026-09-24; NaN where the figure draws no
% bar); the signal's SEM is the quadrature sum of the two arms.
%
% The sampling day within August is not recorded, so beside the mean each
% comparison shows the model's daily range over the month and the
% timing-agnostic "miss": the distance of the observation outside that
% range, zero when some sampling date would match. Its RMS over the five
% years is the score that cannot be moved by the unknown date.
%
% Time axis: PH_bulk_d is the daily series the frame saves (one row per
% day from the start date, as Plot_1D_Then_All_12 reads it), so its dates
% come from Datam(1,1:3); a result whose pH series is as long as the hourly
% Date vector is paired hour by hour instead.
if nargin < 2, label = res_dir; end
O = struct('control',[6.202 6.275 6.310 6.133 5.978; 6.021 6.285 6.316 6.128 6.247], ...
           'basalt', [6.126 6.600 6.663 6.919 6.919; 6.104 6.403 6.566 6.479 6.510]);
SEM = struct('control',[0.063 0.076 0.051 0.073 0.054; 0.113 0.083 0.081 0.071 0.078], ...
             'basalt', [0.074 0.060 0.042 NaN   0.031; 0.071 0.068 0.056 0.048 0.048]);
prof = @(o) -log10((10*10.^(-o(1,:))+20*10.^(-o(2,:)))/30);
rmse = @(e) sqrt(mean(e(isfinite(e)).^2));
tgt = {'0-10 cm','10-30 cm','0-30 H-avg'};
nm = {'control','basalt'}; id = [10 11];
fprintf('\n===== %s =====\n',label);

%% Daily series of both runs on a common axis, then the seasonal statistics.
series = cell(1,2);
for r = 1:2
    d = dir(fullfile(res_dir,sprintf('TC_2D_run%02d_*.mat',id(r))));
    if isempty(d), fprintf('run %02d: no result file\n',id(r)); return; end
    S = load(fullfile(d(1).folder,d(1).name),'PH_bulk_d','Date','Datam');
    ph = S.PH_bulk_d(:);
    if numel(ph) == numel(S.Date)
        dates = S.Date(:);
    else
        if ~isfield(S,'Datam')
            error('phobs:NoDailyAxis','%s has no Datam to date its daily pH.',d(1).name);
        end
        dates = datenum(S.Datam(1,1),S.Datam(1,2),S.Datam(1,3)) + (0:numel(ph)-1)';
    end
    series{r} = struct('ph',ph,'date',dates);
end
n = min(numel(series{1}.ph),numel(series{2}.ph));
if any(abs(series{1}.date(1:n)-series{2}.date(1:n)) > 1e-6)
    error('phobs:AxisMismatch','The two runs are not on the same time axis.');
end
dv = datevec(series{1}.date(1:n));
window = false(n,5);
for k = 1:5
    window(:,k) = dv(:,1)==2015+k & dv(:,2)==8;
end
stat = @(x) [arrayfun(@(k) mean(x(window(:,k))),1:5); ...
             arrayfun(@(k) min(x(window(:,k))),1:5); ...
             arrayfun(@(k) max(x(window(:,k))),1:5)];   % rows: mean, min, max
M = cell(1,2);
for r = 1:2, M{r} = stat(series{r}.ph(1:n)); end
D = stat(series{2}.ph(1:n)-series{1}.ph(1:n));         % daily signal

for r = 1:2
    o = [O.(nm{r}); prof(O.(nm{r}))];
    fprintf('\n%s (run %02d)        2016    2017    2018    2019    2020\n',nm{r},id(r));
    fprintf('  model mean     %7.3f %7.3f %7.3f %7.3f %7.3f\n',M{r}(1,:));
    fprintf('  model min      %7.3f %7.3f %7.3f %7.3f %7.3f\n',M{r}(2,:));
    fprintf('  model max      %7.3f %7.3f %7.3f %7.3f %7.3f\n',M{r}(3,:));
    for t = 1:3
        e = M{r}(1,:)-o(t,:);
        miss = max(0,o(t,:)-M{r}(3,:))+max(0,M{r}(2,:)-o(t,:));
        fprintf('  obs %-10s %7.3f %7.3f %7.3f %7.3f %7.3f\n',tgt{t},o(t,:));
        if t < 3
            fprintf('    SEM          %7.3f %7.3f %7.3f %7.3f %7.3f\n',SEM.(nm{r})(t,:));
        end
        fprintf('    resid        %7.3f %7.3f %7.3f %7.3f %7.3f  RMSE %.3f bias %+.3f\n', ...
            e,rmse(e),mean(e(isfinite(e))));
        fprintf('    miss         %7.3f %7.3f %7.3f %7.3f %7.3f  RMS  %.3f\n',miss,rmse(miss));
    end
    e = [M{r}(1,:)-o(1,:), M{r}(1,:)-o(2,:)];
    miss = [max(0,o(1,:)-M{r}(3,:))+max(0,M{r}(2,:)-o(1,:)), ...
            max(0,o(2,:)-M{r}(3,:))+max(0,M{r}(2,:)-o(2,:))];
    fprintf('    pooled over the two intervals              RMSE %.3f  miss RMS %.3f\n', ...
        rmse(e),rmse(miss));
end

do = [O.basalt-O.control; prof(O.basalt)-prof(O.control)];
fprintf('\nbasalt minus control (daily)\n');
fprintf('  model mean     %7.3f %7.3f %7.3f %7.3f %7.3f\n',D(1,:));
fprintf('  model min      %7.3f %7.3f %7.3f %7.3f %7.3f\n',D(2,:));
fprintf('  model max      %7.3f %7.3f %7.3f %7.3f %7.3f\n',D(3,:));
for t = 1:3
    e = D(1,:)-do(t,:);
    miss = max(0,do(t,:)-D(3,:))+max(0,D(2,:)-do(t,:));
    fprintf('  obs %-10s %7.3f %7.3f %7.3f %7.3f %7.3f\n',tgt{t},do(t,:));
    if t < 3
        fprintf('    SEM          %7.3f %7.3f %7.3f %7.3f %7.3f\n', ...
            sqrt(SEM.basalt(t,:).^2+SEM.control(t,:).^2));
    end
    fprintf('    resid        %7.3f %7.3f %7.3f %7.3f %7.3f  RMSE %.3f  2020 model/obs %.2f\n', ...
        e,rmse(e),D(1,5)/do(t,5));
    fprintf('    miss         %7.3f %7.3f %7.3f %7.3f %7.3f  RMS  %.3f\n',miss,rmse(miss));
end
e = [D(1,:)-do(1,:), D(1,:)-do(2,:)];
miss = [max(0,do(1,:)-D(3,:))+max(0,D(2,:)-do(1,:)), ...
        max(0,do(2,:)-D(3,:))+max(0,D(2,:)-do(2,:))];
fprintf('    pooled over the two intervals              RMSE %.3f  miss RMS %.3f\n', ...
    rmse(e),rmse(miss));
end
