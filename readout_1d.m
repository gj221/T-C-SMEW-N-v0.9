function readout_1d(sim_dir)
%% Read-out of the five-year 1-D runs of this root: the campaign set
%% (Base_1D, Base_1D_case8) and the lower and upper sets of the uncertainty
%% ensemble (Uncertainty_1D, Uncertainty_1D_case8), whichever have finished.
%% For each set it prints the August means of the bulk and pore pH and of
%% the base saturation, the alkalinity the rock released and the solid
%% calcite left, for the control (run 10), the basalt (run 11), the 50 t/ha
%% calcite (run 12) and the 5 t/ha calcite (run 12, workbook case 8), then
%% the basalt-minus-control effect beside the field's. It only reads.
%%
%% August is the month the soil was sampled. The field rows are the 0-30 cm
%% hydrogen-activity average of the two sampled intervals (phobs.m holds
%% the values and their source) and the base-saturation gain read from
%% Beerling et al. (2024), Fig. 1F.
root = fileparts(mfilename('fullpath'));
if nargin < 1 || isempty(sim_dir)
    sim_dir = fullfile(root,'v2','SIM','Subfield_12');
end
years = 2016:2020;
sets = { ...
    'campaign set',fullfile(sim_dir,'Base_1D'),fullfile(sim_dir,'Base_1D_case8'); ...
    'lower set',fullfile(sim_dir,'Uncertainty_1D','lower'),fullfile(sim_dir,'Uncertainty_1D_case8','lower'); ...
    'upper set',fullfile(sim_dir,'Uncertainty_1D','upper'),fullfile(sim_dir,'Uncertainty_1D_case8','upper')};
O = struct('control',[6.202 6.275 6.310 6.133 5.978; 6.021 6.285 6.316 6.128 6.247], ...
           'basalt', [6.126 6.600 6.663 6.919 6.919; 6.104 6.403 6.566 6.479 6.510]);
profile = @(o) -log10((10*10.^(-o(1,:))+20*10.^(-o(2,:)))/30);
field_control = profile(O.control);
field_basalt = profile(O.basalt);
field_gain = [NaN 2.4 4.5 5.9 5.9];        % base-saturation points, basalt minus control
fprintf('Field, 0-30 cm      control pH %s\n',sprintf(' %6.2f',field_control));
fprintf('                    basalt pH  %s\n',sprintf(' %6.2f',field_basalt));
fprintf('                    effect     %s\n',sprintf(' %+6.2f',field_basalt-field_control));
fprintf('                    base saturation gain %s\n',sprintf(' %6.1f',field_gain));
for k = 1:size(sets,1)
    fprintf('\n===== %s =====\n',sets{k,1});
    rows = { ...
        'control',sets{k,2},10; 'basalt',sets{k,2},11; ...
        'calcite 50 t/ha',sets{k,2},12; 'calcite 5 t/ha',sets{k,3},12};
    A = cell(size(rows,1),1);
    for r = 1:size(rows,1)
        A{r} = august(rows{r,2},rows{r,3},years);
        if isempty(A{r})
            fprintf('%-16s not finished (%s)\n',rows{r,1},rows{r,2});
            continue
        end
        fprintf('%-16s alpha %.4g, KNmaxd %.4g, %d days\n',rows{r,1}, ...
            A{r}.alpha,A{r}.knmaxd,A{r}.n_days);
        fprintf('   %-28s %s\n','year',sprintf(' %7d',years));
        fprintf('   %-28s %s\n','bulk pH, August',sprintf(' %7.3f',A{r}.bulk));
        fprintf('   %-28s %s\n','pore pH, August',sprintf(' %7.3f',A{r}.pore));
        fprintf('   %-28s %s\n','base saturation %, August',sprintf(' %7.2f',A{r}.base));
        fprintf('   %-28s %s\n','released by 31 Aug, eq m-2',sprintf(' %7.2f',A{r}.released));
        fprintf('   %-28s %s\n','solid calcite, mol m-2',sprintf(' %7.2f',A{r}.calcite));
        fprintf('   %-28s %s\n','Mg carbonate, mol m-2',sprintf(' %7.2f',A{r}.mgco3));
        fprintf('   bulk pH over the run: min %.3f, max %.3f, all finite %d\n', ...
            A{r}.bulk_min,A{r}.bulk_max,A{r}.finite);
    end
    if ~isempty(A{1}) && ~isempty(A{2})
        fprintf('%-16s\n','basalt effect');
        fprintf('   %-28s %s\n','bulk pH, model',sprintf(' %+7.3f',A{2}.bulk-A{1}.bulk));
        fprintf('   %-28s %s\n','bulk pH, field',sprintf(' %+7.3f',field_basalt-field_control));
        fprintf('   %-28s %s\n','base saturation, model',sprintf(' %7.2f',A{2}.base-A{1}.base));
        fprintf('   %-28s %s\n','base saturation, field',sprintf(' %7.2f',field_gain));
        fprintf('   %-28s %s\n','control pH, model - field',sprintf(' %+7.3f',A{1}.bulk-field_control));
    end
end
end

function a = august(res_dir,run_id,years)
%% August statistics of one result, or [] when the run has no result yet.
a = [];
listing = dir(fullfile(res_dir,sprintf('TC_2D_run%02d_*.mat',run_id)));
listing = listing(~[listing.isdir]);
if isempty(listing)
    return
end
S = load(fullfile(listing(1).folder,listing(1).name), ...
    'PH_bulk_d','PH_pore_d','Datam','ERW','run_provenance');
bulk = S.PH_bulk_d(:);
pore = S.PH_pore_d(:);
n_d = numel(bulk);
E = S.ERW{1};
bg = E.bg;
n_t = numel(bg.pH);
per_day = round(n_t/n_d);
t0 = datenum(S.Datam(1,1),S.Datam(1,2),S.Datam(1,3)); %#ok<DATNM>
day = datevec(t0+(0:n_d-1)');
step = datevec(t0+(1:n_t)'/per_day);
base = 100*(bg.f_Ca(:)+bg.f_Mg(:)+bg.f_K(:)+bg.f_Na(:));
released = cumsum(bg.dAlk_sil(:)+bg.dAlk_carb(:))/1e6;
calcite = max(bg.CaCO3(:),0)/1e6;
names = strtrim(string(E.var.mineral(:)'));
kc = find(strcmpi(names,"calcite"),1);
if ~isempty(kc) && size(bg.M_min,2) >= kc
    calcite = calcite+max(bg.M_min(:,kc),0)/E.const.MM_min(kc)/1e6;
end
mgco3 = max(bg.MgCO3(:),0)/1e6;
a.bulk = nan(1,numel(years)); a.pore = a.bulk; a.base = a.bulk;
a.released = a.bulk; a.calcite = a.bulk; a.mgco3 = a.bulk;
for y = 1:numel(years)
    in_day = day(:,1) == years(y) & day(:,2) == 8;
    in_step = step(:,1) == years(y) & step(:,2) == 8;
    if ~any(in_day) || ~any(in_step)
        continue
    end
    a.bulk(y) = mean(bulk(in_day));
    a.pore(y) = mean(pore(in_day));
    a.base(y) = mean(base(in_step));
    last = find(in_step,1,'last');
    a.released(y) = released(last);
    a.calcite(y) = calcite(last);
    a.mgco3(y) = mgco3(last);
end
a.bulk_min = min(bulk); a.bulk_max = max(bulk);
a.finite = all(isfinite(bulk)) && all(isfinite(pore));
a.n_days = n_d;
a.alpha = NaN; a.knmaxd = NaN;
if isfield(S,'run_provenance')
    p = S.run_provenance;
    if isfield(p,'AppeloAlpha'); a.alpha = p.AppeloAlpha(1); end
    if isfield(p,'KNmaxd'); a.knmaxd = p.KNmaxd(1); end
end
end
