function ok = precheck_1d(root,original_dir,ensemble_dir)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The acceptance checks of the uncertainty ensemble's five-year 1-D runs
%% (UNCERTAINTY_BOUNDS_REVISION_2026-10-09.md, section 3), printed as a
%% table with PASS or FAIL per check, before any 2-D or 35-year array is
%% submitted. It only reads. From the checkout root:
%%
%%   matlab -batch "precheck_1d"
%%   precheck_1d(root)                           % another checkout
%%   precheck_1d(root,original_dir,ensemble_dir) % other result folders
%%
%% The original set is the campaign's own Base_1D result (runs 10 to 12);
%% the lower and upper sets sit under Uncertainty_1D/lower and /upper.
%% For the control (run 10) of each moved set:
%%   N2O 2016 (unfertilised) below 5 kg N/ha; 2016-2020 within 8-50;
%%   peak daily N2O below 150 g N/ha/day;
%%   NO3 leaching below 100 kg N/ha in every year;
%%   NEE 2016 within +/-150 g C/m2 of the original's;
%%   bulk pH: the lower set below and the upper set above the original in
%%   every August mean, neither below 5.4 nor above 6.4;
%%   basalt effect (run 11 - run 10) between +0.1 and +0.8 in August 2020;
%%   calcite 50 t/ha (run 12) at the carbonate ceiling, 7.9-8.3, in
%%   every August after the application.
%% Daily rows follow Summarize_Uncertainty: row r is dated start + r - 3,
%% rows 1 and 2 are the spin-up row and the extra BIOGEO day. Lives at the
%% package root, outside every campaign stamp.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 1 || isempty(root); root = pwd; end
if nargin < 2 || isempty(original_dir)
    original_dir = fullfile(root,'v2','SIM','Subfield_12','Base_1D');
end
if nargin < 3 || isempty(ensemble_dir)
    ensemble_dir = fullfile(root,'v2','SIM','Subfield_12','Uncertainty_1D');
end
sets = {'original',original_dir; 'lower',fullfile(ensemble_dir,'lower'); ...
    'upper',fullfile(ensemble_dir,'upper')};
runs = [10 11 12];
names = {'control','basalt','calcite 50 t/ha'};
years = 2016:2020;

%% Read every result that exists.
R = cell(size(sets,1),numel(runs));
for s = 1:size(sets,1)
    for r = 1:numel(runs)
        R{s,r} = read_result(sets{s,2},runs(r),years);
        if isempty(R{s,r})
            fprintf('%-8s run %d: no result in %s\n',sets{s,1},runs(r),sets{s,2});
        else
            fprintf('%-8s run %d: %d days, alpha %.4g, KNmaxd %.4g (%s)\n',sets{s,1},runs(r), ...
                R{s,r}.n_days,R{s,r}.alpha,R{s,r}.knmaxd,R{s,r}.file);
        end
    end
end

%% Per set and run: the yearly table.
for s = 1:size(sets,1)
    for r = 1:numel(runs)
        a = R{s,r};
        if isempty(a); continue; end
        fprintf('\n===== %s, run %d (%s) =====\n',sets{s,1},runs(r),names{r});
        fprintf('   %-30s %s\n','year',sprintf(' %8d',years));
        fprintf('   %-30s %s\n','N2O, kg N/ha',sprintf(' %8.2f',a.n2o_kg_ha));
        fprintf('   %-30s %s\n','peak daily N2O, g N/ha/day',sprintf(' %8.1f',a.n2o_peak_g_ha_d));
        fprintf('   %-30s %s\n','NO3 leached, kg N/ha',sprintf(' %8.1f',a.no3_kg_ha));
        fprintf('   %-30s %s\n','NEE, g C/m2',sprintf(' %+8.1f',a.nee_g_m2));
        fprintf('   %-30s %s\n','Rh, g C/m2',sprintf(' %8.1f',a.rh_g_m2));
        fprintf('   %-30s %s\n','bulk pH, August mean',sprintf(' %8.3f',a.ph_aug));
        fprintf('   %-30s %s\n','bulk pH, year minimum',sprintf(' %8.3f',a.ph_min));
        fprintf('   %-30s %s\n','bulk pH, year maximum',sprintf(' %8.3f',a.ph_max));
        fprintf('   %-30s %s\n','released by 31 Aug, eq m-2',sprintf(' %8.2f',a.released_aug));
    end
end

%% The checks.
checks = cell(0,4);
for s = 2:3
    tag = sets{s,1};
    if strcmp(tag,'lower'); want = -1; else; want = 1; end
    c = R{s,1}; o = R{1,1}; b = R{s,2}; ob = R{1,2}; k = R{s,3};
    if ~isempty(c)
        y16 = years == 2016;
        checks(end+1,:) = {tag,'N2O 2016 < 5 kg N/ha',c.n2o_kg_ha(y16) < 5, ...
            sprintf('%.2f',c.n2o_kg_ha(y16))}; %#ok<AGROW>
        total = sum(c.n2o_kg_ha,'omitnan');
        checks(end+1,:) = {tag,'N2O 2016-2020 within 8-50 kg N/ha',total >= 8 && total <= 50, ...
            sprintf('%.2f',total)}; %#ok<AGROW>
        peak = max(c.n2o_peak_g_ha_d);
        checks(end+1,:) = {tag,'peak daily N2O < 150 g N/ha/day',peak < 150,sprintf('%.1f',peak)}; %#ok<AGROW>
        checks(end+1,:) = {tag,'NO3 leached < 100 kg N/ha every year',all(c.no3_kg_ha < 100), ...
            sprintf('max %.1f',max(c.no3_kg_ha))}; %#ok<AGROW>
        if ~isempty(o)
            d = c.nee_g_m2(y16)-o.nee_g_m2(y16);
            checks(end+1,:) = {tag,'NEE 2016 within +/-150 g C/m2 of the original',abs(d) <= 150, ...
                sprintf('%+.1f (original %+.1f)',d,o.nee_g_m2(y16))}; %#ok<AGROW>
            side = sign(c.ph_aug-o.ph_aug);
            checks(end+1,:) = {tag,'control August pH on its side of the original every year', ...
                all(side(isfinite(side)) == want),sprintf('%+.3f ',c.ph_aug-o.ph_aug)}; %#ok<AGROW>
        end
        checks(end+1,:) = {tag,'control August pH within 5.4-6.4', ...
            all(c.ph_aug >= 5.4 & c.ph_aug <= 6.4),sprintf('%.3f ',c.ph_aug)}; %#ok<AGROW>
    end
    if ~isempty(c) && ~isempty(b)
        effect = b.ph_aug-c.ph_aug;
        e20 = effect(years == 2020);
        checks(end+1,:) = {tag,'basalt effect +0.1 to +0.8 in August 2020', ...
            e20 >= 0.1 && e20 <= 0.8,sprintf('%+.3f (all: %s)',e20,sprintf('%+.2f ',effect))}; %#ok<AGROW>
        if ~isempty(ob) && ~isempty(o)
            oe = ob.ph_aug-o.ph_aug;
            later = years >= 2017;
            checks(end+1,:) = {tag,'basalt August pH on its side of the original every year after 2016', ...
                all(sign(b.ph_aug(later)-ob.ph_aug(later)) == want), ...
                sprintf('%+.3f ',b.ph_aug-ob.ph_aug)}; %#ok<AGROW>
            fprintf('\nbasalt effect, August: original %s | %s %s\n',sprintf('%+.2f ',oe),tag,sprintf('%+.2f ',effect));
        end
    end
    if ~isempty(k)
        after = years >= 2017;
        checks(end+1,:) = {tag,'calcite 50 t/ha at the carbonate ceiling (7.9-8.3) every August after application', ...
            all(k.ph_aug(after) >= 7.9 & k.ph_aug(after) <= 8.3),sprintf('%.2f ',k.ph_aug)}; %#ok<AGROW>
    end
end
fprintf('\n%-7s %-80s %-5s %s\n','set','check','','value');
ok = true;
for i = 1:size(checks,1)
    verdict = 'FAIL'; if checks{i,3}; verdict = 'PASS'; end
    ok = ok && checks{i,3};
    fprintf('%-7s %-80s %-5s %s\n',checks{i,1},checks{i,2},verdict,checks{i,4});
end
if isempty(checks)
    ok = false;
    fprintf('No check could be made: no moved-set result was found.\n');
elseif ok
    fprintf('\nPRECHECK PASS: the 2-D and 35-year arrays may be submitted.\n');
else
    fprintf('\nPRECHECK FAIL: do not submit the 2-D or 35-year arrays; see the note, section 2c (fallback x/1.5).\n');
end
end

function a = read_result(res_dir,run_id,years)
a = [];
listing = dir(fullfile(res_dir,sprintf('TC_2D_run%02d_*.mat',run_id)));
listing = listing(~[listing.isdir]);
if isempty(listing); return; end
file = fullfile(listing(1).folder,listing(1).name);
S = load(file,'PH_bulk_d','Datam','Domain','R_litter','R_microbe','R_ew', ...
    'NPP_H','NPP_L','Ccrown_t','run_provenance','ERW');
bulk = S.PH_bulk_d(:);
nd = numel(bulk);
start = datenum(S.Datam(1,1),S.Datam(1,2),S.Datam(1,3)); %#ok<DATNM>
Date = start+(0:nd-1).'-2;
dv = datevec(Date);
complete = (1:nd).' >= 3;
rh = double(S.R_litter(:))+double(S.R_microbe(:))+double(S.R_ew(:));
cc = numel(S.NPP_L)/nd;
npp = sum(reshape(double(S.Ccrown_t),nd,cc).*(reshape(double(S.NPP_H),nd,cc)+reshape(double(S.NPP_L),nd,cc)),2);
nee = rh-npp;
n2o = double(S.Domain.N2Oflx(:));      % g N m-2 d-1
no3 = double(S.Domain.LEAK_NO3(:));    % g N m-2 d-1
bg = S.ERW{1}.bg;
n_t = numel(bg.pH);
per_day = round(n_t/nd);
step = datevec(start+(1:n_t)'/per_day-2);
released = cumsum(bg.dAlk_sil(:)+bg.dAlk_carb(:))/1e6;
n = numel(years);
a = struct('file',listing(1).name,'n_days',nd,'alpha',S.run_provenance.AppeloAlpha(1), ...
    'knmaxd',S.run_provenance.KNmaxd,'n2o_kg_ha',nan(1,n),'n2o_peak_g_ha_d',nan(1,n), ...
    'no3_kg_ha',nan(1,n),'nee_g_m2',nan(1,n),'rh_g_m2',nan(1,n),'ph_aug',nan(1,n), ...
    'ph_min',nan(1,n),'ph_max',nan(1,n),'released_aug',nan(1,n));
for y = 1:n
    in = complete & dv(:,1) == years(y);
    if ~any(in); continue; end
    a.n2o_kg_ha(y) = 10*sum(n2o(in),'omitnan');
    a.n2o_peak_g_ha_d(y) = 1e4*max(n2o(in));
    a.no3_kg_ha(y) = 10*sum(no3(in),'omitnan');
    a.nee_g_m2(y) = sum(nee(in),'omitnan');
    a.rh_g_m2(y) = sum(rh(in),'omitnan');
    a.ph_min(y) = min(bulk(in)); a.ph_max(y) = max(bulk(in));
    aug = in & dv(:,2) == 8;
    if any(aug); a.ph_aug(y) = mean(bulk(aug),'omitnan'); end
    in_step = step(:,1) == years(y) & step(:,2) == 8;
    if any(in_step); a.released_aug(y) = released(find(in_step,1,'last')); end
end
end
