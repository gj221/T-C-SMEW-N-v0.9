function A = Kantola_N_Magnitude_Audit(S,scenario,nd)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Kantola_N_Magnitude_Audit
%%
%% Retained states are g N m^-2 and retained daily fluxes are
%% g N m^-2 d^-1.  One g m^-2 equals 10 kg ha^-1.  On drainage days,
%% (g N m^-2 d^-1)/(mm d^-1) equals g N L^-1 because 1 mm m^-2 is 1 L.
%% Min_N is net mineralization, (NH4_in+NO3_in) in BIOGEOCHEMISTRY_DYNAMIC3,
%% and is negative on days of net microbial immobilization; it is an
%% internal transfer, not part of the boundary ledger, and is reported
%% signed.
%% Row 1 of every retained state is written before the first daily block:
%% Ccrown_t starts as ones() and its row 1 is never overwritten, while
%% Nreserve_L(1,:) holds the legacy seed reserve of 1000 g N per crown m^2
%% and crop (MOD_PARAM_Kantola_Maize_Soy). Row-1 plant-side products are
%% therefore per-crown constants under a placeholder crown of one, not
%% ground-area pools. The ledger takes soil pools from row 1 and plant-side
%% pools from row 2, whose retained crown is the campaign's actual crown
%% ([0 0] on 1 January, so the carried reserve is masked); the row-1
%% placeholder product is reported as row1_placeholder_plant_N_gNm2.
%% The daily biogeochemistry block runs at the hour with Datam(i,4) == 1
%% (01:00) and receives sum(Lk(pdind)) with pdind = max(1,i-24):i-1
%% (MAIN_FRAME_BG_ERW_2D), i.e. the 24 hours from 01:00 of the previous
%% day through 00:00 of the current day. The dry-day check and the mobile
%% NO3 concentration rebuild exactly that window from the retained hourly
%% dates; on the first retained day the window holds only the hours before
%% the first block.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

required = {'P','TNIT_H','TNIT_L','Nreserve_H','Nreserve_L', ...
    'FertN_d','DepN_d','BfixN','N2Oflx','N2flx','VOL', ...
    'LEAK_NO3','LEAK_NH4','LEAK_DON','Lk','AgrHarNut','SeedInput', ...
    'Nuptake_H','Nuptake_L','DENITflx','Min_N','Ccrown_t', ...
    'ISOIL_H','ISOIL_L','TexN_H','TexN_L','Date'};
for k = 1:numel(required)
    assert(isfield(S,required{k}), ...
        'Kantola_N_Magnitude_Audit:MissingField', ...
        'Missing retained field %s.',required{k});
end

if nargin < 3 || isempty(nd)
    assert(isfield(S,'i'),'Kantola_N_Magnitude_Audit:MissingDay', ...
        'A completed day or retained iteration is required.');
    nd = min(size(S.P,1)-1,floor(double(S.i)/24));
end
assert(isscalar(nd) && isfinite(nd) && nd == floor(nd) && nd >= 1, ...
    'Kantola_N_Magnitude_Audit:BadDay','Completed day must be a positive integer.');
assert(size(S.P,1) >= nd+1 && size(S.P,2) >= 34, ...
    'Kantola_N_Magnitude_Audit:ShortState', ...
    'The retained pool state does not contain day %d.',nd);

scenario = lower(char(scenario));
assert(ismember(scenario,{'control','basalt','calcite'}), ...
    'Kantola_N_Magnitude_Audit:BadScenario','Unknown scenario %s.',scenario);

P = double(S.P(1:nd+1,:));
assert(all(isfinite(P(:))) && min(P(:,23:34),[],'all') >= -1e-12, ...
    'Kantola_N_Magnitude_Audit:InvalidPool', ...
    'Retained nitrogen pools are nonfinite or negative.');

A = struct();
A.scenario = scenario;
A.day = nd;
A.standing_NH4_gNm2 = P(end,31);
A.standing_NO3_gNm2 = P(end,32);
A.standing_NH4_kgNha = 10*A.standing_NH4_gNm2;
A.standing_NO3_kgNha = 10*A.standing_NO3_gNm2;
A.cumulative_NO3_leak_gNm2 = sum_daily(S.LEAK_NO3,nd,'LEAK_NO3');
A.cumulative_NH4_leak_gNm2 = sum_daily(S.LEAK_NH4,nd,'LEAK_NH4');
A.cumulative_NO3_leak_kgNha = 10*A.cumulative_NO3_leak_gNm2;
A.cumulative_NH4_leak_kgNha = 10*A.cumulative_NH4_leak_gNm2;
A.plant_uptake_gNm2 = sum_weighted_daily( ...
    S.Nuptake_H,S.Ccrown_t,nd,'Nuptake_H') + ...
    sum_weighted_daily(S.Nuptake_L,S.Ccrown_t,nd,'Nuptake_L');
A.denitrification_gNm2 = sum_daily(S.DENITflx,nd,'DENITflx');
A.fertilizer_gNm2 = sum_daily(S.FertN_d,nd,'FertN_d');
net_mineralization = signed_daily_state(S.Min_N,nd,'Min_N');
A.mineralization_gNm2 = sum(net_mineralization);
A.net_immobilization_gNm2 = -sum(net_mineralization(net_mineralization < 0));
A.net_immobilization_days = nnz(net_mineralization < -1e-12);
A.seed_input_gNm2 = sum_daily(S.SeedInput(:,2),nd,'SeedInput(:,2)');

leak_no3 = daily_state(S.LEAK_NO3,nd,'LEAK_NO3');
drainage = hourly_to_daily(S.Lk,nd,'Lk'); %#ok<NASGU> retained daily total
driving_drainage = block_window_sum(S.Lk,S.Date,nd,'Lk');
wet = driving_drainage > 0;
assert(all(leak_no3(~wet) <= 1e-12), ...
    'Kantola_N_Magnitude_Audit:DryLeakage', ...
    'NO3 leakage is positive while the drainage window of its daily block was dry.');
mobile = NaN(nd,1);
mobile(wet) = 1000*leak_no3(wet)./driving_drainage(wet);
if any(wet)
    A.maximum_NO3_mobile_mgNL = max(mobile(wet));
    A.latest_NO3_mobile_mgNL = mobile(find(wet,1,'last'));
else
    A.maximum_NO3_mobile_mgNL = NaN;
    A.latest_NO3_mobile_mgNL = NaN;
end

pending_litter_n = weighted_pending_litter_n( ...
    S.ISOIL_H,S.Ccrown_t,nd+1,'ISOIL_H') + ...
    weighted_pending_litter_n(S.ISOIL_L,S.Ccrown_t,nd+1,'ISOIL_L');
total_transfer_n = weighted_row_sum( ...
    S.TexN_H,S.Ccrown_t,nd+1,'TexN_H') + ...
    weighted_row_sum(S.TexN_L,S.Ccrown_t,nd+1,'TexN_L');
organ_harvest_n = total_transfer_n-pending_litter_n;
assert(all(isfinite(organ_harvest_n)) && min(organ_harvest_n) >= -1e-10, ...
    'Kantola_N_Magnitude_Audit:OrganExport', ...
    'Derived organ harvest N export is nonfinite or negative.');
plant_n = weighted_plant_n(S.TNIT_H,S.Nreserve_H,S.Ccrown_t,nd+1,'H') + ...
    weighted_plant_n(S.TNIT_L,S.Nreserve_L,S.Ccrown_t,nd+1,'L') + ...
    pending_litter_n;
storage = sum(P(:,23:34),2) + plant_n;
storage(1) = sum(P(1,23:34)) + plant_n(2);
A.row1_placeholder_plant_N_gNm2 = plant_n(1)-plant_n(2);
% The row-2 rule is exact only when no plant-side transfer happened on day
% 1, which holds when no crop is present, or when no legacy reserve is
% carried; refuse any other configuration rather than mis-close silently.
assert(all(double(S.Ccrown_t(2,:)) == 0) || ...
    max(double(S.Nreserve_L(2,:)),[],'all') < 1, ...
    'Kantola_N_Magnitude_Audit:LegacyReserveLive', ...
    'A crop is present on day 1 while the legacy seed reserve is still carried.');
inputs = A.fertilizer_gNm2 + A.seed_input_gNm2 + ...
    sum_daily(S.DepN_d,nd,'DepN_d') + sum_daily(S.BfixN,nd,'BfixN');
outputs = sum_daily(S.N2Oflx,nd,'N2Oflx') + ...
    sum_daily(S.N2flx,nd,'N2flx') + sum_daily(S.VOL,nd,'VOL') + ...
    A.cumulative_NO3_leak_gNm2 + A.cumulative_NH4_leak_gNm2 + ...
    sum_daily(S.LEAK_DON,nd,'LEAK_DON') + ...
    sum_daily(S.AgrHarNut(:,1),nd,'AgrHarNut(:,1)') + ...
    sum(organ_harvest_n(2:end));
A.N_storage_change_gNm2 = storage(end)-storage(1);
A.pending_litter_N_gNm2 = pending_litter_n(end);
A.organ_harvest_export_gNm2 = sum(organ_harvest_n(2:end));
A.N_closure_residual_gNm2 = inputs-outputs-A.N_storage_change_gNm2;
assert(abs(A.N_closure_residual_gNm2) <= 1e-8, ...
    'Kantola_N_Magnitude_Audit:NClosure', ...
    'Nitrogen boundary ledger residual is %.12g g N m^-2.', ...
    A.N_closure_residual_gNm2);

if strcmp(scenario,'control')
    A.observed_2017_NO3_leak_kgNha = 67.06;
elseif strcmp(scenario,'basalt')
    A.observed_2017_NO3_leak_kgNha = 77.49;
else
    A.observed_2017_NO3_leak_kgNha = NaN;
end
is_2017 = calendar_year_mask(S.Date,nd,2017);
if nnz(is_2017) == 365
    A.modeled_2017_NO3_leak_kgNha = 10*sum(leak_no3(is_2017));
    A.NO3_leak_2017_ratio = A.modeled_2017_NO3_leak_kgNha ./ ...
        A.observed_2017_NO3_leak_kgNha;
else
    A.modeled_2017_NO3_leak_kgNha = NaN;
    A.NO3_leak_2017_ratio = NaN;
end
end

function value = sum_daily(x,nd,name)
value = sum(daily_state(x,nd,name));
end

function value = daily_state(x,nd,name)
assert(size(x,1) >= nd+1, ...
    'Kantola_N_Magnitude_Audit:ShortFlux','%s lacks the state row for day %d.',name,nd);
x = double(x(2:nd+1,:));
assert(all(isfinite(x(:))) && min(x(:)) >= -1e-12, ...
    'Kantola_N_Magnitude_Audit:InvalidFlux','%s is nonfinite or negative.',name);
value = sum(x,2);
end

function value = signed_daily_state(x,nd,name)
% Net fluxes carry their sign; only finiteness is a defect.
assert(size(x,1) >= nd+1, ...
    'Kantola_N_Magnitude_Audit:ShortFlux','%s lacks the state row for day %d.',name,nd);
x = double(x(2:nd+1,:));
assert(all(isfinite(x(:))), ...
    'Kantola_N_Magnitude_Audit:InvalidFlux','%s is nonfinite.',name);
value = sum(x,2);
end

function value = block_window_sum(x,dates,nd,name)
% Sum x over the model's daily-block window for each retained day: the
% block for retained day d runs at the d-th hour with hour-of-day 1 and
% reads the 24 hourly values before it (MAIN_FRAME_BG_ERW_2D pdind).
assert(size(x,1) >= 24*nd && numel(dates) >= 24*nd, ...
    'Kantola_N_Magnitude_Audit:ShortFlux','%s or Date ends before hour %d.',name,24*nd);
x = double(x(1:24*nd,:));
assert(all(isfinite(x(:))) && min(x(:)) >= -1e-12, ...
    'Kantola_N_Magnitude_Audit:InvalidFlux','%s is nonfinite or negative.',name);
hour = round(mod(double(dates(1:24*nd)),1)*24);
block = find(hour(:) == 1);
assert(numel(block) == nd, ...
    'Kantola_N_Magnitude_Audit:DailyBlockCount', ...
    'The retained dates hold %d daily-block hours for %d days.',numel(block),nd);
hourly = sum(x,2);
value = zeros(nd,1);
for d = 1:nd
    value(d) = sum(hourly(max(1,block(d)-24):block(d)-1));
end
end

function value = hourly_to_daily(x,nd,name)
assert(size(x,1) >= 24*nd, ...
    'Kantola_N_Magnitude_Audit:ShortFlux','%s ends before hour %d.',name,24*nd);
x = double(x(1:24*nd,:));
assert(all(isfinite(x(:))) && min(x(:)) >= -1e-12, ...
    'Kantola_N_Magnitude_Audit:InvalidFlux','%s is nonfinite or negative.',name);
value = sum(reshape(sum(x,2),24,nd),1)';
end

function value = sum_weighted_daily(x,w,nd,name)
value = sum(weighted_daily_state(x,w,nd,name));
end

function value = weighted_daily_state(x,w,nd,name)
assert(size(x,1) >= nd+1 && size(w,1) >= nd+1, ...
    'Kantola_N_Magnitude_Audit:ShortFlux', ...
    '%s or Ccrown_t lacks the state row for day %d.',name,nd);
x = double(x(2:nd+1,:));
w = double(w(2:nd+1,:));
assert(isequal(size(x),size(w)), ...
    'Kantola_N_Magnitude_Audit:CrownShape', ...
    '%s and Ccrown_t have different retained dimensions.',name);
assert(all(isfinite(x(:))) && min(x(:)) >= -1e-12, ...
    'Kantola_N_Magnitude_Audit:InvalidFlux','%s is nonfinite or negative.',name);
assert(all(isfinite(w(:))) && min(w(:)) >= -1e-12 && max(w(:)) <= 1+1e-12, ...
    'Kantola_N_Magnitude_Audit:InvalidCrown', ...
    'Ccrown_t is nonfinite or outside [0,1].');
value = sum(x.*w,2);
end

function value = weighted_row_sum(x,w,n,name)
assert(size(x,1) >= n && size(w,1) >= n, ...
    'Kantola_N_Magnitude_Audit:ShortState', ...
    '%s or Ccrown_t ends before retained row %d.',name,n);
x = double(x(1:n,:));
w = double(w(1:n,:));
assert(isequal(size(x),size(w)), ...
    'Kantola_N_Magnitude_Audit:CrownShape', ...
    '%s and Ccrown_t have different retained dimensions.',name);
assert(all(isfinite(x(:))) && min(x(:)) >= -1e-12, ...
    'Kantola_N_Magnitude_Audit:InvalidState','%s is nonfinite or negative.',name);
assert(all(isfinite(w(:))) && min(w(:)) >= -1e-12 && max(w(:)) <= 1+1e-12, ...
    'Kantola_N_Magnitude_Audit:InvalidCrown', ...
    'Ccrown_t is nonfinite or outside [0,1].');
value = sum(x.*w,2);
end

function value = weighted_plant_n(tissue,reserve,w,n,name)
% A reserve below zero is tissue diluted below its nominal N (legacy T&C),
% so the sign check falls on the plant total, tissue plus reserve.
assert(isequal(size(tissue),size(reserve)), ...
    'Kantola_N_Magnitude_Audit:CrownShape', ...
    'TNIT_%s and Nreserve_%s have different retained dimensions.',name,name);
value = weighted_row_sum(double(tissue)+double(reserve),w,n, ...
    ['TNIT_' name '+Nreserve_' name]);
end

function value = weighted_pending_litter_n(x,w,n,name)
assert(size(x,1) >= n && size(w,1) >= n && ...
        size(x,2) == size(w,2) && size(x,3) >= 12, ...
    'Kantola_N_Magnitude_Audit:PendingShape', ...
    '%s or Ccrown_t has incompatible retained dimensions.',name);
x = double(x(1:n,:,10:12));
w = double(w(1:n,:));
assert(all(isfinite(x(:))) && min(x(:)) >= -1e-12, ...
    'Kantola_N_Magnitude_Audit:InvalidPending', ...
    '%s pending litter N is nonfinite or negative.',name);
assert(all(isfinite(w(:))) && min(w(:)) >= -1e-12 && max(w(:)) <= 1+1e-12, ...
    'Kantola_N_Magnitude_Audit:InvalidCrown', ...
    'Ccrown_t is nonfinite or outside [0,1].');
value = sum(sum(x,3).*w,2);
end

function mask = calendar_year_mask(x,nd,target_year)
assert(numel(x) >= 24*nd, ...
    'Kantola_N_Magnitude_Audit:ShortDate', ...
    'Date ends before the last completed hour of day %d.',nd);
x = double(x(:));
x = x(1:24*nd);
assert(all(isfinite(x)) && all(diff(x) > 0), ...
    'Kantola_N_Magnitude_Audit:InvalidDate', ...
    'Retained hourly dates must be finite and strictly increasing.');
day_date = floor(x(1:24:end));
assert(all(diff(day_date) == 1), ...
    'Kantola_N_Magnitude_Audit:InvalidDate', ...
    'Retained hourly dates do not map to consecutive model days.');
date_parts = datevec(day_date);
mask = date_parts(:,1) == target_year;
end
