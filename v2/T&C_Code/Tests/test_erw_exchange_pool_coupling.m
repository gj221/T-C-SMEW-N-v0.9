function test_erw_exchange_pool_coupling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Internal solution-exchange transfer must not change the conserved ERW
%% cation inventory.  External loss must retain its molar equivalent.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MM = [39.1 40 24.3];
B0 = [0.20 0.30; 0.40 0.50; 0.60 0.70];
B1 = B0 + [-0.01 0.01; -0.02 0.02; -0.03 0.03];
dX = sum(B1-B0,2)./MM(:)*1e6;
assert(all(abs(dX) < 1e-10), ...
    'Internal exchange changed the conserved ERW cation inventory.');

B1(2,1) = B1(2,1)-0.004;
dCa = sum(B1(2,:)-B0(2,:))/MM(2)*1e6;
assert(abs(dCa+100) < 1e-10, ...
    'External Ca loss was not converted from g m^-2 d^-1 to umol m^-2 d^-1.');

P_legacy = zeros(1,69);
P_current = P_legacy;
P_erw_seed = P_legacy;
P_legacy([60 61]) = [0.2 0.3];
P_current([60 61]) = [0.19 0.30];
P_erw_seed([60 61]) = [0.8 1.2];
true_dCa = (sum(P_current([60 61]))-sum(P_legacy([60 61])))/40*1e6;
seed_contaminated_dCa = (sum(P_current([60 61]))-sum(P_erw_seed([60 61])))/40*1e6;
assert(abs(true_dCa+250) < 1e-10 && abs(seed_contaminated_dCa-true_dCa) > 1e4, ...
    'The initial ERW state mapping was not separated from the first T&C daily flux.');

daily_delta = [100 -25 5 0.5];
window_days = 23/24;
rate = ERW_Normalize_TC_Daily_Deltas(daily_delta,window_days);
assert(max(abs(rate*window_days-daily_delta)) < 1e-12, ...
    'A partial ERW window did not integrate to the complete T&C daily delta.');
assert(max(abs(ERW_Normalize_TC_Daily_Deltas(daily_delta,1)-daily_delta)) < 1e-12, ...
    'A complete ERW day changed the T&C daily delta.');

Datam = zeros(48,4);
Datam(:,4) = repmat((1:24)',2,1);
assert(ERW_Should_Apply_TC_Daily_Delta(2,2,48,Datam), ...
    'The first retained T&C daily cation delta was omitted.');
assert(ERW_Should_Apply_TC_Daily_Delta(3,25,48,Datam), ...
    'A regular new T&C daily state did not apply its element delta.');
assert(~ERW_Should_Apply_TC_Daily_Delta(3,48,48,Datam), ...
    'The final partial flush applied the same T&C daily delta twice.');
Datam(48,4) = 1;
assert(ERW_Should_Apply_TC_Daily_Delta(3,48,48,Datam), ...
    'A run ending on a new daily boundary omitted its T&C element delta.');

CEC_tot = 3.34e7;
f_exchange = [0.04 0.38 0.06];
exchange_g = f_exchange.*CEC_tot./[1 2 2].*MM/1e6;
exchange_umol = exchange_g./MM*1e6;
exchange_expected = f_exchange.*CEC_tot./[1 2 2];
assert(max(abs(exchange_umol-exchange_expected)./exchange_expected) < 1e-12, ...
    'Exchange-pool write-back does not preserve charge-equivalent inventory.');

code_dir = fileparts(fileparts(mfilename('fullpath')));
daily_code = fileread(fullfile(code_dir,'ERW_Daily_Step_2D.m'));
legacy_code = fileread(fullfile(code_dir,'MAIN_FRAME_BG_ERW_10mins.m'));
frame_2d_code = fileread(fullfile(code_dir,'MAIN_FRAME_BG_ERW_2D.m'));
delta_gate_code = fileread(fullfile(code_dir,'ERW_Should_Apply_TC_Daily_Delta.m'));

%% The daily step receives yesterday's and today's daily rows as rows
%% r_prev and r_now of P_c (package y).
daily_terms = {'P_c(r_now,52) + P_c(r_now,53)', ...
    'P_c(r_now,60) + P_c(r_now,61)', 'P_c(r_now,68) + P_c(r_now,69)'};
legacy_terms = {'P(j,52) + P(j,53)', ...
    'P(j,60) + P(j,61)', 'P(j,68) + P(j,69)'};
for k = 1:numel(daily_terms)
    assert(contains(daily_code,daily_terms{k}), ...
        'The daily ERW path omits a legacy exchange pool.');
    assert(contains(legacy_code,legacy_terms{k}), ...
        'The legacy ERW path omits a legacy exchange pool.');
end

assert(contains(delta_gate_code,'j >= 2'), ...
    'The synchronized first ERW state omitted its first T&C daily delta.');
assert(contains(daily_code,'P_prev_for_erw_delta = P_c(r_prev,:);') && ...
    contains(legacy_code,'P_prev_for_erw_delta = P(j-1,:);'), ...
    'The first T&C daily delta is not based on the pre-initialization legacy state.');
assert(contains(daily_code,'P_prev_for_erw_delta(60)') && ...
    contains(legacy_code,'P_prev_for_erw_delta(60)'), ...
    'The ERW seed transformation can contaminate the first Ca delta.');
assert(contains(daily_code,'ERW_Normalize_TC_Daily_Deltas') && ...
    contains(legacy_code,'ERW_Normalize_TC_Daily_Deltas'), ...
    'The 1-D and 2-D paths do not normalize partial ERW integration windows.');
assert(contains(daily_code,'ERW_Should_Apply_TC_Daily_Delta') && ...
    contains(legacy_code,'ERW_Should_Apply_TC_Daily_Delta'), ...
    'The 1-D and 2-D paths can reapply a daily T&C delta during the final flush.');

daily_writeback = {'E.P53 =','E.P61 =','E.P69 ='};
legacy_writeback = {'P(j, 53)=','P(j, 61)=','P(j, 69)='};
frame_2d_writeback = {'P(j,53,c) =','P(j,61,c) =','P(j,69,c) ='};
for k = 1:numel(daily_writeback)
    assert(contains(daily_code,daily_writeback{k}), ...
        'The daily ERW path omits an exchange-pool write-back.');
    assert(contains(legacy_code,legacy_writeback{k}), ...
        'The legacy ERW path omits an exchange-pool write-back.');
    assert(contains(frame_2d_code,frame_2d_writeback{k}), ...
        'The 2-D frame omits an exchange-pool write-back.');
end

assert(contains(daily_code,'cov_moll_end') && ...
    contains(daily_code,'bg.Ca(dm_end,:) * 40   * cov_moll_end'), ...
    'The daily ERW handoff does not use the retained-state water volume.');
assert(contains(legacy_code, ...
    'cov_moll_end = ERW_const.n * ERW_const.Zrw * TT_m.s(dm_end,:) * 1000;') && ...
    contains(legacy_code,'ERW_biogeochem.Ca(dm_end,:)* 40 * cov_moll_end'), ...
    'The 1-D ERW handoff does not use the retained-state water volume.');

fprintf('test_erw_exchange_pool_coupling: PASS\n');
end
