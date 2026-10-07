function test_erw_al_leaching_scale
%%% Aluminium leaching is formed on the solver's Al_CEC_scale basis, the
%%% same construction as the sodium and residual-anion leaching terms.
%%%
%%% The legacy TC-SMEW line (Other models/T&C-SMEW, ERW_function.m lines
%%% 263 and 271) carried an extra factor 1e-3 on the mobile-inventory term
%%% of Leak_Al. Al, AlOH4 and Al_tot are all on the conv_Al (Al_CEC_scale)
%%% basis: the Al_w0 seed divides (Al_tot - exchangeable Al*Al_CEC_scale)
%%% by the pore-water volume, AlOH4 closes Al_w by subtraction, and Al_tot
%%% accumulates EW*Al_CEC_scale. So (Al+AlOH4)*water_l*f_Lk_biogeo is
%%% already on Al_tot's basis and the factor put that one term a thousand
%%% times below the other two in the min(). This test pins the corrected
%%% source and replays the construction.

test_dir = fileparts(mfilename('fullpath'));
code_dir = fullfile(test_dir,'..');

erw_source = fileread(fullfile(code_dir,'ERW_function.m'));
assert(~contains(erw_source,'f_Lk_biogeo*1e-3'), ...
    'ERW_function.m must not scale the aluminium leaching term by 1e-3.');
with_drainage = ['Leak_Al = max(0,min([(Al(j_only)+AlOH4(j_only))*n*Zrw*s*1000*f_Lk_biogeo, ' ...
    'Al_tot(j_only), Lw(j_only)*1000*(Al(j_only)+AlOH4(j_only))]));'];
without_drainage = ['Leak_Al = max(0,min([(Al(j_only)+AlOH4(j_only))*n*Zrw*s*1000*f_Lk_biogeo, ' ...
    'Al_tot(j_only)]));'];
assert(contains(erw_source,with_drainage) && contains(erw_source,without_drainage), ...
    'Both Leak_Al branches must form the mobile inventory on the Al_CEC_scale basis.');
assert(contains(erw_source,'Al_w0 = max((Al_tot-(f_Al(j_idx)/3)*CEC_seed*Al_CEC_scale)/(n*Zrw*s*1000), 0);'), ...
    'The Al_w0 seed line that fixes the scale must be unchanged.');
assert(contains(erw_source,'Al_tot = Al_tot(j_only)+(Al_weathering-Leak_Al)*dtw;') && ...
    contains(erw_source,'Al_weathering = sum(min_st(:,5)''.*EW(j_only,:))*Al_CEC_scale;'), ...
    'The Al_tot balance and its Al_CEC_scale weathering source must be unchanged.');

% Replay on a synthetic state on the solver's scales: Al_CEC_scale = 1e3,
% Al_tot in nmol/m2, pore water n*Zrw*s*1000 in L/m2 (line 653's own
% denominator), f_Al = 0 so the seed is the dissolved inventory alone.
Al_CEC_scale = 1e3;
Al_tot = 5.0e5;
f_Al = 0; CEC_seed = 5.07e7;
n = 0.45; Zrw = 0.3; s = 0.7;
water_l = n*Zrw*s*1000;
Al_w0 = max((Al_tot-(f_Al/3)*CEC_seed*Al_CEC_scale)/water_l, 0);
assert(Al_w0 > 0,'The replayed seed must be positive.');
Al_free = 0.3*Al_w0; AlOH4 = 0.7*Al_w0;      % the two mobile species, on Al_w0's units
f_Lk_biogeo = 0.02; Lw = 0.0007;             % 1/d, m/d

leak_fixed = max(0,min([(Al_free+AlOH4)*water_l*f_Lk_biogeo, Al_tot, Lw*1000*(Al_free+AlOH4)]));
leak_legacy = max(0,min([(Al_free+AlOH4)*water_l*f_Lk_biogeo*1e-3, Al_tot, Lw*1000*(Al_free+AlOH4)]));
term_fixed = (Al_free+AlOH4)*water_l*f_Lk_biogeo;
term_legacy = term_fixed*1e-3;

% With every dissolved Al mobile and f_Al = 0, the mobile-inventory term is
% exactly Al_tot*f_Lk_biogeo: the Na/An form pool*fraction on Al_tot's basis.
assert(abs(term_fixed-Al_tot*f_Lk_biogeo) <= 1e-12*Al_tot*f_Lk_biogeo, ...
    'The corrected term must equal Al_tot*f_Lk_biogeo when all dissolved Al is mobile.');
assert(abs(term_fixed/term_legacy-1000) < 1e-9, ...
    'The legacy factor scaled the mobile-inventory term by exactly 1e-3.');
assert(leak_fixed == min(term_fixed,min(Al_tot,Lw*1000*(Al_free+AlOH4))), ...
    'Leak_Al must be the minimum of the three terms on one basis.');
assert(leak_legacy > 0 && leak_fixed/leak_legacy > 10, ...
    'At a slow-drainage state the correction must raise Leak_Al by at least an order of magnitude.');
assert(leak_fixed == Lw*1000*(Al_free+AlOH4), ...
    'At this state the drainage-flux cap, not the flawed term, must bind after the correction.');

% Dimensional closure of the Al_tot balance: the leaching term subtracts
% from Al_tot in the units Al_weathering adds (nmol/m2/d = umol/m2/d *
% Al_CEC_scale), so a mobile inventory drained at f_Lk_biogeo per day must
% drain the same fraction of Al_tot per day when it is the whole inventory.
EW_Al = 12.5;                                 % umol Al / m2 / d from weathering
Al_weathering = EW_Al*Al_CEC_scale;
dtw = 1/144;
Al_tot_next = Al_tot+(Al_weathering-term_fixed)*dtw;
assert(abs(Al_tot_next-(Al_tot*(1-f_Lk_biogeo*dtw)+Al_weathering*dtw)) <= 1e-9*Al_tot, ...
    'The Al_tot balance must drain f_Lk_biogeo of the mobile inventory per day on the weathering basis.');

fprintf('PASS: aluminium leaching is formed on the solver Al scale, matching Na and An.\n');
end
