function members = Uncertainty_Members
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The three parameter sets of the Kantola uncertainty ensemble: the
%% original (calibrated) values and two moved sets, revised on 9 October
%% 2026 (UNCERTAINTY_BOUNDS_REVISION_2026-10-09.md) after the first
%% five-year run of the ensemble.
%%
%% Sixteen parameters move together: the twelve soil-biogeochemistry
%% literals the T&C-SMEW field analysis varied (Zhang et al., GCB, SI
%% section 4 and Table S2), the dissolution factor (0.8-1.2) and the
%% exponent on the ERW pore saturation (s^2 and sqrt(s)) as the legacy
%% analysis moved them, plus the two parameters this model adds, KNmaxd
%% and the Appelo exchange slope alpha.
%%
%% The ranges. The twelve literals move a factor of two around the shipped
%% values (fd, a fraction, by 0.1), clipped to Table S2 where the factor
%% would leave it (pepf's floor 0.0031, Vdb's ceiling 0.06). Table S2's
%% own ends are cross-site priors spanning two to three orders of
%% magnitude; run on pools spun up at the shipped values, they
%% mineralised about a decade of litter in the first year (the first
%% run's upper set), and a factor of two is the bracket within which the
%% spun-up state holds. KNmaxd runs 0.002-0.0064 d^-1 around the
%% calibrated 0.0048: the per-year chamber fits bracket it (2017 modelled
%% 2.4x the chambers, 2018 0.75x: 0.0048/2.4 and 0.0048/0.75), Hou et al.
%% (2024) measure 0.002-0.003 d^-1, and 0.0064 was run in the sweep grid.
%% Alpha runs 1.3-6 around the campaign's 3.4 (the fitted mesocosm value,
%% Appelo's and SCEPTER's default, and the value bracketing the field's
%% pH-versus-base-saturation slope), with the exchange constants slid by
%% the campaign's rule.
%%
%% The assignment. "lower" and "upper" name the ends of the basalt pH
%% response, not of the parameter ranges: each parameter's end sits in
%% the set it pushes that response towards. The lower set is the fast-
%% nitrogen, low-alpha, slow-weathering set (literals x2: more
%% nitrification acidity and nitrate export; KNmaxd 0.002: less
%% denitrification consuming acidity; alpha 1.3; dissolution factor 0.8;
%% s^2). The upper set is the slow-nitrogen, high-alpha, fast-weathering
%% set (literals x0.5; KNmaxd 0.0064; alpha 6; 1.2; sqrt(s)). Other
%% outputs need not follow: N2O pairs fast mineralisation with slow
%% denitrification and vice versa, so its band is narrow.
%%
%% The original is the calibrated campaign run itself (Base_1D, Base_2D and
%% the SSP projections, with their 5 t/ha twins), so the ensemble reuses
%% those results rather than running it again. Its base values are what
%% Biogeochemistry_Parameter.m ships: the SI's original column, except mrb
%% and mrf, which keep the legacy code's 0.010 and 0.005 under which the
%% spin-up was made. test_uncertainty_campaign pins every value here.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

literals = { ...   % name, shipped, lower set (x2), upper set (x0.5)
    'Vm',   100,    200,    50; ...
    'Vpc',  200,    400,    100; ...
    'Km',   0.25,   0.5,    0.125; ...
    'rem',  0.018,  0.036,  0.009; ...
    'rep',  0.018,  0.036,  0.009; ...
    'pepf', 0.006,  0.012,  0.0031; ...   % x0.5 = 0.003, below Table S2's floor
    'mrb',  0.010,  0.020,  0.005; ...
    'mrf',  0.005,  0.010,  0.0025; ...
    'fd',   0.4,    0.5,    0.3; ...      % a fraction: +/-0.1
    'Vpl',  23,     46,     11.5; ...
    'Vdb',  0.04,   0.06,   0.02; ...     % x2 = 0.08, above Table S2's ceiling
    'Vdf',  0.02,   0.04,   0.01};
dissolution_factor = [0.8 1.2];
saturation_exponent = [2 0.5];          % s^2 in the lower set, sqrt(s) in the upper
knmaxd = [0.002 0.0064];                % d^-1
alpha_grid = linspace(1.3,3.4,5);       % formed as Sweep_Cases forms it
alpha = [alpha_grid(1) 6];              % the sweep grid's ends around the campaign's 3.4

members = struct('id',{1 2 3},'tag',{'original','lower','upper'}, ...
    'reused',{true false false},'BiogeoOverrides',[],'DissolutionFactor',[], ...
    'SaturationExponent',[],'KNmaxd',[],'AppeloAlpha',[]);
for L = 1:2
    members(L+1).BiogeoOverrides = cell2struct(literals(:,2+L),literals(:,1),1);
    members(L+1).DissolutionFactor = dissolution_factor(L);
    members(L+1).SaturationExponent = saturation_exponent(L);
    members(L+1).KNmaxd = knmaxd(L);
    members(L+1).AppeloAlpha = alpha(L);
end
end
