function members = Uncertainty_Members
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The three parameter sets of the Kantola uncertainty ensemble, after the
%% T&C-SMEW field analysis (Zhang et al., GCB, SI section 4 and Table S2):
%% the original values, every parameter at the lower end of its range, and
%% every parameter at the upper end.
%%
%% Sixteen parameters move together: the twelve soil-biogeochemistry
%% literals of Table S2, the dissolution factor (0.8-1.2) and the exponent
%% on the ERW pore saturation (s^2 at the low end, sqrt(s) at the high end)
%% as the legacy analysis moved them, plus the two parameters this model
%% adds, KNmaxd (0.002-0.01 d^-1, the soil-model sensitivity bounds) and
%% the Appelo exchange slope alpha (1.3 and 6, the ends of the sweep grid,
%% around the campaign's 3.4, with the exchange constants slid by the
%% campaign's rule). "Lower" and "upper" name the ends of the parameter
%% ranges, not of any output: the sixteen moves do not all push a given
%% output the same way.
%%
%% The original is the calibrated campaign run itself (Base_1D, Base_2D and
%% the SSP projections, with their 5 t/ha twins), so the ensemble reuses
%% those results rather than running it again. Its base values are what
%% Biogeochemistry_Parameter.m ships: the SI's original column, except mrb
%% and mrf, which keep the legacy code's 0.010 and 0.005 under which the
%% spin-up was made.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

literals = { ...   % name, minimum, maximum  (Table S2)
    'Vm',   1.2,    528; ...
    'Vpc',  4.8,    792; ...
    'Km',   0.057,  0.95; ...
    'rem',  0.0075, 0.075; ...
    'rep',  0.0075, 0.075; ...
    'pepf', 0.0031, 0.031; ...
    'mrb',  0.0012, 0.03; ...
    'mrf',  0.0012, 0.03; ...
    'fd',   0.2,    0.8; ...
    'Vpl',  4.8,    792; ...
    'Vdb',  0.0024, 0.06; ...
    'Vdf',  0.0024, 0.06};
dissolution_factor = [0.8 1.2];
saturation_exponent = [2 0.5];          % low end s^2, high end sqrt(s)
knmaxd = [0.002 0.01];                  % d^-1
alpha_grid = linspace(1.3,3.4,5);       % formed as Sweep_Cases forms it
alpha = [alpha_grid(1) 6];              % the sweep grid's ends around the campaign's 3.4

members = struct('id',{1 2 3},'tag',{'original','lower','upper'}, ...
    'reused',{true false false},'BiogeoOverrides',[],'DissolutionFactor',[], ...
    'SaturationExponent',[],'KNmaxd',[],'AppeloAlpha',[]);
for L = 1:2
    members(L+1).BiogeoOverrides = cell2struct(literals(:,1+L),literals(:,1),1);
    members(L+1).DissolutionFactor = dissolution_factor(L);
    members(L+1).SaturationExponent = saturation_exponent(L);
    members(L+1).KNmaxd = knmaxd(L);
    members(L+1).AppeloAlpha = alpha(L);
end
end
