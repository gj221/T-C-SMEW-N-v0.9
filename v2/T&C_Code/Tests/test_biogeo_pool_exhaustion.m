function test_biogeo_pool_exhaustion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A dissolved pool that the day's sinks empty must end at zero, not a
%% rounding error below it: the next day's uptake guard rejects any
%% negative pool. The all-minimum ensemble member stopped that way in
%% mid-September 2016, once its soybean took the soil's last mineral N
%% (-8.7e-27 g g^-1 soil after uptake).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

code_dir = fileparts(fileparts(mfilename('fullpath')));
addpath(code_dir);
unit_source = fileread(fullfile(code_dir,'BIOGEO_UNIT.m'));
assert(contains(unit_source,'= Biogeo_Uptake_Limit(') && ...
    contains(unit_source,'P = Biogeo_Pool_Update(Ptm1,dP,ZBIOG*rsd*1000);'), ...
    'BIOGEO_UNIT must limit uptake and update the pools through the helpers.');

%% Uptake below the headroom is untouched, bit for bit.
[scale,total] = Biogeo_Uptake_Limit([1 2],[0.1 0],[0.3 0.7]);
assert(isequal(scale,[1 1]) && isequal(total,[0.3 0.7]), ...
    'Uptake within the headroom must pass unchanged.');

%% Leaching keeps priority and uptake takes the rest.
[scale,total] = Biogeo_Uptake_Limit(1,0.4,1);
assert(abs(scale-0.6) < 1e-15 && total == 1-0.4, ...
    'Uptake must be scaled to the pool left after leaching.');

%% demand*(headroom/demand) rounds one ulp above the headroom here; the
%% total must stop at the headroom so the pool lands on zero.
headroom = 5.1023098237380445e-05;
demand = 0.0015238160855020927;
assert(demand*(headroom/demand) > headroom, ...
    'The fixture must reproduce the one-ulp overshoot.');
[~,total] = Biogeo_Uptake_Limit(headroom,0,demand);
soil = 3.25e5; %% 0.25 m x 1300 kg m^-3, in g soil m^-2
assert(total == headroom && headroom/soil-total/soil >= 0, ...
    'An exhausting uptake must leave the pool at zero, not below it.');

[~,total] = Biogeo_Uptake_Limit(1,0,NaN);
assert(isnan(total),'A nonfinite uptake must reach the pool guard.');
[scale,total] = Biogeo_Uptake_Limit(0,0,0.2);
assert(scale == 0 && total == 0,'An empty pool supplies no uptake.');

%% Daily update: a pool emptied by its sinks can round a few ulps below
%% zero; that is zero. A real overdraw is still an error. A fertilizer-size
%% gross term (20 g m^-2) must not shrink the rounding band.
Ptm1 = [1e-5 2 0 1e-6];
dP = [-(1e-5+eps(1e-5)) 0.5 0 (20-(20+1e-6))-4*eps(20)];
assert(all(Ptm1([1 4])+dP([1 4]) < 0),'The fixture must round below zero.');
P = Biogeo_Pool_Update(Ptm1,dP,soil);
assert(isequal(P,[0 2.5 0 0]), ...
    'A pool emptied to rounding must end at zero; the others add exactly.');
overdraw_failed = false;
try
    Biogeo_Pool_Update([1 1],[-1.1 0],soil);
catch ME
    overdraw_failed = strcmp(ME.identifier,'BIOGEO_UNIT:NegativePoolState');
end
assert(overdraw_failed,'A real overdraw of a pool must be rejected.');
P = Biogeo_Pool_Update([1 1],[NaN 0],soil);
assert(isnan(P(1)),'A nonfinite update must reach the finite-state guard.');

%% Hourly lateral exchange (2-D runs): a column drained of its nitrate
%% lands on zero, and a real overdraw is still an error.
P0 = zeros(1,79);
P0(32) = 1e-5;
dM = zeros(1,11);
dM(6) = -(1e-5+eps(1e-5)); %% NO3 is the sixth mobile species
assert(P0(32)+dM(6) < 0,'The lateral fixture must round below zero.');
P1 = Biogeo_Lateral_Pools(P0,dM);
assert(P1(32) == 0,'A column drained to rounding must hold no nitrate.');
dM(6) = -2e-5;
lateral_failed = false;
try
    Biogeo_Lateral_Pools(P0,dM);
catch ME
    lateral_failed = strcmp(ME.identifier,'Biogeo_Lateral_Pools:NegativePool');
end
assert(lateral_failed,'A real lateral overdraw must be rejected.');

fprintf('PASS: exhausted soil pools end at zero.\n');
end
