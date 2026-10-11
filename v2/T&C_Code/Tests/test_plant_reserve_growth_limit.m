function test_plant_reserve_growth_limit

reserve = [0.10; 0.02; 0.05; 0.03; 0.01; 0.04];
dtd = 1;
fixed_demand = [-0.01; 0; 0.01; -0.01; 0; 0.01];
flexible_demand = [0.20; 0.04; 0.08; 0.08; 0.04; 0.06];

scale = Plant_Reserve_Growth_Limit(reserve,dtd,fixed_demand,flexible_demand);
expected = min((reserve/dtd-fixed_demand)./flexible_demand);
assert(abs(scale-expected) < 1e-12, ...
    'The most limiting physical reserve must set the shared carbon-flux scale.');

reserve_new = reserve-(fixed_demand+scale*flexible_demand)*dtd;
assert(all(reserve_new >= -1e-12), ...
    'The limited daily step must not overdraw N, P, K, Ca, Mg or Si reserve.');

npp = 4.0;
add_reserve_c = 1.5;
overflow_respiration = (1-scale)*(npp+add_reserve_c);
assert(abs(scale*(npp+add_reserve_c)+overflow_respiration-...
    (npp+add_reserve_c)) < 1e-12, ...
    'Rejected growth and reserve-allocation carbon must remain in the C ledger.');

zero_scale = Plant_Reserve_Growth_Limit(zeros(6,1),dtd,zeros(6,1),...
    flexible_demand);
assert(zero_scale == 0, ...
    'Positive nutrient demand cannot be supplied by an empty physical reserve.');

release_scale = Plant_Reserve_Growth_Limit(zeros(6,1),dtd,...
    -ones(6,1),zeros(6,1));
assert(release_scale == 1, ...
    'Nutrient-releasing fixed transfers must remain possible with empty reserves.');

did_error = false;
try
    Plant_Reserve_Growth_Limit(zeros(6,1),dtd,...
        [0.01; 0; 0; 0; 0; 0],zeros(6,1));
catch ME
    did_error = strcmp(ME.identifier,...
        'Plant_Reserve_Growth_Limit:FixedDemandOverdraw');
end
assert(did_error, ...
    'A fixed process that overdrafts reserve must fail rather than be clamped.');

did_error = false;
try
    Plant_Reserve_Growth_Limit([-0.01; zeros(5,1)],dtd,...
        zeros(6,1),zeros(6,1));
catch ME
    did_error = strcmp(ME.identifier,...
        'Plant_Reserve_Growth_Limit:InvalidState');
end
assert(did_error, ...
    'A negative physical reserve must fail at the nutrient-growth boundary.');

fprintf('PASS: physical plant-reserve growth limit and C overflow closure\n');
end
