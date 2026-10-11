function scale = Plant_Reserve_Growth_Limit(reserve,dtd,fixed_demand,flexible_demand)

reserve = reserve(:);
fixed_demand = fixed_demand(:);
flexible_demand = flexible_demand(:);

if numel(reserve) ~= 6 || numel(fixed_demand) ~= 6 || ...
        numel(flexible_demand) ~= 6 || ~isfinite(dtd) || dtd <= 0 || ...
        any(~isfinite(reserve)) || any(reserve < 0) || ...
        any(~isfinite(fixed_demand)) || any(~isfinite(flexible_demand))
    error('Plant_Reserve_Growth_Limit:InvalidState',...
        'Plant reserves and nutrient demands must be finite and physical.');
end

available = reserve/dtd;
tol = 64*eps(max(1,max(abs([available; fixed_demand; flexible_demand]))));
lower = 0;
upper = 1;

for j = 1:6
    if flexible_demand(j) > tol
        upper = min(upper,(available(j)-fixed_demand(j))/flexible_demand(j));
    elseif flexible_demand(j) < -tol
        lower = max(lower,(available(j)-fixed_demand(j))/flexible_demand(j));
    elseif fixed_demand(j) > available(j)+tol
        error('Plant_Reserve_Growth_Limit:FixedDemandOverdraw',...
            'A fixed plant transfer exceeds the available physical reserve.');
    end
end

lower = max(0,lower);
upper = min(1,upper);
if lower > upper+tol
    error('Plant_Reserve_Growth_Limit:InfeasibleDemand',...
        'No shared carbon-flux scale satisfies every physical reserve.');
end

scale = max(lower,min(1,max(0,upper)));
end
