function [scale,total] = Biogeo_Uptake_Limit(pool,leak,total)
% Joint uptake and leaching limit on the dissolved mineral pools, g m^-2 and
% g m^-2 d^-1 over one day. Leaching keeps priority and uptake is scaled to
% the rest. The scaled total can round one ulp above that headroom, so it is
% set to it and an exhausted pool stops at zero instead of crossing it.

pool = pool(:).';
leak = leak(:).';
total = total(:).';
headroom = max(0,pool-leak);
scale = ones(size(total));
over = total > headroom;
scale(over) = headroom(over)./total(over);
total = total.*scale;
above = total > headroom;
total(above) = headroom(above);
end
