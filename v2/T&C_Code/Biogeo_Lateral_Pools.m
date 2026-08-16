function [P_out, mobile, names, pool_idx, fraction_fields] = Biogeo_Lateral_Pools(P_in, dM)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pack and update the pools treated as dissolved by Biogeo_Leakage.
%% P_in and dM have one row per cell. Pool units remain gX/m2.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

names = {'DOC_bacteria','DOC_fungi','DON','DOP','NH4','NO3','P','K','Ca','Mg','Si'};
pool_idx = {12,13,33,47,31,32,43,52,60,68,76};
fraction_fields = {'aDOC','aDOC','aDON','aDOP','aNH4','aNO3', ...
                   'aP','aK','aCa','aMg','aSi'};

if ~(isnumeric(P_in) && ismatrix(P_in) && size(P_in,2) >= 76)
    error('Biogeo_Lateral_Pools:BadPools', ...
        'P_in must have one row per cell and at least 76 pool columns.');
end

source_idx = unique([pool_idx{:}]);
source = P_in(:,source_idx);
%% A pool driven to exactly zero by the donor limiters can land a
%% rounding step below it. Tolerate and clamp rather than aborting a
%% multi-hour run; a meaningful negative still errors, with its size.
pool_tol = 1e-12*max(1,max(abs(source(:))));
if any(~isfinite(source(:))) || any(source(:) < -pool_tol)
    error('Biogeo_Lateral_Pools:BadPools', ...
        ['Laterally mobile pools must be finite and nonnegative; most ' ...
         'negative is %.6g against a %.3g tolerance.'], ...
        min(source(:)),pool_tol);
end
P_in(:,source_idx) = max(0,source);
P_out = P_in;
mobile = zeros(size(P_in,1),numel(pool_idx));
for s = 1:numel(pool_idx)
    mobile(:,s) = sum(P_in(:,pool_idx{s}),2);
end

if nargin < 2 || isempty(dM)
    return
end
if ~(isnumeric(dM) && isequal(size(dM),size(mobile)) && all(isfinite(dM(:))))
    error('Biogeo_Lateral_Pools:BadChange', ...
        'dM must match the cell-by-species mobile-pool matrix.');
end

for c = 1:size(P_in,1)
    for s = 1:numel(pool_idx)
        idx = pool_idx{s};
        new_total = mobile(c,s) + dM(c,s);
        tol = 1e-12*max(1,mobile(c,s));
        if new_total < -tol
            error('Biogeo_Lateral_Pools:NegativePool', ...
                'Lateral transfer overdraws %s in cell %d.',names{s},c);
        elseif new_total < 0
            new_total = 0;
        end
        if isscalar(idx)
            P_out(c,idx) = new_total;
        elseif mobile(c,s) > 0
            P_out(c,idx) = P_in(c,idx) * (new_total/mobile(c,s));
        else
            P_out(c,idx) = new_total/numel(idx);
        end
    end
end

end
