function P = Biogeo_Pool_Update(Ptm1,dP,soil)
% Daily update of the soil pools, g m^-2; soil is ZBIOG*rsd, g soil m^-2.
% Sinks limited to a pool's content can still leave it a few ulps below
% zero once the day's terms are summed. The donor limits treat 16 eps of a
% pool per gram of soil as nothing, so within that band the pool is empty;
% a larger deficit is an overdraw those limits failed to prevent.

P = Ptm1+dP;
tol = 16*eps(max(1,abs(P)/soil))*soil;
P(P < 0 & P >= -tol) = 0;
if any(P < 0)
    [worst,k] = min(P);
    error('BIOGEO_UNIT:NegativePoolState', ...
        'Biogeochemistry pool %d is %.6g g m^-2 after the daily update.',k,worst);
end
end
