function OH = ERW_Derived_OH(H,k_w)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Derive hydroxide from the retained water ion product and proton pool.
%%
%% H and OH are stored in umol L^-1 and k_w in (umol L^-1)^2, so
%% k_w/H returns OH in umol L^-1 without an additional conversion.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~(isnumeric(H) && isnumeric(k_w) && isequal(size(H),size(k_w)))
    error('ERW_Derived_OH:BadInput', ...
        'H and k_w must be numeric arrays with identical dimensions.');
end
if any(~isfinite(H),'all') || any(H <= 0,'all')
    error('ERW_Derived_OH:BadHydrogen', ...
        'H must contain finite positive concentrations.');
end
if any(~isfinite(k_w),'all') || any(k_w < 0,'all')
    error('ERW_Derived_OH:BadIonProduct', ...
        'k_w must contain finite nonnegative ion products.');
end

OH = k_w./H;
if any(~isfinite(OH),'all') || any(OH < 0,'all')
    error('ERW_Derived_OH:InvalidResult', ...
        'Derived OH is nonfinite or negative.');
end

end
