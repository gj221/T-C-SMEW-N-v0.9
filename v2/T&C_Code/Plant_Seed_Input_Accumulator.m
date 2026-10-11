function total = Plant_Seed_Input_Accumulator(total,seed_per_crown,crown_fraction)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Accumulate an explicit crop-establishment boundary on ground area.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if numel(total)~=7 || numel(seed_per_crown)~=7 || ...
        any(~isfinite(total)) || any(~isfinite(seed_per_crown)) || ...
        ~isscalar(crown_fraction) || ~isfinite(crown_fraction) || ...
        crown_fraction<0 || crown_fraction>1
    error('Plant_Seed_Input_Accumulator:InvalidInput', ...
        'Seed inputs must be finite seven-element vectors and crown fraction must be in [0,1].');
end
if any(seed_per_crown<0)
    error('Plant_Seed_Input_Accumulator:NegativeSeedInput', ...
        'An explicit seed boundary cannot contain a negative element input.');
end

total=reshape(total,1,7)+crown_fraction*reshape(seed_per_crown,1,7);

end