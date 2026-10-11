function pH_mean = pH_Area_Weighted(pH_cell,weights)
%%% Area-weighted hydrogen-activity pH index.  This is exact for the
%%% conservative porewater H activity mean.  For bulk-extract pH it is a
%%% reporting index only; an exact pooled extract requires pooled chemical
%%% inventories followed by the extraction solve.

if ~isnumeric(pH_cell) || ~isreal(pH_cell) || any(~isfinite(pH_cell(:)))
    error('ERW:InvalidPH', 'Cell pH values must be finite and real.');
end
if ~isnumeric(weights) || ~isreal(weights) || any(~isfinite(weights(:))) || ...
        any(weights(:) < 0) || sum(weights(:)) <= 0
    error('ERW:InvalidAreaWeights', ...
        'Area weights must be finite, nonnegative and have a positive sum.');
end

weights = weights(:)/sum(weights(:));
if isvector(pH_cell)
    if numel(weights) == 1
        pH_cell = pH_cell(:);
    elseif numel(pH_cell) == numel(weights)
        pH_cell = reshape(pH_cell,1,[]);
    end
end
if size(pH_cell,2) ~= numel(weights)
    error('ERW:InvalidAreaWeights', ...
        'The number of area weights must match the pH cell columns.');
end

pH_mean = -log10(10.^(-pH_cell)*weights);
end
