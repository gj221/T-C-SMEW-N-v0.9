function rates = Biogeo_Cation_Partition_Rates(rates,opt_ERW,exchange_scheme)
%%% Select legacy or SCEPTER/Appelo cation partition ownership.
%%% rates = [exchange release, exchange uptake, fixation, fixed release], d^-1.

if opt_ERW ~= 1
    return
end
if isstring(exchange_scheme) && isscalar(exchange_scheme)
    exchange_scheme = char(exchange_scheme);
end
if ~ischar(exchange_scheme) || ...
        ~any(strcmpi(exchange_scheme,{'TC_GT','appelo'}))
    error('BIOGEO:InvalidExchangeScheme', ...
        'Exchange scheme must be TC_GT or appelo.');
end
if strcmpi(exchange_scheme,'appelo')
    % Appelo already partitions the conserved solution-exchange inventory.
    rates(:) = 0;
else
    rates(1:2) = 0;
end

end
