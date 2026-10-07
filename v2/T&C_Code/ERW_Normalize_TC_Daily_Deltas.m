function rates = ERW_Normalize_TC_Daily_Deltas(daily_deltas,window_days)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Convert complete T&C daily inventory changes to rates over the ERW
%% integration window so their substep integral preserves each inventory.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~(isscalar(window_days) && isfinite(window_days) && window_days > 0)
    error('ERW:invalidIntegrationWindow', ...
        'The ERW integration window must be positive and finite.');
end
if any(~isfinite(daily_deltas(:)))
    error('ERW:invalidDailyDelta', ...
        'T&C daily inventory changes must be finite.');
end
rates = daily_deltas./window_days;
end
