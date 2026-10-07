function idx = ERW_Application_Indices(time_application,dtw,n_steps)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ERW_Application_Indices
%%
%% ERW row one is the initial state at time zero. Row k is timestamped
%% (k-1)*dtw, so an application at time ta belongs to ta/dtw+1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

assert(isnumeric(time_application) && all(isfinite(time_application(:))), ...
    'ERW:applicationTime','Rock-application times must be finite and numeric.');
assert(isscalar(dtw) && isfinite(dtw) && dtw > 0, ...
    'ERW:applicationTimeStep','The ERW timestep must be finite and positive.');
assert(isscalar(n_steps) && isfinite(n_steps) && ...
    n_steps == floor(n_steps) && n_steps >= 1, ...
    'ERW:applicationStepCount', ...
    'The ERW schedule length must be a positive integer.');

grid_time = time_application/dtw;
if any(abs(grid_time-round(grid_time)) > 1e-10)
    error('ERW:applicationGrid', ...
        'Rock-application times must lie on the ERW timestep grid.');
end
idx = round(grid_time)+1;
idx = idx(idx >= 1 & idx <= n_steps);
end
