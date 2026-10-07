function workers = Resolve_2D_Parfor_Workers(n_columns)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Resolve the bounded worker count for independent daily ERW columns.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~(isnumeric(n_columns) && isscalar(n_columns) && isfinite(n_columns) && ...
        n_columns == fix(n_columns) && n_columns >= 1)
    error('Resolve_2D_Parfor_Workers:BadColumnCount', ...
        'n_columns must be a positive integer.');
end

raw = strtrim(getenv('TC_PARALLEL_WORKERS'));
if isempty(raw)
    workers = 0;
    return
end

workers = str2double(raw);
if ~(isfinite(workers) && workers == fix(workers) && workers >= 0)
    error('Resolve_2D_Parfor_Workers:BadEnvironment', ...
        'TC_PARALLEL_WORKERS must be a nonnegative integer.');
end
workers = min(workers,n_columns);
if workers < 2; workers = 0; end

end
