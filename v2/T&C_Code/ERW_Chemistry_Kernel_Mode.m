function mode = ERW_Chemistry_Kernel_Mode(bulk_pH_params)
%%% Which implementation of the ERW chemistry kernels a step runs.
%%%
%%%   'matlab'  the original code; the default, so a run that sets nothing
%%%             computes exactly as before
%%%   'fast'    the restructured MATLAB kernels (ERW_Appelo_Seed_Kernel,
%%%             ERW_Charge_Root_Kernel) and a one-pass ERW_data; the same
%%%             floating-point operations in the same order, so every output
%%%             is bit-identical (test_erw_chemistry_kernels)
%%%   'mex'     the same kernels compiled with MATLAB Coder
%%%             (Build_ERW_Chemistry_Mex).  MATLAB refuses MEX functions on
%%%             thread-based workers, so a run whose daily ERW columns run on
%%%             a parpool('Threads') pool must use 'fast'; 'mex' serves the
%%%             1-D frame and serial or process-pool 2-D runs.
%%%
%%% bulk_pH_params.erw_kernel decides when it is present and nonempty;
%%% otherwise the environment variable TC_ERW_KERNEL does; otherwise
%%% 'matlab'.  getenv is available on thread-based workers.

if isstruct(bulk_pH_params) && isfield(bulk_pH_params,'erw_kernel') && ...
        ~isempty(bulk_pH_params.erw_kernel)
    value = bulk_pH_params.erw_kernel;
else
    value = getenv('TC_ERW_KERNEL');
end
if isempty(value)
    mode = 'matlab';
    return
end
if isstring(value) && isscalar(value)
    value = char(value);
end
if ~ischar(value)
    error('ERW:invalidKernel', ...
        'The ERW chemistry kernel must be matlab, fast or mex.');
end
mode = lower(strtrim(value));
switch mode
    case {'matlab','fast'}
    case 'mex'
        ERW_Chemistry_Mex_Check();
    otherwise
        error('ERW:invalidKernel', ...
            'The ERW chemistry kernel must be matlab, fast or mex, not %s.',value);
end
end
