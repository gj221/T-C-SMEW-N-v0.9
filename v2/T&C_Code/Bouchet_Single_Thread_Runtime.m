function observed = Bouchet_Single_Thread_Runtime(error_id_prefix)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pin the campaign to single-threaded numerics and report what was found.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Multithreaded BLAS reorders floating-point reductions, so a rerun of the
%% same campaign would not reproduce bit for bit. The submission scripts
%% export OMP_NUM_THREADS, MKL_NUM_THREADS and OPENBLAS_NUM_THREADS as 1,
%% but MATLAB manages its own BLAS threading and does not guarantee that
%% those variables survive startup unchanged: on Bouchet R2025a,
%% MKL_NUM_THREADS does not arrive as 1 even though OMP_NUM_THREADS does.
%%
%% Asserting on the variables therefore tests something MATLAB owns rather
%% than something that governs the numerics. What governs the numerics is
%% maxNumCompThreads: with -singleCompThread MATLAB runs its BLAS on one
%% thread whatever the variables say. So this function
%%
%%   - hard-fails unless maxNumCompThreads is 1, which is the real guarantee
%%   - sets the three variables to 1 rather than asserting them, so any
%%     child process, parpool worker or external library inherits the same
%%     constraint instead of picking up an allocation-sized default
%%   - returns what it observed, for the run provenance
%%
%% INPUT
%%   error_id_prefix  identifier prefix for the calling campaign
%% OUTPUT
%%   observed         struct of the values seen before normalisation

if nargin < 1 || isempty(error_id_prefix)
    error_id_prefix = 'Bouchet';
end
error_id_prefix = char(error_id_prefix);

names = {'OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS'};
observed = struct();
for k = 1:numel(names)
    before = strtrim(getenv(names{k}));
    observed.(names{k}) = before;
    if ~strcmp(before,'1')
        fprintf(['Thread environment: %s was {%s}; setting it to 1 for ' ...
            'this process and its children.\n'],names{k},before);
        setenv(names{k},'1');
    end
end

%% The authoritative check. MATLAB's own computational thread count is what
%% decides whether its BLAS reductions are deterministic.
thread_count = maxNumCompThreads;
observed.maxNumCompThreads = thread_count;
if thread_count ~= 1
    error([error_id_prefix ':MatlabThreads'], ...
        ['MATLAB must run on one computational thread; found %d. Start ' ...
         'it with -singleCompThread.'],thread_count);
end

blas = regexprep(strtrim(evalc('version -blas')),'\s+',' ');
observed.blas = blas;
fprintf('Computational threads: %d\nBLAS: %s\n',thread_count,blas);

end
