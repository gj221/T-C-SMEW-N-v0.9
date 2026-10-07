function Run_Worker_Scaling_Test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Measure what actually limits a 2-D run, before committing to 35-year
%% campaigns.  Three arms, one array job:
%%
%%   1   6 workers, checkpoints every 10 min   (production settings)
%%   2  12 workers, checkpoints every 10 min
%%   3   6 workers, checkpoints disabled
%%
%% Arm 1 against arm 2 is the worker speedup.  The daily ERW chemistry is
%% the only parallel section (parfor over 12 columns in
%% MAIN_FRAME_BG_ERW_2D), and 6 workers means every worker handles two
%% columns, so 12 should halve that section -- but nothing here assumes how
%% large a share of runtime it is.  The completed five-year runs argue it
%% may be small: run 1 is the control, with no mineral dissolution at all,
%% and it finished only 9% faster than run 2.
%%
%% Arm 1 against arm 3 is the cost of checkpointing, which production pays
%% 258 times over a 43-hour run.
%%
%% MaxRSS from all three, set beside the known five-year figures (17.3 GB
%% control, 24.2 GB basalt), gives a second point on the memory curve.  That
%% matters more than the timing: extrapolating linearly puts a 35-year run
%% at 120-170 GB against the 48 G the campaign currently requests, so the
%% production runs would be killed part-way through.  Two points separate
%% fixed cost from per-day growth and size the request on evidence.
%%
%% Run 2 (basalt, whole_field) is used throughout: it is the heavier of the
%% two measured runs, and sizing should follow the worst case.
%%
%% This is a diagnostic and is not part of any campaign: no controller
%% token, no validation stamp, no provenance gate.  Its results go to their
%% own SIM subtree and can be deleted with this directory once the numbers
%% are recorded.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NDAYS = 365;
RUN_ID = 2;

arm = read_positive_integer('SLURM_ARRAY_TASK_ID');
switch arm
    case 1
        workers = 6;  ckpt_minutes = 10;
    case 2
        workers = 12; ckpt_minutes = 10;
    case 3
        workers = 6;  ckpt_minutes = 0;
    otherwise
        error('Worker_Scaling_Test:BadArm', ...
            'SLURM_ARRAY_TASK_ID must be 1, 2 or 3.');
end

test_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(test_dir));
one_d_dir = fullfile(repo_root,'runtime', ...
    'Base_1D');
run_dir = fullfile(repo_root,'v2','run');
code_dir = fullfile(repo_root,'v2','T&C_Code');
addpath(run_dir,code_dir,fullfile(code_dir,'Inputs'), ...
    fullfile(code_dir,'Parameters'),one_d_dir);

matlab_release = version('-release');
if ~strcmp(matlab_release,'2025a')
    error('Worker_Scaling_Test:WrongMatlabRelease', ...
        'Timings are only comparable with production under R2025a; loaded R%s.', ...
        matlab_release);
end
Bouchet_Single_Thread_Runtime('Worker_Scaling_Test');

allocated_cpus = read_positive_integer('SLURM_CPUS_PER_TASK');
if allocated_cpus < 12
    error('Worker_Scaling_Test:InsufficientCpus', ...
        ['Reserve 12 CPUs for every arm, so that worker count is the only ' ...
        'difference between them; got %d.'],allocated_cpus);
end
%% Run_One_Subfield reads TC_CKPT_MINUTES from the environment so the value
%% reaches worker processes. It is set here rather than in the submitter
%% because it varies by arm.
setenv('TC_CKPT_MINUTES',sprintf('%d',ckpt_minutes));
setenv('TC_CAMPAIGN_TAG','Worker_Scaling_Test');

P = Campaign_Parameters;
results_dir = fullfile(repo_root,'v2','SIM','Worker_Scaling_Test', ...
    sprintf('arm%d_w%02d_ckpt%02d',arm,workers,ckpt_minutes));
if ~exist(results_dir,'dir'); mkdir(results_dir); end

fprintf(['Arm %d: %d workers, checkpoints %s\n' ...
    'Run %d, %d days from %s, 12 columns\n' ...
    'Node %s, %d CPUs allocated\n'], ...
    arm,workers, ...
    ternary(ckpt_minutes==0,'disabled',sprintf('every %d min',ckpt_minutes)), ...
    RUN_ID,NDAYS, ...
    char(datetime(P.StartDate,'ConvertFrom','datenum','Format','yyyy-MM-dd')), ...
    getenv('SLURMD_NODENAME'),allocated_cpus);

wall = tic;
result_file = Run_One_Subfield(RUN_ID, ...
    'NDays',NDAYS,'StartDate',P.StartDate, ...
    'GridNx',2,'GridNy',2, ...
    'AppeloAlpha',P.AppeloAlpha,'AppeloLogKH',P.AppeloLogKH, ...
    'KNmaxd',P.KNmaxd,'Kbot',P.Kbot, ...
    'ParallelWorkers',workers,'SaveResults',true, ...
    'ResultDirectory',results_dir);
elapsed = toc(wall);

fprintf('\n===== ARM %d RESULT =====\n',arm);
fprintf('  workers      %d\n',workers);
fprintf('  checkpoints  %s\n', ...
    ternary(ckpt_minutes==0,'disabled',sprintf('%d min',ckpt_minutes)));
fprintf('  elapsed      %.1f s  (%.2f h)\n',elapsed,elapsed/3600);
fprintf('  per day      %.3f s\n',elapsed/NDAYS);
%% The projection horizon is 12784 days. Scaling by days is the honest
%% extrapolation: the cost per day is dominated by the hourly loop and the
%% daily ERW call, neither of which depends on how far the run has got.
fprintf('  -> 12784 d   %.1f h (linear)\n',elapsed/NDAYS*12784/3600);
fprintf('  result       %s\n',result_file);
d = dir(result_file);
if ~isempty(d)
    fprintf('  result size  %.1f MB\n',d(1).bytes/1024^2);
end
fprintf('=========================\n');
fprintf(['Read MaxRSS from sacct; it is not visible from inside MATLAB.\n' ...
    'Compare arms 1 and 2 for the worker speedup, 1 and 3 for the\n' ...
    'checkpoint cost, and diff the arm 1 and arm 2 result MATs: worker\n' ...
    'count must not change the numbers.\n']);
end

function out = ternary(cond,a,b)
if cond; out = a; else; out = b; end
end

function value = read_positive_integer(name)
raw = strtrim(getenv(name));
value = str2double(raw);
if isempty(raw) || ~(isfinite(value) && value == fix(value) && value >= 1)
    error('Worker_Scaling_Test:BadEnvironment', ...
        '%s must be a positive integer.',name);
end
end
