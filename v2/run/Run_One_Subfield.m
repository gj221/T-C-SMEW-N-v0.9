function result_file = Run_One_Subfield(run_id, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run_One_Subfield
%%
%% Run one of the twelve simulations defined by Run_Matrix_12 and save
%% the result. Runs 1 to 9 use the retained sub-field domain. Runs 10 to
%% 12 use one plot-scale column and the original vertical transport only.
%%
%% This is the 2-D counterpart of ERW_run: it performs the same workbook
%% read, spin-up load, forcing load and window slice, then partitions the
%% plot into sub-plot cells and hands the whole set to
%% MAIN_FRAME_BG_ERW_2D instead of the single-column frame.
%%
%% TREATMENT PLACEMENT
%%   Each cell's ColCfg entry carries apply_rock, set from the run's
%%   layout and the cell's region.  Cells that receive no amendment still
%%   run the full carbonate and cation-exchange chemistry, so an
%%   untreated field interior and a treated margin are directly
%%   comparable inside one simulation.
%%
%% Usage
%%   Run_One_Subfield(5)                       % run 5 of the matrix
%%   Run_One_Subfield(11)                      % basalt, one-dimensional
%%   Run_One_Subfield(5,'NDays',3)             % short smoke test
%%   Run_One_Subfield(12,'StartDate','2016-11-09','NDays',3)
%%   Run_One_Subfield(5,'GridNx',1,'GridNy',1) % coarser main grid
%%   Run_One_Subfield(11,'AppeloAlpha',0.8)     % calibration trajectory
%%   Run_One_Subfield(11,'AppeloLogKH',[5.9 4.8 10.47 10.986 16.47])
%%   Run_One_Subfield(10,'Pcla',0.12,'Psan',0.22,'KNmaxd',0.0032)
%%   Run_One_Subfield(2,'ParallelWorkers',6)    % parallel daily ERW columns
%%
%% Returns the path of the saved result file.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ResultDirectory defaults to the usual SIM/Subfield_12, so behaviour is
%% unchanged when it is not given. It exists because two runs of the same
%% run_id launched concurrently would otherwise race on one output path.
opt = struct('NDays',[],'StartDate',[],'GridNx',[],'GridNy',[], ...
    'AppeloAlpha',[],'AppeloLogKH',[],'Pcla',[],'Psan',[], ...
    'KNmaxd',[],'ParallelWorkers',[], ...
    'SaveResults',true,'ResultDirectory','');
if mod(numel(varargin),2) ~= 0
    error('Run_One_Subfield:BadOptions','Options must be name/value pairs.');
end
for k = 1:2:numel(varargin)
    name = char(varargin{k});
    if ~isfield(opt,name)
        error('Run_One_Subfield:BadOptions','Unknown option %s.',name);
    end
    opt.(name) = varargin{k+1};
end
validate_optional_integer(opt.NDays,'NDays');
validate_optional_integer(opt.GridNx,'GridNx');
validate_optional_integer(opt.GridNy,'GridNy');
validate_optional_nonnegative_scalar(opt.AppeloAlpha,'AppeloAlpha');
validate_optional_fraction(opt.Pcla,'Pcla');
validate_optional_fraction(opt.Psan,'Psan');
validate_optional_nonnegative_scalar(opt.KNmaxd,'KNmaxd');
validate_optional_nonnegative_integer(opt.ParallelWorkers,'ParallelWorkers');
if xor(isempty(opt.Pcla),isempty(opt.Psan))
    error('Run_One_Subfield:BadOptions', ...
        'Pcla and Psan must be supplied together.');
end
if ~isempty(opt.AppeloLogKH) && ~(isnumeric(opt.AppeloLogKH) && ...
        isvector(opt.AppeloLogKH) && numel(opt.AppeloLogKH) == 5 && ...
        all(isfinite(opt.AppeloLogKH)))
    error('Run_One_Subfield:BadOptions', ...
        'AppeloLogKH must be five finite values ordered [Na K Ca Mg Al].');
end
opt.AppeloLogKH = reshape(opt.AppeloLogKH,1,[]);
if ~(islogical(opt.SaveResults) && isscalar(opt.SaveResults))
    error('Run_One_Subfield:BadOptions','SaveResults must be scalar logical.');
end
start_datenum = [];
if ~isempty(opt.StartDate)
    try
        start_datenum = datenum(opt.StartDate);
    catch
        error('Run_One_Subfield:BadStartDate', ...
            'StartDate must be a scalar date accepted by datenum.');
    end
    if ~(isnumeric(start_datenum) && isscalar(start_datenum) && isfinite(start_datenum))
        error('Run_One_Subfield:BadStartDate', ...
            'StartDate must be a scalar date accepted by datenum.');
    end
end

RM = Run_Matrix_12();
if ~(isnumeric(run_id) && isscalar(run_id) && isfinite(run_id) && ...
        run_id == fix(run_id) && run_id >= 1 && run_id <= numel(RM))
    error('Run_One_Subfield:BadId','run_id must be 1..%d', numel(RM));
end

R = RM(run_id);

%% This function lives in v2/run; the repository root is one level up.
script_directory = fileparts(mfilename('fullpath'));
if isempty(script_directory); script_directory = pwd; end
current_directory = fileparts(script_directory);
Directory           = fullfile(current_directory,'T&C_Code');
Input_Directory     = fullfile(Directory,'Inputs');
Parameter_Directory = fullfile(Directory,'Parameters');
addpath(Directory, Input_Directory, Parameter_Directory);

CONTROL_CASE = 5;
NN_case_study_ERW_MOD = CONTROL_CASE;
OPT_ERW = 1;                    %% chemistry on; apply_rock places the rock
OPT_SoilBiogeochemistry = 1;

fprintf('\n=================================================\n');
fprintf(' Run %d: %s\n', R.id, R.name);
fprintf('   treatment %s (workbook case %d)\n', R.treatment, R.case_number);
fprintf('   layout    %s\n', R.layout);
fprintf('=================================================\n');

%% ---------------------------------------------------------------------
%% Workbook
%% ---------------------------------------------------------------------
wb = fullfile(Input_Directory,'ERW_MOD_PARAM_parameters.xlsx');
MOD_opts = detectImportOptions(wb,"Sheet",'MOD_PARAM');
MOD_opts.DataRange = "A2:Y6";
MOD_opts.VariableTypes = ["double","string","string","double","double", ...
    "string","double","double","string","double", ...
    "double","double","double","double","string", ...
    "double","double","string","double","double", ...
    "double","double","double","string","string"];
MOD_PARAM_input = readtable(wb,MOD_opts);

ERW_opts = detectImportOptions(wb,"Sheet",'ERW_MOD');
ERW_opts.DataRange = "A2:AI6";
ERW_opts.VariableTypes = {'double','string','string','double','double',...
    'string','double','double','string','double',...
    'string','string','string','string','double', ...
    'double','double','double','double','double',...
    'double','double','double','string','double',...
    'string','string','double','double','double',...
    'double','double','double','double','double'};
ERW_opts = setvaropts(ERW_opts,{'f_CEC_in','rock_filepath','mineral','rock_f_in'},'FillValue','NaN');
ERW_PARAM_input = readtable(wb,ERW_opts);

[MOD_PARAM_input,ERW_PARAM_input,ERW_PARAM_control] = ...
    Select_Subfield_Parameter_Rows( ...
    MOD_PARAM_input,ERW_PARAM_input,R.case_number);
MOD_PARAM_input.Zbio(:) = R.biogeo_depth_mm;
if ~isempty(opt.Pcla)
    texture_sum = opt.Pcla + opt.Psan + MOD_PARAM_input.Porg;
    if opt.Pcla + opt.Psan <= 0 || any(~isfinite(texture_sum)) || ...
            any(texture_sum > 1)
        error('Run_One_Subfield:BadOptions', ...
            ['Pcla and Psan must have a positive sum, and Pcla, Psan ' ...
             'and Porg must sum to no more than one.']);
    end
    MOD_PARAM_input.Pcla(:) = opt.Pcla;
    MOD_PARAM_input.Psan(:) = opt.Psan;
end
id_location = char(MOD_PARAM_input.Location);
fprintf('Active biogeochemical depth: %.0f mm\n',R.biogeo_depth_mm);
fprintf('Texture: Pcla %.4f, Psan %.4f\n', ...
    MOD_PARAM_input.Pcla,MOD_PARAM_input.Psan);

default_biogeo = Biogeochemistry_Parameter(1);
effective_KNmaxd = default_biogeo.KNmaxd;
if ~isempty(opt.KNmaxd)
    effective_KNmaxd = opt.KNmaxd;
end
fprintf('KNmaxd: %.8g d^-1\n',effective_KNmaxd);

%% Spin-up pool state
load(fullfile(Input_Directory,['Data_' id_location '_BG.mat']));  %#ok<LOAD>

%% ---------------------------------------------------------------------
%% Forcing
%% ---------------------------------------------------------------------
NN  = MOD_PARAM_input.NumberOfDays*24;
dt  = 3600;
dth = 1;
eval(append('Zs','=',eval('MOD_PARAM_input.Zs'),';'));
ms = length(Zs)-1;
cc = MOD_PARAM_input.cc;

load(fullfile(Input_Directory,['Data_' id_location '.mat']));     %#ok<LOAD>

x1 = 61351;                     %% 01-Jan-2016
if ~isempty(start_datenum)
    [date_error,x1] = min(abs(Date-start_datenum));
    if isempty(x1) || date_error > 1/48
        error('Run_One_Subfield:BadStartDate', ...
            'StartDate is outside the available hourly forcing record.');
    end
end
NN = 1826*24;
if ~isempty(opt.NDays)
    NN = opt.NDays*24;
end
x2 = x1+NN-1;
if x2 > numel(Date)
    error('Run_One_Subfield:BadWindow', ...
        'Requested StartDate and NDays exceed the forcing record.');
end

Date=Date(x1:x2); Pr=Pr(x1:x2); Ta=Ta(x1:x2); Ws=Ws(x1:x2); ea=ea(x1:x2);
SAD1=SAD1(x1:x2); SAD2=SAD2(x1:x2); SAB1=SAB1(x1:x2); SAB2=SAB2(x1:x2);
Pre=Pre(x1:x2); N=N(x1:x2); Tdew=Tdew(x1:x2); esat=esat(x1:x2);
PARB=PARB(x1:x2); PARD=PARD(x1:x2);

%% Kantola's first amendment event is 10-Nov-2016. A shortened smoke
%% window ending before that date contains no rock application, so those
%% cells must initialise from the control row. The full run retains the
%% treatment row because its future application events are in the window.
treatment_in_window = strcmp(R.treatment,'control') || ...
    Date(end) >= datenum(2016,11,10);

t_bef = -0.67; t_aft = 1.67;
Ds = esat-ea; Ds(Ds<0) = 0;
Ca_file = load(fullfile(Input_Directory,'Ca_Data.mat'));
d1 = find(abs(Ca_file.Date_CO2-Date(1))<1/36);
d2 = find(abs(Ca_file.Date_CO2-Date(end))<1/36);
Ca = Ca_file.Ca(d1:d2);
Oa = 210000;
Ws(Ws<=0) = 0.01;
%% datevec returns outputs shaped like its input, so a row-vector Date
%% would make [YE MO DA HO] a 1 x 4N row instead of an N x 4 matrix.
%% Force columns, as the single-column frame does.
[YE,MO,DA,HO] = datevec(Date);
Datam = [YE(:) MO(:) DA(:) HO(:)];
clear YE MO DA HO
NNd = ceil(NN/24)+1;

%% ---------------------------------------------------------------------
%% Plot partition and per-column configuration
%% ---------------------------------------------------------------------
g = R.geom;
nx = g.main_grid_nx; ny = g.main_grid_ny;
if R.use_subdomains
    if ~isempty(opt.GridNx); nx = opt.GridNx; end
    if ~isempty(opt.GridNy); ny = opt.GridNy; end
    topo = Build_Cell_Topology(g.field_side_m,g.margin_width_m,nx,ny);
else
    topo = Build_Single_Cell_Topology(g.field_side_m^2);
end
parallel_workers = resolve_parallel_workers(opt.ParallelWorkers,topo.N_cells);
run_provenance.matlab_version = version;
run_provenance.matlab_release = version('-release');
run_provenance.computer = computer;
run_provenance.blas = regexprep(strtrim(evalc('version -blas')),'\s+',' ');
run_provenance.max_threads = maxNumCompThreads;
run_provenance.parallel_workers = parallel_workers;
run_provenance.numerical_identity = getenv('TC_NUMERICAL_IDENTITY');
%% The campaign gate that this process passed before it started.  Filename
%% wildcards cannot tell which code produced a MAT, so the stamp travels
%% with the result and a later job can reject anything validated
%% differently instead of mixing two model versions in one campaign.
run_provenance.validation_fingerprint = getenv('TC_VALIDATION_FINGERPRINT');
run_provenance.campaign = getenv('TC_CAMPAIGN_TAG');
run_provenance.slurm_job_id = getenv('SLURM_JOB_ID');
run_provenance.slurm_step_id = getenv('SLURM_STEP_ID');
slurm_restart_count = getenv('SLURM_RESTART_COUNT');
if isempty(slurm_restart_count); slurm_restart_count = '0'; end
run_provenance.slurm_restart_count = slurm_restart_count;
run_provenance.slurm_array_task_id = getenv('SLURM_ARRAY_TASK_ID');
run_provenance.node = getenv('SLURMD_NODENAME');
run_provenance.date_start = Date(1);
run_provenance.date_end = Date(end);
run_provenance.n_days = NN/24;
run_provenance.n_timesteps = NN;
run_provenance.n_columns = topo.N_cells;
run_provenance.result_directory = char(opt.ResultDirectory);
run_provenance.Pcla = double(MOD_PARAM_input.Pcla);
run_provenance.Psan = double(MOD_PARAM_input.Psan);
run_provenance.KNmaxd = effective_KNmaxd;

base_param = fullfile(Parameter_Directory,[char(MOD_PARAM_input.MOD_PARAM) '.m']);
is_main = strcmp(topo.cell_region,'main');
ColCfg = repmat(struct('PARAM_IC',base_param,'PFT_OVERRIDE','', ...
                       'apply_rock',false,'ERW_PARAM',ERW_PARAM_control, ...
                       'appelo_alpha',opt.AppeloAlpha, ...
                       'appelo_logKH',opt.AppeloLogKH, ...
                       'KNmaxd',opt.KNmaxd), ...
                topo.N_cells,1);
for c = 1:topo.N_cells
    if is_main(c)
        ColCfg(c).apply_rock = R.apply_main;
    else
        ColCfg(c).apply_rock = R.apply_margin;
        ColCfg(c).PFT_OVERRIDE = R.margin_pft;   %% '' keeps the crop
    end
    if ColCfg(c).apply_rock && treatment_in_window
        ColCfg(c).ERW_PARAM = ERW_PARAM_input;
    else
        ColCfg(c).ERW_PARAM = ERW_PARAM_control;
    end
end
if R.use_subdomains
    fprintf('Domain: %d columns, %.0f m2 main + %.0f m2 margin (%g m wide)\n', ...
        topo.N_cells,sum(topo.cell_area_m2(is_main)), ...
        sum(topo.cell_area_m2(~is_main)),g.margin_width_m);
    fprintf('Rock on: %d main cell(s), %d margin cell(s)\n', ...
        sum([ColCfg(is_main).apply_rock]),sum([ColCfg(~is_main).apply_rock]));
else
    fprintf('Domain: one %.0f m2 column; vertical transport only\n', ...
        topo.cell_area_m2);
    fprintf('Rock on: %d column(s)\n',sum([ColCfg.apply_rock]));
end
if ~treatment_in_window && ~strcmp(R.treatment,'control')
    fprintf('No amendment event falls inside this shortened window.\n');
end

%% ---------------------------------------------------------------------
%% Run
%% ---------------------------------------------------------------------
%% Rolling checkpoint for this run.  One file per run, overwritten in
%% place, so a long column that aborts late (or a machine that sleeps)
%% resumes instead of replaying from day 1.  MAIN_FRAME refuses to resume
%% a checkpoint whose code fingerprint has moved on unless
%% ERW_ckpt_force_resume is set, so an edit to the chemistry starts clean.
ckpt_dir = char(opt.ResultDirectory);
if isempty(ckpt_dir)
    ckpt_dir = fullfile(current_directory,'SIM','Subfield_12');
end
ckpt_dir = fullfile(ckpt_dir,'checkpoints');
if ~exist(ckpt_dir,'dir'); mkdir(ckpt_dir); end
ckpt_tag = R.name;
if ~isempty(opt.AppeloAlpha)
    ckpt_tag = sprintf('%s__alpha%s',ckpt_tag,number_tag(opt.AppeloAlpha));
end
if ~isempty(opt.AppeloLogKH)
    logkh_parts = arrayfun(@number_tag,opt.AppeloLogKH,'UniformOutput',false);
    ckpt_tag = sprintf('%s__logkh%s',ckpt_tag,strjoin(logkh_parts,'_'));
end
if ~isempty(opt.Pcla)
    ckpt_tag = sprintf('%s__pcla%s__psan%s',ckpt_tag, ...
        number_tag(opt.Pcla),number_tag(opt.Psan));
end
if ~isempty(opt.KNmaxd)
    ckpt_tag = sprintf('%s__knmaxd%s',ckpt_tag,number_tag(opt.KNmaxd));
end
ckpt_tag = sprintf('%s__%dd__from%s',ckpt_tag,NN/24,datestr(Date(1),'yyyymmdd'));
ERW_ckpt_file = fullfile(ckpt_dir, ...
    sprintf('ckpt_run%02d_%s.mat',R.id,ckpt_tag));
%% Interval comes from the environment so it reaches the worker
%% processes the launcher spawns (a launcher-workspace variable does
%% not).  TC_CKPT_MINUTES = 0 disables checkpointing entirely.
if ~exist('ERW_ckpt_minutes','var') || isempty(ERW_ckpt_minutes)
    ERW_ckpt_minutes = str2double(getenv('TC_CKPT_MINUTES'));
    if ~isfinite(ERW_ckpt_minutes) || ERW_ckpt_minutes < 0
        ERW_ckpt_minutes = 10;
    end
end
if ERW_ckpt_minutes == 0
    ERW_ckpt_file = '';
end
%% Debugging escape hatch, also via the environment so it reaches the
%% worker.  Resumes a checkpoint whose code fingerprint has moved on:
%% legitimate when iterating on a LATE-timestep failure whose earlier
%% trajectory the edit cannot touch, and never for production results,
%% because the run then contains two model versions spliced at the
%% resume point.  MAIN_FRAME warns and names that iteration.
if ~exist('ERW_ckpt_force_resume','var') || isempty(ERW_ckpt_force_resume)
    ERW_ckpt_force_resume = strcmpi(strtrim(getenv('TC_CKPT_FORCE_RESUME')), '1');
end

parallel_pool_created = start_parallel_pool(parallel_workers);
old = pwd;
cleanup = onCleanup(@() restore_run_environment(old,parallel_pool_created));
cd(Directory);
%% The deficit counters are persistent, so clear them here: several runs can
%% share one MATLAB process and each is reported separately.
Plant_Carbon_Deficit('reset');
MAIN_FRAME_BG_ERW_2D;
clear cleanup

%% ---------------------------------------------------------------------
%% Carbohydrate reserve overdraws
%% ---------------------------------------------------------------------
%% Reported, not corrected. Legacy carried a negative B(4) forward and the
%% clamped ratio in Plant_Exports absorbed it, so this carbon never appeared
%% in any budget. The campaign conservation checks audit ERW mass and
%% soil-atmosphere exchange, not the plant carbohydrate pool, so without this
%% line the overdraw stays invisible. A single event of order 1e-3 gC/m2 is a
%% boundary artifact of the adaptive step; a count that scales with the run
%% length is a systematic sink and the nitrogen response needs revisiting.
carbon_deficit = Plant_Carbon_Deficit('report');
if carbon_deficit.reserve_events > 0 || carbon_deficit.other_pool_events > 0
    fprintf(['Carbon reserve overdraw: %d events, %.6g gC/m2 total, ' ...
        '%.6g gC/m2 worst; %d other negative pool states.\n'], ...
        carbon_deficit.reserve_events,carbon_deficit.reserve_total_gC_m2, ...
        carbon_deficit.reserve_worst_gC_m2,carbon_deficit.other_pool_events);
else
    fprintf('Carbon reserve overdraw: none.\n');
end

%% ---------------------------------------------------------------------
%% Save
%% ---------------------------------------------------------------------
result_file = '';
if opt.SaveResults
    fl_save = char(opt.ResultDirectory);
    if isempty(fl_save)
        fl_save = fullfile(current_directory,'SIM','Subfield_12');
    end
    if ~exist(fl_save,'dir'); mkdir(fl_save); end
    tag = R.name;
    if R.use_subdomains && (~isempty(opt.GridNx) || ~isempty(opt.GridNy))
        tag = sprintf('%s__grid%dx%d',tag,nx,ny);
    end
    if ~isempty(opt.AppeloAlpha)
        tag = sprintf('%s__alpha%s',tag,number_tag(opt.AppeloAlpha));
    end
    if ~isempty(opt.AppeloLogKH)
        logkh_parts = arrayfun(@number_tag,opt.AppeloLogKH,'UniformOutput',false);
        tag = sprintf('%s__logkh%s',tag,strjoin(logkh_parts,'_'));
    end
    if ~isempty(opt.Pcla)
        tag = sprintf('%s__pcla%s__psan%s',tag, ...
            number_tag(opt.Pcla),number_tag(opt.Psan));
    end
    if ~isempty(opt.KNmaxd)
        tag = sprintf('%s__knmaxd%s',tag,number_tag(opt.KNmaxd));
    end
    if ~isempty(opt.NDays); tag = sprintf('%s__%dd',tag,opt.NDays); end
    if ~isempty(start_datenum); tag = sprintf('%s__from%s',tag,datestr(start_datenum,'yyyymmdd')); end
    result_file = fullfile(fl_save,sprintf('TC_2D_run%02d_%s.mat',R.id,tag));
    tmp_result = [tempname(fl_save) '.mat'];
    try
        save(tmp_result,'-v7.3');
        movefile(tmp_result,result_file,'f');
    catch ME
        if isfile(tmp_result); delete(tmp_result); end
        rethrow(ME)
    end
    fprintf('Saved: %s\n', result_file);
    %% The run finished and its results are on disk, so the rolling
    %% checkpoint is now dead weight.  Drop it: these files are large
    %% and keeping one per completed run would fill the volume.
    if exist('ERW_ckpt_file','var') && ~isempty(ERW_ckpt_file) && isfile(ERW_ckpt_file)
        delete(ERW_ckpt_file);
    end
end

end

function validate_optional_integer(value,name)
if isempty(value); return; end
if ~(isnumeric(value) && isscalar(value) && isfinite(value) && ...
        value == fix(value) && value >= 1)
    error('Run_One_Subfield:BadOptions','%s must be a positive integer.',name);
end
end

function validate_optional_nonnegative_scalar(value,name)
if isempty(value); return; end
if ~(isnumeric(value) && isscalar(value) && isfinite(value) && value >= 0)
    error('Run_One_Subfield:BadOptions', ...
        '%s must be a finite nonnegative scalar.',name);
end
end

function validate_optional_fraction(value,name)
if isempty(value); return; end
if ~(isnumeric(value) && isreal(value) && isscalar(value) && ...
        isfinite(value) && value >= 0 && value <= 1)
    error('Run_One_Subfield:BadOptions', ...
        '%s must be a finite scalar in [0,1].',name);
end
end

function validate_optional_nonnegative_integer(value,name)
if isempty(value); return; end
if ~(isnumeric(value) && isscalar(value) && isfinite(value) && ...
        value == fix(value) && value >= 0)
    error('Run_One_Subfield:BadOptions', ...
        '%s must be a nonnegative integer.',name);
end
end

function workers = resolve_parallel_workers(requested,n_cells)
if isempty(requested)
    text = strtrim(getenv('TC_PARALLEL_WORKERS'));
    if isempty(text)
        requested = 0;
    else
        requested = str2double(text);
    end
end
validate_optional_nonnegative_integer(requested,'ParallelWorkers');
workers = min(requested,n_cells);
if workers < 2; workers = 0; end
end

function created = start_parallel_pool(workers)
created = false;
if workers == 0; return; end
if isempty(ver('parallel')) || ~license('test','Distrib_Computing_Toolbox')
    error('Run_One_Subfield:MissingParallelToolbox', ...
        'ParallelWorkers requires Parallel Computing Toolbox.');
end
pool = gcp('nocreate');
if isempty(pool)
    pool = parpool('Threads',workers);
    created = true;
elseif ~isa(pool,'parallel.ThreadPool')
    error('Run_One_Subfield:WrongParallelPool', ...
        'ParallelWorkers requires a thread-based pool.');
elseif pool.NumWorkers < workers
    error('Run_One_Subfield:SmallParallelPool', ...
        'The existing thread pool has %d workers; %d are required.', ...
        pool.NumWorkers,workers);
end
fprintf('Daily ERW chemistry: %d parallel workers\n',workers);
clear pool
end

function restore_run_environment(old,delete_pool)
cd(old);
if delete_pool
    pool = gcp('nocreate');
    if ~isempty(pool); delete(pool); end
end
end

function tag = number_tag(value)
tag = strrep(sprintf('%.17g',value),'.','p');
tag = strrep(tag,'-','m');
tag = strrep(tag,'+','');
end
