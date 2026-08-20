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
%%   Run_One_Subfield(5,'GridNx',1,'GridNy',1) % coarser main grid
%%
%% Returns the path of the saved result file.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ResultDirectory defaults to the usual SIM/Subfield_12, so behaviour is
%% unchanged when it is not given. It exists because two runs of the same
%% run_id launched concurrently would otherwise race on one output path.
opt = struct('NDays',[],'GridNx',[],'GridNy',[],'SaveResults',true, ...
    'ResultDirectory','');
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
if ~(islogical(opt.SaveResults) && isscalar(opt.SaveResults))
    error('Run_One_Subfield:BadOptions','SaveResults must be scalar logical.');
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
id_location = char(MOD_PARAM_input.Location);
fprintf('Active biogeochemical depth: %.0f mm\n',R.biogeo_depth_mm);

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
NN = 1826*24;
x2 = x1+NN-1;

if ~isempty(opt.NDays)
    NN = opt.NDays*24;
    x2 = x1+NN-1;
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

base_param = fullfile(Parameter_Directory,[char(MOD_PARAM_input.MOD_PARAM) '.m']);
is_main = strcmp(topo.cell_region,'main');
ColCfg = repmat(struct('PARAM_IC',base_param,'PFT_OVERRIDE','', ...
                       'apply_rock',false,'ERW_PARAM',ERW_PARAM_control), ...
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
old = cd(Directory);
cleanup = onCleanup(@() cd(old));
MAIN_FRAME_BG_ERW_2D;
clear cleanup
cd(old);

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
    if ~isempty(opt.NDays); tag = sprintf('%s__%dd',tag,opt.NDays); end
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
end

end

function validate_optional_integer(value,name)
if isempty(value); return; end
if ~(isnumeric(value) && isscalar(value) && isfinite(value) && ...
        value == fix(value) && value >= 1)
    error('Run_One_Subfield:BadOptions','%s must be a positive integer.',name);
end
end
