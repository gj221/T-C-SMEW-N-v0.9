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
%%   Run_One_Subfield(11,'CECpHSensitivity',6)  % empirical TC_GT only
%%   Run_One_Subfield(11,'AppeloLogKH',[5.9 4.8 10.47 10.986 16.47])
%%   Run_One_Subfield(10,'Pcla',0.12,'Psan',0.22,'KNmaxd',0.0032)
%%   Run_One_Subfield(10,'Kbot',NaN)            % explicit free drainage
%%   Run_One_Subfield(2,'ParallelWorkers',6)    % parallel daily ERW columns
%%   Run_One_Subfield(3,'Scenario','ssp245','LastApplicationYear',2050, ...
%%       'ResultDirectory',D)                  % rock every November to 2050
%%   Run_One_Subfield(5,'FieldShapefile',S,'GridCellSizeM',5)
%%                                             % field grid, 5 m margin
%%   Run_One_Subfield(5,'FieldShapefile',S,'GridCellSizeM',5,'FieldDEM',D)
%%                                             % ... with topography
%%   Run_One_Subfield(5,'FieldTopologyFile',F) % a saved field grid
%%
%% FIELD GRIDS (from T&C-SMEW-N v1.0)
%%   Without a field option runs 1 to 9 use the 12-column plot, exactly as
%%   before. FieldShapefile replaces it with v1.0's polygon grid on the real
%%   field: squares of GridCellSizeM metres, kept under the majority rule
%%   (FieldMinCoverage, default 'majority'; each square more than half
%%   inside, on the lattice of least misrepresented area), the width taken
%%   from the terrain (Run_Topographic_Grid_Convergence). Every sub-field
%%   layout keeps its margin as the plot does, as the field's outer band:
%%   squares more than half inside the field shrunk by MarginWidthM metres
%%   (default the plot's 5 m) are main, the rest margin
%%   (Build_Margin_Field_Topology). All nine runs simulate the same squares;
%%   MarginWidthM 0 labels every square main, v1.0's whole-field grid, for
%%   runs 1 to 3 only. FieldDEM attaches topography, which switches the 2-D
%%   frame's OPT_DEM on; CellAgronomicFile gives every cell its own soil
%%   (OPT_CELL_INPUTS). Building a grid needs the Mapping Toolbox, so a
%%   cluster campaign loads one saved by Build_Field_Grid_Topology
%%   (FieldTopologyFile) and every run uses the identical grid.
%%
%% Returns the path of the saved result file.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ResultDirectory defaults to the usual SIM/Subfield_12, so behaviour is
%% unchanged when it is not given. It exists because two runs of the same
%% run_id launched concurrently would otherwise race on one output path.
opt = struct('NDays',[],'StartDate',[],'GridNx',[],'GridNy',[], ...
    'AppeloAlpha',[],'AppeloLogKH',[],'Pcla',[],'Psan',[], ...
    'KNmaxd',[],'Kbot',[],'CECpHSensitivity',[],'ParallelWorkers',[], ...
    'SaveResults',true,'ResultDirectory','','Scenario','', ...
    'CalciteCase',[],'LastApplicationYear',[], ...
    'BiogeoOverrides',[],'DissolutionFactor',[],'SaturationExponent',[], ...
    'FieldShapefile','','GridCellSizeM',[],'FieldMinCoverage',[], ...
    'FieldFeatureIndices',[],'FieldPolygonIndex',[],'MarginWidthM',[], ...
    'FieldDEM','','DEMResolutionM',[],'FieldTopologyFile','', ...
    'CellAgronomicFile','');
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
%% A scenario run reads a different forcing file and must not share a
%% result directory with the validated campaign: Bouchet_Check_Existing_Result
%% globs TC_2D_run%02d_*.mat and errors on more than one match.
if ~isempty(opt.Scenario)
    if ~(ischar(opt.Scenario) || (isstring(opt.Scenario) && isscalar(opt.Scenario)))
        error('Run_One_Subfield:BadOptions', ...
            'Scenario must be a character vector or scalar string.');
    end
    if isempty(opt.ResultDirectory)
        error('Run_One_Subfield:BadOptions', ...
            ['A Scenario run needs an explicit ResultDirectory so its ' ...
            'results do not collide with the validated campaign.']);
    end
end
validate_optional_integer(opt.NDays,'NDays');
validate_optional_integer(opt.GridNx,'GridNx');
validate_optional_integer(opt.GridNy,'GridNy');
validate_optional_nonnegative_scalar(opt.AppeloAlpha,'AppeloAlpha');
%%% Variable-charge sensitivity of the exchanger, cmolc kg-1 per pH unit.
%%% Empty leaves whatever the parameter file declares, which defaults to the
%%% legacy constant CEC; a value here overrides it for every column.
validate_optional_nonnegative_scalar(opt.CECpHSensitivity,'CECpHSensitivity');
validate_optional_fraction(opt.Pcla,'Pcla');
validate_optional_fraction(opt.Psan,'Psan');
validate_optional_nonnegative_scalar(opt.KNmaxd,'KNmaxd');
validate_optional_kbot(opt.Kbot);
%% Uncertainty-ensemble knobs (runtime/Uncertainty_1D): a
%% struct of Biogeochemistry_Parameter literals to replace for this run, a
%% dissolution factor for the feedstock rows of the ERW workbook, and an
%% exponent on the active-layer pore saturation the ERW module receives.
%% Each is recorded in run_provenance and in the result and checkpoint
%% tags; empty leaves the shipped value and the run bit-identical.
validate_optional_positive_scalar(opt.DissolutionFactor,'DissolutionFactor');
validate_optional_positive_scalar(opt.SaturationExponent,'SaturationExponent');
%% Rock goes on every 10 November from 2016 to LastApplicationYear. The
%% default, 2019, is the field record the validation runs reproduce; a
%% projection extends it and so needs its own ResultDirectory.
validate_optional_integer(opt.LastApplicationYear,'LastApplicationYear');
last_application_year = 2019;
if ~isempty(opt.LastApplicationYear)
    if opt.LastApplicationYear < 2016
        error('Run_One_Subfield:BadOptions', ...
            'LastApplicationYear must be 2016 or later, the first field application.');
    end
    if isempty(opt.ResultDirectory)
        error('Run_One_Subfield:BadOptions', ...
            ['A LastApplicationYear run needs an explicit ResultDirectory so its ' ...
            'results do not collide with the field-record campaign.']);
    end
    last_application_year = opt.LastApplicationYear;
end
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
%% Field grids (see FIELD GRIDS above).
opt.FieldShapefile = validate_optional_file_name(opt.FieldShapefile,'FieldShapefile');
opt.FieldTopologyFile = validate_optional_file_name( ...
    opt.FieldTopologyFile,'FieldTopologyFile');
opt.CellAgronomicFile = validate_optional_file_name( ...
    opt.CellAgronomicFile,'CellAgronomicFile');
if ~(ischar(opt.FieldDEM) || (isstring(opt.FieldDEM) && isscalar(opt.FieldDEM)))
    error('Run_One_Subfield:BadOptions', ...
        'FieldDEM must be a GeoTIFF or a directory of tiles.');
end
opt.FieldDEM = char(opt.FieldDEM);
if ~isempty(opt.FieldDEM) && ~(isfile(opt.FieldDEM) || isfolder(opt.FieldDEM))
    error('Run_One_Subfield:MissingDEM','FieldDEM does not exist: %s',opt.FieldDEM);
end
validate_optional_positive_scalar(opt.GridCellSizeM,'GridCellSizeM');
validate_optional_nonnegative_scalar(opt.MarginWidthM,'MarginWidthM');
validate_optional_positive_scalar(opt.DEMResolutionM,'DEMResolutionM');
validate_optional_index_vector(opt.FieldFeatureIndices,'FieldFeatureIndices');
validate_optional_integer(opt.FieldPolygonIndex,'FieldPolygonIndex');
%% FieldMinCoverage: v1.0's coverage rule, 'majority' or a fraction of the
%% full-cell area below 1 (Build_Projected_Field_Topology). This tree
%% defaults to the majority rule.
field_coverage = 'majority';
if ~isempty(opt.FieldMinCoverage)
    if (ischar(opt.FieldMinCoverage) || isstring(opt.FieldMinCoverage)) && ...
            strcmp(char(opt.FieldMinCoverage),'majority')
        opt.FieldMinCoverage = 'majority';
    else
        validate_optional_positive_scalar(opt.FieldMinCoverage,'FieldMinCoverage');
        if opt.FieldMinCoverage >= 1
            error('Run_One_Subfield:BadOptions', ...
                ['FieldMinCoverage is a fraction of full-cell area and must be ' ...
                 'below 1, or ''majority''.']);
        end
    end
    field_coverage = opt.FieldMinCoverage;
end
has_field_shape = ~isempty(opt.FieldShapefile);
has_field_topology = ~isempty(opt.FieldTopologyFile);
uses_field_grid = has_field_shape || has_field_topology;
if has_field_shape && has_field_topology
    error('Run_One_Subfield:ConflictingFieldGrid', ...
        'Give FieldShapefile or FieldTopologyFile, not both.');
end
%% Without a shapefile these would be ignored by the grid but still name
%% the result and the checkpoint, so they are refused instead; a saved
%% grid carries its own.
if ~has_field_shape && (~isempty(opt.GridCellSizeM) || ...
        ~isempty(opt.FieldFeatureIndices) || ~isempty(opt.FieldPolygonIndex) || ...
        ~isempty(opt.FieldDEM) || ~isempty(opt.DEMResolutionM) || ...
        ~isempty(opt.FieldMinCoverage) || ~isempty(opt.MarginWidthM))
    error('Run_One_Subfield:IncompleteFieldGrid', ...
        ['GridCellSizeM, FieldMinCoverage, FieldFeatureIndices, ' ...
         'FieldPolygonIndex, MarginWidthM, FieldDEM and DEMResolutionM ' ...
         'require FieldShapefile.']);
end
if uses_field_grid && (~isempty(opt.GridNx) || ~isempty(opt.GridNy))
    error('Run_One_Subfield:ConflictingGridOptions', ...
        'GridNx and GridNy cannot be combined with a field grid.');
end
%% Configure_Cell_Agronomic_Inputs writes each cell's own sand and clay, so
%% a run-level texture would be overwritten while still naming the result.
if ~isempty(opt.CellAgronomicFile) && ~isempty(opt.Pcla)
    error('Run_One_Subfield:ConflictingTexture', ...
        ['Pcla and Psan set one texture for every column; CellAgronomicFile ' ...
         'sets texture cell by cell. Supply one or the other.']);
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
if uses_field_grid && ~R.use_subdomains
    error('Run_One_Subfield:BadFieldLayout', ...
        'A field grid is for the sub-field runs 1 to 9; run %d is one column.', ...
        run_id);
end
%% Every sub-field layout keeps the plot's margin by default, the field's
%% outer band, as the 12-column plot keeps its outer 5 m ring for all nine
%% runs; a margin layout cannot do without one.
field_margin_m = R.geom.margin_width_m;
if ~isempty(opt.MarginWidthM); field_margin_m = opt.MarginWidthM; end
if has_field_shape && field_margin_m == 0 && ~strcmp(R.layout,'whole_field')
    error('Run_One_Subfield:MarginRequired', ...
        'The %s layout needs a margin: MarginWidthM must be positive.',R.layout);
end

%% The calcite application rate lives in the workbook as M_rock_event, so a
%% second rate is a second ERW row rather than a fourth treatment: the run
%% matrix, its ids, and every test that pins them stay as they are. Only a
%% calcite run has a rate to vary, so asking for one on any other run would
%% duplicate a result that already exists under a misleading tag.
if ~isempty(opt.CalciteCase)
    validateattributes(opt.CalciteCase,{'numeric'}, ...
        {'scalar','real','finite','integer','positive'});
    if ~strcmp(R.treatment,'calcite')
        error('Run_One_Subfield:CalciteCaseOnNonCalciteRun', ...
            'CalciteCase applies only to a calcite run; run %d is %s.', ...
            run_id,R.treatment);
    end
    if isempty(opt.ResultDirectory)
        error('Run_One_Subfield:BadOptions', ...
            ['A CalciteCase run needs an explicit ResultDirectory so its ' ...
            'results do not collide with the default-rate campaign.']);
    end
    R.case_number = opt.CalciteCase;
end

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
ERW_opts.DataRange = "A2:AJ7";
ERW_opts.VariableTypes = {'double','string','string','double','double',...
    'string','double','double','string','double',...
    'string','string','string','string','double', ...
    'double','double','double','double','double',...
    'double','double','double','string','double',...
    'string','string','double','double','double',...
    'double','double','double','double','double',...
    'double'};
ERW_opts = setvaropts(ERW_opts,{'f_CEC_in','rock_filepath','mineral','rock_f_in'},'FillValue','NaN');
ERW_PARAM_input = readtable(wb,ERW_opts);

[MOD_PARAM_input,ERW_PARAM_input,ERW_PARAM_control] = ...
    Select_Subfield_Parameter_Rows( ...
    MOD_PARAM_input,ERW_PARAM_input,R.case_number);
MOD_PARAM_input.Zbio(:) = R.biogeo_depth_mm;
if ~isempty(opt.Kbot)
    MOD_PARAM_input.Kbot(:) = opt.Kbot;
end
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
if ~isempty(opt.DissolutionFactor)
    %% Only feedstock rows carry a dissolution factor; the control row's
    %% NaN (no rock) stays NaN, so a control run records the request and
    %% weathers nothing either way.
    feedstock_rows = isfinite(ERW_PARAM_input.diss_f);
    ERW_PARAM_input.diss_f(feedstock_rows) = opt.DissolutionFactor;
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
if ~isempty(opt.BiogeoOverrides)
    %% The parameter script owns the rules (which literals may move, their
    %% domain); asking it now fails before any forcing is loaded.
    Biogeochemistry_Parameter(1,[],opt.BiogeoOverrides);
    override_names = fieldnames(opt.BiogeoOverrides);
    for k = 1:numel(override_names)
        fprintf('BiogeoPar.%s: %.8g (shipped %.8g)\n',override_names{k}, ...
            opt.BiogeoOverrides.(override_names{k}), ...
            default_biogeo.(override_names{k}));
    end
end
saturation_exponent = 1;
if ~isempty(opt.SaturationExponent)
    saturation_exponent = opt.SaturationExponent;
    fprintf('ERW pore-saturation exponent: %.8g\n',saturation_exponent);
end
if ~isempty(opt.DissolutionFactor)
    fprintf('Dissolution factor: %.8g\n',opt.DissolutionFactor);
end
effective_Kbot = double(MOD_PARAM_input.Kbot);
if isnan(effective_Kbot)
    fprintf('Kbot: free drainage\n');
else
    fprintf('Kbot: %.8g mm h^-1\n',effective_Kbot);
end

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

%% Only the meteorological forcing is scenario-specific. The
%% biogeochemical seed loaded above holds pools, moisture and N
%% availability at the start date, none of which depend on which
%% future follows, so it stays shared across scenarios.
forcing_id = id_location;
if ~isempty(opt.Scenario)
    forcing_id = [id_location '_' char(opt.Scenario)];
end
load(fullfile(Input_Directory,['Data_' forcing_id '.mat']));     %#ok<LOAD>

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

%% A scenario forcing file carries its own Ca on the same hourly axis as
%% Date, because the SSP CO2 pathway must match the climate deltas it was
%% built with. Detect it before Date is trimmed.
forcing_Ca = exist('Ca','var')==1 && numel(Ca)==numel(Date);

Date=Date(x1:x2); Pr=Pr(x1:x2); Ta=Ta(x1:x2); Ws=Ws(x1:x2); ea=ea(x1:x2);
SAD1=SAD1(x1:x2); SAD2=SAD2(x1:x2); SAB1=SAB1(x1:x2); SAB2=SAB2(x1:x2);
Pre=Pre(x1:x2); N=N(x1:x2); Tdew=Tdew(x1:x2); esat=esat(x1:x2);
PARB=PARB(x1:x2); PARD=PARD(x1:x2);
if forcing_Ca; Ca = Ca(x1:x2); end

%% Kantola's first amendment event is 10-Nov-2016. A shortened smoke
%% window ending before that date contains no rock application, so those
%% cells must initialise from the control row. The full run retains the
%% treatment row because its future application events are in the window.
treatment_in_window = strcmp(R.treatment,'control') || ...
    Date(end) >= datenum(2016,11,10);

t_bef = -0.67; t_aft = 1.67;
Ds = esat-ea; Ds(Ds<0) = 0;
%% Ca_Data.mat is the shared historical series and ends in 2023, so it
%% cannot serve a projection run. Use it only when the forcing carries no Ca
%% of its own, and fail loudly rather than returning a short or empty series.
if ~forcing_Ca
    Ca_file = load(fullfile(Input_Directory,'Ca_Data.mat'));
    d1 = find(abs(Ca_file.Date_CO2-Date(1))<1/36);
    d2 = find(abs(Ca_file.Date_CO2-Date(end))<1/36);
    if isempty(d1) || isempty(d2)
        lo = datetime(Date(1),'ConvertFrom','datenum','Format','yyyy-MM-dd');
        hi = datetime(Date(end),'ConvertFrom','datenum','Format','yyyy-MM-dd');
        error('Run_One_Subfield:CO2WindowNotCovered', ...
            ['Ca_Data.mat does not cover %s to %s and the forcing carries ' ...
            'no Ca of its own. Use a forcing file with an embedded SSP CO2 ' ...
            'series, or extend Ca_Data.mat.'],char(lo),char(hi));
    end
    Ca = Ca_file.Ca(d1:d2);
end
if numel(Ca) ~= numel(Date)
    error('Run_One_Subfield:CO2LengthMismatch', ...
        'CO2 series has %d values for %d forcing hours.', ...
        numel(Ca),numel(Date));
end
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
if has_field_shape
    topo = Build_Field_Grid_Topology(opt.FieldShapefile, ...
        'GridCellSizeM',opt.GridCellSizeM,'CoverageRule',field_coverage, ...
        'MarginWidthM',field_margin_m, ...
        'FeatureIndices',opt.FieldFeatureIndices, ...
        'PolygonIndex',opt.FieldPolygonIndex, ...
        'FieldDEM',opt.FieldDEM,'DEMResolutionM',opt.DEMResolutionM);
elseif has_field_topology
    topo = load_field_grid_topology(opt.FieldTopologyFile);
elseif R.use_subdomains
    if ~isempty(opt.GridNx); nx = opt.GridNx; end
    if ~isempty(opt.GridNy); ny = opt.GridNy; end
    topo = Build_Cell_Topology(g.field_side_m,g.margin_width_m,nx,ny);
else
    topo = Build_Single_Cell_Topology(g.field_side_m^2);
end
if uses_field_grid
    nx = topo.main_grid_nx; ny = topo.main_grid_ny;
    if ~strcmp(R.layout,'whole_field') && all(strcmp(topo.cell_region,'main'))
        error('Run_One_Subfield:MarginRequired', ...
            'The %s layout needs a field grid with a margin band.',R.layout);
    end
end
%% The 2-D frame's two field switches follow what this run was given:
%% topography exactly when the grid carries a DEM, per-cell soils exactly
%% when a cell table is supplied. Without either the frame is the flat,
%% uniform frame bit for bit.
OPT_DEM = double(isfield(topo,'topography_attached') && topo.topography_attached);
OPT_CELL_INPUTS = double(~isempty(opt.CellAgronomicFile));
parallel_workers = resolve_parallel_workers(opt.ParallelWorkers,topo.N_cells);
run_provenance.matlab_version = version;
run_provenance.matlab_release = version('-release');
run_provenance.computer = computer;
run_provenance.blas = regexprep(strtrim(evalc('version -blas')),'\s+',' ');
run_provenance.max_threads = maxNumCompThreads;
run_provenance.parallel_workers = parallel_workers;
run_provenance.numerical_identity = getenv('TC_NUMERICAL_IDENTITY');
%% The ERW chemistry kernels the steps run (TC_ERW_KERNEL): matlab, the
%% original code, or fast or mex, which compute the same numbers bit for
%% bit. Checkpoints do not pin it, so the segments before a resume keep
%% their own in execution_segments. MATLAB runs no MEX function on a
%% thread-based worker, so mex with a pool stops here, not at the first day.
run_provenance.erw_kernel = ERW_Chemistry_Kernel_Mode([]);
if strcmp(run_provenance.erw_kernel,'mex') && parallel_workers > 0
    error('Run_One_Subfield:MexOnThreadPool', ...
        ['TC_ERW_KERNEL=mex cannot run on the %d thread-pool workers of ' ...
         'this run; TC_ERW_KERNEL=fast computes the same numbers there.'], ...
        parallel_workers);
end
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
run_provenance.sweep_layout = getenv('TC_SWEEP_LAYOUT');
sweep_case_id = getenv('TC_SWEEP_CASE_ID');
if isempty(sweep_case_id)
    sweep_case_id = run_provenance.slurm_array_task_id;
end
run_provenance.sweep_case_id = sweep_case_id;
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
run_provenance.Kbot = effective_Kbot;
run_provenance.last_application_year = last_application_year;
run_provenance.biogeo_overrides = opt.BiogeoOverrides;
run_provenance.dissolution_factor = opt.DissolutionFactor;
run_provenance.saturation_exponent = saturation_exponent;
%% Recorded only on a field grid or with per-cell soils, so a plot run's
%% provenance keeps its field set.
if uses_field_grid
    run_provenance.field_grid_build = topo.field_grid_build;
    run_provenance.field_topology_file = opt.FieldTopologyFile;
    if has_field_topology
        run_provenance.field_topology_file_sha256 = file_sha256(opt.FieldTopologyFile);
    end
    run_provenance.opt_dem = OPT_DEM;
end
if OPT_CELL_INPUTS == 1
    run_provenance.opt_cell_inputs = OPT_CELL_INPUTS;
    run_provenance.cell_agronomic_file = opt.CellAgronomicFile;
    run_provenance.cell_agronomic_file_sha256 = file_sha256(opt.CellAgronomicFile);
end

base_param = fullfile(Parameter_Directory,[char(MOD_PARAM_input.MOD_PARAM) '.m']);
is_main = strcmp(topo.cell_region,'main');
ColCfg = repmat(struct('PARAM_IC',base_param,'PFT_OVERRIDE','', ...
                       'apply_rock',false,'ERW_PARAM',ERW_PARAM_control, ...
                       'appelo_alpha',opt.AppeloAlpha, ...
                       'appelo_logKH',opt.AppeloLogKH, ...
                       'cec_dpH',opt.CECpHSensitivity, ...
                       'KNmaxd',opt.KNmaxd, ...
                       'last_application_year',last_application_year, ...
                       'biogeo_overrides',opt.BiogeoOverrides, ...
                       'saturation_exponent',saturation_exponent), ...
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
%% Per-cell soil (OPT_CELL_INPUTS). Resolve validates the table against
%% the grid and closes the texture budget; Configure clones each column's
%% MOD and ERW rows and overwrites PH, Psan, Pcla, Porg and CEC_tot_in.
%% Unlike v1.0 no per-cell spin-up state is handed over: every cell starts
%% from the same spin-up pools, with its own soil. The new ColCfg fields
%% exist only on this path, so ColCfg keeps its shape on every other run.
if OPT_CELL_INPUTS == 1
    CellAgronomicInputs = Resolve_Cell_Agronomic_Inputs(topo, ...
        readtable(opt.CellAgronomicFile,'VariableNamingRule','preserve'));
    ERW_cell_rows = cell(topo.N_cells,1);
    for c = 1:topo.N_cells
        ERW_cell_rows{c} = ColCfg(c).ERW_PARAM;
    end
    CellInputs = Configure_Cell_Agronomic_Inputs(MOD_PARAM_input, ...
        vertcat(ERW_cell_rows{:}),CellAgronomicInputs);
    for c = 1:topo.N_cells
        ColCfg(c).MOD_PARAM = CellInputs(c).MOD_PARAM;
        ColCfg(c).ERW_PARAM = CellInputs(c).ERW_PARAM;
        ColCfg(c).agronomic_identity = CellInputs(c).agronomic_identity;
    end
    fprintf('Cell soil inputs: %d columns from %s\n', ...
        topo.N_cells,opt.CellAgronomicFile);
end
if uses_field_grid
    field_area = topo.cell_overlap_area_m2(:);
    field_missed_m2 = 0;
    if isfield(topo,'field_missed_area_m2')
        field_missed_m2 = topo.field_missed_area_m2;
    end
    fprintf(['Domain: %d columns of %.3g m on a field grid (%.2f deg): ' ...
        '%.0f m2 main + %.0f m2 margin (%g m wide)\n'], ...
        topo.N_cells,topo.cell_size_m,topo.grid_rotation_degrees, ...
        sum(field_area(is_main)),sum(field_area(~is_main)),topo.buffer_width_m);
    fprintf(['Grid: %.0f m2 simulated, %.0f m2 of it outside the field; ' ...
        '%.0f m2 of the field in dropped squares\n'], ...
        sum(topo.cell_area_m2),sum(topo.cell_area_m2)-sum(field_area), ...
        field_missed_m2);
    if OPT_DEM == 1
        fprintf(['Topography: relief %.2f m, slope %.1f-%.1f%% (mean ' ...
            '%.1f%%), %d outlet cell(s)\n'], ...
            max(topo.cell_z_m)-min(topo.cell_z_m),100*min(topo.cell_slope), ...
            100*max(topo.cell_slope),100*mean(topo.cell_slope), ...
            sum(topo.cell_is_outlet));
    else
        fprintf('Topography: none (flat field grid)\n');
    end
    fprintf('Rock on: %d main cell(s), %d margin cell(s)\n', ...
        sum([ColCfg(is_main).apply_rock]),sum([ColCfg(~is_main).apply_rock]));
elseif R.use_subdomains
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
if ~isempty(opt.Scenario)
    ckpt_tag = sprintf('%s__%s',ckpt_tag,char(opt.Scenario));
end
if ~isempty(opt.CalciteCase)
    ckpt_tag = sprintf('%s__case%d',ckpt_tag,opt.CalciteCase);
end
%% A field-grid run, and a topographic or per-cell-soil one, never shares a
%% checkpoint or a result file with a plot run or with each other.
if uses_field_grid
    ckpt_tag = sprintf('%s__%s',ckpt_tag,field_run_tag(topo));
end
if OPT_CELL_INPUTS == 1
    ckpt_tag = sprintf('%s__%s',ckpt_tag,cell_soil_tag(opt.CellAgronomicFile));
end
if ~isempty(opt.CECpHSensitivity)
    ckpt_tag = sprintf('%s__dcec%s',ckpt_tag,number_tag(opt.CECpHSensitivity));
end
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
if ~isempty(opt.BiogeoOverrides)
    ckpt_tag = sprintf('%s__bg%s',ckpt_tag,biogeo_override_tag(opt.BiogeoOverrides));
end
if ~isempty(opt.DissolutionFactor)
    ckpt_tag = sprintf('%s__fd%s',ckpt_tag,number_tag(opt.DissolutionFactor));
end
if ~isempty(opt.SaturationExponent)
    ckpt_tag = sprintf('%s__sexp%s',ckpt_tag,number_tag(opt.SaturationExponent));
end
if ~isempty(opt.Kbot)
    ckpt_tag = sprintf('%s__kbot%s',ckpt_tag,kbot_tag(opt.Kbot));
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
cec_dpH_effective = zeros(topo.N_cells,1);
cec_ref_pH_effective = nan(topo.N_cells,1);
for c = 1:topo.N_cells
    cec_params = ERW{c}.const.bulk_pH_params;
    cec_scheme = 'TC_GT';
    if isfield(cec_params,'scheme') && ~isempty(cec_params.scheme)
        cec_scheme = char(string(cec_params.scheme));
    end
    [cec_dpH_effective(c),cec_ref_pH_effective(c)] = ...
        ERW_Resolve_CEC_pH_Sensitivity(cec_params,cec_scheme);
end
run_provenance.CECpHSensitivity = cec_dpH_effective;
run_provenance.CECrefpH = cec_ref_pH_effective;
run_provenance.CEC_cmolc_kg = cellfun(@(x) ...
    double(x.const.bulk_pH_params.CEC_cmolc_kg),ERW(:));
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
    if ~isempty(opt.Scenario)
        tag = sprintf('%s__%s',tag,char(opt.Scenario));
    end
    if ~isempty(opt.CalciteCase)
        tag = sprintf('%s__case%d',tag,opt.CalciteCase);
    end
    if uses_field_grid
        tag = sprintf('%s__%s',tag,field_run_tag(topo));
    elseif R.use_subdomains && (~isempty(opt.GridNx) || ~isempty(opt.GridNy))
        tag = sprintf('%s__grid%dx%d',tag,nx,ny);
    end
    if OPT_CELL_INPUTS == 1
        tag = sprintf('%s__%s',tag,cell_soil_tag(opt.CellAgronomicFile));
    end
    if ~isempty(opt.CECpHSensitivity)
        tag = sprintf('%s__dcec%s',tag,number_tag(opt.CECpHSensitivity));
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
    if ~isempty(opt.BiogeoOverrides)
        tag = sprintf('%s__bg%s',tag,biogeo_override_tag(opt.BiogeoOverrides));
    end
    if ~isempty(opt.DissolutionFactor)
        tag = sprintf('%s__fd%s',tag,number_tag(opt.DissolutionFactor));
    end
    if ~isempty(opt.SaturationExponent)
        tag = sprintf('%s__sexp%s',tag,number_tag(opt.SaturationExponent));
    end
    if ~isempty(opt.Kbot)
        tag = sprintf('%s__kbot%s',tag,kbot_tag(opt.Kbot));
    end
    if ~isempty(opt.NDays); tag = sprintf('%s__%dd',tag,opt.NDays); end
    if ~isempty(start_datenum); tag = sprintf('%s__from%s',tag,datestr(start_datenum,'yyyymmdd')); end
    result_file = fullfile(fl_save,sprintf('TC_2D_run%02d_%s.mat',R.id,tag));
    tmp_result = [tempname(fl_save) '.mat'];
    try
        if FULL_RECORDS
            save(tmp_result,'-v7.3');
        else
            %% A daily-store result (package y): the stores, the hourly field
            %% means, the final state and the run's inputs and metadata. The
            %% hourly windows and the daily arrays outside the stores are in
            %% Final_state, in order; under their own names a reader would
            %% take a window for the whole run.
            result_vars = setdiff(who,[{'cleanup'} HOURLY_NAMES DAILY_STATE_ONLY]);
            save(tmp_result,result_vars{:},'-v7.3');
        end
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

function validate_optional_kbot(value)
if isempty(value); return; end
if ~(isnumeric(value) && isreal(value) && isscalar(value) && ...
        (isnan(value) || (isfinite(value) && value >= 0)))
    error('Run_One_Subfield:BadOptions', ...
        'Kbot must be NaN for free drainage or a finite nonnegative scalar.');
end
end

function validate_optional_positive_scalar(value,name)
if isempty(value); return; end
if ~(isnumeric(value) && isreal(value) && isscalar(value) && ...
        isfinite(value) && value > 0)
    error('Run_One_Subfield:BadOptions', ...
        '%s must be a finite positive scalar.',name);
end
end

function tag = biogeo_override_tag(overrides)
%% One or two replaced literals spell themselves; more take a count and
%% the head of a digest of the sorted name=value list, which
%% run_provenance.biogeo_overrides carries in full.
names = sort(fieldnames(overrides));
parts = cellfun(@(n) sprintf('%s%s',n,number_tag(overrides.(n))), ...
    names,'UniformOutput',false);
if numel(parts) <= 2
    tag = strjoin(parts,'_');
    return
end
digest = javaMethod('getInstance','java.security.MessageDigest','SHA-256');
digest.update(typecast(uint8(strjoin(parts,',')),'int8'));
raw = typecast(digest.digest(),'uint8');
hex = lower(reshape(dec2hex(raw,2).',1,[]));
tag = sprintf('%dfields_%s',numel(parts),hex(1:8));
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
text = sprintf('%.12g',value);
if str2double(text) ~= value
    text = sprintf('%s_h%s',text,num2hex(value));
end
tag = strrep(text,'.','p');
tag = strrep(tag,'-','m');
tag = strrep(tag,'+','');
end

function tag = kbot_tag(value)
if isnan(value)
    tag = 'free';
else
    tag = number_tag(value);
end
end

function file_name = validate_optional_file_name(value,name)
if ~(ischar(value) || (isstring(value) && isscalar(value)))
    error('Run_One_Subfield:BadOptions', ...
        '%s must be a character vector or scalar string.',name);
end
file_name = char(value);
if isempty(file_name)
    return
end
if ~isfile(file_name)
    error('Run_One_Subfield:MissingInputFile', ...
        '%s does not exist: %s',name,file_name);
end
[~,attributes] = fileattrib(file_name);
file_name = attributes.Name;
end

function validate_optional_index_vector(value,name)
if isempty(value); return; end
if ~(isnumeric(value) && isvector(value) && all(isfinite(value)) && ...
        all(value == fix(value)) && all(value >= 1) && ...
        numel(unique(value)) == numel(value))
    error('Run_One_Subfield:BadOptions', ...
        '%s must contain unique positive integers.',name);
end
end

function topo = load_field_grid_topology(file_name)
%% A grid saved by Build_Field_Grid_Topology. Its build record names the
%% result exactly as a run built from the same inputs would be named.
S = load(file_name,'topo');
if ~isfield(S,'topo') || ~isstruct(S.topo) || ~isscalar(S.topo)
    error('Run_One_Subfield:BadFieldTopologyFile', ...
        '%s holds no topology struct named topo.',file_name);
end
topo = S.topo;
required = {'N_cells','N_edges','edge_a','edge_b','cell_area_m2', ...
    'cell_overlap_area_m2','cell_region','cell_size_m','main_grid_nx', ...
    'main_grid_ny','buffer_width_m','grid_rotation_degrees','field_grid_build'};
missing = required(~isfield(topo,required));
if ~isempty(missing)
    error('Run_One_Subfield:BadFieldTopologyFile', ...
        '%s is not a field grid (missing %s); build it with Build_Field_Grid_Topology.', ...
        file_name,strjoin(missing,', '));
end
if topo.N_cells < 2
    error('Run_One_Subfield:BadFieldTopologyFile', ...
        '%s holds a one-cell grid; a sub-field run needs a lateral graph.',file_name);
end
end

function tag = field_run_tag(topo)
%% v1.0's field grid tag (shapefile, cell width, polygon, features), with
%% the margin and the DEM spacing; the polygon only when the shapefile
%% holds more than one field, and the coverage rule only when it is not the
%% majority rule. The shapefile name is cut to 24 characters, since a
%% result name must stay within 255.
build = topo.field_grid_build;
[~,shape_name] = fileparts(build.shapefile);
shape_name = regexprep(shape_name,'[^A-Za-z0-9]+','_');
shape_name = shape_name(1:min(end,24));
tag = sprintf('field_%s__cell%sm',shape_name,number_tag(build.cell_size_m));
if isfield(topo,'setup_count') && topo.setup_count > 1
    tag = sprintf('%s__polygon%02d',tag,build.polygon_index);
end
if ~isempty(build.feature_indices)
    parts = arrayfun(@(x) sprintf('%d',x),build.feature_indices, ...
        'UniformOutput',false);
    tag = sprintf('%s__features%s',tag,strjoin(parts,'_'));
end
if build.margin_width_m > 0
    tag = sprintf('%s__margin%s',tag,number_tag(build.margin_width_m));
end
if ~isempty(build.dem_files)
    tag = sprintf('%s__dem%s',tag,number_tag(build.dem_resolution_m));
end
if ~(ischar(build.coverage_rule) && strcmp(build.coverage_rule,'majority'))
    tag = sprintf('%s__cov%s',tag,strrep(sprintf('%g',build.coverage_rule),'.','p'));
end
end

function tag = cell_soil_tag(file_name)
hex = file_sha256(file_name);
tag = ['cellsoil_' hex(1:12)];
end

function hex = file_sha256(file_name)
fid = fopen(file_name,'rb');
if fid < 0
    error('Run_One_Subfield:FileRead','Cannot read %s',file_name);
end
md = javaMethod('getInstance','java.security.MessageDigest','SHA-256');
while true
    bytes = fread(fid,64*2^20,'*uint8');
    if isempty(bytes); break; end
    md.update(typecast(bytes,'int8'));
end
fclose(fid);
raw = typecast(md.digest(),'uint8');
hex = lower(reshape(dec2hex(raw,2).',1,[]));
end
