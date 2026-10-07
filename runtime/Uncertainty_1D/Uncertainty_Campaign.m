function C = Uncertainty_Campaign(dimension)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The uncertainty campaign an array element belongs to, read from the
%% controller's environment as the Base and SSP campaigns read theirs:
%%   TC_SSP           unset: the five-year validation on the legacy forcing;
%%                    ssp126|ssp245|ssp585: that scenario's projection
%%   TC_HORIZON       5yr (or unset) without TC_SSP, 35yr with it
%%   TC_CALCITE_CASE  unset, or an ERW workbook row for the calcite runs
%%                    only (8, the agronomic 5 t/ha rate)
%% dimension is '1D' (runs 10-12, one column, serial) or '2D' (runs 1-9 on
%% the 2x2 grid, twelve columns). The five-year originals are the Base
%% campaigns and the 35-year ones the SSP projections, so each horizon
%% takes its run length, application years and worker count from the
%% campaign it is compared with. Every run of the arm is run with each set
%% the ensemble runs (member-major), under one tag per dimension, horizon
%% and case.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

here = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(here));
calibrated = fullfile(repo_root,'runtime','Base_1D');
addpath(calibrated);

if ~ismember(dimension,{'1D','2D'})
    error('Kantola_Uncertainty:BadDimension','dimension must be 1D or 2D.');
end
C.dimension = dimension;
ssp = strtrim(getenv('TC_SSP'));
horizon = strtrim(getenv('TC_HORIZON'));
if isempty(ssp)
    if ~isempty(horizon) && ~strcmp(horizon,'5yr')
        error('Kantola_Uncertainty:BadHorizon', ...
            'Without TC_SSP the ensemble is the 5-year validation; TC_HORIZON=%s is not it.',horizon);
    end
    horizon = '5yr';
    C.n_days = 1826;                  % 2016-2020, as Campaign_Parameters
    C.last_application_year = 2019;   % the field record, the driver's default
    scenario_tag = '';
else
    if ~ismember(ssp,{'ssp126','ssp245','ssp585'})
        error('Kantola_Uncertainty:BadScenario', ...
            'TC_SSP must be ssp126, ssp245 or ssp585; got "%s".',ssp);
    end
    if ~strcmp(horizon,'35yr')
        error('Kantola_Uncertainty:BadHorizon', ...
            'A scenario ensemble is the 35-year projection; TC_HORIZON=%s is not it.',horizon);
    end
    C.n_days = 12784;                 % 2016-2050, as the SSP campaigns
    C.last_application_year = 2050;
    scenario_tag = ['_35yr_' ssp];
end
C.scenario = ssp;
C.horizon = horizon;
[C.calcite_case,C.case_suffix] = Calcite_Case;

C.tag = ['Uncertainty_' dimension scenario_tag C.case_suffix];
if isempty(ssp)
    C.original_tag = ['Base_' dimension C.case_suffix];
else
    C.original_tag = ['SSP_' dimension scenario_tag C.case_suffix];
end
C.options = {};
if strcmp(dimension,'1D')
    C.runs = [10 11 12];
    if ~isempty(C.calcite_case); C.runs = 12; end
    C.n_columns = 1;
    C.workers = 0;
else
    C.runs = 1:9;
    if ~isempty(C.calcite_case); C.runs = [3 6 9]; end
    C.n_columns = 12;
    C.options = {'GridNx',2,'GridNy',2};
    %% The worker count of the original's campaign: it is part of the
    %% declared numerical identity, though it moves no number
    %% (test_2d_parfor_equivalence).
    C.workers = 6;
    if ~isempty(ssp); C.workers = 12; end
end
if ~isempty(ssp)
    C.options = [C.options {'Scenario',ssp,'LastApplicationYear',C.last_application_year}];
end
if ~isempty(C.calcite_case)
    C.options = [C.options {'CalciteCase',C.calcite_case}];
end

%% Files outside the task's own directory that decide what it runs.
C.dependencies = {fullfile(calibrated,'Campaign_Parameters.m'), ...
    fullfile(calibrated,'Calcite_Case.m')};
if strcmp(dimension,'2D')
    C.dependencies = [{fullfile(here,'*.m'),fullfile(here,'*.sh'), ...
        fullfile(repo_root,'runtime','Base_2D','Check_Lateral_Conservation.m')} ...
        C.dependencies];
end
if ~isempty(ssp)
    C.dependencies{end+1} = fullfile(repo_root,'v2','T&C_Code','Inputs', ...
        sprintf('Data_Kantola_%s.mat',ssp));
end

%% TC_UC_WITH_ORIGINAL=1 runs the original set here as well, as its own
%% campaign runs it: for when that campaign's results are far off (the
%% week partition's per-user limits hold the 35-year SSP 2-D runs to a few
%% at a time), and as a check that the reused results are this model's.
with_original = strtrim(getenv('TC_UC_WITH_ORIGINAL'));
if ~ismember(with_original,{'','0','1'})
    error('Kantola_Uncertainty:BadWithOriginal', ...
        'TC_UC_WITH_ORIGINAL must be unset, 0 or 1; got "%s".',with_original);
end
members = Uncertainty_Members;
if ~strcmp(with_original,'1')
    members = members(~[members.reused]);
end
[r,m] = ndgrid(1:numel(C.runs),1:numel(members));
C.elements = struct('member',num2cell(reshape(members(m),1,[])), ...
    'run_id',num2cell(reshape(C.runs(r),1,[])));
end
