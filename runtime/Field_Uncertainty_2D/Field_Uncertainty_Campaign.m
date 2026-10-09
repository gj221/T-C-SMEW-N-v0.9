function C = Field_Uncertainty_Campaign
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The field-grid uncertainty campaign an array element belongs to: the
%% lower and upper sets of the Kantola ensemble (Uncertainty_Members) on
%% runs of the 35-year SSP campaign on a field grid, read from the
%% controller's environment:
%%   TC_SSP             ssp126 | ssp245 | ssp585
%%   TC_CALCITE_CASE    unset, or an ERW workbook row for the calcite runs
%%                      only (8, the agronomic 5 t/ha rate)
%%   TC_FIELD_TOPOLOGY  the grid file beside the deployed root
%%   TC_FIELD_RUNS      the runs, separated by commas (the controller
%%                      writes them so); unset is every run of the arm
%%   TC_FIELD_WORKERS   the pool
%% The original set is the field SSP campaign's own result
%% (runtime/Field_SSP_2D, tag original_tag), never run here. Each set runs
%% as that campaign runs: 12784 days, rock to 2050, the scenario forcing,
%% the field grid in place of the 2x2 plot grid. Elements are member-major
%% over the runs, as Uncertainty_Campaign orders them.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

prefix = 'Kantola_FieldUncertainty';
here = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(here));
calibrated = fullfile(repo_root,'runtime','Base_1D');
uncertainty_1d = fullfile(repo_root,'runtime','Uncertainty_1D');
field_ssp = fullfile(repo_root,'runtime','Field_SSP_2D');
addpath(calibrated,uncertainty_1d,field_ssp);

ssp = strtrim(getenv('TC_SSP'));
if ~ismember(ssp,{'ssp126','ssp245','ssp585'})
    error([prefix ':BadScenario'], ...
        'TC_SSP must be ssp126, ssp245 or ssp585; got "%s".',ssp);
end
C.dimension = '2D';
C.scenario = ssp;
C.horizon = '35yr';
C.n_days = 12784;                 % 2016-2050, as the SSP campaigns
C.last_application_year = 2050;
[C.calcite_case,C.case_suffix] = Calcite_Case;
[C.topology_file,C.topology_name] = Field_Campaign_Grid('read',repo_root,prefix);
C.topo = Field_Campaign_Grid('load',C.topology_file,prefix);
C.n_columns = C.topo.N_cells;
C.tag = sprintf('Field_Uncertainty_2D_35yr_%s%s__%s',ssp,C.case_suffix,C.topology_name);
C.original_tag = sprintf('Field_SSP_2D_35yr_%s%s__%s',ssp,C.case_suffix,C.topology_name);

runs_text = strtrim(getenv('TC_FIELD_RUNS'));
if isempty(runs_text)
    C.runs = 1:9;
    if ~isempty(C.calcite_case); C.runs = [3 6 9]; end
else
    C.runs = Field_Campaign_Grid('runs',runs_text,C.calcite_case,prefix);
end
raw_workers = strtrim(getenv('TC_FIELD_WORKERS'));
field_workers = str2double(raw_workers);
if isempty(raw_workers) || ~(isfinite(field_workers) && field_workers == fix(field_workers) && ...
        field_workers >= 2 && field_workers <= 64)
    error([prefix ':BadWorkers'],'TC_FIELD_WORKERS must be a whole number from 2 to 64.');
end
%% The pool of the original's campaign; the count moves no number but is
%% part of the declared identity.
C.workers = min(field_workers,C.n_columns);
C.options = {'FieldTopologyFile',C.topology_file,'Scenario',ssp, ...
    'LastApplicationYear',C.last_application_year};
if ~isempty(C.calcite_case)
    C.options = [C.options {'CalciteCase',C.calcite_case}];
end

%% Files outside the task's own directory that decide what it runs: the
%% member, knob and result rules, the requeue helper, the grid helper and
%% the conservation check, the calibrated set, the forcing and the grid.
C.dependencies = {fullfile(uncertainty_1d,'*.m'),fullfile(uncertainty_1d,'*.sh'), ...
    fullfile(field_ssp,'*.m'), ...
    fullfile(calibrated,'Campaign_Parameters.m'), ...
    fullfile(calibrated,'Calcite_Case.m'), ...
    fullfile(calibrated,'Fingerprint_Digest.m'), ...
    fullfile(repo_root,'v2','T&C_Code','Inputs',sprintf('Data_Kantola_%s.mat',ssp)), ...
    C.topology_file};

members = Uncertainty_Members;
members = members(~[members.reused]);
[r,m] = ndgrid(1:numel(C.runs),1:numel(members));
C.elements = struct('member',num2cell(reshape(members(m),1,[])), ...
    'run_id',num2cell(reshape(C.runs(r),1,[])));
end
