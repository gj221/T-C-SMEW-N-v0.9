function C = Field_Uncertainty_Campaign
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The field-grid uncertainty campaign an array element belongs to: the
%% lower and upper sets of the Kantola ensemble (Uncertainty_Members) on
%% runs of the central campaign on a field grid, read from the
%% controller's environment:
%%   TC_SSP             unset: the five-year field record on the legacy
%%                      forcing; ssp126 | ssp245 | ssp585: that scenario's
%%                      35-year projection (Field_Campaign_Grid 'horizon')
%%   TC_HORIZON         5yr (or unset) without TC_SSP, 35yr (or unset) with it
%%   TC_CALCITE_CASE    unset, or an ERW workbook row for the calcite runs
%%                      only (8, the agronomic 5 t/ha rate)
%%   TC_FIELD_TOPOLOGY  the grid file beside the deployed root
%%   TC_FIELD_RUNS      the runs, separated by commas (the controller
%%                      writes them so); unset is every run of the arm
%%   TC_FIELD_WORKERS   the pool
%% The original set is the central field campaign's own result
%% (runtime/Field_SSP_2D, tag original_tag), never run here. Each set runs
%% as that campaign runs: the horizon's days, forcing and rock schedule,
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

H = Field_Campaign_Grid('horizon',prefix);
C.dimension = '2D';
C.scenario = H.scenario;
C.horizon = H.horizon;
C.n_days = H.n_days;
C.last_application_year = H.last_application_year;
C.forcing = H.forcing;
[C.calcite_case,C.case_suffix] = Calcite_Case;
[C.topology_file,C.topology_name] = Field_Campaign_Grid('read',repo_root,prefix);
C.topo = Field_Campaign_Grid('load',C.topology_file,prefix);
C.n_columns = C.topo.N_cells;
C.tag = sprintf('Field_Uncertainty_2D_%s%s__%s',H.tag,C.case_suffix,C.topology_name);
C.original_tag = sprintf('Field_SSP_2D_%s%s__%s',H.tag,C.case_suffix,C.topology_name);

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
C.options = [{'FieldTopologyFile',C.topology_file} H.options];
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
    fullfile(repo_root,'v2','T&C_Code','Inputs',C.forcing), ...
    C.topology_file};

members = Uncertainty_Members;
members = members(~[members.reused]);
[r,m] = ndgrid(1:numel(C.runs),1:numel(members));
C.elements = struct('member',num2cell(reshape(members(m),1,[])), ...
    'run_id',num2cell(reshape(C.runs(r),1,[])));
end
