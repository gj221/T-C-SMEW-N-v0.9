%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot the vertical 1-D group and the full twelve-scenario set in one
%% run.  Standalone runnable script: no function signature, no arguments.
%%
%% Default behaviour
%%   1) The vertical 1-D weathering and nitrogen figure (runs 10:12).
%%      Base cations, cumulative gaseous N loss and WFPS are panels that
%%      the summary figure does not carry.  Prints the years-2-onward
%%      nitrogen statistics.
%%   2) Separate summary and time-dynamics figures for each metric.  Each
%%      compares control, basalt and calcite in every available 1-D and
%%      2-D configuration.
%%
%% PNG and FIG copies of every figure are written to v2/figures.
%%
%% Stage 2 is skipped with a warning listing the missing runs, so the
%% script stays runnable while the 2-D parallel runs are still in
%% progress.  Stage 1 raises when a 1-D result is missing.
%%
%% Optional workspace variables, set before running the script:
%%   Runs             complete treatment groups        [1:12 with the
%%                    single-group summary and dynamics stage; empty
%%                    means the two-stage default]
%%   ResultDirectory  folder holding the run results     [SIM/Subfield_12]
%%   FigureDirectory  folder the figures are written to           [figures]
%%   NDays            shortened-run filename tag       []
%%   Scenario         scenario filename tag, e.g. 'ssp245'; results are
%%                    then TC_2D_run%02d_<name>__<Scenario>[__%dd].mat  ['']
%%   UncertaintyFile  uncertainty_summary.mat of the 1-D parameter
%%                    ensemble (Summarize_Uncertainty); draws
%%                    the member min-max as bands under the 1-D lines
%%                    and as range bars on the 1-D summary bars       ['']
%%   Visible          figure visibility                ['on']
%%   OutputFile       summary PNG filename prefix      []
%%   TimeOutputFile   dynamics PNG filename prefix     []
%%   DryRun           print the plan without plotting  [false]
%%
%% Examples
%%   matlab -batch "cd('v2/plotting'); Plot_1D_Then_All_12"
%%   matlab -batch "cd('v2/plotting'); NDays=1; Plot_1D_Then_All_12"
%%   matlab -batch "cd('v2/plotting'); Runs=10:12; Plot_1D_Then_All_12"
%%   matlab -batch "cd('v2/plotting'); ResultDirectory='../SIM/Projection'; Scenario='ssp245'; Plot_1D_Then_All_12"
%%
%% Result variables are left in the workspace: summary, six-element fig
%% and time_fig arrays, six-element plot_file and time_plot_file cells,
%% and fig1d, plot_file_1d, stats_1d for the 1-D weathering stage.
%%
%% NEE follows the TC sign convention: positive is carbon released to the
%% atmosphere, normalised by the full domain area.  The comparison figures
%% show all three treatments in each available configuration.  Chamber
%% dates from figure3a_n2o_digitised.csv define the observation-day N2O mean.
%% The observed bulk pH carries the SEM of the block means read from
%% Beerling2024_Fig1E_vector_extraction.csv.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% This script lives in v2/plotting; the repository root is one level up.
%% The launchers in v2/run and the model code are put on the path here so
%% the script runs from any folder.
here = fileparts(mfilename('fullpath'));
if isempty(here); here = pwd; end
root = fileparts(here);
addpath(fullfile(root,'run'),fullfile(root,'T&C_Code'), ...
    fullfile(root,'T&C_Code','Parameters'));
if ~exist('ResultDirectory','var') || isempty(ResultDirectory)
    ResultDirectory = fullfile(root,'SIM','Subfield_12');
end
if ~exist('FigureDirectory','var') || isempty(FigureDirectory)
    FigureDirectory = fullfile(root,'figures');
end
if ~exist('Runs','var'); Runs = []; end
if ~exist('NDays','var'); NDays = []; end
if ~exist('Scenario','var') || isempty(Scenario); Scenario = ''; end
if ~exist('UncertaintyFile','var') || isempty(UncertaintyFile); UncertaintyFile = ''; end
if ~exist('Visible','var') || isempty(Visible); Visible = 'on'; end
if ~exist('OutputFile','var') || isempty(OutputFile); OutputFile = ''; end
if ~exist('TimeOutputFile','var') || isempty(TimeOutputFile); TimeOutputFile = ''; end
if ~exist('DryRun','var') || isempty(DryRun); DryRun = false; end

validate_options(NDays,Scenario,UncertaintyFile,Visible,DryRun,ResultDirectory);
if ~isfolder(FigureDirectory); mkdir(FigureDirectory); end
opt = struct('ResultDirectory',ResultDirectory, ...
    'FigureDirectory',FigureDirectory,'NDays',NDays, ...
    'Scenario',char(Scenario),'Uncertainty',load_uncertainty(UncertaintyFile), ...
    'Visible',Visible,'OutputFile',OutputFile, ...
    'TimeOutputFile',TimeOutputFile,'Root',root);

summary = []; fig = []; plot_file = []; time_fig = []; time_plot_file = [];
fig1d = []; plot_file_1d = []; stats_1d = [];

if isempty(Runs)
    %% Two-stage default: 1-D group first, then all twelve.
    one_d_runs = 10:12;
    if DryRun
        fprintf('(dry run) stage 1: 1-D weathering figure, runs 10:12\n');
        fprintf('(dry run) stage 2: summary and dynamics, runs 1:12\n');
        return
    end
    missing_1d = missing_results(one_d_runs,opt);
    if ~isempty(missing_1d)
        error('Plot_1D_Then_All_12:MissingResults', ...
            'Cannot plot the 1-D group without results for run(s) %s.', ...
            mat2str(missing_1d));
    end
    [fig1d,plot_file_1d,stats_1d] = plot_group_figure(one_d_runs,opt);
    missing_12 = missing_results(1:12,opt);
    if ~isempty(missing_12)
        fprintf(['\nSkipping the all-12 plot: missing results for run(s) %s.\n', ...
            'Plot again once the parallel runs finish.\n'],mat2str(missing_12));
        return
    end
    [summary,fig,plot_file,time_fig,time_plot_file] = ...
        plot_subfield(1:12,opt);
    return
end

%% Single group: summary and time dynamics for the requested runs.
run_ids = validate_runs(Runs);
if DryRun
    fprintf('(dry run) summary and dynamics, runs %s\n',mat2str(run_ids));
    return
end
[summary,fig,plot_file,time_fig,time_plot_file] = plot_subfield(run_ids,opt);

%% =======================================================================
function run_ids = validate_runs(Runs)
run_ids = reshape(Runs,1,[]);
if ~(isnumeric(run_ids) && isreal(run_ids) && ~isempty(run_ids) && ...
        all(isfinite(run_ids)) && all(run_ids == fix(run_ids)) && ...
        all(run_ids >= 1 & run_ids <= 12))
    error('Plot_1D_Then_All_12:BadRuns', ...
        'Runs must contain complete treatment groups from 1 to 12.');
end
group_ids = unique(ceil(run_ids/3),'stable');
expected_ids = reshape(3*group_ids+(-2:0).',1,[]);
if ~isequal(run_ids,expected_ids)
    error('Plot_1D_Then_All_12:BadRuns', ...
        'Runs must contain complete treatment groups in matrix order.');
end
end

%% =======================================================================
function validate_options(NDays,Scenario,UncertaintyFile,Visible,DryRun,ResultDirectory)
if ~isempty(NDays) && ~(isnumeric(NDays) && isscalar(NDays) && ...
        isfinite(NDays) && NDays == fix(NDays) && NDays >= 1)
    error('Plot_1D_Then_All_12:BadOptions','NDays must be a positive integer.');
end
if ~(ischar(Scenario) || (isstring(Scenario) && isscalar(Scenario)))
    error('Plot_1D_Then_All_12:BadOptions', ...
        'Scenario must be a character vector or scalar string.');
end
if ~(ischar(UncertaintyFile) || (isstring(UncertaintyFile) && isscalar(UncertaintyFile)))
    error('Plot_1D_Then_All_12:BadOptions', ...
        'UncertaintyFile must be a character vector or scalar string.');
end
if ~isempty(UncertaintyFile) && ~isfile(UncertaintyFile)
    error('Plot_1D_Then_All_12:MissingUncertainty', ...
        'Uncertainty summary not found: %s',UncertaintyFile);
end
if ~(ischar(Visible) && ismember(char(Visible),{'on','off'}))
    error('Plot_1D_Then_All_12:BadOptions','Visible must be on or off.');
end
if ~(islogical(DryRun) && isscalar(DryRun))
    error('Plot_1D_Then_All_12:BadOptions','DryRun must be scalar logical.');
end
if ~isfolder(ResultDirectory)
    error('Plot_1D_Then_All_12:MissingDirectory', ...
        'Result directory not found: %s',ResultDirectory);
end
end

%% =======================================================================
function [fig,plot_file,stats] = plot_group_figure(run_ids,opt)
%% Weathering signal and nitrogen partitioning for one treatment group.
%%
%% Former Plot_1D_Group.  This figure complements the summary figure
%% rather than repeating it: base cations Ca, Mg and Si, cumulative
%% gaseous N loss (N2 + N2O), and water-filled pore space live here.
%% The denitrification N2O fraction is reported from the saved total N2O,
%% nitrification N2O and N2 fluxes. It is not plotted as a daily curve.
%%
%% Colours, line weights, axis treatment, application markers and legend
%% placement follow the summary figure so the two figures read as one set.
if numel(run_ids) ~= 3
    error('Plot_1D_Then_All_12:BadRuns', ...
        'Expected one complete treatment group of three runs, got %d.', ...
        numel(run_ids));
end
RM = Run_Matrix_12();

treatment_labels = {'Control','Basalt','Calcite'};
treatment_colors = [0 0 0; 0 0.25 0.60; 0.35 0.70 0.90];
application_color = [0.45 0.45 0.45];

%% Pool column indices in P.
iCa = 60; iMg = 68; iSi = 76;

suffix = filename_suffix(opt);

Date_d = datetime.empty(0,1);
for i = 1:3
    id = run_ids(i);
    result_file = result_filename(opt,id,RM(id).name);
    if ~isfile(result_file)
        error('Plot_1D_Then_All_12:MissingResult','Missing result: %s',result_file);
    end
    S = load(result_file);
    P = S.P; if ndims(P) == 3; P = P(:,:,1); end
    nd = size(P,1);
    if i == 1
        %% Daily row j holds the BIOGEO day running from 01:00 of
        %% Date(1)+(j-3); row 2 is the lone first hour.
        Date_d = datetime(S.Datam(1,1),S.Datam(1,2),S.Datam(1,3)) + ...
            caldays((0:nd-1).'-2);
        ca = zeros(nd,3); mg = zeros(nd,3); si = zeros(nd,3);
        gaseous_n = zeros(nd,3); n2o = zeros(nd,3); wfps = zeros(nd,3);
        n2o_nit = zeros(nd,3); n2 = zeros(nd,3);
    end
    ca(:,i) = P(:,iCa); mg(:,i) = P(:,iMg); si(:,i) = P(:,iSi);
    n2o(:,i) = required_daily(S,'N2Oflx',nd);
    n2o_nit(:,i) = required_daily(S,'N2Oflx_nit',nd);
    n2(:,i)  = required_daily(S,'N2flx',nd);
    gaseous_n(:,i) = n2(:,i) + n2o(:,i);
    wfps(:,i) = Daily_WFPS_From_Hourly(S.O,S.p,nd);
end
check_uncertainty_axis(opt.Uncertainty,Date_d);

%% The saved N2O flux includes nitrification. Subtract that source before
%% calculating the N2O/N2 split for denitrification.
nd_all = size(gaseous_n,1);
eq = min(nd_all,round(365.25)+1):nd_all;
rN2O_fw = zeros(1,3);
for i = 1:3
    rN2O_fw(i) = denitrification_n2o_share( ...
        n2o(eq,i),n2o_nit(eq,i),n2(eq,i));
end

time_value = {ca,mg,si,cumsum(gaseous_n,1),wfps};
time_label = {'Ca solution','Mg solution','Si solution', ...
    'Cumulative gaseous N loss (N_2 + N_2O)','Water-filled pore space'};
time_units = {'g Ca m^{-2}','g Mg m^{-2}','g Si m^{-2}', ...
    'g N m^{-2}','%'};
time_tag = {'CaSolution','MgSolution','SiSolution', ...
    'CumulativeGaseousNLoss','WFPS'};

application_dates = datetime([2016 11 10; 2017 11 10; ...
    2018 11 10; 2019 11 10]);
application_dates = application_dates( ...
    application_dates >= Date_d(1) & application_dates <= Date_d(end));
management = kantola_management(Date_d);

group_label = strrep(RM(run_ids(1)).layout,'_','\_');
fig = figure('Color','w','Visible',char(opt.Visible), ...
    'Name','T&C-SMEW-N one-group weathering and nitrogen partitioning');
fig.Position(3:4) = [1150 720];
%% Two rows: the three base cations that carry the weathering signal, then
%% the two nitrogen panels. A six-column grid with spans lets a row of three
%% sit above a row of two without leaving an empty tile.
tl = tiledlayout(fig,2,6,'TileSpacing','compact','Padding','compact');
tile_span = {[1 2],[1 2],[1 2],[1 3],[1 3]};

time_axes = gobjects(5,1);
legend_lines = gobjects(3,1);
application_legend = gobjects(1);
fertilizer_legend = gobjects(1);
n_band = 0;
for measurement = 1:5
    ax = nexttile(tl,tile_span{measurement});
    time_axes(measurement) = ax;
    ax.Tag = sprintf('TimeMetric%d',measurement);
    n_band = n_band+draw_uncertainty_bands(ax,Date_d,opt.Uncertainty, ...
        time_tag{measurement},run_ids,treatment_colors);
    ln = plot(ax,Date_d,time_value{measurement},'LineWidth',1.15);
    for treatment = 1:3
        ln(treatment).Color = treatment_colors(treatment,:);
    end
    set(ln,'Tag',time_tag{measurement});
    if ~isempty(application_dates)
        marker = xline(ax,application_dates,':','Color',application_color, ...
            'LineWidth',1.1,'Tag','FeedstockApplication');
        if measurement == 1; application_legend = marker(1); end
    end
    if measurement >= 4
        fertilizer_marker = draw_kantola_management(ax,management);
        if measurement == 4 && isgraphics(fertilizer_marker)
            fertilizer_legend = fertilizer_marker;
        end
    end
    format_time_axis(ax);
    ylabel(ax,time_units{measurement},'Color','k');
    if measurement >= 4; xlabel(ax,'Date','Color','k'); end
    label_y = 0.96;
    if measurement >= 4; label_y = 0.88; end
    text(ax,0.02,label_y,time_label{measurement},'Units','normalized', ...
        'VerticalAlignment','top','Color','k','FontWeight','bold');
    if measurement == 1; legend_lines = ln; end
end
linkaxes(time_axes,'x');
set(time_axes(1:3),'XTickLabel',[]);

legend_handles = legend_lines;
legend_labels = treatment_labels;
if isgraphics(application_legend)
    legend_handles(end+1) = application_legend;
    legend_labels{end+1} = 'Feedstock application';
end
if isgraphics(fertilizer_legend)
    legend_handles(end+1) = fertilizer_legend;
    legend_labels{end+1} = 'Maize fertilizer';
end
lgd = legend(legend_handles,legend_labels,'Orientation','horizontal', ...
    'Box','off','TextColor','k','Color','w','EdgeColor','none');
lgd.Layout.Tile = 'south';
title(tl,sprintf('%s: weathering signal and nitrogen partitioning', ...
    group_label),'Color','k');
if n_band > 0; annotate_uncertainty(tl,opt.Uncertainty,'Shaded'); end

plot_file = fullfile(opt.FigureDirectory, ...
    sprintf('%s_weathering_nitrogen%s.png',RM(run_ids(1)).layout,suffix));
exportgraphics(fig,plot_file,'Resolution',200);
[plot_folder,plot_name] = fileparts(plot_file);
savefig(fig,fullfile(plot_folder,[plot_name '.fig']));
fprintf('Saved: %s\n',plot_file);
fprintf('Saved: %s\n',fullfile(plot_folder,[plot_name '.fig']));

%% Equilibrium statistics, excluding year one. The spin-up starts NO3 near
%% 0.737 gN/m2 and roughly doubles before settling, so a first-year mean
%% reports the transient rather than the model's own steady state.
nd = size(gaseous_n,1);
sl = min(nd,round(365.25)+1):nd;
f = 10/(numel(sl)/365.25);          %% gN/m2 over the window -> kg N/ha/yr
stats = struct();
fprintf('\nyears 2 onward\n');
fprintf('N fluxes: kg N/ha/yr; rN2O and WFPS: %%; Ca end: g Ca/m2\n');
fprintf('%-8s %7s %8s %10s %8s %9s\n', ...
    'run','N2O','gas N','rN2O','WFPS','Ca end');
for i = 1:3
    st = struct('treatment',treatment_labels{i}, ...
        'N2O',sum(n2o(sl,i),'omitnan')*f, ...
        'gaseous_N_loss',sum(gaseous_n(sl,i),'omitnan')*f, ...
        'rN2O',rN2O_fw(i), ...
        'wfps',mean(wfps(sl,i),'omitnan'), ...
        'Ca_end',ca(end,i),'Mg_end',mg(end,i),'Si_end',si(end,i));
    stats.(lower(treatment_labels{i})) = st;
    fprintf('%-8s %7.1f %8.1f %9.1f%% %7.1f%% %9.2f\n', ...
        st.treatment,st.N2O,st.gaseous_N_loss,st.rN2O,st.wfps,st.Ca_end);
end

end

%% =======================================================================
function [summary,fig,plot_file,time_fig,time_plot_file] = plot_subfield( ...
        run_ids,opt)
%% Twelve-scenario summary and daily dynamics for complete treatment
%% groups.  Former Plot_Subfield_12.  The 1-D group (runs 10:12) can be
%% plotted before the 2-D runs finish; the full set needs all 12 results.
RM_all = Run_Matrix_12();
RM = RM_all(run_ids);
n = numel(RM);
n2o = nan(n,1);
no3_leach = nan(n,1);
no3_final = nan(n,1);
ph_final = nan(n,1);
nee = nan(n,1);
crop_harvest = nan(n,1);
Date_d = datetime.empty(0,1);
n2o_daily = [];
no3_leach_daily = [];
no3_pool_daily = [];
nee_daily = [];
crop_harvest_daily = [];
has_bulk_pore_ph = true;
carbon_fields = {'NPP_H','NPP_L','Ccrown_t','R_litter','R_microbe','R_ew'};
result_fields = {'R','Domain','Datam','RB_L','topo','ColCfg'};
for i = 1:n
    run_id = RM(i).id;
    result_file = result_filename(opt,run_id,RM(i).name);
    if ~isfile(result_file)
        error('Plot_1D_Then_All_12:MissingResult', ...
            'Result not found: %s',result_file);
    end
    saved = whos('-file',result_file);
    saved_names = {saved.name};
    if ~all(ismember(carbon_fields,saved_names))
        error('Plot_1D_Then_All_12:IncompleteCarbonResult', ...
            'Carbon state is incomplete in %s.',result_file);
    end
    if ~all(ismember(result_fields,saved_names))
        error('Plot_1D_Then_All_12:IncompleteResult', ...
            'Required result state is incomplete in %s.',result_file);
    end
    S = load(result_file,result_fields{:},carbon_fields{:});
    if ~isfield(S,'R') || ~isfield(S,'Domain') || ...
            S.R.id ~= run_id || ~strcmp(S.R.name,RM(i).name)
        error('Plot_1D_Then_All_12:WrongResult', ...
            'Result does not match run %d: %s',run_id,result_file);
    end
    %% PH_bulk_d and PH_pore_d get the same completeness, shape, length and
    %% finiteness checks as the other Domain series whenever a run carries
    %% them, instead of reaching format_ph_axis unchecked: an empty or
    %% mismatched pair used to fail as a generic assignment error or a
    %% NaN ylim rather than one of this function's own error identifiers.
    run_has_bulk_pore_ph = isfield(S.Domain,'PH_bulk_d') && isfield(S.Domain,'PH_pore_d');
    required = {'N2Oflx','LEAK_NO3','NO3_pool','PH_d'};
    if run_has_bulk_pore_ph
        required = [required {'PH_bulk_d','PH_pore_d'}];
    end
    if ~all(isfield(S.Domain,required)) || ~isfield(S,'Datam') || ...
            isempty(S.Domain.NO3_pool) || isempty(S.Domain.PH_d) || ...
            (run_has_bulk_pore_ph && (isempty(S.Domain.PH_bulk_d) || ...
            isempty(S.Domain.PH_pore_d)))
        error('Plot_1D_Then_All_12:IncompleteResult', ...
            'Domain summary is incomplete in %s.',result_file);
    end
    measurement = {S.Domain.N2Oflx,S.Domain.LEAK_NO3, ...
        S.Domain.NO3_pool,S.Domain.PH_d};
    if run_has_bulk_pore_ph
        measurement = [measurement {S.Domain.PH_bulk_d,S.Domain.PH_pore_d}];
    end
    if ~all(cellfun(@(x) isnumeric(x) && isvector(x),measurement))
        error('Plot_1D_Then_All_12:IncompleteResult', ...
            'Domain measurements must be numeric vectors in %s.',result_file);
    end
    n_time = cellfun(@numel,measurement);
    if any(n_time ~= n_time(1))
        error('Plot_1D_Then_All_12:TimeLengthMismatch', ...
            'Domain measurements have different daily lengths in %s.',result_file);
    end
    if any(cellfun(@(x) any(~isfinite(x(:))),measurement))
        error('Plot_1D_Then_All_12:NonfiniteMeasurement', ...
            'Run %d does not contain finite values for every plotted measurement.',run_id);
    end
    if ~isnumeric(S.Datam) || isempty(S.Datam) || size(S.Datam,2) < 3 || ...
            any(~isfinite(S.Datam(1,1:3)))
        error('Plot_1D_Then_All_12:IncompleteResult', ...
            'Daily start date is missing or invalid in %s.',result_file);
    end
    run_date = datetime(S.Datam(1,1),S.Datam(1,2),S.Datam(1,3)) + ...
        caldays((0:n_time(1)-1).'-2);
    [run_nee,run_crop_harvest] = carbon_dynamics( ...
        S,n_time(1),result_file);
    has_bulk_pore_ph = has_bulk_pore_ph && run_has_bulk_pore_ph;
    if i == 1
        Date_d = run_date;
        n2o_daily = zeros(n_time(1),n);
        no3_leach_daily = zeros(n_time(1),n);
        no3_pool_daily = zeros(n_time(1),n);
        ph_bulk_daily = zeros(n_time(1),n);
        ph_pore_daily = zeros(n_time(1),n);
        nee_daily = zeros(n_time(1),n);
        crop_harvest_daily = zeros(n_time(1),n);
    elseif ~isequal(run_date,Date_d)
        error('Plot_1D_Then_All_12:TimeAxisMismatch', ...
            'Run %d does not use the same daily time axis as run %d.', ...
            run_id,RM(1).id);
    end
    n2o_daily(:,i) = S.Domain.N2Oflx(:); %#ok<AGROW>
    no3_leach_daily(:,i) = S.Domain.LEAK_NO3(:); %#ok<AGROW>
    no3_pool_daily(:,i) = S.Domain.NO3_pool(:); %#ok<AGROW>
    if run_has_bulk_pore_ph
        ph_bulk_daily(:,i) = S.Domain.PH_bulk_d(:);
        ph_pore_daily(:,i) = S.Domain.PH_pore_d(:);
    else
        ph_bulk_daily(:,i) = S.Domain.PH_d(:);
        ph_pore_daily(:,i) = S.Domain.PH_d(:);
    end
    nee_daily(:,i) = run_nee; %#ok<AGROW>
    crop_harvest_daily(:,i) = run_crop_harvest; %#ok<AGROW>
    n2o(i) = sum(n2o_daily(:,i));
    no3_leach(i) = sum(no3_leach_daily(:,i));
    no3_final(i) = no3_pool_daily(end,i);
    ph_final(i) = ph_bulk_daily(end,i);
    nee(i) = sum(nee_daily(:,i));
    crop_harvest(i) = sum(crop_harvest_daily(:,i));
end
check_uncertainty_axis(opt.Uncertainty,Date_d);

summary = table([RM.id].',string({RM.layout}).',string({RM.treatment}).', ...
    n2o,no3_leach,no3_final,ph_final,nee,crop_harvest, ...
    'VariableNames',{'RunID','ScenarioGroup','Treatment', ...
    'CumulativeN2O_gN_m2','CumulativeNO3Leaching_gN_m2', ...
    'FinalNO3Pool_gN_m2','FinalBulkPH', ...
    'CumulativeNEE_gC_m2','CumulativeCropYield_gC_m2'});

%% Digitised Val Martin et al. chamber measurements for maize 2016-2018.
%% The field data are reported in g N ha^-1 d^-1, so the modelled flux is
%% carried in those units wherever the two are shown together.
n2o_obs = n2o_observations(opt.Root);
n2o_obs = clip_observations(n2o_obs,Date_d);
%% Beerling bulk-pH observations by year, no date component: these compare
%% against the pH time-dynamics figure as a categorical points-and-years
%% panel rather than against a specific day, so clip_observations (which
%% keys on obs.date) does not apply.
bulk_ph_obs = bulk_pH_observations(opt.Root);
n2o_ha_daily = 1e4*n2o_daily;
[n2o_obs_day_mean,~,n_obs_day] = ...
    paired_n2o_means(Date_d,n2o_ha_daily,n2o_obs);
fprintf(['Crop-matched N2O observation days inside the run window: ' ...
    '%d\n'],n_obs_day);

all_group_labels = {'Whole field','Margin with crop', ...
    'Margin with prairie','Vertical 1-D'};
all_group_stems = {'Whole_Field','Margin_Crop','Margin_Prairie','Vertical_1D'};
group_ids = unique(ceil(run_ids/3),'stable');
group_labels = all_group_labels(group_ids);
if isequal(run_ids,1:12)
    output_stem = 'Subfield_12';
    summary_name = 'T&C-SMEW-N twelve-scenario summary';
    time_name = 'T&C-SMEW-N twelve-scenario time dynamics';
elseif isequal(run_ids,1:9)
    output_stem = 'Subfield_2D';
    summary_name = 'T&C-SMEW-N 2-D scenario summary';
    time_name = 'T&C-SMEW-N 2-D scenario time dynamics';
elseif isscalar(group_ids)
    output_stem = all_group_stems{group_ids};
    summary_name = sprintf('T&C-SMEW-N %s summary',group_labels{1});
    time_name = sprintf('T&C-SMEW-N %s time dynamics',group_labels{1});
else
    run_tag = sprintf('%02d_',run_ids);
    output_stem = ['Subfield_runs_' run_tag(1:end-1)];
    summary_name = 'T&C-SMEW-N selected-scenario summary';
    time_name = 'T&C-SMEW-N selected-scenario time dynamics';
end
suffix = filename_suffix(opt);
comparison = treatment_configurations(RM);
uncertainty_range = summary_uncertainty_ranges(opt.Uncertainty,Date_d,n2o_obs);

summary_value = {n2o,no3_leach,no3_final, ...
    ph_final,nee,n2o_obs_day_mean};
summary_label = {'Cumulative N_2O','Cumulative NO_3 leaching', ...
    'Final soil NO_3 pool','Final bulk pH', ...
    'Cumulative NEE (+ to atmosphere)', ...
    'Mean N_2O on observation days'};
summary_units = {'g N m^{-2}','g N m^{-2}','g N m^{-2}', ...
    'pH','g C m^{-2}','g N ha^{-1} d^{-1}'};
summary_tag = {'CumulativeN2O','CumulativeNO3Leaching', ...
    'NO3Pool','FinalBulkPH','CumulativeNEE', ...
    'ObservationDayN2O'};
summary_stem = {'cumulative_n2o','cumulative_no3_leaching', ...
    'soil_no3_pool','final_bulk_ph','cumulative_nee', ...
    'observation_day_n2o'};
fig = gobjects(6,1);
plot_file = cell(6,1);
for measurement = 1:6
    fig(measurement) = figure('Color','w','Visible',char(opt.Visible), ...
        'Name',sprintf('%s: %s',summary_name,summary_label{measurement}));
    fig(measurement).Position(3:4) = [1050 700];
    ax = axes('Parent',fig(measurement));
    value = summary_value{measurement}(comparison.columns);
    b = draw_treatment_bars(ax,value,comparison.treatment_colors, ...
        comparison.run_ids,summary_tag{measurement});
    n_range = draw_uncertainty_ranges(ax,b,comparison.run_ids, ...
        opt.Uncertainty,uncertainty_range{measurement});
    format_configuration_axis(ax,comparison.configuration_labels);
    title(ax,summary_label{measurement},'Color','k');
    if n_range > 0; annotate_uncertainty(ax,opt.Uncertainty,'Range bars'); end
    ylabel(ax,summary_units{measurement},'Color','k');
    lgd = legend(ax,b,comparison.treatment_labels, ...
        'Location','southoutside','Orientation','horizontal', ...
        'Box','off','TextColor','k');
    lgd.AutoUpdate = 'off';
    plot_file{measurement} = metric_output_file(opt.OutputFile, ...
        opt.FigureDirectory,[output_stem '_summary'], ...
        summary_stem{measurement},suffix);
    export_figure_pair(fig(measurement),plot_file{measurement});
end

time_value = {cumsum(n2o_daily,1),cumsum(no3_leach_daily,1), ...
    no3_pool_daily,ph_bulk_daily,nee_daily,n2o_ha_daily};
time_label = {'Cumulative N_2O','Cumulative NO_3 leaching', ...
    'Soil NO_3 pool','pH','NEE (+ to atmosphere)','N_2O'};
time_units = {'g N m^{-2}','g N m^{-2}','g N m^{-2}','pH', ...
    'g C m^{-2} d^{-1}','g N ha^{-1} d^{-1}'};
time_tag = {'CumulativeN2O','CumulativeNO3Leaching','NO3Pool', ...
    'BulkPH','NEE','DailyN2O'};
envelope_field = {'CumulativeN2O','CumulativeNO3Leaching','NO3Pool', ...
    'BulkPH','NEE','DailyN2O_gN_ha_d'};
time_stem = {'cumulative_n2o','cumulative_no3_leaching', ...
    'soil_no3_pool','ph','nee','daily_n2o'};
trace_colors = comparison.treatment_colors(comparison.trace_treatment,:);
application_dates = datetime([2016 11 10; 2017 11 10; ...
    2018 11 10; 2019 11 10]);
application_dates = application_dates( ...
    application_dates >= Date_d(1) & application_dates <= Date_d(end));
management = kantola_management(Date_d);
application_color = [0.45 0.45 0.45];
time_fig = gobjects(6,1);
time_plot_file = cell(6,1);
for measurement = 1:6
    time_fig(measurement) = figure('Color','w', ...
        'Visible',char(opt.Visible), ...
        'Name',sprintf('%s: %s',time_name,time_label{measurement}));
    time_fig(measurement).Position(3:4) = [1400 900];
    %% The pH figure carries bulk and pore-water pH in two stacked panels,
    %% each scaled to its own data range, so the 0.4-0.5 unit offset between
    %% them no longer leaves half of a shared axis empty. pH and daily N2O
    %% additionally carry a bottom panel comparing modelled years against
    %% the digitised field observations: one point per year, never a line,
    %% since a chamber campaign or a yearly Beerling sample is a
    %% measurement, not a trend to interpolate through.
    split_ph = (measurement == 4) && has_bulk_pore_ph;
    show_ph_obs = (measurement == 4) && bulk_ph_obs.available;
    %% The chamber N2O record is the exception: it is dated at sub-weekly
    %% resolution against the same axis the model is drawn on, so it goes
    %% straight onto the modelled traces rather than into a panel of its
    %% own, and adds no row here.
    show_n2o_obs = (measurement == 6) && ~isempty(n2o_obs.date);
    n_rows = 1+double(split_ph)+double(show_ph_obs);
    row = 0;
    if split_ph
        %% subplot, not tiledlayout: a tiledlayout's OuterPosition is not
        %% honoured under 'Padding','compact' in this MATLAB release (the
        %% child tiles keep rendering flush to the figure edge regardless
        %% of the OuterPosition set on the layout), which left no real
        %% room for the legend band added below. subplot axes are plain
        %% axes with the ordinary auto-margin behaviour the legend
        %% placement below relies on.
        row = row+1; ax = subplot(n_rows,1,row,'Parent',time_fig(measurement));
        row = row+1; ax_pore = subplot(n_rows,1,row,'Parent',time_fig(measurement));
    elseif n_rows > 1
        row = row+1; ax = subplot(n_rows,1,row,'Parent',time_fig(measurement));
    else
        row = row+1; ax = axes('Parent',time_fig(measurement));
    end
    n_band = draw_uncertainty_bands(ax,Date_d,opt.Uncertainty, ...
        envelope_field{measurement},comparison.trace_ids,trace_colors);
    ln = plot(ax,Date_d,time_value{measurement}(:,comparison.trace_columns), ...
        'LineWidth',1.35);
    hold(ax,'on');
    for trace = 1:numel(comparison.trace_columns)
        treatment = comparison.trace_treatment(trace);
        configuration = comparison.trace_configuration(trace);
        ln(trace).Color = comparison.treatment_colors(treatment,:);
        ln(trace).LineStyle = comparison.configuration_styles{configuration};
        ln(trace).Tag = time_tag{measurement};
        ln(trace).UserData = comparison.trace_ids(trace);
    end
    panel_axes = ax;
    if split_ph
        draw_uncertainty_bands(ax_pore,Date_d,opt.Uncertainty, ...
            'PoreWaterPH',comparison.trace_ids,trace_colors);
        pore_lines = plot(ax_pore,Date_d, ...
            ph_pore_daily(:,comparison.trace_columns), ...
            'LineWidth',1.35,'Tag','PoreWaterPH');
        hold(ax_pore,'on');
        for trace = 1:numel(comparison.trace_columns)
            treatment = comparison.trace_treatment(trace);
            configuration = comparison.trace_configuration(trace);
            pore_lines(trace).Color = comparison.treatment_colors(treatment,:);
            pore_lines(trace).LineStyle = ...
                comparison.configuration_styles{configuration};
            pore_lines(trace).UserData = comparison.trace_ids(trace);
        end
        panel_axes = [ax ax_pore];
        %% Tight y-limits from each panel's own data (bands included), set
        %% before the management overlay adds its label headroom.
        format_ph_axis(ax,with_uncertainty(ph_bulk_daily(:,comparison.trace_columns), ...
            opt.Uncertainty,'BulkPH',comparison.trace_ids));
        format_ph_axis(ax_pore,with_uncertainty(ph_pore_daily(:,comparison.trace_columns), ...
            opt.Uncertainty,'PoreWaterPH',comparison.trace_ids));
    elseif measurement == 4
        format_ph_axis(ax,with_uncertainty(ph_bulk_daily(:,comparison.trace_columns), ...
            opt.Uncertainty,'BulkPH',comparison.trace_ids));
    end
    for panel = panel_axes
        if ~isempty(application_dates)
            marker = xline(panel,application_dates,':', ...
                'Color',application_color,'LineWidth',1.1, ...
                'Tag','FeedstockApplication');
            set(marker,'HandleVisibility','off');
        end
        fertilizer_marker = draw_kantola_management(panel,management);
        if isgraphics(fertilizer_marker)
            fertilizer_marker.HandleVisibility = 'off';
        end
        format_time_axis(panel);
        ylabel(panel,time_units{measurement},'Color','k');
    end
    if split_ph
        ax.XTickLabel = {};
        title(ax,'Bulk pH','Color','k');
        title(ax_pore,'Pore-water pH','Color','k');
        xlabel(ax_pore,'Date','Color','k');
        %% The two panels share Date_d but keep independent y-limits;
        %% link only x so a zoom or pan on one follows through on the
        %% other instead of leaving the blanked top ticks stale.
        linkaxes(panel_axes,'x');
    else
        xlabel(ax,'Date','Color','k');
        title(ax,time_label{measurement},'Color','k');
    end
    if n_band > 0; annotate_uncertainty(ax,opt.Uncertainty,'Shaded'); end
    if show_n2o_obs
        %% Overlaid on the modelled traces, not stacked beneath them: the
        %% comparison being made is per-day, and two panels sharing an x
        %% axis still leave the reader estimating vertical distance by eye
        %% across a gap. On one axes an observed spike and the modelled
        %% value for that day are read off the same gridline.
        overlay_n2o_observations(ax,n2o_obs,comparison.treatment_colors);
    end
    if show_ph_obs
        row = row+1; ax_obs = subplot(n_rows,1,row,'Parent',time_fig(measurement));
        draw_pH_comparison(ax_obs,[],[],bulk_ph_obs,comparison.treatment_colors);
    end
    %% Panels in top-to-bottom order for the layout pass below.
    all_panels = ax;
    if split_ph; all_panels = [all_panels ax_pore]; end
    if show_ph_obs; all_panels = [all_panels ax_obs]; end
    %% A single legend spelling out all twelve "configuration treatment"
    %% combinations rendered cleanly in local testing but overlapped
    %% illegibly on Bouchet's MATLAB/Linux stack regardless of
    %% Orientation, NumColumns or Location -- most likely a font-metrics
    %% or legend-layout difference between MATLAB releases/platforms that
    %% only shows up once the auto-column layout has to fit this many
    %% long labels. Two short single-row legends -- colour for treatment,
    %% line style for configuration -- convey the same mapping without
    %% ever exercising that multi-column path, and are positioned
    %% explicitly rather than via 'Location', which is what triggered the
    %% original layout.
    %% Legends are built first and the panels laid out afterwards. The band
    %% the legends need is measured from their rendered row height rather
    %% than reserved by a fixed guess, so the panels can then expand into
    %% exactly the space left above it -- no empty strip, and no legend
    %% crowded against the axis it sits under.
    obs_labels = {}; obs_colors = zeros(0,3);
    if show_n2o_obs
        obs_labels = {'Observed control','Observed basalt'};
        obs_colors = comparison.treatment_colors(1:2,:);
    end
    [lgd,lgd_configuration,lgd_observation,band_top] = add_split_legends( ...
        time_fig(measurement),ax, ...
        ln(1:numel(comparison.treatment_labels)),comparison.configuration_labels, ...
        comparison.configuration_styles,comparison.treatment_labels, ...
        comparison.treatment_colors,obs_labels,obs_colors);
    lgd.AutoUpdate = 'off';
    lgd_configuration.AutoUpdate = 'off';
    if isgraphics(lgd_observation); lgd_observation.AutoUpdate = 'off'; end
    layout_time_panels(all_panels,band_top);
    if split_ph
        hold(ax_pore,'off');
    end
    hold(ax,'off');
    time_plot_file{measurement} = metric_output_file( ...
        opt.TimeOutputFile,opt.FigureDirectory, ...
        [output_stem '_time_dynamics'],time_stem{measurement},suffix);
    export_figure_pair(time_fig(measurement),time_plot_file{measurement});
end

end

%% =======================================================================
function comparison = treatment_configurations(RM)
wanted_ids = [10 11 12; 1 2 3; 4 5 6; 7 8 9];
wanted_labels = {'1-D','2-D whole field','2-D buffer only', ...
    '2-D native vegetation buffer'};
wanted_styles = {'-','--',':','-.'};
available = all(ismember(wanted_ids,[RM.id]),2);
comparison.run_ids = wanted_ids(available,:);
comparison.configuration_labels = wanted_labels(available);
comparison.configuration_styles = wanted_styles(available);
comparison.treatment_labels = {'Control','Basalt','Calcite'};
comparison.treatment_colors = [0 0 0; 0 0.25 0.60; 0.35 0.70 0.90];
comparison.columns = zeros(size(comparison.run_ids));
for configuration = 1:size(comparison.run_ids,1)
    for treatment = 1:3
        comparison.columns(configuration,treatment) = find( ...
            [RM.id] == comparison.run_ids(configuration,treatment),1);
    end
end

n_trace = numel(comparison.run_ids);
comparison.trace_columns = zeros(1,n_trace);
comparison.trace_ids = zeros(1,n_trace);
comparison.trace_labels = cell(1,n_trace);
comparison.trace_configuration = zeros(1,n_trace);
comparison.trace_treatment = zeros(1,n_trace);
trace = 0;
for configuration = 1:size(comparison.run_ids,1)
    for treatment = 1:3
        trace = trace+1;
        comparison.trace_columns(trace) = ...
            comparison.columns(configuration,treatment);
        comparison.trace_ids(trace) = ...
            comparison.run_ids(configuration,treatment);
        comparison.trace_labels{trace} = sprintf('%s %s', ...
            comparison.configuration_labels{configuration}, ...
            lower(comparison.treatment_labels{treatment}));
        comparison.trace_configuration(trace) = configuration;
        comparison.trace_treatment(trace) = treatment;
    end
end
end

%% =======================================================================
function b = draw_treatment_bars(ax,value,colors,run_ids,tag)
if size(run_ids,1) == 1
    b = gobjects(3,1);
    position = 1+[-0.24 0 0.24];
    hold(ax,'on');
    for treatment = 1:3
        b(treatment) = bar(ax,position(treatment),value(treatment),0.22, ...
            'FaceColor',colors(treatment,:),'Tag',tag);
        b(treatment).UserData = run_ids(:,treatment).';
    end
    hold(ax,'off');
    xlim(ax,[0.5 1.5]);
else
    b = bar(ax,value,'grouped','Tag',tag);
    for treatment = 1:3
        b(treatment).FaceColor = colors(treatment,:);
        b(treatment).UserData = run_ids(:,treatment).';
    end
end
end

%% =======================================================================
function format_configuration_axis(ax,labels)
ax.Color = 'w';
ax.XColor = 'k';
ax.YColor = 'k';
ax.XTick = 1:numel(labels);
ax.XTickLabel = labels;
ax.XTickLabelRotation = 15;
grid(ax,'off');
end

%% =======================================================================
function output_file = metric_output_file(requested,output_dir,base,metric,suffix)
if isempty(requested)
    output_file = fullfile(output_dir, ...
        sprintf('%s_%s%s.png',base,metric,suffix));
    return
end
[folder,name] = fileparts(char(requested));
output_file = fullfile(folder,sprintf('%s_%s.png',name,metric));
end

%% =======================================================================
function export_figure_pair(fig,png_file)
exportgraphics(fig,png_file,'Resolution',200);
[folder,name] = fileparts(png_file);
fig_file = fullfile(folder,[name '.fig']);
savefig(fig,fig_file);
fprintf('Saved: %s\n',png_file);
fprintf('Saved: %s\n',fig_file);
end

%% =======================================================================
function missing = missing_results(run_ids,opt)
RM = Run_Matrix_12();
absent = arrayfun(@(id) ~isfile(result_filename(opt,id,RM(id).name)),run_ids);
missing = reshape(run_ids(absent),[],1);
end

%% =======================================================================
function result_file = result_filename(opt,run_id,name)
%% Results are TC_2D_run%02d_<name>[__<Scenario>][__%dd].mat, the tag
%% order Run_One_Subfield writes them with.
result_file = fullfile(opt.ResultDirectory, ...
    sprintf('TC_2D_run%02d_%s%s.mat',run_id,name,filename_suffix(opt)));
end

%% =======================================================================
function suffix = filename_suffix(opt)
suffix = '';
if ~isempty(opt.Scenario); suffix = sprintf('__%s',opt.Scenario); end
if ~isempty(opt.NDays); suffix = sprintf('%s__%dd',suffix,opt.NDays); end
end

%% =======================================================================
function U = load_uncertainty(file)
%% The 1-D parameter ensemble summary Summarize_Uncertainty
%% writes, or [] when none is requested.  Members are the second
%% dimension of the daily arrays; run_ids map the third to Run_Matrix_12.
U = [];
if isempty(file); return; end
U = load(char(file),'Date','run_ids','daily','envelope','summary_range','in_envelope');
if ~all(isfield(U,{'Date','run_ids','daily','envelope','summary_range','in_envelope'})) || ...
        ~isfield(U.daily,'DailyN2O_gN_ha_d') || ~islogical(U.in_envelope) || ...
        numel(U.in_envelope) ~= size(U.daily.DailyN2O_gN_ha_d,2) || ~U.in_envelope(1)
    error('Plot_1D_Then_All_12:BadUncertainty', ...
        'Uncertainty summary is incomplete: %s',file);
end
%% The bands and range bars span the members flagged for the envelope; the
%% sensitivity-only members are in the daily arrays but not in the spread.
U.in_envelope = U.in_envelope(:).';
U.n_members = nnz(U.in_envelope);
U.file = char(file);
end

%% =======================================================================
function check_uncertainty_axis(U,Date_d)
%% The ensemble is summarised on its own daily axis; the bands are only
%% meaningful under lines drawn on exactly the same days.
if isempty(U); return; end
ensemble_date = datetime(U.Date(:),'ConvertFrom','datenum');
if numel(ensemble_date) ~= numel(Date_d) || ...
        any(abs(days(ensemble_date-Date_d(:))) > 1e-6)
    error('Plot_1D_Then_All_12:UncertaintyAxis', ...
        'The uncertainty summary is not on the plotted daily axis: %s',U.file);
end
end

%% =======================================================================
function n_band = draw_uncertainty_bands(ax,Date_d,U,metric,run_ids,colors)
%% Member min-max envelope of every run the ensemble holds, as a
%% translucent band in the run's treatment colour, drawn before the lines
%% so they stay on top.  Days whose envelope is NaN are left out.
n_band = 0;
if isempty(U) || ~isfield(U.envelope,metric); return; end
env = U.envelope.(metric);
for k = 1:numel(run_ids)
    r = find(U.run_ids == run_ids(k),1);
    if isempty(r); continue; end
    low = env.min(:,r); high = env.max(:,r);
    ok = isfinite(low) & isfinite(high);
    if ~any(ok); continue; end
    hold(ax,'on');
    fill(ax,[Date_d(ok); flipud(Date_d(ok))],[low(ok); flipud(high(ok))], ...
        colors(k,:),'FaceAlpha',0.18,'EdgeColor','none', ...
        'Tag','UncertaintyBand','HandleVisibility','off');
    n_band = n_band+1;
end
end

%% =======================================================================
function value = with_uncertainty(value,U,metric,run_ids)
%% The envelope columns of the runs the ensemble holds appended to the
%% plotted values, so an axis scaled from the data also covers the bands.
if isempty(U) || ~isfield(U.envelope,metric); return; end
r = find(ismember(U.run_ids,run_ids));
value = [value U.envelope.(metric).min(:,r) U.envelope.(metric).max(:,r)];
end

%% =======================================================================
function range = summary_uncertainty_ranges(U,Date_d,n2o_obs)
%% Member min-max of the six summary quantities, in the summary figures'
%% order.  Five are stored; the observation-day N2O mean is derived per
%% member from the daily flux with the same day pairing as the bars.
range = cell(6,1);
if isempty(U); return; end
stored = {'CumulativeN2O_gN_m2','CumulativeNO3Leaching_gN_m2', ...
    'FinalNO3Pool_gN_m2','FinalBulkPH','CumulativeNEE_gC_m2'};
for k = 1:numel(stored)
    if isfield(U.summary_range,stored{k})
        range{k} = U.summary_range.(stored{k});
    end
end
n_runs = numel(U.run_ids);
obs_day = nan(U.n_members,n_runs);
for r = 1:n_runs
    obs_day(:,r) = paired_n2o_means(Date_d, ...
        U.daily.DailyN2O_gN_ha_d(:,U.in_envelope,r),n2o_obs);
end
range{6} = struct('min',min(obs_day,[],1),'max',max(obs_day,[],1));
end

%% =======================================================================
function n_range = draw_uncertainty_ranges(ax,b,run_ids,U,range)
%% Member min-max on the bars of the runs the ensemble holds, as a thin
%% black bar spanning exactly that range whatever the plotted value is.
n_range = 0;
if isempty(U) || isempty(range); return; end
was_held = ishold(ax);
hold(ax,'on');
for configuration = 1:size(run_ids,1)
    for treatment = 1:size(run_ids,2)
        r = find(U.run_ids == run_ids(configuration,treatment),1);
        if isempty(r); continue; end
        low = range.min(r); high = range.max(r);
        if ~(isfinite(low) && isfinite(high)); continue; end
        errorbar(ax,b(treatment).XEndPoints(configuration),(low+high)/2, ...
            (high-low)/2,'LineStyle','none','Color','k','LineWidth',0.9, ...
            'CapSize',6,'Tag','UncertaintyRange','HandleVisibility','off');
        n_range = n_range+1;
    end
end
if ~was_held; hold(ax,'off'); end
end

%% =======================================================================
function annotate_uncertainty(target,U,what)
%% One short note per figure that carries the ensemble: what is drawn
%% ('Shaded' bands or 'Range bars') and how many members set the range.
%% A subtitle, so it is neither a legend entry nor a panel label.
subtitle(target,sprintf('%s: parameter-uncertainty range (%d members)', ...
    what,U.n_members),'FontSize',8,'FontWeight','normal','Color',[0.3 0.3 0.3]);
end

%% =======================================================================
function [nee,crop_harvest] = carbon_dynamics( ...
        S,n_time,result_file)
required = {'NPP_H','NPP_L','Ccrown_t','R_litter','R_microbe', ...
    'R_ew','RB_L','topo','ColCfg'};
if ~all(isfield(S,required)) || ~isstruct(S.topo) || ...
        ~isfield(S.topo,'cell_area_m2') || ~isstruct(S.ColCfg)
    error('Plot_1D_Then_All_12:IncompleteCarbonResult', ...
        'Carbon state is incomplete in %s.',result_file);
end

A_cell = S.topo.cell_area_m2(:);
n_col = numel(A_cell);
if n_col < 1 || numel(S.ColCfg) ~= n_col || ...
        any(~isfinite(A_cell)) || any(A_cell <= 0)
    error('Plot_1D_Then_All_12:IncompleteCarbonResult', ...
        'Carbon cell areas or configuration are invalid in %s.',result_file);
end

crop_cell = false(n_col,1);
for c = 1:n_col
    if ~isfield(S.ColCfg(c),'PFT_OVERRIDE') || ...
            ~(ischar(S.ColCfg(c).PFT_OVERRIDE) || ...
            (isstring(S.ColCfg(c).PFT_OVERRIDE) && ...
            isscalar(S.ColCfg(c).PFT_OVERRIDE)))
        error('Plot_1D_Then_All_12:IncompleteCarbonResult', ...
            'Cell PFT configuration is incomplete in %s.',result_file);
    end
    crop_cell(c) = isempty(char(S.ColCfg(c).PFT_OVERRIDE));
end

if ~(isnumeric(S.NPP_H) && isnumeric(S.NPP_L) && ...
        isnumeric(S.Ccrown_t) && isnumeric(S.R_litter) && ...
        isnumeric(S.R_microbe) && isnumeric(S.R_ew))
    error('Plot_1D_Then_All_12:CarbonTimeLengthMismatch', ...
        'NEE arrays must be numeric in %s.',result_file);
end
if ~isnumeric(S.RB_L)
    error('Plot_1D_Then_All_12:CropTimeLengthMismatch', ...
        'Crop-removal arrays must be numeric in %s.',result_file);
end
cc = size(S.NPP_L,2);
if cc < 1 || ...
        size(S.NPP_H,1) ~= n_time || size(S.NPP_H,2) ~= cc || ...
        size(S.NPP_H,3) ~= n_col || ...
        size(S.NPP_L,1) ~= n_time || size(S.NPP_L,2) ~= cc || ...
        size(S.NPP_L,3) ~= n_col || ...
        size(S.Ccrown_t,1) ~= n_time || size(S.Ccrown_t,2) ~= cc || ...
        size(S.Ccrown_t,3) ~= n_col || ...
        size(S.R_litter,1) ~= n_time || size(S.R_litter,2) ~= 1 || ...
        size(S.R_litter,3) ~= n_col || ...
        size(S.R_microbe,1) ~= n_time || size(S.R_microbe,2) ~= 1 || ...
        size(S.R_microbe,3) ~= n_col || ...
        size(S.R_ew,1) ~= n_time || size(S.R_ew,2) ~= 1 || ...
        size(S.R_ew,3) ~= n_col
    error('Plot_1D_Then_All_12:CarbonTimeLengthMismatch', ...
        'NEE arrays do not match the daily domain dimensions in %s.',result_file);
end
if size(S.RB_L,1) ~= n_time || size(S.RB_L,2) ~= cc || ...
        size(S.RB_L,3) ~= 7 || size(S.RB_L,4) ~= n_col
    error('Plot_1D_Then_All_12:CropTimeLengthMismatch', ...
        'Crop-removal arrays do not match the daily domain dimensions in %s.',result_file);
end
if any(~isfinite(S.NPP_H(:))) || any(~isfinite(S.NPP_L(:))) || ...
        any(~isfinite(S.Ccrown_t(:))) || any(~isfinite(S.R_litter(:))) || ...
        any(~isfinite(S.R_microbe(:))) || any(~isfinite(S.R_ew(:)))
    error('Plot_1D_Then_All_12:NonfiniteCarbonMeasurement', ...
        'NEE arrays contain nonfinite values in %s.',result_file);
end
if any(~isfinite(S.RB_L(:)))
    error('Plot_1D_Then_All_12:NonfiniteCropMeasurement', ...
        'Crop-removal arrays contain nonfinite values in %s.',result_file);
end
if any(S.Ccrown_t(:) < 0) || any(S.Ccrown_t(:) > 1)
    error('Plot_1D_Then_All_12:UnphysicalCarbonMeasurement', ...
        'Crown fractions must lie between zero and one in %s.',result_file);
end

NPP_H = reshape(S.NPP_H,n_time,cc,n_col);
NPP_L = reshape(S.NPP_L,n_time,cc,n_col);
Ccrown_t = reshape(S.Ccrown_t,n_time,cc,n_col);
R_litter = reshape(S.R_litter,n_time,n_col);
R_microbe = reshape(S.R_microbe,n_time,n_col);
R_ew = reshape(S.R_ew,n_time,n_col);
RB_L = reshape(S.RB_L,n_time,cc,7,n_col);
nee_cell = zeros(n_time,n_col);
crop_harvest_cell = zeros(n_time,n_col);
for c = 1:n_col
    ecosystem_npp = sum(Ccrown_t(:,:,c).* ...
        (NPP_H(:,:,c)+NPP_L(:,:,c)),2);
    nee_cell(:,c) = R_litter(:,c)+R_microbe(:,c)+R_ew(:,c)-ecosystem_npp;
    if ~crop_cell(c); continue; end
    %% Both Kantola crops set fract_left_fr=0, so all B5 removal is yield.
    for k = 1:cc
        crown = Ccrown_t(:,k,c);
        fruit_removed = max(0,reshape(RB_L(:,k,5,c),n_time,1));
        crop_harvest_cell(:,c) = crop_harvest_cell(:,c) + ...
            crown.*fruit_removed;
    end
end
w_col = A_cell/sum(A_cell);
nee = nee_cell*w_col;
crop_harvest = crop_harvest_cell*w_col;
end

%% =======================================================================
function r = denitrification_n2o_share(n2o_total,n2o_nit,n2_flux)
n2o_denit = n2o_total-n2o_nit;
values = [n2o_total(:); n2o_nit(:); n2_flux(:)];
if any(~isfinite(values))
    error('Plot_1D_Then_All_12:NonfiniteNitrogenFlux', ...
        'The saved N2O or N2 flux contains a nonfinite value.');
end
tol = 100*eps(max([1; abs(values)]));
if any(n2o_total < -tol) || any(n2o_nit < -tol) || ...
        any(n2o_denit < -tol) || any(n2_flux < -tol)
    error('Plot_1D_Then_All_12:InvalidNitrogenFlux', ...
        'The saved fluxes imply a negative denitrification N2O or N2 flux.');
end
n2o_denit(abs(n2o_denit) <= tol) = 0;
n2_flux(abs(n2_flux) <= tol) = 0;
denit = n2o_denit+n2_flux;
denit_total = sum(denit);
if denit_total <= tol
    r = NaN;
else
    r = 100*sum(n2o_denit)/denit_total;
end
end

%% =======================================================================
function v = required_daily(S,name,nd)
if ~isfield(S,name)
    error('Plot_1D_Then_All_12:MissingDailyFlux', ...
        'The saved result does not contain %s.',name);
end
v = S.(name);
if ~isnumeric(v)
    error('Plot_1D_Then_All_12:BadDailyFlux', ...
        '%s must be numeric.',name);
end
v = v(:);
if numel(v) < nd
    error('Plot_1D_Then_All_12:ShortDailyFlux', ...
        '%s has %d values, but the result contains %d saved days.', ...
        name,numel(v),nd);
end
v = v(1:nd);
end

function obs = n2o_observations(root)
%% Digitised Val Martin et al. chamber N2O measurements from the Energy
%% Farm maize plots, 2016-2018, in g N ha^-1 d^-1.  The same file
%% also carries the traced CLM5 curves from that figure; only the measured
%% points are returned here.  A missing file is a warning rather than an
%% error, so the script stays runnable without the digitised data.
obs = struct('date',datetime.empty(0,1),'flux',zeros(0,1), ...
    'treatment',strings(0,1),'file','');
obs.file = fullfile(root,'Experimental observations', ...
    'Web plot digitizer','figure3a_n2o_digitised.csv');
if ~isfile(obs.file)
    warning('Plot_1D_Then_All_12:MissingObservations', ...
        'Digitised N2O observations not found: %s',obs.file);
    return
end
T = readtable(obs.file,'TextType','string');
required = {'series','source','year','day_of_year','n2o_flux_gN_ha_day'};
if ~all(ismember(required,T.Properties.VariableNames))
    error('Plot_1D_Then_All_12:BadObservations', ...
        'Digitised observations are missing a required column: %s',obs.file);
end
T = T(T.source == "Energy Farm measurement",:);
%% Kantola has soybean in 2016; only the shared maize years are comparable.
T = T(ismember(T.year,[2017 2018]),:);
%% The digitiser resolves the x axis to a fraction of a day; the field
%% campaign sampled on whole days, so the nearest day is the sample date.
sample_date = datetime(T.year,1,1) + caldays(round(T.day_of_year)-1);
treatment = strings(height(T),1);
treatment(T.series == "EF Cont") = "control";
treatment(T.series == "EF ERW") = "basalt";
keep = treatment ~= "";
[obs.date,order] = sort(sample_date(keep));
flux = T.n2o_flux_gN_ha_day(keep);
treatment = treatment(keep);
obs.flux = flux(order);
obs.treatment = treatment(order);
end

%% =======================================================================
function obs = bulk_pH_observations(root)
%% Beerling bulk-pH observations.  The 0-10 and 10-30 cm measurements stay
%% separate because the model has one 0-30 cm bulk-pH state.
obs = struct('available',false,'labels',strings(0,1));
obs_file = fullfile(root,'Experimental observations', ...
    'Web plot digitizer','Data_Beerling.mat');
if ~isfile(obs_file)
    warning('Plot_1D_Then_All_12:MissingObservations', ...
        'Digitised bulk-pH observations not found: %s',obs_file);
    return
end
S = load(obs_file,'pH_Exp');
required = {'X','control_0','control_10','basalt_0','basalt_10'};
if ~isfield(S,'pH_Exp') || ~isstruct(S.pH_Exp) || ...
        ~all(isfield(S.pH_Exp,required))
    error('Plot_1D_Then_All_12:BadPHObservations', ...
        'Digitised bulk-pH structure is incomplete in %s.',obs_file);
end
obs.labels = strip(string(S.pH_Exp.X(:)));
for k = 2:numel(required)
    field_name = required{k};
    value = S.pH_Exp.(field_name);
    if ~isnumeric(value) || ~isvector(value) || ...
            numel(value) ~= numel(obs.labels) || any(~isfinite(value(:)))
        error('Plot_1D_Then_All_12:BadPHObservations', ...
            'Digitised bulk-pH series %s is invalid in %s.', ...
            field_name,obs_file);
    end
    obs.(field_name) = double(value(:));
end
obs.sem = bulk_pH_sem(root,obs.labels);
obs.available = true;
end

%% =======================================================================
function sem = bulk_pH_sem(root,labels)
%% Standard errors of the Beerling block means, read from the vector
%% extraction of Fig. 1E and aligned to the digitised labels by period
%% ('pre', then the year).  An empty SEM in the file is a point whose bar
%% the figure does not draw; it stays NaN and gets no bar here either.
series = {'control_0','control','0-10'; 'control_10','control','10-30'; ...
    'basalt_0','basalt','0-10'; 'basalt_10','basalt','10-30'};
sem = struct();
for k = 1:size(series,1); sem.(series{k,1}) = nan(numel(labels),1); end
csv_file = fullfile(root,'Experimental observations', ...
    'Web plot digitizer','Beerling2024_Fig1E_vector_extraction.csv');
if ~isfile(csv_file)
    warning('Plot_1D_Then_All_12:MissingObservations', ...
        'Bulk-pH standard errors not found: %s',csv_file);
    return
end
opts = detectImportOptions(csv_file,'CommentStyle','#','Delimiter',',');
opts = setvartype(opts,{'treatment','depth_cm','period'},'string');
opts = setvartype(opts,{'pH','SEM'},'double');
T = readtable(csv_file,opts);
for k = 1:size(series,1)
    rows = T.treatment == series{k,2} & T.depth_cm == series{k,3};
    for j = 1:numel(labels)
        pick = find(rows & T.period == labels(j));
        if ~isscalar(pick)
            error('Plot_1D_Then_All_12:BadPHObservations', ...
                '%s has %d rows for %s %s cm %s.',csv_file,numel(pick), ...
                series{k,2},series{k,3},labels(j));
        end
        sem.(series{k,1})(j) = T.SEM(pick);
    end
end
end

%% =======================================================================
function obs = clip_observations(obs,Date_d)
%% Keep only the samples inside the simulated window.  The time panels
%% share one x axis, so a point outside it would stretch every panel.
if isempty(obs.date) || isempty(Date_d); return; end
keep = obs.date >= Date_d(1) & obs.date <= Date_d(end);
obs.date = obs.date(keep);
obs.flux = obs.flux(keep);
obs.treatment = obs.treatment(keep);
end

%% =======================================================================
function [model_mean,treatment_mean,n_day] = paired_n2o_means(Date_d,flux,obs)
%% Model and field means restricted to the days the chambers were read.
%% Comparing a growing-season model mean against a chamber mean would put
%% the two on different day sets: the campaign sampled around fertiliser
%% and rainfall events, so its mean sits well above an all-day mean.  The
%% model side therefore uses only the sampled days.  The two field series
%% were not always read on the same dates, so each keeps its own days.
model_mean = nan(size(flux,2),1);
treatment_mean = struct('control',NaN,'basalt',NaN);
n_day = 0;
if isempty(obs.date); return; end
is_sample_day = ismember(Date_d,unique(obs.date));
n_day = sum(is_sample_day);
if n_day > 0
    model_mean = mean(flux(is_sample_day,:),1).';
end
for name = ["control" "basalt"]
    pick = obs.treatment == name;
    if any(pick)
        treatment_mean.(name) = mean(obs.flux(pick));
    end
end
end

%% =======================================================================
function draw_pH_comparison(ax,~,~,obs,colors)
%% Show the observed bulk pH of both sampled depth intervals beneath the
%% simulated time series. The model is one lumped 0-30 cm state, so the two
%% intervals are drawn as measured -- filled 0-10 cm, open 10-30 cm -- and
%% bracket where a lumped value can fall; they are not merged into a depth
%% average.
format_comparison_axis(ax);
if ~obs.available
    text(ax,0.5,0.5,'Bulk-pH observations unavailable', ...
        'Units','normalized','HorizontalAlignment','center','Color','k');
    return
end
labels = obs.labels;
x = (1:numel(labels)).';
hold(ax,'on');
treatments = {'control','basalt'};
depth = struct('field',{'_0','_10'},'tag',{'0_10','10_30'}, ...
    'marker',{'o','s'},'offset',{-0.08,0.08});
for treatment = 1:2
    for d = 1:2
        value = obs.([treatments{treatment} depth(d).field]);
        if d == 1; face = colors(treatment,:); else; face = 'none'; end
        %% SEM first, so the marker sits on top of its bar.
        errorbar(ax,x+depth(d).offset,value, ...
            obs.sem.([treatments{treatment} depth(d).field]), ...
            'LineStyle','none','Color',colors(treatment,:), ...
            'LineWidth',0.75,'CapSize',4,'HandleVisibility','off', ...
            'Tag',sprintf('Observed_%s_pH_%s_SEM',treatments{treatment}, ...
            depth(d).tag));
        scatter(ax,x+depth(d).offset,value,36,colors(treatment,:), ...
            depth(d).marker,'MarkerFaceColor',face, ...
            'MarkerEdgeColor',colors(treatment,:),'LineWidth',0.75, ...
            'Tag',sprintf('Observed_%s_pH_%s',treatments{treatment}, ...
            depth(d).tag));
    end
end
hold(ax,'off');
ax.XLim = [0.5 numel(labels)+0.5];
ax.XTick = x;
ax.XTickLabel = labels;
ylabel(ax,'pH','Color','k');
title(ax,'Observed bulk pH: filled 0-10 cm, open 10-30 cm', ...
    'Color','k','FontSize',8);
end

%% =======================================================================
function overlay_n2o_observations(ax,obs,colors)
%% Chamber N2O drawn onto the modelled series at the dates it was sampled.
%% Markers only, never a connecting line: consecutive campaigns are days to
%% weeks apart, and joining them would draw flux on days nobody measured --
%% exactly the fabricated detail that makes a measurement look like a model.
%% n2o_observations already restricts obs to the maize years (2017, 2018)
%% and clip_observations to the model window, so every point has a modelled
%% day beneath it.
%%
%% Each point takes its treatment's colour and a white edge. The colour is
%% what makes the comparison readable -- an observation sits in the same
%% hue as the traces it tests -- and the edge is what keeps it visible when
%% it lands on top of one of them, which is precisely where the points that
%% matter fall.
if isempty(obs.date); return; end
was_held = ishold(ax);
hold(ax,'on');
treatments = {'control','basalt'};
for treatment = 1:2
    pick = obs.treatment == treatments{treatment};
    if ~any(pick); continue; end
    %% HandleVisibility off keeps these out of the treatment legend, which
    %% is built from the modelled lines and already carries this colour key.
    plot(ax,obs.date(pick),obs.flux(pick),'o', ...
        'MarkerFaceColor',colors(treatment,:),'MarkerEdgeColor','w', ...
        'MarkerSize',5,'LineWidth',0.5,'LineStyle','none', ...
        'HandleVisibility','off', ...
        'Tag',sprintf('Observed_%s_N2O_series',treatments{treatment}));
end
if ~was_held; hold(ax,'off'); end
end

%% =======================================================================
function format_comparison_axis(ax)
ax.Color = 'w';
ax.XColor = 'k';
ax.YColor = 'k';
ax.XTickLabelRotation = 0;
grid(ax,'off');
end

%% =======================================================================
function [lgd_treatment,lgd_configuration,lgd_observation,band_top] = ...
    add_split_legends(fig,treatment_axes,treatment_lines, ...
    configuration_labels,configuration_styles,treatment_labels, ...
    treatment_colors,observation_labels,observation_colors)
%% Two short single-row legends stacked at an explicit figure-normalized
%% position: treatment_axes' own bottom margin must already have been
%% shrunk by the caller to leave room for both rows below y_bottom plus
%% their own height. A second, invisible full-figure axes hosts the
%% configuration legend, since one axes can carry only one legend and
%% treatment_axes already carries the treatment legend.
PAD_BELOW = 0.015;   %% clear space beneath the bottom legend row
GAP_ABOVE = 0.030;   %% clear space between the top legend row and the axes

lgd_treatment = legend(treatment_axes,treatment_lines,treatment_labels, ...
    'Orientation','horizontal','Box','off','TextColor','k');
lgd_treatment.Units = 'normalized';

overlay = axes('Parent',fig,'Position',[0 0 1 1],'Visible','off', ...
    'HitTest','off','HandleVisibility','off');
hold(overlay,'on');
configuration_proxy = gobjects(numel(configuration_labels),1);
for c = 1:numel(configuration_labels)
    configuration_proxy(c) = plot(overlay,NaN,NaN,configuration_styles{c}, ...
        'Color',[0.3 0.3 0.3],'LineWidth',1.35);
end
lgd_configuration = legend(overlay,configuration_proxy,configuration_labels, ...
    'Orientation','horizontal','Box','off','TextColor','k');
lgd_configuration.Units = 'normalized';

%% Third row, only where measurements are overlaid: the marker proxies
%% carry the same face colour and white edge as the plotted points, so the
%% key matches what is actually on the axes rather than approximating it.
lgd_observation = gobjects(1);
if nargin >= 9 && ~isempty(observation_labels)
    overlay_obs = axes('Parent',fig,'Position',[0 0 1 1],'Visible','off', ...
        'HitTest','off','HandleVisibility','off');
    hold(overlay_obs,'on');
    observation_proxy = gobjects(numel(observation_labels),1);
    for c = 1:numel(observation_labels)
        observation_proxy(c) = plot(overlay_obs,NaN,NaN,'o', ...
            'MarkerFaceColor',observation_colors(c,:),'MarkerEdgeColor','w', ...
            'MarkerSize',5,'LineWidth',0.5,'LineStyle','none');
    end
    lgd_observation = legend(overlay_obs,observation_proxy,observation_labels, ...
        'Orientation','horizontal','Box','off','TextColor','k');
    lgd_observation.Units = 'normalized';
end

%% One drawnow for every legend, then stack them on a common row height so
%% the rows sit evenly regardless of which happens to render tallest.
drawnow;
row_h = max([lgd_treatment.Position(4) lgd_configuration.Position(4)]);
n_rows = 2;
if isgraphics(lgd_observation)
    row_h = max(row_h,lgd_observation.Position(4));
    n_rows = 3;
end
place_legend_row(lgd_treatment,PAD_BELOW,0.26);
place_legend_row(lgd_configuration,PAD_BELOW+row_h,0.60);
if isgraphics(lgd_observation)
    place_legend_row(lgd_observation,PAD_BELOW+2*row_h,0.40);
end
band_top = PAD_BELOW+n_rows*row_h+GAP_ABOVE;
end

%% =======================================================================
function layout_time_panels(panels,band_top)
%% Stack the panels between the legend band and the top of the figure,
%% filling the space instead of leaving subplot's default margins. Panels
%% arrive top-to-bottom.
%%
%% The gap below each panel is sized to what that panel actually renders:
%% a panel whose tick labels have been blanked (the upper pH panels, which
%% hand their date axis to the one below) needs only enough room for the
%% next panel's title, while one still carrying rotated date labels and an
%% axis label needs roughly half again as much. Using a single worst-case
%% gap everywhere is what leaves the empty strips.
LEFT = 0.10; WIDTH = 0.87;
n = numel(panels);
%% An axes' Position is its plot box only. Tick labels and the axis label
%% render BELOW it and the title ABOVE it, so a layout built from Position
%% alone prints those decorations across whatever sits next to the panel --
%% the "Date" label over the legend rows, or over the title of the panel
%% beneath. TightInset gives what those decorations actually occupy, in the
%% same normalized units, so the spacing is measured rather than guessed.
%%
%% Two passes: the first uses rough spacing only to get the panels to a
%% representative size, the second reads each panel's insets and sets the
%% real layout. A single constant cannot serve, because within one figure
%% the panels differ -- the pH figure's middle panel carries rotated dates
%% and a "Date" label while the one above it has both suppressed.
gaps = repmat(0.085,1,max(0,n-1));
bottom_margin = 0.06;
top = 0.945;
for pass = 1:2
    h = (top-band_top-bottom_margin-sum(gaps))/n;
    if h <= 0; return; end   %% pathological figure size: keep subplot's layout
    y = top;
    for k = 1:n
        y = y-h;
        panels(k).Position = [LEFT y WIDTH h];
        if k <= numel(gaps); y = y-gaps(k); end
    end
    if pass == 1
        drawnow;
        for k = 1:n-1
            %% clear the upper panel's labels and the lower panel's title
            gaps(k) = panels(k).TightInset(2)+panels(k+1).TightInset(4)+0.008;
        end
        bottom_margin = panels(end).TightInset(2)+0.012;
        top = 0.995-panels(1).TightInset(4);
    end
end
end

%% =======================================================================
function row_h = place_legend_row(lgd,y,min_width)
%% Centre one horizontal legend row at an explicit figure-normalized height
%% and return its height so the caller can stack the next row above it.
%%
%% The width is floored rather than taken as given. MATLAB's auto-computed
%% legend width, read immediately after creation, is occasionally stale --
%% narrower than its own entries need -- on Bouchet's MATLAB/Linux stack;
%% freezing the position then locks that undersized box and the entries lay
%% out overlapping inside it. This is why the four-entry configuration row
%% rendered correctly on the pH figure, which does enough drawing beforehand
%% for the width to settle, and illegibly on the single-panel figures.
%%
%% In practice the floor binds rather than merely rescuing: the auto width
%% comes back below it locally too, so these values, not MATLAB's, set the
%% row width in the normal case. They are sized from a correctly rendered
%% Bouchet figure, where the configuration row spanned about 0.64 of the
%% figure. That asymmetry is deliberate -- a row wider than its entries need
%% only spreads them out, whereas one too narrow is what makes them
%% unreadable -- but it does mean changing the label text warrants
%% re-checking these numbers rather than trusting them to adapt.
drawnow;
w = max(lgd.Position(3),min_width);
lgd.Position = [0.5-w/2, y, w, lgd.Position(4)];
drawnow;
row_h = lgd.Position(4);
end

%% =======================================================================
function format_time_axis(ax)
ax.Color = 'w';
ax.XColor = 'k';
ax.YColor = 'k';
ax.XTickLabelRotation = 20;
grid(ax,'off');
end

%% =======================================================================
function format_ph_axis(ax,value)
limits = [min(value,[],'all') max(value,[],'all')];
span = diff(limits);
if span < 0.02
    centre = mean(limits);
    limits = centre+[-0.05 0.05];
else
    limits = limits+0.05*span*[-1 1];
end
ylim(ax,limits);
ytickformat(ax,'%.2f');
end

%% =======================================================================
function management = kantola_management(Date_d)
maize_years = [2017 2018 2020];
soybean_years = [2016 2019];
fertilizer_doy = [81 100];

year_start = datetime(maize_years,ones(size(maize_years)), ...
    ones(size(maize_years)));
fertilizer_start = year_start+caldays(fertilizer_doy(1)-1);
fertilizer_end = year_start+caldays(fertilizer_doy(2)-1);
keep = fertilizer_end >= Date_d(1) & fertilizer_start <= Date_d(end);
management.fertilizer_start = fertilizer_start(keep);
management.fertilizer_end = fertilizer_end(keep);

crop_years = [soybean_years maize_years];
crop_names = [repmat({'Soybean'},size(soybean_years)), ...
    repmat({'Maize'},size(maize_years))];
[crop_years,order] = sort(crop_years);
crop_names = crop_names(order);
crop_start = datetime(crop_years,ones(size(crop_years)), ...
    ones(size(crop_years)));
crop_end = datetime(crop_years,12*ones(size(crop_years)), ...
    31*ones(size(crop_years)));
keep = crop_end >= Date_d(1) & crop_start <= Date_d(end);
management.crop_years = crop_years(keep);
management.crop_names = crop_names(keep);
management.crop_start = crop_start(keep);
management.crop_end = crop_end(keep);
management.date_start = Date_d(1);
management.date_end = Date_d(end);
end

%% =======================================================================
function fertilizer_legend = draw_kantola_management(ax,management)
fertilizer_legend = gobjects(1);
fertilizer_color = [0.75 0.30 0.10];
if ~isempty(management.fertilizer_start)
    marker = xline(ax,management.fertilizer_start,'--', ...
        'Color',fertilizer_color,'LineWidth',1.0, ...
        'Tag','KantolaFertilizerStart');
    stop_marker = xline(ax,management.fertilizer_end,':', ...
        'Color',fertilizer_color,'LineWidth',1.0, ...
        'Tag','KantolaFertilizerEnd');
    set(stop_marker,'HandleVisibility','off');
    fertilizer_legend = marker(1);
end

data_ylim = ylim(ax);
data_span = diff(data_ylim);
if all(isfinite(data_ylim)) && data_span > 0
    ylim(ax,[data_ylim(1) data_ylim(2)+0.12*data_span]);
end
date_span = days(management.date_end-management.date_start);
for k = 1:numel(management.crop_years)
    left = management.crop_start(k);
    right = management.crop_end(k);
    if left < management.date_start; left = management.date_start; end
    if right > management.date_end; right = management.date_end; end
    if date_span > 0
        x = days(left+(right-left)/2-management.date_start)/date_span;
    else
        x = 0.5;
    end
    text(ax,x,0.97,sprintf('%d %s',management.crop_years(k), ...
        management.crop_names{k}),'Units','normalized', ...
        'HorizontalAlignment','center','VerticalAlignment','top', ...
        'Color',[0.25 0.25 0.25],'FontSize',8, ...
        'BackgroundColor','w','Margin',1, ...
        'Tag','KantolaCropLabel');
end
end
