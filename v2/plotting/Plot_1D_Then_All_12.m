%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot the vertical 1-D group and the full twelve-scenario set in one
%% run.  Standalone runnable script: no function signature, no arguments.
%%
%% Default behaviour
%%   1) The vertical 1-D weathering and nitrogen figure (runs 10:12).
%%      Base cations, cumulative denitrification and WFPS are panels that
%%      the summary figure does not carry.  Prints the years-2-onward
%%      nitrogen statistics.
%%   2) The twelve-scenario summary and time dynamics (runs 1:12),
%%      including the associated FIG files.
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
%%   Visible          figure visibility                ['on']
%%   OutputFile       summary PNG (single-group mode)  []
%%   TimeOutputFile   dynamics PNG (single-group mode) []
%%   DryRun           print the plan without plotting  [false]
%%
%% Examples
%%   matlab -batch "cd('v2/plotting'); Plot_1D_Then_All_12"
%%   matlab -batch "cd('v2/plotting'); NDays=1; Plot_1D_Then_All_12"
%%   matlab -batch "cd('v2/plotting'); Runs=10:12; Plot_1D_Then_All_12"
%%
%% Result variables are left in the workspace: summary, fig,
%% time_fig, plot_file, time_plot_file for the summary stage, and
%% fig1d, plot_file_1d, stats_1d for the 1-D weathering stage.
%%
%% NEE follows the TC sign convention: positive is carbon released to the
%% atmosphere.  Yield uses the Kantola low crop PFTs only.  Both are
%% normalised by the full domain area.
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
if ~exist('Visible','var') || isempty(Visible); Visible = 'on'; end
if ~exist('OutputFile','var') || isempty(OutputFile); OutputFile = ''; end
if ~exist('TimeOutputFile','var') || isempty(TimeOutputFile); TimeOutputFile = ''; end
if ~exist('DryRun','var') || isempty(DryRun); DryRun = false; end

validate_options(NDays,Visible,DryRun,ResultDirectory);
if ~isfolder(FigureDirectory); mkdir(FigureDirectory); end
opt = struct('ResultDirectory',ResultDirectory, ...
    'FigureDirectory',FigureDirectory,'NDays',NDays, ...
    'Visible',Visible,'OutputFile',OutputFile, ...
    'TimeOutputFile',TimeOutputFile);

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
function validate_options(NDays,Visible,DryRun,ResultDirectory)
if ~isempty(NDays) && ~(isnumeric(NDays) && isscalar(NDays) && ...
        isfinite(NDays) && NDays == fix(NDays) && NDays >= 1)
    error('Plot_1D_Then_All_12:BadOptions','NDays must be a positive integer.');
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
%% denitrification (N2 + N2O), and water-filled pore space live here.
%% The denitrification N2O fraction rN2O is computed and reported
%% flux-weighted in the printed statistics, but is not plotted: as a daily
%% curve it is just the partition function re-evaluated on that day's pH,
%% which reads as a larger treatment effect than the run actually produces.
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

suffix = '';
if ~isempty(opt.NDays); suffix = sprintf('__%dd',opt.NDays); end

Date_d = datetime.empty(0,1);
for i = 1:3
    id = run_ids(i);
    result_file = fullfile(opt.ResultDirectory, ...
        sprintf('TC_2D_run%02d_%s%s.mat',id,RM(id).name,suffix));
    if ~isfile(result_file)
        error('Plot_1D_Then_All_12:MissingResult','Missing result: %s',result_file);
    end
    S = load(result_file);
    P = S.P; if ndims(P) == 3; P = P(:,:,1); end
    nd = size(P,1);
    if i == 1
        Date_d = datetime(S.Datam(1,1),S.Datam(1,2),S.Datam(1,3)) + ...
            caldays((0:nd-1).');
        ca = zeros(nd,3); mg = zeros(nd,3); si = zeros(nd,3);
        den = zeros(nd,3); n2o = zeros(nd,3); wfps = zeros(nd,3);
        ph = zeros(nd,3); n2 = zeros(nd,3);
    end
    ca(:,i) = P(:,iCa); mg(:,i) = P(:,iMg); si(:,i) = P(:,iSi);
    n2o(:,i) = daily(S,'N2Oflx',nd);
    n2(:,i)  = daily(S,'N2flx',nd);
    den(:,i) = n2(:,i) + n2o(:,i);
    wfps(:,i) = water_filled_pore_space(S,nd);
    ph(:,i) = daily(S,'PH_d',nd);
end

%% Denitrification N2O fraction, evaluated from the model's own partition.
%%
%% Do not compute this as N2Oflx/(N2Oflx+N2flx) from the saved fluxes. That
%% ratio is not the denitrification partition, because N2Oflx carries a
%% second source:
%%
%%   N2Oflx = f_N2O_nit*NO3flx  +  rN2O*DENITflx
%%
%% The first term is nitrification N2O and has nothing to do with the N2O
%% against N2 split.  Reconstruct rN2O directly from daily pH using the
%% same expression BIOGEOCHEMISTRY_DYNAMIC3 uses, defined on every day.
BiogeoPar = Biogeochemistry_Parameter(1);
rN2O = 100./(1+exp(BiogeoPar.k_rN2O*(ph-BiogeoPar.pH0_rN2O)));
rN2O(ph <= 0) = NaN;

%% Denitrification-weighted rN2O over the equilibrium years, used to label
%% panel 5.  What the nitrogen actually experiences is the daily share
%% weighted by the denitrification flux, which differs from the plain mean
%% because denitrification concentrates on wet days that do not carry the
%% mean pH.
nd_all = size(den,1);
eq = min(nd_all,round(365.25)+1):nd_all;
rN2O_fw = zeros(1,3);
for i = 1:3
    rN2O_fw(i) = rN2O_fluxweighted(rN2O(eq,i),n2(eq,i));
end

time_value = {ca,mg,si,cumsum(den,1),wfps};
time_label = {'Ca solution','Mg solution','Si solution', ...
    'Cumulative denitrification (N_2 + N_2O)','Water-filled pore space'};
time_units = {'\mu mol','\mu mol','\mu mol','g N m^{-2}','%'};
time_tag = {'CaSolution','MgSolution','SiSolution', ...
    'CumulativeDenitrification','WFPS'};

application_dates = datetime([2016 11 10; 2017 11 10; ...
    2018 11 10; 2019 11 10]);
application_dates = application_dates( ...
    application_dates >= Date_d(1) & application_dates <= Date_d(end));

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
for measurement = 1:5
    ax = nexttile(tl,tile_span{measurement});
    time_axes(measurement) = ax;
    ax.Tag = sprintf('TimeMetric%d',measurement);
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
    format_time_axis(ax);
    ylabel(ax,time_units{measurement},'Color','k');
    if measurement >= 4; xlabel(ax,'Date','Color','k'); end
    text(ax,0.02,0.96,time_label{measurement},'Units','normalized', ...
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
lgd = legend(legend_handles,legend_labels,'Orientation','horizontal', ...
    'Box','off','TextColor','k','Color','w','EdgeColor','none');
lgd.Layout.Tile = 'south';
title(tl,sprintf('%s: weathering signal and nitrogen partitioning', ...
    group_label),'Color','k');

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
nd = size(den,1);
sl = min(nd,round(365.25)+1):nd;
f = 10/(numel(sl)/365.25);          %% gN/m2 over the window -> kg N/ha/yr
stats = struct();
fprintf('\nyears 2 onward, kg N/ha/yr\n');
fprintf('%-8s %7s %8s %10s %8s %9s\n', ...
    'run','N2O','denit','rN2O','WFPS','Ca end');
for i = 1:3
    st = struct('treatment',treatment_labels{i}, ...
        'N2O',sum(n2o(sl,i),'omitnan')*f, ...
        'denit',sum(den(sl,i),'omitnan')*f, ...
        'rN2O',rN2O_fw(i), ...
        'wfps',mean(wfps(sl,i),'omitnan'), ...
        'Ca_end',ca(end,i),'Mg_end',mg(end,i),'Si_end',si(end,i));
    stats.(lower(treatment_labels{i})) = st;
    fprintf('%-8s %7.1f %8.1f %9.1f%% %7.1f%% %9.2f\n', ...
        st.treatment,st.N2O,st.denit,st.rN2O,st.wfps,st.Ca_end);
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
    tag = RM(i).name;
    if ~isempty(opt.NDays); tag = sprintf('%s__%dd',tag,opt.NDays); end
    result_file = fullfile(opt.ResultDirectory, ...
        sprintf('TC_2D_run%02d_%s.mat',run_id,tag));
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
    required = {'N2Oflx','LEAK_NO3','NO3_pool','PH_d'};
    if ~all(isfield(S.Domain,required)) || ~isfield(S,'Datam') || ...
            isempty(S.Domain.NO3_pool) || isempty(S.Domain.PH_d)
        error('Plot_1D_Then_All_12:IncompleteResult', ...
            'Domain summary is incomplete in %s.',result_file);
    end
    measurement = {S.Domain.N2Oflx,S.Domain.LEAK_NO3, ...
        S.Domain.NO3_pool,S.Domain.PH_d};
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
        caldays((0:n_time(1)-1).');
    [run_nee,run_crop_harvest] = carbon_dynamics( ...
        S,n_time(1),result_file);
    run_has_bulk_pore_ph = isfield(S.Domain,'PH_bulk_d') && isfield(S.Domain,'PH_pore_d');
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
        ph_bulk_daily(:,i) = S.Domain.PH_bulk_d(:); %#ok<AGROW>
        ph_pore_daily(:,i) = S.Domain.PH_pore_d(:); %#ok<AGROW>
    else
        ph_bulk_daily(:,i) = S.Domain.PH_d(:); %#ok<AGROW>
        ph_pore_daily(:,i) = S.Domain.PH_d(:); %#ok<AGROW>
    end
    nee_daily(:,i) = run_nee; %#ok<AGROW>
    crop_harvest_daily(:,i) = run_crop_harvest; %#ok<AGROW>
    n2o(i) = sum(n2o_daily(:,i));
    no3_leach(i) = sum(no3_leach_daily(:,i));
    no3_final(i) = no3_pool_daily(end,i);
    ph_final(i) = ph_pore_daily(end,i);
    nee(i) = sum(nee_daily(:,i));
    crop_harvest(i) = sum(crop_harvest_daily(:,i));
end

summary = table([RM.id].',string({RM.layout}).',string({RM.treatment}).', ...
    n2o,no3_leach,no3_final,ph_final,nee,crop_harvest, ...
    'VariableNames',{'RunID','ScenarioGroup','Treatment', ...
    'CumulativeN2O_gN_m2','CumulativeNO3Leaching_gN_m2', ...
    'FinalNO3Pool_gN_m2','FinalPoreWaterPH', ...
    'CumulativeNEE_gC_m2','CumulativeCropYield_gC_m2'});

all_group_labels = {'Whole field','Margin with crop', ...
    'Margin with prairie','Vertical 1-D'};
all_group_stems = {'Whole_Field','Margin_Crop','Margin_Prairie','Vertical_1D'};
group_ids = unique(ceil(run_ids/3),'stable');
group_labels = all_group_labels(group_ids);
n_group = numel(group_labels);
treatment_labels = {'Control','Basalt','Calcite'};
treatment_colors = [0 0 0; 0 0.25 0.60; 0.35 0.70 0.90];
if isequal(run_ids,1:12)
    output_stem = 'Subfield_12';
    summary_name = 'T&C-SMEW-N twelve-scenario summary';
    time_name = 'T&C-SMEW-N twelve-scenario time dynamics';
    summary_title = 'Treatment response across four scenario groups';
    time_title = 'Treatment dynamics across four scenario groups';
elseif isequal(run_ids,1:9)
    output_stem = 'Subfield_2D';
    summary_name = 'T&C-SMEW-N 2-D scenario summary';
    time_name = 'T&C-SMEW-N 2-D scenario time dynamics';
    summary_title = 'Treatment response across three 2-D scenario groups';
    time_title = 'Treatment dynamics across three 2-D scenario groups';
elseif n_group == 1
    output_stem = all_group_stems{group_ids};
    summary_name = sprintf('T&C-SMEW-N %s summary',group_labels{1});
    time_name = sprintf('T&C-SMEW-N %s time dynamics',group_labels{1});
    summary_title = sprintf('Treatment response: %s',group_labels{1});
    time_title = sprintf('Treatment dynamics: %s',group_labels{1});
else
    run_tag = sprintf('%02d_',run_ids);
    output_stem = ['Subfield_runs_' run_tag(1:end-1)];
    summary_name = 'T&C-SMEW-N selected-scenario summary';
    time_name = 'T&C-SMEW-N selected-scenario time dynamics';
    summary_title = 'Treatment response across selected scenario groups';
    time_title = 'Treatment dynamics across selected scenario groups';
end
fig = figure('Color','w','Visible',char(opt.Visible), ...
    'Name',summary_name);
fig.Position(3:4) = [1400 1050];
tl = tiledlayout(fig,3,2,'TileSpacing','compact','Padding','compact');
summary_axes = gobjects(6,1);
summary_axes(1) = nexttile(tl);
b = draw_panel(summary_axes(1),n2o,'Cumulative N_2O','g N m^{-2}', ...
    group_labels,treatment_colors);
summary_axes(2) = nexttile(tl);
draw_panel(summary_axes(2),no3_leach,'Cumulative NO_3 leaching','g N m^{-2}', ...
    group_labels,treatment_colors);
summary_axes(3) = nexttile(tl);
draw_panel(summary_axes(3),no3_final,'Final soil NO_3 pool','g N m^{-2}', ...
    group_labels,treatment_colors);
summary_axes(4) = nexttile(tl);
draw_panel(summary_axes(4),ph_final,'Final pore-water pH','pH', ...
    group_labels,treatment_colors);
summary_axes(5) = nexttile(tl);
draw_panel(summary_axes(5),nee,'Cumulative NEE (+ to atmosphere)', ...
    'g C m^{-2} domain', ...
    group_labels,treatment_colors);
summary_axes(6) = nexttile(tl);
draw_panel(summary_axes(6),crop_harvest,'Cumulative crop yield', ...
    'g C m^{-2} domain',group_labels,treatment_colors);
for iax = 1:6
    summary_axes(iax).Tag = sprintf('SummaryRow%dColumn%d', ...
        ceil(iax/2),mod(iax-1,2)+1);
end
set(summary_axes(1:4),'XTickLabel',[]);
lgd = legend(b,treatment_labels,'Orientation','horizontal','Box','off', ...
    'TextColor','k','Color','w','EdgeColor','none');
lgd.Layout.Tile = 'south';
title(tl,summary_title,'Color','k');

suffix = '';
if ~isempty(opt.NDays); suffix = sprintf('__%dd',opt.NDays); end
if isempty(opt.OutputFile)
    plot_file = fullfile(opt.FigureDirectory, ...
        [output_stem '_summary' suffix '.png']);
else
    plot_file = char(opt.OutputFile);
end
exportgraphics(fig,plot_file,'Resolution',200);
[plot_folder,plot_name] = fileparts(plot_file);
fig_file = fullfile(plot_folder,[plot_name '.fig']);
savefig(fig,fig_file);
fprintf('Saved: %s\n',plot_file);
fprintf('Saved: %s\n',fig_file);

time_fig = figure('Color','w','Visible',char(opt.Visible), ...
    'Name',time_name);
time_fig.Position(3:4) = [max(1200,700*n_group) 1050];
n_time_columns = 2*n_group;
time_tl = tiledlayout(time_fig,3,n_time_columns, ...
    'TileSpacing','compact','Padding','compact');
time_value = {cumsum(n2o_daily,1),cumsum(no3_leach_daily,1), ...
    no3_pool_daily,ph_pore_daily,nee_daily, ...
    cumsum(crop_harvest_daily,1)};
time_label = {'Cumulative N_2O','Cumulative NO_3 leaching', ...
    'Soil NO_3 pool','Pore-water pH','NEE (+ to atmosphere)', ...
    'Cumulative crop yield'};
time_units = {'g N m^{-2}','g N m^{-2}','g N m^{-2}','pH', ...
    'g C m^{-2} domain d^{-1}','g C m^{-2} domain'};
time_tag = {'CumulativeN2O','CumulativeNO3Leaching','NO3Pool', ...
    'PoreWaterPH','NEE','CumulativeCropYield'};
application_dates = datetime([2016 11 10; 2017 11 10; ...
    2018 11 10; 2019 11 10]);
application_dates = application_dates( ...
    application_dates >= Date_d(1) & application_dates <= Date_d(end));
application_color = [0.45 0.45 0.45];
time_axes = gobjects(6,n_group);
legend_lines = gobjects(3,1);
application_legend = gobjects(1);
for plot_row = 1:3
    for group = 1:n_group
        for pair_column = 1:2
            measurement = 2*(plot_row-1)+pair_column;
            tile_column = 2*(group-1)+pair_column;
            tile_id = (plot_row-1)*n_time_columns+tile_column;
            ax = nexttile(time_tl,tile_id);
            time_axes(measurement,group) = ax;
            ax.Tag = sprintf('TimeMetric%dGroup%d',measurement,group);
            ids = (group-1)*3+(1:3);
            ln = plot(ax,Date_d,time_value{measurement}(:,ids),'LineWidth',1.15);
            for treatment = 1:3
                ln(treatment).Color = treatment_colors(treatment,:);
            end
            set(ln,'Tag',time_tag{measurement});
            if measurement == 4 && has_bulk_pore_ph
                ln_bulk = plot(ax,Date_d,ph_bulk_daily(:,ids), ...
                    'LineStyle','--','LineWidth',1.15);
                for treatment = 1:3
                    ln_bulk(treatment).Color = treatment_colors(treatment,:);
                end
                set(ln_bulk,'Tag','BulkPH');
            end
            if ~isempty(application_dates)
                marker = xline(ax,application_dates,':','Color',application_color, ...
                    'LineWidth',1.1,'Tag','FeedstockApplication');
                if measurement == 1 && group == 1
                    application_legend = marker(1);
                end
            end
            format_time_axis(ax);
            if plot_row == 1
                title(ax,group_labels{group},'Color','k');
            end
            if group == 1
                ylabel(ax,time_units{measurement},'Color','k');
            end
            if plot_row == 3; xlabel(ax,'Date','Color','k'); end
            text(ax,0.02,0.96,time_label{measurement},'Units','normalized', ...
                'VerticalAlignment','top','Color','k','FontWeight','bold');
            if measurement == 1 && group == 1; legend_lines = ln; end
        end
    end
end
linkaxes(time_axes(:),'x');
set(time_axes(1:4,:),'XTickLabel',[]);
for measurement = 1:6
    measurement_axes = time_axes(measurement,:);
    measurement_lim = vertcat(measurement_axes.YLim);
    set(measurement_axes,'YLim', ...
        [min(measurement_lim(:,1)) max(measurement_lim(:,2))]);
end
legend_handles = legend_lines;
legend_labels = treatment_labels;
if isgraphics(application_legend)
    legend_handles(end+1) = application_legend;
    legend_labels{end+1} = 'Feedstock application';
end
time_lgd = legend(legend_handles,legend_labels,'Orientation','horizontal', ...
    'Box','off','TextColor','k','Color','w','EdgeColor','none');
time_lgd.Layout.Tile = 'south';
title(time_tl,time_title,'Color','k');

if isempty(opt.TimeOutputFile)
    time_plot_file = fullfile(opt.FigureDirectory, ...
        [output_stem '_time_dynamics' suffix '.png']);
else
    time_plot_file = char(opt.TimeOutputFile);
end
exportgraphics(time_fig,time_plot_file,'Resolution',200);
[time_plot_folder,time_plot_name] = fileparts(time_plot_file);
time_fig_file = fullfile(time_plot_folder,[time_plot_name '.fig']);
savefig(time_fig,time_fig_file);
fprintf('Saved: %s\n',time_plot_file);
fprintf('Saved: %s\n',time_fig_file);

end

%% =======================================================================
function missing = missing_results(run_ids,opt)
missing = zeros(0,1);
RM = Run_Matrix_12();
for id = run_ids
    tag = RM(id).name;
    if ~isempty(opt.NDays); tag = sprintf('%s__%dd',tag,opt.NDays); end
    result_file = fullfile(opt.ResultDirectory, ...
        sprintf('TC_2D_run%02d_%s.mat',id,tag));
    if ~isfile(result_file); missing(end+1,1) = id; end %#ok<AGROW>
end
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
function r = rN2O_fluxweighted(share_pct,n2_flux)
%% Denitrification-weighted mean of the N2O fraction [%].
%%
%% The plain mean over days is not what the nitrogen actually experiences.
%% Denitrification is concentrated on wet days, and those days do not carry
%% the mean pH, so weighting by the flux changes the answer materially.
%% Recover total denitrification from N2flx = (1-rN2O)*DENITflx.
f = share_pct/100;
den = n2_flux./max(1-f,eps);
r = 100*sum(f.*den,'omitnan')/max(sum(den,'omitnan'),eps);
end

%% =======================================================================
function v = daily(S,name,nd)
%% Daily flux as a column, or zeros when the run did not save it.
if isfield(S,name); v = S.(name); v = v(:); v = v(1:nd); else; v = zeros(nd,1); end
end

%% =======================================================================
function w = water_filled_pore_space(S,nd)
%% Depth-weighted WFPS of the biogeochemical layer [%], daily mean.
p = S.p;
O = S.O; if ndims(O) == 3; O = O(:,:,1); end
sph = size(O,1)/nd;
Osat = sum(p.Bio_Zs.*p.Soil_Param.Osat);
w = zeros(nd,1);
for j = 1:nd
    rows = round((j-1)*sph)+1:round(j*sph);
    w(j) = sum(p.Bio_Zs.*mean(O(rows,:),1))/Osat*100;
end
end

%% =======================================================================
function b = draw_panel(ax,value,panel_title,y_label,group_labels,colors)
n_group = numel(group_labels);
if n_group == 1
    b = gobjects(3,1);
    hold(ax,'on');
    for treatment = 1:3
        b(treatment) = bar(ax,treatment,value(treatment),0.8, ...
            'FaceColor',colors(treatment,:));
    end
    hold(ax,'off');
else
    b = bar(ax,reshape(value,3,n_group).','grouped');
    for treatment = 1:3
        b(treatment).FaceColor = colors(treatment,:);
    end
end
ax.Color = 'w';
ax.XColor = 'k';
ax.YColor = 'k';
if n_group == 1
    ax.XTick = 2;
else
    ax.XTick = 1:n_group;
end
ax.XTickLabel = group_labels;
ax.XTickLabelRotation = 18;
grid(ax,'off');
title(ax,panel_title,'Color','k');
ylabel(ax,y_label,'Color','k');
end

%% =======================================================================
function format_time_axis(ax)
ax.Color = 'w';
ax.XColor = 'k';
ax.YColor = 'k';
ax.XTickLabelRotation = 20;
grid(ax,'off');
end