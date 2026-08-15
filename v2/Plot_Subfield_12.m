function [summary,fig,plot_file,time_fig,time_plot_file] = Plot_Subfield_12(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot summary values and daily dynamics for complete treatment groups.
%% NEE follows the TC sign convention: positive is carbon released to the
%% atmosphere. Yield uses the Kantola low crop PFTs only. Both are
%% normalised by the full domain area.
%%
%% Options
%%   ResultDirectory  result folder              [SIM/Subfield_12]
%%   Runs             complete treatment groups  [1:12]
%%   NDays            shortened-run filename tag []
%%   Visible          figure visibility          ['on']
%%   OutputFile       output PNG and companion FIG [summary in result folder]
%%   TimeOutputFile   output PNG and companion FIG [time dynamics in result folder]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

here = fileparts(mfilename('fullpath'));
opt = struct('ResultDirectory',fullfile(here,'SIM','Subfield_12'), ...
    'Runs',1:12,'NDays',[],'Visible','on','OutputFile','', ...
    'TimeOutputFile','');
if mod(numel(varargin),2) ~= 0
    error('Plot_Subfield_12:BadOptions','Options must be name/value pairs.');
end
for k = 1:2:numel(varargin)
    name = char(varargin{k});
    if ~isfield(opt,name)
        error('Plot_Subfield_12:BadOptions','Unknown option %s.',name);
    end
    opt.(name) = varargin{k+1};
end
if ~isempty(opt.NDays) && ~(isnumeric(opt.NDays) && isscalar(opt.NDays) && ...
        isfinite(opt.NDays) && opt.NDays == fix(opt.NDays) && opt.NDays >= 1)
    error('Plot_Subfield_12:BadOptions','NDays must be a positive integer.');
end
if ~ismember(char(opt.Visible),{'on','off'})
    error('Plot_Subfield_12:BadOptions','Visible must be on or off.');
end
run_ids = reshape(opt.Runs,1,[]);
if ~(isnumeric(run_ids) && isreal(run_ids) && ~isempty(run_ids) && ...
        all(isfinite(run_ids)) && all(run_ids == fix(run_ids)) && ...
        all(run_ids >= 1 & run_ids <= 12))
    error('Plot_Subfield_12:BadRuns', ...
        'Runs must contain complete treatment groups from 1 to 12.');
end
group_ids = unique(ceil(run_ids/3),'stable');
expected_ids = reshape(3*group_ids+(-2:0).',1,[]);
if ~isequal(run_ids,expected_ids)
    error('Plot_Subfield_12:BadRuns', ...
        'Runs must contain complete treatment groups in matrix order.');
end
if ~isfolder(opt.ResultDirectory)
    error('Plot_Subfield_12:MissingDirectory', ...
        'Result directory not found: %s',opt.ResultDirectory);
end

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
ph_daily = [];
nee_daily = [];
crop_harvest_daily = [];
carbon_fields = {'NPP_H','NPP_L','Ccrown_t','R_litter','R_microbe','R_ew'};
result_fields = {'R','Domain','Datam','RB_L','topo','ColCfg'};
for i = 1:n
    run_id = RM(i).id;
    tag = RM(i).name;
    if ~isempty(opt.NDays); tag = sprintf('%s__%dd',tag,opt.NDays); end
    result_file = fullfile(opt.ResultDirectory, ...
        sprintf('TC_2D_run%02d_%s.mat',run_id,tag));
    if ~isfile(result_file)
        error('Plot_Subfield_12:MissingResult','Result not found: %s',result_file);
    end
    saved = whos('-file',result_file);
    saved_names = {saved.name};
    if ~all(ismember(carbon_fields,saved_names))
        error('Plot_Subfield_12:IncompleteCarbonResult', ...
            'Carbon state is incomplete in %s.',result_file);
    end
    if ~all(ismember(result_fields,saved_names))
        error('Plot_Subfield_12:IncompleteResult', ...
            'Required result state is incomplete in %s.',result_file);
    end
    S = load(result_file,result_fields{:},carbon_fields{:});
    if ~isfield(S,'R') || ~isfield(S,'Domain') || ...
            S.R.id ~= run_id || ~strcmp(S.R.name,RM(i).name)
        error('Plot_Subfield_12:WrongResult', ...
            'Result does not match run %d: %s',run_id,result_file);
    end
    required = {'N2Oflx','LEAK_NO3','NO3_pool','PH_d'};
    if ~all(isfield(S.Domain,required)) || ~isfield(S,'Datam') || ...
            isempty(S.Domain.NO3_pool) || isempty(S.Domain.PH_d)
        error('Plot_Subfield_12:IncompleteResult', ...
            'Domain summary is incomplete in %s.',result_file);
    end
    measurement = {S.Domain.N2Oflx,S.Domain.LEAK_NO3, ...
        S.Domain.NO3_pool,S.Domain.PH_d};
    if ~all(cellfun(@(x) isnumeric(x) && isvector(x),measurement))
        error('Plot_Subfield_12:IncompleteResult', ...
            'Domain measurements must be numeric vectors in %s.',result_file);
    end
    n_time = cellfun(@numel,measurement);
    if any(n_time ~= n_time(1))
        error('Plot_Subfield_12:TimeLengthMismatch', ...
            'Domain measurements have different daily lengths in %s.',result_file);
    end
    if any(cellfun(@(x) any(~isfinite(x(:))),measurement))
        error('Plot_Subfield_12:NonfiniteMeasurement', ...
            'Run %d does not contain finite values for every plotted measurement.',run_id);
    end
    if ~isnumeric(S.Datam) || isempty(S.Datam) || size(S.Datam,2) < 3 || ...
            any(~isfinite(S.Datam(1,1:3)))
        error('Plot_Subfield_12:IncompleteResult', ...
            'Daily start date is missing or invalid in %s.',result_file);
    end
    run_date = datetime(S.Datam(1,1),S.Datam(1,2),S.Datam(1,3)) + ...
        caldays((0:n_time(1)-1).');
    [run_nee,run_crop_harvest] = carbon_dynamics( ...
        S,n_time(1),result_file);
    if i == 1
        Date_d = run_date;
        n2o_daily = zeros(n_time(1),n);
        no3_leach_daily = zeros(n_time(1),n);
        no3_pool_daily = zeros(n_time(1),n);
        ph_daily = zeros(n_time(1),n);
        nee_daily = zeros(n_time(1),n);
        crop_harvest_daily = zeros(n_time(1),n);
    elseif ~isequal(run_date,Date_d)
        error('Plot_Subfield_12:TimeAxisMismatch', ...
            'Run %d does not use the same daily time axis as run %d.', ...
            run_id,RM(1).id);
    end
    n2o_daily(:,i) = S.Domain.N2Oflx(:); %#ok<AGROW>
    no3_leach_daily(:,i) = S.Domain.LEAK_NO3(:); %#ok<AGROW>
    no3_pool_daily(:,i) = S.Domain.NO3_pool(:); %#ok<AGROW>
    ph_daily(:,i) = S.Domain.PH_d(:); %#ok<AGROW>
    nee_daily(:,i) = run_nee; %#ok<AGROW>
    crop_harvest_daily(:,i) = run_crop_harvest; %#ok<AGROW>
    n2o(i) = sum(n2o_daily(:,i));
    no3_leach(i) = sum(no3_leach_daily(:,i));
    no3_final(i) = no3_pool_daily(end,i);
    ph_final(i) = ph_daily(end,i);
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
    plot_file = fullfile(opt.ResultDirectory, ...
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
    no3_pool_daily,ph_daily,nee_daily, ...
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
    time_plot_file = fullfile(opt.ResultDirectory, ...
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

function [nee,crop_harvest] = carbon_dynamics( ...
        S,n_time,result_file)
required = {'NPP_H','NPP_L','Ccrown_t','R_litter','R_microbe', ...
    'R_ew','RB_L','topo','ColCfg'};
if ~all(isfield(S,required)) || ~isstruct(S.topo) || ...
        ~isfield(S.topo,'cell_area_m2') || ~isstruct(S.ColCfg)
    error('Plot_Subfield_12:IncompleteCarbonResult', ...
        'Carbon state is incomplete in %s.',result_file);
end

A_cell = S.topo.cell_area_m2(:);
n_col = numel(A_cell);
if n_col < 1 || numel(S.ColCfg) ~= n_col || ...
        any(~isfinite(A_cell)) || any(A_cell <= 0)
    error('Plot_Subfield_12:IncompleteCarbonResult', ...
        'Carbon cell areas or configuration are invalid in %s.',result_file);
end

crop_cell = false(n_col,1);
for c = 1:n_col
    if ~isfield(S.ColCfg(c),'PFT_OVERRIDE') || ...
            ~(ischar(S.ColCfg(c).PFT_OVERRIDE) || ...
            (isstring(S.ColCfg(c).PFT_OVERRIDE) && ...
            isscalar(S.ColCfg(c).PFT_OVERRIDE)))
        error('Plot_Subfield_12:IncompleteCarbonResult', ...
            'Cell PFT configuration is incomplete in %s.',result_file);
    end
    crop_cell(c) = isempty(char(S.ColCfg(c).PFT_OVERRIDE));
end

if ~(isnumeric(S.NPP_H) && isnumeric(S.NPP_L) && ...
        isnumeric(S.Ccrown_t) && isnumeric(S.R_litter) && ...
        isnumeric(S.R_microbe) && isnumeric(S.R_ew))
    error('Plot_Subfield_12:CarbonTimeLengthMismatch', ...
        'NEE arrays must be numeric in %s.',result_file);
end
if ~isnumeric(S.RB_L)
    error('Plot_Subfield_12:CropTimeLengthMismatch', ...
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
    error('Plot_Subfield_12:CarbonTimeLengthMismatch', ...
        'NEE arrays do not match the daily domain dimensions in %s.',result_file);
end
if size(S.RB_L,1) ~= n_time || size(S.RB_L,2) ~= cc || ...
        size(S.RB_L,3) ~= 7 || size(S.RB_L,4) ~= n_col
    error('Plot_Subfield_12:CropTimeLengthMismatch', ...
        'Crop-removal arrays do not match the daily domain dimensions in %s.',result_file);
end
if any(~isfinite(S.NPP_H(:))) || any(~isfinite(S.NPP_L(:))) || ...
        any(~isfinite(S.Ccrown_t(:))) || any(~isfinite(S.R_litter(:))) || ...
        any(~isfinite(S.R_microbe(:))) || any(~isfinite(S.R_ew(:)))
    error('Plot_Subfield_12:NonfiniteCarbonMeasurement', ...
        'NEE arrays contain nonfinite values in %s.',result_file);
end
if any(~isfinite(S.RB_L(:)))
    error('Plot_Subfield_12:NonfiniteCropMeasurement', ...
        'Crop-removal arrays contain nonfinite values in %s.',result_file);
end
if any(S.Ccrown_t(:) < 0) || any(S.Ccrown_t(:) > 1)
    error('Plot_Subfield_12:UnphysicalCarbonMeasurement', ...
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

function format_time_axis(ax)
ax.Color = 'w';
ax.XColor = 'k';
ax.YColor = 'k';
ax.XTickLabelRotation = 20;
grid(ax,'off');
end
