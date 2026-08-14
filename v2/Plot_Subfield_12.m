function [summary,fig,plot_file,time_fig,time_plot_file] = Plot_Subfield_12(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot summary values and daily dynamics for all twelve scenarios.
%%
%% Options
%%   ResultDirectory  result folder              [SIM/Subfield_12]
%%   NDays            shortened-run filename tag []
%%   Visible          figure visibility          ['on']
%%   OutputFile       output PNG                  [summary in result folder]
%%   TimeOutputFile   time-dynamics PNG           [time dynamics in result folder]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

here = fileparts(mfilename('fullpath'));
opt = struct('ResultDirectory',fullfile(here,'SIM','Subfield_12'), ...
    'NDays',[],'Visible','on','OutputFile','','TimeOutputFile','');
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
if ~isfolder(opt.ResultDirectory)
    error('Plot_Subfield_12:MissingDirectory', ...
        'Result directory not found: %s',opt.ResultDirectory);
end

RM = Run_Matrix_12();
n = numel(RM);
n2o = nan(n,1);
no3_leach = nan(n,1);
no3_final = nan(n,1);
ph_final = nan(n,1);
Date_d = datetime.empty(0,1);
n2o_daily = [];
no3_leach_daily = [];
no3_pool_daily = [];
ph_daily = [];
for k = 1:n
    tag = RM(k).name;
    if ~isempty(opt.NDays); tag = sprintf('%s__%dd',tag,opt.NDays); end
    result_file = fullfile(opt.ResultDirectory, ...
        sprintf('TC_2D_run%02d_%s.mat',k,tag));
    if ~isfile(result_file)
        error('Plot_Subfield_12:MissingResult','Result not found: %s',result_file);
    end
    S = load(result_file,'R','Domain','Datam');
    if ~isfield(S,'R') || ~isfield(S,'Domain') || ...
            S.R.id ~= k || ~strcmp(S.R.name,RM(k).name)
        error('Plot_Subfield_12:WrongResult', ...
            'Result does not match run %d: %s',k,result_file);
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
            'Run %d does not contain finite values for every plotted measurement.',k);
    end
    if ~isnumeric(S.Datam) || isempty(S.Datam) || size(S.Datam,2) < 3 || ...
            any(~isfinite(S.Datam(1,1:3)))
        error('Plot_Subfield_12:IncompleteResult', ...
            'Daily start date is missing or invalid in %s.',result_file);
    end
    run_date = datetime(S.Datam(1,1),S.Datam(1,2),S.Datam(1,3)) + ...
        caldays((0:n_time(1)-1).');
    if k == 1
        Date_d = run_date;
        n2o_daily = zeros(n_time(1),n);
        no3_leach_daily = zeros(n_time(1),n);
        no3_pool_daily = zeros(n_time(1),n);
        ph_daily = zeros(n_time(1),n);
    elseif ~isequal(run_date,Date_d)
        error('Plot_Subfield_12:TimeAxisMismatch', ...
            'Run %d does not use the same daily time axis as run 1.',k);
    end
    n2o_daily(:,k) = S.Domain.N2Oflx(:); %#ok<AGROW>
    no3_leach_daily(:,k) = S.Domain.LEAK_NO3(:); %#ok<AGROW>
    no3_pool_daily(:,k) = S.Domain.NO3_pool(:); %#ok<AGROW>
    ph_daily(:,k) = S.Domain.PH_d(:); %#ok<AGROW>
    n2o(k) = sum(n2o_daily(:,k));
    no3_leach(k) = sum(no3_leach_daily(:,k));
    no3_final(k) = no3_pool_daily(end,k);
    ph_final(k) = ph_daily(end,k);
end

summary = table([RM.id].',string({RM.layout}).',string({RM.treatment}).', ...
    n2o,no3_leach,no3_final,ph_final, ...
    'VariableNames',{'RunID','ScenarioGroup','Treatment', ...
    'CumulativeN2O_gN_m2','CumulativeNO3Leaching_gN_m2', ...
    'FinalNO3Pool_gN_m2','FinalSoilPH'});

group_labels = {'Whole field','Margin with crop','Margin with prairie','Vertical 1-D'};
treatment_labels = {'Control','Basalt','Calcite'};
treatment_colors = [0 0 0; 0 0.25 0.60; 0.35 0.70 0.90];
fig = figure('Color','w','Visible',char(opt.Visible), ...
    'Name','T&C-SMEW-N twelve-scenario summary');
fig.Position(3:4) = [1400 800];
tl = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
b = draw_panel(nexttile(tl),n2o,'Cumulative N_2O','g N m^{-2}', ...
    group_labels,treatment_colors);
draw_panel(nexttile(tl),no3_leach,'Cumulative NO_3 leaching','g N m^{-2}', ...
    group_labels,treatment_colors);
draw_panel(nexttile(tl),no3_final,'Final soil NO_3 pool','g N m^{-2}', ...
    group_labels,treatment_colors);
draw_panel(nexttile(tl),ph_final,'Final soil pH','pH', ...
    group_labels,treatment_colors);
lgd = legend(b,treatment_labels,'Orientation','horizontal','Box','off', ...
    'TextColor','k','Color','w');
lgd.Layout.Tile = 'north';
title(tl,'Treatment response across four scenario groups','Color','k');

suffix = '';
if ~isempty(opt.NDays); suffix = sprintf('__%dd',opt.NDays); end
if isempty(opt.OutputFile)
    plot_file = fullfile(opt.ResultDirectory, ...
        ['Subfield_12_summary' suffix '.png']);
else
    plot_file = char(opt.OutputFile);
end
exportgraphics(fig,plot_file,'Resolution',200);
fprintf('Saved: %s\n',plot_file);

time_fig = figure('Color','w','Visible',char(opt.Visible), ...
    'Name','T&C-SMEW-N twelve-scenario time dynamics');
time_fig.Position(3:4) = [1700 1150];
time_tl = tiledlayout(time_fig,4,4,'TileSpacing','compact','Padding','compact');
time_value = {cumsum(n2o_daily,1),cumsum(no3_leach_daily,1), ...
    no3_pool_daily,ph_daily};
time_label = {'Cumulative N_2O','Cumulative NO_3 leaching', ...
    'Soil NO_3 pool','Soil pH'};
time_units = {'g N m^{-2}','g N m^{-2}','g N m^{-2}','pH'};
time_tag = {'CumulativeN2O','CumulativeNO3Leaching','NO3Pool','SoilPH'};
time_axes = gobjects(16,1);
legend_lines = gobjects(3,1);
iax = 0;
for row = 1:4
    for group = 1:4
        iax = iax+1;
        ax = nexttile(time_tl);
        time_axes(iax) = ax;
        ids = (group-1)*3+(1:3);
        ln = plot(ax,Date_d,time_value{row}(:,ids),'LineWidth',1.15);
        for treatment = 1:3
            ln(treatment).Color = treatment_colors(treatment,:);
        end
        set(ln,'Tag',time_tag{row});
        format_time_axis(ax);
        if row == 1; title(ax,group_labels{group},'Color','k'); end
        if group == 1; ylabel(ax,time_units{row},'Color','k'); end
        if row == 4; xlabel(ax,'Date','Color','k'); end
        if group == 1
            text(ax,0.02,0.96,time_label{row},'Units','normalized', ...
                'VerticalAlignment','top','Color','k','FontWeight','bold');
        end
        if row == 1 && group == 1; legend_lines = ln; end
    end
end
linkaxes(time_axes,'x');
for row = 1:4
    ids = (row-1)*4+(1:4);
    row_lim = vertcat(time_axes(ids).YLim);
    set(time_axes(ids),'YLim',[min(row_lim(:,1)) max(row_lim(:,2))]);
end
time_lgd = legend(legend_lines,treatment_labels,'Orientation','horizontal', ...
    'Box','off','TextColor','k','Color','w');
time_lgd.Layout.Tile = 'north';
title(time_tl,'Daily treatment dynamics across four scenario groups','Color','k');

if isempty(opt.TimeOutputFile)
    time_plot_file = fullfile(opt.ResultDirectory, ...
        ['Subfield_12_time_dynamics' suffix '.png']);
else
    time_plot_file = char(opt.TimeOutputFile);
end
exportgraphics(time_fig,time_plot_file,'Resolution',200);
fprintf('Saved: %s\n',time_plot_file);

end

function b = draw_panel(ax,value,panel_title,y_label,group_labels,colors)
b = bar(ax,reshape(value,3,4).','grouped');
for treatment = 1:3; b(treatment).FaceColor = colors(treatment,:); end
ax.Color = 'w';
ax.XColor = 'k';
ax.YColor = 'k';
ax.GridColor = [0.75 0.75 0.75];
ax.XTick = 1:4;
ax.XTickLabel = group_labels;
ax.XTickLabelRotation = 18;
grid(ax,'on');
title(ax,panel_title,'Color','k');
ylabel(ax,y_label,'Color','k');
end

function format_time_axis(ax)
ax.Color = 'w';
ax.XColor = 'k';
ax.YColor = 'k';
ax.GridColor = [0.75 0.75 0.75];
ax.XTickLabelRotation = 20;
grid(ax,'on');
end
