function [fig,plot_file] = Plot_Uncertainty(summary_file,figure_dir,visible)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Per-year member spread of the Kantola parameter-uncertainty ensemble.
%%
%% The legacy T&C-SMEW uncertainty figure, drawn from the summary file
%% Summarize_Uncertainty writes (see the campaign README for the
%% contract) with MATLAB's own boxchart: one box per treatment and year
%% over the envelope members, the base parameter set as a filled marker,
%% the sensitivity-only members (KNmaxd at 0.02 d^-1, alpha at the next two
%% sweep points) as open markers of their own shape, and the field
%% observations in black with their standard errors.
%%
%%   [fig,plot_file] = Plot_Uncertainty(summary_file,figure_dir,visible)
%%
%% Three figures, PNG and FIG each, written to figure_dir:
%%   1) growing-season bulk pH per treatment, and the basalt-minus-control
%%      signal; Beerling et al. (2024) pH by depth interval (filled 0-10 cm,
%%      open 10-30 cm), the observed signal with its propagated SEM;
%%   2) calendar-year cumulative N2O and NO3 leaching per treatment;
%%   3) cumulative feedstock Ca and Mg release and the potential CDR of the
%%      amended treatments; Beerling et al. (2024) Dataset S06 in black.
%% Colours follow Plot_1D_Then_All_12: control black, basalt dark blue,
%% calcite light blue.  visible is 'on' (default) or 'off'.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 3 || isempty(visible); visible = 'on'; end
if ~(ischar(visible) || (isstring(visible) && isscalar(visible))) || ...
        ~ismember(char(visible),{'on','off'})
    error('Plot_Uncertainty:BadOptions','visible must be on or off.');
end
if ~isfile(summary_file)
    error('Plot_Uncertainty:MissingSummary', ...
        'Uncertainty summary not found: %s',summary_file);
end
S = load(summary_file,'run_ids','treatments','yearly','obs','in_envelope','members');
%% A summary may name what it covers (Summarize_Uncertainty: the plot, a
%% field grid, the scenario); the figures' titles say it.
label = 'Kantola 1-D';
contents = whos('-file',summary_file);
if any(strcmp({contents.name},'label'))
    L = load(summary_file,'label');
    label = char(L.label);
end
required = {'run_ids','treatments','yearly','obs','in_envelope','members'};
yearly_fields = {'years','pH_season_mean','N2O_cumulative_gN_m2', ...
    'NO3_leaching_cumulative_gN_m2','Ca_release_mol_m2','Mg_release_mol_m2', ...
    'CDR_potential_tCO2_ha'};
if ~all(isfield(S,required)) || ~all(isfield(S.yearly,yearly_fields)) || ...
        ~islogical(S.in_envelope) || ...
        numel(S.in_envelope) ~= size(S.yearly.pH_season_mean,1) || ~S.in_envelope(1) || ...
        numel(S.members) ~= numel(S.in_envelope) || ~isfield(S.members,'parameter')
    error('Plot_Uncertainty:BadSummary', ...
        'Uncertainty summary is incomplete: %s',summary_file);
end
if ~isfolder(figure_dir); mkdir(figure_dir); end

treatments = cellstr(S.treatments(:).');
years = S.yearly.years(:).';
in_envelope = S.in_envelope(:).';
sensitivity = find(~in_envelope);
sensitivity_names = arrayfun(@sensitivity_label,S.members(sensitivity),'UniformOutput',false);
palette = struct('control',[0 0 0],'basalt',[0 0.25 0.60],'calcite',[0.35 0.70 0.90]);
colors = zeros(numel(treatments),3);
for r = 1:numel(treatments)
    colors(r,:) = palette.(treatments{r});
end
labels = treatments;
for r = 1:numel(labels); labels{r}(1) = upper(labels{r}(1)); end
%% One slot per treatment around each year tick; amended-only panels
%% keep the same slots so a treatment sits at the same offset everywhere.
offset = 0.26*((1:numel(treatments))-(numel(treatments)+1)/2);
%% Observations sit just left of their treatment's box centre, so the base
%% marker on the box stays visible beside them.
obs_shift = -0.07;
obs = S.obs;
obs_years = obs.years(:).';
member_note = sprintf(['Boxes: %d envelope members, whiskers 1.5 IQR; ' ...
    'filled marker: base parameter set; open markers: sensitivity-only members'], ...
    nnz(in_envelope));
fig = gobjects(3,1);
plot_file = cell(3,1);

%% ---- 1) Growing-season bulk pH and the treatment signal -----------------
fig(1) = figure('Color','w','Visible',char(visible), ...
    'Name','Kantola uncertainty: growing-season bulk pH');
fig(1).Position(3:4) = [1100 520];
tl = tiledlayout(fig(1),1,2,'TileSpacing','compact','Padding','compact');
ax = nexttile(tl);
[handles,sensitivity_handles] = treatment_boxes(ax,years,offset, ...
    S.yearly.pH_season_mean,colors,1:numel(treatments),in_envelope,sensitivity);
depths = {'0_10','10_30'};
obs_handles = gobjects(1,2);
obs_names = {'Observed 0-10 cm','Observed 10-30 cm'};
for r = 1:numel(treatments)
    for d = 1:2
        name = sprintf('%s_%s',treatments{r},depths{d});
        if ~isfield(obs.pH,name); continue; end
        eb = observed_points(ax,obs_years+offset(r)+obs_shift,obs.pH.(name).value, ...
            obs.pH.(name).sem,depths{d},sprintf('Observed_%s_pH_%s',treatments{r},depths{d}));
        if ~isgraphics(obs_handles(d)); obs_handles(d) = eb; end
    end
end
drawn = isgraphics(obs_handles);
obs_handles = obs_handles(drawn); obs_names = obs_names(drawn);
format_year_axis(ax,years);
ylabel(ax,'Bulk pH (1 May - 30 Sep mean)','Color','k');
title(ax,'Growing-season bulk pH','Color','k');

ax = nexttile(tl);
bas = find(strcmp(treatments,'basalt'),1);
if isfield(S.yearly,'pH_season_signal') && ~isempty(bas)
    member_boxes(ax,years,S.yearly.pH_season_signal(in_envelope,:),colors(bas,:));
    base_marker(ax,years,S.yearly.pH_season_signal(1,:),colors(bas,:));
    sensitivity_markers(ax,years,S.yearly.pH_season_signal(sensitivity,:),colors(bas,:));
    %% Observed signal: basalt minus control per depth interval, SEM
    %% propagated as the root sum of squares of the two block SEMs.
    for depth = {'0_10','10_30'}
        c = sprintf('control_%s',depth{1}); b = sprintf('basalt_%s',depth{1});
        if ~isfield(obs.pH,c) || ~isfield(obs.pH,b); continue; end
        observed_points(ax,obs_years+obs_shift,obs.pH.(b).value-obs.pH.(c).value, ...
            hypot(obs.pH.(b).sem,obs.pH.(c).sem),depth{1}, ...
            sprintf('Observed_signal_pH_%s',depth{1}));
    end
    yline(ax,0,'-','Color',[0.6 0.6 0.6],'HandleVisibility','off');
else
    text(ax,0.5,0.5,'No basalt-minus-control signal in the summary', ...
        'Units','normalized','HorizontalAlignment','center','Color','k');
end
format_year_axis(ax,years);
ylabel(ax,'\DeltapH (basalt - control)','Color','k');
title(ax,'Basalt pH signal','Color','k');
lgd = legend([handles obs_handles sensitivity_handles], ...
    [labels obs_names sensitivity_names],'NumColumns',4,'Box','off','TextColor','k');
lgd.Layout.Tile = 'south';
title(tl,[label ': parameter-uncertainty ensemble, bulk pH'],'Color','k');
subtitle(tl,member_note,'Color',[0.3 0.3 0.3],'FontSize',9);
plot_file{1} = save_pair(fig(1),fullfile(figure_dir,'Kantola_uncertainty_ph.png'));

%% ---- 2) Nitrogen: calendar-year cumulative N2O and NO3 leaching --------
fig(2) = figure('Color','w','Visible',char(visible), ...
    'Name','Kantola uncertainty: nitrogen losses');
fig(2).Position(3:4) = [1100 520];
tl = tiledlayout(fig(2),1,2,'TileSpacing','compact','Padding','compact');
metric = {'N2O_cumulative_gN_m2','NO3_leaching_cumulative_gN_m2'};
metric_title = {'Cumulative N_2O per calendar year','Cumulative NO_3 leaching per calendar year'};
metric_units = {'g N m^{-2} yr^{-1}','g N m^{-2} yr^{-1}'};
for m = 1:2
    ax = nexttile(tl);
    [handles,sensitivity_handles] = treatment_boxes(ax,years,offset, ...
        S.yearly.(metric{m}),colors,1:numel(treatments),in_envelope,sensitivity);
    format_year_axis(ax,years);
    ylabel(ax,metric_units{m},'Color','k');
    title(ax,metric_title{m},'Color','k');
end
lgd = legend([handles sensitivity_handles],[labels sensitivity_names], ...
    'NumColumns',4,'Box','off','TextColor','k');
lgd.Layout.Tile = 'south';
title(tl,[label ': parameter-uncertainty ensemble, nitrogen losses'],'Color','k');
subtitle(tl,member_note,'Color',[0.3 0.3 0.3],'FontSize',9);
plot_file{2} = save_pair(fig(2),fullfile(figure_dir,'Kantola_uncertainty_nitrogen.png'));

%% ---- 3) Feedstock release and potential CDR of the amended treatments --
fig(3) = figure('Color','w','Visible',char(visible), ...
    'Name','Kantola uncertainty: feedstock release and CDR potential');
fig(3).Position(3:4) = [1100 820];
tl = tiledlayout(fig(3),2,2,'TileSpacing','compact','Padding','compact');
amended = find(~strcmp(treatments,'control'));
metric = {'Ca_release_mol_m2','Mg_release_mol_m2','CDR_potential_tCO2_ha'};
metric_title = {'Cumulative Ca release','Cumulative Mg release', ...
    'Potential CDR, 2 x (Ca + Mg) as CO_2'};
metric_units = {'mol Ca m^{-2}','mol Mg m^{-2}','t CO_2 ha^{-1}'};
%% Dataset S06 release is in mmol m-2; the model arrays are mol m-2.
obs_series = {scaled_obs(obs,'Ca_release_mmol_m2',1e-3), ...
    scaled_obs(obs,'Mg_release_mmol_m2',1e-3), ...
    scaled_obs(obs,'CDR_potential_tCO2_ha',1)};
tile_span = {[1 1],[1 1],[1 2]};
obs_handle = gobjects(1);
for m = 1:3
    ax = nexttile(tl,tile_span{m});
    [handles,sensitivity_handles] = treatment_boxes(ax,years,offset, ...
        S.yearly.(metric{m}),colors,amended,in_envelope,sensitivity);
    if ~isempty(obs_series{m}) && ~isempty(bas)
        eb = observed_points(ax,obs_years+offset(bas)+obs_shift,obs_series{m}.value, ...
            obs_series{m}.sem,'0_10',sprintf('Observed_%s',metric{m}));
        if m == 3; obs_handle = eb; end
    end
    format_year_axis(ax,years);
    ylabel(ax,metric_units{m},'Color','k');
    title(ax,metric_title{m},'Color','k');
end
legend_handles = [handles sensitivity_handles];
legend_names = [labels(amended) sensitivity_names];
if isgraphics(obs_handle)
    legend_handles = [handles obs_handle sensitivity_handles];
    legend_names = [labels(amended) {'Observed basalt (Beerling et al. 2024)'} sensitivity_names];
end
lgd = legend(legend_handles,legend_names,'NumColumns',4,'Box','off','TextColor','k');
lgd.Layout.Tile = 'south';
title(tl,[label ': parameter-uncertainty ensemble, weathering'],'Color','k');
subtitle(tl,member_note,'Color',[0.3 0.3 0.3],'FontSize',9);
plot_file{3} = save_pair(fig(3),fullfile(figure_dir,'Kantola_uncertainty_weathering.png'));
end

%% =======================================================================
function [handles,sensitivity_handles] = treatment_boxes(ax,years,offset, ...
        value,colors,runs,in_envelope,sensitivity)
%% Boxes over the envelope members for the selected runs of one yearly
%% array [n_members x n_years x n_runs], each in its treatment slot, then
%% the base member and the sensitivity-only members on top. The markers of
%% the first run are returned for the legend (one entry per member).
handles = gobjects(1,numel(runs));
for k = 1:numel(runs)
    handles(k) = member_boxes(ax,years+offset(runs(k)), ...
        value(in_envelope,:,runs(k)),colors(runs(k),:));
end
sensitivity_handles = gobjects(1,0);
for k = 1:numel(runs)
    base_marker(ax,years+offset(runs(k)),value(1,:,runs(k)),colors(runs(k),:));
    h = sensitivity_markers(ax,years+offset(runs(k)),value(sensitivity,:,runs(k)), ...
        colors(runs(k),:));
    if k == 1; sensitivity_handles = h; end
end
end

%% =======================================================================
function h = sensitivity_markers(ax,x,value,color)
%% The sensitivity-only members as open markers in the treatment colour,
%% one shape per member: value is [n_sensitivity x n_years].
shapes = 'dv^sph';
h = gobjects(1,size(value,1));
for j = 1:size(value,1)
    h(j) = plot(ax,x,value(j,:),shapes(mod(j-1,numel(shapes))+1), ...
        'MarkerSize',6,'MarkerFaceColor','w','MarkerEdgeColor',color, ...
        'LineWidth',0.9,'LineStyle','none','Tag','SensitivityMember', ...
        'HandleVisibility','off');
end
end

%% =======================================================================
function label = sensitivity_label(member)
%% Legend text of a sensitivity-only member: the parameter and its value.
switch member.parameter
    case 'KNmaxd'
        label = sprintf('KNmaxd = %.3g d^{-1} (%s)',member.KNmaxd,member.level);
    case 'AppeloAlpha'
        label = sprintf('\\alpha = %.4g',member.AppeloAlpha);
    otherwise
        label = sprintf('%s (%s)',member.parameter,member.level);
end
end

%% =======================================================================
function bx = member_boxes(ax,x,value,color)
%% One box per year over the members: value is [n_members x n_years].
hold(ax,'on');
xg = repmat(x(:).',size(value,1),1);
bx = boxchart(ax,xg(:),value(:),'BoxFaceColor',color,'BoxEdgeColor',color, ...
    'WhiskerLineColor',color,'MarkerStyle','none','LineWidth',0.75, ...
    'BoxWidth',0.2,'BoxFaceAlpha',0.35,'Tag','MemberBox');
end

%% =======================================================================
function base_marker(ax,x,value,color)
plot(ax,x,value,'o','MarkerSize',5,'MarkerFaceColor',color, ...
    'MarkerEdgeColor','w','LineWidth',0.5,'LineStyle','none', ...
    'Tag','BaseMember','HandleVisibility','off');
end

%% =======================================================================
function eb = observed_points(ax,x,value,sem,depth,tag)
%% Observations in black: filled 0-10 cm, open 10-30 cm, SEM as the bar.
face = 'k';
if strcmp(depth,'10_30'); face = 'w'; end
keep = isfinite(value);
eb = errorbar(ax,x(keep),value(keep),sem(keep),'o','Color','k', ...
    'MarkerSize',5,'MarkerFaceColor',face,'LineStyle','none', ...
    'LineWidth',0.9,'CapSize',4,'Tag',tag,'HandleVisibility','off');
end

%% =======================================================================
function series = scaled_obs(obs,name,scale)
series = [];
if ~isfield(obs,name); return; end
series = struct('value',scale*obs.(name).value,'sem',scale*obs.(name).sem);
end

%% =======================================================================
function format_year_axis(ax,years)
ax.Color = 'w'; ax.XColor = 'k'; ax.YColor = 'k';
step = max(1,ceil(numel(years)/10));
ax.XTick = years(1:step:end);
ax.XTickLabel = arrayfun(@(y) sprintf('%d',y),ax.XTick,'UniformOutput',false);
ax.XLim = [years(1)-0.5 years(end)+0.5];
if numel(years) > 1 && numel(years) <= 10
    xline(ax,years(1:end-1)+0.5,':','Color',[0.8 0.8 0.8],'HandleVisibility','off');
end
xlabel(ax,'Year','Color','k');
grid(ax,'off');
hold(ax,'off');
end

%% =======================================================================
function png_file = save_pair(fig,png_file)
exportgraphics(fig,png_file,'Resolution',200);
[folder,name] = fileparts(png_file);
savefig(fig,fullfile(folder,[name '.fig']));
fprintf('Saved: %s\n',png_file);
fprintf('Saved: %s\n',fullfile(folder,[name '.fig']));
end
