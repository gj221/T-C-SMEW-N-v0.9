function out = Plot_Field_Maps(result_file,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot_Field_Maps
%%
%% Maps of a 2-D run, cell by cell, over a sequence of time windows: how
%% N2O, NO3 leaching, NEE, pH or any other Field_Cell_Variables quantity
%% varies across the plot and how that pattern changes in time. One
%% figure per variable: the field-mean daily series with the cell range
%% and the window boundaries on top, a map per window below, one colour
%% scale for all windows. Every value drawn is also written to a CSV.
%%
%%   Plot_Field_Maps(result_file)              N2O, NO3 leaching, NEE, pH by year
%%   Plot_Field_Maps(f,'Variables',{'N2O','pH_bulk'},'Interval','season')
%%   Plot_Field_Maps(f,'Variables','NEE','Interval',30,'Period',{'2017-01-01','2018-12-31'})
%%   Plot_Field_Maps(f,'Variables','pH_bulk','Snapshots',datetime(2017:2020,8,1),'Statistic','end')
%%   Plot_Field_Maps(f,'Reference',control_file)        run minus a run on the same grid
%%   Plot_Field_Maps(f,'Interval','month','Static',false,'Animation',true)
%%   Field_Cell_Variables('list')                      every variable key
%%
%% Options
%%   'Variables'    a key, a cell of keys, or Field_Cell_Variables rows
%%                  [N2O NO3_leaching NEE pH_bulk]
%%   'Interval'     'all', 'year', 'growing_season', 'season', 'month',
%%                  'week', 'day', N days, or a calendarDuration   ['year']
%%   'Windows'      explicit [start stop] rows; replaces Interval
%%   'Snapshots'    one-day windows at these dates; replaces Interval
%%   'Period'       [t0 t1]: clip the windows
%%   'YearStartMonth', 'GrowingSeasonMonths'  see Field_Cell_Windows
%%   'Statistic'    overrides the window statistic of every variable (one
%%                  name) or of each (a cell, '' keeps the default):
%%                  sum, mean, hmean, end, start, change, max, min
%%   'Reference'    a result on the same grid and dates; maps show this
%%                  run minus the reference
%%   'Anomaly'      subtract each window's field mean            [false]
%%   'ColorLimits'  'shared' (one scale per variable), 'window' (one per
%%                  map), [lo hi], or a cell of these per variable ['shared']
%%   'Layout'       'per_variable', or 'combined' (one figure, a row per
%%                  variable, a column per window)        ['per_variable']
%%   'TimeSeries'   draw the field-mean series above the maps     [true]
%%   'TraceCells'   cells whose own series are drawn and marked     [[]]
%%   'Coordinates'  'auto', 'local' or 'projected'              ['auto']
%%   'ClipToField'  draw the in-field part of perimeter squares   [true]
%%   'Static'       write the window figures                       [true]
%%   'Animation'    also write a GIF stepping through the windows [false]
%%   'FrameDelay'   seconds per GIF frame                           [0.8]
%%   'MaxPanels'    most windows on one static figure                [36]
%%   'FigureDirectory'  [v2/figures/field_maps/<result name>]
%%   'Export'       write PNG, FIG, CSV, MAT (and GIF)             [true]
%%   'Visible'      'on' or 'off'                                  ['on']
%%   'Resolution'   PNG dpi                                         [250]
%%
%% Window statistics: fluxes are totals over the window (kg N ha-1 for
%% N2O over a year), pools are read at the window end, pH is the mean
%% hydrogen-ion activity. Field means weigh cells by their in-field area
%% (Field_Aggregation_Areas), and pH by hydrogen-ion activity
%% (pH_Area_Weighted), as Domain does.
%%
%% out holds the windows, and per variable the values [windows x cells],
%% the field means, the colour limits and the files; out.table is every
%% value in long form (cell, zone, vegetation, feedstock, position).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
p.addParameter('Variables',{'N2O','NO3_leaching','NEE','pH_bulk'});
p.addParameter('Interval','year');
p.addParameter('Windows',[]);
p.addParameter('Snapshots',[]);
p.addParameter('Period',[]);
p.addParameter('YearStartMonth',1);
p.addParameter('GrowingSeasonMonths',[5 10]);
p.addParameter('Statistic','');
p.addParameter('Reference','');
p.addParameter('Anomaly',false);
p.addParameter('ColorLimits','shared');
p.addParameter('Layout','per_variable',@(x) any(strcmpi(x,{'per_variable','combined'})));
p.addParameter('TimeSeries',true);
p.addParameter('TraceCells',[]);
p.addParameter('Coordinates','auto');
p.addParameter('ClipToField',true);
p.addParameter('Static',true);
p.addParameter('Animation',false);
p.addParameter('FrameDelay',0.8);
p.addParameter('MaxPanels',36);
p.addParameter('FigureDirectory','');
p.addParameter('Export',true);
p.addParameter('Visible','on',@(x) any(strcmpi(x,{'on','off'})));
p.addParameter('Resolution',250);
p.parse(varargin{:});
o = p.Results;

here = fileparts(mfilename('fullpath'));
v2 = fileparts(here);
res = Field_Cell_Result(result_file);
if ~res.has_series
    error('Plot_Field_Maps:NotAResult','%s holds a grid, not a model result.',res.name);
end
N = res.N_cells;
specs = variable_rows(o.Variables,o.Statistic);
n_var = numel(specs);
[win,tag,words] = Field_Cell_Windows(res.day,o.Interval,'Windows',o.Windows, ...
    'Snapshots',o.Snapshots,'Period',o.Period,'YearStartMonth',o.YearStartMonth, ...
    'GrowingSeasonMonths',o.GrowingSeasonMonths);
n_win = numel(win);
if o.Static && n_win > o.MaxPanels
    error('Plot_Field_Maps:TooManyWindows', ...
        ['%d windows exceed MaxPanels (%d) for one figure: choose a longer ' ...
        'Interval or a Period, raise MaxPanels, or set Static false and ' ...
        'Animation true.'],n_win,o.MaxPanels);
end
trace = o.TraceCells(:).';
if any(trace < 1 | trace > N | trace ~= round(trace))
    error('Plot_Field_Maps:BadTraceCells','TraceCells are cell numbers from 1 to %d.',N);
end
ref = [];
if ~isempty(o.Reference)
    ref = Field_Cell_Result(o.Reference);
    check_reference(res,ref);
end
limit_mode = per_variable(o.ColorLimits,n_var);

G = Field_Cell_Geometry(res.topo,'Coordinates',o.Coordinates,'ClipToField',o.ClipToField);
in_field = res.weight(:).' > 0;

%% Values: per window and cell, the field means, and the daily series.
V = struct('key',{},'label',{},'units',{},'stat',{},'values',{},'field_mean',{}, ...
    'run_values',{},'reference_values',{},'daily_mean',{},'daily_low',{}, ...
    'daily_high',{},'daily_reference',{},'daily_trace',{},'clim',{}, ...
    'diverging',{},'png','','fig','','gif','');
for v = 1:n_var
    s = specs(v);
    hmean = strcmp(s.stat,'hmean');
    daily = Field_Cell_Series(res,s);
    values = s.scale*window_values(daily,win,s.stat);
    fm = weighted_mean(values,res.weight,hmean);
    run_values = values;
    reference_values = [];
    daily_reference = [];
    if ~isempty(ref)
        daily_ref = Field_Cell_Series(ref,s);
        reference_values = s.scale*window_values(daily_ref,win,s.stat);
        values = values - reference_values;
        fm = fm - weighted_mean(reference_values,ref.weight,hmean);
        daily_reference = s.scale*weighted_mean(daily_ref,ref.weight,hmean);
    end
    if o.Anomaly
        values = values - weighted_mean(values,res.weight,hmean && isempty(ref));
    end
    entry.key = s.key;
    entry.label = s.label;
    entry.units = s.units;
    entry.stat = s.stat;
    entry.values = values;
    entry.field_mean = fm;
    entry.run_values = run_values;
    entry.reference_values = reference_values;
    entry.daily_mean = s.scale*weighted_mean(daily,res.weight,hmean);
    entry.daily_low = s.scale*min(daily(:,in_field),[],2);
    entry.daily_high = s.scale*max(daily(:,in_field),[],2);
    if s.scale < 0
        [entry.daily_low,entry.daily_high] = deal(entry.daily_high,entry.daily_low);
    end
    entry.daily_reference = daily_reference;
    entry.daily_trace = s.scale*daily(:,trace);
    diverging = s.diverging || ~isempty(ref) || o.Anomaly;
    entry.clim = colour_limits(values(:,in_field),limit_mode{v},diverging);
    entry.diverging = diverging;
    entry.png = ''; entry.fig = ''; entry.gif = '';
    V(v) = entry;
end

%% Where the files go.
figure_dir = o.FigureDirectory;
if isempty(figure_dir)
    figure_dir = fullfile(v2,'figures','field_maps',res.short_name);
end
suffix = tag;
if ~isempty(ref)
    suffix = sprintf('%s_minus_%s',suffix,ref.short_name(1:min(end,24)));
end
if o.Anomaly
    suffix = [suffix '_anomaly'];
end
if o.Export && ~isfolder(figure_dir); mkdir(figure_dir); end
title_note = '';
if ~isempty(ref)
    title_note = sprintf(', minus %s',tex_escape(ref.run_label));
end
if o.Anomaly
    title_note = [title_note ', departure from the field mean'];
end

figures = gobjects(0);
if o.Static
    if strcmpi(o.Layout,'combined')
        fig = combined_figure(res,G,V,win,words,title_note,o);
        figures(end+1) = fig;
        if o.Export
            [png,fig_file] = export_pair(fig,fullfile(figure_dir,['maps_' suffix]),o.Resolution);
            [V.png] = deal(png); [V.fig] = deal(fig_file);
        end
    else
        for v = 1:n_var
            fig = variable_figure(res,G,V(v),win,words,title_note,trace,o);
            figures(end+1) = fig; %#ok<AGROW>
            if o.Export
                [V(v).png,V(v).fig] = export_pair(fig, ...
                    fullfile(figure_dir,[V(v).key '_' suffix]),o.Resolution);
            end
        end
    end
end
if o.Animation && o.Export
    for v = 1:n_var
        V(v).gif = write_animation(res,G,V(v),win,title_note, ...
            fullfile(figure_dir,[V(v).key '_' suffix '.gif']),o);
    end
end

%% Every value drawn, in long form.
T = value_table(res,G,V,win,~isempty(ref));
out = struct('result',res.file,'reference','','windows',win,'interval',tag, ...
    'variables',V,'table',T,'figures',figures,'figure_directory',figure_dir, ...
    'csv','','mat','');
if ~isempty(ref); out.reference = ref.file; end
if o.Export
    out.csv = fullfile(figure_dir,['cell_values_' suffix '.csv']);
    writetable(T,out.csv);
    out.mat = fullfile(figure_dir,['cell_values_' suffix '.mat']);
    cells = cell_table(res,G);
    windows = rmfield(win,'rows');
    variables = rmfield(V,{'png','fig','gif'});
    save(out.mat,'cells','windows','variables','-v7');
    fprintf('Saved: %s\n',out.csv);
end
end

%% =======================================================================
function specs = variable_rows(variables,statistic)
if ischar(variables) || (isstring(variables) && isscalar(variables))
    variables = {char(variables)};
elseif isstring(variables)
    variables = cellstr(variables);
elseif isstruct(variables)
    variables = num2cell(variables);
end
specs = cellfun(@Field_Cell_Variables,variables(:).','UniformOutput',false);
specs = [specs{:}];
stats = {'sum','mean','hmean','end','start','change','max','min'};
if ~isempty(statistic)
    if ischar(statistic) || isstring(statistic)
        statistic = repmat(cellstr(statistic),1,numel(specs));
    end
    if numel(statistic) ~= numel(specs)
        error('Plot_Field_Maps:BadStatistic','Give one statistic, or one per variable.');
    end
    for k = 1:numel(specs)
        if isempty(statistic{k}); continue; end
        if ~ismember(statistic{k},stats)
            error('Plot_Field_Maps:BadStatistic','Statistic must be one of: %s.', ...
                strjoin(stats,', '));
        end
        specs(k).stat = char(statistic{k});
    end
end
end

function c = per_variable(x,n)
if iscell(x)
    if numel(x) ~= n
        error('Plot_Field_Maps:BadColorLimits','Give one colour-limit setting per variable.');
    end
    c = x(:).';
else
    c = repmat({x},1,n);
end
end

function check_reference(res,ref)
if ~ref.has_series
    error('Plot_Field_Maps:ReferenceGrid','The reference is not a model result.');
end
same = ref.N_cells == res.N_cells && ...
    max(abs(ref.topo.cell_x_lo(:)-res.topo.cell_x_lo(:))) < 1e-6 && ...
    max(abs(ref.topo.cell_y_lo(:)-res.topo.cell_y_lo(:))) < 1e-6;
if ~same
    error('Plot_Field_Maps:ReferenceGrid','The reference run uses a different grid.');
end
if ~isequal(ref.Date_d,res.Date_d)
    error('Plot_Field_Maps:ReferencePeriod','The reference run covers different dates.');
end
end

%% =======================================================================
function v = window_values(X,win,stat)
v = nan(numel(win),size(X,2));
for k = 1:numel(win)
    seg = X(win(k).rows,:);
    switch stat
        case 'sum'
            v(k,:) = sum(seg,1,'omitnan');
            v(k,all(isnan(seg),1)) = NaN;
        case 'mean'
            v(k,:) = mean(seg,1,'omitnan');
        case 'hmean'
            v(k,:) = -log10(mean(10.^(-seg),1,'omitnan'));
        case 'end'
            v(k,:) = X(win(k).last_row,:);
        case 'start'
            v(k,:) = X(win(k).prev_row,:);
        case 'change'
            v(k,:) = X(win(k).last_row,:) - X(win(k).prev_row,:);
        case 'max'
            v(k,:) = max(seg,[],1,'omitnan');
        case 'min'
            v(k,:) = min(seg,[],1,'omitnan');
        otherwise
            error('Plot_Field_Maps:BadStatistic','Unknown statistic "%s".',stat);
    end
end
end

function m = weighted_mean(X,weight,hmean)
%% Area-weighted mean over the cells of each row, cells without a value
%% left out; hydrogen-ion activity for pH, as pH_Area_Weighted.
ok = isfinite(X);
W = ok.*weight(:).';
total = sum(W,2);
W = W./total;
X(~ok) = 0;
if hmean
    m = -log10(sum(W.*10.^(-X),2));
else
    m = sum(W.*X,2);
end
m(total == 0) = NaN;
end

function lim = colour_limits(values,mode,diverging)
if isnumeric(mode) && numel(mode) == 2
    lim = sort(double(mode(:).'));
    return
end
if (ischar(mode) || isstring(mode)) && strcmpi(mode,'window')
    lim = zeros(size(values,1),2);
    for k = 1:size(values,1)
        lim(k,:) = one_range(values(k,:),diverging);
    end
    return
end
lim = one_range(values(:),diverging);
end

function lim = one_range(v,diverging)
v = v(isfinite(v));
if isempty(v)
    lim = [0 1];
elseif diverging
    m = max(abs(v));
    if m == 0; m = 1; end
    lim = [-m m];
else
    lim = [min(v) max(v)];
    if lim(2) <= lim(1)
        lim = lim + [-1 1]*max(abs(lim(1))*1e-3,1e-9);
    end
end
end

%% =======================================================================
function fig = variable_figure(res,G,entry,win,words,title_note,trace,o)
n_win = numel(win);
nc = columns_for(n_win);
nr = ceil(n_win/nc);
strip = logical(o.TimeSeries);
fig = figure('Color','w','Visible',o.Visible,'Name',entry.key);
fig.Position(3:4) = [min(2400,90+nc*270+110) min(2600,90+strip*230+nr*265)];
tl = tiledlayout(fig,double(strip)+2*nr,nc,'TileSpacing','compact','Padding','compact');
title(tl,{tex_escape(res.run_label), variable_title(entry,words,title_note)}, ...
    'Color','k','FontWeight','bold','FontSize',12);
if strip
    ax_s = nexttile(tl,1,[1 nc]);
    draw_strip(ax_s,res,entry,win,trace);
end
cmap = colour_map(entry.diverging);
per_window = size(entry.clim,1) > 1;
ax = gobjects(n_win,1);
for k = 1:n_win
    row = double(strip) + 2*floor((k-1)/nc) + 1;
    col = mod(k-1,nc)+1;
    ax(k) = nexttile(tl,(row-1)*nc+col,[2 1]);
    Draw_Field_Cells(ax(k),G,entry.values(k,:),'ScaleBar',k == 1,'NorthArrow',k == 1);
    colormap(ax(k),cmap);
    set(ax(k),'CLim',entry.clim(min(k,size(entry.clim,1)),:));
    mark_cells(ax(k),G,trace);
    title(ax(k),win(k).label,'Color','k','FontSize',10);
    subtitle(ax(k),sprintf('%s | mean %s',win(k).short_range, ...
        format_value(entry.field_mean(k),entry)),'Color',[0.25 0.25 0.25],'FontSize',8);
    if per_window
        cb = colorbar(ax(k),'southoutside');
        cb.Color = 'k';
    end
end
if ~per_window
    cb = colorbar(ax(end));
    cb.Layout.Tile = 'east';
    cb.Color = 'k';
    cb.Label.String = sprintf('%s [%s]',entry.label,entry.units);
end
end

function fig = combined_figure(res,G,V,win,words,title_note,o)
%% A row per variable (a nested layout carrying its own colour bar and
%% label), a column per window.
n_var = numel(V);
n_win = numel(win);
fig = figure('Color','w','Visible',o.Visible,'Name','Field maps');
fig.Position(3:4) = [min(2600,140+n_win*250+130) min(2600,110+n_var*270)];
outer = tiledlayout(fig,n_var,1,'TileSpacing','compact','Padding','compact');
title(outer,{tex_escape(res.run_label), sprintf('Cell maps %s%s',words,title_note)}, ...
    'Color','k','FontWeight','bold','FontSize',12);
for v = 1:n_var
    inner = tiledlayout(outer,1,n_win,'TileSpacing','compact','Padding','compact');
    inner.Layout.Tile = v;
    cmap = colour_map(V(v).diverging);
    for k = 1:n_win
        ax = nexttile(inner);
        Draw_Field_Cells(ax,G,V(v).values(k,:),'ScaleBar',v == 1 && k == 1, ...
            'NorthArrow',v == 1 && k == 1);
        colormap(ax,cmap);
        set(ax,'CLim',V(v).clim(min(k,size(V(v).clim,1)),:));
        if v == 1
            title(ax,win(k).label,'Color','k','FontSize',10);
        end
        text(ax,0.5,-0.01,sprintf('mean %s',format_value(V(v).field_mean(k),V(v))), ...
            'Units','normalized','HorizontalAlignment','center', ...
            'VerticalAlignment','top','FontSize',8,'Color',[0.25 0.25 0.25]);
    end
    cb = colorbar(ax);
    cb.Layout.Tile = 'east';
    cb.Color = 'k';
    cb.Label.String = sprintf('[%s]',V(v).units);
    ylabel(inner,sprintf('%s (%s)',V(v).label,stat_words(V(v).stat)), ...
        'FontSize',10,'FontWeight','bold','Color','k');
end
end

function draw_strip(ax,res,entry,win,trace)
%% Field-mean daily series, the range over the cells, any traced cells,
%% the reference, and the window boundaries. Each day's value sits at its
%% midday, inside the day the window edges bound.
hold(ax,'on');
rows = 2:res.NNd;
x = res.day(rows) + hours(12);
h = gobjects(0); names = {};
mark = {};
if numel(rows) <= 31
    mark = {'Marker','.','MarkerSize',12};      %% a short run: show each day
end
lo = entry.daily_low(rows); hi = entry.daily_high(rows);
if any(isfinite(lo)) && any(abs(hi-lo) > 0)
    band = area(ax,x,[lo hi-lo],'LineStyle','none');
    band(1).FaceColor = 'none';
    band(2).FaceColor = [0.70 0.80 0.92];
    band(2).FaceAlpha = 0.6;
    h(end+1) = band(2); names{end+1} = 'Range over the cells';
end
colours = lines(max(1,numel(trace)));
for k = 1:numel(trace)
    h(end+1) = plot(ax,x,entry.daily_trace(rows,k),'-','Color',colours(k,:), ...
        'LineWidth',0.8,mark{:}); %#ok<AGROW>
    if res.is_main(trace(k)); zone = 'main'; else; zone = 'margin'; end
    names{end+1} = sprintf('Cell %d (%s)',trace(k),zone); %#ok<AGROW>
end
if ~isempty(entry.daily_reference)
    h(end+1) = plot(ax,x,entry.daily_reference(rows),'--','Color',[0.45 0.45 0.45], ...
        'LineWidth',1,mark{:});
    names{end+1} = 'Reference field mean';
end
h(end+1) = plot(ax,x,entry.daily_mean(rows),'-','Color',[0.10 0.35 0.65],'LineWidth',1.3,mark{:});
names{end+1} = 'Field mean';
edges = unique([[win.start] [win.stop]+caldays(1)]);
xline(ax,edges,':','Color',[0.4 0.4 0.4],'LineWidth',0.6,'HandleVisibility','off');
xlim(ax,[res.day(2) res.day(end)+caldays(1)]);
shown = [lo; hi; entry.daily_mean(rows); reshape(entry.daily_trace(rows,:),[],1)];
if ~isempty(entry.daily_reference)
    shown = [shown; entry.daily_reference(rows)];
end
ylim(ax,strip_limits(shown));
ax.Color = 'w'; ax.XColor = 'k'; ax.YColor = 'k'; ax.Box = 'on';
grid(ax,'on');
units = entry.units;
if strcmp(entry.stat,'sum') || strcmp(entry.stat,'change')
    units = [units ' d^{-1}'];
end
ylabel(ax,units,'Color','k');
legend(ax,h,names,'Location','northoutside','Box','off','FontSize',8, ...
    'TextColor','k','Orientation','horizontal','AutoUpdate','off');
end

function lim = strip_limits(v)
%% The data range, keeping zero in view for a positive series that comes
%% near it (a flux), tight otherwise (pH, pools).
v = v(isfinite(v));
if isempty(v); lim = [0 1]; return; end
lo = min(v); hi = max(v);
if lo >= 0 && lo < 0.25*hi; lo = 0; end
span = hi-lo;
if span <= 0; span = max(abs(hi)*0.05,1e-9); end
lim = [lo-0.05*span*(lo ~= 0) hi+0.08*span];
end

function mark_cells(ax,G,trace)
if isempty(trace); return; end
colours = lines(numel(trace));
for k = 1:numel(trace)
    c = trace(k);
    plot(ax,G.centre_x(c),G.centre_y(c),'o','MarkerSize',7,'LineWidth',1.5, ...
        'MarkerEdgeColor',colours(k,:),'MarkerFaceColor','w','Tag','TracedCell');
end
end

function text_out = variable_title(entry,words,title_note)
text_out = sprintf('%s, %s %s [%s]%s',entry.label,stat_words(entry.stat),words, ...
    entry.units,title_note);
end

function w = stat_words(stat)
switch stat
    case 'sum';    w = 'total';
    case {'mean','hmean'}; w = 'mean';
    case 'end';    w = 'at the end of the window';
    case 'start';  w = 'at the start of the window';
    case 'change'; w = 'change';
    case 'max';    w = 'maximum';
    case 'min';    w = 'minimum';
    otherwise;     w = stat;
end
end

function text_out = format_value(v,entry)
if ~isfinite(v)
    text_out = 'n/a';
elseif strncmp(entry.units,'pH',2)
    text_out = sprintf('%.2f',v);
elseif abs(v) >= 100
    text_out = sprintf('%.0f',v);
elseif abs(v) >= 10
    text_out = sprintf('%.1f',v);
else
    text_out = sprintf('%.3g',v);
end
end

function nc = columns_for(n)
if n <= 4
    nc = n;
elseif n <= 6
    nc = 3;
elseif n <= 12
    nc = 4;
elseif n <= 20
    nc = 5;
else
    nc = 6;
end
end

function cmap = colour_map(diverging)
if diverging
    anchors = [0.13 0.32 0.62; 0.54 0.72 0.86; 0.97 0.97 0.97; 0.94 0.62 0.48; 0.70 0.12 0.15];
    cmap = interp1(linspace(0,1,5),anchors,linspace(0,1,256));
else
    cmap = parula(256);
end
end

function text = tex_escape(text)
text = strrep(char(text),'_',[char(92) '_']);
end

%% =======================================================================
function gif_file = write_animation(res,G,entry,win,title_note,gif_file,o)
%% One frame per window: the strip with the window shaded, and its map.
fig = figure('Color','w','Visible','off','Name',[entry.key ' animation']);
fig.Position(3:4) = [760 760];
tl = tiledlayout(fig,4,1,'TileSpacing','compact','Padding','compact');
ax_s = nexttile(tl,1);
draw_strip(ax_s,res,entry,win,[]);
yl = ylim(ax_s);
ax = nexttile(tl,2,[3 1]);
cmap = colour_map(entry.diverging);
frame_file = [tempname '.png'];
cleanup = onCleanup(@() delete_if(frame_file,fig));
map = [];
for k = 1:numel(win)
    cla(ax);
    Draw_Field_Cells(ax,G,entry.values(k,:),'ScaleBar',true,'NorthArrow',true);
    colormap(ax,cmap);
    set(ax,'CLim',entry.clim(min(k,size(entry.clim,1)),:));
    if k == 1
        cb = colorbar(ax,'eastoutside');
        cb.Color = 'k';
        cb.Label.String = sprintf('%s [%s]',entry.label,entry.units);
    end
    title(ax,sprintf('%s: %s | mean %s',win(k).label,win(k).short_range, ...
        format_value(entry.field_mean(k),entry)),'Color','k','FontSize',11);
    delete(findobj(ax_s,'Tag','CurrentWindow'));
    patch_x = [win(k).start win(k).stop+caldays(1)];
    r = xregion(ax_s,patch_x(1),patch_x(2),'FaceColor',[0.95 0.55 0.15], ...
        'FaceAlpha',0.25);
    r.Tag = 'CurrentWindow';
    ylim(ax_s,yl);
    title(tl,{tex_escape(res.run_label),sprintf('%s, %s%s',entry.label, ...
        stat_words(entry.stat),title_note)},'Color','k','FontWeight','bold','FontSize',11);
    print(fig,frame_file,'-dpng','-r110');
    rgb = imread(frame_file);
    if isempty(map)
        [ind,map] = rgb2ind(rgb,256,'nodither');
        imwrite(ind,map,gif_file,'gif','LoopCount',Inf,'DelayTime',o.FrameDelay);
    else
        ind = rgb2ind(rgb,map,'nodither');
        imwrite(ind,map,gif_file,'gif','WriteMode','append','DelayTime',o.FrameDelay);
    end
end
fprintf('Saved: %s\n',gif_file);
end

function delete_if(file,fig)
if isfile(file); delete(file); end
if isgraphics(fig); close(fig); end
end

function [png,fig_file] = export_pair(fig,stem,resolution)
png = [stem '.png'];
fig_file = [stem '.fig'];
exportgraphics(fig,png,'Resolution',resolution);
savefig(fig,fig_file);
fprintf('Saved: %s\n',png);
end

%% =======================================================================
function T = cell_table(res,G)
N = res.N_cells;
z = nan(N,1); slope = nan(N,1);
if res.has_topography
    z = res.topo.cell_z_m(:);
    slope = 100*res.topo.cell_slope(:);
end
pft = res.pft(:); rock = res.rock(:);
pft(cellfun(@isempty,pft)) = {'unknown'};
rock(cellfun(@isempty,rock)) = {'unknown'};
T = table((1:N).',res.region(:),pft,rock,G.centre_x(:),G.centre_y(:),z,slope, ...
    res.area_m2(:),'VariableNames',{'cell','zone','vegetation','feedstock', ...
    'x_m','y_m','elevation_m','slope_percent','field_area_m2'});
end

function T = value_table(res,G,V,win,has_reference)
cells = cell_table(res,G);
N = res.N_cells;
n_win = numel(win);
parts = cell(numel(V),1);
names = month_free_dates(win);
for v = 1:numel(V)
    w_idx = repelem((1:n_win).',N,1);
    c_idx = repmat((1:N).',n_win,1);
    values = V(v).values.';
    T = cells(c_idx,:);
    T = addvars(T,repmat(string(V(v).key),height(T),1), ...
        repmat(string(plain(V(v).units)),height(T),1), ...
        repmat(string(V(v).stat),height(T),1),w_idx, ...
        string({win(w_idx).label}).',names.start(w_idx),names.stop(w_idx), ...
        'Before',1,'NewVariableNames',{'variable','units','statistic','window', ...
        'window_label','window_start','window_stop'});
    T.value = values(:);
    fm = V(v).field_mean;
    T.field_mean = fm(w_idx);
    if has_reference
        run_values = V(v).run_values.';
        reference_values = V(v).reference_values.';
        T.run_value = run_values(:);
        T.reference_value = reference_values(:);
    end
    parts{v} = T;
end
T = vertcat(parts{:});
end

function names = month_free_dates(win)
names.start = string(arrayfun(@(w) char(w.start,'yyyy-MM-dd'),win,'UniformOutput',false)).';
names.stop = string(arrayfun(@(w) char(w.stop,'yyyy-MM-dd'),win,'UniformOutput',false)).';
end

function text = plain(text)
text = strrep(strrep(text,'^{',''),'}','');
end
