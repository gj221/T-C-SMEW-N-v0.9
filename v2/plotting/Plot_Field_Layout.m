function [fig,files] = Plot_Field_Layout(source,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot_Field_Layout
%%
%% The cells of a 2-D run as the model sees them:
%%   (a) the grid, its main cells and its margin band, the field outline
%%   (b) the vegetation of each cell (crop rotation or a PFT override)
%%   (c) the feedstock each cell receives
%%   (d) with a DEM: cell elevation, contours of the DEM, the overland
%%       flow paths (Attach_Topography_To_Topology's multiple-flow-
%%       direction weights) and the outlet cells
%%   (e) slope and (f) topographic wetness index
%% (b) and (c) need a model result; (d) to (f) need topography.
%%
%%   Plot_Field_Layout(result_file)
%%   Plot_Field_Layout(topology_file)       a saved field grid
%%   Plot_Field_Layout(topo)
%%   [fig,files] = Plot_Field_Layout(...,'Visible','off')
%%
%% Options
%%   'FigureDirectory'  [v2/figures/field_maps/<result name>]
%%   'Name'             file stem                         ['field_layout']
%%   'Visible'          'on' or 'off'                               ['on']
%%   'Export'           write PNG and FIG                           [true]
%%   'Resolution'       PNG dpi                                      [250]
%%   'Coordinates'      'auto', 'local' or 'projected'            ['auto']
%%   'ClipToField'      draw the in-field part of perimeter squares [true]
%%   'CellNumbers'      true, false or 'auto' (up to 60 cells)    ['auto']
%%   'FlowArrows'       true, false or 'auto' (up to 4000 cells)  ['auto']
%%
%% Colours follow Plot_Grid_Setup: main cells dark blue, margin cells
%% light blue, black text on white.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
p.addParameter('FigureDirectory','');
p.addParameter('Name','field_layout');
p.addParameter('Visible','on',@(x) any(strcmpi(x,{'on','off'})));
p.addParameter('Export',true);
p.addParameter('Resolution',250);
p.addParameter('Coordinates','auto');
p.addParameter('ClipToField',true);
p.addParameter('CellNumbers','auto');
p.addParameter('FlowArrows','auto');
p.parse(varargin{:});
o = p.Results;

here = fileparts(mfilename('fullpath'));
v2 = fileparts(here);
res = Field_Cell_Result(source);
topo = res.topo;
G = Field_Cell_Geometry(topo,'Coordinates',o.Coordinates,'ClipToField',o.ClipToField);
main_c = [0.10 0.35 0.65];
buff_c = [0.82 0.91 0.98];

has_config = ~isempty(res.ColCfg);
has_dem = res.has_topography;
has_twi = has_dem && isfield(topo,'cell_twi') && any(isfinite(topo.cell_twi));
panels = {'zones'};
if has_config; panels = [panels {'pft','rock'}]; end
if has_dem; panels = [panels {'elevation','slope'}]; end
if has_twi; panels = [panels {'twi'}]; end
n_panel = numel(panels);
n_col = min(3,n_panel);
n_row = ceil(n_panel/n_col);

fig = figure('Color','w','Visible',o.Visible,'Name','Field layout');
fig.Position(3:4) = [430*n_col 400*n_row+70];
tl = tiledlayout(fig,n_row,n_col,'TileSpacing','compact','Padding','compact');
title(tl,{tex_escape(res.run_label), grid_summary(res,G)}, ...
    'Color','k','FontWeight','bold','FontSize',12);

letters = 'abcdef';
for k = 1:n_panel
    ax = nexttile(tl);
    switch panels{k}
        case 'zones'
            values = 2-double(res.is_main);
            Draw_Field_Cells(ax,G,values,'Axes','on','ScaleBar',false,'NorthArrow',true);
            colormap(ax,[main_c; buff_c]);
            set(ax,'CLim',[0.5 2.5]);
            draw_cell_numbers(ax,G,res,o.CellNumbers);
            h = [patch(ax,NaN,NaN,main_c) patch(ax,NaN,NaN,buff_c)];
            names = {'Main cell','Margin cell'};
            if G.clipped
                h(end+1) = patch(ax,NaN,NaN,[0.96 0.96 0.96],'EdgeColor',[0.78 0.78 0.78]); %#ok<AGROW>
                names{end+1} = 'Grid square outside the field'; %#ok<AGROW>
            end
            h(end+1) = plot(ax,NaN,NaN,'k-','LineWidth',1.4); %#ok<AGROW>
            names{end+1} = 'Field boundary'; %#ok<AGROW>
            if ~isempty(G.zone_x)
                h(end+1) = plot(ax,NaN,NaN,'k--','LineWidth',0.9); %#ok<AGROW>
                names{end+1} = 'Main / margin boundary'; %#ok<AGROW>
            end
            place_legend(ax,h,names);
            title(ax,sprintf('(%s) Grid cells and margin zone',letters(k)),'Color','k');
        case 'pft'
            categorical_panel(ax,G,res.pft,@pft_colour);
            title(ax,sprintf('(%s) Vegetation (PFT)',letters(k)),'Color','k');
        case 'rock'
            categorical_panel(ax,G,res.rock,@rock_colour);
            title(ax,sprintf('(%s) Feedstock applied',letters(k)),'Color','k');
        case 'elevation'
            z = topo.cell_z_m(:);
            Draw_Field_Cells(ax,G,z);
            colormap(ax,terrain_map(256));
            set(ax,'CLim',limits(z));
            draw_contours(ax,G,topo);
            n_arrow = draw_flow(ax,G,topo,o.FlowArrows);
            outlet = find(topo.cell_is_outlet(:));
            h_out = plot(ax,G.centre_x(outlet),G.centre_y(outlet),'o','MarkerSize',6, ...
                'MarkerFaceColor',[0.85 0.10 0.10],'MarkerEdgeColor','k','Tag','Outlet');
            cb = colorbar(ax,'eastoutside');
            cb.Label.String = 'Elevation [m]';
            cb.Color = 'k';
            names = {sprintf('Outlet cell (%d)',numel(outlet))};
            h = h_out;
            if n_arrow > 0
                h(end+1) = plot(ax,NaN,NaN,'-','Color',[0.1 0.1 0.1],'LineWidth',1); %#ok<AGROW>
                names{end+1} = 'Overland flow direction'; %#ok<AGROW>
            end
            place_legend(ax,h,names);
            title(ax,sprintf('(%s) Elevation and overland flow (relief %.2f m)', ...
                letters(k),max(z)-min(z)),'Color','k');
        case 'slope'
            s = 100*topo.cell_slope(:);
            Draw_Field_Cells(ax,G,s);
            colormap(ax,sequential_map(256,[1.00 0.97 0.80; 0.99 0.73 0.33; 0.80 0.36 0.07; 0.40 0.15 0.03]));
            set(ax,'CLim',limits(s));
            cb = colorbar(ax,'eastoutside');
            cb.Label.String = 'Slope [%]';
            cb.Color = 'k';
            title(ax,sprintf('(%s) Slope (mean %.2f %%)',letters(k),mean(s)),'Color','k');
        case 'twi'
            t = topo.cell_twi(:);
            Draw_Field_Cells(ax,G,t);
            colormap(ax,sequential_map(256,[0.97 0.98 1.00; 0.62 0.79 0.88; 0.19 0.51 0.74; 0.03 0.19 0.42]));
            set(ax,'CLim',limits(t));
            cb = colorbar(ax,'eastoutside');
            cb.Label.String = 'TWI [-]';
            cb.Color = 'k';
            title(ax,sprintf('(%s) Topographic wetness index',letters(k)),'Color','k');
    end
end

files = struct('png','','fig','');
if o.Export
    figure_dir = o.FigureDirectory;
    if isempty(figure_dir)
        figure_dir = fullfile(v2,'figures','field_maps',res.short_name);
    end
    if ~isfolder(figure_dir); mkdir(figure_dir); end
    files.png = fullfile(figure_dir,[o.Name '.png']);
    files.fig = fullfile(figure_dir,[o.Name '.fig']);
    exportgraphics(fig,files.png,'Resolution',o.Resolution);
    savefig(fig,files.fig);
    fprintf('Saved: %s\n',files.png);
end
end

%% =======================================================================
function text = grid_summary(res,G)
topo = res.topo;
n_main = sum(res.is_main);
if isfield(topo,'cell_size_m') && ~isempty(topo.cell_size_m)
    text = sprintf('%d cells of %g m (%d main, %d margin)', ...
        res.N_cells,topo.cell_size_m,n_main,res.N_cells-n_main);
else
    text = sprintf('%d cells (%d main, %d margin)',res.N_cells,n_main,res.N_cells-n_main);
end
text = sprintf('%s, field %.0f m^2',text,sum(res.area_m2));
if res.has_topography
    text = [text ', DEM attached'];
else
    text = [text ', flat'];
end
if strcmp(G.coordinates,'projected') && isfield(topo,'grid_rotation_degrees')
    text = sprintf('%s, grid rotated %.1f deg',text,topo.grid_rotation_degrees);
end
end

function text = tex_escape(text)
text = strrep(char(text),'_',[char(92) '_']);
end

%% =======================================================================
function categorical_panel(ax,G,labels,colour_fn)
[names,~,idx] = unique(labels(:),'stable');
colours = zeros(numel(names),3);
for k = 1:numel(names)
    colours(k,:) = colour_fn(names{k},k);
end
Draw_Field_Cells(ax,G,idx);
colormap(ax,colours);
set(ax,'CLim',[0.5 numel(names)+0.5]);
h = gobjects(numel(names),1);
for k = 1:numel(names)
    h(k) = patch(ax,NaN,NaN,colours(k,:),'EdgeColor',[0.25 0.25 0.25]);
end
counts = accumarray(idx,1);
legend_names = arrayfun(@(k) sprintf('%s (%d)',names{k},counts(k)),1:numel(names), ...
    'UniformOutput',false);
place_legend(ax,h,legend_names);
end

function c = pft_colour(name,k)
if startsWith(name,'Crop')
    c = [0.93 0.76 0.29];
elseif contains(lower(name),'prairie') || contains(lower(name),'grass')
    c = [0.36 0.62 0.29];
else
    palette = [0.55 0.35 0.65; 0.80 0.45 0.45; 0.40 0.70 0.75; 0.60 0.60 0.30];
    c = palette(mod(k-1,size(palette,1))+1,:);
end
end

function c = rock_colour(name,~)
switch lower(name)
    case 'none';    c = [0.97 0.97 0.97];
    case 'basalt';  c = [0.28 0.28 0.31];
    case 'calcite'; c = [0.87 0.80 0.62];
    otherwise;      c = [0.60 0.60 0.60];
end
end

function place_legend(ax,h,names)
lgd = legend(ax,h,names,'Location','southoutside','Orientation','vertical', ...
    'Box','off','TextColor','k','FontSize',8,'AutoUpdate','off');
if numel(names) > 3
    lgd.NumColumns = 2;
end
end

%% =======================================================================
function draw_cell_numbers(ax,G,res,mode)
if ischar(mode) && strcmp(mode,'auto')
    mode = res.N_cells <= 60;
end
if ~mode; return; end
for c = 1:res.N_cells
    if res.is_main(c); colour = 'w'; else; colour = 'k'; end
    text(ax,G.centre_x(c),G.centre_y(c),sprintf('%d',c),'Color',colour, ...
        'FontSize',7,'HorizontalAlignment','center','Clipping','on');
end
end

function draw_contours(ax,G,topo)
if ~isfield(topo,'topography_grid') || ~isstruct(topo.topography_grid) || ...
        ~all(isfield(topo.topography_grid,{'x','y','Z'}))
    return
end
grid = topo.topography_grid;
z = grid.Z;
if all(~isfinite(z(:))); return; end
relief = max(z(:))-min(z(:));
step = nice_step(relief/10);
levels = ceil(min(z(:))/step)*step:step:max(z(:));
if numel(levels) < 2; return; end
C = contourc(double(grid.x(:).'),double(grid.y(:).'),double(z),levels);
k = 1;
xs = []; ys = [];
while k < size(C,2)
    n = C(2,k);
    xs = [xs C(1,k+1:k+n) NaN]; %#ok<AGROW>
    ys = [ys C(2,k+1:k+n) NaN]; %#ok<AGROW>
    k = k+n+1;
end
[px,py] = G.dem_to_plot(xs,ys);
outside = ~inpolygon(px,py,G.outline_x,G.outline_y);
px(outside) = NaN; py(outside) = NaN;           %% the field only
%% Opaque and lifted above the cell faces: a translucent line is depth-
%% sorted behind the faces and vanishes where it crosses them.
plot3(ax,px,py,ones(size(px)),'-','Color',[0.30 0.30 0.30],'LineWidth',0.5, ...
    'Tag','DEMContours');
end

function n_arrow = draw_flow(ax,G,topo,mode)
n_arrow = 0;
if ~isfield(topo,'edge_overland_weight_ab'); return; end
if ischar(mode) && strcmp(mode,'auto')
    mode = topo.N_cells <= 4000;
end
if ~mode; return; end
%% One arrow per cell: the weight-averaged direction of its receivers (the
%% multiple-flow-direction split), 0.4 of a square long when all of the
%% cell's runoff leaves in one direction, shorter when it spreads.
a = topo.edge_a(:); b = topo.edge_b(:);
w_ab = topo.edge_overland_weight_ab(:); w_ba = topo.edge_overland_weight_ba(:);
from = [a; b]; to = [b; a]; weight = [w_ab; w_ba];
keep = weight > 0;
from = from(keep); to = to(keep); weight = weight(keep);
if isempty(from); return; end
dx = G.centre_x(to)-G.centre_x(from);
dy = G.centre_y(to)-G.centre_y(from);
len = hypot(dx,dy);
N = topo.N_cells;
ux = accumarray(from,weight.*dx./len,[N 1]);
uy = accumarray(from,weight.*dy./len,[N 1]);
moving = hypot(ux,uy) > 0.05;
L = 0.4*G.cell_size_m;
quiver(ax,G.centre_x(moving),G.centre_y(moving),L*ux(moving),L*uy(moving),0, ...
    'Color',[0.1 0.1 0.1],'LineWidth',0.9,'MaxHeadSize',0.6,'Tag','OverlandFlow');
n_arrow = sum(moving);
end

%% =======================================================================
function lim = limits(v)
v = v(isfinite(v));
if isempty(v); lim = [0 1]; return; end
lim = [min(v) max(v)];
if lim(2) <= lim(1)
    lim = lim + [-0.5 0.5]*max(abs(lim(1))*1e-3,1e-6);
end
end

function s = nice_step(raw)
if ~(raw > 0); s = 1; return; end
e = 10^floor(log10(raw));
f = raw/e;
choices = [1 2 2.5 5 10];
s = choices(find(choices >= f,1))*e;
end

function cmap = terrain_map(n)
anchors = [0.20 0.45 0.30; 0.55 0.72 0.42; 0.90 0.86 0.58; 0.72 0.55 0.37; 0.50 0.38 0.30];
cmap = interp1(linspace(0,1,size(anchors,1)),anchors,linspace(0,1,n));
end

function cmap = sequential_map(n,anchors)
cmap = interp1(linspace(0,1,size(anchors,1)),anchors,linspace(0,1,n));
end
