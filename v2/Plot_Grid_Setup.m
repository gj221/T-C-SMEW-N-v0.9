function [fig,plot_file] = Plot_Grid_Setup(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot_Grid_Setup — 1-D and 2-D computational grid of T&C-SMEW-N v0.9
%%
%% One figure, two panels:
%%
%%   panel A   1-D soil column: 20 layers on the variable-thickness mesh
%%             Zs = [0 10 50 100 150 200 300 400 600 800 1000 1250 1500
%%             1750 2000 2500 3000 3500 4000 4500 5000] mm (5 m deep),
%%             with the active biogeochemical depth Zbio = 300 mm marked
%%             and a top-metre inset.  This column is the domain of the
%%             classic 1-D driver (MAIN_FRAME_BG_ERW_10mins) and of the
%%             vertical_1d group of Run_Matrix_12.
%%
%%   panel B   2-D field domain of Build_Cell_Topology: 60 m x 60 m field,
%%             w = 5 m margin ring, 2 x 2 main grid (25 m cells) plus
%%             4 side strips and 4 corners = 12 cells, with shared-face
%%             lateral coupling drawn between adjacent cell centres.
%%
%% Geometry is read from Run_Matrix_12 and Build_Cell_Topology so the
%% figure always matches the runs.  The layer mesh matches the workbook
%% value (ERW_MOD_PARAM_parameters.xlsx, Kantola rows) and
%% MOD_PARAM_Kantola_Maize_Soy.m.
%%
%% Options
%%   OutputFile       output PNG (companion FIG is saved next to it)
%%                    [this folder / Grid_Setup.png]
%%   Visible          figure visibility          ['on']
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

here = fileparts(mfilename('fullpath'));
opt = struct('OutputFile',fullfile(here,'Grid_Setup.png'), ...
    'Visible','on');
if mod(numel(varargin),2) ~= 0
    error('Plot_Grid_Setup:BadOptions','Options must be name/value pairs.');
end
for k = 1:2:numel(varargin)
    name = char(varargin{k});
    if ~isfield(opt,name)
        error('Plot_Grid_Setup:BadOptions','Unknown option %s.',name);
    end
    opt.(name) = varargin{k+1};
end

addpath(here,fullfile(here,'T&C_Code'));
RM = Run_Matrix_12();
geom = RM(1).geom;
Zbio_mm = RM(1).biogeo_depth_mm;

Zs_mm = [0 10 50 100 150 200 300 400 600 800 1000 1250 1500 1750 ...
    2000 2500 3000 3500 4000 4500 5000];
dz_mm = diff(Zs_mm);
nL = numel(dz_mm);

edge_col = [0.25 0.25 0.25];
soil_fill = {[0.93 0.88 0.80],[1 1 1]};
active_fill = [0.62 0.78 0.95];
main_fill = [0.90 0.95 0.99];
margin_fill = [0.86 0.93 0.83];

fig = figure('Color','w','Visible',char(opt.Visible), ...
    'Name','T&C-SMEW-N grid setup');
fig.Position(3:4) = [1500 800];
tl = tiledlayout(fig,1,2,'TileSpacing','loose','Padding','compact');
title(tl,'T&C-SMEW-N computational grid: 1-D column and 2-D field domain', ...
    'Color','k','FontWeight','bold');

%% ================= panel A: 1-D soil column =================
axA = nexttile(tl);
hold(axA,'on');

rectangle(axA,'Position',[0 0 1 Zbio_mm/1000],'FaceColor',active_fill, ...
    'EdgeColor','none','FaceAlpha',0.30);
for i = 1:nL
    rectangle(axA,'Position',[0 Zs_mm(i)/1000 1 dz_mm(i)/1000], ...
        'FaceColor',soil_fill{mod(i,2)+1},'EdgeColor',edge_col,'LineWidth',1);
end
for i = 1:nL
    if dz_mm(i) >= 150
        text(axA,0.30,(Zs_mm(i)+dz_mm(i)/2)/1000,sprintf('%d',i), ...
            'FontSize',9,'HorizontalAlignment','center', ...
            'VerticalAlignment','middle');
        text(axA,1.08,(Zs_mm(i)+dz_mm(i)/2)/1000,sprintf('%d',dz_mm(i)), ...
            'FontSize',8,'HorizontalAlignment','left', ...
            'VerticalAlignment','middle');
    end
end
text(axA,0.50,Zbio_mm/2000, ...
    sprintf('Z_{bio} = %d mm\n(active biogeochemistry)',Zbio_mm), ...
    'FontSize',8,'HorizontalAlignment','center', ...
    'VerticalAlignment','middle','Color',[0.10 0.30 0.60], ...
    'FontWeight','bold');
rectangle(axA,'Position',[0 0 1 1],'EdgeColor',[0.55 0.60 0.65], ...
    'LineWidth',1,'LineStyle','--');
text(axA,1.02,0.50,'0-1 m detail to the right','FontSize',7, ...
    'Color',[0.4 0.4 0.4],'Rotation',90,'VerticalAlignment','bottom');

axA.YDir = 'reverse';
axA.YTick = 0:1:5;
axA.YTickLabel = {'0','1','2','3','4','5'};
axA.XLim = [0 1.35];
axA.YLim = [0 5];
axA.XTick = [];
axA.Box = 'on';
ylabel(axA,'Depth below surface [m]','FontSize',10);
title(axA,sprintf('1-D soil column: %d layers, 5 m deep',nL), ...
    'FontSize',11);
text(axA,0.02,5.35, ...
    sprintf('Same column in every 2-D cell and in the vertical-1-D runs (%d layers, variable thickness).',nL), ...
    'FontSize',8,'Color',[0.3 0.3 0.3]);

pax = axA.Position;
ins_pos = [pax(1)+0.42*pax(3) pax(2)+0.03*pax(4) 0.30*pax(3) 0.42*pax(4)];
axI = axes('Parent',fig,'Position',ins_pos,'Color','w','Box','on', ...
    'YDir','reverse','XLim',[-0.20 1.10],'YLim',[0 1]);
hold(axI,'on');
rectangle(axI,'Position',[0 0 1 Zbio_mm/1000],'FaceColor',active_fill, ...
    'EdgeColor','none','FaceAlpha',0.30);
for i = 1:9
    rectangle(axI,'Position',[0 Zs_mm(i)/1000 1 dz_mm(i)/1000], ...
        'FaceColor',soil_fill{mod(i,2)+1},'EdgeColor',edge_col,'LineWidth',0.8);
    text(axI,0.45,(Zs_mm(i)+dz_mm(i)/2)/1000,sprintf('%d',i), ...
        'FontSize',7,'HorizontalAlignment','center', ...
        'VerticalAlignment','middle');
    text(axI,0.82,(Zs_mm(i)+dz_mm(i)/2)/1000,sprintf('%d',dz_mm(i)), ...
        'FontSize',6,'HorizontalAlignment','left', ...
        'VerticalAlignment','middle');
end
for i = 1:10
    text(axI,-0.02,Zs_mm(i)/1000,sprintf('%d',Zs_mm(i)), ...
        'FontSize',6,'HorizontalAlignment','right', ...
        'VerticalAlignment','middle');
end
axI.XTick = [];
axI.YTick = [];
title(axI,'Top metre detail (mm)','FontSize',8);

%% ================= panel B: 2-D field domain =================
topo = Build_Cell_Topology(geom.field_side_m,geom.margin_width_m, ...
    geom.main_grid_nx,geom.main_grid_ny);
L = topo.field_side_m;
w = topo.buffer_width_m;

axB = nexttile(tl);
hold(axB,'on');

for c = 1:topo.N_cells
    x0 = topo.cell_x_lo(c); x1 = topo.cell_x_hi(c);
    y0 = topo.cell_y_lo(c); y1 = topo.cell_y_hi(c);
    if strcmp(topo.cell_region{c},'main')
        fill_col = main_fill;
    else
        fill_col = margin_fill;
    end
    rectangle(axB,'Position',[x0 y0 x1-x0 y1-y0],'FaceColor',fill_col, ...
        'EdgeColor',edge_col,'LineWidth',1.2);
end
for e = 1:topo.N_edges
    a = topo.edge_a(e); b = topo.edge_b(e);
    plot(axB,[topo.cell_xc(a) topo.cell_xc(b)], ...
        [topo.cell_yc(a) topo.cell_yc(b)],'Color',[0.70 0.70 0.70], ...
        'LineWidth',1.4);
end
for c = 1:topo.N_cells
    xc = topo.cell_xc(c); yc = topo.cell_yc(c);
    if strcmp(topo.cell_region{c},'main')
        text(axB,xc,yc,sprintf('%d',c),'FontSize',13,'FontWeight','bold', ...
            'HorizontalAlignment','center','VerticalAlignment','middle');
    else
        nm = strrep(strrep(topo.cell_name{c},'Cell_Buf_',''),'Cell_Cor_','');
        text(axB,xc,yc,sprintf('%s\n%d',nm,c),'FontSize',9, ...
            'HorizontalAlignment','center','VerticalAlignment','middle');
    end
end

h_main = patch(axB,NaN,NaN,main_fill,'EdgeColor','none');
h_margin = patch(axB,NaN,NaN,margin_fill,'EdgeColor','none');
h_edge = plot(axB,NaN,NaN,'Color',[0.70 0.70 0.70],'LineWidth',1.4);
legend(axB,[h_main h_margin h_edge], ...
    {'Main cells (1-4)','Margin cells (5-12)','Lateral coupling (shared face)'}, ...
    'Orientation','horizontal','Box','off','Location','northoutside', ...
    'FontSize',9);

x0_ = 0; y_ = 61.5;
plot(axB,[x0_ L],[y_ y_],'k-','LineWidth',0.8);
plot(axB,[x0_ x0_],[y_-0.8 y_+0.8],'k-','LineWidth',0.8);
plot(axB,[L L],[y_-0.8 y_+0.8],'k-','LineWidth',0.8);
text(axB,L/2,y_+0.4,sprintf('%d m',L),'FontSize',9, ...
    'HorizontalAlignment','center','VerticalAlignment','bottom');

x_ = -2.0;
plot(axB,[x_ x_],[0 L],'k-','LineWidth',0.8);
plot(axB,[x_-0.8 x_+0.8],[0 0],'k-','LineWidth',0.8);
plot(axB,[x_-0.8 x_+0.8],[L L],'k-','LineWidth',0.8);
text(axB,x_-0.4,L/2,sprintf('%d m',L),'FontSize',9,'Rotation',90, ...
    'HorizontalAlignment','center','VerticalAlignment','bottom');

y_ = 63.0;
plot(axB,[L-w L],[y_ y_],'k-','LineWidth',0.8);
plot(axB,[L-w L-w],[y_-0.8 y_+0.8],'k-','LineWidth',0.8);
plot(axB,[L L],[y_-0.8 y_+0.8],'k-','LineWidth',0.8);
text(axB,L-w/2,y_+0.4,sprintf('w = %d m margin',w),'FontSize',9, ...
    'HorizontalAlignment','center','VerticalAlignment','bottom');

y_ = -2.0;
plot(axB,[w L-w],[y_ y_],'k-','LineWidth',0.8);
plot(axB,[w w],[y_-0.8 y_+0.8],'k-','LineWidth',0.8);
plot(axB,[L-w L-w],[y_-0.8 y_+0.8],'k-','LineWidth',0.8);
text(axB,L/2,y_-0.7,sprintf('%d m',L-2*w),'FontSize',9, ...
    'HorizontalAlignment','center','VerticalAlignment','top');

y_ = -3.3;
x_cell = w + (L-2*w)/geom.main_grid_nx;
plot(axB,[w x_cell],[y_ y_],'k-','LineWidth',0.8);
plot(axB,[w w],[y_-0.8 y_+0.8],'k-','LineWidth',0.8);
plot(axB,[x_cell x_cell],[y_-0.8 y_+0.8],'k-','LineWidth',0.8);
text(axB,(w+x_cell)/2,y_-0.7,sprintf('%d m',x_cell-w),'FontSize',9, ...
    'HorizontalAlignment','center','VerticalAlignment','top');

text(axB,2.5,62.0,'N','FontSize',13,'FontWeight','bold', ...
    'HorizontalAlignment','center');
text(axB,L/2,-5.0, ...
    sprintf('%d cells = %d main (%d x %d, %d m each) + 8 margin cells (4 side strips + 4 corners).', ...
    topo.N_cells,geom.main_grid_nx*geom.main_grid_ny,geom.main_grid_nx, ...
    geom.main_grid_ny,(L-2*w)/geom.main_grid_nx), ...
    'FontSize',8.5,'Color',[0.3 0.3 0.3],'HorizontalAlignment','center');
text(axB,L/2,-6.2, ...
    'whole_field: amendment on all cells | margin_crop / margin_prairie: margin only | vertical_1d: single column', ...
    'FontSize',8.5,'Color',[0.3 0.3 0.3],'HorizontalAlignment','center');

axB.XLim = [-4.5 64.5];
axB.YLim = [-7.2 65.5];
axB.Box = 'on';
axB.XTick = 0:10:60;
axB.YTick = 0:10:60;
axB.DataAspectRatio = [1 1 1];
axB.XGrid = 'on';
axB.YGrid = 'on';
axB.GridColor = [0.85 0.85 0.85];
xlabel(axB,'Easting [m]','FontSize',10);
ylabel(axB,'Northing [m]','FontSize',10);
title(axB, ...
    sprintf('2-D domain: %d m x %d m field, %d m margin ring, %d coupled cells', ...
    L,L,w,topo.N_cells),'FontSize',11);

%% ================= save =================
exportgraphics(fig,opt.OutputFile,'Resolution',200);
[plot_folder,plot_name] = fileparts(opt.OutputFile);
fig_file = fullfile(plot_folder,[plot_name '.fig']);
savefig(fig,fig_file);
fprintf('Saved: %s\n',opt.OutputFile);
fprintf('Saved: %s\n',fig_file);
plot_file = opt.OutputFile;

end
