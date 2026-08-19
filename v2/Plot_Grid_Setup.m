%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot_Grid_Setup — simple schematic of the T&C-SMEW-N v0.9 grid
%%
%% Standalone script; run it directly (e.g.  run Plot_Grid_Setup).
%%
%% Four panels:
%%
%%   (a)  1-D top-down:    a single 60 m x 60 m column (3600 m2)
%%   (b)  2-D top-down:    60 m x 60 m field, 5 m margin ring, 12 cells
%%                         (4 dark-blue main cells + 8 light-blue buffer
%%                         cells), grid lines shown
%%   (c)  1-D horizontal:  the 20-layer, 5 m soil column - a single
%%                         column has no lateral exchange
%%   (d)  2-D horizontal:  an east-west section of the field - adjacent
%%                         cells exchange water, heat and solutes
%%                         laterally across their shared faces
%%
%% Buffer (margin) cells are light blue and main cells dark blue in every
%% panel.  Outputs: Grid_Setup.png and Grid_Setup.fig next to this file.
%%
%% Geometry matches Run_Matrix_12 and Build_Cell_Topology:
%% L = 60 m field, w = 5 m margin, 2 x 2 main grid (25 m cells).
%% The layer mesh matches ERW_MOD_PARAM_parameters.xlsx (Kantola rows).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(here,'T&C_Code'));

%% geometry
L  = 60;                 % field side [m]
w  = 5;                  % margin width [m]
nx = 2;  ny = 2;         % main grid
Zs_mm = [0 10 50 100 150 200 300 400 600 800 1000 1250 1500 1750 ...
    2000 2500 3000 3500 4000 4500 5000];
Zd = Zs_mm/1000;         % layer boundaries [m]
nL = numel(Zs_mm)-1;     % number of soil layers

%% colours
main_c  = [0.10 0.35 0.65];   % dark blue: main cells
buff_c  = [0.78 0.89 0.98];   % light blue: buffer (margin) cells
edge_c  = [0.15 0.15 0.15];
layer_c = [0.60 0.63 0.65];
arrow_c = [0.85 0.33 0.10];

%% figure
fig = figure('Color','w','Name','T&C-SMEW-N 1-D and 2-D grid setup');
fig.Position(3:4) = [1250 850];
tl = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
title(tl,'T&C-SMEW-N grid setup: 1-D and 2-D domains', ...
    'Color','k','FontWeight','bold','FontSize',13);

%% ================= (a) 1-D top-down =================
ax = nexttile(tl);
hold(ax,'on'); axis(ax,'equal');
rectangle(ax,'Position',[0 0 L L],'FaceColor',main_c, ...
    'EdgeColor',edge_c,'LineWidth',1.5);
text(ax,L/2,L/2,{'1-D single column','60 m x 60 m (3600 m^{2})'}, ...
    'Color','w','FontSize',10,'FontWeight','bold', ...
    'HorizontalAlignment','center');
text(ax,L/2,-3,'60 m','FontSize',9,'HorizontalAlignment','center');
text(ax,-3,L/2,'60 m','FontSize',9,'Rotation',90, ...
    'HorizontalAlignment','center');
ax.XLim = [-8 68]; ax.YLim = [-8 68];
ax.XTick = []; ax.YTick = []; ax.Box = 'on';
title(ax,'(a) 1-D top-down','FontSize',11);

%% ================= (b) 2-D top-down =================
ax = nexttile(tl);
hold(ax,'on'); axis(ax,'equal');

buff_xy = [0 0 L w; 0 L-w L w; 0 0 w L; L-w 0 w L; ...
    0 L-w w w; L-w L-w w w; 0 0 w w; L-w 0 w w];
for k = 1:size(buff_xy,1)
    rectangle(ax,'Position',buff_xy(k,:),'FaceColor',buff_c, ...
        'EdgeColor',edge_c,'LineWidth',1.2);
end
main_xy = [w w L-2*w (L-2*w)/nx; w+(L-2*w)/nx w L-2*w (L-2*w)/nx; ...
    w w+(L-2*w)/ny L-2*w (L-2*w)/ny; w+(L-2*w)/nx w+(L-2*w)/ny L-2*w (L-2*w)/ny];
cell_xc = [w+(L-2*w)/4 w+3*(L-2*w)/4 w+(L-2*w)/4 w+3*(L-2*w)/4];
cell_yc = [w+(L-2*w)/4 w+(L-2*w)/4 w+3*(L-2*w)/4 w+3*(L-2*w)/4];
for k = 1:4
    rectangle(ax,'Position',main_xy(k,:),'FaceColor',main_c, ...
        'EdgeColor','w','LineWidth',1.5);
    text(ax,cell_xc(k),cell_yc(k),sprintf('%d',k),'Color','w', ...
        'FontSize',13,'FontWeight','bold', ...
        'HorizontalAlignment','center');
end

text(ax,L/2,L-2.5,'5 m (margin)','FontSize',8.5, ...
    'HorizontalAlignment','center');
text(ax,L/2,-3,'60 m','FontSize',9,'HorizontalAlignment','center');
text(ax,-3,L/2,'60 m','FontSize',9,'Rotation',90, ...
    'HorizontalAlignment','center');
text(ax,w+(L-2*w)/4,w+2.5,'25 m','FontSize',8.5,'Color','w', ...
    'HorizontalAlignment','center');

ax.XLim = [-8 68]; ax.YLim = [-8 68];
ax.XTick = []; ax.YTick = []; ax.Box = 'on';
title(ax,'(b) 2-D top-down: 12 cells (4 main + 8 buffer)', ...
    'FontSize',11);

%% legend
h_main  = patch(ax,NaN,NaN,main_c,'EdgeColor','none');
h_buff  = patch(ax,NaN,NaN,buff_c,'EdgeColor','none');
h_arrow = plot(ax,NaN,NaN,'Color',arrow_c,'LineWidth',1.8);
lgd = legend([h_main h_buff h_arrow], ...
    {'Main cells','Buffer (margin) cells','Lateral coupling'}, ...
    'Orientation','horizontal','Box','off','FontSize',10);
lgd.Layout.Tile = 'south';

%% ================= (c) 1-D horizontal =================
ax = nexttile(tl);
hold(ax,'on');
rectangle(ax,'Position',[0 -5 1 5],'FaceColor',[0.93 0.95 0.97], ...
    'EdgeColor',edge_c,'LineWidth',1.5);
for i = 2:nL+1
    plot(ax,[0 1],[-Zd(i) -Zd(i)],'Color',layer_c,'LineWidth',0.5);
end
plot(ax,[0 1],[0 0],'Color',edge_c,'LineWidth',2);
text(ax,0.05,-0.25,'soil surface','FontSize',8,'Color',edge_c);
text(ax,0.5,-2.4,sprintf('%d soil layers\n5 m deep',nL), ...
    'FontSize',9.5,'HorizontalAlignment','center', ...
    'Color',[0.1 0.1 0.1]);
plot(ax,[0.45 0.55],[-4.1 -4.5],'Color',[0.78 0.15 0.15],'LineWidth',2);
plot(ax,[0.55 0.45],[-4.1 -4.5],'Color',[0.78 0.15 0.15],'LineWidth',2);
text(ax,0.5,-4.75,'no lateral coupling','FontSize',8.5, ...
    'HorizontalAlignment','center','Color',[0.78 0.15 0.15], ...
    'FontWeight','bold');
ax.XLim = [-0.15 1.25]; ax.YLim = [-5.9 0.6];
ax.XTick = []; ax.YTick = -5:0;
ax.YTickLabel = {'0','1','2','3','4','5'};
ax.Box = 'on';
ylabel(ax,'Depth [m]','FontSize',10);
title(ax,'(c) 1-D horizontal: single column','FontSize',11);

%% ================= (d) 2-D horizontal =================
ax = nexttile(tl);
hold(ax,'on');
col_x = 0:3;
col_fill = {buff_c,main_c,main_c,buff_c};
col_name = {'W margin','Main 1','Main 2','E margin'};
for k = 1:4
    rectangle(ax,'Position',[col_x(k) -5 1 5],'FaceColor',col_fill{k}, ...
        'EdgeColor',edge_c,'LineWidth',1.2);
    for i = 2:nL+1
        plot(ax,col_x(k)+[0 1],[-Zd(i) -Zd(i)],'Color',layer_c, ...
            'LineWidth',0.5);
    end
    text(ax,col_x(k)+0.5,-5.7,col_name{k},'FontSize',8.5, ...
        'HorizontalAlignment','center');
end
plot(ax,[0 4],[0 0],'Color',edge_c,'LineWidth',2);
for k = 1:3
    gx = col_x(k)+1;
    quiver(ax,gx-0.40,-0.25,0.22,0,'Color',arrow_c,'LineWidth',1.6, ...
        'MaxHeadSize',1.2,'AutoScale','off');
    quiver(ax,gx+0.40,-0.25,-0.22,0,'Color',arrow_c,'LineWidth',1.6, ...
        'MaxHeadSize',1.2,'AutoScale','off');
end
text(ax,2,0.75,'lateral coupling (water, heat, solutes)', ...
    'FontSize',9,'HorizontalAlignment','center','Color',arrow_c, ...
    'FontWeight','bold');
text(ax,2,-6.7,'every cell carries the same 20-layer column', ...
    'FontSize',8,'HorizontalAlignment','center','Color',[0.4 0.4 0.4]);
ax.XLim = [-0.35 4.35]; ax.YLim = [-7.4 1.3];
ax.XTick = []; ax.YTick = -5:0;
ax.YTickLabel = {'0','1','2','3','4','5'};
ax.Box = 'on';
ylabel(ax,'Depth [m]','FontSize',10);
title(ax,'(d) 2-D horizontal: lateral exchange between cells', ...
    'FontSize',11);

%% ================= save =================
plot_file = fullfile(here,'Grid_Setup.png');
exportgraphics(fig,plot_file,'Resolution',250);
[plot_folder,plot_name] = fileparts(plot_file);
fig_file = fullfile(plot_folder,[plot_name '.fig']);
savefig(fig,fig_file);
fprintf('Saved: %s\n',plot_file);
fprintf('Saved: %s\n',fig_file);