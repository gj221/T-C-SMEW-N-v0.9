%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot_Grid_Setup â€” simple schematic of the T&C-SMEW-N v0.9 grid
%%
%% Standalone script; run it directly (e.g.  run Plot_Grid_Setup).
%%
%% Four panels, white backgrounds, black axis and legend text:
%%
%%   (a)  1-D top-down:    a single 60 m x 60 m column (3600 m2)
%%   (b)  2-D top-down:    60 m x 60 m field, 5 m margin ring, 12 cells
%%                         (4 dark-blue main cells + 8 light-blue buffer
%%                         cells); orange arrows show the lateral
%%                         coupling across every shared face
%%   (c)  1-D horizontal:  the 20-layer, 5 m soil column - light brown
%%                         active zone (0-300 mm) over darker brown
%%                         deeper soil, layer lines kept visible; a
%%                         single column has no lateral transport
%%   (d)  2-D horizontal:  an east-west section of the field - light
%%                         brown active zone over darker brown soil in
%%                         every column, blue caps (dark blue main,
%%                         light blue margin); red lateral transport
%%                         arrows run through the 0-300 mm active zone
%%
%% Buffer (margin) cells are light blue and main cells dark blue in every
%% panel.  Cell geometry is taken from Build_Cell_Topology(L, w, nx, ny)
%% so the top-down drawing always matches the runs (L = 60 m, w = 5 m,
%% 2 x 2 main grid of 25 m cells).  The layer mesh matches
%% ERW_MOD_PARAM_parameters.xlsx (Kantola rows).
%%
%% Outputs: Grid_Setup.png and Grid_Setup.fig next to this file.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

here = fileparts(mfilename('fullpath'));
addpath(here,fullfile(here,'T&C_Code'));

%% geometry
L  = 60;                 % field side [m]
w  = 5;                  % margin width [m]
nx = 2;  ny = 2;         % main grid
Zs_mm = [0 10 50 100 150 200 300 400 600 800 1000 1250 1500 1750 ...
    2000 2500 3000 3500 4000 4500 5000];
Zd = Zs_mm/1000;         % layer boundaries [m]
nL = numel(Zs_mm)-1;     % number of soil layers
topo = Build_Cell_Topology(L,w,nx,ny);

%% colours
main_c  = [0.10 0.35 0.65];   % dark blue: main cells
buff_c  = [0.82 0.91 0.98];   % light blue: buffer (margin) cells
edge_c  = [0.15 0.15 0.15];
soil_light_c = [0.86 0.75 0.58];  % light brown: active zone (0-300 mm)
soil_deep_c  = [0.70 0.56 0.38];  % darker brown: soil below the active zone
soil_alt_c   = [0.75 0.62 0.44];  % alternating shade in the deep zone
soil_line_c  = [0.42 0.30 0.18];  % brown layer lines
arrow_c      = [0.85 0.10 0.10];  % red: lateral coupling
kw_txt = {'Color','k'};           % black text
kw_bg  = {'BackgroundColor','w','EdgeColor','none'};  % text on white box

%% figure
fig = figure('Color','w','Name','T&C-SMEW-N 1-D and 2-D grid setup');
fig.Position(3:4) = [1400 1400];
tl = tiledlayout(fig,4,2,'TileSpacing','tight','Padding','compact');
title(tl,'T&C-SMEW-N grid setup: 1-D and 2-D domains', ...
    'Color','k','FontWeight','bold','FontSize',13);

%% ================= (a) 1-D top-down =================
ax = nexttile(tl);
hold(ax,'on'); axis(ax,'equal'); ax.Color = 'w';
rectangle(ax,'Position',[0 0 L L],'FaceColor',main_c, ...
    'EdgeColor',edge_c,'LineWidth',1.5);
text(ax,L/2,L/2,{'1-D single column','60 m x 60 m (3600 m^{2})'}, ...
    'Color','w','FontSize',10,'FontWeight','bold', ...
    'HorizontalAlignment','center');
plot(ax,[0 L],[-4.5 -4.5],'Color','k','LineWidth',1);
plot(ax,[0 0],[-5.3 -3.7],'Color','k','LineWidth',1);
plot(ax,[L L],[-5.3 -3.7],'Color','k','LineWidth',1);
text(ax,L/2,-5.6,'60 m','FontSize',9,'HorizontalAlignment','center',kw_txt{:});
plot(ax,[-4.5 -4.5],[0 L],'Color','k','LineWidth',1);
plot(ax,[-5.3 -3.7],[0 0],'Color','k','LineWidth',1);
plot(ax,[-5.3 -3.7],[L L],'Color','k','LineWidth',1);
text(ax,-6.5,L/2,'60 m','FontSize',9,'Rotation',90, ...
    'HorizontalAlignment','center',kw_txt{:});
ax.XLim = [-9.5 69.5]; ax.YLim = [-9.5 69.5];
ax.XTick = []; ax.YTick = [];
ax.XColor = 'k'; ax.YColor = 'k';
ax.Box = 'on';
title(ax,'(a) 1-D top-down: single column','FontSize',11,'Color','k');

%% ================= (b) 2-D top-down =================
ax = nexttile(tl);
hold(ax,'on'); axis(ax,'equal'); ax.Color = 'w';

for c = 1:topo.N_cells
    if strcmp(topo.cell_region{c},'main')
        fill_c = main_c;
    else
        fill_c = buff_c;
    end
    rectangle(ax,'Position',[topo.cell_x_lo(c) topo.cell_y_lo(c) ...
        topo.cell_dx(c) topo.cell_dy(c)],'FaceColor',fill_c, ...
        'EdgeColor',edge_c,'LineWidth',1.2);
end
tol = 1e-9*L;
for e = 1:topo.N_edges
    a = topo.edge_a(e); b = topo.edge_b(e);
    if abs(topo.cell_x_hi(a)-topo.cell_x_lo(b)) < tol || ...
            abs(topo.cell_x_hi(b)-topo.cell_x_lo(a)) < tol
        mx = min(topo.cell_x_hi(a),topo.cell_x_hi(b));
        my = 0.5*(max(topo.cell_y_lo(a),topo.cell_y_lo(b)) + ...
            min(topo.cell_y_hi(a),topo.cell_y_hi(b)));
        quiver(ax,mx-1.6,my,1.6,0,'Color',arrow_c,'LineWidth',2.5, ...
            'MaxHeadSize',0.5,'AutoScale','off');
        quiver(ax,mx+1.6,my,-1.6,0,'Color',arrow_c,'LineWidth',2.5, ...
            'MaxHeadSize',0.5,'AutoScale','off');
    else
        my = min(topo.cell_y_hi(a),topo.cell_y_hi(b));
        mx = 0.5*(max(topo.cell_x_lo(a),topo.cell_x_lo(b)) + ...
            min(topo.cell_x_hi(a),topo.cell_x_hi(b)));
        quiver(ax,mx,my-1.6,0,1.6,'Color',arrow_c,'LineWidth',2.5, ...
            'MaxHeadSize',0.5,'AutoScale','off');
        quiver(ax,mx,my+1.6,0,-1.6,'Color',arrow_c,'LineWidth',2.5, ...
            'MaxHeadSize',0.5,'AutoScale','off');
    end
end
for c = 1:topo.N_cells
    if strcmp(topo.cell_region{c},'main')
        text(ax,topo.cell_xc(c),topo.cell_yc(c),sprintf('%d',c), ...
            'Color','w','FontSize',13,'FontWeight','bold', ...
            'HorizontalAlignment','center');
    end
end

plot(ax,[L-w L],[-4.5 -4.5],'Color','k','LineWidth',1);
plot(ax,[L-w L-w],[-5.3 -3.7],'Color','k','LineWidth',1);
plot(ax,[L L],[-5.3 -3.7],'Color','k','LineWidth',1);
text(ax,L-2.5,-5.8,'5 m margin','FontSize',8.5, ...
    'HorizontalAlignment','center',kw_txt{:});
plot(ax,[0 L],[-6.7 -6.7],'Color','k','LineWidth',1);
plot(ax,[0 0],[-7.5 -5.9],'Color','k','LineWidth',1);
plot(ax,[L L],[-7.5 -5.9],'Color','k','LineWidth',1);
text(ax,L/2,-7.8,'60 m','FontSize',9,'HorizontalAlignment','center',kw_txt{:});
plot(ax,[-4.5 -4.5],[0 L],'Color','k','LineWidth',1);
plot(ax,[-5.3 -3.7],[0 0],'Color','k','LineWidth',1);
plot(ax,[-5.3 -3.7],[L L],'Color','k','LineWidth',1);
text(ax,-6.5,L/2,'60 m','FontSize',9,'Rotation',90, ...
    'HorizontalAlignment','center',kw_txt{:});
plot(ax,[w w+(L-2*w)/nx],[w+7.0 w+7.0],'Color','w','LineWidth',1);
plot(ax,[w w],[w+6.2 w+7.8],'Color','w','LineWidth',1);
plot(ax,[w+(L-2*w)/nx w+(L-2*w)/nx],[w+6.2 w+7.8],'Color','w','LineWidth',1);
text(ax,w+(L-2*w)/(2*nx),w+6.2,'25 m','FontSize',8.5,'Color','w', ...
    'HorizontalAlignment','center');

ax.XLim = [-9.5 69.5]; ax.YLim = [-9.5 69.5];
ax.XTick = []; ax.YTick = [];
ax.XColor = 'k'; ax.YColor = 'k';
ax.Box = 'on';
title(ax,'(b) 2-D top-down: 12 cells, coupling across shared faces', ...
    'FontSize',11,'Color','k');

%% legend
h_main  = patch(ax,NaN,NaN,main_c,'EdgeColor','none');
h_buff  = patch(ax,NaN,NaN,buff_c,'EdgeColor','none');
h_soil  = patch(ax,NaN,NaN,soil_deep_c,'EdgeColor',soil_line_c);
h_zone  = patch(ax,NaN,NaN,soil_light_c,'EdgeColor',soil_line_c);
h_arrow = plot(ax,NaN,NaN,'Color',arrow_c,'LineWidth',2);
lgd = legend([h_main h_buff h_soil h_zone h_arrow], ...
    {'Main cell','Margin cell','Soil column (20 layers)','Active zone (0-300 mm)','Lateral coupling'}, ...
    'Orientation','horizontal','Box','off','FontSize',9, ...
    'TextColor','k','Color','w','EdgeColor','none');
lgd.Layout.Tile = 'south';

%% ================= (c) 1-D horizontal =================
ax = nexttile(tl,[3 1]);
hold(ax,'on'); ax.Color = 'w';
for i = 1:nL
    if Zd(i+1) <= 0.3
        fill_c = soil_light_c;
    else
        fill_c = soil_deep_c;
        if mod(i,2)==0; fill_c = soil_alt_c; end
    end
    rectangle(ax,'Position',[0 -Zd(i) 1 Zd(i+1)-Zd(i)], ...
        'FaceColor',fill_c,'EdgeColor',soil_line_c,'LineWidth',0.6);
end
plot(ax,[0 1],[-0.3 -0.3],'--','Color',soil_line_c,'LineWidth',1);
plot(ax,[0 1],[0 0],'Color',edge_c,'LineWidth',2);
text(ax,0.03,0.12,'soil surface','FontSize',8, ...
    'HorizontalAlignment','left',kw_txt{:},kw_bg{:});
text(ax,0.5,-0.15,{'active zone','0-300 mm'},'FontSize',8, ...
    'HorizontalAlignment','center','VerticalAlignment','middle', ...
    kw_txt{:},kw_bg{:});
text(ax,0.5,-2.4,{'20 soil layers','5 m deep'}, ...
    'FontSize',9.5,'HorizontalAlignment','center','FontWeight','bold', ...
    kw_txt{:},kw_bg{:});
plot(ax,[0.45 0.55],[-4.2 -4.6],'Color',[0.78 0.15 0.15],'LineWidth',2);
plot(ax,[0.55 0.45],[-4.2 -4.6],'Color',[0.78 0.15 0.15],'LineWidth',2);
text(ax,0.5,-4.85,'no lateral transport','FontSize',8.5, ...
    'HorizontalAlignment','center','FontWeight','bold',kw_txt{:},kw_bg{:});
ax.XLim = [-0.05 1.25]; ax.YLim = [-5.15 0.5];
ax.XTick = []; ax.YTick = -5:0;
ax.YTickLabel = {'0','1','2','3','4','5'};
ax.XColor = 'k'; ax.YColor = 'k';
ax.Box = 'on';
ylabel(ax,'Depth [m]','FontSize',10,'Color','k');
title(ax,'(c) 1-D horizontal: single soil column','FontSize',11,'Color','k');

%% ================= (d) 2-D horizontal =================
ax = nexttile(tl,[3 1]);
hold(ax,'on'); ax.Color = 'w';
col_x = 0:3;
col_cap = {buff_c,main_c,main_c,buff_c};
col_name = {'W margin','Main 1','Main 2','E margin'};
for k = 1:4
    for i = 1:nL
        if Zd(i+1) <= 0.3
            fill_c = soil_light_c;
        else
            fill_c = soil_deep_c;
            if mod(i,2)==0; fill_c = soil_alt_c; end
        end
        rectangle(ax,'Position',[col_x(k) -Zd(i) 1 Zd(i+1)-Zd(i)], ...
            'FaceColor',fill_c,'EdgeColor',soil_line_c,'LineWidth',0.6);
    end
    text(ax,col_x(k)+0.5,-5.6,col_name{k},'FontSize',8.5, ...
        'HorizontalAlignment','center',kw_txt{:},kw_bg{:});
end
plot(ax,[0 4],[-0.3 -0.3],'--','Color',soil_line_c,'LineWidth',1);
for k = 1:4
    rectangle(ax,'Position',[col_x(k) -0.12 1 0.12], ...
        'FaceColor',col_cap{k},'EdgeColor',edge_c,'LineWidth',1);
end
plot(ax,[0 4],[0 0],'Color',edge_c,'LineWidth',2);
text(ax,0.03,0.12,'soil surface','FontSize',8, ...
    'HorizontalAlignment','left',kw_txt{:},kw_bg{:});
for k = 1:3
    gx = col_x(k)+1;
    quiver(ax,gx-0.42,-0.18,0.24,0,'Color',arrow_c,'LineWidth',2.5, ...
        'MaxHeadSize',0.8,'AutoScale','off');
    quiver(ax,gx+0.42,-0.18,-0.24,0,'Color',arrow_c,'LineWidth',2.5, ...
        'MaxHeadSize',0.8,'AutoScale','off');
end
text(ax,2,1.0,'lateral transport between adjacent cells', ...
    'FontSize',10,'HorizontalAlignment','center','FontWeight','bold', ...
    kw_txt{:},kw_bg{:});
text(ax,2,0.5,'water, heat and solutes','FontSize',8.5, ...
    'HorizontalAlignment','center',kw_txt{:},kw_bg{:});
plot(ax,[4.16 4.16],[0 -0.3],'Color','k','LineWidth',1.2);
plot(ax,[4.11 4.21],[0 0],'Color','k','LineWidth',1.2);
plot(ax,[4.11 4.21],[-0.3 -0.3],'Color','k','LineWidth',1.2);
text(ax,4.24,-0.15,{'active zone','0-300 mm'},'FontSize',8, ...
    'HorizontalAlignment','left','VerticalAlignment','middle', ...
    kw_txt{:},kw_bg{:});
text(ax,2,-6.15,'every cell carries the same 20-layer soil column', ...
    'FontSize',8,'HorizontalAlignment','center',kw_txt{:},kw_bg{:});
ax.XLim = [-0.3 5.15]; ax.YLim = [-6.3 1.5];
ax.XTick = []; ax.YTick = -5:0;
ax.YTickLabel = {'0','1','2','3','4','5'};
ax.XColor = 'k'; ax.YColor = 'k';
ax.Box = 'on';
ylabel(ax,'Depth [m]','FontSize',10,'Color','k');
title(ax,'(d) 2-D horizontal: lateral transport through the active zone', ...
    'FontSize',11,'Color','k');

%% ================= save =================
plot_file = fullfile(here,'Grid_Setup.png');
exportgraphics(fig,plot_file,'Resolution',250);
[plot_folder,plot_name] = fileparts(plot_file);
fig_file = fullfile(plot_folder,[plot_name '.fig']);
savefig(fig,fig_file);
fprintf('Saved: %s\n',plot_file);
fprintf('Saved: %s\n',fig_file);