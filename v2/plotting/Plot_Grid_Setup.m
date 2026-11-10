%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot_Grid_Setup and schematic of the T&C-SMEW-N v0.9 1-D and 2-D grids
%%
%% Standalone script; run it directly (e.g.  run Plot_Grid_Setup).
%%
%% Publication figure, 18 cm wide, Arial 6-8 pt:
%%
%%   (a)  1-D plan view:   one 60 m x 60 m column (3600 m2) on a planar
%%                         slope, contours every 0.2 m, section A-A'
%%   (b)  2-D plan view:   the same plot as 12 cells (4 main cells of
%%                         25 m and 8 margin cells in the 5 m ring);
%%                         arrows mark the lateral exchange across each
%%                         of the 20 shared faces
%%   (c)  1-D section A-A': the 20-layer, 5 m column under the slope,
%%                         closed sides and free drainage at the base;
%%                         the inset draws the top 0.4 m of the mesh and
%%                         the 0-300 mm active zone at true scale
%%   (d)  2-D section A-A': cells 5, 4, 2 and 6 at their true widths;
%%                         water and heat cross each face in every layer,
%%                         aqueous solutes and soil CO2 in the active zone
%%
%% The slope is the mean slope of the Kantola basalt field (2.3 %,
%% runtime/Field_2D/README.md; the control field is 1.0 %), drawn as a
%% plane for illustration: the v0.9 runs use no DEM (OPT_DEM = 0) and
%% Slo_top = 0.  Sections carry a vertical exaggeration VE so the slope
%% and the active zone are visible; tick labels give true elevation.
%%
%% Cell geometry is taken from Build_Cell_Topology(L, w, nx, ny), as in
%% Run_Matrix_12 (L = 60 m, w = 5 m, 2 x 2 main grid of 25 m cells).  The
%% layer mesh matches ERW_MOD_PARAM_parameters.xlsx (Kantola rows,
%% Zbio = 300 mm).
%%
%% Outputs: Grid_Setup.png (600 dpi), Grid_Setup.pdf (vector) and
%% Grid_Setup.fig in v2/figures.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% This script lives in v2/plotting; the repository root is one level up.
here = fileparts(mfilename('fullpath'));
if isempty(here); here = pwd; end
root = fileparts(here);
addpath(fullfile(root,'T&C_Code'));
figure_dir = fullfile(root,'figures');
if ~isfolder(figure_dir); mkdir(figure_dir); end

%% geometry
L  = 60;                 % field side [m]
w  = 5;                  % margin width [m]
nx = 2;  ny = 2;         % main grid
slope = 0.023;           % planar slope [m/m], Kantola basalt field mean
VE = 5;                  % vertical exaggeration of the sections [-]
x_sec = 47;              % plan position of section A-A' [m]
Zs_mm = [0 10 50 100 150 200 300 400 600 800 1000 1250 1500 1750 ...
    2000 2500 3000 3500 4000 4500 5000];
Zd = Zs_mm/1000;         % layer boundaries [m]
Zbio = 0.3;              % active biogeochemical depth [m]
topo = Build_Cell_Topology(L,w,nx,ny);
z_surf = @(s) slope*(L - s);          % surface along A-A', s = 0 at A (upslope)
sec_cells = find(topo.cell_x_lo(:) <= x_sec & topo.cell_x_hi(:) > x_sec);
[~,o] = sort(topo.cell_yc(sec_cells),'descend');
sec_cells = sec_cells(o).';          % cells crossed by A-A', upslope first
s_face = L - [topo.cell_y_hi(sec_cells(:)).' topo.cell_y_lo(sec_cells(end))];  % [m]
z_cont = 0.2:0.2:slope*L;             % contour levels [m]
pr = '′'; nd = '–';             % prime and en dash
sec_name = ['A' nd 'A' pr];

%% colours
main_c  = [0.21 0.45 0.70];   % blue: main cells
buff_c  = [0.80 0.89 0.97];   % light blue: margin cells
edge_c  = [0.15 0.15 0.15];
cont_c  = [0.35 0.35 0.35];   % contours
soil_light_c = [0.88 0.78 0.60];  % light brown: active zone (0-300 mm)
soil_deep_c  = [0.72 0.58 0.40];  % darker brown: soil below the active zone
soil_alt_c   = [0.77 0.64 0.46];  % alternating shade in the deep zone
soil_line_c  = [0.45 0.33 0.20];  % brown layer lines
water_c  = [0.84 0.27 0.00];      % vermillion: water and heat, every layer
solute_c = [0.00 0.55 0.40];      % bluish green: solutes and CO2, active zone
drain_c  = [0.30 0.30 0.30];
font = {'FontName','Arial','FontSize',7};

%% figure (cm); axes are placed by hand so the equal-aspect panels fill them
fig_w = 18; fig_h = 16.6;
fig = figure('Color','w','Name','T&C-SMEW-N 1-D and 2-D grid setup', ...
    'Units','centimeters','Position',[2 2 fig_w fig_h],'PaperPositionMode','auto');
pos_cm = @(x,y,wd,ht) [x/fig_w y/fig_h wd/fig_w ht/fig_h];
plan_xlim = [-9 67]; plan_ylim = [-9 68];
sec_xlim = [-2.5 62.5]; sec_ylim = VE*[-5.6 slope*L+0.75];

%% ================= (a) 1-D plan view =================
ax = axes('Parent',fig,'Position',pos_cm(1.25,8.35,7.6,7.7));
hold(ax,'on'); axis(ax,'equal'); axis(ax,'off');
patch(ax,[0 L L 0],[0 0 L L],main_c,'EdgeColor',edge_c,'LineWidth',0.8);
plot_contours(ax,z_cont,slope,L,cont_c,font);
text(ax,28,39.1,'Single column, 3600 m^{2}','Color','w', ...
    'HorizontalAlignment','center',font{:},'FontSize',8,'FontWeight','bold');
draw_arrow(ax,20,27,20,14,'w',1.0,2.2,2.0,false);
text(ax,22.5,13.0,sprintf('%.1f %% slope (%.1f°)',100*slope,atand(slope)), ...
    'Color','w',font{:},'FontSize',6.5);
plot_section_line(ax,x_sec,L,pr,font);
plot_dimension(ax,[0 L],-4.5,'60 m','below',font);
plot_dimension(ax,[0 L],-4.5,'60 m','left',font);
set(ax,'XLim',plan_xlim,'YLim',plan_ylim);
panel_title(ax,'a','1-D, plan view',font);

%% ================= (b) 2-D plan view =================
ax = axes('Parent',fig,'Position',pos_cm(10.0,8.35,7.6,7.7));
hold(ax,'on'); axis(ax,'equal'); axis(ax,'off');
for c = 1:topo.N_cells
    if strcmp(topo.cell_region{c},'main'); fill_c = main_c; else; fill_c = buff_c; end
    patch(ax,topo.cell_x_lo(c)+[0 1 1 0]*topo.cell_dx(c), ...
        topo.cell_y_lo(c)+[0 0 1 1]*topo.cell_dy(c),fill_c, ...
        'EdgeColor',edge_c,'LineWidth',0.8);
end
plot_contours(ax,z_cont,slope,L,cont_c,font);
tol = 1e-9*L;
for e = 1:topo.N_edges
    a = topo.edge_a(e); b = topo.edge_b(e);
    if abs(topo.cell_x_hi(a)-topo.cell_x_lo(b)) < tol || ...
            abs(topo.cell_x_hi(b)-topo.cell_x_lo(a)) < tol
        mx = min(topo.cell_x_hi(a),topo.cell_x_hi(b));
        my = 0.5*(max(topo.cell_y_lo(a),topo.cell_y_lo(b)) + ...
            min(topo.cell_y_hi(a),topo.cell_y_hi(b)));
        draw_arrow(ax,mx-1.6,my,mx+1.6,my,water_c,1.0,1.2,1.4,true);
    else
        my = min(topo.cell_y_hi(a),topo.cell_y_hi(b));
        mx = 0.5*(max(topo.cell_x_lo(a),topo.cell_x_lo(b)) + ...
            min(topo.cell_x_hi(a),topo.cell_x_hi(b)));
        draw_arrow(ax,mx,my-1.6,mx,my+1.6,water_c,1.0,1.2,1.4,true);
    end
end
for c = 1:topo.N_cells
    if strcmp(topo.cell_region{c},'main')
        text(ax,topo.cell_xc(c),topo.cell_yc(c),sprintf('%d',c),'Color','w', ...
            'HorizontalAlignment','center',font{:},'FontSize',9,'FontWeight','bold');
    else
        %% corner numbers sit toward the outer corner, clear of the arrows
        off = 0.7*(topo.cell_dx(c) == w && topo.cell_dy(c) == w)* ...
            [sign(topo.cell_xc(c)-L/2) sign(topo.cell_yc(c)-L/2)];
        text(ax,topo.cell_xc(c)+off(1),topo.cell_yc(c)+off(2),sprintf('%d',c), ...
            'Color','k','HorizontalAlignment','center',font{:},'FontSize',5.5);
    end
end
plot_section_line(ax,x_sec,L,pr,font);
plot_dimension(ax,[0 L],-4.5,'60 m','below',font);
plot_dimension(ax,[w w+(L-2*w)/nx],L+3.2,'25 m','above',font);
plot_dimension(ax,[L-w L],L+3.2,'5 m','above',font);
set(ax,'XLim',plan_xlim,'YLim',plan_ylim);
panel_title(ax,'b','2-D, plan view',font);

%% ================= (c) 1-D section A-A' =================
ax = axes('Parent',fig,'Position',pos_cm(1.25,2.8,7.6,4.5));
hold(ax,'on'); axis(ax,'equal');
draw_section_column(ax,0,L,z_surf,Zd,Zbio,VE,soil_light_c,soil_deep_c,soil_alt_c,soil_line_c);
draw_cap(ax,0,L,z_surf,VE,main_c,edge_c);
draw_sides(ax,z_surf,Zd(end),L,VE);
for sx = [10 30 50]
    draw_arrow(ax,sx,VE*(z_surf(sx)-Zd(end))-0.4,sx,VE*(z_surf(sx)-Zd(end))-3.4, ...
        drain_c,0.9,1.2,1.3,false);
end
%% inset: top 0.4 m of the mesh at true scale, called out near A
z_in = 0.4; sc = 40;                       % depth shown [m], inset units per m
box = [19 46 -19.9 0.6];                   % x0 x1 y0 y1 [axes units]
sx0 = 5; sx1 = 9;                          % called-out strip along A-A'
patch(ax,[sx0 sx1 sx1 sx0],VE*[z_surf(sx0) z_surf(sx1) z_surf(sx1)-z_in z_surf(sx0)-z_in], ...
    'w','FaceColor','none','EdgeColor','k','LineWidth',0.6);
plot(ax,[sx1 box(1)],[VE*z_surf(sx1) box(4)],':k','LineWidth',0.6);
plot(ax,[sx1 box(1)],[VE*(z_surf(sx1)-z_in) box(3)],':k','LineWidth',0.6);
patch(ax,box([1 2 2 1]),box([3 3 4 4]),'w','EdgeColor','k','LineWidth',0.6);
y0 = box(4) - 3.3;
for i = 1:find(Zd < z_in,1,'last')
    if Zd(i+1) <= Zbio; fill_c = soil_light_c; elseif mod(i,2) == 0; fill_c = soil_alt_c; else; fill_c = soil_deep_c; end
    patch(ax,[31 37 37 31],y0-sc*[Zd(i) Zd(i) min(Zd(i+1),z_in) min(Zd(i+1),z_in)], ...
        fill_c,'EdgeColor',soil_line_c,'LineWidth',0.25);
end
for d = Zd(Zd <= z_in & Zd ~= 0.01)
    text(ax,30.2,y0-sc*d,sprintf('%.2f',d),'HorizontalAlignment','right', ...
        font{:},'FontSize',5.5);
end
plot(ax,[38 38.8 38.8 38],y0-sc*[0 0 Zbio Zbio],'k','LineWidth',0.5);
text(ax,39.6,y0-sc*Zbio/2,{'Active','zone'},font{:},'FontSize',5.5);
text(ax,24.2,y0-sc*z_in/2,'Depth (m)','Rotation',90,'HorizontalAlignment','center', ...
    font{:},'FontSize',5.5);
text(ax,mean(box(1:2)),box(4)-0.6,'Upper 0.4 m of 20 layers','HorizontalAlignment','center', ...
    'VerticalAlignment','top',font{:},'FontSize',5.5);
plot_section_axes(ax,sec_xlim,sec_ylim,VE,0:20:L,sec_name,font);
panel_title(ax,'c',['1-D, section ' sec_name],font);

%% ================= (d) 2-D section A-A' =================
ax = axes('Parent',fig,'Position',pos_cm(10.0,2.8,7.6,4.5));
hold(ax,'on'); axis(ax,'equal');
for k = 1:numel(sec_cells)
    s0 = s_face(k); s1 = s_face(k+1);
    draw_section_column(ax,s0,s1,z_surf,Zd,Zbio,VE,soil_light_c,soil_deep_c,soil_alt_c,soil_line_c);
    if strcmp(topo.cell_region{sec_cells(k)},'main'); fill_c = main_c; else; fill_c = buff_c; end
    draw_cap(ax,s0,s1,z_surf,VE,fill_c,edge_c);
    text(ax,0.5*(s0+s1),VE*(z_surf(0.5*(s0+s1))+0.22)+0.5,sprintf('%d',sec_cells(k)), ...
        'HorizontalAlignment','center','VerticalAlignment','bottom',font{:},'FontWeight','bold');
end
for f = s_face(2:end-1)
    plot(ax,[f f],VE*[z_surf(f)-Zd(end) z_surf(f)],'Color',edge_c,'LineWidth',0.6);
    for d = 1:4
        draw_arrow(ax,f-2.6,VE*(z_surf(f)-d),f+2.6,VE*(z_surf(f)-d),water_c,1.0,1.2,1.3,true);
    end
    draw_arrow(ax,f-2.6,VE*(z_surf(f)-Zbio/2),f+2.6,VE*(z_surf(f)-Zbio/2),solute_c,1.0,1.2,1.3,true);
end
draw_sides(ax,z_surf,Zd(end),L,VE);
for sx = 0.5*(s_face(1:end-1)+s_face(2:end))
    draw_arrow(ax,sx,VE*(z_surf(sx)-Zd(end))-0.4,sx,VE*(z_surf(sx)-Zd(end))-3.4, ...
        drain_c,0.9,1.2,1.3,false);
end
plot_section_axes(ax,sec_xlim,sec_ylim,VE,s_face,sec_name,font);
panel_title(ax,'d',['2-D, section ' sec_name],font);

%% ================= legend =================
ax = axes('Parent',fig,'Position',pos_cm(0.4,0.1,17.2,1.5));
hold(ax,'on'); axis(ax,'equal'); axis(ax,'off');
set(ax,'XLim',[0 172],'YLim',[0 15]);
row = [12.5 7.5 2.5]; col = [2 58 114];
patch(ax,col(1)+[0 5 5 0],row(1)+[-1.4 -1.4 1.4 1.4],main_c,'EdgeColor',edge_c,'LineWidth',0.5);
text(ax,col(1)+6.5,row(1),'Main cell',font{:});
patch(ax,col(1)+[0 5 5 0],row(2)+[-1.4 -1.4 1.4 1.4],buff_c,'EdgeColor',edge_c,'LineWidth',0.5);
text(ax,col(1)+6.5,row(2),'Margin cell',font{:});
plot(ax,col(1)+[0 5],row(3)*[1 1],'--','Color',cont_c,'LineWidth',0.5);
text(ax,col(1)+6.5,row(3),'Elevation contour (0.2 m)',font{:});
patch(ax,col(2)+[0 5 5 0],row(1)+[-1.4 -1.4 1.4 1.4],soil_light_c,'EdgeColor',soil_line_c,'LineWidth',0.5);
text(ax,col(2)+6.5,row(1),['Active zone (0' nd '0.3 m)'],font{:});
patch(ax,col(2)+[0 5 5 0],row(2)+[-1.4 -1.4 1.4 1.4],soil_deep_c,'EdgeColor',soil_line_c,'LineWidth',0.5);
text(ax,col(2)+6.5,row(2),['Deeper soil (0.3' nd '5 m)'],font{:});
plot(ax,col(2)+[2.5 2.5],row(3)+[-1.6 1.6],'k','LineWidth',1.8);
text(ax,col(2)+6.5,row(3),'Closed side, no lateral flux',font{:});
draw_arrow(ax,col(3),row(1),col(3)+6,row(1),water_c,1.0,1.4,1.6,true);
text(ax,col(3)+7.5,row(1),'Water and heat, every layer',font{:});
draw_arrow(ax,col(3),row(2),col(3)+6,row(2),solute_c,1.0,1.4,1.6,true);
text(ax,col(3)+7.5,row(2),'Aqueous solutes and soil CO_{2}, active zone',font{:});
draw_arrow(ax,col(3)+3,row(3)+1.7,col(3)+3,row(3)-1.7,drain_c,0.9,1.3,1.5,false);
text(ax,col(3)+7.5,row(3),'Free drainage',font{:});

%% ================= save =================
plot_file = fullfile(figure_dir,'Grid_Setup.png');
exportgraphics(fig,plot_file,'Resolution',600);
[plot_folder,plot_name] = fileparts(plot_file);
pdf_file = fullfile(plot_folder,[plot_name '.pdf']);
exportgraphics(fig,pdf_file,'ContentType','vector');
fig_file = fullfile(plot_folder,[plot_name '.fig']);
savefig(fig,fig_file);
fprintf('Saved: %s\n',plot_file);
fprintf('Saved: %s\n',pdf_file);
fprintf('Saved: %s\n',fig_file);

%% ================= local functions =================
function draw_arrow(ax,x0,y0,x1,y1,c,lw,hl,hw,both)
%% Arrow with filled heads sized in data units (the axes are equal-aspect).
u = [x1-x0, y1-y0]/hypot(x1-x0,y1-y0);
n = [-u(2) u(1)];
p0 = [x0 y0]; p1 = [x1 y1];
q0 = p0 + both*hl*u; q1 = p1 - hl*u;
plot(ax,[q0(1) q1(1)],[q0(2) q1(2)],'Color',c,'LineWidth',lw);
tri = [p1; p1-hl*u+0.5*hw*n; p1-hl*u-0.5*hw*n];
patch(ax,tri(:,1),tri(:,2),c,'EdgeColor','none');
if both
    tri = [p0; p0+hl*u+0.5*hw*n; p0+hl*u-0.5*hw*n];
    patch(ax,tri(:,1),tri(:,2),c,'EdgeColor','none');
end
end

function draw_section_column(ax,s0,s1,z_surf,Zd,Zbio,VE,c_act,c_deep,c_alt,c_line)
%% Soil layers between s0 and s1, each at a fixed depth below the slope.
for i = 1:numel(Zd)-1
    if Zd(i+1) <= Zbio
        fill_c = c_act;
    elseif mod(i,2) == 0
        fill_c = c_alt;
    else
        fill_c = c_deep;
    end
    patch(ax,[s0 s1 s1 s0], ...
        VE*[z_surf(s0)-Zd(i) z_surf(s1)-Zd(i) z_surf(s1)-Zd(i+1) z_surf(s0)-Zd(i+1)], ...
        fill_c,'EdgeColor',c_line,'LineWidth',0.25);
end
plot(ax,[s0 s1],VE*[z_surf(s0)-Zbio z_surf(s1)-Zbio],'--','Color',c_line,'LineWidth',0.6);
end

function draw_cap(ax,s0,s1,z_surf,VE,c,edge_c)
%% Land-cover band above the soil surface, coloured by cell type.
patch(ax,[s0 s1 s1 s0],VE*[z_surf(s0) z_surf(s1) z_surf(s1)+0.22 z_surf(s0)+0.22],c, ...
    'EdgeColor',edge_c,'LineWidth',0.6);
end

function draw_sides(ax,z_surf,Zmax,L,VE)
%% Closed lateral boundary at both ends of the section.
plot(ax,[0 0],VE*[z_surf(0)-Zmax z_surf(0)+0.22],'k','LineWidth',1.8);
plot(ax,[L L],VE*[z_surf(L)-Zmax z_surf(L)+0.22],'k','LineWidth',1.8);
end

function plot_contours(ax,z_cont,slope,L,c,font)
%% Contours of the plane; elevation above the downslope edge.
for k = 1:numel(z_cont)
    yk = z_cont(k)/slope;
    plot(ax,[0 L],[yk yk],'--','Color',c,'LineWidth',0.5);
    text(ax,L+0.8,yk,sprintf('%.1f',z_cont(k)),'Color',c,font{:},'FontSize',5.5);
end
end

function plot_section_axes(ax,xl,yl,VE,xt,sec_name,font)
%% Tick labels give true elevation; the y data are scaled by VE.
set(ax,'XLim',xl,'YLim',yl,'XTick',xt,'YTick',VE*(-5:1), ...
    'YTickLabel',arrayfun(@(v) sprintf('%d',v),-5:1,'UniformOutput',false), ...
    'TickDir','out','TickLength',[0.012 0.012],'XColor','k','YColor','k', ...
    'Box','off','LineWidth',0.5,font{:});
xlabel(ax,['Distance along ' sec_name ' (m)'],font{:});
ylabel(ax,'Elevation (m)',font{:});
text(ax,xl(2),yl(2),sprintf('VE \\times%d',VE),'HorizontalAlignment','right', ...
    'VerticalAlignment','top',font{:},'FontSize',6);
end

function plot_section_line(ax,x_sec,L,pr,font)
plot(ax,[x_sec x_sec],[L+1.5 -1.5],'-.','Color','k','LineWidth',0.8);
text(ax,x_sec,L+3.4,'A','HorizontalAlignment','center',font{:},'FontWeight','bold');
text(ax,x_sec,-3.0,['A' pr],'HorizontalAlignment','center',font{:},'FontWeight','bold');
end

function plot_dimension(ax,span,pos,label,kind,font)
%% Dimension bar with end ticks, placed 'below', 'left' or 'above'.
t = 0.8;
switch kind
    case 'below'
        plot(ax,span,[pos pos],'k','LineWidth',0.5);
        plot(ax,[span(1) span(1)],pos+[-t t],'k',[span(2) span(2)],pos+[-t t],'k','LineWidth',0.5);
        text(ax,mean(span),pos-1.8,label,'HorizontalAlignment','center', ...
            'VerticalAlignment','top',font{:});
    case 'left'
        plot(ax,[pos pos],span,'k','LineWidth',0.5);
        plot(ax,pos+[-t t],[span(1) span(1)],'k',pos+[-t t],[span(2) span(2)],'k','LineWidth',0.5);
        text(ax,pos-1.8,mean(span),label,'HorizontalAlignment','center', ...
            'VerticalAlignment','bottom','Rotation',90,font{:});
    case 'above'
        plot(ax,span,[pos pos],'k','LineWidth',0.5);
        plot(ax,[span(1) span(1)],pos+[-t t],'k',[span(2) span(2)],pos+[-t t],'k','LineWidth',0.5);
        text(ax,mean(span),pos+1.4,label,'HorizontalAlignment','center', ...
            'VerticalAlignment','bottom',font{:});
end
end

function panel_title(ax,letter,txt,font)
title(ax,['\bf' letter '\rm   ' txt],font{:},'FontSize',8,'FontWeight','normal', ...
    'HorizontalAlignment','left','Units','normalized','Position',[0 1.0 0]);
end
