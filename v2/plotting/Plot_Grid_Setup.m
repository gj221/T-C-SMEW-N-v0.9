%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot_Grid_Setup and schematic of the T&C-SMEW-N v0.9 1-D and 2-D grids
%%
%% Standalone script; run it directly (e.g.  run Plot_Grid_Setup).
%%
%% Publication figure, 18 cm wide, Arial 7 pt:
%%
%%   (a)  1-D block:  one 186 m x 186 m column (3.46 ha), 5 m deep, under a
%%                    planar slope; section line A-A'
%%   (b)  2-D block:  the same plot as 16 cells (4 main cells of 73 m, the
%%                    20 m ring split at the grid lines into 8 strip pieces
%%                    and 4 corners) with its perimeter open: a boundary
%%                    column of untreated land (Open_Perimeter_Topology's
%%                    ghost) beyond each of the 16 exposed sides; arrows
%%                    mark the lateral exchange across the 24 shared faces
%%                    and the 16 boundary faces
%%   (c)  1-D section A-A':  one column, drainage at the base
%%   (d)  2-D section A-A':  the cells crossed by A-A' at their true
%%                    widths and the ghosts at both ends; water and heat
%%                    cross each face in every layer, aqueous solutes and
%%                    soil CO2 in the active zone
%%   (e)  upper 0.4 m of the 20-layer mesh and the 0-300 mm active zone
%%
%% The plot is plan 4a-C of 6 November 2026: a square of the basalt
%% field's area (3.46 ha, outline about 185 m x 186 m) with a 20 m margin,
%% run flat and on the plane fitted to the field's DEM, 2.17 % downhill
%% toward 223 degrees, drawn here.  Elevations are above the lowest
%% (south-west) corner.  The slope cannot be seen at true scale, so (a) to
%% (d) are stretched vertically by VE through the axes' data aspect
%% ratio; every elevation and tick stays in true metres.
%%
%% Cell and ghost geometry is taken from Build_Cell_Topology(L, w, nx, ny,
%% true) and Open_Perimeter_Topology (runtime/open_perimeter_20261110), so
%% the script runs on a tree carrying that patch (L = 186 m, w = 20 m,
%% 2 x 2 main grid of 73 m cells).  The
%% layer mesh matches ERW_MOD_PARAM_parameters.xlsx (Kantola rows,
%% Zbio = 300 mm).
%%
%% Outputs: Grid_Setup.png (600 dpi) and Grid_Setup.fig in v2/figures.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% This script lives in v2/plotting; the repository root is one level up.
here = fileparts(mfilename('fullpath'));
if isempty(here); here = pwd; end
root = fileparts(here);
addpath(fullfile(root,'T&C_Code'));
figure_dir = fullfile(root,'figures');
if ~isfolder(figure_dir); mkdir(figure_dir); end

%% geometry
L  = 186;                % field side [m], the basalt field's 3.46 ha as a square
w  = 20;                 % margin width [m], one lime-spreader pass
nx = 2;  ny = 2;         % main grid
slope = 0.0217;          % planar slope [m/m], plane fitted to the basalt field's DEM
aspect = 223;            % downhill direction [degrees from north]
VE = 10;                 % vertical exaggeration of (a) to (d) [-]
x_sec = 146;             % plan position of section A-A' [m]
sc = L/60;               % arrow and label scale against the 60 m plot they were sized on
Zs_mm = [0 10 50 100 150 200 300 400 600 800 1000 1250 1500 1750 ...
    2000 2500 3000 3500 4000 4500 5000];
Zd = Zs_mm/1000;         % layer boundaries [m]
Zbio = 0.3;              % active biogeochemical depth [m]
topo = Build_Cell_Topology(L,w,nx,ny,true);
if ~isfield(topo,'split_strips')
    error('Plot_Grid_Setup:Tree', ...
        'This tree''s Build_Cell_Topology does not split the strips; apply the open-perimeter patch.');
end
topo_x = Open_Perimeter_Topology(topo);
ghost = topo.N_cells + (1:topo_x.N_ghost);
host = topo_x.ghost_host(:).';
g_off = [topo_x.cell_xc(ghost)-topo.cell_xc(host) topo_x.cell_yc(ghost)-topo.cell_yc(host)];  % ghost minus host centre [m]
gx = -slope*sind(aspect); gy = -slope*cosd(aspect);   % elevation gradient, as Open_Perimeter_Topology
z_plan = @(x,y) gx*x + gy*y;          % surface elevation above the south-west corner [m]
z_surf = @(s) z_plan(x_sec,L - s);    % the same along A-A', s = 0 at A (north)
sec_cells = find(topo.cell_x_lo(:) <= x_sec & topo.cell_x_hi(:) > x_sec);
[~,o] = sort(topo.cell_yc(sec_cells),'descend');
sec_cells = sec_cells(o).';          % cells crossed by A-A', upslope first
s_face = L - [topo.cell_y_hi(sec_cells(:)).' topo.cell_y_lo(sec_cells(end))];  % [m]
pr = '′'; nd = '–';                  % prime and en dash
sec_name = ['A' nd 'A' pr];

%% colours
main_c  = [0.21 0.45 0.70];   % blue: main cells
buff_c  = [0.80 0.89 0.97];   % light blue: margin cells
edge_c  = [0.15 0.15 0.15];
soil_c  = [0.88 0.78 0.60; 0.72 0.58 0.40; 0.77 0.64 0.46];  % active zone, deeper soil, alternate shade
soil_line_c = [0.45 0.33 0.20];   % brown layer lines
water_c  = [0.84 0.27 0.00];      % vermillion: water and heat, every layer
solute_c = [0.00 0.55 0.40];      % bluish green: solutes and CO2, active zone
drain_c  = [0.30 0.30 0.30];
ghost_c  = [0.80 0.80 0.80];      % grey: ghost columns
ghost_soil_c = [0.88 0.88 0.88; 0.76 0.76 0.76; 0.82 0.82 0.82];
font = {'FontName','Arial','FontSize',7};

%% figure (cm); axes are placed by hand so the fixed-aspect panels fill them
sec_xlim = [-w-1.5*sc L+w+1.5*sc];
z_end = z_surf([-w L+w]);
sec_ylim = [floor(min(z_end)-Zd(end)-1.6) ceil(max(z_end)+1.4)];
sec_h = 6.3*VE*diff(sec_ylim)/diff(sec_xlim);   % section axes height [cm] at equal x scale
y_top = 2.1 + sec_h + 0.55;                     % bottom of the block panels [cm]
fig_w = 18; fig_h = y_top + 5.7;
fig = figure('Color','w','Name','T&C-SMEW-N 1-D and 2-D grid setup', ...
    'Units','centimeters','Position',[2 2 fig_w fig_h],'PaperPositionMode','auto');
pos_cm = @(x,y,wd,ht) [x/fig_w y/fig_h wd/fig_w ht/fig_h];

%% ================= (a) 1-D block =================
ax = axes('Parent',fig,'Position',pos_cm(0.1,y_top,8.9,5.6));
hold(ax,'on');
draw_block(ax,L,z_plan,Zd,Zbio,soil_c,soil_line_c,edge_c,[]);
patch(ax,[0 L L 0],[0 0 L L],z_plan([0 L L 0],[0 0 L L]),main_c, ...
    'EdgeColor',edge_c,'LineWidth',0.6);
draw_section_mark(ax,x_sec,L,z_plan,Zd(end),pr,font);
text(ax,-2*sc,0,z_plan(0,0)+0.3,'0 m','HorizontalAlignment','right',font{:});
text(ax,L+3*sc,0,z_plan(L,0)+0.3,sprintf('%.2f m',z_plan(L,0)),font{:});
text(ax,L+2*sc,L,z_plan(L,L),sprintf('%.2f m',z_plan(L,L)),font{:});
text(ax,L+2*sc,0,z_plan(L,0)-Zd(end),sprintf('%.2f m',z_plan(L,0)-Zd(end)),font{:});
text(ax,L/2,-4*sc,z_plan(L/2,0)-Zd(end),sprintf('%d m',L),'HorizontalAlignment','center','VerticalAlignment','top',font{:});
text(ax,L+4*sc,L/2,z_plan(L,L/2)-Zd(end),sprintf('%d m',L),'VerticalAlignment','top',font{:});
setup_block_axes(ax,VE,L,w,z_plan,Zd(end));
panel_letter(fig,pos_cm(0.2,y_top+5.3,0.5,0.4),'a',font);

%% ================= (b) 2-D block =================
ax = axes('Parent',fig,'Position',pos_cm(9.0,y_top,8.9,5.6));
hold(ax,'on');
draw_block(ax,L,z_plan,Zd,Zbio,soil_c,soil_line_c,edge_c,w+(0:ny)*(L-2*w)/ny);
for c = 1:topo.N_cells
    if strcmp(topo.cell_region{c},'main'); fill_c = main_c; else; fill_c = buff_c; end
    xs = topo.cell_x_lo(c)+[0 1 1 0]*topo.cell_dx(c);
    ys = topo.cell_y_lo(c)+[0 0 1 1]*topo.cell_dy(c);
    patch(ax,xs,ys,z_plan(xs,ys),fill_c,'EdgeColor',edge_c,'LineWidth',0.6);
end
for g = 1:numel(ghost)
    xs = topo.cell_x_lo(host(g))+[0 1 1 0]*topo.cell_dx(host(g))+g_off(g,1);
    ys = topo.cell_y_lo(host(g))+[0 0 1 1]*topo.cell_dy(host(g))+g_off(g,2);
    patch(ax,xs,ys,z_plan(xs,ys),ghost_c,'FaceAlpha',0.35,'EdgeColor',edge_c,'LineStyle','--','LineWidth',0.5);
    if g_off(g,2) < 0 || g_off(g,1) > 0   % outer faces seen from the south-east
        p0 = [xs(1) ys(1)] + (g_off(g,1) > 0)*[topo.cell_dx(host(g)) 0];
        p1 = p0 + (g_off(g,2) < 0)*[topo.cell_dx(host(g)) 0] + (g_off(g,1) > 0)*[0 topo.cell_dy(host(g))];
        patch(ax,[p0(1) p1(1) p1(1) p0(1)],[p0(2) p1(2) p1(2) p0(2)], ...
            [z_plan(p0(1),p0(2)) z_plan(p1(1),p1(2)) z_plan(p1(1),p1(2))-Zd(end) z_plan(p0(1),p0(2))-Zd(end)], ...
            ghost_soil_c(2,:),'FaceAlpha',0.35,'EdgeColor',edge_c,'LineStyle','--','LineWidth',0.5);
    end
end
z_top = @(x,y) z_plan(x,y) + 0.02;   % arrows sit just above the surface
for g = 1:numel(ghost)
    m = [topo.cell_xc(host(g)) topo.cell_yc(host(g))] + g_off(g,:)/2;   % face centre
    u = 2.4*sc*g_off(g,:)/norm(g_off(g,:));
    draw_arrow(ax,m(1)-u(1),m(2)-u(2),m(1)+u(1),m(2)+u(2),water_c,1.0,1.8*sc,2.0*sc,true,z_top,1);
end
tol = 1e-9*L;
for e = 1:topo.N_edges
    a = topo.edge_a(e); b = topo.edge_b(e);
    if abs(topo.cell_x_hi(a)-topo.cell_x_lo(b)) < tol || ...
            abs(topo.cell_x_hi(b)-topo.cell_x_lo(a)) < tol
        mx = min(topo.cell_x_hi(a),topo.cell_x_hi(b));
        my = 0.5*(max(topo.cell_y_lo(a),topo.cell_y_lo(b)) + ...
            min(topo.cell_y_hi(a),topo.cell_y_hi(b)));
        draw_arrow(ax,mx-2.4*sc,my,mx+2.4*sc,my,water_c,1.0,1.8*sc,2.0*sc,true,z_top,1);
    else
        my = min(topo.cell_y_hi(a),topo.cell_y_hi(b));
        mx = 0.5*(max(topo.cell_x_lo(a),topo.cell_x_lo(b)) + ...
            min(topo.cell_x_hi(a),topo.cell_x_hi(b)));
        draw_arrow(ax,mx,my-2.4*sc,mx,my+2.4*sc,water_c,1.0,1.8*sc,2.0*sc,true,z_top,1);
    end
end
for c = find(strcmp(topo.cell_region,'main')).'
    text(ax,topo.cell_xc(c),topo.cell_yc(c),z_plan(topo.cell_xc(c),topo.cell_yc(c))+0.4,sprintf('%d',c), ...
        'Color','w','HorizontalAlignment','center',font{:},'FontSize',8,'FontWeight','bold');
end
draw_section_mark(ax,x_sec,L,z_plan,Zd(end),pr,font);
setup_block_axes(ax,VE,L,w,z_plan,Zd(end));
panel_letter(fig,pos_cm(9.0,y_top+5.3,0.5,0.4),'b',font);

%% ================= (c) 1-D section A-A' =================
ax = axes('Parent',fig,'Position',pos_cm(1.15,2.1,6.3,sec_h));
hold(ax,'on');
draw_layers(ax,0,L,z_surf(0),z_surf(L),Zd,Zbio,soil_c,soil_line_c);
draw_cap(ax,0,L,z_surf,main_c,edge_c);
for sx = [10 30 50]*sc
    draw_arrow(ax,sx,z_surf(sx)-Zd(end)-0.15,sx,z_surf(sx)-Zd(end)-1.15,drain_c,0.8,1.0*sc,1.2*sc,false,[],VE);
end
setup_section_axes(ax,sec_xlim,sec_ylim,[0 L/2 L],VE,sec_name,font);
panel_letter(fig,pos_cm(0.1,2.1+sec_h+0.4,0.5,0.4),'c',font);

%% ================= (d) 2-D section A-A' =================
ax = axes('Parent',fig,'Position',pos_cm(8.45,2.1,6.3,sec_h));
hold(ax,'on');
for s0 = [-w L]                      % ghosts beyond A and A'
    draw_layers(ax,s0,s0+w,z_surf(s0),z_surf(s0+w),Zd,Zbio,ghost_soil_c,ghost_c*0.7);
    draw_cap(ax,s0,s0+w,z_surf,ghost_c,edge_c);
    plot(ax,[s0 s0+w s0+w s0 s0],[z_surf(s0)+0.3 z_surf(s0+w)+0.3 z_surf(s0+w)-Zd(end) z_surf(s0)-Zd(end) z_surf(s0)+0.3], ...
        '--','Color',edge_c,'LineWidth',0.5);
end
for k = 1:numel(sec_cells)
    s0 = s_face(k); s1 = s_face(k+1);
    draw_layers(ax,s0,s1,z_surf(s0),z_surf(s1),Zd,Zbio,soil_c,soil_line_c);
    if strcmp(topo.cell_region{sec_cells(k)},'main'); fill_c = main_c; else; fill_c = buff_c; end
    draw_cap(ax,s0,s1,z_surf,fill_c,edge_c);
    if strcmp(topo.cell_region{sec_cells(k)},'main')
        text(ax,0.5*(s0+s1),z_surf(0.5*(s0+s1))+0.45,sprintf('%d',sec_cells(k)), ...
            'HorizontalAlignment','center','VerticalAlignment','bottom',font{:},'FontWeight','bold');
    end
end
for f = s_face
    plot(ax,[f f],[z_surf(f)-Zd(end) z_surf(f)],'Color',edge_c,'LineWidth',0.6);
    for d = 1:4
        draw_arrow(ax,f-2.6*sc,z_surf(f)-d,f+2.6*sc,z_surf(f)-d,water_c,1.0,1.1*sc,1.2*sc,true,[],VE);
    end
    draw_arrow(ax,f-2.6*sc,z_surf(f)-Zbio/2,f+2.6*sc,z_surf(f)-Zbio/2,solute_c,1.0,1.1*sc,1.2*sc,true,[],VE);
end
for sx = 0.5*(s_face(1:end-1)+s_face(2:end))
    draw_arrow(ax,sx,z_surf(sx)-Zd(end)-0.15,sx,z_surf(sx)-Zd(end)-1.15,drain_c,0.8,1.0*sc,1.2*sc,false,[],VE);
end
setup_section_axes(ax,sec_xlim,sec_ylim,s_face,VE,sec_name,font);
panel_letter(fig,pos_cm(7.6,2.1+sec_h+0.4,0.5,0.4),'d',font);

%% ================= (e) upper 0.4 m of the mesh =================
z_in = 0.4;              % depth shown [m]
ax = axes('Parent',fig,'Position',pos_cm(15.75,2.1,0.7,sec_h));
hold(ax,'on');
for i = 1:find(Zd < z_in,1,'last')
    patch(ax,[0 1 1 0],-[Zd(i) Zd(i) min(Zd(i+1),z_in) min(Zd(i+1),z_in)], ...
        layer_colour(i,Zd,Zbio,soil_c),'EdgeColor',soil_line_c,'LineWidth',0.4);
end
plot(ax,[-0.15 -0.35 -0.35 -0.15],-[0 0 Zbio Zbio],'k','LineWidth',0.5);
text(ax,-0.5,-Zbio/2,'Active zone','Rotation',90,'HorizontalAlignment','center', ...
    'VerticalAlignment','bottom',font{:});
d_lab = Zd(Zd <= z_in & Zd ~= 0.01);   % 10 mm sits too close to 0 to label
set(ax,'XLim',[0 1],'YLim',[-z_in 0],'XTick',[],'YTick',fliplr(-d_lab), ...
    'YTickLabel',arrayfun(@(v) sprintf('%.2f',v),fliplr(d_lab),'UniformOutput',false), ...
    'YAxisLocation','right','TickDir','out','TickLength',[0.04 0.04], ...
    'XColor','none','YColor','k','Box','off','LineWidth',0.5,'Clipping','off',font{:});
ylabel(ax,'Depth (m)',font{:});
panel_letter(fig,pos_cm(14.95,2.1+sec_h+0.4,0.5,0.4),'e',font);

%% ================= legend =================
ax = axes('Parent',fig,'Position',pos_cm(0.4,0.05,17.2,1.0));
hold(ax,'on'); axis(ax,'equal'); axis(ax,'off');
set(ax,'XLim',[0 172],'YLim',[0 10]);
row = [7.2 2.5]; col = [2 30 66 124];
legend_patch(ax,col(1),row(1),main_c,edge_c,'Main cell',font);
legend_patch(ax,col(1),row(2),buff_c,edge_c,'Margin cell',font);
legend_patch(ax,col(2),row(1),soil_c(1,:),soil_line_c,['Active zone (0' nd '0.3 m)'],font);
legend_patch(ax,col(2),row(2),soil_c(2,:),soil_line_c,['Deeper soil (0.3' nd '5 m)'],font);
draw_arrow(ax,col(3),row(1),col(3)+6,row(1),water_c,1.0,1.4,1.6,true,[],1);
text(ax,col(3)+7.5,row(1),'Water and heat, every layer',font{:});
draw_arrow(ax,col(3),row(2),col(3)+6,row(2),solute_c,1.0,1.4,1.6,true,[],1);
text(ax,col(3)+7.5,row(2),'Solutes and soil CO_{2}, active zone',font{:});
patch(ax,col(4)+[0 5 5 0],row(1)+[-1.4 -1.4 1.4 1.4],ghost_c,'EdgeColor',edge_c,'LineStyle','--','LineWidth',0.5);
text(ax,col(4)+6.5,row(1),'Boundary column',font{:});
draw_arrow(ax,col(4)+2.5,row(2)+1.7,col(4)+2.5,row(2)-1.7,drain_c,0.9,1.3,1.5,false,[],1);
text(ax,col(4)+6.5,row(2),'Drainage',font{:});

%% ================= save =================
plot_file = fullfile(figure_dir,'Grid_Setup.png');
exportgraphics(fig,plot_file,'Resolution',600);
[plot_folder,plot_name] = fileparts(plot_file);
fig_file = fullfile(plot_folder,[plot_name '.fig']);
savefig(fig,fig_file);
fprintf('Saved: %s\n',plot_file);
fprintf('Saved: %s\n',fig_file);

%% ================= local functions =================
function draw_arrow(ax,x0,y0,x1,y1,c,lw,hl,hw,both,zfun,ar)
%% Arrow with filled heads; hl and hw are in x units.  ar is the axes'
%% y exaggeration, so the heads are built in display space and keep their
%% shape.  With zfun the arrow is lifted onto the surface z = zfun(x,y).
p0 = [x0 y0*ar]; p1 = [x1 y1*ar];
u = (p1-p0)/norm(p1-p0);
n = [-u(2) u(1)];
q0 = p0 + both*hl*u; q1 = p1 - hl*u;
tri = {[p1; p1-hl*u+0.5*hw*n; p1-hl*u-0.5*hw*n]};
if both; tri{2} = [p0; p0+hl*u+0.5*hw*n; p0+hl*u-0.5*hw*n]; end
line_x = [q0(1) q1(1)]; line_y = [q0(2) q1(2)]/ar;
if isempty(zfun)
    plot(ax,line_x,line_y,'Color',c,'LineWidth',lw);
    for k = 1:numel(tri)
        patch(ax,tri{k}(:,1),tri{k}(:,2)/ar,c,'EdgeColor','none');
    end
else
    plot3(ax,line_x,line_y,zfun(line_x,line_y),'Color',c,'LineWidth',lw);
    for k = 1:numel(tri)
        patch(ax,tri{k}(:,1),tri{k}(:,2)/ar,zfun(tri{k}(:,1),tri{k}(:,2)/ar),c,'EdgeColor','none');
    end
end
end

function c = layer_colour(i,Zd,Zbio,soil_c)
%% Active zone, then the deeper layers in two alternating shades.
if Zd(i+1) <= Zbio
    c = soil_c(1,:);
else
    c = soil_c(2+mod(i,2),:);
end
end

function draw_layers(ax,p0,p1,zt0,zt1,Zd,Zbio,soil_c,c_line)
%% Soil layers on a vertical face from p0 to p1 (s along a section, or
%% [x y] on a block face), each at a fixed depth below the surface whose
%% elevation is zt0 at p0 and zt1 at p1.
for i = 1:numel(Zd)-1
    zs = [zt0-Zd(i) zt1-Zd(i) zt1-Zd(i+1) zt0-Zd(i+1)];
    if isscalar(p0)
        patch(ax,[p0 p1 p1 p0],zs,layer_colour(i,Zd,Zbio,soil_c), ...
            'EdgeColor',c_line,'LineWidth',0.2);
    else
        patch(ax,[p0(1) p1(1) p1(1) p0(1)],[p0(2) p1(2) p1(2) p0(2)],zs, ...
            layer_colour(i,Zd,Zbio,soil_c),'EdgeColor',c_line,'LineWidth',0.2);
    end
end
end

function draw_block(ax,L,z_plan,Zd,Zbio,soil_c,c_line,edge_c,cuts)
%% Soil block under the plane: the downslope (y = 0) and east (x = L)
%% faces carry the layer mesh; cuts are the cell boundaries drawn on both.
draw_layers(ax,[0 0],[L 0],z_plan(0,0),z_plan(L,0),Zd,Zbio,soil_c,c_line);
draw_layers(ax,[L 0],[L L],z_plan(L,0),z_plan(L,L),Zd,Zbio,soil_c,c_line);
for q = cuts
    plot3(ax,[q q],[0 0],z_plan(q,0)-[0 Zd(end)],'Color',edge_c,'LineWidth',0.5);
    plot3(ax,[L L],[q q],z_plan(L,q)-[0 Zd(end)],'Color',edge_c,'LineWidth',0.5);
end
plot3(ax,[0 L L L],[0 0 0 L],[z_plan(0,0) z_plan(L,0) z_plan(L,0) z_plan(L,L)]-Zd(end), ...
    'Color',edge_c,'LineWidth',0.6);
plot3(ax,[0 0],[0 0],z_plan(0,0)-[0 Zd(end)],'Color',edge_c,'LineWidth',0.6);
plot3(ax,[L L],[L L],z_plan(L,L)-[0 Zd(end)],'Color',edge_c,'LineWidth',0.6);
plot3(ax,[L L],[0 0],z_plan(L,0)-[0 Zd(end)],'Color',edge_c,'LineWidth',0.6);
end

function draw_section_mark(ax,x_sec,L,z_plan,Zmax,pr,font)
%% Section line on the surface and down the downslope face.
plot3(ax,[x_sec x_sec],[L 0],z_plan([x_sec x_sec],[L 0])+0.03,'-.k','LineWidth',0.8);
plot3(ax,[x_sec x_sec],[0 0],z_plan(x_sec,0)-[0 Zmax],':k','LineWidth',0.8);
text(ax,x_sec,L+L/5,z_plan(x_sec,L),'A','HorizontalAlignment','center',font{:},'FontWeight','bold');
text(ax,x_sec,-L/5,z_plan(x_sec,0),['A' pr],'HorizontalAlignment','center',font{:},'FontWeight','bold');
end

function draw_cap(ax,s0,s1,z_surf,c,edge_c)
%% Land-cover band above the soil surface, coloured by cell type.
patch(ax,[s0 s1 s1 s0],[z_surf(s0) z_surf(s1) z_surf(s1)+0.3 z_surf(s0)+0.3],c, ...
    'EdgeColor',edge_c,'LineWidth',0.5);
end

function setup_block_axes(ax,VE,L,w,z_plan,Zmax)
%% Oblique view from the south-east corner, z stretched by VE; the limits
%% hold the plot and its boundary columns.
axis(ax,'off');
zc = z_plan([0 L 0 L],[0 0 L L]);
set(ax,'DataAspectRatio',[1 1 1/VE],'XLim',[-w L+w]+[-0.05 0.15]*L,'YLim',[-w L+w]+[-0.1 0.1]*L, ...
    'ZLim',[min(zc)-Zmax-1 max(zc)+1],'Clipping','off');
view(ax,62,18);
end

function setup_section_axes(ax,xl,yl,xt,VE,sec_name,font)
set(ax,'DataAspectRatio',[1 1/VE 1],'XLim',xl,'YLim',yl,'XTick',xt,'YTick',2*ceil(yl(1)/2):2:yl(2), ...
    'TickDir','out','TickLength',[0.012 0.012],'XColor','k','YColor','k', ...
    'Box','off','LineWidth',0.5,font{:});
xlabel(ax,['Distance along ' sec_name ' (m)'],font{:});
ylabel(ax,'Elevation (m)',font{:});
end

function legend_patch(ax,x,y,c,edge_c,label,font)
patch(ax,x+[0 5 5 0],y+[-1.4 -1.4 1.4 1.4],c,'EdgeColor',edge_c,'LineWidth',0.5);
text(ax,x+6.5,y,label,font{:});
end

function panel_letter(fig,pos,letter,font)
annotation(fig,'textbox',pos,'String',letter,'LineStyle','none', ...
    'Margin',1,'VerticalAlignment','top',font{:},'FontSize',9,'FontWeight','bold');
end
