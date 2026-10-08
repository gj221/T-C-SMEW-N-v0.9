function h = Draw_Field_Cells(ax,G,value,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Draw_Field_Cells
%%
%% Draw one value per cell on ax through the axes colormap: the cell faces
%% of Field_Cell_Geometry, cells without a value in NaNColor, the part of
%% perimeter squares outside the field in pale grey, the faces between
%% main and margin cells dashed, and the field outline.
%%
%%   h = Draw_Field_Cells(ax,G,value)
%%
%% Options
%%   'EdgeColor'   cell edges; 'auto' draws them up to 400 cells
%%   'ShowOutside' the simulated part of perimeter squares outside the field
%%   'ShowZones'   the main/margin faces (dashed)
%%   'ShowOutline' the field outline
%%   'NaNColor'    cells with no value
%%   'Axes'        'off' (no ticks: map panels) or 'on' (metres)
%%   'ScaleBar'    draw a scale bar in the lower left
%%   'NorthArrow'  draw the north arrow in the upper right
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
p.addParameter('EdgeColor','auto');
p.addParameter('ShowOutside',true);
p.addParameter('ShowZones',true);
p.addParameter('ShowOutline',true);
p.addParameter('NaNColor',[0.86 0.86 0.86]);
p.addParameter('Axes','off');
p.addParameter('ScaleBar',false);
p.addParameter('NorthArrow',false);
p.parse(varargin{:});
o = p.Results;

value = double(value(:));
if numel(value) ~= G.N_cells
    error('Draw_Field_Cells:BadValues','Expected %d values, one per cell.',G.N_cells);
end
edge = o.EdgeColor;
edge_width = 0.3;
if ischar(edge) && strcmp(edge,'auto')
    if G.N_cells <= 400
        edge = [0.25 0.25 0.25];
    elseif G.N_cells <= 3000
        edge = [0.60 0.60 0.60];         %% thin and light: the squares stay visible
        edge_width = 0.15;
    else
        edge = 'none';
    end
end
hold(ax,'on');
h = struct('outside',[],'cells',[],'nan',[],'zones',[],'outline',[]);
if o.ShowOutside && G.clipped
    h.outside = patch(ax,'Faces',G.square_faces,'Vertices',G.square_vertices, ...
        'FaceColor',[0.96 0.96 0.96],'EdgeColor',[0.78 0.78 0.78],'LineWidth',0.3, ...
        'Tag','OutsideField');
end
face_value = value(G.face_cell);
ok = isfinite(face_value);
if any(ok)
    h.cells = patch(ax,'Faces',G.faces(ok,:),'Vertices',G.vertices, ...
        'FaceVertexCData',face_value(ok),'FaceColor','flat', ...
        'EdgeColor',edge,'LineWidth',edge_width,'Tag','CellValues');
end
if any(~ok)
    h.nan = patch(ax,'Faces',G.faces(~ok,:),'Vertices',G.vertices, ...
        'FaceColor',o.NaNColor,'EdgeColor',edge,'LineWidth',edge_width,'Tag','CellNoValue');
end
if o.ShowZones && ~isempty(G.zone_x)
    h.zones = plot(ax,G.zone_x,G.zone_y,'--','Color',[0.05 0.05 0.05], ...
        'LineWidth',0.9,'Tag','MarginBoundary');
end
if o.ShowOutline
    h.outline = plot(ax,G.outline_x,G.outline_y,'-','Color','k', ...
        'LineWidth',1.4,'Tag','FieldOutline');
end
axis(ax,'equal');
pad = 0.03*max(diff(G.extent(1:2)),diff(G.extent(3:4)));
xlim(ax,G.extent(1:2)+[-pad pad]);
ylim(ax,G.extent(3:4)+[-pad pad]);
ax.Color = 'w';
if strcmpi(o.Axes,'off')
    ax.XTick = []; ax.YTick = [];
    ax.XColor = 'none'; ax.YColor = 'none';
else
    ax.XColor = 'k'; ax.YColor = 'k';
    ax.Box = 'on';
    xlabel(ax,G.xlabel,'Color','k');
    ylabel(ax,G.ylabel,'Color','k');
end
if o.ScaleBar
    draw_scale_bar(ax,G,pad);
end
if o.NorthArrow && ~isempty(G.north)
    draw_north_arrow(ax,G,pad);
end
end

%% =======================================================================
function draw_scale_bar(ax,G,pad)
width = diff(G.extent(1:2));
nice = [1 2 5 10 20 25 50 100 200 250 500 1000];
len = nice(find(nice <= 0.3*width,1,'last'));
if isempty(len); len = nice(1); end
x0 = G.extent(1); y0 = G.extent(3)-0.6*pad;
plot(ax,[x0 x0+len],[y0 y0],'-','Color','k','LineWidth',2,'Tag','ScaleBar');
text(ax,x0+len/2,y0,sprintf('%g m',len),'HorizontalAlignment','center', ...
    'VerticalAlignment','top','FontSize',8,'Color','k');
end

function draw_north_arrow(ax,G,pad)
L = 0.08*max(diff(G.extent(1:2)),diff(G.extent(3:4)));
x0 = G.extent(2)-0.2*pad; y0 = G.extent(4)-L-0.2*pad;
u = G.north(:).'/norm(G.north);
quiver(ax,x0,y0,u(1)*L,u(2)*L,0,'Color','k','LineWidth',1.2, ...
    'MaxHeadSize',0.6,'Tag','NorthArrow');
text(ax,x0+u(1)*L*1.25,y0+u(2)*L*1.25,'N','FontSize',8,'FontWeight','bold', ...
    'HorizontalAlignment','center','Color','k');
end
