function grid = Field_Topographic_Attributes(dem)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Field_Topographic_Attributes
%%
%% Derive the terrain attributes used by the topographic lateral/vertical
%% flow code from a gridded elevation model, on the DEM grid itself.
%% Attach_Topography_To_Topology samples these at the cell centres;
%% Run_Field_DEM plots them.
%%
%% INPUT
%%   dem   struct with
%%           .Z  [ny x nx]  elevation [m] (NaN outside the mapped area)
%%         and one of
%%           .x [1 x nx], .y [1 x ny]   projected coordinates of grid
%%                                      column / row centres [m], or
%%           .R                         a map raster reference object
%%                                      (readgeoraster).
%%
%% OUTPUT (all [ny x nx] unless noted, on an ascending-x, ascending-y grid)
%%   .Z            elevation [m]
%%   .x .y         resolved ascending coordinate vectors [m]
%%   .cellsize     nominal grid spacing 0.5*(dx+dy) [m]
%%   .slope        steepest slope [-] (fraction, = tan of the slope angle)
%%   .aspect_deg   aspect [deg from north, clockwise]
%%   .dzdx .dzdy   surface gradient components [-] (x east, y north)
%%   .curv_prof    profile (down-slope) curvature [1/m], convex-up positive
%%   .curv_plan    plan (contour) curvature [1/m], divergent positive
%%   .twi          topographic wetness index ln(a / tan beta) [-]
%%   .contrib_area specific catchment area a [m] (MFD, Quinn 1991) =
%%                 upslope_area_m2 / grid spacing
%%   .upslope_area_m2  total upslope (contributing) area [m2]
%%   .svf          sky-view factor [-] in [0,1]
%%   .flow_i .flow_j   D-infinity downslope unit vector components [-]
%%
%% References
%%   Zevenbergen & Thorne (1987); Quinn et al. (1991); Tarboton (1997);
%%   Beven & Kirkby (1979).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isstruct(dem) || ~isscalar(dem) || ~isfield(dem,'Z')
    error('Field_Topographic_Attributes:BadInput', ...
        'dem must be a scalar struct with a .Z elevation matrix.');
end
Z = double(dem.Z);
if ~(ismatrix(Z) && ~isvector(Z))
    error('Field_Topographic_Attributes:BadInput', ...
        'dem.Z must be a 2-D elevation matrix.');
end
[ny,nx] = size(Z);

if isfield(dem,'x') && isfield(dem,'y')
    x = double(dem.x(:)).'; y = double(dem.y(:)).';
elseif isfield(dem,'R')
    R = dem.R;
    x = R.XWorldLimits(1)+((1:nx)-0.5)*R.CellExtentInWorldX;
    y = R.YWorldLimits(1)+((1:ny)-0.5)*R.CellExtentInWorldY;
    if ~isprop(R,'ColumnsStartFrom') || strcmp(R.ColumnsStartFrom,'north')
        Z = flipud(Z);   % Z rows run north->south; match ascending y
    end
else
    error('Field_Topographic_Attributes:BadInput', ...
        'dem needs .x/.y coordinate vectors or a raster reference .R.');
end
if numel(x) ~= nx || numel(y) ~= ny
    error('Field_Topographic_Attributes:BadInput', ...
        'dem.x/.y lengths must match size(dem.Z).');
end
if numel(y) > 1 && y(2) < y(1); y = fliplr(y); Z = flipud(Z); end
if numel(x) > 1 && x(2) < x(1); x = fliplr(x); Z = fliplr(Z); end
dx = median(diff(x)); dy = median(diff(y));
if ~(isfinite(dx) && dx > 0 && isfinite(dy) && dy > 0)
    error('Field_Topographic_Attributes:BadInput', ...
        'dem.x/.y must be monotone with a positive spacing.');
end
cellsize = 0.5*(dx+dy);

%% Fill small interior NaN gaps so the derivatives are defined; keep a
%% record of where elevation was genuinely missing.
missing = ~isfinite(Z);
if any(missing(:)) && ~all(missing(:))
    Z = fill_small_gaps(Z);
end

%% Slope [-] and aspect [deg from N].  Slope_Aspect_indexes('grad') returns
%% SLO = tan(steepest slope angle) and ASP in degrees from north.
[slope,aspect_deg] = Slope_Aspect_indexes(Z,cellsize,'grad');
slope = max(0,slope);

%% Curvature (Zevenbergen & Thorne 1987), convex-up / divergent positive.
[gx,gy]   = gradient(Z,dx,dy);
[gxx,gxy] = gradient(gx,dx,dy);
[~,gyy]   = gradient(gy,dx,dy);
gmag2 = gx.^2 + gy.^2;
curv_prof = -(gxx.*gx.^2 + 2*gxy.*gx.*gy + gyy.*gy.^2) ./ ...
    (max(gmag2,eps).*(1+gmag2).^1.5);
curv_plan = -(gxx.*gy.^2 - 2*gxy.*gx.*gy + gyy.*gx.^2) ./ ...
    max(gmag2,eps).^1.5;
flat = gmag2 < 1e-10;
curv_prof(flat) = 0; curv_plan(flat) = 0;

%% D-infinity downslope direction (Tarboton 1997) for display + overland.
try
    [Rdir,~] = dem_flow(Z,dx,dy);
    flow_i = -sin(Rdir);   % grid row increases downward (north up); flip y
    flow_j =  cos(Rdir);
    flow_i(~isfinite(flow_i)) = 0;
    flow_j(~isfinite(flow_j)) = 0;
catch ME
    warning('Field_Topographic_Attributes:NoFlowDirection', ...
        'dem_flow failed (%s); D-infinity display vectors set to zero.', ...
        ME.message);
    flow_i = zeros(ny,nx); flow_j = zeros(ny,nx);
end

%% Contributing area (MFD, Quinn 1991) and TWI.
[upslope_area_m2,specific_area] = mfd_contributing_area(Z,dx,dy);
beta = atan(max(slope,1e-4));
twi = log(specific_area./tan(beta));

%% Sky-view factor (terrain openness).
try
    [HZ,Zasp] = Horizon_Angle(Z,cellsize);
    svf = Sky_View_Factor(Z,atan(slope)*180/pi,aspect_deg,HZ,Zasp);
    %% Guard BEFORE clipping: max/min drop NaN silently, which would turn
    %% an undefined sky view into total obstruction (0).
    bad = ~isfinite(svf) & isfinite(Z);
    if any(bad(:))
        warning('Field_Topographic_Attributes:SkyViewNaN', ...
            '%d sky-view values undefined; set to open sky (1).',nnz(bad));
        svf(bad) = 1;
    end
    svf = min(1,max(0,svf));
    svf(~isfinite(Z)) = NaN;   % genuinely missing nodes stay undefined
catch ME
    %% Never silently: an unobstructed sky (svf = 1) is a physics change
    %% (no diffuse shading, no terrain-reflected term), so say so.
    warning('Field_Topographic_Attributes:NoSkyView', ...
        ['Sky-view factor could not be computed (%s); using svf = 1 ' ...
         '(open sky) everywhere.'], ME.message);
    svf = ones(ny,nx);
end

Z(missing) = NaN;   % restore genuine gaps in the returned elevation

grid = struct('Z',Z,'x',x,'y',y,'cellsize',cellsize, ...
    'slope',slope,'aspect_deg',aspect_deg,'dzdx',gx,'dzdy',gy, ...
    'curv_prof',curv_prof,'curv_plan',curv_plan, ...
    'twi',twi,'contrib_area',specific_area, ...
    'upslope_area_m2',upslope_area_m2,'svf',svf, ...
    'flow_i',flow_i,'flow_j',flow_j,'missing',missing);
end

%% =======================================================================
function [a,specific_area] = mfd_contributing_area(Z,dx,dy)
[ny,nx] = size(Z);
cellarea = dx*dy;
a = cellarea*ones(ny,nx);
[~,order] = sort(Z(:),'descend');   % NaNs sort first and are skipped
[ii,jj] = ndgrid(1:ny,1:nx);
di = [-1 -1 -1 0 0 1 1 1];
dj = [-1 0 1 -1 1 -1 0 1];
ln = sqrt((di*dy).^2+(dj*dx).^2);
p = 1.1;
for node = order.'
    if ~isfinite(Z(node)); continue; end
    i = ii(node); j = jj(node);
    drop = zeros(1,8); valid = false(1,8);
    for m = 1:8
        ni = i+di(m); nj = j+dj(m);
        if ni>=1 && ni<=ny && nj>=1 && nj<=nx && isfinite(Z(ni,nj))
            d = (Z(i,j)-Z(ni,nj))/ln(m);
            if d > 0; drop(m) = d; valid(m) = true; end
        end
    end
    if ~any(valid); continue; end
    wsum = sum(drop(valid).^p);
    for m = find(valid)
        a(i+di(m),j+dj(m)) = a(i+di(m),j+dj(m)) + a(i,j)*drop(m).^p/wsum;
    end
end
specific_area = a/(0.5*(dx+dy));
end

%% =======================================================================
function Zf = fill_small_gaps(Z)
%% Nearest-neighbour fill of NaN cells (iterative dilation).  Adequate for
%% the sliver gaps a clipped lidar DEM leaves; large voids should be
%% pre-filled before this function is called.
Zf = Z;
bad = ~isfinite(Zf);
guard = 0;
while any(bad(:)) && guard < 500
    guard = guard + 1;
    [ii,jj] = find(bad);
    for k = 1:numel(ii)
        i = ii(k); j = jj(k);
        w = Zf(max(1,i-1):min(end,i+1),max(1,j-1):min(end,j+1));
        w = w(isfinite(w));
        if ~isempty(w); Zf(i,j) = mean(w); end
    end
    bad = ~isfinite(Zf);
end
end
