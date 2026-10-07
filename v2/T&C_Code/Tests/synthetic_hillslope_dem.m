function dem = synthetic_hillslope_dem(kind,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% synthetic_hillslope_dem   (TEST FIXTURE -- not for production use)
%%
%% Build an idealised elevation grid for verifying the topographic lateral
%% and vertical flow code without any downloaded DEM.  The output plugs
%% straight into Attach_Topography_To_Topology(topo, dem).
%%
%%   kind = 'plane'   uniform-slope plane, steepest descent toward -y
%%          'hollow'  concave (convergent) hillslope, valley along x = 0
%%          'nose'    convex (divergent) spur, ridge along x = 0
%%
%% Name/value options
%%   'Extent'   [Lx Ly]  domain size [m]            (default [60 60])
%%   'CellSize' scalar   grid spacing [m]           (default 1)
%%   'Slope'    scalar   mean down-y gradient [-]   (default 0.05)
%%   'Relief'   scalar   cross-slope relief [m] for hollow/nose (default 3)
%%   'Base'     scalar   elevation at (x=0,y=Ly) [m](default 100)
%%
%% Returns a struct with .Z [ny x nx], .x [1 x nx], .y [1 x ny] (ascending),
%% plus .kind and the resolved options for provenance.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 1 || isempty(kind); kind = 'plane'; end
opt = struct('Extent',[60 60],'CellSize',1,'Slope',0.05,'Relief',3,'Base',100);
for k = 1:2:numel(varargin)
    opt.(varargin{k}) = varargin{k+1};
end
Lx = opt.Extent(1); Ly = opt.Extent(2); h = opt.CellSize;

x = 0:h:Lx;
y = 0:h:Ly;
[X,Y] = meshgrid(x,y);

% Base plane: elevation decreases with y (steepest descent toward -y).
Z = opt.Base + opt.Slope*(Ly - Y);

switch lower(kind)
    case 'plane'
        % nothing further
    case 'hollow'
        % add a concave cross-section: low along the centre line x = Lx/2
        xc = (X - Lx/2)/(Lx/2);
        Z = Z + opt.Relief*(xc.^2 - mean(xc(:).^2));
    case 'nose'
        xc = (X - Lx/2)/(Lx/2);
        Z = Z - opt.Relief*(xc.^2 - mean(xc(:).^2));
    otherwise
        error('synthetic_hillslope_dem:BadKind', ...
            'kind must be plane, hollow or nose.');
end

dem = struct('Z',Z,'x',x,'y',y,'kind',lower(kind),'options',opt);
end
