function test_total_head_darcy()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% test_total_head_darcy
%%
%% Physics checks on the elevation-head (gravity) term added to
%% Darcy_Lateral_Topology and on Attach_Topography_To_Topology:
%%
%%   1  a flat surface (cell_z_m = const) reproduces the no-topography
%%      Darcy result bit for bit
%%   2  with equal matric head, water flows from the HIGH cell to the LOW
%%      cell (pure gravity drainage)
%%   3  the closed perimeter still conserves water exactly with gravity on
%%   4  a large elevation drop overrides an opposing matric-head gradient
%%   5  a steeper face carries more flux (monotone in slope)
%%   6  the slope (down-dip) face distance raises the flow resistance
%%   7  uniform recharge on a tilted plane builds a downslope-increasing
%%      storage profile and saturates the toe first
%%   8  theta stays within [Ohy, Osat] under a strong tilted-plane gradient
%%   9  Attach_Topography_To_Topology defaults (no elevation) leave the
%%      topology flat and the solver unchanged
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n=== test_total_head_darcy ===\n');
here = fileparts(mfilename('fullpath'));
code = fileparts(here);
addpath(code,fullfile(code,'Inputs'),fullfile(code,'Parameters'),here);
n_pass = 0; n_fail = 0;

    function check(name,cond,detail)
        if nargin < 3; detail = ''; end
        if cond
            n_pass = n_pass + 1; fprintf('  PASS  %s\n',name);
        else
            n_fail = n_fail + 1; fprintf('  FAIL  %s  %s\n',name,detail);
        end
    end

%% --- soil (Saxton-Rawls, SPAR = 2) ----------------------------------
ms = 4;
Zlayer_m = [0.05 0.10 0.15 0.20];
[Osat,Lslope,Pe,Ks,O33] = Soil_parameters(0.30,0.30,0.02);
Ohy_s = 0.02;
sp = struct('SPAR',2,'Ks',Ks*ones(1,ms),'Osat',Osat*ones(1,ms), ...
    'Ohy',Ohy_s*ones(1,ms),'L',Lslope*ones(1,ms),'Pe',Pe*ones(1,ms), ...
    'O33',O33*ones(1,ms),'alpVG',1e-3*ones(1,ms),'nVG',2*ones(1,ms), ...
    'lVG',0.5*ones(1,ms),'Ks_mac',zeros(1,ms),'Omac',zeros(1,ms), ...
    'alpVGM',1e-3*ones(1,ms),'nVGM',2*ones(1,ms),'lVGM',0.5*ones(1,ms), ...
    'Phy1',10000,'s_SVG',zeros(1,ms),'bVG',zeros(1,ms));

%% ---------------------------------------------------------------------
%% 1. Flat surface == no topography, bit for bit
%% ---------------------------------------------------------------------
topo = Build_Cell_Topology(60,5,2,2);
soil = repmat(sp,topo.N_cells,1);
theta = 0.20 + 0.05*mod((1:topo.N_cells)',3);
theta = repmat(theta,1,ms);

[Q0,dth0] = Darcy_Lateral_Topology(topo,theta,Zlayer_m,soil,1);
topo_flat = Attach_Topography_To_Topology(topo, ...
    struct('cell_z_m',37*ones(topo.N_cells,1)));       % constant elevation
[Q1,dth1] = Darcy_Lateral_Topology(topo_flat,theta,Zlayer_m,soil,1);
check('flat surface reproduces the matric-head-only flux exactly', ...
    isequal(Q0,Q1) && isequal(dth0,dth1), ...
    sprintf('max|dQ| = %.3e',max(abs(Q0(:)-Q1(:)))));
check('Attach_Topography_To_Topology marks the topology attached', ...
    isfield(topo_flat,'topography_attached') && topo_flat.topography_attached);
check('constant elevation gives zero edge drop', ...
    max(abs(topo_flat.edge_dz_m)) < 1e-12);

%% ---------------------------------------------------------------------
%% 2-3. Pure gravity: equal matric head, cell 1 above cell 2
%% ---------------------------------------------------------------------
row = row_topology(2,1);
row = Attach_Topography_To_Topology(row,struct('cell_z_m',[1.0;0.0]));
soil2 = repmat(sp,2,1);
th_eq = repmat(0.30,2,ms);
[Qg,dthg] = Darcy_Lateral_Topology(row,th_eq,Zlayer_m,soil2,1);
check('equal matric head + elevation drop drives flow high -> low', ...
    all(Qg(:) > 0), sprintf('Q = %s',mat2str(Qg,3)));
dVg = sum(dthg.*repmat(Zlayer_m,2,1),2).*row.cell_area_m2(:);
check('gravity flux still conserves water exactly', ...
    abs(sum(dVg)) < 1e-9*sum(abs(dVg))+1e-15, ...
    sprintf('net dV = %.3e',sum(dVg)));
check('the high cell loses and the low cell gains', ...
    dVg(1) < 0 && dVg(2) > 0);

%% ---------------------------------------------------------------------
%% 4. A big elevation drop overrides an opposing matric-head gradient.
%%    Cell 1 is HIGHER but DRIER; without gravity water would climb to it.
%% ---------------------------------------------------------------------
row4 = row_topology(2,1);
row4 = Attach_Topography_To_Topology(row4,struct('cell_z_m',[8.0;0.0]));
th4 = [repmat(0.30,1,ms); repmat(0.35,1,ms)];   % cell 1 drier(high), cell 2 wetter(low)
Q_noZ = Darcy_Lateral_Topology(row_topology(2,1),th4,Zlayer_m,soil2,1);
Q_Z   = Darcy_Lateral_Topology(row4,th4,Zlayer_m,soil2,1);
check('without gravity, water flows uphill toward the drier high cell', ...
    all(Q_noZ(:) < 0));
check('an 8 m drop reverses that: net flux is downslope', ...
    sum(Q_Z(:)) > 0, sprintf('Q_Z = %s',mat2str(Q_Z,3)));

%% ---------------------------------------------------------------------
%% 5. Steeper face -> more gravity flux (equal soil and moisture)
%% ---------------------------------------------------------------------
flux_vs_slope = zeros(1,4);
slopes = [0.01 0.05 0.15 0.35];
for s = 1:numel(slopes)
    r = row_topology(2,10);
    r = Attach_Topography_To_Topology(r,struct('cell_z_m',[10*slopes(s);0]));
    Qs = Darcy_Lateral_Topology(r,repmat(0.28,2,ms),Zlayer_m,soil2,1);
    flux_vs_slope(s) = sum(Qs(:));
end
check('gravity flux increases monotonically with face slope', ...
    all(diff(flux_vs_slope) > 0), mat2str(flux_vs_slope,3));

%% ---------------------------------------------------------------------
%% 6. The slope (down-dip) face distance raises the resistance, so a
%%    tilted face carries less flux than the same drop over a planar
%%    distance would.
%% ---------------------------------------------------------------------
r6 = row_topology(2,2);                            % 2 m centre spacing
r6 = Attach_Topography_To_Topology(r6,struct('cell_z_m',[3.0;0.0])); % steep
check('slope face distance exceeds the planar face distance', ...
    r6.edge_face_distance_surface_m > r6.edge_face_distance_m + 1e-9);
% same head difference, artificially planar distance:
r6_planar = r6;
r6_planar = rmfield(r6_planar,{'edge_face_distance_surface_a_m', ...
    'edge_face_distance_surface_b_m','edge_face_distance_surface_m'});
Q6_surf   = Darcy_Lateral_Topology(r6,repmat(0.30,2,ms),Zlayer_m,soil2,1);
Q6_planar = Darcy_Lateral_Topology(r6_planar,repmat(0.30,2,ms),Zlayer_m,soil2,1);
check('using the slope distance reduces the flux vs the planar distance', ...
    sum(Q6_surf(:)) < sum(Q6_planar(:)));

%% ---------------------------------------------------------------------
%% 7. Uniform recharge on a tilted plane with a fixed-head (draining) toe
%%    -> downslope-increasing storage in the interior.  The toe cell is
%%    held at a fixed drained moisture each step (Dirichlet boundary),
%%    emulating the outlet that WP5 will route explicitly.
%% ---------------------------------------------------------------------
n = 12; h = 5;
plane = row_topology(n,h);
z = ((n:-1:1)' - 1)*h*0.06;                        % 6% slope toward cell n
plane = Attach_Topography_To_Topology(plane,struct('cell_z_m',z));
soil_n = repmat(sp,n,1);
theta_p = repmat(0.24,n,ms);
theta_toe = repmat(0.24,1,ms);                     % held fixed each step
recharge_frac = 2e-4;                              % per step, per layer
for step = 1:1500
    theta_p(1:n-1,:) = min(theta_p(1:n-1,:) + recharge_frac, Osat-1e-6);
    theta_p(n,:) = theta_toe;
    [~,dth] = Darcy_Lateral_Topology(plane,theta_p,Zlayer_m,soil_n,1);
    theta_p = theta_p + dth;
    theta_p(n,:) = theta_toe;                      % re-impose the boundary
end
storage = sum(theta_p(1:n-1,:).*repmat(Zlayer_m,n-1,1),2);
check('interior storage increases monotonically downslope', ...
    all(diff(storage) > -1e-6), mat2str(storage.',3));
check('the foot of the slope holds more water than the crest', ...
    storage(end) > storage(1) + 1e-6, ...
    sprintf('foot %.5f vs crest %.5f m',storage(end),storage(1)));

%% ---------------------------------------------------------------------
%% 8. theta stays physical under a strong tilted-plane gradient
%% ---------------------------------------------------------------------
theta_x = repmat(linspace(Osat-1e-3,Ohy_s+1e-3,n)',1,ms);
[~,dth_x] = Darcy_Lateral_Topology(plane,theta_x,Zlayer_m,soil_n,6);
theta_end = theta_x + dth_x;
check('theta stays within [Ohy, Osat] over a 6 h strong gradient', ...
    all(theta_end(:) >= Ohy_s-1e-9) && all(theta_end(:) <= Osat+1e-9), ...
    sprintf('range [%.5f %.5f]',min(theta_end(:)),max(theta_end(:))));

%% ---------------------------------------------------------------------
%% 9. No-elevation attach leaves the topology flat
%% ---------------------------------------------------------------------
r9 = row_topology(3,1);
r9 = Attach_Topography_To_Topology(r9,struct('cell_z_m',zeros(3,1)));
check('zero elevation -> zero slope, zero edge drop', ...
    max(r9.cell_slope) < 1e-12 && max(abs(r9.edge_dz_m)) < 1e-12);
Qz = Darcy_Lateral_Topology(r9,repmat(0.30,3,ms),Zlayer_m,repmat(sp,3,1),1);
check('flat attach: uniform moisture -> zero flux', max(abs(Qz(:))) < 1e-12);

%% ---------------------------------------------------------------------
fprintf('  ---------------------------------------\n');
fprintf('  %d passed, %d failed\n\n',n_pass,n_fail);
if n_fail > 0
    error('test_total_head_darcy: %d check(s) failed',n_fail);
end
end

%% =======================================================================
function topo = row_topology(n,h)
%% A 1-D chain of n square cells of side h, edges between neighbours.
topo.N_cells = n;
topo.N_edges = n-1;
topo.edge_a = (1:n-1)';
topo.edge_b = (2:n)';
topo.edge_face_width_m      = repmat(h,n-1,1);
topo.edge_face_distance_m   = repmat(h,n-1,1);
topo.edge_face_distance_a_m = repmat(h/2,n-1,1);
topo.edge_face_distance_b_m = repmat(h/2,n-1,1);
topo.edge_active = ones(n-1,1);
topo.cell_area_m2 = repmat(h^2,n,1);
topo.cell_xc = ((1:n)'-0.5)*h;
topo.cell_yc = repmat(h/2,n,1);
topo.outer_face_width_m = repmat(2*h,n,1);
topo.outer_face_width_m([1 n]) = 3*h;
topo.adj = sparse([topo.edge_a;topo.edge_b],[topo.edge_b;topo.edge_a],1,n,n) > 0;
topo.cell_region = repmat({'main'},n,1);
end
