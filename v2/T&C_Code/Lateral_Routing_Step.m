function [theta_cell, dM_cell, diag] = Lateral_Routing_Step( ...
    topo, theta_cell, C_aq_cell, soilpar_cell, Zlayer_m, dt_h)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lateral_Routing_Step
%%
%% One hydrologic-timestep lateral routing call on the sub-plot cell
%% topology.  Operator-split from the vertical solver: call once per
%% timestep, after every column has taken its vertical step.
%%
%%   1. Darcy_Lateral_Topology     lateral water exchange (closed perimeter)
%%   2. Solute_Lateral_Transport   upwind advection of the aqueous phase
%%
%% Inputs
%%   topo          topology struct from Build_Cell_Topology
%%   theta_cell    [N_cells x ms] volumetric water content per cell-layer
%%   C_aq_cell     [N_cells x ms x N_species] aqueous concentration
%%                 [mass/m3 water].  The caller owns the pool mapping and
%%                 the re-aggregation back into model state.
%%   soilpar_cell  [N_cells x 1] struct array (see Soil_Param_To_Routing);
%%                 per-cell entries may differ, so cells can carry
%%                 different soil hydraulic properties.
%%   Zlayer_m      [1 x ms] soil layer thickness [m]
%%   dt_h          hydrologic timestep [h]
%%
%% Outputs
%%   theta_cell    updated [N_cells x ms]
%%   dM_cell       [N_cells x ms x N_species] solute mass change [mass]
%%   diag          .sub_steps  adaptive sub-steps taken by the Darcy solver
%%                 .Q_edge     [N_edges x ms] per-edge water flux [m3/h]
%%                 .buf_to_main.water  net water buffer->main this step [m3]
%%                 .buf_to_main.net    [N_species x 1] net mass buffer->main
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N_cells = topo.N_cells;
ms      = numel(Zlayer_m);
N_sp    = size(C_aq_cell, 3);

assert(size(theta_cell,1) == N_cells, 'theta_cell row count must equal topo.N_cells');
assert(size(theta_cell,2) == ms, 'theta_cell column count must equal numel(Zlayer_m)');
assert(size(C_aq_cell,1) == N_cells, 'C_aq_cell row count must equal topo.N_cells');
assert(size(C_aq_cell,2) == ms, 'C_aq_cell layer count must equal numel(Zlayer_m)');
assert(numel(soilpar_cell) == N_cells, 'soilpar_cell must have topo.N_cells entries');

%% --- 1. Lateral Buckingham-Darcy water flux ---------------------------
theta_before = theta_cell;
[Q_edge, dtheta_lat, sub_steps] = Darcy_Lateral_Topology( ...
    topo, theta_cell, Zlayer_m, soilpar_cell, dt_h);
theta_cell = theta_cell + dtheta_lat;

%% --- 2. Upwind advective solute redistribution ------------------------
Vw_layer = theta_before .* (topo.cell_area_m2(:)*Zlayer_m(:).');
M_aq_cell = C_aq_cell .* reshape(Vw_layer,N_cells,ms,1);
[dM_cell, dM_edge] = Solute_Lateral_Transport( ...
    topo,Q_edge,C_aq_cell,M_aq_cell,dt_h);

%% --- 3. Buffer <-> main diagnostic ------------------------------------
%% Sign each buffer-main edge so that positive means buffer -> main.
is_main_a = strcmp(topo.cell_region(topo.edge_a), 'main');
is_main_b = strcmp(topo.cell_region(topo.edge_b), 'main');
bm_sign = zeros(topo.N_edges, 1);
bm_sign(~is_main_a &  is_main_b) = +1;
bm_sign( is_main_a & ~is_main_b) = -1;

diag.sub_steps = sub_steps;
diag.Q_edge    = Q_edge;
diag.buf_to_main.water = sum(bm_sign .* sum(Q_edge, 2)) * dt_h;     %% m3
net_sp = zeros(N_sp, 1);
for e = find(bm_sign ~= 0).'
    net_sp = net_sp + bm_sign(e) * reshape(sum(dM_edge(e,:,:), 2), N_sp, 1);
end
diag.buf_to_main.net = net_sp;

end
