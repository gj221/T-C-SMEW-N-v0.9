function [dM_cell, dM_edge] = Solute_Lateral_Transport( ...
    topo, Q_edge, C_aq_cell, dt_h)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solute_Lateral_Transport
%%
%% First-order upwind advection on the cell topology, driven by the
%% per-edge per-layer volumetric water fluxes Q_edge produced by
%% Darcy_Lateral_Topology.  Transport is bulk advection only; the donor
%% cell is selected by the sign of Q_edge (Patankar 1980 upwind scheme).
%%
%% Hydrodynamic dispersion and molecular diffusion are deliberately NOT
%% included.  At sub-plot cell separations (tens of metres) with
%% unsaturated lateral fluxes this small, advection dominates by orders
%% of magnitude, and a dispersivity would be a free parameter with
%% nothing in the experiment to constrain it.
%%
%% The domain perimeter is closed (zero-flux Neumann), matching
%% Darcy_Lateral_Topology, so no advective mass leaves the domain
%% laterally.  Vertical solute export is handled by the biogeochemistry's
%% own leaching terms (Biogeo_Leakage / LEAK_*), not here.
%%
%% Coupling note
%%   This routine moves only the AQUEOUS mass.  The caller must derive
%%   aqueous concentration from totals before the call and re-run the ERW
%%   cation-exchange equilibrium after it, so the lateral step stays
%%   consistent with the upstream Gaines-Thomas convention without
%%   duplicating any equilibrium code.
%%
%% Inputs
%%   topo        topology struct (Build_Cell_Topology / lumped)
%%   Q_edge      [N_edges x ms] per-edge volumetric water flux [m3/h],
%%               positive cell_a -> cell_b
%%   C_aq_cell   [N_cells x ms x N_species] aqueous concentration
%%               [mass / m3 water]; units arbitrary but consistent
%%   dt_h        timestep [h]
%%
%% Outputs
%%   dM_cell     [N_cells x ms x N_species] change in solute mass [mass]
%%   dM_edge     [N_edges x ms x N_species] signed per-edge mass flux over
%%               dt_h, positive from topo.edge_a to topo.edge_b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N_cells = topo.N_cells;
N_edges = topo.N_edges;
[~, ms, N_species] = size(C_aq_cell);

assert(size(Q_edge,1) == N_edges, 'Q_edge row count must equal N_edges');
assert(size(Q_edge,2) == ms,      'Q_edge col count must equal ms');
assert(size(C_aq_cell,1) == N_cells, 'C_aq_cell row count must equal N_cells');

dM_cell = zeros(N_cells, ms, N_species);
dM_edge = zeros(N_edges, ms, N_species);

for e = 1:N_edges
    a = topo.edge_a(e);
    b = topo.edge_b(e);
    for k = 1:ms
        Q = Q_edge(e,k);                              % m3/h
        if Q > 0
            donor = a;
        elseif Q < 0
            donor = b;
        else
            continue;
        end
        for s = 1:N_species
            dM = Q * C_aq_cell(donor,k,s) * dt_h;     % mass over dt_h
            dM_cell(a,k,s) = dM_cell(a,k,s) - dM;
            dM_cell(b,k,s) = dM_cell(b,k,s) + dM;
            dM_edge(e,k,s) = dM;                      % positive a -> b
        end
    end
end

end
