function topo = Build_Single_Cell_Topology(area_m2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Build_Single_Cell_Topology
%%
%% Degenerate one-cell topology: a single column covering the whole plot,
%% with no edges and therefore no lateral exchange.
%%
%% Its purpose is verification.  Running MAIN_FRAME_BG_ERW_2D on this
%% topology must reproduce the single-column MAIN_FRAME_BG_ERW_10mins
%% exactly, because the lateral block is skipped when N_cells == 1.  Any
%% divergence is a transcription error in the multi-column driver rather
%% than a difference in physics, which is what makes this the sharpest
%% available test of the port.
%%
%% Input
%%   area_m2  plot area [m2] (default 3600, the 60 m square)
%%
%% Output
%%   topo     struct with the same field names Build_Cell_Topology emits
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 1 || isempty(area_m2); area_m2 = 3600; end
if ~isnumeric(area_m2) || ~isscalar(area_m2) || ...
        ~isfinite(area_m2) || area_m2 <= 0
    error('Build_Single_Cell_Topology:BadArea', ...
        'Single-cell plot area must be finite and positive.');
end
L = sqrt(area_m2);

topo.field_side_m         = L;
topo.buffer_width_m       = 0;
topo.N_cells              = 1;
topo.cell_id              = 1;
topo.cell_name            = {'Cell_Single'};
topo.cell_region          = {'main'};
topo.cell_subtype         = {'main'};
topo.cell_x_lo = 0;  topo.cell_x_hi = L;
topo.cell_y_lo = 0;  topo.cell_y_hi = L;
topo.cell_xc   = L/2; topo.cell_yc  = L/2;
topo.cell_dx   = L;  topo.cell_dy   = L;
topo.cell_area_m2         = area_m2;
topo.N_edges              = 0;
topo.edge_a               = zeros(0,1);
topo.edge_b               = zeros(0,1);
topo.edge_face_width_m    = zeros(0,1);
topo.edge_face_distance_m = zeros(0,1);
topo.edge_face_distance_a_m = zeros(0,1);
topo.edge_face_distance_b_m = zeros(0,1);
topo.edge_active          = zeros(0,1);
topo.adj                  = sparse(false(1,1));
topo.outer_face_width_m   = 4*L;

end
