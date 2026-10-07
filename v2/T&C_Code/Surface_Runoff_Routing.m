function [q_deliver_mm, C_deliver, exit_water_m3, exit_mass, diag, ...
    exit_by_cell_m3, exit_mass_by_cell] = ...
    Surface_Runoff_Routing(topo, runoff_mm, C_runoff, dt_h)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Surface_Runoff_Routing
%%
%% Route this step's ponded overland runoff down the DEM topology,
%% carrying its dissolved load.  Returns (a) the depth + concentration
%% delivered to each cell as its NEXT vertical step's re-infiltration
%% opportunity (q_runon) and (b) the water + mass that leaves the domain
%% at the topographic outlet(s).
%%
%% MODEL -- runoff / run-on connectivity
%%   Ponded water is partitioned among a cell's downslope neighbours by
%%   the Quinn (1991) multiple-flow-direction weights (Attach_Topography_
%%   To_Topology) and cascaded downslope within the step, in descending
%%   elevation order.  What happens at a receiving cell depends on what
%%   the vertical solver just said about it:
%%
%%     TRANSMITTING  the cell produced runoff itself this step (runoff_mm
%%                   > 0): it is at infiltration capacity (Hortonian) or
%%                   saturated (Dunne), so run-on cannot infiltrate there
%%                   -- it passes straight through, downslope, in the
%%                   same step, together with the cell's own runoff.
%%     ABSORBING     the cell produced no runoff (spare infiltration
%%                   capacity): run-on stops here and is offered to its
%%                   next vertical step as q_runon, where the infiltration
%%                   solver decides how much re-infiltrates; any excess
%%                   becomes that cell's runoff and is routed next step.
%%     OUTLET        a TRANSMITTING outlet cell passes water and mass out
%%                   of the domain; a dry outlet cell absorbs like any
%%                   other dry cell (it is a full cell of field, not a
%%                   channel) and exports only once it is at capacity.
%%     PIT           an interior cell with no lower neighbour keeps what
%%                   reaches it (offered back as q_runon).
%%
%%   On a hillslope that is generating runoff everywhere (a storm) this
%%   gives same-step delivery to the outlet, as overland travel times
%%   (minutes) demand; a one-hop-per-step rule would instead walk the
%%   water down the slope over N hours and offer it re-infiltration at
%%   every hop after the rain has stopped -- a systematic bias toward
%%   retention.  Where runoff meets dry ground, the run-on stops and is
%%   allowed to infiltrate, the classic runoff-runon mechanism.  Solute
%%   moves with the water (upwind, fully mixed at each confluence).  The
%%   scheme is mass-exact; the guard at the end enforces it.
%%
%%   No topography attached -> no weights, no outlet: every cell keeps its
%%   own runoff (q_deliver == runoff), the pre-topography behaviour.
%%
%% INPUTS
%%   topo        topology; uses cell_z_m, edge_overland_weight_ab/ba,
%%               cell_is_outlet, cell_area_m2, edge_a/edge_b, N_edges
%%   runoff_mm   [N_cells x 1] ponded runoff depth this step [mm]
%%   C_runoff    [N_cells x N_species] its aqueous concentration
%%               [mass / m3 water].  [] => no solute tracking.
%%   dt_h        hydrologic timestep [h] (interface symmetry / future use)
%%
%% OUTPUTS
%%   q_deliver_mm      [N_cells x 1] depth delivered to each cell [mm]
%%   C_deliver         [N_cells x N_species] its concentration [mass/m3]
%%   exit_water_m3     scalar water leaving at outlets [m3]
%%   exit_mass         [N_species x 1] solute mass leaving at outlets
%%   diag              .edge_flux_m3 [N_edges x 1] overland water per edge
%%                     (a->b positive), .pit_cells, .transmitting (logical)
%%   exit_by_cell_m3   [N_cells x 1] outlet water attributed per cell
%%   exit_mass_by_cell [N_cells x N_species] outlet mass per cell
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N_cells = topo.N_cells;
runoff_mm = runoff_mm(:);
if numel(runoff_mm) ~= N_cells || any(~isfinite(runoff_mm)) || any(runoff_mm < -1e-9)
    error('Surface_Runoff_Routing:BadRunoff', ...
        'runoff_mm must hold one finite nonnegative depth per cell.');
end
runoff_mm = max(0,runoff_mm);
track_solute = nargin >= 3 && ~isempty(C_runoff);
if track_solute
    N_sp = size(C_runoff,2);
    if size(C_runoff,1) ~= N_cells || any(~isfinite(C_runoff(:))) || any(C_runoff(:) < 0)
        error('Surface_Runoff_Routing:BadConc', ...
            'C_runoff must be [N_cells x N_species], finite and nonnegative.');
    end
else
    N_sp = 0;
end
if nargin < 4; dt_h = 1; end %#ok<NASGU>

A = topo.cell_area_m2(:);
has_topo = isfield(topo,'edge_overland_weight_ab') && ...
    isfield(topo,'cell_z_m') && ~isempty(topo.cell_z_m);

diag = struct('edge_flux_m3',zeros(topo.N_edges,1),'pit_cells',[], ...
    'transmitting',false(N_cells,1));
exit_water_m3      = 0;
exit_by_cell_m3    = zeros(N_cells,1);
exit_mass_by_cell  = zeros(N_cells,max(N_sp,1));

if ~has_topo || all(runoff_mm == 0)
    q_deliver_mm = runoff_mm;
    if track_solute
        C_deliver = C_runoff; exit_mass = zeros(N_sp,1);
        exit_mass_by_cell = exit_mass_by_cell(:,1:N_sp);
    else
        C_deliver = zeros(N_cells,0); exit_mass = zeros(0,1);
        exit_mass_by_cell = zeros(N_cells,0);
    end
    return
end

is_outlet = false(N_cells,1);
if isfield(topo,'cell_is_outlet'); is_outlet = logical(topo.cell_is_outlet(:)); end
w_ab = topo.edge_overland_weight_ab(:);
if isfield(topo,'edge_overland_weight_ba')
    w_ba = topo.edge_overland_weight_ba(:);
else
    w_ba = zeros(topo.N_edges,1);
end
z = topo.cell_z_m(:);

pond_m3 = (runoff_mm/1000).*A;                       % this step's own runoff
mass    = zeros(N_cells,max(N_sp,1));
if track_solute; mass(:,1:N_sp) = pond_m3 .* C_runoff; end
%% A cell that produced runoff is at capacity this step and transmits
%% run-on; a tolerance keeps round-off residue from flagging a dry cell.
transmitting = runoff_mm > 1e-9;

inflow_m3   = zeros(N_cells,1);                      % run-on arriving this step
inflow_mass = zeros(N_cells,max(N_sp,1));
delivered_m3   = zeros(N_cells,1);
delivered_mass = zeros(N_cells,max(N_sp,1));
pit = false(N_cells,1);

%% Process every contributor before the cell it drains to.  The order is
%% a topological sort of the drainage digraph the WEIGHTS define, not of
%% elevation: with the nodal drop the two coincide, but user-supplied
%% weights or outlet overrides could make them differ, and water pushed
%% onto an already-swept cell would be stranded (the mass-balance guard
%% below would then abort the run).  Descending elevation is the
%% fallback if the weight graph is somehow cyclic.
src = [topo.edge_a(w_ab > 0); topo.edge_b(w_ba > 0)];
dst = [topo.edge_b(w_ab > 0); topo.edge_a(w_ba > 0)];
order = [];
if ~isempty(src)
    G = digraph(src(:),dst(:),ones(numel(src),1),N_cells);
    if isdag(G)
        order = toposort(G);
    else
        warning('Surface_Runoff_Routing:CyclicWeights', ...
            ['The overland weight graph contains a cycle; sweeping in ' ...
             'descending elevation instead.']);
    end
end
if isempty(order); [~,order] = sort(z,'descend'); end
order = order(:).';
for c = order
    Vc = pond_m3(c) + inflow_m3(c);
    Mc = mass(c,:) + inflow_mass(c,:);
    if Vc <= 0; continue; end
    if ~transmitting(c)
        %% Spare infiltration capacity: run-on stops here (own runoff is
        %% zero by definition of ~transmitting, up to the tolerance).
        %% This holds at an outlet cell too -- it is a full cell of field,
        %% not a channel, and only exports once it is itself at capacity.
        delivered_m3(c)     = delivered_m3(c) + Vc;
        delivered_mass(c,:) = delivered_mass(c,:) + Mc;
        continue
    end
    if is_outlet(c)
        exit_water_m3      = exit_water_m3 + Vc;
        exit_by_cell_m3(c) = exit_by_cell_m3(c) + Vc;
        exit_mass_by_cell(c,:) = exit_mass_by_cell(c,:) + Mc;
        continue
    end
    ea_out = find(topo.edge_a == c & w_ab > 0);
    eb_out = find(topo.edge_b == c & w_ba > 0);
    targets = [topo.edge_b(ea_out); topo.edge_a(eb_out)];
    weights = [w_ab(ea_out);        w_ba(eb_out)];
    edges   = [ea_out;              eb_out];
    esign   = [ones(numel(ea_out),1); -ones(numel(eb_out),1)];   % a->b positive
    wsum = sum(weights);
    if wsum <= 0
        pit(c) = true;
        delivered_m3(c)     = delivered_m3(c) + Vc;
        delivered_mass(c,:) = delivered_mass(c,:) + Mc;
        continue
    end
    weights = weights/wsum;
    Cc = Mc/Vc;                                      % mixed concentration
    for t = 1:numel(targets)
        nt = targets(t);
        dV = Vc*weights(t);
        diag.edge_flux_m3(edges(t)) = diag.edge_flux_m3(edges(t)) + esign(t)*dV;
        inflow_m3(nt)     = inflow_m3(nt) + dV;
        inflow_mass(nt,:) = inflow_mass(nt,:) + dV*Cc;
    end
end

q_deliver_mm = 1000*delivered_m3./A;
if track_solute
    C_deliver = zeros(N_cells,N_sp);
    nz = delivered_m3 > 0;
    C_deliver(nz,:) = delivered_mass(nz,1:N_sp)./delivered_m3(nz);
    exit_mass         = sum(exit_mass_by_cell(:,1:N_sp),1).';
    exit_mass_by_cell = exit_mass_by_cell(:,1:N_sp);
else
    C_deliver = zeros(N_cells,0);
    exit_mass = zeros(0,1);
    exit_mass_by_cell = zeros(N_cells,0);
end
diag.pit_cells = find(pit);
diag.transmitting = transmitting;

%% Mass-balance guard.
in_water  = sum((runoff_mm/1000).*A);
out_water = sum(delivered_m3) + exit_water_m3;
if abs(in_water-out_water) > 1e-8*max(1,in_water)
    error('Surface_Runoff_Routing:WaterImbalance', ...
        'Routed water (%.6g m3) != input (%.6g m3).',out_water,in_water);
end
if track_solute
    in_mass  = sum(pond_m3.*C_runoff,1).';
    out_mass = sum(delivered_mass(:,1:N_sp),1).' + exit_mass;
    bad = abs(in_mass-out_mass) > 1e-7*max(1,max(abs(in_mass)));
    if any(bad)
        error('Surface_Runoff_Routing:MassImbalance', ...
            'Routed solute species %d does not conserve.',find(bad,1));
    end
end
end
