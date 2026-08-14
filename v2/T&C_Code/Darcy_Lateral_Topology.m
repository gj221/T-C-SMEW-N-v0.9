function [Q_edge, dtheta_lat, sub_steps] = Darcy_Lateral_Topology( ...
    topo, theta_cell, Zlayer_m, soilpar_cell, dt_h, Oice_cell)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Darcy_Lateral_Topology
%
% Unsaturated Buckingham-Darcy lateral flow on the cell adjacency graph
% built by Build_Cell_Topology (optionally reduced by
% Lump_Topology_By_Class).  Operator-split: call once per hydrologic
% timestep AFTER the per-column vertical solver.
%
% BOUNDARY CONDITION
%   The domain perimeter is a zero-flux Neumann boundary: no water
%   crosses it laterally.  Water still leaves each column vertically
%   through the existing bottom-leakage physics (Leakage_Bottom), which
%   the caller handles.  The closed perimeter is a treatment-comparison
%   assumption: there is no lateral exchange with soil outside the 60 m
%   domain.
%
% CONVENTION NOTE (mechanistic correction over the source port)
%   Conductivity_Suction returns a convention-dependent second output:
%   SPAR=2 (Saxton-Rawls) returns POSITIVE suction [mm], while
%   SPAR=1/3/4 (van Genuchten variants) return NEGATIVE matric head
%   [mm].  The flux law below is written on matric HEAD, so Po is
%   converted per SPAR.  Assuming a single convention reverses the flow
%   direction for half the soil parameterisations.
%
% References
%   Buckingham (1907); Richards (1931); Bear (1972) ch. 9;
%   Freeze & Cherry (1979); Patankar (1980).
%
% Inputs
%   topo          topology struct (Build_Cell_Topology / lumped)
%   theta_cell    [N_cells x ms]  volumetric water content [m3/m3]
%   Zlayer_m      [1 x ms]        soil layer thickness [m]
%   soilpar_cell  [N_cells x 1] struct array of Conductivity_Suction
%                 parameters (see Soil_Param_To_Routing); each field
%                 scalar or [1 x ms]
%   dt_h          hydrologic timestep [h]
%   Oice_cell     optional [N_cells x ms] volumetric ice content [m3/m3]
%
% Outputs
%   Q_edge        [N_edges x ms] mean volumetric flux over dt_h [m3/h],
%                 positive from cell_a -> cell_b
%   dtheta_lat    [N_cells x ms] change in theta over dt_h [m3/m3]
%   sub_steps     number of adaptive sub-steps actually taken
%
% Coupling note
%   On a flat field the elevation head is identical in every cell, so
%   gravity cancels horizontally and lateral flux is purely
%   matric-head-driven.  Mass conservation is exact at the finite-volume
%   level (every face is shared by exactly two cells); the sub-step is
%   adaptively rescaled to keep theta within [0, Osat] without ever
%   clipping a volume.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N_cells = topo.N_cells;
N_edges = topo.N_edges;
ms      = numel(Zlayer_m);

if ~(isnumeric(dt_h) && isscalar(dt_h) && isfinite(dt_h) && dt_h > 0)
    error('Darcy_Lateral_Topology:BadTimeStep', ...
        'dt_h must be a positive finite scalar.');
end
if ~(isnumeric(Zlayer_m) && isvector(Zlayer_m) && ...
        all(isfinite(Zlayer_m)) && all(Zlayer_m > 0))
    error('Darcy_Lateral_Topology:BadLayerDepth', ...
        'Zlayer_m must contain positive finite layer thicknesses.');
end
assert(size(theta_cell,1) == N_cells, 'theta_cell row count must equal N_cells');
assert(size(theta_cell,2) == ms,      'theta_cell col count must equal ms');
if any(~isfinite(theta_cell(:)))
    error('Darcy_Lateral_Topology:BadTheta', ...
        'Liquid water content must be finite.');
end
if nargin < 6 || isempty(Oice_cell)
    Oice_cell = zeros(N_cells,ms);
end
if ~isequal(size(Oice_cell),size(theta_cell)) || ...
        any(~isfinite(Oice_cell(:))) || any(Oice_cell(:) < 0)
    error('Darcy_Lateral_Topology:BadIce', ...
        'Oice_cell must match theta_cell and be finite and nonnegative.');
end
if N_edges == 0
    Q_edge = zeros(0,ms);
    dtheta_lat = zeros(N_cells,ms);
    sub_steps = 0;
    return
end

%% --- 1. Per-cell K(theta), head(theta), Osat for every layer ----------
[K_cell, h_cell, Osat_cell, Ohy_cell] = eval_soil_curves( ...
    topo, theta_cell, soilpar_cell, ms);
theta_upper = max(Ohy_cell,Osat_cell-Oice_cell);
tol = 1e-10;
if any(theta_cell(:) < Ohy_cell(:)-tol) || any(theta_cell(:) > theta_upper(:)+tol)
    error('Darcy_Lateral_Topology:BadTheta', ...
        'Liquid water must lie between residual content and ice-limited saturation.');
end

%% --- 2. Volume-based nominal sub-step ---------------------------------
Q_initial = darcy_flux_rates(topo,Zlayer_m,K_cell,h_cell);
incident_rate = zeros(N_cells,ms);
for e = 1:N_edges
    incident_rate(topo.edge_a(e),:) = incident_rate(topo.edge_a(e),:) + abs(Q_initial(e,:));
    incident_rate(topo.edge_b(e),:) = incident_rate(topo.edge_b(e),:) + abs(Q_initial(e,:));
end
mobile_volume = (topo.cell_area_m2(:)*Zlayer_m(:).') .* ...
    max(theta_upper-Ohy_cell,0);
exchange_ratio = max(incident_rate(:)*dt_h ./ max(mobile_volume(:),eps));
dt_sub_h = dt_h / max(1,ceil(2*exchange_ratio));

%% --- 3. Adaptive sub-stepped explicit Darcy update --------------------
theta_w   = theta_cell;
Q_edge    = zeros(N_edges, ms);                % volume accumulator [m3]
elapsed_h = 0;
sub_steps = 0;
dt_sub_h_min = eps*max(dt_h,1);
while elapsed_h < dt_h * (1 - 1e-12)
    dt_try_h = min(dt_sub_h, dt_h - elapsed_h);
    [Q_step, theta_new, scale] = darcy_step( ...
        topo, theta_w, Zlayer_m, K_cell, h_cell, dt_try_h, ...
        Ohy_cell, theta_upper);
    dt_used_h = dt_try_h * scale;
    if dt_used_h < dt_sub_h_min
        error(['Darcy_Lateral_Topology: substep collapsed below safety ' ...
               'floor (dt_used=%.3g h, floor=%.3g h). Check forcing, ' ...
               'soil state, or topology cell sizes.'], dt_used_h, dt_sub_h_min);
    end
    Q_edge    = Q_edge + Q_step * dt_used_h;   % m3
    theta_w   = theta_new;
    elapsed_h = elapsed_h + dt_used_h;
    sub_steps = sub_steps + 1;
    if sub_steps > 100000
        error('Darcy_Lateral_Topology:TooManySubsteps', ...
            'Adaptive lateral flow exceeded 100000 substeps.');
    end
    [K_cell, h_cell] = eval_soil_curves(topo, theta_w, soilpar_cell, ms);
end
%% Time-weighted mean rate over the hydrologic step
Q_edge     = Q_edge / dt_h;
dtheta_lat = theta_w - theta_cell;

end

%% =======================================================================
function [K_cell, h_cell, Osat_cell, Ohy_cell] = eval_soil_curves(topo, theta, soilpar_cell, ms)
% K [mm/h] and matric HEAD h [mm, negative when unsaturated] per cell-layer.
N_cells = topo.N_cells;
K_cell    = zeros(N_cells, ms);
h_cell    = zeros(N_cells, ms);
Osat_cell = zeros(N_cells, ms);
Ohy_cell  = zeros(N_cells, ms);
for c = 1:N_cells
    sp = soilpar_cell(c);
    for k = 1:ms
        scal = @(x) get_layer_par(x, k);
        [Ko, Po] = Conductivity_Suction( ...
            sp.SPAR, scal(sp.Ks), scal(sp.Osat), scal(sp.Ohy), ...
            scal(sp.L), scal(sp.Pe), scal(sp.O33), ...
            scal(sp.alpVG), scal(sp.nVG), scal(sp.lVG), ...
            scal(sp.Ks_mac), scal(sp.Omac), scal(sp.alpVGM), ...
            scal(sp.nVGM), scal(sp.lVGM), scal(sp.Phy1), ...
            scal(sp.s_SVG), scal(sp.bVG), theta(c,k));
        K_cell(c,k) = Ko;
        if sp.SPAR == 2
            h_cell(c,k) = -Po;    % Saxton-Rawls returns suction -> head
        else
            h_cell(c,k) = Po;     % van Genuchten already returns head
        end
        Osat_cell(c,k) = scal(sp.Osat);
        Ohy_cell(c,k)  = scal(sp.Ohy);
    end
end
end

%% =======================================================================
function [Q_step, theta_new, scale] = darcy_step( ...
    topo, theta, Zlayer_m, K_cell, h_cell, dt_h, theta_lower, theta_upper)
% scale in (0,1] is the fraction of dt_h actually applied.  scale<1 means
% the candidate fluxes would have pushed a (cell,layer) outside
% [0, Osat]; Q_step and dV are rescaled uniformly to land exactly on the
% bound.  Mass-conservative: the unused fraction of dt_h is recovered by
% the next iteration of the outer adaptive loop.
N_cells = topo.N_cells;
N_edges = topo.N_edges;
ms      = numel(Zlayer_m);

Q_step = darcy_flux_rates(topo,Zlayer_m,K_cell,h_cell);
bound_tol = 1e-12;
for e = 1:N_edges
    a = topo.edge_a(e); b = topo.edge_b(e);
    for k = 1:ms
        if Q_step(e,k) > 0 && ...
                (theta(a,k) <= theta_lower(a,k)+bound_tol || ...
                 theta(b,k) >= theta_upper(b,k)-bound_tol)
            Q_step(e,k) = 0;
        elseif Q_step(e,k) < 0 && ...
                (theta(b,k) <= theta_lower(b,k)+bound_tol || ...
                 theta(a,k) >= theta_upper(a,k)-bound_tol)
            Q_step(e,k) = 0;
        end
    end
end

% Per-cell volumetric in/out over dt_h
dV_cell = zeros(N_cells, ms);                 % m3
for e = 1:N_edges
    a = topo.edge_a(e); b = topo.edge_b(e);
    for k = 1:ms
        Q = Q_step(e,k) * dt_h;               % m3
        dV_cell(a,k) = dV_cell(a,k) - Q;      % positive Q leaves a
        dV_cell(b,k) = dV_cell(b,k) + Q;
    end
end

% Largest scale in (0,1] that keeps every cell-layer inside the liquid
% residual and ice-limited saturation bounds.
scale = 1;
for c = 1:N_cells
    A_c = topo.cell_area_m2(c);
    for k = 1:ms
        V_layer = A_c * Zlayer_m(k);                       % m3
        dV      = dV_cell(c,k);                            % m3
        if dV > 0
            headroom_V = (theta_upper(c,k) - theta(c,k)) * V_layer;
            if headroom_V <= 0
                scale = 0; break;
            elseif dV > headroom_V
                scale = min(scale, headroom_V / dV);
            end
        elseif dV < 0
            removable_V = (theta(c,k) - theta_lower(c,k)) * V_layer;
            if removable_V <= 0
                scale = 0; break;
            elseif -dV > removable_V
                scale = min(scale, -removable_V / dV);
            end
        end
    end
    if scale == 0; break; end
end

if scale < 1
    dV_cell = dV_cell * scale;
end

theta_new = theta;
for c = 1:N_cells
    A_c = topo.cell_area_m2(c);
    for k = 1:ms
        theta_new(c,k) = theta(c,k) + dV_cell(c,k) / (A_c * Zlayer_m(k));
    end
end

end

%% =======================================================================
function Q = darcy_flux_rates(topo,Zlayer_m,K_cell,h_cell)
N_edges = topo.N_edges;
ms = numel(Zlayer_m);
Q = zeros(N_edges,ms);
for e = 1:N_edges
    if isfield(topo,'edge_active') && ~topo.edge_active(e)
        continue
    end
    a = topo.edge_a(e);
    b = topo.edge_b(e);
    for k = 1:ms
        Ka = K_cell(a,k);
        Kb = K_cell(b,k);
        if Ka <= 0 || Kb <= 0
            continue
        end
        [da,db] = edge_half_distances(topo,e);
        resistance_h = da/(Ka/1000) + db/(Kb/1000); % h
        q = ((h_cell(a,k)-h_cell(b,k))/1000) / resistance_h; % m/h
        A_face = topo.edge_face_width_m(e) * Zlayer_m(k); % m2
        Q(e,k) = q * A_face;                  % m3/h
    end
end
end

%% =======================================================================
function [da,db] = edge_half_distances(topo,e)
if isfield(topo,'edge_face_distance_a_m') && ...
        isfield(topo,'edge_face_distance_b_m')
    da = topo.edge_face_distance_a_m(e);
    db = topo.edge_face_distance_b_m(e);
else
    da = topo.edge_face_distance_m(e)/2;
    db = da;
end
if ~(isfinite(da) && isfinite(db) && da > 0 && db > 0)
    error('Darcy_Lateral_Topology:BadGeometry', ...
        'Every active edge must have positive finite half-cell distances.');
end
end

%% =======================================================================
function v = get_layer_par(x, k)
% Resolve a scalar or per-layer soil parameter to its layer-k value.
if isscalar(x)
    v = x;
elseif numel(x) >= k
    v = x(k);
else
    v = x(end);
end
end
