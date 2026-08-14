function [T_cell, diag] = Soil_Lateral_Heat_Transport( ...
    topo, Q_edge, theta_before, theta_after, T_cell, lambda_cell, ...
    cv_cell, Zlayer_m, dt_h)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lateral soil conduction and liquid-water sensible heat advection.
%% Energy is referenced to 0 degC. The closed perimeter has zero heat flux.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N_cells = topo.N_cells;
N_edges = topo.N_edges;
ms = numel(Zlayer_m);
state_size = [N_cells ms];

if ~isequal(size(theta_before),state_size) || ...
        ~isequal(size(theta_after),state_size) || ...
        ~isequal(size(T_cell),state_size) || ...
        ~isequal(size(lambda_cell),state_size) || ...
        ~isequal(size(cv_cell),state_size) || ...
        ~isequal(size(Q_edge),[N_edges ms])
    error('Soil_Lateral_Heat_Transport:BadSize', ...
        'State arrays must be N_cells by N_layers and Q_edge must match the topology.');
end
if ~(isscalar(dt_h) && isfinite(dt_h) && dt_h > 0) || ...
        any(~isfinite(theta_before(:))) || any(~isfinite(theta_after(:))) || ...
        any(~isfinite(T_cell(:))) || any(~isfinite(lambda_cell(:))) || ...
        any(lambda_cell(:) < 0) || any(~isfinite(cv_cell(:))) || any(cv_cell(:) <= 0)
    error('Soil_Lateral_Heat_Transport:BadState', ...
        'Thermal states must be finite, with nonnegative conductivity and positive heat capacity.');
end

diag.advective_edge_J = zeros(N_edges,ms);
diag.conductive_edge_J = zeros(N_edges,ms);
diag.sub_steps = 0;
diag.energy_residual_J = 0;
if N_edges == 0
    return
end

rho_cw = 1000*4186;                                      %% J/m3/K
bulk_volume = topo.cell_area_m2(:)*Zlayer_m(:).';        %% m3
C_before = cv_cell.*bulk_volume;                         %% J/K
dV = (theta_after-theta_before).*bulk_volume;            %% m3
C_after = C_before + rho_cw*dV;                          %% J/K
if any(C_after(:) <= 0)
    error('Soil_Lateral_Heat_Transport:BadCapacity', ...
        'Lateral water exchange produced nonpositive layer heat capacity.');
end

T_initial = T_cell;
dE_adv = zeros(N_cells,ms);
for e = 1:N_edges
    if isfield(topo,'edge_active') && ~topo.edge_active(e); continue; end
    a = topo.edge_a(e); b = topo.edge_b(e);
    for k = 1:ms
        Q = Q_edge(e,k);
        if Q > 0
            donor_T = T_initial(a,k);
        elseif Q < 0
            donor_T = T_initial(b,k);
        else
            continue
        end
        dE = rho_cw*Q*dt_h*donor_T;
        diag.advective_edge_J(e,k) = dE;
        dE_adv(a,k) = dE_adv(a,k)-dE;
        dE_adv(b,k) = dE_adv(b,k)+dE;
    end
end
T_cell = (C_before.*T_initial+dE_adv)./C_after;

G_edge = zeros(N_edges,ms);                              %% W/K
incident_G = zeros(N_cells,ms);
for e = 1:N_edges
    if isfield(topo,'edge_active') && ~topo.edge_active(e); continue; end
    a = topo.edge_a(e); b = topo.edge_b(e);
    [da,db] = edge_half_distances(topo,e);
    for k = 1:ms
        la = lambda_cell(a,k); lb = lambda_cell(b,k);
        if la <= 0 || lb <= 0; continue; end
        A_face = topo.edge_face_width_m(e)*Zlayer_m(k);
        G_edge(e,k) = A_face/(da/la+db/lb);
        incident_G(a,k) = incident_G(a,k)+G_edge(e,k);
        incident_G(b,k) = incident_G(b,k)+G_edge(e,k);
    end
end

dt_s = dt_h*3600;
thermal_ratio = max(incident_G(:)*dt_s./C_after(:));
diag.sub_steps = max(1,ceil(thermal_ratio/0.45));
dt_sub_s = dt_s/diag.sub_steps;
for it = 1:diag.sub_steps
    dE_cond = zeros(N_cells,ms);
    for e = 1:N_edges
        a = topo.edge_a(e); b = topo.edge_b(e);
        for k = 1:ms
            dE = G_edge(e,k)*(T_cell(a,k)-T_cell(b,k))*dt_sub_s;
            diag.conductive_edge_J(e,k) = diag.conductive_edge_J(e,k)+dE;
            dE_cond(a,k) = dE_cond(a,k)-dE;
            dE_cond(b,k) = dE_cond(b,k)+dE;
        end
    end
    T_cell = T_cell+dE_cond./C_after;
end

diag.energy_residual_J = sum(C_after(:).*T_cell(:)) - ...
    sum(C_before(:).*T_initial(:));

end

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
    error('Soil_Lateral_Heat_Transport:BadGeometry', ...
        'Every active edge must have positive finite half-cell distances.');
end
end
