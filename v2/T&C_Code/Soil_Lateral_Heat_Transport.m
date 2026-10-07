function [T_cell, diag, Oliq_cell, Oice_cell] = Soil_Lateral_Heat_Transport( ...
    topo, Q_edge, theta_before, theta_after, T_cell, lambda_cell, ...
    capacity_cell, Zlayer_m, dt_h, phase)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lateral soil conduction and liquid-water sensible heat advection.
%% Energy is referenced to liquid water and soil at 0 degC. The closed
%% perimeter has zero heat flux.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N_cells = topo.N_cells;
N_edges = topo.N_edges;
ms = numel(Zlayer_m);
state_size = [N_cells ms];

if ~isequal(size(theta_before),state_size) || ...
        ~isequal(size(theta_after),state_size) || ...
        ~isequal(size(T_cell),state_size) || ...
        ~isequal(size(lambda_cell),state_size) || ...
        ~isequal(size(capacity_cell),state_size) || ...
        ~isequal(size(Q_edge),[N_edges ms]) || ...
        numel(topo.cell_area_m2) ~= N_cells
    error('Soil_Lateral_Heat_Transport:BadSize', ...
        'State arrays must be N_cells by N_layers and Q_edge must match the topology.');
end
if ~(isscalar(dt_h) && isfinite(dt_h) && dt_h > 0) || ...
        any(~isfinite(Zlayer_m(:))) || any(Zlayer_m(:) <= 0) || ...
        any(~isfinite(topo.cell_area_m2(:))) || any(topo.cell_area_m2(:) <= 0) || ...
        any(~isfinite(Q_edge(:))) || ...
        any(~isfinite(theta_before(:))) || any(theta_before(:) < 0) || ...
        any(~isfinite(theta_after(:))) || any(theta_after(:) < 0) || ...
        any(~isfinite(T_cell(:))) || any(T_cell(:) <= -273.15) || ...
        any(~isfinite(lambda_cell(:))) || any(lambda_cell(:) < 0) || ...
        any(~isfinite(capacity_cell(:))) || any(capacity_cell(:) <= 0)
    error('Soil_Lateral_Heat_Transport:BadState', ...
        'Thermal states must be finite, with physical water, temperature, conductivity and capacity.');
end

rho_cw = 1000*4186;                                     %% J/m3/K
rho_i_Lf = 916.2*333700;                                %% J/m3 ice
bulk_volume = topo.cell_area_m2(:)*Zlayer_m(:).';       %% m3
T_initial = T_cell;

diag.advective_edge_J = zeros(N_edges,ms);
diag.conductive_edge_J = zeros(N_edges,ms);
diag.sub_steps = 0;
diag.energy_residual_J = 0;
diag.energy_scale_J = 1;
diag.energy_residual_relative = 0;
Oliq_cell = theta_after;
if nargin < 10 || isempty(phase)
    Oice_cell = zeros(N_cells,ms);
else
    Oice_cell = phase.Oice_cell;
end
if N_edges == 0
    return
end

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

use_phase = nargin >= 10 && ~isempty(phase);
if use_phase
    required = {'Oice_cell','parameter_cell','Pre','OPT_FR_SOIL'};
    if ~isstruct(phase) || ~isscalar(phase) || ...
            ~all(isfield(phase,required)) || ...
            ~isequal(size(phase.Oice_cell),state_size) || ...
            any(~isfinite(phase.Oice_cell(:))) || any(phase.Oice_cell(:) < 0) || ...
            ~iscell(phase.parameter_cell) || numel(phase.parameter_cell) ~= N_cells || ...
            ~(isscalar(phase.Pre) && isfinite(phase.Pre) && phase.Pre > 0) || ...
            ~(isscalar(phase.OPT_FR_SOIL) && isfinite(phase.OPT_FR_SOIL))
        error('Soil_Lateral_Heat_Transport:BadPhaseState', ...
            'Phase transport requires ice, parameters, pressure and the frozen-soil option.');
    end
    %% T&C owns the phase split, and its state is deliberately NOT on the
    %% local freezing curve: HYDROLOGIC_UNIT solves Richards flow and then
    %% Cryosuction_stabilizer moves liquid water BETWEEN layers while
    %% holding Oice fixed. So (T,theta,Oice) cannot satisfy a local
    %% equilibrium relation, at any tolerance, and re-deriving the split
    %% here would re-partition roughly 0.03 m3/m3 every step and inject
    %% spurious latent heat.
    %%
    %% This lateral call is an operator-split sub-step: it redistributes
    %% heat and leaves phase change to the vertical solver. Ice is held at
    %% T&C's prognostic value and the latent term rides in the apparent
    %% heat capacity, which already embeds dOliq/dT along the freezing
    %% curve, so temperature still responds with the right thermal inertia
    %% near the freezing front.
    %% Instead, carry that departure from equilibrium as a constant offset
    %% across the lateral step. Because the offset is fixed, the enthalpy
    %% curve is merely shifted by -rho_i_Lf*offset, so the same inversion
    %% applies to a shifted target and latent heat is still resolved. When
    %% the incoming state IS on the curve the offset is zero and this
    %% reduces exactly to the pure equilibrium formulation.
    Otot_before = theta_before+phase.Oice_cell;
    ice_offset = zeros(N_cells,ms);
    H_initial = zeros(N_cells,ms);
    for c = 1:N_cells
        [~,~,ice_eq,~,cv] = Soil_Lateral_Thermal_State( ...
            T_initial(c,:),Otot_before(c,:),phase.Pre, ...
            phase.parameter_cell{c},phase.OPT_FR_SOIL);
        ice_offset(c,:) = phase.Oice_cell(c,:)-ice_eq;
        H_initial(c,:) = cv.*T_initial(c,:)-rho_i_Lf*phase.Oice_cell(c,:);
    end
    Otot_after = theta_after+phase.Oice_cell;
    E_cell = H_initial.*bulk_volume+dE_adv;
    T_guess = T_initial+dE_adv./(capacity_cell.*bulk_volume);
    [T_cell,Oliq_cell,Oice_cell,capp_cell] = invert_enthalpy_offset( ...
        E_cell./bulk_volume,T_guess,Otot_after,phase,rho_i_Lf,ice_offset);
    stability_capacity_cell = capp_cell;
else
    cv_before = capacity_cell;
    cv_after = cv_before+rho_cw*(theta_after-theta_before);
    if any(cv_after(:) <= 0)
        error('Soil_Lateral_Heat_Transport:BadCapacity', ...
            'Lateral water exchange produced nonpositive sensible heat capacity.');
    end
    H_initial = cv_before.*T_initial;
    E_cell = H_initial.*bulk_volume+dE_adv;
    T_cell = E_cell./(cv_after.*bulk_volume);
    sensible_capacity_cell = cv_after;
    stability_capacity_cell = sensible_capacity_cell;
end

G_edge = zeros(N_edges,ms);                             %% W/K
incident_G = zeros(N_cells,ms);
for e = 1:N_edges
    if isfield(topo,'edge_active') && ~topo.edge_active(e); continue; end
    a = topo.edge_a(e); b = topo.edge_b(e);
    [da,db] = edge_half_distances(topo,e);
    for k = 1:ms
        la = lambda_cell(a,k); lb = lambda_cell(b,k);
        if la <= 0 || lb <= 0; continue; end
        A_face = topo.edge_face_width_m(e)*Zlayer_m(k);
        if ~(isfinite(A_face) && A_face > 0)
            error('Soil_Lateral_Heat_Transport:BadGeometry', ...
                'Every active edge must have positive finite face area.');
        end
        G_edge(e,k) = A_face/(da/la+db/lb);
        incident_G(a,k) = incident_G(a,k)+G_edge(e,k);
        incident_G(b,k) = incident_G(b,k)+G_edge(e,k);
    end
end

remaining_s = dt_h*3600;
while remaining_s > 0
    rate = max(incident_G(:)./(stability_capacity_cell(:).*bulk_volume(:)));
    if rate > 0
        dt_sub_s = min(remaining_s,0.45/rate);
    else
        dt_sub_s = remaining_s;
    end
    if remaining_s-dt_sub_s <= 16*eps(remaining_s)
        dt_sub_s = remaining_s;
    end
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
    E_cell = E_cell+dE_cond;
    if use_phase
        [T_cell,Oliq_cell,Oice_cell,capp_cell] = invert_enthalpy_offset( ...
            E_cell./bulk_volume,T_cell,Otot_after,phase,rho_i_Lf,ice_offset);
        stability_capacity_cell = capp_cell;
    else
        T_cell = E_cell./(sensible_capacity_cell.*bulk_volume);
    end
    remaining_s = remaining_s-dt_sub_s;
    diag.sub_steps = diag.sub_steps+1;
    if diag.sub_steps > 100000
        error('Soil_Lateral_Heat_Transport:TooManySubsteps', ...
            'Lateral heat transport exceeded the stability sub-step limit.');
    end
end

if use_phase
    %% Final enthalpy uses the same offset ice the inversion returned, so
    %% the ledger closes on the state actually handed back to T&C.
    H_final = zeros(N_cells,ms);
    for c = 1:N_cells
        [~,~,~,~,cv] = Soil_Lateral_Thermal_State( ...
            T_cell(c,:),Otot_after(c,:),phase.Pre, ...
            phase.parameter_cell{c},phase.OPT_FR_SOIL);
        H_final(c,:) = cv.*T_cell(c,:)-rho_i_Lf*Oice_cell(c,:);
    end
    if max(abs(Oliq_cell(:)+Oice_cell(:)-Otot_after(:))) > 1e-10
        error('Soil_Lateral_Heat_Transport:WaterPhaseResidual', ...
            'Post-transport liquid and ice do not conserve total soil water.');
    end
else
    H_final = sensible_capacity_cell.*T_cell;
end
E_initial = H_initial.*bulk_volume;
E_final = H_final.*bulk_volume;
diag.energy_residual_J = sum(E_final(:))-sum(E_initial(:));
diag.energy_scale_J = max([sum(abs(E_initial(:))),sum(abs(E_final(:))), ...
    sum(abs(diag.advective_edge_J(:))), ...
    sum(abs(diag.conductive_edge_J(:))),1]);
diag.energy_residual_relative = abs(diag.energy_residual_J)/diag.energy_scale_J;

end

function [T,Oliq,Oice,capp] = invert_enthalpy_offset( ...
    H_target,T_guess,Otot,phase,rho_i_Lf,ice_offset)
%% Invert enthalpy on the equilibrium curve shifted by a fixed ice offset.
%% H_tc(T) = H_eq(T) - rho_i_Lf*offset, so solving H_tc(T) = H_target is
%% solving H_eq(T) = H_target + rho_i_Lf*offset. The offset is then added
%% back to the equilibrium ice and the liquid follows from total water.
[T,~,Oice_eq,capp] = invert_enthalpy_cells( ...
    H_target+rho_i_Lf*ice_offset,T_guess,Otot,phase,rho_i_Lf);
Oice = min(Otot,max(0,Oice_eq+ice_offset));
Oliq = Otot-Oice;
if any(~isfinite(Oice(:))) || any(Oliq(:) < -1e-12)
    error('Soil_Lateral_Heat_Transport:BadOffsetPhase', ...
        'Offset phase update produced a nonphysical liquid or ice content.');
end
end

function [T,Oliq,Oice,capp] = invert_enthalpy_cells( ...
    H_target,T_guess,Otot,phase,rho_i_Lf)

[N_cells,ms] = size(H_target);
T = T_guess;
Oliq = zeros(N_cells,ms);
Oice = zeros(N_cells,ms);
capp = zeros(N_cells,ms);
for c = 1:N_cells
    [T(c,:),Oliq(c,:),Oice(c,:),capp(c,:)] = invert_enthalpy_column( ...
        H_target(c,:),T_guess(c,:),Otot(c,:),phase.Pre, ...
        phase.parameter_cell{c},phase.OPT_FR_SOIL,rho_i_Lf);
end

end

function [T,Oliq,Oice,capp] = invert_enthalpy_column( ...
    H_target,T_guess,Otot,Pre,p,OPT_FR_SOIL,rho_i_Lf)

T = T_guess;
[H,capp,Oice,Oliq] = enthalpy_state( ...
    T,Otot,Pre,p,OPT_FR_SOIL,rho_i_Lf);
scale = max(max(abs(H_target),rho_i_Lf*Otot),1);
residual = H-H_target;
if ~any(abs(residual) > 1e-12*scale); return; end

T_lo = max(-273.14,T_guess-200);
T_hi = T_guess+200;
[H_lo,~,~,~] = enthalpy_state(T_lo,Otot,Pre,p,OPT_FR_SOIL,rho_i_Lf);
[H_hi,~,~,~] = enthalpy_state(T_hi,Otot,Pre,p,OPT_FR_SOIL,rho_i_Lf);
if any(H_target < H_lo-1e-12*scale) || any(H_target > H_hi+1e-12*scale)
    error('Soil_Lateral_Heat_Transport:UnbracketedEnthalpy', ...
        'The lateral energy state is outside the physical temperature bracket.');
end

for it = 1:60
    active = abs(residual) > 1e-12*scale;
    if ~any(active); return; end
    above = residual > 0;
    T_hi(above & active) = T(above & active);
    T_lo(~above & active) = T(~above & active);
    T_new = T-residual./capp;
    use_bisection = ~isfinite(T_new) | T_new <= T_lo | T_new >= T_hi;
    T_new(use_bisection) = 0.5*(T_lo(use_bisection)+T_hi(use_bisection));
    T(active) = T_new(active);
    [H,capp,Oice,Oliq] = enthalpy_state( ...
        T,Otot,Pre,p,OPT_FR_SOIL,rho_i_Lf);
    residual = H-H_target;
end

if any(abs(residual) > 1e-10*scale)
    error('Soil_Lateral_Heat_Transport:EnthalpyConvergence', ...
        'The phase-equilibrium enthalpy inversion did not converge.');
end

end

function [H,capp,Oice,Oliq] = enthalpy_state( ...
    T,Otot,Pre,p,OPT_FR_SOIL,rho_i_Lf)

[~,capp,Oice,Oliq,cv] = Soil_Lateral_Thermal_State( ...
    T,Otot,Pre,p,OPT_FR_SOIL);
H = cv.*T-rho_i_Lf*Oice;

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
