function [ERW, diag] = ERW_Lateral_Transport( ...
    ERW, topo, Q_edge, theta_cell, Zlayer_m, active_fraction_cell, dt_h)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lateral advection of ERW-only aqueous stores and soil-gas CO2 diffusion.
%% Totals are stored by ERW in umol/m2; concentrations are umol/L.
%% theta_cell is the start-of-step state; Q_edge is the mean Darcy rate.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N_cells = topo.N_cells;
N_edges = topo.N_edges;
ms = numel(Zlayer_m);
names = {'Na','Al','An','DIC'};
total_fields = {'Na_tot','Al_tot','An_tot','IC_tot'};
N_sp = numel(names);

if ~iscell(ERW) || numel(ERW) ~= N_cells || ...
        ~isequal(size(Q_edge),[N_edges ms]) || ...
        ~isequal(size(theta_cell),[N_cells ms]) || ...
        ~isequal(size(active_fraction_cell),[N_cells ms])
    error('ERW_Lateral_Transport:BadSize', ...
        'ERW and state arrays must match the topology and soil layers.');
end
if ~(isscalar(dt_h) && isfinite(dt_h) && dt_h > 0) || ...
        any(~isfinite(theta_cell(:))) || any(theta_cell(:) < 0) || ...
        any(~isfinite(active_fraction_cell(:))) || ...
        any(active_fraction_cell(:) < 0) || any(active_fraction_cell(:) > 1)
    error('ERW_Lateral_Transport:BadState', ...
        'Water state and active-layer fractions must be finite and physical.');
end

diag.aqueous_names = names;
diag.aqueous_edge_umol = zeros(N_edges,N_sp);
diag.CO2_edge_umol = zeros(N_edges,1);
diag.buf_to_main.aqueous_umol = zeros(N_sp,1);
diag.buf_to_main.CO2_umol = 0;
diag.mass_residual_umol = zeros(N_sp,1);
if N_edges == 0
    return
end
if any(cellfun(@(x) ~isfield(x,'init') || ~x.init,ERW))
    error('ERW_Lateral_Transport:NotInitialized', ...
        'Every ERW column must be initialized before lateral transport.');
end

A_cell = topo.cell_area_m2(:);
[Q_active,~] = Restrict_Lateral_Flux_To_Active_Depth( ...
    topo,Q_edge,active_fraction_cell);

C_aq = zeros(N_cells,ms,N_sp);                         %% umol/m3
M_aq = zeros(N_cells,ms,N_sp);                         %% umol
water_volume = zeros(N_cells,ms);                      %% m3
for c = 1:N_cells
    E = ERW{c}; row = E.dm_end;
    concentration = [E.bg.Na(row), E.bg.Al(row)+E.bg.AlOH4(row), ...
                     E.bg.An(row), E.bg.DIC(row)]*1000;
    water_volume(c,:) = theta_cell(c,:).*Zlayer_m(:).'.* ...
        active_fraction_cell(c,:)*A_cell(c);
    for s = 1:N_sp
        C_aq(c,:,s) = concentration(s);
        aqueous_mass = concentration(s)*water_volume(c,:);
        total_mass = max(0,E.bg.(total_fields{s})(row)*A_cell(c));
        available = min(total_mass,sum(aqueous_mass));
        if available > 0 && sum(aqueous_mass) > 0
            M_aq(c,:,s) = aqueous_mass*(available/sum(aqueous_mass));
        end
    end
end

[dM_aq,dM_edge] = Solute_Lateral_Transport( ...
    topo,Q_active,C_aq,M_aq,dt_h);
diag.aqueous_edge_umol = reshape(sum(dM_edge,2),N_edges,N_sp);
dV_cell = zeros(N_cells,ms);
for e = 1:N_edges
    dV = Q_active(e,:)*dt_h;
    a = topo.edge_a(e); b = topo.edge_b(e);
    dV_cell(a,:) = dV_cell(a,:)-dV;
    dV_cell(b,:) = dV_cell(b,:)+dV;
end
water_volume_after = water_volume+dV_cell;
volume_tol = 1e-12*max(1,max(water_volume(:)));
if any(water_volume_after(:) < -volume_tol)
    error('ERW_Lateral_Transport:NegativeWater', ...
        'Lateral transport overdraws active-layer water.');
end
water_volume_after = max(0,water_volume_after);
M_aq_after = M_aq+dM_aq;
mass_tol = 1e-12*max(1,max(M_aq(:)));
if any(M_aq_after(:) < -mass_tol)
    error('ERW_Lateral_Transport:NegativeAqueousMass', ...
        'Lateral transport overdraws an aqueous ERW pool.');
end
M_aq_after = max(0,M_aq_after);
dM_area = reshape(sum(dM_aq,2),N_cells,N_sp)./A_cell;
for c = 1:N_cells
    row = ERW{c}.dm_end;
    for s = 1:N_sp
        value = ERW{c}.bg.(total_fields{s})(row)+dM_area(c,s);
        if value < -1e-9
            error('ERW_Lateral_Transport:NegativePool', ...
                'Lateral transfer overdraws ERW %s in cell %d.',names{s},c);
        end
        ERW{c}.bg.(total_fields{s})(row) = max(0,value);
    end
    if isfield(ERW{c}.bg,'Alk_tot')
        ERW{c}.bg.Alk_tot(row) = ERW{c}.bg.Alk_tot(row) + ...
            dM_area(c,1)-dM_area(c,3);
    end
    volume = sum(water_volume_after(c,:));
    dissolved = reshape(sum(M_aq_after(c,:,:),2),1,N_sp);
    if volume == 0
        if any(dissolved > mass_tol)
            error('ERW_Lateral_Transport:DryAqueousPool', ...
                'Aqueous ERW mass remains in a dry active layer.');
        end
        concentration = zeros(1,N_sp);
    else
        concentration = dissolved/(1000*volume);          %% umol/L
    end
    al_total = ERW{c}.bg.Al(row)+ERW{c}.bg.AlOH4(row);
    if al_total > 0
        al_fraction = ERW{c}.bg.Al(row)/al_total;
    else
        al_fraction = 1;
    end
    ERW{c}.bg.Na(row) = concentration(1);
    ERW{c}.bg.Al(row) = concentration(2)*al_fraction;
    ERW{c}.bg.AlOH4(row) = concentration(2)*(1-al_fraction);
    ERW{c}.bg.An(row) = concentration(3);
    ERW{c}.bg.DIC(row) = concentration(4);
end

%% Millington-Quirk gas diffusivity already held in E.bg.D [m2/d].
requested = zeros(N_edges,1);
donor = zeros(N_edges,1);
available_gas = zeros(N_cells,1);
air_volume = zeros(N_cells,1);
for c = 1:N_cells
    E = ERW{c}; row = E.dm_end;
    s = min(1,max(0,E.current_s));
    air_volume(c) = E.const.n*(1-s)*E.const.Zrw*A_cell(c);
    gas_mass = max(0,E.bg.CO2_air(row)*1000*air_volume(c));
    available_gas(c) = min(gas_mass,max(0,E.bg.IC_tot(row)*A_cell(c)));
end
for e = 1:N_edges
    if isfield(topo,'edge_active') && ~topo.edge_active(e); continue; end
    a = topo.edge_a(e); b = topo.edge_b(e);
    if air_volume(a) <= 0 || air_volume(b) <= 0; continue; end
    Ea = ERW{a}; Eb = ERW{b};
    ra = Ea.dm_end; rb = Eb.dm_end;
    Da = Ea.bg.D(ra); Db = Eb.bg.D(rb);
    if Da <= 0 || Db <= 0; continue; end
    [da,db] = edge_half_distances(topo,e);
    face_area = topo.edge_face_width_m(e)*min(Ea.const.Zrw,Eb.const.Zrw);
    conductance = face_area/(da/Da+db/Db);                %% m3/d
    Ca = Ea.bg.CO2_air(ra)*1000;                          %% umol/m3
    Cb = Eb.bg.CO2_air(rb)*1000;
    requested(e) = conductance*(Ca-Cb)*(dt_h/24);         %% umol
    if requested(e) > 0
        donor(e) = a;
    elseif requested(e) < 0
        donor(e) = b;
    end
end

requested_out = zeros(N_cells,1);
for e = 1:N_edges
    if donor(e) > 0
        requested_out(donor(e)) = requested_out(donor(e))+abs(requested(e));
    end
end
scale = ones(N_cells,1);
limited = requested_out > available_gas;
scale(limited) = available_gas(limited)./requested_out(limited);
dIC = zeros(N_cells,1);
for e = 1:N_edges
    if donor(e) == 0; continue; end
    dM = requested(e)*scale(donor(e));
    a = topo.edge_a(e); b = topo.edge_b(e);
    diag.CO2_edge_umol(e) = dM;
    dIC(a) = dIC(a)-dM;
    dIC(b) = dIC(b)+dM;
end
for c = 1:N_cells
    row = ERW{c}.dm_end;
    ERW{c}.bg.IC_tot(row) = ERW{c}.bg.IC_tot(row)+dIC(c)/A_cell(c);
    if air_volume(c) > 0
        CO2_air = ERW{c}.bg.CO2_air(row)+dIC(c)/(1000*air_volume(c));
        if CO2_air < -1e-9
            error('ERW_Lateral_Transport:NegativeGas', ...
                'Lateral diffusion overdraws soil-gas CO2 in cell %d.',c);
        end
        ERW{c}.bg.CO2_air(row) = max(0,CO2_air);
    elseif abs(dIC(c)) > 1e-9
        error('ERW_Lateral_Transport:DryGasVolume', ...
            'CO2 cannot diffuse into a cell with zero soil-air volume.');
    end
end

bm_sign = edge_sign_buffer_to_main(topo);
diag.buf_to_main.aqueous_umol = ...
    (bm_sign.'*diag.aqueous_edge_umol).';
diag.buf_to_main.CO2_umol = bm_sign.'*diag.CO2_edge_umol;
diag.mass_residual_umol = reshape(sum(dM_aq,[1 2]),N_sp,1);
diag.mass_residual_umol(4) = diag.mass_residual_umol(4)+sum(dIC);

end

function bm_sign = edge_sign_buffer_to_main(topo)
is_main_a = strcmp(topo.cell_region(topo.edge_a),'main');
is_main_b = strcmp(topo.cell_region(topo.edge_b),'main');
bm_sign = zeros(topo.N_edges,1);
bm_sign(~is_main_a & is_main_b) = 1;
bm_sign(is_main_a & ~is_main_b) = -1;
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
    error('ERW_Lateral_Transport:BadGeometry', ...
        'Every active edge must have positive finite half-cell distances.');
end
end
