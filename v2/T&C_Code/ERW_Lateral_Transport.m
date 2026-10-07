function [ERW, diag] = ERW_Lateral_Transport( ...
    ERW, topo, Q_edge, theta_cell, Zlayer_m, active_fraction_cell, dt_h)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lateral advection of ERW-only aqueous stores and soil-gas CO2 diffusion.
%% Transport amounts are umol. ERW stores aqueous Al in nmol units.
%% theta_cell is the pre-lateral state; Q_edge is the mean Darcy rate.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N_cells = topo.N_cells;
N_edges = topo.N_edges;
ms = numel(Zlayer_m);
names = {'Na','Al','An','DIC'};
total_fields = {'Na_tot','Al_tot','An_tot','IC_tot'};
N_sp = numel(names);
concentration_to_umol_m3 = [1000 1 1000 1000];
total_to_umol = [1 1e-3 1 1];
umol_to_total = 1./total_to_umol;

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
if any(cellfun(@(x) ~isfield(x,'init') || ~x.init,ERW))
    error('ERW_Lateral_Transport:NotInitialized', ...
        'Every ERW column must be initialized before lateral transport.');
end

A_cell = topo.cell_area_m2(:);

C_aq = zeros(N_cells,ms,N_sp);                         %% umol/m3
M_aq = zeros(N_cells,ms,N_sp);                         %% umol
M_full = zeros(N_cells,N_sp);                          %% umol, dissolved
water_volume = zeros(N_cells,ms);                      %% m3
V_before = zeros(N_cells,1);                           %% m3
V_vertical = zeros(N_cells,1);                         %% m3
pore_volume = zeros(N_cells,1);                        %% m3
erw_fraction = zeros(N_cells,ms);
for c = 1:N_cells
    erw_fraction(c,:) = min(active_fraction_cell(c,:), ...
        erw_layer_fraction(ERW{c}.const.Zrw,Zlayer_m));
end
[Q_active,~] = Restrict_Lateral_Flux_To_Active_Depth( ...
    topo,Q_edge,erw_fraction);
for c = 1:N_cells
    E = ERW{c}; row = E.dm_end;
    concentration_stored = [E.bg.Na(row), E.bg.Al(row)+E.bg.AlOH(row)+ ...
                     E.bg.AlOH2(row)+E.bg.AlOH3(row)+E.bg.AlOH4(row), ...
                     E.bg.An(row), E.bg.DIC(row)];
    concentration = concentration_stored.*concentration_to_umol_m3;
    pore_volume(c) = E.const.n*E.const.Zrw*A_cell(c);
    V_before(c) = E.current_s*pore_volume(c);
    water_volume(c,:) = theta_cell(c,:).*Zlayer_m(:).'.* ...
        erw_fraction(c,:)*A_cell(c);
    V_vertical(c) = sum(water_volume(c,:));
    if V_before(c) <= 0 || V_vertical(c) <= 0
        continue
    end
    layer_weight = water_volume(c,:)/V_vertical(c);
    for s = 1:N_sp
        dissolved = concentration(s)*V_before(c);
        total_mass = E.bg.(total_fields{s})(row)*A_cell(c)*total_to_umol(s);
        if ~(isfinite(total_mass) && total_mass >= 0)
            error('ERW_Lateral_Transport:NegativeInventory', ...
                'ERW %s inventory is negative in cell %d.',names{s},c);
        end
        available = min(total_mass,dissolved);
        M_full(c,s) = available;
        C_aq(c,:,s) = (available/V_vertical(c)).* ...
            (erw_fraction(c,:) > 0);
        if available > 0
            M_aq(c,:,s) = available*layer_weight;
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
V_after = V_before+sum(dV_cell,2);
sat_tol = 1e-12*max(1,max(pore_volume));
if any(V_before < -sat_tol) || any(V_before-pore_volume > sat_tol) || ...
        any(V_vertical < -sat_tol) || any(V_vertical-pore_volume > sat_tol) || ...
        any(V_after < -sat_tol) || any(V_after-pore_volume > sat_tol)
    error('ERW_Lateral_Transport:BadSaturation', ...
        'ERW water volume lies outside the weathering-zone pore volume.');
end
V_before = min(pore_volume,max(0,V_before));
V_after = min(pore_volume,max(0,V_after));
M_aq_after = M_aq+dM_aq;
mass_tol = 1e-12*max(1,max(M_aq(:)));
if any(M_aq_after(:) < -mass_tol)
    error('ERW_Lateral_Transport:NegativeAqueousMass', ...
        'Lateral transport overdraws an aqueous ERW pool.');
end
dM_area = reshape(sum(dM_aq,2),N_cells,N_sp)./A_cell;
for c = 1:N_cells
    row = ERW{c}.dm_end;
    for s = 1:N_sp
        value = ERW{c}.bg.(total_fields{s})(row)+ ...
            dM_area(c,s)*umol_to_total(s);
        if value < -1e-9*umol_to_total(s)
            error('ERW_Lateral_Transport:NegativePool', ...
                'Lateral transfer overdraws ERW %s in cell %d.',names{s},c);
        end
        if value < 0
            error('ERW_Lateral_Transport:NegativeInventory', ...
                'Lateral transfer leaves ERW %s negative in cell %d.',names{s},c);
        end
        ERW{c}.bg.(total_fields{s})(row) = value;
    end
    if isfield(ERW{c}.bg,'Alk_tot')
        ERW{c}.bg.Alk_tot(row) = ERW{c}.bg.Alk_tot(row) + ...
            dM_area(c,1)-dM_area(c,3);
    end
    dissolved = M_full(c,:) + reshape(sum(dM_aq(c,:,:),2),1,N_sp);
    if any(dissolved < 0)
        error('ERW_Lateral_Transport:NegativeInventory', ...
            'Lateral transfer leaves a dissolved ERW inventory negative in cell %d.',c);
    end
    if V_after(c) <= 0
        if any(dissolved > mass_tol)
            error('ERW_Lateral_Transport:DryAqueousPool', ...
                'Aqueous ERW mass remains in a dry weathering zone.');
        end
        concentration = zeros(1,N_sp);
    else
        concentration = dissolved./ ...
            (V_after(c)*concentration_to_umol_m3);
    end
    al_fields = {'Al','AlOH','AlOH2','AlOH3','AlOH4'};
    al_state = zeros(1,numel(al_fields));
    for k = 1:numel(al_fields)
        al_state(k) = ERW{c}.bg.(al_fields{k})(row);
    end
    al_total = sum(al_state);
    if al_total > 0
        al_fraction = al_state/al_total;
    else
        al_fraction = [1 0 0 0 0];
    end
    ERW{c}.bg.Na(row) = concentration(1);
    for k = 1:numel(al_fields)
        ERW{c}.bg.(al_fields{k})(row) = concentration(2)*al_fraction(k);
    end
    ERW{c}.bg.An(row) = concentration(3);
    ERW{c}.bg.DIC(row) = concentration(4);
end

%% Millington-Quirk gas diffusivity already held in E.bg.D [m2/d].
requested = zeros(N_edges,1);
donor = zeros(N_edges,1);
available_gas = zeros(N_cells,1);
gas_mass_cell = zeros(N_cells,1);
air_volume = zeros(N_cells,1);
gas_concentration = zeros(N_cells,1);
for c = 1:N_cells
    E = ERW{c}; row = E.dm_end;
    air_volume(c) = pore_volume(c)-V_after(c);
    if air_volume(c) < -1e-12*max(1,pore_volume(c))
        error('ERW_Lateral_Transport:NegativeAirVolume', ...
            'Soil-air volume cannot be negative.');
    end
    air_volume(c) = max(0,air_volume(c));
    %% Air pushed out or drawn in by the lateral water is booked once, by
    %% the ADV term of the next ERW step, which sees the same saturation
    %% change; here the soil air keeps its concentration.
    gas_mass = E.bg.CO2_air(row)*1000*air_volume(c);
    if ~(isfinite(gas_mass) && gas_mass >= 0)
        error('ERW_Lateral_Transport:NegativeGas', ...
            'Soil-gas CO2 inventory is negative in cell %d.',c);
    end
    gas_mass_cell(c) = gas_mass;
    if air_volume(c) > 0
        gas_concentration(c) = gas_mass/air_volume(c);
    elseif gas_mass > 1e-12*max(1,gas_mass)
        error('ERW_Lateral_Transport:DryGasVolume', ...
            'Soil-gas CO2 remains in a saturated weathering zone.');
    end
    IC_available = ERW{c}.bg.IC_tot(row)*A_cell(c);
    if ~(isfinite(IC_available) && IC_available >= 0)
        error('ERW_Lateral_Transport:NegativeGas', ...
            'Total inorganic carbon is negative in cell %d.',c);
    end
    available_gas(c) = min(gas_mass,IC_available);
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
    Ca = gas_concentration(a);                            %% umol/m3
    Cb = gas_concentration(b);
    requested(e) = conductance*(Ca-Cb)*(dt_h/24);         %% umol
    %% Explicit pairwise exchange: never transfer more than the amount
    %% that equalizes the two soil-air concentrations, otherwise a large
    %% conductance-to-air-volume ratio overshoots and oscillates.
    equalizing = (Ca-Cb)/(1/air_volume(a)+1/air_volume(b));
    if abs(requested(e)) > abs(equalizing)
        requested(e) = equalizing;
    end
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
edge_transfer = zeros(N_edges,1);
for e = 1:N_edges
    if donor(e) == 0; continue; end
    dM = requested(e)*scale(donor(e));
    edge_transfer(e) = dM;
    a = topo.edge_a(e); b = topo.edge_b(e);
    dIC(a) = dIC(a)-dM;
    dIC(b) = dIC(b)+dM;
end
open_air = air_volume > 0;
if any(open_air)
    Cmin = min(gas_concentration(open_air));
    Cmax = max(gas_concentration(open_air));
    dC = zeros(N_cells,1);
    dC(open_air) = dIC(open_air)./air_volume(open_air);
    bounds = 1;
    for c = find(open_air).'
        if dC(c) > 0
            bounds = min(bounds,(Cmax-gas_concentration(c))/dC(c));
        elseif dC(c) < 0
            bounds = min(bounds,(gas_concentration(c)-Cmin)/(-dC(c)));
        end
    end
    bounds = min(1,max(0,bounds));
    edge_transfer = edge_transfer*bounds;
    dIC = dIC*bounds;
end
diag.CO2_edge_umol = edge_transfer;
for c = 1:N_cells
    row = ERW{c}.dm_end;
    ERW{c}.bg.IC_tot(row) = ERW{c}.bg.IC_tot(row)+dIC(c)/A_cell(c);
    if air_volume(c) > 0
        CO2_air = (gas_mass_cell(c)+dIC(c))/(1000*air_volume(c));
        if CO2_air < -1e-9
            error('ERW_Lateral_Transport:NegativeGas', ...
                'Lateral diffusion overdraws soil-gas CO2 in cell %d.',c);
        end
        if CO2_air < 0
            error('ERW_Lateral_Transport:NegativeGas', ...
                'Lateral diffusion leaves soil-gas CO2 negative in cell %d.',c);
        end
        ERW{c}.bg.CO2_air(row) = CO2_air;
    elseif abs(dIC(c)) > 1e-9
        error('ERW_Lateral_Transport:DryGasVolume', ...
            'CO2 cannot diffuse into a cell with zero soil-air volume.');
    else
        ERW{c}.bg.CO2_air(row) = 0;
    end
    ERW{c}.current_s = V_after(c)/pore_volume(c);
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

function frac = erw_layer_fraction(Zrw,Zlayer_m)
%% Fraction of each soil layer lying inside the ERW weathering depth [m].
if ~(isscalar(Zrw) && isfinite(Zrw) && Zrw > 0)
    error('ERW_Lateral_Transport:BadZrw', ...
        'The ERW weathering depth must be a positive finite scalar.');
end
thickness = Zlayer_m(:).';
top = [0 cumsum(thickness)];
frac = max(0,min(Zrw,top(2:end))-top(1:end-1))./thickness;
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
    error('ERW_Lateral_Transport:BadGeometry', ...
        'Every active edge must have positive finite half-cell distances.');
end
end
