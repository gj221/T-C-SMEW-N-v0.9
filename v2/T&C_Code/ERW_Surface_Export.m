function [ERW, diag] = ERW_Surface_Export(ERW, topo, runoff_mm, Rd_mm, out_frac, ...
    pond_mass_prev_umol, pond_field_prev_umol, pond_in_water_mm, field_fraction, dt_h)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ERW_Surface_Export
%%
%% Overland export and run-on delivery of the ERW-only aqueous stores
%% (Na, total dissolved Al, residual strong anion An, DIC) with this
%% step's routed surface runoff -- the surface counterpart of
%% ERW_Lateral_Transport, which moves the same stores with the subsurface
%% Darcy fluxes.  Without it, water leaving the field at the outlet would
%% take the T&C cations (Ca, Mg, K, ...) but leave every mole of DIC /
%% alkalinity behind: the exported solution would be charge-imbalanced
%% and the alkalinity-export ledger (the CDR quantity) would miss the
%% overland pathway.
%%
%% TWO SOURCES, exactly as the driver treats the T&C pools:
%%   saturation-excess runoff (Rd_mm) is weathering-zone pore water and
%%     carries the stored aqueous concentration, capped at the dissolved
%%     inventory ERW_Lateral_Transport would allow to move
%%     (min(total, dissolved));
%%   infiltration-excess runoff re-mobilises the PONDED STORE -- the
%%     run-on delivered last step that has not infiltrated -- in the
%%     share out_frac the caller computed against the water that actually
%%     reached the surface (WIS).  out_frac <= 1 by construction, so the
%%     pond limb can never draw on the weathering-zone inventory.
%% The remainder of the pond (1 - out_frac) infiltrated this step and is
%% credited to the stores now; what the router delivers this step goes to
%% the NEXT step's pond store, not into the aqueous inventory, because its
%% water has not infiltrated yet.  Store updates follow ERW_Lateral_
%% Transport's conventions (umol transport, Al stores in nmol, Alk_tot +=
%% dNa - dAn, concentration re-derived over the water volume INCLUDING
%% the water that entered).
%%
%% FIELD-ORIGIN ATTRIBUTION rides along as passive tracer columns and is
%% carried through the pond store, so a stop-over on a partial perimeter
%% cell does not re-label the load with that cell's overlap.
%%
%% INPUTS
%%   runoff_mm            [N x 1] this step's ponded runoff (Rh + Rd) [mm]
%%   Rd_mm                [N x 1] saturation-excess part of it [mm]
%%   out_frac             [N x 1] share of the pond store leaving with Rh
%%   pond_mass_prev_umol  [N x 4] pond store (last step's delivery) [umol/cell]
%%   pond_field_prev_umol [N x 4] its field-origin part
%%   pond_in_water_mm     [N x 1] last step's run-on depth that infiltrated
%%                        this step [mm]
%%   field_fraction       [N x 1] polygon overlap fraction of each cell
%%   dt_h                 hydrologic step [h]
%%
%% OUTPUT diag
%%   .names {'Na','Al','An','DIC'};  .exit_umol, .exit_field_umol [4 x 1]
%%   .delivered_umol, .delivered_field_umol [N x 4]  -> next step's pond
%%   .removed_umol [N x 4] taken from the stores (the Rd limb)
%%   .exit_by_cell_umol [N x 4]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = topo.N_cells;
names = {'Na','Al','An','DIC'};
total_fields = {'Na_tot','Al_tot','An_tot','IC_tot'};
N_sp = numel(names);
concentration_to_umol_m3 = [1000 1 1000 1000];   % stored -> umol/m3
total_to_umol = [1 1e-3 1 1];                     % stored total*A -> umol
umol_to_total = 1./total_to_umol;
al_fields = {'Al','AlOH','AlOH2','AlOH3','AlOH4'};

if ~iscell(ERW) || numel(ERW) ~= N || ...
        any(cellfun(@(x) ~isfield(x,'init') || ~x.init,ERW))
    error('ERW_Surface_Export:NotInitialized', ...
        'Every ERW column must be initialized before surface export.');
end
runoff_mm = max(0,runoff_mm(:)); Rd_mm = max(0,Rd_mm(:));
out_frac = min(1,max(0,out_frac(:)));
field_fraction = field_fraction(:); pond_in_water_mm = max(0,pond_in_water_mm(:));
if ~isequal(numel(runoff_mm),N,numel(Rd_mm),numel(out_frac), ...
        numel(field_fraction),numel(pond_in_water_mm)) || ...
        ~isequal(size(pond_mass_prev_umol),[N N_sp]) || ...
        ~isequal(size(pond_field_prev_umol),[N N_sp])
    error('ERW_Surface_Export:BadSize','Inputs must match the topology.');
end
A = topo.cell_area_m2(:);

diag = struct('names',{names},'exit_umol',zeros(N_sp,1), ...
    'exit_field_umol',zeros(N_sp,1),'delivered_umol',zeros(N,N_sp), ...
    'delivered_field_umol',zeros(N,N_sp),'removed_umol',zeros(N,N_sp), ...
    'exit_by_cell_umol',zeros(N,N_sp));

%% --- dissolved inventory and concentration per cell -----------------
V = zeros(N,1); conc = zeros(N,N_sp); available = zeros(N,N_sp);
pore_volume = zeros(N,1);
al_fraction = zeros(N,numel(al_fields));
for c = 1:N
    E = ERW{c}; row = E.dm_end;
    stored = [E.bg.Na(row), E.bg.Al(row)+E.bg.AlOH(row)+E.bg.AlOH2(row)+ ...
              E.bg.AlOH3(row)+E.bg.AlOH4(row), E.bg.An(row), E.bg.DIC(row)];
    conc(c,:) = stored.*concentration_to_umol_m3;            % umol/m3
    pore_volume(c) = E.const.n*E.const.Zrw*A(c);              % m3
    V(c) = E.current_s*pore_volume(c);                        % m3
    for s = 1:N_sp
        total_mass = E.bg.(total_fields{s})(row)*A(c)*total_to_umol(s);
        if ~(isfinite(total_mass) && total_mass >= 0)
            error('ERW_Surface_Export:NegativeInventory', ...
                'ERW %s inventory is negative in cell %d.',names{s},c);
        end
        available(c,s) = max(0,min(total_mass,conc(c,s)*max(V(c),0)));
    end
    al_state = zeros(1,numel(al_fields));
    for k = 1:numel(al_fields); al_state(k) = E.bg.(al_fields{k})(row); end
    if sum(al_state) > 0
        al_fraction(c,:) = al_state/sum(al_state);
    else
        al_fraction(c,:) = [1 0 0 0 0];
    end
end

%% --- the two limbs [umol per cell] -------------------------------------
soil_out   = min((Rd_mm/1000).*A.*conc, available);       % from the stores
soil_out(~isfinite(soil_out) | soil_out < 0) = 0;
pond_out   = pond_mass_prev_umol.*out_frac;              % leaves again
pond_in    = pond_mass_prev_umol - pond_out;             % infiltrated here
pond_field_out = pond_field_prev_umol.*out_frac;
routed     = soil_out + pond_out;
routed_field = soil_out.*field_fraction + pond_field_out;

if all(runoff_mm <= 0) && all(pond_in(:) <= 0)
    return
end

depth_m = runoff_mm/1000;
C_surf  = routed./max(depth_m.*A,eps);                    % umol/m3
C_field = routed_field./max(depth_m.*A,eps);
C_surf(~isfinite(C_surf) | C_surf < 0) = 0;
C_field(~isfinite(C_field) | C_field < 0) = 0;

%% --- route with the same water as the T&C pools ----------------------
if any(runoff_mm > 0)
    [q_deliver_mm,C_deliver_ext,~,exit_mass_ext,~,~,exit_by_cell_ext] = ...
        Surface_Runoff_Routing(topo,runoff_mm,[C_surf, C_field],dt_h);
    delivered       = (q_deliver_mm/1000).*A.*C_deliver_ext(:,1:N_sp);
    delivered_field = (q_deliver_mm/1000).*A.*C_deliver_ext(:,N_sp+1:end);
    diag.exit_umol         = exit_mass_ext(1:N_sp);
    diag.exit_field_umol   = exit_mass_ext(N_sp+1:end);
    diag.exit_by_cell_umol = exit_by_cell_ext(:,1:N_sp);
else
    delivered = zeros(N,N_sp); delivered_field = zeros(N,N_sp);
end

%% --- update the stores as ERW_Lateral_Transport does ----------------
%% Net change to the inventory this step: what infiltrated out of the
%% pond minus what the Rd limb took.  The routed delivery is NOT added
%% here -- its water has not infiltrated; it becomes the next pond store.
dM = pond_in - soil_out;
mass_tol = 1e-12*max(1,max(available(:)));
for c = 1:N
    row = ERW{c}.dm_end;
    for s = 1:N_sp
        value = ERW{c}.bg.(total_fields{s})(row) + dM(c,s)/A(c)*umol_to_total(s);
        if value < -1e-9*umol_to_total(s)
            error('ERW_Surface_Export:NegativePool', ...
                'Surface export overdraws ERW %s in cell %d.',names{s},c);
        end
        ERW{c}.bg.(total_fields{s})(row) = max(0,value);
    end
    if isfield(ERW{c}.bg,'Alk_tot')
        ERW{c}.bg.Alk_tot(row) = ERW{c}.bg.Alk_tot(row) + (dM(c,1)-dM(c,3))/A(c);
    end
    dissolved = available(c,:) + dM(c,:);
    if any(dissolved < -mass_tol)
        error('ERW_Surface_Export:NegativeInventory', ...
            'Surface export leaves a dissolved ERW inventory negative in cell %d.',c);
    end
    dissolved = max(0,dissolved);
    %% Concentration over the water volume INCLUDING the run-on that
    %% infiltrated this step (ERW_Lateral_Transport's V_after convention);
    %% the saturation follows, capped at the pore volume.
    V_eff = V(c) + pond_in_water_mm(c)/1000*A(c);
    if V_eff <= 0
        %% A dry weathering zone holds the mass in its *_tot inventory
        %% (conserved) with nothing to dissolve it in until the next
        %% vertical step wets the zone.
        new_conc = zeros(1,N_sp);
    else
        new_conc = dissolved./(V_eff*concentration_to_umol_m3);
        ERW{c}.current_s = min(1, V_eff/pore_volume(c));
    end
    ERW{c}.bg.Na(row) = new_conc(1);
    for k = 1:numel(al_fields)
        ERW{c}.bg.(al_fields{k})(row) = new_conc(2)*al_fraction(c,k);
    end
    ERW{c}.bg.An(row) = new_conc(3);
    ERW{c}.bg.DIC(row) = new_conc(4);
end

diag.delivered_umol       = delivered;
diag.delivered_field_umol = delivered_field;
diag.removed_umol         = soil_out;
end
