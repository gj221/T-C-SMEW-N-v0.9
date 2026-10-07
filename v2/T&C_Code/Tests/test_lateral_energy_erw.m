function test_lateral_energy_erw()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Unequal-cell resistance, physical water bounds, heat, and ERW transport.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n=== test_lateral_energy_erw ===\n');
here = fileparts(mfilename('fullpath'));
code = fileparts(here);
addpath(code,fullfile(code,'Inputs'),fullfile(code,'Parameters'));
n_pass = 0; n_fail = 0;

    function check(name, cond, detail)
        if nargin < 3; detail = ''; end
        if cond
            n_pass = n_pass + 1;
            fprintf('  PASS  %s\n',name);
        else
            n_fail = n_fail + 1;
            fprintf('  FAIL  %s  %s\n',name,detail);
        end
    end

    function check_error(name, fcn, expected_id)
        try
            fcn();
            check(name,false,'expected an error');
        catch ME
            check(name,strcmp(ME.identifier,expected_id),ME.identifier);
        end
    end

%% Unequal cells retain each half-cell distance at their shared face.
topo = Build_Cell_Topology(60,5,2,2);
check('edge half distances sum to the centre distance', ...
    max(abs(topo.edge_face_distance_a_m + topo.edge_face_distance_b_m - ...
    topo.edge_face_distance_m)) < 1e-12);
is_bm = xor(strcmp(topo.cell_region(topo.edge_a),'main'), ...
            strcmp(topo.cell_region(topo.edge_b),'main'));
bm_halves = sort([topo.edge_face_distance_a_m(is_bm), ...
                  topo.edge_face_distance_b_m(is_bm)],2);
check('retained 25 m and 5 m cells use 12.5 m and 2.5 m resistances', ...
    all(abs(bm_halves(:,1)-2.5) < 1e-12) && ...
    all(abs(bm_halves(:,2)-12.5) < 1e-12));

%% Darcy transmissibility uses the two unequal half-cell resistances.
[sp, Osat, Ohy] = soil_parameters(1);
sp2 = sp;
sp2.Ks = 0.35*sp.Ks;
soil = [sp;sp2];
two = two_cell_topology(2,8,1,[20;80]);
theta = [0.36;0.24];
dt_h = 1e-9;
Q = Darcy_Lateral_Topology(two,theta,0.2,soil,dt_h);
[Ka,ha] = soil_state(sp,theta(1));
[Kb,hb] = soil_state(sp2,theta(2));
Q_expected = (two.edge_face_width_m*0.2) * ((ha-hb)/1000) / ...
    (two.edge_face_distance_a_m/(Ka/1000) + ...
     two.edge_face_distance_b_m/(Kb/1000));
check('Darcy flux uses distance-weighted hydraulic resistance', ...
    abs(Q-Q_expected) < 1e-6*max(abs(Q_expected),eps), ...
    sprintf('got %.12g expected %.12g m3/h',Q,Q_expected));
check_error('Darcy rejects a zero timestep', ...
    @() Darcy_Lateral_Topology(two,theta,0.2,soil,0), ...
    'Darcy_Lateral_Topology:BadTimeStep');

C_bad = ones(2,1,1);
M_bad = ones(2,1,1);
check_error('solute transport rejects a nonfinite water flux', ...
    @() Solute_Lateral_Transport(two,NaN,C_bad,M_bad,1), ...
    'Solute_Lateral_Transport:BadFlux');

%% Liquid water cannot drain below residual or fill ice-occupied pores.
two_fast = two_cell_topology(0.005,0.005,1,[0.01;0.01]);
soil_fast = repmat(sp,2,1);
theta_floor = [Osat-1e-6;Ohy+1e-6];
[~,dtheta_floor] = Darcy_Lateral_Topology( ...
    two_fast,theta_floor,0.1,soil_fast,24);
theta_floor = theta_floor+dtheta_floor;
check('Darcy donor stays at or above residual water content', ...
    all(theta_floor >= Ohy-1e-12),sprintf('minimum theta %.12g',min(theta_floor)));

Oice = [0;0.20];
theta_ice = [Osat-1e-6;Osat-Oice(2)-1e-6];
[~,dtheta_ice] = Darcy_Lateral_Topology( ...
    two_fast,theta_ice,0.1,soil_fast,24,Oice);
theta_ice = theta_ice+dtheta_ice;
check('ice-filled pore space cannot accept liquid water', ...
    all(theta_ice <= Osat-Oice+1e-12), ...
    sprintf('theta %s, upper %s',mat2str(theta_ice),mat2str(Osat-Oice)));

%% A blocked frozen receiver must not force zero-flux substeps.
two_blocked = two_cell_topology(12.5,2.5,25,[625;250]);
delta_liquid = 0.02;
theta_blocked = [Ohy+delta_liquid;Ohy];
ice_blocked = [Osat-Ohy-delta_liquid;Osat-Ohy];
[Q_blocked,dtheta_blocked,n_blocked] = Darcy_Lateral_Topology( ...
    two_blocked,theta_blocked,0.04,soil_fast,1,ice_blocked);
check('ice-blocked zero flux completes in one substep', ...
    Q_blocked == 0 && all(dtheta_blocked == 0) && n_blocked == 1, ...
    sprintf('Q %.3e, max dtheta %.3e, substeps %d', ...
    Q_blocked,max(abs(dtheta_blocked)),n_blocked));

%% Empirical frozen-pore impedance is applied before edge resistance.
two_impedance = two_cell_topology(0.5,0.5,1,[1e12;1e12]);
mobile_range = Osat-Ohy;
theta_impedance = Ohy+mobile_range*[0.50;0.40];
ice_fraction = 0.25;
ice_impedance = ice_fraction/(1-ice_fraction)* ...
    (theta_impedance-Ohy);
Q_unfrozen = Darcy_Lateral_Topology( ...
    two_impedance,theta_impedance,0.1,soil_fast,1);
Q_frozen = Darcy_Lateral_Topology( ...
    two_impedance,theta_impedance,0.1,soil_fast,1,ice_impedance);
expected_ratio = 10^(-5*ice_fraction);
check('lateral frozen-pore impedance matches the vertical legacy Omega=5', ...
    abs(Q_frozen/Q_unfrozen-expected_ratio) < 1e-10, ...
    sprintf('ratio %.12g expected %.12g',Q_frozen/Q_unfrozen,expected_ratio));
%% Lateral flow uses the horizontal conductivity K*aR, as T&C's own does.
soil_aniso = soil_fast;
[soil_aniso.aR] = deal(0.5);
Q_aniso = Darcy_Lateral_Topology( ...
    two_impedance,theta_impedance,0.1,soil_aniso,1);
check('the anisotropy ratio scales lateral conductivity', ...
    abs(Q_aniso/Q_unfrozen-0.5) < 1e-10, ...
    sprintf('ratio %.12g expected 0.5',Q_aniso/Q_unfrozen));
check_error('Darcy rejects ice above mobile pore volume', ...
    @() Darcy_Lateral_Topology(two_impedance,theta_impedance,0.1, ...
    soil_fast,1,1.01*mobile_range*ones(2,1)), ...
    'Darcy_Lateral_Topology:BadIce');
check_error('Darcy rejects combined liquid and ice above pore volume', ...
    @() Darcy_Lateral_Topology(two_impedance,theta_impedance,0.1, ...
    soil_fast,1,0.60*mobile_range*ones(2,1)), ...
    'Darcy_Lateral_Topology:BadIce');

%% Soil conduction is directional, bounded, and energy conserving.
two_heat = two_cell_topology(0.5,0.5,1,[1;1]);
Q_zero = 0;
theta0 = [0.25;0.25];
T0 = [20;0];
lambda = [2;2];
cv = [2e6;2e6];
[T_cond,heat_diag] = Soil_Lateral_Heat_Transport( ...
    two_heat,Q_zero,theta0,theta0,T0,lambda,cv,1,1);
E0 = sum(cv.*two_heat.cell_area_m2.*T0);
E1 = sum(cv.*two_heat.cell_area_m2.*T_cond);
check('soil heat conducts from warm to cold', ...
    T_cond(1) < T0(1) && T_cond(2) > T0(2));
check('soil conduction does not overshoot neighbouring temperatures', ...
    all(T_cond >= min(T0)-1e-12) && all(T_cond <= max(T0)+1e-12));
check('closed lateral conduction conserves sensible energy', ...
    abs(E1-E0) < 1e-11*max(abs(E0),1),sprintf('residual %.6g J',E1-E0));
check('heat diagnostic closes its internal energy ledger', ...
    heat_diag.energy_residual_relative < 1e-12);

%% Moving liquid water carries sensible heat and changes heat capacity.
Q_adv = 0.1;
theta1 = [0.30;0.10];
theta2 = [0.20;0.20];
T1 = [20;0];
[T_adv,adv_diag] = Soil_Lateral_Heat_Transport( ...
    two_heat,Q_adv,theta1,theta2,T1,zeros(2,1),cv,1,1);
rho_cw = 1000*4186;
C_after = cv.*two_heat.cell_area_m2 + rho_cw*(theta2-theta1).*two_heat.cell_area_m2;
E_after = sum(C_after.*T_adv);
E_before = sum(cv.*two_heat.cell_area_m2.*T1);
check('advected water keeps the donor temperature and warms the receiver', ...
    abs(T_adv(1)-20) < 1e-12 && T_adv(2) > 0 && T_adv(2) < 20);
check('water-carried heat and changed capacity conserve sensible energy', ...
    abs(E_after-E_before) < 1e-11*max(abs(E_before),1));
check('advective heat is positive from warm cell a to cold cell b', ...
    adv_diag.advective_edge_J > 0);

%% Frozen lateral heat uses apparent capacity only to invert physical enthalpy.
p_thermal = thermal_parameters;
T_freeze = [0.2;-0.2];
Otot_freeze = [0.35;0.35];
lambda_freeze = zeros(2,1); capp_freeze = zeros(2,1);
cv_freeze = zeros(2,1); ice_freeze = zeros(2,1); liquid_freeze = zeros(2,1);
for c = 1:2
    [lambda_freeze(c),capp_freeze(c),ice_freeze(c), ...
        liquid_freeze(c),cv_freeze(c)] = Soil_Lateral_Thermal_State( ...
        T_freeze(c),Otot_freeze(c),1013,p_thermal,1);
end
check('freezing adds latent apparent heat capacity', ...
    capp_freeze(2) > cv_freeze(2) && ice_freeze(2) > 0 && ...
    abs(liquid_freeze(2)+ice_freeze(2)-Otot_freeze(2)) < 1e-12);
check_error('thermal state rejects water below residual content', ...
    @() Soil_Lateral_Thermal_State(-1, ...
    p_thermal.Soil_Param.Ohy-1e-6,1013,p_thermal,1), ...
    'Soil_Lateral_Thermal_State:BadWaterContent');
check_error('thermal state rejects water above porosity', ...
    @() Soil_Lateral_Thermal_State(-1, ...
    p_thermal.Soil_Param.Osat+1e-6,1013,p_thermal,1), ...
    'Soil_Lateral_Thermal_State:BadWaterContent');
enthalpy_step = 1e-4;
[~,capp_mid] = Soil_Lateral_Thermal_State( ...
    -0.2,0.35,1013,p_thermal,1);
H_minus = physical_soil_enthalpy( ...
    -0.2-enthalpy_step,0.35,1013,p_thermal,1);
H_plus = physical_soil_enthalpy( ...
    -0.2+enthalpy_step,0.35,1013,p_thermal,1);
capp_fd = (H_plus-H_minus)/(2*enthalpy_step);
check('apparent capacity is the derivative of physical enthalpy', ...
    abs(capp_mid-capp_fd) < 1e-5*max(1,abs(capp_fd)), ...
    sprintf('capp %.9g, finite difference %.9g',capp_mid,capp_fd));
phase_freeze = phase_context(ice_freeze,p_thermal,2);
[T_latent,latent_diag,liquid_after,ice_after] = ...
    Soil_Lateral_Heat_Transport(two_heat,0,liquid_freeze, ...
    liquid_freeze,T_freeze,lambda_freeze,capp_freeze,1,1,phase_freeze);
[T_sensible,~] = Soil_Lateral_Heat_Transport(two_heat,0, ...
    liquid_freeze,liquid_freeze,T_freeze,lambda_freeze,cv_freeze,1,1);
check('latent capacity slows warming of frozen soil', ...
    T_latent(2)-T_freeze(2) < T_sensible(2)-T_freeze(2));
check('post-transport liquid and ice remain in phase equilibrium', ...
    ice_after(2) <= ice_freeze(2) && ...
    max(abs(liquid_after+ice_after-Otot_freeze)) < 1e-12);
check('frozen conduction closes the scaled physical-energy ledger', ...
    latent_diag.energy_residual_relative < 1e-10, ...
    sprintf('scaled residual %.6g',latent_diag.energy_residual_relative));

%% A stiff frozen edge remains bounded under the apparent-capacity limiter.
two_stiff = two_cell_topology(1e-4,1e-4,1,[1;1]);
T_stiff = [2;-2]; Otot_stiff = [0.35;0.35];
lambda_stiff = zeros(2,1); capp_stiff = zeros(2,1);
ice_stiff = zeros(2,1); liquid_stiff = zeros(2,1);
for c = 1:2
    [lambda_stiff(c),capp_stiff(c),ice_stiff(c),liquid_stiff(c)] = ...
        Soil_Lateral_Thermal_State( ...
        T_stiff(c),Otot_stiff(c),1013,p_thermal,1);
end
[T_stiff_after,stiff_diag,liquid_stiff_after,ice_stiff_after] = ...
    Soil_Lateral_Heat_Transport(two_stiff,0,liquid_stiff, ...
    liquid_stiff,T_stiff,lambda_stiff,capp_stiff,1,1, ...
    phase_context(ice_stiff,p_thermal,2));
check('high-conductance frozen heat remains within neighbour bounds', ...
    stiff_diag.sub_steps > 1 && all(T_stiff_after >= min(T_stiff)-1e-10) && ...
    all(T_stiff_after <= max(T_stiff)+1e-10));
check('stiff frozen heat conserves water and physical enthalpy', ...
    max(abs(liquid_stiff_after+ice_stiff_after-Otot_stiff)) < 1e-12 && ...
    stiff_diag.energy_residual_relative < 1e-10);

%% Apparent capacity is dH/dT, not stored enthalpy divided by temperature.
T_phase = [5;-2];
Otot_phase = [0.35;p_thermal.Soil_Param.Ohy];
lambda_phase = zeros(2,1); capp_phase = zeros(2,1);
ice_phase = zeros(2,1); liquid_phase = zeros(2,1);
H_phase0 = zeros(2,1);
for c = 1:2
    [lambda_phase(c),capp_phase(c),ice_phase(c),liquid_phase(c)] = ...
        Soil_Lateral_Thermal_State(T_phase(c),Otot_phase(c),1013,p_thermal,1);
    H_phase0(c) = physical_soil_enthalpy( ...
        T_phase(c),Otot_phase(c),1013,p_thermal,1);
end
theta_phase_after = liquid_phase+[-0.08;0.08];
[T_phase_after,phase_diag,liquid_phase_after,ice_phase_after] = ...
    Soil_Lateral_Heat_Transport(two_heat,0.08,liquid_phase, ...
    theta_phase_after,T_phase,lambda_phase,capp_phase,1,1, ...
    phase_context(ice_phase,p_thermal,2));
Otot_phase_after = theta_phase_after+ice_phase;
H_phase1 = zeros(2,1);
liquid_phase_check = zeros(2,1);
ice_phase_check = zeros(2,1);
for c = 1:2
    [~,~,ice_phase_check(c),liquid_phase_check(c)] = ...
        Soil_Lateral_Thermal_State(T_phase_after(c), ...
        Otot_phase_after(c),1013,p_thermal,1);
    H_phase1(c) = physical_soil_enthalpy( ...
        T_phase_after(c),Otot_phase_after(c),1013,p_thermal,1);
end
phase_residual = sum(H_phase1)-sum(H_phase0);
check('freezing incoming water conserves physical enthalpy', ...
    ice_phase_after(2) > 0 && ...
    abs(phase_residual) < 1e-10*max(1,sum(abs(H_phase0))), ...
    sprintf('residual %.6g J',phase_residual));
check('cold-receiver phase update returns all transported water', ...
    max(abs(liquid_phase_after+ice_phase_after-Otot_phase_after)) < 1e-12);
check('returned cold-receiver phases match thermal equilibrium', ...
    max(abs(liquid_phase_after-liquid_phase_check)) < 1e-12 && ...
    max(abs(ice_phase_after-ice_phase_check)) < 1e-12);
check('cold-receiver diagnostic uses a scaled physical residual', ...
    phase_diag.energy_residual_relative < 1e-10, ...
    sprintf('scaled residual %.6g',phase_diag.energy_residual_relative));
%% T&C never hands over a state on the local freezing curve: Richards flow
%% and Cryosuction_stabilizer move liquid water between layers after the
%% phase split is fixed. The lateral step must therefore ACCEPT an
%% off-curve state, carry the departure, and still close its ledgers.
phase_off = phase_context(ice_phase,p_thermal,2);
ice_bump = 1e-4;
phase_off.Oice_cell(1) = phase_off.Oice_cell(1)+ice_bump;
Otot_off = theta_phase_after+phase_off.Oice_cell;
[T_off,diag_off,liquid_off,ice_off] = Soil_Lateral_Heat_Transport( ...
    two_heat,0.08,liquid_phase,theta_phase_after,T_phase, ...
    lambda_phase,capp_phase,1,1,phase_off);
check('lateral heat accepts an off-equilibrium initial state', ...
    all(isfinite(T_off)) && all(ice_off >= 0) && all(liquid_off >= 0));
check('off-equilibrium water still closes', ...
    max(abs(liquid_off+ice_off-Otot_off)) < 1e-12);
check('off-equilibrium energy ledger still closes', ...
    diag_off.energy_residual_relative < 1e-10, ...
    sprintf('scaled residual %.6g',diag_off.energy_residual_relative));
ice_eq_off = zeros(2,1);
for c = 1:2
    [~,~,ice_eq_off(c),~] = Soil_Lateral_Thermal_State( ...
        T_off(c),Otot_off(c),1013,p_thermal,1);
end
check('the departure from equilibrium is carried, not erased', ...
    abs((ice_off(1)-ice_eq_off(1))-ice_bump) < 1e-12);

%% One cell is an exact no-op, including the new thermal path.
one = Build_Single_Cell_Topology(1);
[sp_one,~,~] = soil_parameters(2);
thermal = struct('T_cell',[7 8],'lambda_cell',[1.2 1.3], ...
    'capp_cell',[1.8e6 1.9e6],'Oice_cell',[0 0]);
thermal.parameter_cell = {p_thermal};
thermal.Pre = 1013;
thermal.OPT_FR_SOIL = 1;
[theta_one,dM_one,~,T_one,ice_one] = Lateral_Routing_Step( ...
    one,[0.2 0.3],zeros(1,2,1),sp_one,[0.4 0.6],1,thermal);
check('single-cell lateral water, mass, and heat paths are exact no-ops', ...
    isequal(theta_one,[0.2 0.3]) && isequal(dM_one,zeros(1,2,1)) && ...
    isequal(T_one,[7 8]) && isequal(ice_one,[0 0]));

%% T&C dissolved pools occupy only the biogeochemical active depth.
area_bio = [20;80];
theta_bio = [0.20 0.30 0.40; 0.25 0.30 0.35];
Z_bio = [0.10 0.20 0.30];
active_bio = [1 0.5 0; 1 0.5 0];
mobile_bio = [2 4; 3 6];
C_bio = Biogeo_Lateral_Concentration( ...
    mobile_bio,area_bio,theta_bio,Z_bio,active_bio);
water_bio = theta_bio.*(area_bio*Z_bio).*active_bio;
mass_bio = reshape(sum(C_bio.*reshape(water_bio,2,3,1),2),2,2);
check('active-depth concentration reconstructs every mobile T&C pool', ...
    max(abs(mass_bio-mobile_bio.*area_bio),[],'all') < 1e-12);
check('T&C solute concentration is zero below the active depth', ...
    all(C_bio(:,3,:)==0,'all'));
Q_deep = zeros(1,3); Q_deep(3) = 0.1;
[Q_deep,~] = Restrict_Lateral_Flux_To_Active_Depth( ...
    two_cell_topology(0.5,0.5,1,area_bio),Q_deep,active_bio);
[dM_deep,~] = Solute_Lateral_Transport(two_cell_topology(0.5,0.5,1,area_bio), ...
    Q_deep,C_bio,C_bio.*reshape(water_bio,2,3,1),1);
check('deep lateral water cannot advect active-layer T&C pools', ...
    all(dM_deep(:)==0));

%% ERW-only aqueous stores and soil-gas CO2 move conservatively.
two_erw = two_cell_topology(0.5,0.5,1,[1;1]);
E = synthetic_erw_state;
Q_erw = 1e-3;
theta_erw = [0.27;0.27];
active = ones(2,1);
tot0 = erw_totals(E,two_erw.cell_area_m2);
[E_aq,erw_diag] = ERW_Lateral_Transport( ...
    E,two_erw,Q_erw,theta_erw,0.3,active,1);
tot1 = erw_totals(E_aq,two_erw.cell_area_m2);
%% Soil air pushed out or drawn in by lateral water is booked once, by the
%% ADV term of the next ERW step, which sees the same saturation change;
%% the lateral step itself closes inorganic carbon.
check('ERW Na, Al, residual anion, DIC and IC lateral transfer conserves mass', ...
    max(abs(tot1-tot0)) < 1e-9*max(max(abs(tot0)),1) && ...
    ~isfield(erw_diag,'CO2_atmosphere_umol_cell') && ...
    ~isfield(erw_diag,'IC_atmosphere_residual_umol'));
erw_source = fileread(fullfile(fileparts(fileparts(mfilename('fullpath'))),'ERW_function.m'));
check('the ERW step books the soil-air displacement', ...
    contains(erw_source,'ADV = n*Zrw*1000*(s-s_pre)*CO2_air;') && ...
    contains(erw_source,'ADV = n*Zrw*1000*(s-s_pre)*CO2_atm;'));
check('ERW aqueous species move with water from cell a to cell b', ...
    E_aq{1}.bg.Na_tot(1) < E{1}.bg.Na_tot(1) && ...
    E_aq{2}.bg.Na_tot(1) > E{2}.bg.Na_tot(1) && ...
    all(erw_diag.aqueous_edge_umol > 0));
check('ERW aqueous concentrations reflect receiver mixing', ...
    E_aq{2}.bg.Na(1) > E{2}.bg.Na(1) && ...
    E_aq{2}.bg.DIC(1) > E{2}.bg.DIC(1));
al_mobile = E{1}.bg.Al(1)+E{1}.bg.AlOH(1)+E{1}.bg.AlOH2(1)+ ...
    E{1}.bg.AlOH3(1)+E{1}.bg.AlOH4(1);
check('all dissolved aluminum species advect', ...
    abs(erw_diag.aqueous_edge_umol(2)-al_mobile*Q_erw) < 1e-9);

E_uniform = synthetic_erw_state;
aq_fields = {'Na','Al','AlOH','AlOH2','AlOH3','AlOH4','An','DIC'};
for k = 1:numel(aq_fields)
    E_uniform{2}.bg.(aq_fields{k})(1) = E_uniform{1}.bg.(aq_fields{k})(1);
end
theta_uniform = repmat(E_uniform{1}.const.n*E_uniform{1}.current_s,2,1);
s0 = E_uniform{1}.current_s;
[E_uniform_after,uniform_diag] = ERW_Lateral_Transport( ...
    E_uniform,two_erw,Q_erw,theta_uniform,0.3,active,1);
for k = 1:numel(aq_fields)
    check(['uniform ' aq_fields{k} ' survives advection'], ...
        abs(E_uniform_after{1}.bg.(aq_fields{k})(1)- ...
        E_uniform{1}.bg.(aq_fields{k})(1)) < 1e-12 && ...
        abs(E_uniform_after{2}.bg.(aq_fields{k})(1)- ...
        E_uniform{2}.bg.(aq_fields{k})(1)) < 1e-12);
end
ds = Q_erw/(E_uniform{1}.const.n*E_uniform{1}.const.Zrw);
check('ERW saturation follows lateral water transfer', ...
    abs(E_uniform_after{1}.current_s-(s0-ds)) < 1e-12 && ...
    abs(E_uniform_after{2}.current_s-(s0+ds)) < 1e-12);
check('uniform ERW transfer remains conservative', ...
    all(abs(uniform_diag.mass_residual_umol) < 1e-9));

E_gas = synthetic_erw_state;
E_gas{1}.bg.DIC(:) = 0; E_gas{2}.bg.DIC(:) = 0;
gas0 = erw_totals(E_gas,two_erw.cell_area_m2);
[E_gas,gas_diag] = ERW_Lateral_Transport( ...
    E_gas,two_erw,0,theta_erw,0.3,active,1);
gas1 = erw_totals(E_gas,two_erw.cell_area_m2);
check('soil-gas CO2 diffuses from high to low concentration', ...
    gas_diag.CO2_edge_umol > 0 && ...
    E_gas{1}.bg.IC_tot(1) < E{1}.bg.IC_tot(1));
check('soil-gas CO2 concentrations relax after diffusion', ...
    E_gas{1}.bg.CO2_air(1) < E{1}.bg.CO2_air(1) && ...
    E_gas{2}.bg.CO2_air(1) > E{2}.bg.CO2_air(1));
check('closed lateral soil-gas diffusion conserves inorganic carbon', ...
    abs(gas1(4)-gas0(4)) < 1e-9*max(abs(gas0(4)),1));

star = star_topology;
E_star = synthetic_erw_state(4);
for c = 1:4
    E_star{c}.bg.D = 50;
    E_star{c}.bg.CO2_air = 2000;
end
E_star{1}.bg.CO2_air = 400;
star_min = min(cellfun(@(x) x.bg.CO2_air(1),E_star));
star_max = max(cellfun(@(x) x.bg.CO2_air(1),E_star));
star0 = erw_totals(E_star,star.cell_area_m2);
[E_star,star_diag] = ERW_Lateral_Transport( ...
    E_star,star,zeros(3,1),0.27*ones(4,1),0.3,ones(4,1),1);
star1 = erw_totals(E_star,star.cell_area_m2);
star_after = cellfun(@(x) x.bg.CO2_air(1),E_star);
check('multi-edge gas diffusion preserves concentration bounds', ...
    all(star_after >= star_min-1e-12) && all(star_after <= star_max+1e-12));
check('multi-edge gas diffusion remains conservative', ...
    abs(star1(4)-star0(4)) < 1e-9*max(abs(star0(4)),1) && ...
    abs(star_diag.mass_residual_umol(4)) < 1e-9);

E_sat = synthetic_erw_state;
E_sat{2}.current_s = 1;
sat0 = erw_totals(E_sat,two_erw.cell_area_m2);
[E_sat,sat_diag] = ERW_Lateral_Transport( ...
    E_sat,two_erw,0,[0.27;0.45],0.3,active,1);
sat1 = erw_totals(E_sat,two_erw.cell_area_m2);
check('a saturated cell closes the lateral soil-gas face', ...
    sat_diag.CO2_edge_umol == 0 && isequal(sat1,sat0));

%% Water-volume changes leave the soil-gas concentration of each cell as it
%% was; the displaced or added air volume is the ERW step's ADV to book.
E_wet = synthetic_erw_state;
E_wet{1}.bg.D(:) = 0; E_wet{2}.bg.D(:) = 0;
pore_volume = E_wet{2}.const.n*E_wet{2}.const.Zrw;
Q_fill = (1-E_wet{2}.current_s)*pore_volume;
wet0 = erw_totals(E_wet,two_erw.cell_area_m2);
[E_wet,wet_diag] = ERW_Lateral_Transport( ...
    E_wet,two_erw,Q_fill,theta_erw,0.3,active,1);
wet1 = erw_totals(E_wet,two_erw.cell_area_m2);
check('lateral wetting to saturation leaves no soil-air volume', ...
    E_wet{2}.current_s == 1 && E_wet{2}.bg.CO2_air(1) == 0);
check('lateral wetting and drying close inorganic carbon', ...
    abs(wet1(4)-wet0(4)) < 1e-10*max(1,abs(wet0(4))) && ...
    abs(wet_diag.mass_residual_umol(4)) < 1e-8);
check('lateral drying keeps the donor soil-gas concentration', ...
    abs(E_wet{1}.bg.CO2_air(1)-E{1}.bg.CO2_air(1)) < 1e-12*E{1}.bg.CO2_air(1));

fprintf('  ---------------------------------------\n');
fprintf('  %d passed, %d failed\n\n',n_pass,n_fail);
if n_fail > 0
    error('test_lateral_energy_erw: %d check(s) failed',n_fail);
end

end

function [sp,Osat,Ohy] = soil_parameters(ms)
[Osat,Lslope,Pe,Ks,O33] = Soil_parameters(0.15,0.30,0.03);
Ohy = 0.02;
sp = struct('SPAR',2,'Ks',Ks*ones(1,ms),'Osat',Osat*ones(1,ms), ...
    'Ohy',Ohy*ones(1,ms),'L',Lslope*ones(1,ms),'Pe',Pe*ones(1,ms), ...
    'O33',O33*ones(1,ms),'alpVG',1e-3*ones(1,ms),'nVG',2*ones(1,ms), ...
    'lVG',0.5*ones(1,ms),'Ks_mac',zeros(1,ms),'Omac',zeros(1,ms), ...
    'alpVGM',1e-3*ones(1,ms),'nVGM',2*ones(1,ms),'lVGM',0.5*ones(1,ms), ...
    'Phy1',10000,'s_SVG',zeros(1,ms),'bVG',zeros(1,ms),'aR',1);
end

function [K,h] = soil_state(sp,theta)
[K,p] = Conductivity_Suction(sp.SPAR,sp.Ks,sp.Osat,sp.Ohy,sp.L,sp.Pe, ...
    sp.O33,sp.alpVG,sp.nVG,sp.lVG,sp.Ks_mac,sp.Omac,sp.alpVGM, ...
    sp.nVGM,sp.lVGM,sp.Phy1,sp.s_SVG,sp.bVG,theta);
h = -p;
end

function topo = two_cell_topology(da,db,width,area)
topo.N_cells = 2;
topo.N_edges = 1;
topo.edge_a = 1;
topo.edge_b = 2;
topo.edge_face_width_m = width;
topo.edge_face_distance_a_m = da;
topo.edge_face_distance_b_m = db;
topo.edge_face_distance_m = da+db;
topo.edge_active = 1;
topo.cell_area_m2 = area;
topo.cell_region = {'main';'buffer'};
end

function p = thermal_parameters
[Osat,L,Pe,Ks,O33,rsd,lan_dry,lan_s,cv_s] = ...
    Soil_parameters(0.15,0.30,0.03);
pvg = 3+2/L;
mvg = 2/(pvg-1);
nVG = 1/(1-mvg);
alpVG = (((-101.9368*Pe)*(2*pvg*(pvg-1))/(pvg+3))* ...
    ((55.6+7.4*pvg+pvg^2)/(147.8+8.1*pvg+0.092*pvg^2)))^-1;
[~,~,~,Ohy] = Soil_parametersII(1,Osat,L,Pe,Ks,O33,nVG,alpVG, ...
    0.2,1,1,10000);
p.SPAR = 2;
p.Phy = 10000;
p.Soil_Param = struct('rsd',rsd,'lan_dry',lan_dry,'lan_s',lan_s, ...
    'cv_s',cv_s,'L',L,'Pe',Pe,'O33',O33,'alpVG',alpVG,'nVG',nVG, ...
    's_SVG',0,'bVG',0,'Osat',Osat,'Ohy',Ohy);
end

function H = physical_soil_enthalpy(T,Otot,Pre,p,opt_freeze)
[~,~,Oice,~,cv] = Soil_Lateral_Thermal_State(T,Otot,Pre,p,opt_freeze);
H = cv.*T-916.2*333700*Oice;
end

function phase = phase_context(Oice,p,N_cells)
phase = struct('Oice_cell',Oice,'Pre',1013,'OPT_FR_SOIL',1);
phase.parameter_cell = repmat({p},N_cells,1);
end

function topo = star_topology
topo.N_cells = 4;
topo.N_edges = 3;
topo.edge_a = [2;3;4];
topo.edge_b = [1;1;1];
topo.edge_face_width_m = ones(3,1);
topo.edge_face_distance_a_m = 0.5*ones(3,1);
topo.edge_face_distance_b_m = 0.5*ones(3,1);
topo.edge_face_distance_m = ones(3,1);
topo.edge_active = ones(3,1);
topo.cell_area_m2 = ones(4,1);
topo.cell_region = {'main';'buffer';'buffer';'buffer'};
end

function E = synthetic_erw_state(N)
if nargin < 1; N = 2; end
E = cell(N,1);
for c = 1:N
    E{c}.init = true;
    E{c}.dm_end = 1;
    E{c}.current_s = 0.60;
    E{c}.const.n = 0.45;
    E{c}.const.Zrw = 0.30;
    E{c}.bg.Na = 100*(c==1) + 10*(c==2);
    E{c}.bg.Al = 20*(c==1) + 2*(c==2);
    E{c}.bg.AlOH = 4*(c==1) + 0.4*(c~=1);
    E{c}.bg.AlOH2 = 3*(c==1) + 0.3*(c~=1);
    E{c}.bg.AlOH3 = 2*(c==1) + 0.2*(c~=1);
    E{c}.bg.AlOH4 = 5*(c==1) + 0.5*(c==2);
    E{c}.bg.An = 50*(c==1) + 5*(c==2);
    E{c}.bg.DIC = 500*(c==1) + 50*(c==2);
    E{c}.bg.CO2_air = 2000*(c==1) + 400*(c==2);
    E{c}.bg.D = 0.50;
    E{c}.bg.Na_tot = 1e5;
    E{c}.bg.Al_tot = 5e4;
    E{c}.bg.An_tot = 5e4;
    E{c}.bg.IC_tot = 5e5;
    E{c}.bg.Alk_tot = 5e4;
end
end

function totals = erw_totals(E,area)
fields = {'Na_tot','Al_tot','An_tot','IC_tot'};
to_umol = [1 1e-3 1 1];
totals = zeros(1,numel(fields));
for c = 1:numel(E)
    for s = 1:numel(fields)
        totals(s) = totals(s) + E{c}.bg.(fields{s})(E{c}.dm_end)* ...
            area(c)*to_umol(s);
    end
end
end
