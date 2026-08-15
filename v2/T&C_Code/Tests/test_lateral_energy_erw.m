function test_lateral_energy_erw()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Unequal-cell resistance, physical water bounds, heat, and ERW transport.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n=== test_lateral_energy_erw ===\n');
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

%% Hansson frozen-pore impedance is applied before edge resistance.
two_impedance = two_cell_topology(0.5,0.5,1,[1e12;1e12]);
mobile_range = Osat-Ohy;
theta_impedance = Ohy+mobile_range*[0.50;0.40];
ice_fraction = 0.25;
ice_impedance = ice_fraction*mobile_range*ones(2,1);
Q_unfrozen = Darcy_Lateral_Topology( ...
    two_impedance,theta_impedance,0.1,soil_fast,1);
Q_frozen = Darcy_Lateral_Topology( ...
    two_impedance,theta_impedance,0.1,soil_fast,1,ice_impedance);
expected_ratio = 10^(-7*ice_fraction);
check('frozen-pore impedance reduces liquid conductivity', ...
    abs(Q_frozen/Q_unfrozen-expected_ratio) < 1e-10, ...
    sprintf('ratio %.12g expected %.12g',Q_frozen/Q_unfrozen,expected_ratio));

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
    abs(heat_diag.energy_residual_J) < 1e-9);

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

%% One cell is an exact no-op, including the new thermal path.
one = Build_Single_Cell_Topology(1);
[sp_one,~,~] = soil_parameters(2);
thermal = struct('T_cell',[7 8],'lambda_cell',[1.2 1.3], ...
    'cv_cell',[1.8e6 1.9e6],'Oice_cell',[0 0]);
[theta_one,dM_one,~,T_one] = Lateral_Routing_Step( ...
    one,[0.2 0.3],zeros(1,2,1),sp_one,[0.4 0.6],1,thermal);
check('single-cell lateral water, mass, and heat paths are exact no-ops', ...
    isequal(theta_one,[0.2 0.3]) && isequal(dM_one,zeros(1,2,1)) && ...
    isequal(T_one,[7 8]));

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
theta_erw = [0.30;0.30];
active = ones(2,1);
tot0 = erw_totals(E,two_erw.cell_area_m2);
[E_aq,erw_diag] = ERW_Lateral_Transport( ...
    E,two_erw,Q_erw,theta_erw,0.3,active,1);
tot1 = erw_totals(E_aq,two_erw.cell_area_m2);
check('ERW Na, Al, residual anion, and DIC advection conserves mass', ...
    max(abs(tot1-tot0)) < 1e-9*max(max(abs(tot0)),1));
check('ERW aqueous species move with water from cell a to cell b', ...
    E_aq{1}.bg.Na_tot(1) < E{1}.bg.Na_tot(1) && ...
    E_aq{2}.bg.Na_tot(1) > E{2}.bg.Na_tot(1) && ...
    all(erw_diag.aqueous_edge_umol > 0));
check('ERW aqueous concentrations reflect receiver mixing', ...
    E_aq{2}.bg.Na(1) > E{2}.bg.Na(1) && ...
    E_aq{2}.bg.DIC(1) > E{2}.bg.DIC(1));

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

E_sat = synthetic_erw_state;
E_sat{2}.current_s = 1;
sat0 = erw_totals(E_sat,two_erw.cell_area_m2);
[E_sat,sat_diag] = ERW_Lateral_Transport( ...
    E_sat,two_erw,0,theta_erw,0.3,active,1);
sat1 = erw_totals(E_sat,two_erw.cell_area_m2);
check('a saturated cell closes the lateral soil-gas face', ...
    sat_diag.CO2_edge_umol == 0 && isequal(sat1,sat0));

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
    'Phy1',10000,'s_SVG',zeros(1,ms),'bVG',zeros(1,ms));
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

function E = synthetic_erw_state
E = cell(2,1);
for c = 1:2
    E{c}.init = true;
    E{c}.dm_end = 1;
    E{c}.current_s = 0.60;
    E{c}.const.n = 0.45;
    E{c}.const.Zrw = 0.30;
    E{c}.bg.Na = 100*(c==1) + 10*(c==2);
    E{c}.bg.Al = 20*(c==1) + 2*(c==2);
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
totals = zeros(1,numel(fields));
for c = 1:numel(E)
    for s = 1:numel(fields)
        totals(s) = totals(s) + E{c}.bg.(fields{s})(E{c}.dm_end)*area(c);
    end
end
end
