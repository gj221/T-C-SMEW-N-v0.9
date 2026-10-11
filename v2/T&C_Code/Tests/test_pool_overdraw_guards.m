function test_pool_overdraw_guards
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Guards against the pool-overdraw failures that stopped the 2-D runs.
%% Covers lateral transport, uptake/leakage, reactions, decomposition,
%% stoichiometric transfers and exact-zero states.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

here = fileparts(mfilename('fullpath'));
code = fileparts(here);
addpath(code,fullfile(code,'Inputs'),fullfile(code,'Parameters'));
fprintf('\n=== test_pool_overdraw_guards ===\n');
npass = 0; nfail = 0;

%% --- 1/2. ERW zero-flux invariant -------------------------------------
topo = Build_Cell_Topology(60,5,2,2);
ms = 4;
Zlayer_m = [0.10 0.20 0.40 0.80];
N = topo.N_cells;
theta = repmat([0.30 0.30 0.28 0.26],N,1);
active = ones(N,ms);
E = erw_columns(N);

Q_zero = zeros(topo.N_edges,ms);
[E_out,diag] = ERW_Lateral_Transport(E,topo,Q_zero,theta,Zlayer_m,active,1,18);

max_drift = 0;
fields = {'Na','Al','AlOH','AlOH2','AlOH3','AlOH4','An','DIC','CO2_air', ...
          'Na_tot','Al_tot','An_tot','IC_tot','Alk_tot'};
for c = 1:N
    for k = 1:numel(fields)
        a = E{c}.bg.(fields{k})(E{c}.dm_end);
        b = E_out{c}.bg.(fields{k})(E_out{c}.dm_end);
        max_drift = max(max_drift,abs(b-a)/max(1,abs(a)));
    end
end
[npass,nfail] = check(max_drift < 1e-12, ...
    sprintf('zero water flux leaves every ERW store unchanged (drift %.3g)', ...
    max_drift),npass,nfail);
[npass,nfail] = check(all(abs(diag.mass_residual_umol) < 1e-6), ...
    'zero-flux mass residual is closed',npass,nfail);

%% Concentration must invert against the ERW pore volume, not T&C water.
c1 = E{1};
V_erw = c1.const.n*c1.current_s*c1.const.Zrw*topo.cell_area_m2(1);
V_tc = sum(theta(1,:).*Zlayer_m)*topo.cell_area_m2(1);
[npass,nfail] = check(abs(V_erw-V_tc)/V_tc > 0.1, ...
    'the two candidate water volumes really do differ',npass,nfail);
implied = c1.bg.Na(1)*1000*V_erw;          %% umol from ERW pore water
[npass,nfail] = check(implied <= c1.bg.Na_tot(1)*topo.cell_area_m2(1)*(1+1e-9), ...
    'ERW pore-water Na does not exceed its own total',npass,nfail);

%% --- 3. Immobilization deficit path -----------------------------------
F = load(fullfile(here,'Fixtures','no3_immobilization_20190906.mat'),'fixture');
F = F.fixture;
F.args{end+1} = 0.8;
%% Close the system: zero every external N source and sink so the only
%% remaining N exchanges are the internal transfers plus gaseous loss.
%% Arg order: 12 FertN, 13 DepN, 14 BfixN, 25 NH4_Uptake, 26 NO3_Uptake,
%% 33 LEAK_NH4, 34 LEAK_NO3, 37 LEAK_DON.
args = F.args;
for k = [12 13 14 25 26 33 34 37]
    args{k} = 0;
end
%% Same forcing, two SOM states: one where SOM N easily covers the
%% organic draw (deficit path inactive) and one where it cannot (active).
%% Litter inputs and gaseous losses are identical in form, so the net N
%% change must agree between them. Any disagreement is N created or
%% destroyed by the withholding logic.
args{2}([31 32]) = 0;      %% no mineral N at all in either case

rich = args;  rich{2}(26) = 1e3;    %% SOM N far above any demand
poor = args;  poor{2}(26) = 1e-9;   %% SOM N essentially absent

[net_rich,P_rich] = closed_N_change(rich);
[net_poor,P_poor] = closed_N_change(poor);

[npass,nfail] = check(all(isfinite(P_poor(23:34))), ...
    'depleted mineral and SOM nitrogen stays finite',npass,nfail);
[npass,nfail] = check(P_poor(31) >= 0 && P_poor(32) >= 0, ...
    sprintf('mineral N stays nonnegative (NH4 %.3g, NO3 %.3g)', ...
    P_poor(31),P_poor(32)),npass,nfail);
[npass,nfail] = check(P_poor(26) >= -1e-12, ...
    sprintf('SOM nitrogen is not overdrawn (%.6g)',P_poor(26)),npass,nfail);
%% The conservation check above is only meaningful if the two cases
%% actually diverge, i.e. the withholding branch really ran.
microbe_rich = sum(P_rich(27:30));
microbe_poor = sum(P_poor(27:30));
[npass,nfail] = check(microbe_poor < microbe_rich - 1e-12, ...
    sprintf('depleted SOM withholds microbial N (%.6g < %.6g)', ...
    microbe_poor,microbe_rich),npass,nfail);
[npass,nfail] = check(abs(net_poor-net_rich) < 1e-9*max(1,abs(net_rich)), ...
    sprintf('withheld immobilization conserves N (%.3g vs %.3g)', ...
    net_poor,net_rich),npass,nfail);
[~,~,~,~,~,~,~,Min_N_poor] = BIOGEOCHEMISTRY_DYNAMIC3(poor{:});
[npass,nfail] = check(abs(Min_N_poor) < 1e-12, ...
    sprintf('reported mineral-N transfer is realized, not requested (%.6g)', ...
    Min_N_poor),npass,nfail);

som = args;
som{2}([9 11 22 27:30 33]) = 0;
som{2}(10) = 0.025;
som{2}(26) = 0.005;
som{50}.lambda_N = 1;
[dN,~,~,~,~,~,~,~,~,~,~,~,~,~,~] = ...
    BIOGEOCHEMISTRY_DYNAMIC3(som{:});
som_N_new = som{2}(26)+dN(26);
[npass,nfail] = check(som_N_new >= -1e-12, ...
    sprintf('other daily sinks cannot overdraw SOM N (%.6g)', ...
    som_N_new),npass,nfail);

%% A dissolved pool depleted by leaching cannot also feed fixation.
K_args = args;
K_args{2}(52:55) = [1e6 0 0 0];
K_args{17} = 0; K_args{18} = 0;
K_args{28} = 0; K_args{36} = K_args{2}(52); K_args{43} = 0;
[dK,~,~,~,~,~,~,~,~,~,~,~,~,~,~] = ...
    BIOGEOCHEMISTRY_DYNAMIC3(K_args{:});
K_new = K_args{2}+dK;
[npass,nfail] = check(K_new(52) >= -1e-12, ...
    sprintf('leaching plus fixation does not overdraw dissolved K (%.6g)', ...
    K_new(52)),npass,nfail);

%% Mineral N reactions share what remains after uptake and leaching.
no3 = F.args;
no3{2}(1:34) = 0; no3{2}(32) = 1;
no3{5}(:) = 0; no3{52}(:) = 0;
for k = [12:16 25:27 32:33 35:38]; no3{k} = 0; end
no3{34} = 1;
dNO3 = BIOGEOCHEMISTRY_DYNAMIC3(no3{:});
NO3_new = no3{2}+dNO3;
[npass,nfail] = check(NO3_new(32) >= -1e-12, ...
    sprintf('NO3 leakage plus denitrification stays within the pool (%.6g)', ...
    NO3_new(32)),npass,nfail);

nh4 = no3;
nh4{2}(32) = 0; nh4{2}(31) = 1;
nh4{34} = 0; nh4{25} = 0.95; nh4{33} = 0.05;
dNH4 = BIOGEOCHEMISTRY_DYNAMIC3(nh4{:});
NH4_new = nh4{2}+dNH4;
[npass,nfail] = check(NH4_new(31) >= -1e-12, ...
    sprintf('NH4 uptake, leakage and reactions stay within the pool (%.6g)', ...
    NH4_new(31)),npass,nfail);

%% Microbial P demand cannot borrow phosphorus from an empty SOM pool.
p = F.args;
p{2}(35:47) = 0; p{5}(13:15) = 0; p{52}(5) = 0;
for k = [15 16 27 35 38 42]; p{k} = 0; end
dP = BIOGEOCHEMISTRY_DYNAMIC3(p{:});
P_new = p{2}+dP;
[npass,nfail] = check(all(P_new(35:47) >= 0), ...
    sprintf('empty organic and mineral P pools remain nonnegative (min %.6g)', ...
    min(P_new(35:47))),npass,nfail);
[npass,nfail] = check(abs(sum(dP(35:47))) < 1e-10, ...
    sprintf('withheld P immobilization conserves phosphorus (%.6g)', ...
    sum(dP(35:47))),npass,nfail);
P_handoff_ok = true;
try
    Biogeo_Lateral_Pools(P_new);
catch
    P_handoff_ok = false;
end
[npass,nfail] = check(P_handoff_ok, ...
    'exact-zero P state is accepted by the lateral pool handoff', ...
    npass,nfail);

p_leak = F.args;
p_leak{2}([38 39:42 47]) = 0; p_leak{2}(43) = 0.05;
p_leak{5}(13:15) = 0; p_leak{52}(5) = 0;
for k = [15 16 27 38 42]; p_leak{k} = 0; end
p_leak{35} = 1e-4;
p_leak_ok = true; P_leak_new = nan(1,79); P_leak_transfer = nan;
try
    [dP_leak,~,~,~,~,~,~,~,P_leak_transfer] = ...
        BIOGEOCHEMISTRY_DYNAMIC3(p_leak{:});
    P_leak_new = p_leak{2}+dP_leak;
catch
    p_leak_ok = false;
end
[npass,nfail] = check(p_leak_ok && all(P_leak_new(35:47) >= -1e-12), ...
    'P leakage and immobilization share the available mineral and SOM pools', ...
    npass,nfail);
[npass,nfail] = check(P_leak_transfer < -0.1, ...
    sprintf('secondary-P release is available to microbial demand (%.6g)', ...
    P_leak_transfer),npass,nfail);

%% Full DON and DOP leakage leaves nothing for simultaneous stabilization.
organic = F.args;
organic{2}(:) = 0; organic{2}(33) = 1; organic{2}(47) = 1;
organic{5}(:) = 0; organic{52}(:) = 0;
for k = 12:46; organic{k} = 0; end
organic{37} = 1; organic{38} = 1;
dOrganic = BIOGEOCHEMISTRY_DYNAMIC3(organic{:});
organic_new = organic{2}+dOrganic;
[npass,nfail] = check(all(organic_new([33 47]) >= -1e-12), ...
    sprintf('leakage plus stabilization does not overdraw DON or DOP (%.6g)', ...
    min(organic_new([33 47]))),npass,nfail);

%% Enzyme demand cannot withdraw more POC or MOC than is available.
carbon = F.args;
carbon{2}(:) = 0; carbon{2}(9:11) = 1; carbon{2}(14:17) = 1e6;
carbon{5}(:) = 0; carbon{52}(:) = 0;
for k = 12:46; carbon{k} = 0; end
[dCarbon,R_litter_carbon,R_microbe_carbon,~,R_ew_carbon] = ...
    BIOGEOCHEMISTRY_DYNAMIC3(carbon{:});
carbon_new = carbon{2}+dCarbon;
carbon_residual = sum(dCarbon(1:22))+R_litter_carbon+ ...
    R_microbe_carbon+R_ew_carbon;
carbon_element_residual = [sum(dCarbon(23:34)) sum(dCarbon(35:47)) ...
    sum(dCarbon(48:55)) sum(dCarbon(56:63)) ...
    sum(dCarbon(64:71)) sum(dCarbon(72:79))];
[npass,nfail] = check(all(carbon_new >= 0), ...
    sprintf('enzyme sinks leave every coupled elemental pool nonnegative (min %.6g)', ...
    min(carbon_new)),npass,nfail);
[npass,nfail] = check(abs(carbon_residual) < 1e-10 && ...
    all(abs(carbon_element_residual) < 1e-10), ...
    sprintf('POC and MOC donor limiting conserves every element (%.6g)', ...
    max(abs([carbon_residual carbon_element_residual]))),npass,nfail);

%% Earthworm loss cannot release nitrogen that is not present.
earthworm = F.args;
earthworm{2}(:) = 0; earthworm{2}(22) = 1;
earthworm{5}(:) = 0; earthworm{52}(:) = 0;
for k = 12:46; earthworm{k} = 0; end
[dEarthworm,~,~,~,R_ew_loss] = BIOGEOCHEMISTRY_DYNAMIC3(earthworm{:});
earthworm_new = earthworm{2}+dEarthworm;
earthworm_C_residual = sum(dEarthworm(1:22))+R_ew_loss;
[npass,nfail] = check(all(earthworm_new >= 0), ...
    sprintf('earthworm turnover leaves C and N pools nonnegative (min %.6g)', ...
    min(earthworm_new)),npass,nfail);
[npass,nfail] = check(abs(earthworm_C_residual) < 1e-10 && ...
    abs(sum(dEarthworm(23:34))) < 1e-10, ...
    'earthworm donor limiting conserves carbon and nitrogen',npass,nfail);

%% Carbon-only litter does not create pseudo elemental withdrawals.
litter = F.args;
litter{2}(:) = 0; litter{2}(1:8) = 1;
litter{5}(:) = 0; litter{52}(:) = 0;
for k = 12:46; litter{k} = 0; end
[dLitter,R_litter_only,R_microbe_litter,~,R_ew_litter] = ...
    BIOGEOCHEMISTRY_DYNAMIC3(litter{:});
litter_new = litter{2}+dLitter;
litter_C_residual = sum(dLitter(1:22))+R_litter_only+ ...
    R_microbe_litter+R_ew_litter;
[npass,nfail] = check(all(litter_new >= 0), ...
    sprintf('carbon-only litter leaves all elemental pools nonnegative (min %.6g)', ...
    min(litter_new)),npass,nfail);
[npass,nfail] = check(abs(litter_C_residual) < 1e-10, ...
    sprintf('carbon-only litter decomposition conserves carbon (%.6g)', ...
    litter_C_residual),npass,nfail);

%% Carbon-only microbes have no nutrient inventory to mineralize on death.
microbe = F.args;
microbe{2}(:) = 0; microbe{2}(18:21) = 1;
microbe{5}(:) = 0; microbe{52}(:) = 0; microbe{47} = 0.5;
for k = 12:46; microbe{k} = 0; end
[dMicrobe,~,R_microbe_only,~,R_ew_microbe] = ...
    BIOGEOCHEMISTRY_DYNAMIC3(microbe{:});
microbe_new = microbe{2}+dMicrobe;
microbe_C_residual = sum(dMicrobe(1:22))+R_microbe_only+R_ew_microbe;
microbe_element_residual = [sum(dMicrobe(23:34)) sum(dMicrobe(35:47)) ...
    sum(dMicrobe(48:55)) sum(dMicrobe(56:63)) ...
    sum(dMicrobe(64:71)) sum(dMicrobe(72:79))];
[npass,nfail] = check(all(microbe_new >= 0), ...
    sprintf('carbon-only microbial turnover leaves all pools nonnegative (min %.6g)', ...
    min(microbe_new)),npass,nfail);
[npass,nfail] = check(abs(microbe_C_residual) < 1e-10 && ...
    all(abs(microbe_element_residual) < 1e-10), ...
    'carbon-only microbial turnover conserves every element',npass,nfail);

%% Mineralization and mortality share microbial N and P inventories.
stoich = F.args;
stoich{2}(:) = 0; stoich{2}(18:21) = 1;
stoich{2}(27:30) = 1e6; stoich{2}(39:42) = 1e6;
stoich{5}(:) = 0; stoich{52}(:) = 0; stoich{47} = 0.5;
for k = 12:46; stoich{k} = 0; end
[dStoich,~,~,~,~,VOL_stoich,N2_stoich,~,~,~,~,~,~,~,N2O_stoich] = ...
    BIOGEOCHEMISTRY_DYNAMIC3(stoich{:});
stoich_new = stoich{2}+dStoich;
stoich_N_residual = sum(dStoich(23:34))+VOL_stoich+N2_stoich+N2O_stoich;
stoich_P_residual = sum(dStoich(35:47));
[npass,nfail] = check(all(stoich_new([27:30 39:42]) >= 0), ...
    sprintf('mineralization plus mortality does not overdraw microbial N or P (min %.6g)', ...
    min(stoich_new([27:30 39:42]))),npass,nfail);
%% Scale the tolerance to the state: this case loads the microbial pools
%% to 1e6, so a fixed 1e-10 absolute bound sits below machine epsilon for
%% the magnitudes involved and fails on rounding alone.
stoich_tol = 1e-10*max(1,max(abs(stoich{2}([23:34 35:47]))));
[npass,nfail] = check(abs(stoich_N_residual) < stoich_tol && ...
    abs(stoich_P_residual) < stoich_tol, ...
    sprintf('microbial donor limiting conserves N and P (%.6g, %.6g)', ...
    stoich_N_residual,stoich_P_residual),npass,nfail);

%% An empty state has no inventory available for any process rate.
zero = F.args;
zero{2}(:) = 0; zero{5}(:) = 0; zero{52}(:) = 0;
for k = 12:46; zero{k} = 0; end
[dZero,R_litter_zero,R_microbe_zero,R_litter_sur_zero,R_ew_zero, ...
    VOL_zero,N2_zero,MinN_zero,MinP_zero,Rb_zero,Ram_zero,Rem_zero, ...
    ProdB_zero,ProdF_zero,N2O_zero] = BIOGEOCHEMISTRY_DYNAMIC3(zero{:});
zero_rates = [R_litter_zero R_microbe_zero R_litter_sur_zero R_ew_zero ...
    VOL_zero N2_zero MinN_zero MinP_zero Rb_zero Ram_zero Rem_zero ...
    ProdB_zero ProdF_zero N2O_zero];
zero_new = zero{2}+dZero;
[npass,nfail] = check(all(isfinite([dZero(:); zero_rates(:)])), ...
    'empty-state rates and pool changes remain finite',npass,nfail);
[npass,nfail] = check(all(zero_new >= 0), ...
    sprintf('empty-state update leaves all 79 pools nonnegative (min %.6g)', ...
    min(zero_new)),npass,nfail);
zero_handoff_ok = true;
try
    Biogeo_Lateral_Pools(zero_new);
catch
    zero_handoff_ok = false;
end
[npass,nfail] = check(zero_handoff_ok, ...
    'empty-state update is accepted by the lateral pool handoff', ...
    npass,nfail);

%% No water volume and no drainage is a strict zero-leakage state.
b = zeros(1,79);
b([31 32 43 52 12 13 33 47 60 68 76]) = 1;
bp = Biogeochemistry_Parameter(1);
zero_leak = cell(1,10);
[zero_leak{:}] = Biogeo_Leakage(b,0,0,bp);
[npass,nfail] = check(all([zero_leak{:}] == 0), ...
    'a dry zero-drainage layer exports no dissolved pools',npass,nfail);
bad_water_rejected = false;
try
    [zero_leak{:}] = Biogeo_Leakage(b,1,0,bp);
catch ME
    bad_water_rejected = strcmp(ME.identifier,'Biogeo_Leakage:NoWaterStorage');
end
[npass,nfail] = check(bad_water_rejected, ...
    'positive drainage with zero water storage is rejected',npass,nfail);

%% DOC leakage and microbial consumption share each DOC donor pool.
doc = F.args;
doc{2}(:) = 0; doc{2}(12) = 1; doc{2}(18) = 10000;
doc{5}(:) = 0; doc{52}(:) = 0; doc{51}(1) = 0;
for k = 12:46; doc{k} = 0; end
doc{32} = 0.04;
[dDOC,R_litter_DOC,R_microbe_DOC,~,R_ew_DOC] = ...
    BIOGEOCHEMISTRY_DYNAMIC3(doc{:});
DOC_new = doc{2}+dDOC;
DOC_C_residual = sum(dDOC(1:22))+R_litter_DOC+R_microbe_DOC+ ...
    R_ew_DOC+doc{32};
[npass,nfail] = check(DOC_new(12) >= -1e-12, ...
    sprintf('DOC leakage plus microbial use stays within the donor (%.6g)', ...
    DOC_new(12)),npass,nfail);
[npass,nfail] = check(abs(DOC_C_residual) < 1e-10, ...
    sprintf('DOC limiting conserves carbon (%.6g)',DOC_C_residual), ...
    npass,nfail);

%% Physical-state interfaces reject negative inventory rather than hiding it.
tiny_negative = F.args;
tiny_negative{2}(:) = 0;
for k = 25:46; tiny_negative{k} = 0; end
tiny_negative{2}(31) = -eps;
[npass,nfail] = check(rejects(@() BIOGEOCHEMISTRY_DYNAMIC3(tiny_negative{:}), ...
    'BIOGEOCHEMISTRY_DYNAMIC3:InvalidPoolState'), ...
    'biogeochemistry rejects a negative physical pool without clamping', ...
    npass,nfail);
[npass,nfail] = check(rejects(@() Biogeo_Lateral_Pools( ...
    [zeros(1,30) -eps zeros(1,48)]), 'Biogeo_Lateral_Pools:BadPools'), ...
    'lateral handoff rejects a negative physical pool without clamping', ...
    npass,nfail);

bad_env = F.args;
bad_env{55} = 1.01;
[npass,nfail] = check(rejects(@() BIOGEOCHEMISTRY_DYNAMIC3(bad_env{:}), ...
    'BIOGEOCHEMISTRY_DYNAMIC3:InvalidEnvironment'), ...
    'anaerobic fraction outside [0,1] is rejected',npass,nfail);
bad_env = F.args;
bad_env{56} = 1.01;
[npass,nfail] = check(rejects(@() BIOGEOCHEMISTRY_DYNAMIC3(bad_env{:}), ...
    'BIOGEOCHEMISTRY_DYNAMIC3:InvalidEnvironment'), ...
    'water-filled pore space outside [0,1] is rejected',npass,nfail);
bad_env = F.args;
bad_env{10} = 1.01;
[npass,nfail] = check(rejects(@() BIOGEOCHEMISTRY_DYNAMIC3(bad_env{:}), ...
    'BIOGEOCHEMISTRY_DYNAMIC3:InvalidEnvironment'), ...
    'effective saturation outside [0,1] is rejected',npass,nfail);
bad_env = F.args;
bad_env{6} = 70;
[npass,nfail] = check(rejects(@() BIOGEOCHEMISTRY_DYNAMIC3(bad_env{:}), ...
    'BIOGEOCHEMISTRY_DYNAMIC3:InvalidEnvironment'), ...
    'temperature outside the empirical response domain is rejected', ...
    npass,nfail);
bad_env = F.args;
bad_env{9} = NaN;
[npass,nfail] = check(rejects(@() BIOGEOCHEMISTRY_DYNAMIC3(bad_env{:}), ...
    'BIOGEOCHEMISTRY_DYNAMIC3:InvalidEnvironment'), ...
    'non-finite microbial bulk-soil pH is rejected',npass,nfail);

fprintf('  ---------------------------------------\n');
fprintf('  %d passed, %d failed\n',npass,nfail);
if nfail > 0
    error('test_pool_overdraw_guards:Failed','%d check(s) failed.',nfail);
end
end

%% =======================================================================
function [net,Pnew] = closed_N_change(args)
%% Net N change of the soil N pools including gaseous loss, measured
%% against the starting state. External sources and sinks are already
%% zeroed, so this equals the litter N input alone.
[dP,~,~,~,~,VOL,N2flx,~,~,~,~,~,~,~,N2Oflx] = ...
    BIOGEOCHEMISTRY_DYNAMIC3(args{:});
Pnew = args{2} + dP;
net = sum(Pnew(23:34)) - sum(args{2}(23:34)) + VOL + N2flx + N2Oflx;
end

%% =======================================================================
function [npass,nfail] = check(ok,label,npass,nfail)
if ok
    fprintf('  PASS  %s\n',label); npass = npass+1;
else
    fprintf(2,'  FAIL  %s\n',label); nfail = nfail+1;
end
end

%% =======================================================================
function tf = rejects(fn,id)
tf = false;
try
    fn();
catch ME
    tf = strcmp(ME.identifier,id);
end
end

%% =======================================================================
function E = erw_columns(N)
E = cell(N,1);
for c = 1:N
    E{c}.init = true;
    E{c}.dm_end = 1;
    E{c}.current_s = 0.60;
    E{c}.const.n = 0.45;
    E{c}.const.Zrw = 0.30;
    E{c}.bg.Na = 100 - 5*c;
    E{c}.bg.Al = 20 - c;
    E{c}.bg.AlOH = 2;
    E{c}.bg.AlOH2 = 1;
    E{c}.bg.AlOH3 = 0.5;
    E{c}.bg.AlOH4 = 5;
    E{c}.bg.An = 50 - 2*c;
    E{c}.bg.DIC = 500 - 10*c;
    %% Uniform soil-gas CO2: gas diffusion is driven by its own gradient,
    %% not by water, so it would legitimately move even at zero water flux.
    E{c}.bg.CO2_air = 2000;
    E{c}.bg.D = 0.50;
    E{c}.bg.Na_tot = 1e5;
    E{c}.bg.Al_tot = 5e4;
    E{c}.bg.An_tot = 5e4;
    E{c}.bg.IC_tot = 5e5;
    E{c}.bg.Alk_tot = 5e4;
end
end
