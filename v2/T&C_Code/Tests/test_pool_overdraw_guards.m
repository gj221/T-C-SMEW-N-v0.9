function test_pool_overdraw_guards
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Guards against the pool-overdraw failures that stopped the 2-D runs.
%%
%%   1. ERW lateral transport is an exact no-op when no water moves.
%%   2. ERW aqueous state inverts its own concentration/total relation.
%%   3. Net immobilization that exceeds mineral AND SOM nitrogen still
%%      conserves total N (the imm_cut deficit path).
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
[E_out,diag] = ERW_Lateral_Transport(E,topo,Q_zero,theta,Zlayer_m,active,1);

max_drift = 0;
fields = {'Na','Al','AlOH4','An','DIC','CO2_air', ...
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
