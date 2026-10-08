function [state,status] = ERW_Appelo_Seed_Kernel(H_scan,Alk_tot,IC_tot, ...
    cation_tot,soil_water_l,bulk_volume_l,s,carbonate,CEC_tot,appelo,gas, ...
    al3_K_AlOH3,al3_Al_CEC_scale,finite_al,has_finite_al) %#codegen
%%% ERW_Appelo_Mass_Conserving_State for the call ERW_function makes, with
%%% the same arithmetic in the same order and none of its call overhead.
%%%
%%% ERW_function calls the reduced Appelo solve on every ten-minute step
%%% under the Al(OH)3 pin.  The original locates the charge root by
%%% scanning 57 pH points and bisecting, and every trial H runs a fixed
%%% 60-step bisection on f_H through three nested functions, so one call
%%% makes about 5,700 inner evaluations and 17,000 function calls.  This
%%% kernel evaluates exactly the same floating-point expressions on the
%%% same operands; it only removes work whose result the original throws
%%% away or repeats:
%%%   - the inner bisection needs only the site gap, so the Al(OH)3
%%%     inventory split (ERW_Al_Equilibrium_Target) runs once per trial H,
%%%     at the final f_H, instead of at all 61 f_H values;
%%%   - terms that depend on H alone are formed once per trial H, as the
%%%     leading operands they already are;
%%%   - the residual at the lower bracket end is the scan's value at that H;
%%%   - the 57 scan points run their inner bisections side by side, element
%%%     by element the same operations (MATLAB's elementwise +, -, *, /,
%%%     .^ and 10.^ round exactly as its scalar forms).
%%% Every result is therefore bit-identical to the original, which
%%% test_erw_chemistry_kernels checks with isequaln over the model's range,
%%% for this MATLAB source and for the MEX built from it
%%% (Build_ERW_Chemistry_Mex).
%%%
%%% Inputs are the original's, with the two function handles replaced by
%%% the forms ERW_function passes:
%%%   H_scan            10.^(-[1e-8:0.25:13.75,14-1e-8])*conv_mol, formed by
%%%                     the caller with MATLAB's own colon and power
%%%   al3_K_AlOH3,      al3_fn(H) = 10^(al3_K_AlOH3-12)*max(H,1e-12)^3 ...
%%%   al3_Al_CEC_scale      *al3_Al_CEC_scale/al3_Al_CEC_scale
%%%   finite_al         [Al_tot Al_CEC_scale K_AlOH3 K1 K2 K3 K4]
%%%   has_finite_al     false for the original's empty finite_al
%%% The exchanger capacity is constant (cec_fn(H) = CEC_tot), the only form
%%% the Appelo scheme admits.
%%%
%%% status 1 means state is the original's result.  status 0 means the
%%% scan found no sign change, where the original stops with
%%% MATLAB:badsubscript on bracket([true,diff(bracket) > 1]); the caller
%%% reruns the original to raise exactly that error.

state = struct('is_valid',false,'H',NaN,'CO2_w',NaN,'Alk',NaN, ...
    'R_alk',NaN,'Al3',NaN,'fraction',NaN(1,6), ...
    'concentration',NaN(1,4),'IC_tot',NaN,'Fs',NaN,'charge_residual',Inf, ...
    'AlOH3_solid',NaN,'Al_residual',NaN);
status = 1;

if numel(cation_tot) ~= 4 || numel(carbonate) ~= 5 || numel(appelo) ~= 6 || numel(gas) ~= 3 || ...
        any(~isfinite([Alk_tot,IC_tot,cation_tot(:).',soil_water_l, ...
        bulk_volume_l,s,carbonate(:).',CEC_tot,appelo(:).',gas(:).'])) || ...
        IC_tot < 0 || any(cation_tot < 0) || soil_water_l <= 0 || bulk_volume_l <= 0 || ...
        s <= 0 || s > 1 || CEC_tot <= 0 || carbonate(4) <= 0 || ...
        gas(1) < 0 || gas(2) < 0 || gas(3) <= 0
    return
end

P = seed_parameters(Alk_tot,IC_tot,cation_tot,soil_water_l,bulk_volume_l, ...
    s,carbonate,CEC_tot,appelo,gas,al3_K_AlOH3,al3_Al_CEC_scale, ...
    finite_al,has_finite_al);

%% The scan: every point's f_H bisection at once, then its residual.
n_scan = numel(H_scan);
[lo_scan,hi_scan,stopped] = seed_bisect_scan(H_scan,P.CEC_tot,P);
residual = NaN(1,n_scan);
for ii = 1:n_scan
    residual(ii) = seed_charge_residual(H_scan(ii),P, ...
        lo_scan(ii),hi_scan(ii),stopped(ii));
end
%% find(...) followed by bracket([true,diff(bracket) > 1]): the first index
%% of every run of consecutive sign changes.
n_runs = 0;
ib = 0;
previous = -1;
for k = 1:n_scan-1
    r0 = residual(k);
    r1 = residual(k+1);
    if isfinite(r0) && isfinite(r1) && ...
            (r0 == 0 || r1 == 0 || sign(r0) ~= sign(r1))
        if n_runs == 0 || k - previous > 1
            n_runs = n_runs + 1;
            if n_runs == 1
                ib = k;
            end
        end
        previous = k;
    end
end
if n_runs == 0
    status = 0;
    return
end
if n_runs ~= 1
    return
end

a = H_scan(ib);
b = H_scan(ib+1);
ra = residual(ib);
for it = 1:60
    H_mid = 0.5*(a+b);
    [lo,hi,stop] = seed_bisect(H_mid,P.CEC_tot,P);
    r_mid = seed_charge_residual(H_mid,P,lo,hi,stop);
    if abs(r_mid) <= 1e-8
        break
    elseif sign(r_mid) == sign(ra)
        a = H_mid;
        ra = r_mid;
    else
        b = H_mid;
    end
end
H = 0.5*(a+b);
[lo,hi,stop] = seed_bisect(H,P.CEC_tot,P);
[r,CO2_w,Fs,Alk,R_alk,Al3,AlOH3_solid,Al_residual,f_vec,conc_vec, ...
    dic_factor] = seed_charge_residual(H,P,lo,hi,stop);
%% state_at_H is the initial state with the trial-H fields filled in.
state.H = H;
state.CO2_w = CO2_w;
state.Alk = Alk;
state.R_alk = R_alk;
state.Al3 = Al3;
state.AlOH3_solid = AlOH3_solid;
state.Al_residual = Al_residual;
state.fraction = f_vec;
state.concentration = conc_vec;
state.IC_tot = CO2_w*dic_factor;
state.Fs = Fs;
state.charge_residual = r;
state.is_valid = all(isfinite([state.H,state.CO2_w,state.R_alk,state.Al3, ...
    state.fraction,state.concentration,state.IC_tot,state.Fs, ...
    state.AlOH3_solid,state.Al_residual])) && ...
    state.H > 0 && state.CO2_w >= 0 && state.IC_tot >= 0 && ...
    state.R_alk >= 0 && all(state.fraction >= 0) && ...
    all(state.concentration >= 0) && abs(sum(state.fraction)-1) <= 1e-10 && ...
    abs(r) <= 1e-6;
end

function P = seed_parameters(Alk_tot,IC_tot,cation_tot,soil_water_l, ...
    bulk_volume_l,s,carbonate,CEC_tot,appelo,gas,al3_K_AlOH3, ...
    al3_Al_CEC_scale,finite_al,has_finite_al)
%%% Everything the trial-H evaluation reads that does not depend on H.
X = cation_tot(:).'/soil_water_l;             %% Ca Mg K Na [umol L-1]
P = struct('Alk_tot',Alk_tot,'IC_tot',IC_tot, ...
    'cation_ok',numel(cation_tot) == 4 && ~any(~isfinite(cation_tot)) && ...
        ~any(cation_tot < 0), ...
    'X1',X(1),'X2',X(2),'X3',X(3),'X4',X(4), ...
    'soil_water_l',soil_water_l,'bulk_volume_l',bulk_volume_l,'s',s, ...
    'k1',carbonate(1),'k2',carbonate(2),'k_w',carbonate(3),'k_H',carbonate(4), ...
    'CEC_tot',CEC_tot, ...
    'K_Na_H',appelo(1),'K_K_H',appelo(2),'K_Ca_H2',appelo(3), ...
    'K_Mg_H2',appelo(4),'K_Al_H3',appelo(5),'alpha',appelo(6), ...
    'gas1',gas(1),'gas2',gas(2),'gas3',gas(3), ...
    'al3_K_AlOH3',al3_K_AlOH3,'al3_Al_CEC_scale',al3_Al_CEC_scale, ...
    'has_finite_al',logical(has_finite_al), ...
    'finite_ok',false,'Al_tot',0,'Al_CEC_scale',0,'K_AlOH3',0, ...
    'K1',0,'K2',0,'K3',0,'K4',0);
if P.has_finite_al && numel(finite_al) == 7
    %% local_valid_finite_al, which depends on finite_al alone.
    P.finite_ok = all(isfinite(finite_al)) && ...
        finite_al(1) >= 0 && finite_al(2) > 0 && all(finite_al(4:7) >= 0);
    P.Al_tot = finite_al(1);
    P.Al_CEC_scale = finite_al(2);
    P.K_AlOH3 = finite_al(3);
    P.K1 = finite_al(4);
    P.K2 = finite_al(5);
    P.K3 = finite_al(6);
    P.K4 = finite_al(7);
end
end

function Al3 = seed_al3(H,P)
%%% al3_fn(H) = al3_pin_fn(H)/Al_CEC_scale as ERW_function defines them.
Al3 = 10^(P.al3_K_AlOH3-12) * max(real(H), 1e-12)^3 * P.al3_Al_CEC_scale / ...
    P.al3_Al_CEC_scale;
end

function [r,CO2_w,Fs,Alk,R_alk,Al3,AlOH3_solid,Al_residual,f_vec, ...
    conc_vec,dic_factor] = seed_charge_residual(H,P,lo,hi,stopped)
%%% charge_residual of ERW_Appelo_Mass_Conserving_State, given the trial
%%% H's finished f_H bisection.
k1 = P.k1; k2 = P.k2;
k_w = P.k_w; k_H = P.k_H;
dic_factor = P.bulk_volume_l*(P.s*(1+k1/H+k1*k2/H^2)+(1-P.s)/k_H);
CO2_w = (P.IC_tot+P.gas3*P.gas1*P.gas2)/(dic_factor+P.gas3*P.gas1/k_H);
Fs = P.gas1*(CO2_w/k_H-P.gas2);
cec_H = P.CEC_tot;
[f_vec,conc_vec,Al3,AlOH3_solid,Al_residual] = ...
    seed_state_final(H,cec_H,P,lo,hi,stopped);
R_alk = sum(f_vec(1:4))*cec_H;
Al_charge = 3*Al3;
if P.has_finite_al
    Al_charge = Al3*(3+2*P.K1/H+P.K1*P.K2/H^2- ...
        P.K1*P.K2*P.K3*P.K4/H^4);
end
Alk = (P.Alk_tot-R_alk)/P.soil_water_l+Al_charge;
r = k1*CO2_w/H+2*k1*k2*CO2_w/H^2-H+k_w/H-Alk;
end

function ok = seed_state_admissible(H,CEC_tot,P)
%%% The early returns of ERW_Appelo_State_At_H that precede its bisection.
ok = (isfinite(H) && H > 0 && isfinite(P.soil_water_l) && ...
    P.soil_water_l > 0 && isfinite(CEC_tot) && CEC_tot > 0) && P.cation_ok;
if ~ok
    return
end
if P.has_finite_al
    ok = P.finite_ok;
else
    Al3 = seed_al3(H,P);
    ok = isfinite(Al3) && Al3 >= 0;
end
end

function [lo,hi,stopped] = seed_bisect(H,CEC_tot,P)
%%% The 60-step f_H bisection of ERW_Appelo_State_At_H at one trial H.
%%% stopped marks the original's early returns: inadmissible inputs or a
%%% non-finite site gap.
lo = 0;
hi = 1;
stopped = ~seed_state_admissible(H,CEC_tot,P);
if stopped
    return
end
soil_water_l = P.soil_water_l;
cec_per_L = CEC_tot/soil_water_l;
alpha = P.alpha;
K_Ca_H2 = P.K_Ca_H2; K_Mg_H2 = P.K_Mg_H2; K_K_H = P.K_K_H; K_Na_H = P.K_Na_H;
K_Al_H3 = P.K_Al_H3;
X1 = P.X1; X2 = P.X2; X3 = P.X3; X4 = P.X4;
has_finite_al = P.has_finite_al;
Al_tot = P.Al_tot;
Al_CEC_scale = P.Al_CEC_scale;
%% local_finite_al terms that depend on H alone.
Al3 = 0; Al3_sat_scaled = 0; Al_sat_fixed = 0; Al_den_fixed = 0;
Al_cap = 0; Al_none = true;
if has_finite_al
    Al_factor = 1 + P.K1/H + P.K1*P.K2/H^2 + ...
        P.K1*P.K2*P.K3/H^3 + ...
        P.K1*P.K2*P.K3*P.K4/H^4;
    Al3_sat_scaled = 10^(P.K_AlOH3-12)*H^3*Al_CEC_scale;
    Al_sat_fixed = Al3_sat_scaled*Al_factor*soil_water_l;
    Al_den_fixed = Al_CEC_scale*Al_factor*soil_water_l;
    Al_cap = CEC_tot*Al_CEC_scale/3;
    Al_none = Al_tot <= 0;
else
    Al3 = seed_al3(H,P);
end
for ii = 1:60
    fh = 0.5*(lo+hi);
    %% local_site_gap
    rho = (fh/H)*10^(alpha*fh);
    rho2 = rho^2;
    A1 = K_Ca_H2*rho2; A2 = K_Mg_H2*rho2; A3 = K_K_H*rho; A4 = K_Na_H*rho;
    c1 = X1/(1+A1*cec_per_L/2);
    c2 = X2/(1+A2*cec_per_L/2);
    c3 = X3/(1+A3*cec_per_L/1);
    c4 = X4/(1+A4*cec_per_L/1);
    if has_finite_al
        %% local_finite_al
        exchange_factor = K_Al_H3*rho^3;
        f_sat = exchange_factor*Al3_sat_scaled/Al_CEC_scale;
        Al_sat = Al_sat_fixed + f_sat*CEC_tot*Al_CEC_scale/3;
        if Al_none || exchange_factor <= 0
            f_Al = 0;
        elseif isfinite(Al_sat) && Al_sat <= Al_tot
            f_Al = f_sat;
        else
            f_Al = Al_tot/(Al_den_fixed/exchange_factor + Al_cap);
        end
    else
        f_Al = K_Al_H3*Al3*rho^3;
    end
    gap = A1*c1+A2*c2+A3*c3+A4*c4+f_Al+fh-1;
    if ~isfinite(gap)
        stopped = true;
        return
    elseif gap > 0
        hi = fh;
    else
        lo = fh;
    end
end
end

function [lo,hi,stopped] = seed_bisect_scan(H,CEC_tot,P)
%%% seed_bisect at every scan H at once.  Each element goes through the
%%% scalar loop's operations; one whose site gap turns non-finite keeps
%%% its bracket from then on, as the scalar loop returns there.
n = numel(H);
lo = zeros(1,n);
hi = ones(1,n);
stopped = false(1,n);
for k = 1:n
    stopped(k) = ~seed_state_admissible(H(k),CEC_tot,P);
end
soil_water_l = P.soil_water_l;
cec_per_L = CEC_tot/soil_water_l;
alpha = P.alpha;
K_Ca_H2 = P.K_Ca_H2; K_Mg_H2 = P.K_Mg_H2; K_K_H = P.K_K_H; K_Na_H = P.K_Na_H;
K_Al_H3 = P.K_Al_H3;
X1 = P.X1; X2 = P.X2; X3 = P.X3; X4 = P.X4;
has_finite_al = P.has_finite_al;
Al_tot = P.Al_tot;
Al_CEC_scale = P.Al_CEC_scale;
Al3 = zeros(1,n); Al3_sat_scaled = zeros(1,n); Al_sat_fixed = zeros(1,n);
Al_den_fixed = zeros(1,n); Al_cap = 0; Al_none = true;
if has_finite_al
    Al_factor = 1 + P.K1./H + P.K1*P.K2./H.^2 + ...
        P.K1*P.K2*P.K3./H.^3 + ...
        P.K1*P.K2*P.K3*P.K4./H.^4;
    Al3_sat_scaled = 10^(P.K_AlOH3-12)*H.^3*Al_CEC_scale;
    Al_sat_fixed = Al3_sat_scaled.*Al_factor*soil_water_l;
    Al_den_fixed = Al_CEC_scale*Al_factor*soil_water_l;
    Al_cap = CEC_tot*Al_CEC_scale/3;
    Al_none = Al_tot <= 0;
else
    for k = 1:n
        Al3(k) = seed_al3(H(k),P);
    end
end
for ii = 1:60
    fh = 0.5*(lo+hi);
    rho = (fh./H).*10.^(alpha*fh);
    rho2 = rho.^2;
    A1 = K_Ca_H2*rho2; A2 = K_Mg_H2*rho2; A3 = K_K_H*rho; A4 = K_Na_H*rho;
    c1 = X1./(1+A1*cec_per_L/2);
    c2 = X2./(1+A2*cec_per_L/2);
    c3 = X3./(1+A3*cec_per_L/1);
    c4 = X4./(1+A4*cec_per_L/1);
    if has_finite_al
        exchange_factor = K_Al_H3*rho.^3;
        f_sat = exchange_factor.*Al3_sat_scaled/Al_CEC_scale;
        Al_sat = Al_sat_fixed + f_sat*CEC_tot*Al_CEC_scale/3;
        %% The scalar branches, last-assigned first.
        f_Al = Al_tot./(Al_den_fixed./exchange_factor + Al_cap);
        use_sat = isfinite(Al_sat) & Al_sat <= Al_tot;
        f_Al(use_sat) = f_sat(use_sat);
        use_zero = Al_none | exchange_factor <= 0;
        f_Al(use_zero) = 0;
    else
        f_Al = K_Al_H3*Al3.*rho.^3;
    end
    gap = A1.*c1+A2.*c2+A3.*c3+A4.*c4+f_Al+fh-1;
    stopped = stopped | ~isfinite(gap);
    to_hi = ~stopped & gap > 0;
    to_lo = ~stopped & ~(gap > 0);
    hi(to_hi) = fh(to_hi);
    lo(to_lo) = fh(to_lo);
end
end

function [f_vec,conc_vec,Al3,AlOH3_solid,Al_residual] = ...
    seed_state_final(H,CEC_tot,P,lo,hi,stopped)
%%% ERW_Appelo_State_At_H after its bisection: the final local_site_gap at
%%% f_H = 0.5*(lo+hi), with the Al(OH)3 partition.  Its Al3 input is
%%% al3_fn(H), which the finite-Al partition replaces; it is formed here
%%% only where the original returns it unchanged.
f_vec = NaN(1,6);
conc_vec = NaN(1,4);
AlOH3_solid = NaN;
Al_residual = NaN;
if stopped
    Al3 = seed_al3(H,P);
    return
end
soil_water_l = P.soil_water_l;
cec_per_L = CEC_tot/soil_water_l;
fh = 0.5*(lo+hi);
rho = (fh/H)*10^(P.alpha*fh);
rho2 = rho^2;
A1 = P.K_Ca_H2*rho2; A2 = P.K_Mg_H2*rho2; A3 = P.K_K_H*rho; A4 = P.K_Na_H*rho;
c1 = P.X1/(1+A1*cec_per_L/2);
c2 = P.X2/(1+A2*cec_per_L/2);
c3 = P.X3/(1+A3*cec_per_L/1);
c4 = P.X4/(1+A4*cec_per_L/1);
if P.has_finite_al
    Al_factor = 1 + P.K1/H + P.K1*P.K2/H^2 + ...
        P.K1*P.K2*P.K3/H^3 + ...
        P.K1*P.K2*P.K3*P.K4/H^4;
    Al3_sat_scaled = 10^(P.K_AlOH3-12)*H^3*P.Al_CEC_scale;
    exchange_factor = P.K_Al_H3*rho^3;
    f_sat = exchange_factor*Al3_sat_scaled/P.Al_CEC_scale;
    Al_sat = Al3_sat_scaled*Al_factor*soil_water_l + ...
        f_sat*CEC_tot*P.Al_CEC_scale/3;
    if P.Al_tot <= 0 || exchange_factor <= 0
        f_Al = 0;
    elseif isfinite(Al_sat) && Al_sat <= P.Al_tot
        f_Al = f_sat;
    else
        f_Al = P.Al_tot/(P.Al_CEC_scale*Al_factor* ...
            soil_water_l/exchange_factor + CEC_tot*P.Al_CEC_scale/3);
    end
    [Al3_scaled,AlOH3_solid,Al_residual] = ERW_Al_Equilibrium_Target( ...
        H,f_Al,P.Al_tot,soil_water_l,CEC_tot, ...
        P.Al_CEC_scale,P.K_AlOH3,P.K1, ...
        P.K2,P.K3,P.K4);
    Al3 = Al3_scaled/P.Al_CEC_scale;
else
    Al3 = seed_al3(H,P);
    f_Al = P.K_Al_H3*Al3*rho^3;
    AlOH3_solid = 0;
    Al_residual = 0;
end
f_vec = [A1*c1,A2*c2,A3*c3,A4*c4,f_Al,fh];
conc_vec = [c1,c2,c3,c4];
gap = sum(f_vec)-1;
if ~isfinite(gap) || abs(gap) > 1e-12
    f_vec(:) = NaN;
    conc_vec(:) = NaN;
end
end
