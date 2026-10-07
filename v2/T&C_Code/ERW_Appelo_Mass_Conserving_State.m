function state = ERW_Appelo_Mass_Conserving_State(Alk_tot,IC_tot,cation_tot, ...
    soil_water_l,bulk_volume_l,s,carbonate,CEC_tot,appelo,gas,al3_fn,~,finite_al, ...
    cec_fn)
%%% Solve Appelo exchange and carbonate charge from conserved pools.

state = struct('is_valid',false,'H',NaN,'CO2_w',NaN,'Alk',NaN, ...
    'R_alk',NaN,'Al3',NaN,'fraction',NaN(1,6), ...
    'concentration',NaN(1,4),'IC_tot',NaN,'Fs',NaN,'charge_residual',Inf, ...
    'AlOH3_solid',NaN,'Al_residual',NaN);

if nargin < 13
    finite_al = [];
end
%%% Optional pH-dependent site density.  Omitted or empty keeps the constant
%%% CEC_tot, so the legacy call signature is unchanged.  This routine locates
%%% the charge-balance root by scanning H, so it has to see the same site
%%% density at each trial H as the residual system it seeds, otherwise the
%%% seed it returns violates cation mass balance by the CEC difference.
if nargin < 14 || isempty(cec_fn)
    cec_fn = @(H) CEC_tot;
end
if ~isa(cec_fn,'function_handle')
    error('ERW:invalidParam','cec_fn must be a function handle of H.');
end

if numel(cation_tot) ~= 4 || numel(carbonate) ~= 5 || numel(appelo) ~= 6 || numel(gas) ~= 3 || ...
        any(~isfinite([Alk_tot,IC_tot,cation_tot(:).',soil_water_l, ...
        bulk_volume_l,s,carbonate(:).',CEC_tot,appelo(:).',gas(:).'])) || ...
        IC_tot < 0 || any(cation_tot < 0) || soil_water_l <= 0 || bulk_volume_l <= 0 || ...
        s <= 0 || s > 1 || CEC_tot <= 0 || carbonate(4) <= 0 || ...
        gas(1) < 0 || gas(2) < 0 || gas(3) <= 0
    return
end

conv_mol = carbonate(5);
pH_scan = [1e-8:0.25:13.75,14-1e-8];
H_scan = 10.^(-pH_scan)*conv_mol;
residual = NaN(size(H_scan));
for ii = 1:numel(H_scan)
    residual(ii) = charge_residual(H_scan(ii));
end
bracket = find(isfinite(residual(1:end-1)) & isfinite(residual(2:end)) & ...
    (residual(1:end-1) == 0 | residual(2:end) == 0 | ...
    sign(residual(1:end-1)) ~= sign(residual(2:end))));
bracket = bracket([true,diff(bracket) > 1]);
if numel(bracket) ~= 1
    return
end

ib = bracket(1);
a = H_scan(ib);
b = H_scan(ib+1);
ra = charge_residual(a);
for it = 1:60
    H_mid = 0.5*(a+b);
    r_mid = charge_residual(H_mid);
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
[r,state_at_H] = charge_residual(H);
state = state_at_H;
state.charge_residual = r;
state.is_valid = all(isfinite([state.H,state.CO2_w,state.R_alk,state.Al3, ...
    state.fraction,state.concentration,state.IC_tot,state.Fs, ...
    state.AlOH3_solid,state.Al_residual])) && ...
    state.H > 0 && state.CO2_w >= 0 && state.IC_tot >= 0 && ...
    state.R_alk >= 0 && all(state.fraction >= 0) && ...
    all(state.concentration >= 0) && abs(sum(state.fraction)-1) <= 1e-10 && ...
    abs(r) <= 1e-6;

    function [r,out] = charge_residual(H)
        out = state;
        k1 = carbonate(1); k2 = carbonate(2);
        k_w = carbonate(3); k_H = carbonate(4);
        K_Na_H = appelo(1); K_K_H = appelo(2);
        K_Ca_H2 = appelo(3); K_Mg_H2 = appelo(4);
        K_Al_H3 = appelo(5); alpha = appelo(6);
        dic_factor = bulk_volume_l*(s*(1+k1/H+k1*k2/H^2)+(1-s)/k_H);
        CO2_w = (IC_tot+gas(3)*gas(1)*gas(2))/(dic_factor+gas(3)*gas(1)/k_H);
        Fs = gas(1)*(CO2_w/k_H-gas(2));
        Al3 = al3_fn(H);
        cec_H = cec_fn(H);
        [f_vec,conc_vec,Al3,AlOH3_solid,Al_residual] = ...
            ERW_Appelo_State_At_H(H,cation_tot,Al3,cec_H,soil_water_l, ...
            alpha,K_Na_H,K_K_H,K_Ca_H2,K_Mg_H2,K_Al_H3,finite_al);
        R_alk = sum(f_vec(1:4))*cec_H;
        Al_charge = 3*Al3;
        if ~isempty(finite_al)
            Al_charge = Al3*(3+2*finite_al.K1/H+finite_al.K1*finite_al.K2/H^2- ...
                finite_al.K1*finite_al.K2*finite_al.K3*finite_al.K4/H^4);
        end
        Alk = (Alk_tot-R_alk)/soil_water_l+Al_charge;
        r = k1*CO2_w/H+2*k1*k2*CO2_w/H^2-H+k_w/H-Alk;
        out.H = H;
        out.CO2_w = CO2_w;
        out.Alk = Alk;
        out.R_alk = R_alk;
        out.Al3 = Al3;
        out.AlOH3_solid = AlOH3_solid;
        out.Al_residual = Al_residual;
        out.fraction = f_vec;
        out.concentration = conc_vec;
        out.IC_tot = CO2_w*dic_factor;
        out.Fs = Fs;
    end
end
