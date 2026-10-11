function [f_vec, conc_vec, Al3, AlOH3_solid, Al_residual] = ...
    ERW_Appelo_State_At_H(H, cation_tot, Al3, CEC_tot, soil_water_l, ...
    alpha, K_Na_H, K_K_H, K_Ca_H2, K_Mg_H2, K_Al_H3, finite_al)
%%% Appelo competitive exchange partition at fixed pore-water H.

f_vec = NaN(1,6);
conc_vec = NaN(1,4);
AlOH3_solid = NaN;
Al_residual = NaN;
if nargin < 12
    finite_al = [];
end
if ~(isfinite(H) && H > 0 && isfinite(soil_water_l) && ...
        soil_water_l > 0 && isfinite(CEC_tot) && CEC_tot > 0) || ...
        numel(cation_tot) ~= 4 || any(~isfinite(cation_tot)) || ...
        any(cation_tot < 0)
    return
end
if isempty(finite_al) && ~(isfinite(Al3) && Al3 >= 0)
    return
end
if ~isempty(finite_al) && ~local_valid_finite_al(finite_al)
    return
end

cec_per_L = CEC_tot/soil_water_l;
X = cation_tot(:).'/soil_water_l;             %% Ca Mg K Na [umol L-1]
lo = 0;
hi = 1;
for ii = 1:60
    fh = 0.5*(lo+hi);
    [gap,f_try,c_try,Al3_try,solid_try,residual_try] = ...
        local_site_gap(fh,H,X,Al3,cec_per_L,CEC_tot,soil_water_l, ...
        alpha,K_Na_H,K_K_H,K_Ca_H2,K_Mg_H2,K_Al_H3,finite_al); %#ok<ASGLU>
    if ~isfinite(gap)
        return
    elseif gap > 0
        hi = fh;
    else
        lo = fh;
    end
end

[gap,f_vec,conc_vec,Al3,AlOH3_solid,Al_residual] = ...
    local_site_gap(0.5*(lo+hi),H,X,Al3,cec_per_L,CEC_tot,soil_water_l, ...
    alpha,K_Na_H,K_K_H,K_Ca_H2,K_Mg_H2,K_Al_H3,finite_al);
if ~isfinite(gap) || abs(gap) > 1e-12
    f_vec(:) = NaN;
    conc_vec(:) = NaN;
end
end

function [gap,f_vec,conc_vec,Al3,AlOH3_solid,Al_residual] = ...
    local_site_gap(fh,H,X,Al3,cec_per_L,CEC_tot,soil_water_l,alpha, ...
    K_Na_H,K_K_H,K_Ca_H2,K_Mg_H2,K_Al_H3,finite_al)

rho = (fh/H)*10^(alpha*fh);
A = [K_Ca_H2*rho^2,K_Mg_H2*rho^2,K_K_H*rho,K_Na_H*rho];
z = [2 2 1 1];
conc_vec = X./(1+A.*cec_per_L./z);
f4 = A.*conc_vec;
if isempty(finite_al)
    f_Al = K_Al_H3*Al3*rho^3;
    AlOH3_solid = 0;
    Al_residual = 0;
else
    [f_Al,Al3,AlOH3_solid,Al_residual] = ...
        local_finite_al(H,rho,CEC_tot,soil_water_l,K_Al_H3,finite_al);
end
f_vec = [f4,f_Al,fh];
gap = sum(f_vec)-1;
end

function [f_Al,Al3,AlOH3_solid,Al_residual] = ...
    local_finite_al(H,rho,CEC_tot,soil_water_l,K_Al_H3,finite_al)

Al_factor = 1 + finite_al.K1/H + finite_al.K1*finite_al.K2/H^2 + ...
    finite_al.K1*finite_al.K2*finite_al.K3/H^3 + ...
    finite_al.K1*finite_al.K2*finite_al.K3*finite_al.K4/H^4;
Al3_sat_scaled = 10^(finite_al.K_AlOH3-12)*H^3*finite_al.Al_CEC_scale;
exchange_factor = K_Al_H3*rho^3;
f_sat = exchange_factor*Al3_sat_scaled/finite_al.Al_CEC_scale;
Al_sat = Al3_sat_scaled*Al_factor*soil_water_l + ...
    f_sat*CEC_tot*finite_al.Al_CEC_scale/3;
if finite_al.Al_tot <= 0 || exchange_factor <= 0
    f_Al = 0;
elseif isfinite(Al_sat) && Al_sat <= finite_al.Al_tot
    f_Al = f_sat;
else
    f_Al = finite_al.Al_tot/(finite_al.Al_CEC_scale*Al_factor* ...
        soil_water_l/exchange_factor + CEC_tot*finite_al.Al_CEC_scale/3);
end
[Al3_scaled,AlOH3_solid,Al_residual] = ERW_Al_Equilibrium_Target( ...
    H,f_Al,finite_al.Al_tot,soil_water_l,CEC_tot, ...
    finite_al.Al_CEC_scale,finite_al.K_AlOH3,finite_al.K1, ...
    finite_al.K2,finite_al.K3,finite_al.K4);
Al3 = Al3_scaled/finite_al.Al_CEC_scale;
end

function valid = local_valid_finite_al(finite_al)
required = {'Al_tot','Al_CEC_scale','K_AlOH3','K1','K2','K3','K4'};
valid = isstruct(finite_al) && all(isfield(finite_al,required));
if valid
    values = [finite_al.Al_tot,finite_al.Al_CEC_scale,finite_al.K_AlOH3, ...
        finite_al.K1,finite_al.K2,finite_al.K3,finite_al.K4];
    valid = all(isfinite(values)) && finite_al.Al_tot >= 0 && ...
        finite_al.Al_CEC_scale > 0 && all(values(4:7) >= 0);
end
end
