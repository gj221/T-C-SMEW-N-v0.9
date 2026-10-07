function [Al3_target,AlOH3_solid,inventory_residual] = ERW_Al_Equilibrium_Target( ...
    H,f_Al,Al_tot,soil_water_l,CEC_tot,Al_CEC_scale,K_AlOH3,K1,K2,K3,K4)
% Finite Al(OH)3 equilibrium on the ERW scaled-Al storage basis.

H = max(real(H),1e-12);
f_Al = max(real(f_Al),0);
Al_exchange = f_Al*CEC_tot*Al_CEC_scale/3;
Al_available = max(real(Al_tot)-Al_exchange,0);
Al_den = H^4+H^3*K1+H^2*K1*K2+H*K1*K2*K3+K1*K2*K3*K4;
Al3_fraction = H^4/Al_den;
Al3_saturation = 10^(K_AlOH3-12)*H^3*Al_CEC_scale;
Al3_inventory = Al3_fraction*Al_available/soil_water_l;
Al3_target = min(Al3_saturation,Al3_inventory);
Al_dissolved = Al3_target/max(Al3_fraction,realmin)*soil_water_l;
AlOH3_solid = max(real(Al_tot)-Al_exchange-Al_dissolved,0);
inventory_residual = Al_exchange+Al_dissolved+AlOH3_solid-real(Al_tot);
end
