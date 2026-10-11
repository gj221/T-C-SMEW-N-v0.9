function P_row = ERW_Initialize_TC_Cation_Pools(P_row, bg, cst, s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Synchronize the T&C cation pools with the finite ERW seed inventory.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

water_l = cst.n*cst.Zrw*s*1000;
P_row(60) = bg.Ca(1)*40*water_l/1e6;
P_row(61) = bg.f_Ca(1)*bg.CEC_tot(1)/2*40/1e6;
P_row(68) = bg.Mg(1)*24.3*water_l/1e6;
P_row(69) = bg.f_Mg(1)*bg.CEC_tot(1)/2*24.3/1e6;
P_row(52) = bg.K(1)*39.1*water_l/1e6;
P_row(53) = bg.f_K(1)*bg.CEC_tot(1)*39.1/1e6;

end
