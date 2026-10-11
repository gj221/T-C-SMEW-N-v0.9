function [Omega, Wr, EW] = ERW_Calcite_Kinetics(Ca, CO3, K_CaCO3, H, ...
    CO2_air, conv_mol, s, diss_f, k_H, k_w, k_CO2, n_H, n_CO2, SA, ...
    rock_fraction, M_calcite, MM_calcite, dtw)

values = [Ca,CO3,K_CaCO3,H,CO2_air,conv_mol,s,diss_f,k_H,k_w,k_CO2, ...
    n_H,n_CO2,SA,rock_fraction,M_calcite,MM_calcite,dtw];
if any(~isfinite(values)) || any(values([1:11,14:18]) < 0) || ...
        (H == 0 && n_H < 0) || (CO2_air == 0 && n_CO2 < 0) || ...
        K_CaCO3 == 0 || conv_mol == 0 || MM_calcite == 0 || dtw == 0
    value_args = num2cell(values);
    error('ERW:invalidCalciteKinetics', ...
        ['Invalid calcite kinetics: Ca=%g CO3=%g Ksp=%g H=%g CO2=%g ' ...
        'conv=%g s=%g diss=%g kH=%g kw=%g kCO2=%g nH=%g nCO2=%g ' ...
        'SA=%g fraction=%g mass=%g MM=%g dt=%g.'],value_args{:});
end

Omega = Ca*CO3/K_CaCO3;
affinity = max(0,1-Omega);
a_H = H/conv_mol;
pCO2 = CO2_air*22.41/conv_mol;
Wr = s*diss_f*(k_H*a_H^n_H + k_w + k_CO2*pCO2^n_CO2)*affinity;
area_fraction = SA*rock_fraction;
EW = Wr*area_fraction;

EW_max = M_calcite/MM_calcite/dtw;
if EW > EW_max
    EW = EW_max;
    if area_fraction > 0
        Wr = EW/area_fraction;
    else
        Wr = 0;
    end
end

end
