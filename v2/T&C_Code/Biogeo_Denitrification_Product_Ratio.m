function RN2N2O = Biogeo_Denitrification_Product_Ratio(fpH,WFPS,NO3_ppm,Resp_kg_C_ha_day)
% DayCent N2:N2O ratio applied to the existing gross denitrification flux.

if ~(isnumeric(fpH) && isreal(fpH) && isscalar(fpH) && isfinite(fpH) && fpH >= 0)
    error('Biogeo:InvalidPHResponse', ...
        'fpH must be a finite nonnegative scalar.');
end
if ~(isnumeric(WFPS) && isreal(WFPS) && isscalar(WFPS) && isfinite(WFPS) && WFPS >= 0 && WFPS <= 1)
    error('Biogeo:InvalidWFPS','WFPS must be a finite scalar between 0 and 1.');
end
if ~(isnumeric(NO3_ppm) && isreal(NO3_ppm) && isscalar(NO3_ppm) && isfinite(NO3_ppm) && NO3_ppm >= 0)
    error('Biogeo:InvalidNO3','NO3_ppm must be a finite nonnegative scalar.');
end
if ~(isnumeric(Resp_kg_C_ha_day) && isreal(Resp_kg_C_ha_day) && isscalar(Resp_kg_C_ha_day) && isfinite(Resp_kg_C_ha_day) && Resp_kg_C_ha_day >= 0)
    error('Biogeo:InvalidRespiration', ...
        'Resp_kg_C_ha_day must be a finite nonnegative scalar.');
end

fNO3 = (1 - (0.5 + atan(pi*0.01*(NO3_ppm - 190))/pi))*25;
fResp = 13 + 30.78*atan(pi*0.07*(Resp_kg_C_ha_day - 13))/pi;
fWFPS = Biogeo_Denitrification_WFPS_Response(WFPS);
%%% Blanc-Betes et al. (2021), Eq. 6, floors the complete N2:N2O product at
%%% 0.1, not any single factor. Flooring fpH alone held the pH term up below
%%% pH 4.1 and still let a dry or nitrate-rich cell drive the product to
%%% zero, which returned denitrification as pure N2O. The nitrate,
%%% respiration and WFPS terms remain the legacy TC-SMEW-DayCent forms
%%% (Parton et al., 1996; Del Grosso et al., 2000), not the Blanc-Betes
%%% Eqs. 7 and 8 substitutions.
RN2N2O = max(0.1,min(fNO3,fResp)*fWFPS*fpH);

if ~(isfinite(RN2N2O) && RN2N2O >= 0)
    error('Biogeo:InvalidProductRatio', ...
        'The DayCent N2:N2O product ratio must be finite and nonnegative.');
end
end
