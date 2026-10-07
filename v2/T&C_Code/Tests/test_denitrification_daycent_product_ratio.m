function test_denitrification_daycent_product_ratio()
% Verify the DayCent controls used only to partition gross denitrification.

pH = 6.0;
WFPS = 0.80;
NO3_ppm = 50.0;
respiration_kg_C_ha_day = 20.0;

fpH = 0.0016*exp(1.006*pH);
fNO3 = (1.0 - (0.5 + atan(pi*0.01*(NO3_ppm - 190.0))/pi))*25.0;
fResp = 13.0 + 30.78*atan(pi*0.07*(respiration_kg_C_ha_day - 13.0))/pi;
pDEN = 13.0^(2.2*WFPS);
fWFPS = 1.4/13.0^(17.0/pDEN);
expected = max(0.1,min(fNO3,fResp)*fWFPS*fpH);

actual = Biogeo_Denitrification_Product_Ratio( ...
    fpH,WFPS,NO3_ppm,respiration_kg_C_ha_day);
assert(abs(actual - expected) < 1e-12, ...
    'DayCent N2:N2O product ratio does not reproduce the source equations.');

low_NO3 = Biogeo_Denitrification_Product_Ratio(fpH,WFPS,10.0,respiration_kg_C_ha_day);
high_NO3 = Biogeo_Denitrification_Product_Ratio(fpH,WFPS,500.0,respiration_kg_C_ha_day);
assert(high_NO3 < low_NO3, ...
    'Increasing NO3 must lower the DayCent N2:N2O ratio.');

dry = Biogeo_Denitrification_Product_Ratio(fpH,0.30,NO3_ppm,respiration_kg_C_ha_day);
wet = Biogeo_Denitrification_Product_Ratio(fpH,0.95,NO3_ppm,respiration_kg_C_ha_day);
assert(wet > dry, ...
    'Increasing WFPS must raise the DayCent N2:N2O ratio.');

floored = Biogeo_Denitrification_Product_Ratio(fpH,0.05,NO3_ppm, ...
    respiration_kg_C_ha_day);
pDEN_dry = 13.0^(2.2*0.05);
unfloored = min(fNO3,fResp)*(1.4/13.0^(17.0/pDEN_dry))*fpH;
assert(unfloored < 0.1 && floored == 0.1, ...
    'A dry cell must return the Blanc-Betes Eq. 6 floor, not pure N2O.');

gross_denitrification = 0.73;
N2O = gross_denitrification/(1.0 + actual);
N2 = gross_denitrification - N2O;
assert(abs(N2 + N2O - gross_denitrification) < 1e-14, ...
    'The product split must conserve gross denitrified N.');

invalid_inputs = {[-0.1 WFPS NO3_ppm respiration_kg_C_ha_day], ...
    [fpH -0.1 NO3_ppm respiration_kg_C_ha_day], ...
    [fpH 1.1 NO3_ppm respiration_kg_C_ha_day], ...
    [fpH WFPS -1.0 respiration_kg_C_ha_day], ...
    [fpH WFPS NO3_ppm -1.0]};
for i = 1:numel(invalid_inputs)
    values = num2cell(invalid_inputs{i});
    failed = false;
    try
        Biogeo_Denitrification_Product_Ratio(values{:});
    catch
        failed = true;
    end
    assert(failed,'Invalid DayCent product-ratio input was accepted.');
end

fprintf('test_denitrification_daycent_product_ratio passed.\n');
end
