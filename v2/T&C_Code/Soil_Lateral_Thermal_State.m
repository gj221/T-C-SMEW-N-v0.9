function [lambda,capp,Oice,Oliq,cv] = Soil_Lateral_Thermal_State( ...
    T,Otot,Pre,p,OPT_FR_SOIL)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Thermal conductivity, apparent heat capacity and equilibrium soil ice.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isequal(size(T),size(Otot)) || any(~isfinite(T(:))) || ...
        any(T(:) <= -273.15) || any(~isfinite(Otot(:))) || ...
        any(Otot(:) < 0) || ~(isscalar(Pre) && isfinite(Pre) && Pre > 0)
    error('Soil_Lateral_Thermal_State:BadState', ...
        'Temperature, water content and pressure must be finite and physical.');
end
if ~isstruct(p) || ~isfield(p,'Soil_Param') || ...
        ~isfield(p,'SPAR') || ~isfield(p,'Phy')
    error('Soil_Lateral_Thermal_State:BadParameters', ...
        'The column parameters must include Soil_Param, SPAR and Phy.');
end

sp = p.Soil_Param;
if ~isfield(sp,'Ohy') || ~isfield(sp,'Osat') || ...
        ~(isscalar(sp.Ohy) || numel(sp.Ohy) == numel(Otot)) || ...
        ~(isscalar(sp.Osat) || numel(sp.Osat) == numel(Otot))
    error('Soil_Lateral_Thermal_State:BadParameters', ...
        'Residual water and porosity must match the thermal state.');
end
Ohy = sp.Ohy;
Osat = sp.Osat;
if ~isscalar(Ohy); Ohy = reshape(Ohy,size(Otot)); end
if ~isscalar(Osat); Osat = reshape(Osat,size(Otot)); end
if any(Otot(:) < Ohy(:)) || any(Otot(:) > Osat(:))
    error('Soil_Lateral_Thermal_State:BadWaterContent', ...
        'Total water must lie between residual content and porosity.');
end
[lambda,cv,~,Oice,Oliq,Cw] = Soil_Thermal_properties_FT( ...
    T,Pre,sp.rsd,sp.lan_dry,sp.lan_s,sp.cv_s,p.SPAR,sp.L,sp.Pe, ...
    sp.O33,sp.alpVG,sp.nVG,p.Phy,sp.s_SVG,sp.bVG,sp.Osat,sp.Ohy, ...
    Otot,OPT_FR_SOIL);

rho_i = 916.2;                                           %% kg/m3
Lf = 333700;                                             %% J/kg
g = 9.81;                                                %% m/s2
cv_w = 4186000;                                          %% J/m3/K
cv_i = rho_i*2117;                                       %% J/m3/K
T_state = reshape(T,size(Cw));
dOliq_dT = Lf*(1000*Cw)/(g*273.15);                      %% 1/K
dOliq_dT = dOliq_dT.*(OPT_FR_SOIL ~= 0);
capp = cv + (rho_i*Lf+T_state.*(cv_w-cv_i)).*dOliq_dT; %% J/m3/K
bad = ~isfinite(capp) | capp <= 0;
capp(bad) = cv(bad);
if any(~isfinite(lambda(:))) || any(lambda(:) < 0) || ...
        any(~isfinite(cv(:))) || any(cv(:) <= 0) || ...
        any(~isfinite(Oice(:))) || any(Oice(:) < 0) || ...
        any(~isfinite(Oliq(:))) || any(Oliq(:) < 0)
    error('Soil_Lateral_Thermal_State:BadProperties', ...
        'Soil thermal properties must be finite and physical.');
end
water_scale = max(1,max(abs(Otot(:))));
if max(abs(Oliq(:)+Oice(:)-Otot(:))) > 1e-12*water_scale
    error('Soil_Lateral_Thermal_State:WaterPhaseResidual', ...
        'Liquid and ice must sum to the supplied total water.');
end

end
