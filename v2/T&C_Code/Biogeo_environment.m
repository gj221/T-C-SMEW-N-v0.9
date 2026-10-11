%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Temperature Humidity Biogeochemistry Zone
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[Se_bio,Se_fc,Psi_bio,Tdp_bio,Vtot,VT,f_anaer_bio,WFPS_bio]=Biogeo_environment(Tdp,O,V,Soil_Param,Phy,SPAR,Bio_Zs,Zbio)
O=mean(O,1); Tdp=mean(Tdp,1);  V=mean(V,1); 
%%%%%%%
Ks_Zs=Soil_Param.Ks_Zs; 
Osat=Soil_Param.Osat;
Ohy =Soil_Param.Ohy; 
L=Soil_Param.L;
Pe=Soil_Param.Pe; 
O33=Soil_Param.O33; 
alpVG=Soil_Param.alpVG; 
nVG=Soil_Param.nVG; 
lVG=Soil_Param.lVG; 
Ofc=Soil_Param.Ofc; 
Ks_mac=Soil_Param.Ks_mac; 
Omac=Soil_Param.Omac; 
alpVGM=Soil_Param.alpVGM;
nVGM =Soil_Param.nVGM;
lVGM =Soil_Param.lVGM;
s_SVG =Soil_Param.s_SVG;
bVG =Soil_Param.bVG;
%%%%%%%%
%%%% Biogeochemistry Layer Water Content
Obio = sum(Bio_Zs.*O);
Ohy = sum(Bio_Zs.*Ohy); 
Osat = sum(Bio_Zs.*Osat); 
Ofc = sum(Bio_Zs.*Ofc); 
WFPS_bio = Biogeo_WFPS(Obio,Osat); %% [-] theta/theta_sat over the active layer
Obio(Obio<=Ohy)= Ohy + 1e-5; 
Obio(Obio>=Osat)= Osat - 1e-5; 
Se_bio = (Obio-Ohy)./(Osat-Ohy);%%
Se_fc = (Ofc-Ohy)./(Osat-Ohy);%%
[Ko,Psi_bio]=Conductivity_Suction(SPAR,sum(Bio_Zs.*Ks_Zs),Osat,Ohy,...
    sum(Bio_Zs.*L),sum(Bio_Zs.*Pe),sum(Bio_Zs.*O33),sum(Bio_Zs.*alpVG),sum(Bio_Zs.*nVG),sum(Bio_Zs.*lVG),...
     sum(Bio_Zs.*Ks_mac),sum(Bio_Zs.*Omac), sum(Bio_Zs.*alpVGM),sum(Bio_Zs.*nVGM),sum(Bio_Zs.*lVGM),Phy,sum(Bio_Zs.*s_SVG),sum(Bio_Zs.*bVG),Obio);
%%% Anaerobic fraction for denitrification: air-entry suction over the
%%% current suction, Pe/|Psi|, both in kPa. Conductivity_Suction returns
%%% Psi_bio in mm of head, so convert with rho*g = 9.81 kPa per m; Pe is
%%% already kPa (Soil_parameters). This is the quantity Dickinson et al.
%%% (2002) mean by Se^(1/lambda) under Brooks-Corey, computed here from
%%% whatever retention curve SPAR actually selects instead of assuming its
%%% shape: under SPAR = 2 (Saxton-Rawls) the curve is linear in O above O33
%%% and Se^(1/L) is NOT the suction ratio, L there being a log-log slope
%%% rather than a pore-size index. The ratio is 1 at saturation (Psi -> Pe
%%% on the Saxton-Rawls wet branch) and falls as the profile drains.
%%% Formed from the UNCLAMPED suction: BIOGEOCHEMISTRY_DYNAMIC3 clamps its
%%% own Psi_s to -10.1 kPa for the Moyano/Manzoni respiration functions,
%%% which would cap this ratio at Pe/10 and never let it reach 1.
Pe_bio = sum(Bio_Zs.*Pe); %% [kPa] air-entry suction over the active layer
Psi_kPa = max(abs(Psi_bio)*9.81/1000, eps); %% [kPa] current suction, unclamped
f_anaer_bio = min(1, max(0, Pe_bio/Psi_kPa)); %% [-] anaerobic fraction
Tdp_bio = sum(Bio_Zs.*Tdp); %% [�C]
Psi_bio=-(Psi_bio/1000)*1000*9.81/1e+6; %%[MPa]
Vtot = sum(V); %% [mm] Volume of entire depth
%%% Total water standing in the biogeochemically active layer, which is
%%% what the passive uptake term needs (supplementary S15: "VT is the total
%%% water volume within the biogeochemically active layer").
%%% sum(Bio_Zs.*V) does not compute that. Bio_Zs are weights normalised to
%%% sum to one, so with V_i = (O_i-Ohy)*dz_i it evaluates to
%%% sum(dz_i^2*(O_i-Ohy))/Zbio, which changes if the same depth is sliced
%%% into different layers: 12.1 mm on the current [10 40 50 50 50 100] mesh,
%%% 9.4 mm at 6x50, 1.9 mm at 30x10. The physical answer is Zbio*(Obio-Ohy)
%%% = 56.5 mm on all three. Obio and Ohy are already the depth-weighted
%%% values computed above, so this is mesh independent by construction.
VT = Biogeo_Mobile_Water_Volume(Zbio,Obio,Ohy); %% [mm] Water volume in the Biogeochemistry layer
end
