%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Temperature Humidity Biogeochemistry Zone
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[Se_bio,Se_fc,Psi_bio,Tdp_bio,Vtot,VT,L_bio]=Biogeo_environment(Tdp,O,V,Soil_Param,Phy,SPAR,Bio_Zs,Zbio)
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
Obio(Obio<=Ohy)= Ohy + 1e-5; 
Obio(Obio>=Osat)= Osat - 1e-5; 
Se_bio = (Obio-Ohy)./(Osat-Ohy);%%
Se_fc = (Ofc-Ohy)./(Osat-Ohy);%%
[Ko,Psi_bio]=Conductivity_Suction(SPAR,sum(Bio_Zs.*Ks_Zs),Osat,Ohy,...
    sum(Bio_Zs.*L),sum(Bio_Zs.*Pe),sum(Bio_Zs.*O33),sum(Bio_Zs.*alpVG),sum(Bio_Zs.*nVG),sum(Bio_Zs.*lVG),...
     sum(Bio_Zs.*Ks_mac),sum(Bio_Zs.*Omac), sum(Bio_Zs.*alpVGM),sum(Bio_Zs.*nVGM),sum(Bio_Zs.*lVGM),Phy,sum(Bio_Zs.*s_SVG),sum(Bio_Zs.*bVG),Obio);
L_bio = sum(Bio_Zs.*L); %% [-] pore size distribution index over the same depth
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
VT = Zbio*(Obio-Ohy); %% [mm] Water volume in the Biogeochemistry layer
end
