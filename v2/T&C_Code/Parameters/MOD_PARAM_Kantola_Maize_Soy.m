%%%%%%%%%%%%%%%%%%% PARAMETERS AND INITIAL CONDITION %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% adapted from previous case study: Zea mays  Glycine max 

%%% Maize–soybean rotations
%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SOIL AND HYDROLOGICAL PARAMETER
%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Rainfall Disaggregation
a_dis = NaN ;
pow_dis = NaN;
%%%%%%%%%%%%%%%%%%%%
fpr=1;
SvF=1; %% Sky View Factor
SN=0; %% Stream Identifier
Slo_top=0;  %% [fraction dy/dx]
Slo_pot=zeros(1,ms); %% [fraction dy/dx]
Asur = 1./cos(atan(Slo_top)); %% Real Area/Projected Area [m^2/m^2]
Ared = 1; %% 1-Frock
aR =1; %%% anisotropy ratio
%Kh=Ks*aR;
%%%%%%%%%%%%%%%%%%%%%%%%%%%
cellsize=1; %%[m^2];
aTop = 1000*cellsize^2./cellsize; %% [mm] Ratio betweeen Area/ContourLenght
Kbot = MOD_PARAM_input.Kbot; %% [mm/h] Conductivity at the bedrock layer
Krock = NaN; %% [mm/h] Conductivity of Fractured Rock
zatm = 8.0; %% Reference Height
%%%%%%%%%%%%%%%%%%
%%%%%%% VEG. SPECIES  --- Low Grasses -- High Decidous
%%%% LAND COVER PARTITION
Cwat = 0.0; Curb = 0.0 ; Crock = 0.0;
Cbare = MOD_PARAM_input.Cbare; 
% Ccrown = [MOD_PARAM_input.Ccrown MOD_PARAM_input.Ccrown];
eval(append('Ccrown','=',eval('MOD_PARAM_input.Ccrown'),';'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% SOIL INPUT %%%%  silty clay loam
Pcla= MOD_PARAM_input.Pcla; 
Psan= MOD_PARAM_input.Psan;
Porg= MOD_PARAM_input.Porg;
Color_Class = 0;
%%%%%%%%%%%%%%%%%%%
[Osat,L,Pe,Ks,O33,rsd,lan_dry,lan_s,cv_s,K_usle]=Soil_parameters(Psan,Pcla,Porg);
%%%%%%%%%%%%
rsd=rsd*ones(1,ms);
lan_dry=lan_dry*ones(1,ms);
lan_s =lan_s*ones(1,ms);
cv_s = cv_s*ones(1,ms);
%%%%%%%%%%%%%%%
SPAR=2; %%% SOIL PARAMETER TYPE 1-VanGenuchten 2-Saxton-Rawls
%nVG=L+1;
%alpVG = 1/(-101.9368*Pe); %%[1/mm]%;
p=3+2/L;
m=2/(p-1); nVG= 1/(1-m);
alpVG=(((-101.9368*Pe)*(2*p*(p-1))/(p+3))*((55.6+7.4*p+p^2)/(147.8+8.1*p+0.092*p^2)))^-1; %%[1/mm]%;
%%%
Osat=Osat*ones(1,ms);
%Ohy = Ohy*ones(1,ms) ; %% [-]
L=L*ones(1,ms);
Pe = Pe*ones(1,ms);
O33 = O33*ones(1,ms);
alpVG= alpVG*ones(1,ms); %% [1/mm]
nVG= nVG*ones(1,ms); %% [-]
Ks_Zs= Ks*ones(1,ms); %%[mm/h]
%%%%%%%%%%% Matric Potential
Kfc = 0.2; %% [mm/h]
Phy = 10000; %% [kPa]
%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
[Ofc,Oss_Lp,Owp_Lp,Ohy]=Soil_parametersII(ms,Osat,L,Pe,Ks_Zs,O33,nVG,alpVG,Kfc,1,1,Phy);
%%%%%%%%%%%%%%
clear Oss_Lp Owp_Lp
%%%%%%%%%%%%%%
Oice = 0;
%%%%%%%%%%%%%
%Zs = [0 10 50 100 150 200 300 400 600 800 1000 1250 1500 1750 2000 2500 3000 3500 4000 4500 5000]; %%% [ms+1]
eval(append('Zs','=',eval('MOD_PARAM_input.Zs'),';'))
Zdes = 10;
Zinf=  10;
Zbio = MOD_PARAM_input.Zbio;
if  not(length(Zs)==ms+1)
    disp('SOIL LAYER MESH INCONSISTENT')
    return
end
[EvL_Zs]=Evaporation_layers(Zs,Zdes); %%% Evaporation Layer fraction
[Inf_Zs]=Evaporation_layers(Zs,Zinf); %%% Infiltration Depth Layer fraction
[Bio_Zs]=Evaporation_layers(Zs,Zbio); %%% Infiltration Depth Layer fraction
dz= diff(Zs); %%%% [mm]  Thickness of the Layers
Dz=zeros(1,ms);
for i = 1:ms
    if i>1
        Dz(i)= (dz(i)+ dz(i-1))/2; %%% Delta Depth Between Middle Layer  [mm]
    else
        Dz(i)=dz(1)/2; %%% Delta Depth Between First Middle Layer and soil surface [mm]
    end
end
%%%%%%%%%%%%%%%%% OTHER PARAMETER
In_max_urb=5;
In_max_rock=0.1; %% [mm]
%%%%%%%%%%%%% SNOW PARAMETER
TminS=-0.8;%% Threshold temperature snow
TmaxS= 2.8;%% Threshold temperature snow
ros_max1=580; %600; %%% [kg/m^3]
ros_max2=300; %450; %%% [kg/m^3]
Th_Pr_sno = 8.0; %%% [mm/day] Threshold Intensity of snow to consider a New SnowFall
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% ICE Parameter
Ice_wc_sp =0.01; %% [-] Specific Maximum water content ice
ros_Ice_thr = 500 ; %% [kg/m^3] Density Thrshold to transform snow into ice
Aice = 0.28; %% [-] Ice albedo
WatFreez_Th = -8; %% [°C] Threshold for freezing lake water
dz_ice = 0.45; %% [mm / h] Water Freezing Layer progression without snow-layer
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%
ExEM = 0.0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% cc -- number of crown area
%%% Root Depth
CASE_ROOT=1;  %%% Type of Root Profile
%%%%%
ZR95_H = [0 0]; %% [mm]
ZR95_L = [920 700];%[900 1200]; %% [mm] ?
ZR50_H = [NaN NaN];
ZR50_L = [NaN NaN];
ZRmax_H = [NaN NaN];
ZRmax_L = [NaN NaN];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Kct=0.75; %%% Factor Vegetation Cover --- for throughfall
%5 Interception Parameter
gcI=3.7; %%% [1/mm]
KcI=0.06; %%%% [mm] -- Mahfouf and Jacquemin 1989
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Interception Parameter
Sp_SN_In= 5.9; %% [mm/LAI]
Sp_LAI_L_In= [0.2 0.2]; %%[mm/LAI]
Sp_LAI_H_In= [0.2 0.2]; %%[mm/LAI]
%%%%%%%%%%% Leaf Dimension
d_leaf_H= [2.0 2.0]; %%[cm]
d_leaf_L= [5.0 4.0];  %% [cm]
%%%%%%%% Biochemical parameter
KnitH=[0.2 0.2]; %%% Canopy Nitrogen Decay
KnitL=[0.15 0.15]; %%% Canopy Nitrogen Decay
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mSl_H = [0.0 0.0];%% [m2 PFT /gC]  Linear increase in Sla with LAI
mSl_L = [0.0 0.0]; % 0.001; %% [m2 LAI/gC * m2 PFT / m2 LAI]  0.0 - 0.004  Brod. Dec. Tree
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%  Photosynthesis Parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------
FI_H=[0.081 0.081];% Intrinsec quantum Efficiency [umolCO2/umolPhotons]
Do_H=[1000 1000]; %%[Pa]
a1_H=[7 7];
go_H=[0.01 0.01];% [mol / s m^2] minimum Stomatal Conductance
CT_H=[4 4]; %%--> 'CT' == 3  'CT' ==  4  %% Photosyntesis Typology for Plants
DSE_H =[0.649 0.649];  %% [kJ/mol] Activation Energy - Plant Dependent
Ha_H =[72 72]; %% [kJ / mol K]  entropy factor - Plant Dependent
gmes_H=[Inf Inf]; %% [mol CO2 / s m^2 ];  mesophyll conductance
rjv_H=[1.97 1.97]; %%% Scaling Jmax - Vmax  [umol electrons / umolCO2 ]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------
FI_L=[0.045 0.075];%[0.055 0.081];% Intrinsec quantum Efficiency [umolCO2/umolPhotons] ????
Do_L=[800 1000]; %%[Pa]
a1_L=[4 9];
go_L=[0.01 0.01];% % [mol / s m^2] minimum Stomatal Conductance
CT_L=[4 3];  %%--> 'CT' == 3  'CT' ==  4  %% Photosyntesis Typology for Plants
DSE_L = [0.649 0.649];%[0.649 0.649];%[0.649 0.649];  %% [kJ/mol] Activation Energy - Plant Dependent ?
Ha_L =[72 72];%%[72 72]; %% [kJ / mol K]  entropy factor - Plant Dependent???
gmes_L=[65535 65535];%[Inf Inf ]; ????
rjv_L= [2.4 2.4];
%%%%%%%%%%%%%%%maiz%%%%%%%%%%%%%%%
%Pss_H = [800]; Pwp_H = [3000]; %%% [kPa]
%Pss_L = [500]; Pwp_L = [3500]; %%% [kPa]
Psi_sto_50_H =  [-2.0 -2.0] ;%% [MPa]  Water Potential at 50% loss conductivity
Psi_sto_00_H =  [-0.5 -0.5]; %% [MPa]  Water Potential at 2% loss conductivity
%%% Leaf
PsiL00_H =  [-2.5 -2.5] ;%%[MPa]  Water Potential at 50% loss conductivity
PsiL50_H =  [-3.5 -3.5]; %% [MPa]  Water Potential at 2% loss conductivity
Kleaf_max_H = [ 5 5] ; %%  %%%  [mmolH20 m^2 leaf s /MPa]
Cl_H  = [1200 1200];  %%%  [500 - 3000]%  Leaf capacitance [mmolH20 / m^2 leaf MPa]
%%% Xylem
Axyl_H = [15.0 15.0] ; %% [cm^2 stem /m^2 PFT]
Kx_max_H = [80000 80000];  %%5550-555550 [mmolH20 /m s MPa]  Xylem Conductivity specific for water;
PsiX50_H = [-3.5 -3.5]; %%[MPa]  Water Potential at 50% loss conductivity
Cx_H= [150 150]; %%% [kg / m^3 sapwood MPa]
%%------------------------
%
Psi_sto_00_L = [-1.8 -1.5];%[-0.8 -0.5];%  %% [MPa]  Water Potential at PLCs loss conductivity
Psi_sto_50_L = [-2.8 -2.8];%  %% [MPa]  Water Potential at 50% loss conductivity
%%% Leaf
PsiL00_L =  [-1.2 -1.2]; %% [MPa]  Water Potential at PLCs% loss conductivity
PsiL50_L =  [-3.5 -3.5];%%[MPa]  Water Potential at 50% loss conductivity
Kleaf_max_L = [5 5] ; %%  %%%  [mmolH20 m^2 leaf s /MPa]
Cl_L  = [1200 1200];  %%%  [500 - 3000]%  [mmolH20 / m^2 leaf MPa]
%%% Xylem
Axyl_L = [0.0 0.0] ; %% [cm^2 stem /m^2 PFT]
Kx_max_L = [80000 80000];  %%5550-555550 [mmolH20 /m s MPa]  Xylem Conductivity specific for water;
PsiX50_L = [-4.5 -4.5]; %%[MPa]  Water Potential at 50% loss conductivity
Cx_L= [150 150]; %%% [kg / m^3 sapwood MPa]
%%%%%%%%%%%%%%%% Root Parameters
[RfH_Zs,RfL_Zs]=Root_Fraction_General(Zs,CASE_ROOT,ZR95_H,ZR50_H,ZR95_L,ZR50_L,ZRmax_H,ZRmax_L); 

%%%% Growth Parameters
PsiG50_H= [-0.45 -0.45];  %%[MPa]
PsiG99_H= [-1.2 -1.2];  %%[MPa]
gcoef_H = [3.5 3.5]; % [gC/m2 day]
%%------
PsiG50_L= [-1.2 -1.2];
PsiG99_L= [-3.5 -3.5];
gcoef_L = [3.5 3.5]; % [gC/m2 day]

%%%%%%%% Vegetation Optical Parameter
[PFT_opt_H(1)]=Veg_Optical_Parameter(0);
[PFT_opt_L(1)]=Veg_Optical_Parameter(16);
[PFT_opt_H(2)]=Veg_Optical_Parameter(0);
[PFT_opt_L(2)]=Veg_Optical_Parameter(16);
OM_H=[1 1];
OM_L=[1 1];
%%%%%%
Sllit = 2 ; %%% Litter Specific Leaf area [m2 Litter / kg DM]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% VEGETATION PART %%%%%%
%%% HIGH VEGETATION
%%%%%%%%%%%%%%%%%%%%%%%%%
Sl_H = [0.016 0.016]; % 0.018 0.05 -0.005 [m^2 gC] specific leaf area of  biomass [m^2 /gC]
Nl_H= [30 30]; %[gC/gN ] Leaf Carbon-Nitrogen ratio
[Stoich_H(1)]=Veg_Stoichiometric_Parameter(Nl_H(1));
[Stoich_H(2)]=Veg_Stoichiometric_Parameter(Nl_H(2));
%PLNR_H = 0.033; %% Percentage of Leaf N in Rubisco  [kgNRubisco / kgN]
r_H = [0.030 0.030];  %% [0.066 -0.011]respiration rate at 10° [gC/gN d ]
gR_H= [0.28 0.28]; % [0.22 - 0.28] growth respiration  [] -- [Rg/(GPP-Rm)]
%LAI_max_H= [2.0 10.0];
aSE_H= [1 1]; %%% Plant Type -- 1 Seasonal Plant --  0 Evergreen  -- 2 Grass species -- 3 Crops
dd_max_H= [1/365 1/365]; %%%0.005  [1/d]  0.0250 -- 0.005-0.025 death maximum for drought
dc_C_H =  [2/365 2/365]; %% [1/ d°C] -- [Factor of increasing mortality]
Tcold_H = [7 7]; %% [°C] Cold Leaf Shed
drn_H=  [1/1095 1/1095]; %% turnover root  [1/d]
dsn_H= [1/365 1/365]; % normal transfer rate sapwood [1/d]
age_cr_H= [150 150]; %% [day] Critical Leaf Age
Bfac_lo_H= [0.95 0.95]; %% Leaf Onset Water Stress
Bfac_ls_H= [NaN NaN]; %% Leaf Shed Water Stress [0-1]
Tlo_H = [12.9 12.9]; %% Mean Temperature for Leaf onset
Tls_H = [NaN NaN]; %% Mean Temperature for Leaf Shed
PAR_th_H= [NaN NaN]; 
dmg_H= [35 35]; %%% Tree 30 Grasses Day of Max Growth
LAI_min_H = [0.01 0.01];
Trr_H = [3.5 3.5]; %% Translocation rate [gC /m^2 d]
mjDay_H = [180 180]; %% Maximum Julian day for leaf onset
LDay_min_H =[12.58 12.58]; %% Minimum Day duration for leaf onset
LtR_H = [0.8 0.8]; %%% Leaf to Root ratio maximum
Mf_H= [1/50 1/50]; %% fruit maturation turnover [1/d]
Wm_H= [0 0 ] ; % wood turnover coefficient [1/d]
eps_ac_H = [1 1]; %% Allocation to reserve parameter [0-1]
LDay_cr_H = [12.30 12.3]; %%%  Threshold for senescence day light [h]
Klf_H =[1/15 1/15]; %% Dead Leaves fall turnover [1/d]
fab_H = [0.74 0.74]; %% fraction above-ground sapwood and reserve
fbe_H = [0.26 0.26]; %% fraction below-ground sapwood and reserve
ff_r_H= [0.1 0.1]; %% Reference allocation to Fruit and reproduction
soCrop_H = [NaN NaN]; 
Sl_emecrop_H = [NaN NaN];
MHcrop_H =[NaN NaN]; %% Maximum height crop  
[ParEx_H(1)]=Exudation_Parameter(0); 
[ParEx_H(2)]=Exudation_Parameter(0); 
[Mpar_H(1)]=Vegetation_Management_Parameter;
[Mpar_H(2)]=Vegetation_Management_Parameter;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% LOW VEGETATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sl_L = [0.034 0.033];%[0.034 0.033]; % 0.016-0.020; %%  0.05 -0.005 [m^2 gC] specific leaf area of  biomass [m^2 /gC]
Nl_L= [15 12]; %[gC/gN ] Leaf Carbon-Nitrogen ratio
[Stoich_L(1)]=Veg_Stoichiometric_Parameter(Nl_L(1),"Maize");
[Stoich_L(2)]=Veg_Stoichiometric_Parameter(Nl_L(2),"Soybean");
%PLNR_L = 0.033; %% Percentage of Leaf N in Rubisco  [kgNRubisco / kgN]
r_L = [0.066 0.066];  %[0.025 0.025];  %% [0.066 -0.011]respiration rate at 10° [gC/gN d ]
gR_L= [0.25 0.25]; % [0.22 - 0.28] growth respiration  [] -- [Rg/(GPP-Rm)]
%LAI_max_L= [2.0 10.0];2
aSE_L= [5 5]; %%% Plant Type -- 1 Seasonal Plant --  0 Evergreen  -- 2 Grass species - ?
dd_max_L= [1/50 1/50];%%%0.005  [1/d]  0.0250 -- 0.005-0.025 death maximum for drought
dc_C_L = [7/365 7/365]; %% [1/ d°C] -- [Factor of increasing mortality]
Tcold_L = [-5 -5];%[6.0 6.0]; %% [°C] Cold Leaf SLed ?
drn_L=  [1/365 1/365]; %% turnover root  [1/d]
dsn_L= [0 0]; % normal transfer rate sapwood [1/d]
age_cr_L= [80 100];%[80 80];%[70 75]; %% [day] Critical Leaf Age ?
Bfac_lo_L= [0.95 0.95]; %[0.99 0.99]; %% Leaf Onset Water Stress
Bfac_ls_L= [NaN NaN]; %% Leaf Shed Water Stress [0-1]
Tlo_L = [15 15];%[7.0 7.0]; %% Mean Temperature for Leaf onset
Tls_L = [7 7]; %% Mean Temperature for Leaf Shed
PAR_th_L= [NaN NaN]; %%% Stem allocation for Crops  
dmg_L= [100 100];%[45 40]; %%% Tree 30 Grasses Day of Max Growth ?
LAI_min_L = [0.1 0.1];
Trr_L = [4.5 1.0];%[4.5 6.5]; %% Translocation rate [gC /m^2 d] ?
mjDay_L = [300 280];%[200 200]; %% Maximum Julian day for leaf onset ?
LDay_min_L = [12.0 12.0];%[13.0 13.0]; %% Minimum Day duration for leaf onset ?
LtR_L = [1.2 1.4];%[2.5 2.4];%[1.2 1.4]; %%% Leaf to Root ratio maximum ?
Mf_L= [0 0]; %% fruit maturation turnover [1/d]
Wm_L= [0 0] ; % wood turnover coefficient [1/d]
eps_ac_L = [0.2 0.2]; %% Allocation to reserve parameter [0-1]
LDay_cr_L = [9 9];%[13.4 13.4]; %%%  Threshold for senescence day light [h] ?
Klf_L =  [1/20 1/20] ; % [1/83]; %% Dead Leaves fall turnover [1/d]
fab_L = [1.0 1.0]; %% fraction above-ground sapwood and reserve ?
fbe_L =[0.0 0.0]; %% fraction below-ground sapwood and reserve
ff_r_L= [0.50 0.20];
soCrop_L = [0.05 0.20]; 
Sl_emecrop_L = [0.03 0.001];%[0.015 0.015]; %%% Additional SLA at emergence ?
MHcrop_L =[2.5 1.0]; %%[m] maximum crop height  
[ParEx_L(1)]=Exudation_Parameter(0); 
% [ParEx_L(2)]=Exudation_Parameter(0); 
[ParEx_L(2)]=Exudation_Parameter(1); % soybean N fixation
[Mpar_L(1)]=Vegetation_Management_Parameter;
[Mpar_L(2)]=Vegetation_Management_Parameter;

%%%%%%%%%%%%%%%
% FROM MEAD_MAIZE_SOY.m
% pretreatment: 2009-2016 maize-soy-maize-maize-soy-maize-maize
% Kantola rotation: -soy-maize-maize-soy-maize for 2016-2020
% spring cultivation operations in April
% maize
Mpar_L(1).Date_sowing = datenum(2009:1:2020,5,15,10,0,0) ; %% Date of Sowing
if NN_case_study_ERW_MOD == 21 % pretreatment
    Mpar_L(1).Date_harvesting = datenum(2009:1:2020,08,27,10,0,0) ; %%% Fully Harvested
else
    Mpar_L(1).Date_harvesting = datenum(2009:1:2020,09,17,10,0,0) ; %%% Fully Harvested
end
Mpar_L(1).Crop_B=[20 30];
Mpar_L(1).Crop_crown =[1 0 1 1 0 1 1 0 1 1 0 1]; 

% soy
Mpar_L(2).Date_sowing = datenum(2009:1:2020,5,25,10,0,0) ; %% Date of Sowing
Mpar_L(2).Date_harvesting = datenum(2009:1:2020,10,10,10,0,0); %%% Fully Harvested
Mpar_L(2).Date_harvesting = datenum(2009:1:2020,9,10,17,0,0); %%% Fully Harvested
Mpar_L(2).Crop_B=[20 30];
Mpar_L(2).Crop_crown =[0 1 0 0 1 0 0 1 0 0 1 0]; 

%%%%%%%%%%%%%
%%%%
Mpar_L(1).fract_left=1;
Mpar_L(1).fract_left_fr=0;
Mpar_L(1).fract_left_AB=0.5; 
Mpar_L(1).fract_left_BG=1;

Mpar_L(2).fract_left=1;
Mpar_L(2).fract_left_fr=0;
Mpar_L(2).fract_left_AB=1;
Mpar_L(2).fract_left_BG=1; 
%%%%%%

% from MOD_PARAM_MAIZE.m
Mpar_L(1).NPK_res_ini = [20 3 20 0.9 1.5 7]; 
Mpar_L(2).NPK_res_ini = [15 3 20 5.8 2   5.5];

%%%%%%%%%%%% PRODUCTIVITY
Vmax_H = [0 0]; 
Vmax_L = [78 88];
%%%%%%%%%%%%%%%%%%%%%%
%%%%
L_day=zeros(NNd,1);
for j=2:24:NN
    [h_S,delta_S,zeta_S,T_sunrise,T_sunset,L_day(j)]= SetSunVariables(Datam(j,:),DeltaGMT,Lon,Lat,t_bef,t_aft);
end
Lmax_day = max(L_day);
clear('h_S','delta_S','zeta_S','T_sunrise','T_sunset','L_day')
%%%%%%
%%%%%%%%% DORMANT 1 - MAX GROWTH 2 - NORMAL GROWTH 3 - SENESCENCE 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ci_sunL(1,:) = [Ca(1)]; % %% [umolCO2/mol]
Ci_sunH(1,:) = [Ca(1)]; %% [umolCO2/mol]
Ci_shdL(1,:) = [Ca(1)]; % %% [umolCO2/mol]
Ci_shdH(1,:) = [Ca(1)]; %% [umolCO2/mol]
%%%%%%%
LAI_H(1,:)=[0.0 0.0]; %
B_H(1,:,:)= [0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0]; %%
Rrootl_H(1,:) = [0 0] ;
PHE_S_H(1,:)=[0 0];
dflo_H(1,:)=[0 0];
AgeL_H(1,:)=[0 0];
e_rel_H(1,:)=[0 0];
hc_H(1,:) =[0 0]; %%
SAI_H(1,:) = [0 0]; %%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1 Live Leaves/ 2 Sapwood/ 3 Fine Roots / 4 Carbohydrate Reserve /5 Fruit and Flower /6 Heartwood - Dead Sapwood / 7 Standing Dead Leaves
LAI_L(1,:)=[0.0 0.0]; %
B_L(1,:,:)= [0 0 0 0 0 0 0 0 ; 0 0 0 0 0 0 0 0]; %%    
Rrootl_L(1,:) = [0 0] ;
% Rrootl_L(1,:) = [4250 4250] ;

PHE_S_L(1,:)=[1 1]; %[3];
dflo_L(1,:)=[0 0];
AgeL_L(1,:)=[0 0] ;%[730];
e_rel_L(1,:)=[1 1];
hc_L(1,:) =[2 0.2]; %%
SAI_L(1,:) = [0.001 0.001]; %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BLit(1)=0.0;  %% [kg DM /m2 PFT] Litter Biomass
%%%%%%%%%%
Nreserve_H(1,:)= [0 0];
Preserve_H(1,:)= [0 0];
Kreserve_H(1,:)= [0 0];
Careserve_H(1,:)= [1000 1000];
Mgreserve_H(1,:)= [1000 1000];
Sireserve_H(1,:)= [1000 1000];
FNC_H(1,:)=[1 1];
NupI_H(1,:,:)= [0 0 0 0 0 0; 0 0 0 0 0 0];
Nreserve_L(1,:)= [1000 1000];
Preserve_L(1,:)= [100 100];
Kreserve_L(1,:)= [1000 1000];
FNC_L(1,:)=[1 1];
NupI_L(1,:,:)= [0 0 0 0 0 0; 0 0 0 0 0 0];
RexmyI(1,:)= [0 0 0]; 
Careserve_L(1,:)= [1000 1000];
Mgreserve_L(1,:)= [1000 1000];
Sireserve_L(1,:)= [1000 1000];
%%%%%%%%%

if OPT_SoilBiogeochemistry == 1

    % p pools and relevant values from spin-up
    if exist("B_out",'var')
        % outputs from BG spinup 
        P(1,:)=  B_out;

        if exist("RexmyI_out",'var')
            RexmyI(1,:)= RexmyI_out;
            NavlI(1,:)=NavlI_out;
            NupI_L(1,:,:) = NupI_L_out;

            O(1,:) = O_out;

        else

            FNC_L(1,:)=1;
        
            NupI_L(1,:,:)= [ 0.100826 0.01101  0.043436947 0.043436947 0.043436947 0.043436947
                0.142 0.0154 0.0702 0.0702 0.0702 0.0702];% from MAIZE.m
            RexmyI(1,:,:) = [ 0.033995 0.1972202  0];% from MAIZE.m
            NavlI(1,:,:) =[  0.244   0.0115  0.0618 0.0618 0.0618 0.0618];% from MAIZE

        end

    else
        Nreserve_L(1,:)= 0 ; %
        Preserve_L(1,:)= 0; %
        Kreserve_L(1,:)= 0 ; %
        FNC_L(1,:)=1;
    
        % NupI_L(1,:,:)= [ 0.144 0.0192 0.0802 0.0802 0.0802 0.0802
        %     0.142 0.0154 0.0702 0.0702 0.0702 0.0702];
        NupI_L(1,:,:)= [ 0.100826 0.01101  0.043436947 0.043436947 0.043436947 0.043436947
            0.142 0.0154 0.0702 0.0702 0.0702 0.0702];% from MAIZE.m
    
        % RexmyI(1,:)= [ 0.002161627596424   0.031874606960811     0];
        RexmyI(1,:,:) = [ 0.033995 0.1972202  0];% from MAIZE.m
    
        % NavlI(1,:)=[10 10 10 10 10 10];
        NavlI(1,:,:) =[  0.244   0.0115  0.0618 0.0618 0.0618 0.0618];% from MAIZE
    
        %%%%%%
    
       P(1,:)= 1000*[     0.000012410591140   0.000800526607103   0.000141269401254  -0.000000000001652  -0.000000000001652   0.000000776268590,...
           0.006703898603793   0.000744877622644   0.226606508676659   0.365236699866817   2.755486852392532   0.013357831027381,...
           0.011253189397630   0.000031750405418   0.000020301337334   0.000013229335591   0.000033835562223   0.005008342341374,...
           0.017281420404693   0.004900096307971  -0.000000000064899   0.000925141432637   0.000012546771363  -0.000000000001652,...
           0.000097031039541   0.551558213756260   0.001000680893127   0.002671209613892   0.000272480142860  -0.000000000031930,...
           0.002889844814344   0.000244299889924   0.000004953828542   0.000092514143263   0.000008418169974  -0.000000000001652,...
           0.000025623189046   0.230170746027678   0.000315926281660   0.000437368108832   0.000040872021429  -0.000000000031930,...
           0.119115782977356   0.150000000000000   2.510785111813128   0.020616211158395   0.000000137302977   0.000021470098324,...
          -0.000000000001652   0.000017481042041   0.008820644116622   0.002690221023989   0.002686977865475   0.134670712696155,...
           0.501030497641396,...
          ...
          0.001598511730988  -0.000000000006909 0.000150160488437   0.009450287825968   0.000048677659335   0.000034578403829   0.000778005941413   0.502217651617589,...
          0.001598511730988  -0.000000000006909 0.000150160488437   0.009450287825968   0.000048677659335   0.000034578403829   0.000778005941413   0.502217651617589,...
          0.001598511730988  -0.000000000006909 0.000150160488437   0.009450287825968   0.000048677659335  0.000034578403829   0.000778005941413   0.502217651617589,...
             ];
        
      
    end

    % Fertiliser application annual schedule
    FertN=zeros(1,366);
    FertN(81:100) = 0.84;
    FertP=zeros(1,366);
    FertK=zeros(1,366);
    FertCa=zeros(1,366);
    FertMg=zeros(1,366);
    FertSi=zeros(1,366);

    Upl=0.01;
    HIST=0;

    [B_IO]=Biogeochemistry_IO(Zs(ms+1),Lat,Lon,Upl,HIST,FertN,FertP,FertK,FertCa,FertMg,FertSi);
    B_IO.SC_par=[1 1 1 1];
    
    B_IO.N_Man=1e-6;
    B_IO.P_Man=1e-6;
    B_IO.K_Man=1e-6;
    B_IO.Ca_Man=1e-6;
    B_IO.Mg_Man=1e-6;
    B_IO.Si_Man=1e-6;
    
    B_IO.Tup_P = 0;
    B_IO.Tup_Ca = 0.003;
    B_IO.Tup_Mg = 0.002;
    B_IO.Tup_Si = 0.0001;

end
PH(1) = MOD_PARAM_input.PH;
%%%
%%%
TBio_L=[1 1];  %%[ton DM / ha ]
TBio_H=[0 0];  %[ton DM / ha ]
%%%%%%%%%%%%%%%%%
Vx_H=[0 0];  %% [mm/ m2 PFT];
Vl_H=[0 0];  %% [mm/ m2 PFT];
Vx_L=[0 0];   %% [mm/ m2 PFT];
% Vl_L=[1 1];   %% [mm/ m2 PFT];
Vl_L=[100 1]; % from MAIZE
%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%% Initial Conditions
%%%%%%%%%%%%%%%%%%%% Initial Conditions
SWE(1)=0; %% [mm]
SND(1)=0;
Ts(1) = Ta(1)-2;
Tdamp(1) =Ta(1); 
Tdp(1,:)= Ta(1)*ones(1,ms);
Tdp_H(1,:)=Tdp(1,1);
Tdp_L(1,:)=Tdp(1,1);
TdpI_H(1,:)=Tdp(1,1);
TdpI_L(1,:)=Tdp(1,1);

%%% Snow_alb = soil_alb initial
snow_alb.dir_vis = 0.2;
snow_alb.dif_vis = 0.2;
snow_alb.dir_nir = 0.2;
snow_alb.dif_nir = 0.2;
In_L(1,:)=0; In_H(1,:)=0;
In_urb(1)=0; In_rock(1)= 0;
In_Litter(1)=0;
SP_wc(1)=0 ; %%[mm]
In_SWE(1)= 0;
ros(1)= 0;
t_sls(1)= 0;
e_sno(1) = 0.97;
tau_sno(1) = 0;
EK(1)=0;
WAT(1) = 0;
ICE(1) = 0;
IP_wc(1)=0;
ICE_D(1)= 0;
FROCK(1)=0;
Ws_under(1)=1;
%%%%%%%%%%%%%% Volume [mm]
if sum(O(1,:))==0
    O(1,:)= [      0.1670    0.1782    0.2320    0.2386    0.2369    0.2646    0.3003    0.3223    0.3290    0.3296    0.3295    0.3303    0.3315    0.3326,...
          0.3341    0.3356    0.3369    0.3383    0.3395    0.3401]; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V(1,:) = (O(1,:) -Ohy).*dz;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%
