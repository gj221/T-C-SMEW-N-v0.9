%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MOD_PARAM generalisation to use parameters for specified PFT(s) 
% default values for other parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loading example PFT parameters
PFT_info_opts = detectImportOptions("Table_Sites_PFT_LOOKUP_latest.xlsx");
PFT_info_opts.DataRange = "A2:BS260";
PFT_info_opts.VariableTypes = ["double","string","string","double","double",...
    "double","string","string","string",repmat("double",1,62)];

PFT_info = readtable("Table_Sites_PFT_LOOKUP_latest.xlsx",PFT_info_opts);
par_names = PFT_info.Properties.VariableNames(1,10:end);

%% INPUTS (defined by users)
% %%%%%%%%%% PFTs selection 
% PFT_type = [1 2]; % row number from excel database
% stature_type = ["H" "H"]; % H - high; L - low

eval(append('PFT_row','=',eval('MOD_PARAM_input.PFT_type'),';'))
eval(append('stature_type','=',eval('MOD_PARAM_input.stature_type'),';'))

%% Load PFT parameters based on selections
% prepare selected PFT parameters
PFT_info = PFT_info(ismember(PFT_info.Number,PFT_row),:);

MOD_PARAM_PFT_pars = PFT_selection_parameters(PFT_info, cc, stature_type,id_location);
fields = fieldnames(MOD_PARAM_PFT_pars);
for ind_vars =1:length(fields)
    eval([fields{ind_vars},'=MOD_PARAM_PFT_pars.',fields{ind_vars},';'])
end

% other veg relevant parameters
CASE_ROOT= 1 ;  %%% Type of Root Profile 
Kct=0.75; %%% Factor Vegetation Cover --- for throughfall
gcI=3.7; %%% [1/mm]
KcI=0.06; %%%% [mm] -- Mahfouf and Jacquemin 1989
%%%%%%%%% Interception Parameter
Sp_SN_In= 5.9; %% [mm/LAI]
Sllit = 2 ; %%% Litter Specific Leaf area [m2 Litter / kg DM]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if  not(exist('PFT_opt_H','var'))
    for len_veg = 1:cc
        [PFT_opt_H(len_veg)]=Veg_Optical_Parameter(0); 
    end 
end

if  not(exist('PFT_opt_L','var'))
    for len_veg = 1:cc
        [PFT_opt_L(len_veg)]=Veg_Optical_Parameter(0); 
    end 
end

%% Other Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
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
Ared = 1; 
aR =1; %%% anisotropy ratio
%%%%%%%%%%%%%%%%%%%%%%%%%%%
cellsize=1; %%[m^2];
aTop = 1000*cellsize^2./cellsize; %% [mm] Ratio betweeen Area/ContourLenght
Kbot =  MOD_PARAM_input.Kbot; %% [mm/h] Conductivity at the bedrock layer
Krock = NaN; %% [mm/h] Conductivity of Fractured Rock 
zatm = 5; %% Reference Height
%%%%%%%%%%%%%%%%%%
%%%% LAND COVER PARTITION
Cwat = 0.0; Curb = 0.0 ; Crock = 0.0;
Cbare = MOD_PARAM_input.Cbare; 
% Ccrown = MOD_PARAM_input.Ccrown;
eval(append('Ccrown','=',eval('MOD_PARAM_input.Ccrown'),';'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% SOIL INPUT %%% Loamy sand
Pcla= MOD_PARAM_input.Pcla; 
Psan= MOD_PARAM_input.Psan;
Porg= MOD_PARAM_input.Porg; 
Color_Class = 0;
%%%%%%%%%%%%%%%%%%%
[Osat,L,Pe,Ks,O33,rsd,lan_dry,lan_s,cv_s,K_usle]=Soil_parameters(Psan,Pcla,Porg);

% compare rsd value
if abs((rsd-MOD_PARAM_input.rsd)/rsd)>1e-2
    disp(rsd)
    % error('check rsd value')
    rsd = MOD_PARAM_input.rsd;
    
end

%%%%%%%%%%%%%%%
rsd=rsd*ones(1,ms); 
lan_dry=lan_dry*ones(1,ms); 
lan_s =lan_s*ones(1,ms); 
cv_s = cv_s*ones(1,ms);
%%%%
SPAR=2; %%% SOIL PARAMETER TYPE 1-VanGenuchten 2-Saxton-Rawls
%nVG=L+1;
%alpVG = 1/(-101.9368*Pe); %%[1/mm]%;
p=3+2/L;
m=2/(p-1); nVG= 1/(1-m);
alpVG=(((-101.9368*Pe)*(2*p*(p-1))/(p+3))*((55.6+7.4*p+p^2)/(147.8+8.1*p+0.092*p^2)))^-1; %%[1/mm]%;
%%%
Osat=Osat*ones(1,ms);
L=L*ones(1,ms);
Pe = Pe*ones(1,ms);
O33 = O33*ones(1,ms);
alpVG= alpVG*ones(1,ms); %% [1/mm]
nVG= nVG*ones(1,ms); %% [-]
Ks_Zs= Ks*ones(1,ms); %%[mm/h]
%%%%
%%%%%%%%%%% Matric Potential
Kfc = 0.2; %% [mm/h]
Phy = 10000; %% [kPa]
%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
[Ofc,Oss_Lp,Owp_Lp,Ohy]=Soil_parametersII(ms,Osat,L,Pe,Ks_Zs,O33,nVG,alpVG,Kfc,1,1,Phy);

if sum(O(1,:)) ==0
    O(1,:) = Ohy+0.2;
end

clear Oss_Lp Owp_Lp 
Oice = 0;
%%%%%%%%%%%%%
Zdes = 10;
Zinf=  10;
Zbio = MOD_PARAM_input.Zbio;

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
%%%% ICE Parameter 
Ice_wc_sp =0.01; %% [-] Specific Maximum water content ice 
ros_Ice_thr = 500 ; %% [kg/m^3] Density Thrshold to transform snow into ice 
Aice = 0.28; %% [-] Ice albedo
WatFreez_Th = -8; %% [°C] Threshold for freezing lake water
dz_ice = 0.45; %% [mm / h] Water Freezing Layer progression without snow-layer
%%%%%%%%%%%%%%%%%%%%
ExEM = 0.0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Root Parameters 
[RfH_Zs,RfL_Zs]=Root_Fraction_General(Zs,CASE_ROOT,ZR95_H,ZR50_H,ZR95_L,ZR50_L,ZRmax_H,ZRmax_L); 
%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%
L_day=zeros(NNd,1);
for j=2:24:NN
 [h_S,delta_S,zeta_S,T_sunrise,T_sunset,L_day(j)]= SetSunVariables(Datam(j,:),DeltaGMT,Lon,Lat,t_bef,t_aft);
end 
Lmax_day = max(L_day); 
if strcmp(id_location,'Kelland')
    Lmax_day = 18;
end

clear('h_S','delta_S','zeta_S','T_sunrise','T_sunset','L_day')


%% Initialisation
%%%%%%%%%%%%%%%%%%%%%%
Ci_sunL(1,:) = [Ca(1)]; % %% [umolCO2/mol]
Ci_sunH(1,:) = [Ca(1)]; %% [umolCO2/mol]
Ci_shdL(1,:) = [Ca(1)]; % %% [umolCO2/mol]
Ci_shdH(1,:) = [Ca(1)]; %% [umolCO2/mol]
%%%%%%%%%%%%%%%%
LAI_H(1,:)=[0.01]*ones(1,cc); %
B_H(1,:,:)= [304 461 507 309 35 0 76 0 ].*ones(cc,1); 
Rrootl_H(1,:)= [5408]*ones(1,cc) ;
PHE_S_H(1,:)=[3 ]*ones(1,cc);
dflo_H(1,:)=[115 ]*ones(1,cc);
AgeL_H(1,:)=[226 ]*ones(1,cc);
e_rel_H(1,:)=[1 ]*ones(1,cc);
hc_H(1,:) =[11.0]*ones(1,cc); %% 
SAI_H(1,:) = [0.1 ]*ones(1,cc); %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
LAI_L(1,:)=[0.66]*ones(1,cc);
B_L(1,:,:)= [66 0 165 118 10 0 2 0 ].*ones(cc,1);
Rrootl_L(1,:)= [2087]*ones(1,cc) ;
PHE_S_L(1,:)=[3]*ones(1,cc);
dflo_L(1,:)=[72 ]*ones(1,cc);
AgeL_L(1,:)=[41]*ones(1,cc);
e_rel_L(1,:)=[1]*ones(1,cc);
hc_L(1,:) =[0.2]*ones(1,cc);
SAI_L(1,:) = [0.001]*ones(1,cc);

if strcmp(id_location,'Kelland')

    Rrootl_L(1,:) = [12] ;
    PHE_S_L(1,:)=[1];
    dflo_L(1,:)=0;
    AgeL_L(1,:)=[14];
    e_rel_L(1,:)=[0];
    
    dmg_L = [120];

    LAI_L(1,:)=[0.1];
    B_L(1,:,:)= [8 0 45 300 5.68 0 16.3 0];

    ff_r_L = 0.5; %%% allocation to flower and fruit to reproduction
    LtR_L = 2; %%% Leaf to Root ratio maximum
    Mf_L = 1/200; %% fruit maturation turnover [1/d]
    Trr_L = 0.1; %% Translocation rate [gC /m^2 d]

    age_cr_L = 200;
    LDay_cr_L = 8;
    Sl_L = 0.02;

    O(1,:) = 0.27;   
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BLit(1)=0.0;  %% [kg DM /m2 PFT] Litter Biomass
%%%%%%%%%%
%%%%
PARI_H(1,:,:)= [129.86 129.62 0.7377].*ones(cc,1);
NBLI_H(1,:)=  [ 1.697 ]*ones(1,cc);
PARI_L(1,:,:)=[0 0 0].*ones(cc,1);
NBLI_L(1,:)=[0 ]*ones(1,cc);

Nreserve_H(1,:)= [100 ]*ones(1,cc);
Preserve_H(1,:)= [100 ]*ones(1,cc);
Kreserve_H(1,:)= [100 ]*ones(1,cc);
Careserve_H(1,:)= [100 ]*ones(1,cc);
Mgreserve_H(1,:)= [100 ]*ones(1,cc);
Sireserve_H(1,:)= [100 ]*ones(1,cc);

FNC_H(1,:)=[1 ]*ones(1,cc);
NupI_H(1,:,:)= [0 0 0 0 0 0].*ones(cc,1);

Nreserve_L(1,:)= [100 ]*ones(1,cc);
Preserve_L(1,:)= [100 ]*ones(1,cc);
Kreserve_L(1,:)= (100 )*ones(1,cc);
Careserve_L(1,:)= [100 ]*ones(1,cc);
Mgreserve_L(1,:)= [100 ]*ones(1,cc);
Sireserve_L(1,:)= [100 ]*ones(1,cc);

FNC_L(1,:)=[1 ]*ones(1,cc);
NupI_L(1,:,:)= [0 0 0 0 0 0].*ones(cc,1);
RexmyI(1,:)= [0 0 0];

%%%
if OPT_SoilBiogeochemistry == 1
    
    % p pools and relevant values from spin-up
    if exist("B_out",'var')
        % outputs from BG spinup 
        P(1,:)=  B_out;

        RexmyI(1,:)= RexmyI_out;
        NavlI(1,:)=NavlI_out;
        NupI_L(1,:,:) = NupI_L_out;

        Nreserve_L(1,:)= 5;
        Preserve_L(1,:)= 1;
        Kreserve_L(1,:)= 1;
        Careserve_L(1,:)= 3;
        Mgreserve_L(1,:)= 1;
        Sireserve_L(1,:)= 5;

    else
        RexmyI(1,:)= [0.020519937248942   0.121476923926910   0 ];
        NavlI(1,:)=[0.286 0.0036 0.048 0.048 0.048 0.048];

        P(1,:)=  [0.0133276436864202	2.13110765064771	0.357464305487915	-0.00407166864001520	-0.00407166864001520	0.00217674015777382	0.525603427859736	0.0584003808733038,...
                613.502506368818	231.466128209168	5170.14541700210	1.79029952549020	7.08410832387938e-05	1.39768842499527e-11	1.08438335366089e-17	6.16625917420274e-12,...
        	    1.80730558943482e-17	-0.296863560015633	6.95636740000000e-316	1.07074456000000e-315	-0.144050930151499	1.00194115466345	0.721982208853361	-0.00291458664098715,...
        	    0.144162785073935	22308631228.0340	0.00613917902104300	0.000310781579205794	-9721962836.01395	-0.0708730576454231	0.557772000567600	0.0281257365905966,...
        	    0.000146315041968597	0.100194115465843	0.645228867937500	-0.00291458664098715	0.125375348189311	12813287742.4868	3.47192733320546e-05	8.08141187070576e-11,...
        	    -9721962836.01395	-0.0708730576454231	1.05431478083063	150	21.3503548114726	24923678015.0349	8.50225016331410e-06	0.738634467466490	-0.00291458664098715,...
        	    0.0302855144241957	239914228.791479	2.99951662237008	3.02124081900651	32.0315217585881	503.424657524299	0.759488137929802	-0.00291458664098715,...
        	    0.0314455147316317	359906109.321036	3.02206968620830	3.04396971412664	32.1183962970903	503.424657524299	0.710233754407748	-0.00291458664098715,...
        	    0.0287057044836347	239882634.876848	2.96924764035710	2.99073527177415	31.9210275225038	503.424657524299	0.738634467466490	-0.00291458664098715,...
        	    0.0302855144241957	239914228.791479	2.99951662237008	3.02124081900651	32.0315217585881	503.424657524299];
    end
    %%%%%%%%%%%%%%%%%
    
    if strcmp(id_location,'Kelland')

        %%%
        Upl=0.00; %% Soil production [mm/yr]
        HIST=0;

        [B_IO]=Biogeochemistry_IO(Zs(ms+1),Lat,Lon,Upl);
        B_IO.FertN = zeros(1,366);
        B_IO.FertN([10:5:115])=0.82; 
        B_IO.DepN = 1e-6;
        B_IO.DepP = 1e-6;
        B_IO.DepK = 1e-6;
        B_IO.DepCa = 1e-6;
        B_IO.DepMg = 1e-6;
        B_IO.DepSi = 1e-6;
        B_IO.N_Man = 1e-6;
        B_IO.P_Man = 1e-6;
        B_IO.K_Man = 1e-6;
        B_IO.Ca_Man = 1e-6;
        B_IO.Mg_Man = 1e-6;
        B_IO.Si_Man = 1e-6;
        B_IO.Lig_fr_Man = 1e-6;

    else
        %%%
        Upl=0.01; %% Soil production [mm/yr]
        HIST=0;
            
        FertN=zeros(1,366);
        FertP=zeros(1,366);
        FertK=zeros(1,366);
        FertCa=zeros(1,366);
        FertMg=zeros(1,366);
        FertSi=zeros(1,366);

        [B_IO]=Biogeochemistry_IO(Zs(ms+1),Lat,Lon,Upl,HIST,FertN,FertP,FertK,FertCa,FertMg,FertSi);
    end

    B_IO.SC_par=[1 1 1 1];
    PH(1) = MOD_PARAM_input.PH;
    
end
%%%%
TBio_L=[1]*ones(1,cc);  %%[ton DM / ha ]
TBio_H=[300 ]*ones(1,cc); %[ton DM / ha ]
%%%%%%%%%%%%%%%%%
Vx_H=[10 ]*ones(1,cc);  %% [mm/ m2 PFT];
Vl_H=[10 ]*ones(1,cc);  %% [mm/ m2 PFT];
Vx_L=[10 ]*ones(1,cc);   %% [mm/ m2 PFT];
Vl_L=[100 ]*ones(1,cc);   %% [mm/ m2 PFT];
%%%%%%%%%%%%%%%%%%%% Initial Conditions
SWE(1)=0; %% [mm]
SND(1)=0;
Tdamp(1) =Ta(1); 

if ~isnan(MOD_PARAM_input.soil_temp)
    Ts(1) = MOD_PARAM_input.soil_temp;
    Tdp(1,:) = MOD_PARAM_input.soil_temp*ones(1,ms);
else
    Ts(1) = Ta(1);
    Tdp(1,:) = Ta(1);
end

Sdp(1,:) = (0.001*dz).*cv_s.*(Tdp(1,:)); %% [J °C/m^2 K]
TdpI_H(1,:)=Tdp(1,1);
TdpI_L(1,:)=Tdp(1,1);
% for stablising r_aut(1:2)
Tdp_H(1,:)=Tdp(1,1);
Tdp_L(1,:)=Tdp(1,1);

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
Ws_under(1,:)=1; 
%%%%%%%%%%%%%% Volume [mm]
%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V(1,:) = (O(1,:) -Ohy).*dz;
%%%%%%%%%%%%%%%%%
