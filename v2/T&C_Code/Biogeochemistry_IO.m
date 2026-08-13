%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction Biogeochemistry_Input_output  %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[B_IO]=Biogeochemistry_IO(Zs,Lat,Lon,Upl,HIST,FertN,FertP,FertK,FertCa,FertMg,FertSi,ManF,N_Man,P_Man,K_Man,Lig_fr_Man,DepN,DepP,DepK)
%%%%%%%%%%%%%%
if nargin <= 4
    %%% External_Input_Output
    [DepN,DepP,DepK,...
        DepCa,DepMg,DepSi,...
        FertN,FertP,FertK,...
        FertCa,FertMg,FertSi,...
        ManF,N_Man,P_Man,K_Man,...
        Ca_Man,Mg_Man,Si_Man,...
        Lig_fr_Man]=Biogeochemistry_IP(Lat,Lon,0);
end
if nargin == 5
    [DepN,DepP,DepK,...
        DepCa,DepMg,DepSi,...
        FertN,FertP,FertK,...
        FertCa,FertMg,FertSi,...
        ManF,N_Man,P_Man,K_Man,...
        Ca_Man,Mg_Man,Si_Man,...
        Lig_fr_Man]=Biogeochemistry_IP(Lat,Lon,HIST);
else
    if  nargin <= 8+6 && nargin > 4
        [DepN,DepP,DepK,...
            DepCa,DepMg,DepSi,...
            ~,~,~,...
            ~,~,~,...
            ManF,N_Man,P_Man,K_Man,...
            Ca_Man,Mg_Man,Si_Man,...
            Lig_fr_Man]=Biogeochemistry_IP(Lat,Lon,HIST);
    end
    if nargin <= 13+6 && nargin > 4+6
        %%% External_Input_Output
        [DepN,DepP,DepK,...
            DepCa,DepMg,DepSi]=Biogeochemistry_IP(Lat,Lon,HIST);
    end
end
%%% External_ Soil Production
[Tup_P,Tup_K,Tup_Ca,Tup_Mg,Tup_Si]=Biogeochemistry_Soil_Prod(Zs,Upl);
%%%%%%%%%%%%%%%
B_IO.FertN = FertN;
B_IO.FertP = FertP;
B_IO.FertK = FertK;

B_IO.FertCa = FertCa;
B_IO.FertMg = FertMg;
B_IO.FertSi = FertSi;

B_IO.DepN = DepN;
B_IO.DepP =DepP;
B_IO.DepK= DepK;

B_IO.DepCa= DepCa;
B_IO.DepMg= DepMg;
B_IO.DepSi= DepSi;
%%%%%%%%%%%
B_IO.Tup_P=Tup_P;
B_IO.Tup_K=Tup_K;

B_IO.Tup_Ca=Tup_Ca;
B_IO.Tup_Mg=Tup_Mg;
B_IO.Tup_Si=Tup_Si;
%%%%% Scaling parameter
B_IO.SC_par=1;
%%%%% Manure parameters 
B_IO.ManF = ManF;  %%% [gC /m2 day]
B_IO.N_Man = N_Man; %% Manure  [gC/gN]
B_IO.P_Man = P_Man; % Manure  [gC/gP]
B_IO.K_Man = K_Man;% Manure  [gC/gK]

B_IO.Ca_Man = Ca_Man;% Manure  [gC/gCa]
B_IO.Mg_Man = Mg_Man;% Manure  [gC/gMg]
B_IO.Si_Man = Si_Man;% Manure  [gC/gSi]

B_IO.Lig_fr_Man = Lig_fr_Man ; %% Lignin fraction in Manure  [g Lignin / g DM]
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction Biogeochemistry_Input_output  %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[DepN,DepP,DepK,...
    DepCa,DepMg,DepSi,...
    FertN,FertP,FertK,...
    FertCa,FertMg,FertSi,...
    ManF,N_Man,P_Man,K_Man,...
    Ca_Man,Mg_Man,Si_Man,...
    Lig_fr_Man]=Biogeochemistry_IP(Lat,Lon,HIST)
%%%%%%%%%%%%%%
%%%% External Inputs [gX/m2 day]
load("Inputs\All_deposition_data.mat")
if HIST == 1 %%%% historical pre-industrial depositions
    DepN = griddata(LonN_Gl,LatN_Gl,Ndep1860',Lon,Lat,'nearest'); DepN=DepN/1000/365;
    DepP = interp2((LonP),(LatP),(DepP_hist'),Lon,Lat,'nearest'); DepP=DepP/1000/365;
else
    DepN = griddata(LonN,LatN,Ndep_tot,Lon,Lat,'nearest'); DepN=DepN/10/365;
    DepP = interp2((LonP),(LatP),(DepP'),Lon,Lat,'nearest'); DepP=DepP/1000/365;
end
DepK =  griddata(LonCC{3},LatCC{3},CC{3},Lon,Lat,'nearest'); DepK=DepK/10/365;

DepCa = DepK;
DepMg = DepK;
DepSi = 0;

%%%%%
FertN=0*ones(1,366);
FertP=0*ones(1,366);
FertK=0*ones(1,366);

FertCa=0*ones(1,366);
FertMg=0*ones(1,366);
FertSi=0*ones(1,366);
%%%% Erosion and export of Elements [gX/m2 d]
Export_K=0;
Export_P=0;
Export_N=0;
Export_C=0;
%%%%%%
ManF=0*ones(1,366);  %%% [gC /m2 day]
N_Man=1; %% Manure  [gC/gN]
P_Man=1; % Manure  [gC/gP]
K_Man=1;% Manure  [gC/gK]

Ca_Man=1;% Manure  [gC/gCa]
Mg_Man=1;% Manure  [gC/gMg]
Si_Man=1;% Manure  [gC/gSi]

Lig_fr_Man=0.01; %% Lignin fraction in Manure  [g Lignin / g DM]
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction Soil_production              %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[Tup_P,Tup_K,Tup_Ca,Tup_Mg,Tup_Si]=Biogeochemistry_Soil_Prod(Zs,Upl)
%%%% Tectoni Uplift
%Upl=0.05; %% Soil production [mm/yr]
rckd = 2500; %%[kg/m3]
Cel_K = 0.021 ; % 0.021;[%]
Cel_P=  0.0005; %% [%] Yang et al 2013 Biogeosciences 

Cel_Ca = 0.021 ; % 0.021;[%]
Cel_Mg = 0.021 ; % 0.021;[%]
Cel_Si = 0.021 ; % 0.021;[%]

%Tup = Upl*Cel*rckd; [gX/m2 yr]
Tup_P= Upl*Cel_P*rckd/365; %% [gP/m2 d]
Tup_K= Upl*Cel_K*rckd/365; %% [gK/m2 d]

Tup_Ca= Upl*Cel_Ca*rckd/365; %% [gCa/m2 d]
Tup_Mg= Upl*Cel_Mg*rckd/365; %% [gMg/m2 d]
Tup_Si= Upl*Cel_Si*rckd/365; %% [gSi/m2 d]

return