%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction BIOGEOCHEMISTRY_DYNAMIC   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[NH4_Uptake,NO3_Uptake,P_Uptake,K_Uptake,...
    Ca_Uptake,Mg_Uptake,Si_Uptake]= Biogeo_uptake(B,Broot,Ts,T,VT,Ccrown,EM,BiogeoPar)
%%%%%%%%%% Plant Nutrient Uptake
%fT= 1.5^((Ts-Tref)/10); Tref=25;
%fT=power((70-Ts)/(70-38),12)*exp(12*(Ts-38)/(70-38)); %%% Xu-Ri and Prentice 2008
fT=1; %% no Temperature depedence on nutrient uptake
%%%%
T(T>VT)=VT;
%%%%
%%%%%%
BmycAM=B(20); BmycAM(BmycAM<0)=0; %% [gC/m2] AM-Mycorrhizal - C
BmycEM=B(21); BmycEM(BmycEM<0)=0; %% [gC/m2] EM-Mycorrhizal - C
%k=100; %% Hinsiger et al 2011 PLSO (Based also on Hypahe lenght vs Root length)
%Broot = Ccrown*(Broot + k*Bmyc); %% [gC/m2]
Broot = Ccrown*(Broot);
%%%%
%%%%%%%%%% Solubility coeffients
aNH4 =BiogeoPar.aNH4;  %%0.1 CLM4
aNO3 = BiogeoPar.aNO3;
aP = BiogeoPar.aP;
aK = BiogeoPar.aK;

aCa = BiogeoPar.aCa;
aMg = BiogeoPar.aMg;
aSi = BiogeoPar.aSi;

%%%%%%%%%%%%
OPT_UPT=2; 
%%%%%%%%%%%%
if OPT_UPT==1
    %%%% Active Uptake Root
    vNH4 = BiogeoPar.vr; %% 0.09148; % [gN /gC day]
    vNO3 = BiogeoPar.vr;
    vP = BiogeoPar.vr/8;
    vK = BiogeoPar.vr/2.5;
    
    vCa = BiogeoPar.vr/2.5;
    vMg = BiogeoPar.vr/2.5;
    vSi = BiogeoPar.vr;

    KnNH4 = BiogeoPar.Kvr ; %%[gN /gC ]
    KnNO3 = BiogeoPar.Kvr ; %%[gN /gC ]
    KnP=   BiogeoPar.Kvr/8 ; %%[gP /gC ]
    KnK = BiogeoPar.Kvr/2.5 ; %%[gK /gC ]

    KnCa = BiogeoPar.Kvr/2.5 ; %%[gCa /gC ]
    KnMg = BiogeoPar.Kvr/2.5 ; %%[gMg /gC ]
    KnSi = BiogeoPar.Kvr/2.5 ; %%[gSi /gC ]

    NH4_Uptake_r = vNH4*B(31)/(B(31)+KnNH4*Broot)*fT*Broot; %% [gN/m2 d]
    NO3_Uptake_r =vNO3*B(32)/(B(32)+KnNO3*Broot)*fT*Broot; %% [gN/m2 d]
    P_Uptake_r = vP*B(43)/(B(43)+KnP*Broot)*fT*Broot; %% [gP/m2 d]
    K_Uptake_r = vK*B(52)/(B(52)+KnK*Broot)*fT*Broot; %% [gK/m2 d]

    Ca_Uptake_r = vCa*B(60)/(B(60)+KnCa*Broot)*fT*Broot; %% [gCa/m2 d]
    Mg_Uptake_r = vMg*B(68)/(B(68)+KnMg*Broot)*fT*Broot; %% [gMg/m2 d]
    Si_Uptake_r = vSi*B(76)/(B(76)+KnSi*Broot)*fT*Broot; %% [gSi/m2 d]


    %%%% Active Uptake EM Mycorrhizal
    vem = BiogeoPar.vem ; % 7.3340e-006 ;  %%[gN /gC^2 day]
    KnNH4 = BiogeoPar.Kvem; % 16.67e-005 ; %%[gN /gC ]
    KnNO3 = BiogeoPar.Kvem ; %%[gN /gC ]
    KnP=   BiogeoPar.Kvem/8 ; %%[gP /gC ]
    KnK =  BiogeoPar.Kvem/2.5 ; %%[gK /gC ]

    KnCa =  BiogeoPar.Kvem/2.5 ; %%[gCa /gC ]
    KnMg =  BiogeoPar.Kvem/2.5 ; %%[gMg /gC ]
    KnSi =  BiogeoPar.Kvem/2.5 ; %%[gSi /gC ]

    NH4_Uptake_em = vem*B(31)/(B(31)+KnNH4*Broot)*fT*Broot*BmycEM; %% [gN/m2 d]
    NO3_Uptake_em = vem*B(32)/(B(32)+KnNO3*Broot)*fT*Broot*BmycEM; %% [gN/m2 d]
    P_Uptake_em = vem*B(43)/(B(43)+KnP*Broot)*fT*Broot*BmycEM; %% [gP/m2 d]
    K_Uptake_em = vem*B(52)/(B(52)+KnK*Broot)*fT*Broot*BmycEM; %% [gK/m2 d]

    Ca_Uptake_em = vem*B(60)/(B(60)+KnCa*Broot)*fT*Broot*BmycEM; %% [gCa/m2 d]
    Mg_Uptake_em = vem*B(68)/(B(68)+KnMg*Broot)*fT*Broot*BmycEM; %% [gMg/m2 d]
    Si_Uptake_em = vem*B(76)/(B(76)+KnSi*Broot)*fT*Broot*BmycEM; %% [gSi/m2 d]

    %%%% Active Uptake AM Mycorrhizal
    vam =  BiogeoPar.vam ; %% 4.4000e-005 ;  %%[gN /gC^2 day]
    KnNH4 = BiogeoPar.Kvam; %%  0.002 ; %%[gN /gC ]
    KnNO3 =   BiogeoPar.Kvam  ; %%[gN /gC ]
    KnP=   BiogeoPar.Kvam/8 ; %%[gP /gC ]
    KnK =  BiogeoPar.Kvam/2.5 ; %%[gK /gC ]

    KnCa =  BiogeoPar.Kvam/2.5 ; %%[gCa /gC ]
    KnMg =  BiogeoPar.Kvam/2.5 ; %%[gMg /gC ]
    KnSi =  BiogeoPar.Kvam/2.5 ; %%[gSi /gC ]

    NH4_Uptake_am = vam*B(31)/(B(31)+KnNH4*Broot)*fT*Broot*BmycAM; %% [gN/m2 d]
    NO3_Uptake_am = vam*B(32)/(B(32)+KnNO3*Broot)*fT*Broot*BmycAM; %% [gN/m2 d]
    P_Uptake_am = vam*B(43)/(B(43)+KnP*Broot)*fT*Broot*BmycAM; %% [gP/m2 d]
    K_Uptake_am = vam*B(52)/(B(52)+KnK*Broot)*fT*Broot*BmycAM; %% [gK/m2 d]

    Ca_Uptake_am = vam*B(60)/(B(60)+KnCa*Broot)*fT*Broot*BmycAM; %% [gCa/m2 d]
    Mg_Uptake_am = vam*B(68)/(B(68)+KnMg*Broot)*fT*Broot*BmycAM; %% [gMg/m2 d]
    Si_Uptake_am = vam*B(76)/(B(76)+KnSi*Broot)*fT*Broot*BmycAM; %% [gSi/m2 d]

end
if OPT_UPT==2
    vr= BiogeoPar.vr;%% [m^2/day]
    vem =   BiogeoPar.vem;%% [m^2/day]
    vam =   BiogeoPar.vam;%% [m^2/day]
    g_r = BiogeoPar.g_r;
    g_em =  BiogeoPar.g_em;
    g_am =  BiogeoPar.g_am;
    %%%%%%%%%%%%%%%%%
    NH4_Uptake_r =(vr*B(31)*(Broot.^1.5))./g_r; %% [gN/m2 d]
    NO3_Uptake_r =(vr*B(32)*(Broot.^1.5))./g_r; %% [gN/m2 d]
    P_Uptake_r = (vr*B(43)*(Broot.^1.5))./g_r; %% [gP/m2 d]
    K_Uptake_r = (vr*B(52)*(Broot.^1.5))./g_r; %% [gK/m2 d]

    Ca_Uptake_r = (vr*B(60)*(Broot.^1.5))./g_r; %% [gCa/m2 d]
    Mg_Uptake_r = (vr*B(68)*(Broot.^1.5))./g_r; %% [gMg/m2 d]
    Si_Uptake_r = (vr*B(76)*(Broot.^1.5))./g_r; %% [gSi/m2 d]

    %%%%%
    NH4_Uptake_em =(vem*B(31)*(BmycEM.^1.5))./g_em; %% [gN/m2 d]
    NO3_Uptake_em =(vem*B(32)*(BmycEM.^1.5))./g_em; %% [gN/m2 d]
    P_Uptake_em = (vem*B(43)*(BmycEM.^1.5))./g_em; %% [gP/m2 d]
    K_Uptake_em = (vem*B(52)*(BmycEM.^1.5))./g_em; %% [gK/m2 d]

    Ca_Uptake_em = (vem*B(60)*(BmycEM.^1.5))./g_em; %% [gCa/m2 d]
    Mg_Uptake_em = (vem*B(68)*(BmycEM.^1.5))./g_em; %% [gMg/m2 d]
    Si_Uptake_em = (vem*B(76)*(BmycEM.^1.5))./g_em; %% [gSi/m2 d]

    %%%%%
    NH4_Uptake_am = (vam*B(31)*(BmycAM.^1.5))./g_am; %% [gN/m2 d]
    NO3_Uptake_am = (vam*B(32)*(BmycAM.^1.5))./g_am; %% [gN/m2 d]
    P_Uptake_am = (vam*B(43)*(BmycAM.^1.5))./g_am; %% [gP/m2 d]
    K_Uptake_am = (vam*B(52)*(BmycAM.^1.5))./g_am; %% [gK/m2 d]

    Ca_Uptake_am = (vam*B(60)*(BmycAM.^1.5))./g_am; %% [gCa/m2 d]
    Mg_Uptake_am = (vam*B(68)*(BmycAM.^1.5))./g_am; %% [gMg/m2 d]
    Si_Uptake_am = (vam*B(76)*(BmycAM.^1.5))./g_am; %% [gSi/m2 d]
end
%%%%%%%%%%%%%
%%%%%%%%%% Active uptake
NH4_Uptake_a = NH4_Uptake_em*EM + NH4_Uptake_am*(1-EM) + NH4_Uptake_r ;
NO3_Uptake_a = NO3_Uptake_em*EM + NO3_Uptake_am*(1-EM) + NO3_Uptake_r ;
P_Uptake_a = P_Uptake_em*EM + P_Uptake_am*(1-EM) + P_Uptake_r ;
K_Uptake_a = K_Uptake_em*EM + K_Uptake_am*(1-EM) + K_Uptake_r ;

Ca_Uptake_a = Ca_Uptake_em*EM + Ca_Uptake_am*(1-EM) + Ca_Uptake_r ;
Mg_Uptake_a = Mg_Uptake_em*EM + Mg_Uptake_am*(1-EM) + Mg_Uptake_r ;
Si_Uptake_a = Si_Uptake_em*EM + Si_Uptake_am*(1-EM) + Si_Uptake_r ;

NH4_Uptake_a =  max(0,min(Ccrown*B(31),NH4_Uptake_a)); %%
NO3_Uptake_a =  max(0,min(Ccrown*B(32),NO3_Uptake_a));
P_Uptake_a =   max(0,min(Ccrown*B(43),P_Uptake_a));
K_Uptake_a =  max(0,min(Ccrown*B(52),K_Uptake_a));

Ca_Uptake_a =  max(0,min(Ccrown*B(60),Ca_Uptake_a));
Mg_Uptake_a =  max(0,min(Ccrown*B(68),Mg_Uptake_a));
Si_Uptake_a =  max(0,min(Ccrown*B(76),Si_Uptake_a));

%%%%%%%%%%%%%%%%% OPT. Passive
if T>0
    NH4_Uptake_p = max(0,min(Ccrown*B(31),aNH4*B(31)*(T)/VT)); %% [gN/m2 d]
    NO3_Uptake_p = max(0,min(Ccrown*B(32),aNO3*B(32)*(T)/VT)); %% [gN/m2 d]
    P_Uptake_p = max(0,min(Ccrown*B(43),aP*B(43)*(T)/VT)); %% [gP/m2 d]
    K_Uptake_p = max(0,min(Ccrown*B(52),aK*B(52)*(T)/VT)); %% [gK/m2 d]

    Ca_Uptake_p = max(0,min(Ccrown*B(60),aCa*B(60)*(T)/VT)); %% [gCa/m2 d]
    Mg_Uptake_p = max(0,min(Ccrown*B(68),aMg*B(68)*(T)/VT)); %% [gMg/m2 d]
    Si_Uptake_p = max(0,min(Ccrown*B(76),aSi*B(76)*(T)/VT)); %% [gSi/m2 d]

else
    NH4_Uptake_p = 0; NO3_Uptake_p =0 ; P_Uptake_p=0;  K_Uptake_p =0;

    Ca_Uptake_p =0;Mg_Uptake_p=0;Si_Uptake_p=0;
end
%%%%%%%%%%%%
NH4_Uptake= max(NH4_Uptake_a,NH4_Uptake_p);
NO3_Uptake=max(NO3_Uptake_a,NO3_Uptake_p);
P_Uptake=max(P_Uptake_a,P_Uptake_p);
K_Uptake=max(K_Uptake_a,K_Uptake_p);

Ca_Uptake=max(Ca_Uptake_a,Ca_Uptake_p);
Mg_Uptake=max(Mg_Uptake_a,Mg_Uptake_p);
Si_Uptake=max(Si_Uptake_a,Si_Uptake_p);
%%%%%
if Broot == 0
    NH4_Uptake = 0;
    NO3_Uptake =0;
    P_Uptake = 0;
    K_Uptake = 0;
    Ca_Uptake = 0;
    Mg_Uptake = 0;
    Si_Uptake = 0;
end
%%%%
return




