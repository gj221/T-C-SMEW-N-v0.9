%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction BIOGEOCHEMISTRY_DYNAMIC   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[LEAK_NH4,LEAK_NO3,LEAK_P,LEAK_K,LEAK_DOC,LEAK_DON,LEAK_DOP,...
    LEAK_Ca,LEAK_Mg,LEAK_Si...
    ]= Biogeo_Leakage(B,Lk,V,BiogeoPar)  
%%%%%%%%%%%%%% Leakage of Nitrogen Compound 
%%%% Proportional to leakage from the nitrogen zone with fraction of
%%%% dissolved ammonium and nitrate 
%%% For simplicity leakage only at the base of soil column 
%%%
Lk(Lk>V)=V; 
%%%% Porporato et al., 2003 
aNH4 =BiogeoPar.aNH4;  
aNO3 = BiogeoPar.aNO3;  
aP = BiogeoPar.aP;
aK = BiogeoPar.aK; 
aDOC = BiogeoPar.aDOC; 
aDON= BiogeoPar.aDON;
aDOP= BiogeoPar.aDOP;

aCa = BiogeoPar.aCa;
aMg = BiogeoPar.aMg;
aSi = BiogeoPar.aSi;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
LEAK_NH4 = max(0,min(B(31),aNH4*B(31)*(Lk)/V)); %% [gN/m2 d]
LEAK_NO3 = max(0,min(B(32),aNO3*B(32)*(Lk)/V)); %% [gN/m2 d]
LEAK_P = max(0,min(B(43),aP*B(43)*Lk/V)); %% [gP/m2 d]
LEAK_K = max(0,min(B(52),aK*B(52)*Lk/V)); %% [gK/m2 d]
LEAK_DOC = max(0,min(B(12)+B(13),aDOC*(B(12)+B(13))*Lk/V)); %% [gC/m2 d]
LEAK_DON = max(0,min(B(33),aDON*B(33)*Lk/V)); %% [gN/m2 d]
LEAK_DOP = max(0,min(B(47),aDOP*B(47)*Lk/V)); %% [gP/m2 d]

LEAK_Ca = max(0,min(B(60),aCa*B(60)*Lk/V)); %% [gCa/m2 d]
LEAK_Mg = max(0,min(B(68),aMg*B(68)*Lk/V)); %% [gMg/m2 d]
LEAK_Si = max(0,min(B(76),aSi*B(76)*Lk/V)); %% [gSi/m2 d]


%%%%%%%%%%%%%%%%%%
return 
