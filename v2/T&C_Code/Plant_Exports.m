%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  PLANT EXPORTS   %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[TexC,TexN,TexP,TexK,...
    TexCa,TexMg,TexSi,...
    TNIT,TPHO,TPOT,...
    TCAL,TMAG,TSIL,...
    NuLIT,NreserveM,PreserveM,KreserveM,...
    CareserveM,MgreserveM,SireserveM,...
    SupN,SupP,SupK,...
    SupCa,SupMg,SupSi,...
    IS]= Plant_Exports(B,Btm1,...
    NuLITtm1,Slf,Sfr,Swm,Sll,Sr,...
    Rexmy,St,Mpar,fab,fbe,RB,Nreserve,Preserve,Kreserve,...
    Careserve,Mgreserve,Sireserve,...
    rNc,rPc,rKc,...
    rCac,rMgc,rSic,...
    ManI,OPT_SB,id_location)  
%%%%%%%%%
%%%%-> All computations are for /m2 PFT not ground 
%%%% INPUT
dtd = 1 ; %% [day] 
%%%%%%%%%%%%%% Stochiometry  - Concentration -
Nl = St.Nl;  %%%  [gC/gN]
Ns = St.Ns;  %%%  [gC/gN]
Nr = St.Nr ; %%%  [gC/gN]
Nc= Ns; %%%  [gC/gN] Carbohydrate Reserve Carbon Nitrogen
Nf = St.Nf; %%%  [gC/gN]
Nh = St.Nh;  %%%  [gC/gN]
%%%
Phol = St.Phol;  %%%  [gC/gP]
Phos = St.Phos;  %%%  [gC/gP]
Phor = St.Phor ; %%%  [gC/gP]
Phoc= Phos; %%%  [gC/gP] Carbohydrate Reserve Carbon Phosophorus
Phof = St.Phof; %%%  [gC/gP]
Phoh = St.Phoh;  %%%  [gC/gP]
%%%
Kpotl = St.Kpotl;  %%%  [gC/gK]
Kpots = St.Kpots;  %%%  [gC/gK]
Kpotr = St.Kpotr ; %%%  [gC/gK]
Kpotc= Kpots; %%%  [gC/gK] Carbohydrate Reserve Carbon Potassium
Kpotf = St.Kpotf; %%%  [gC/gK]
Kpoth = St.Kpoth;  %%%  [gC/gK]
ftransL = St.ftransL;
ftransR = St.ftransR;

%%% Cal
Call = St.Call;  %%%  [gC/gCa]
Cals = St.Cals;  %%%  [gC/gCa]
Calr = St.Calr ; %%%  [gC/gCa]
Calc= Cals; %%%  [gC/gCa] Carbohydrate Reserve Carbon Calcium
Calf = St.Calf; %%%  [gC/gCa]
Calh = St.Calh;  %%%  [gC/gCa]

%%% Mag
Magl = St.Magl;  %%%  [gC/gMg]
Mags = St.Mags;  %%%  [gC/gMg]
Magr = St.Magr ; %%%  [gC/gMg]
Magc= Mags; %%%  [gC/gMg] Carbohydrate Reserve Carbon Magnesium
Magf = St.Magf; %%%  [gC/gMg]
Magh = St.Magh;  %%%  [gC/gMg]


%%% Sil
Sill = St.Sill;  %%%  [gC/gSi]
Sils = St.Sils;  %%%  [gC/gSi]
Silr = St.Silr ; %%%  [gC/gSi]
Silc= Sils; %%%  [gC/gSi] Carbohydrate Reserve Carbon Silicon 
Silf = St.Silf; %%%  [gC/gSi]
Silh = St.Silh;  %%%  [gC/gSi]

%%%%% Lignin Fraction 
Lig_fr_l=St.Lig_fr_l; %% Lignin fraction in leaves  [g Lignin / g DM] 
Lig_fr_fr=St.Lig_fr_fr; %% Lignin fraction in fruit [g Lignin / g DM] 
Lig_fr_h=St.Lig_fr_h; %% Lignin fraction in wood [g Lignin / g DM] 
Lig_fr_r=St.Lig_fr_r; %% Lignin fraction in fine roots [g Lignin / g DM] 
%%%%%% 
%%%%%%%% Fraction of material left in the field 
fract_left =Mpar.fract_left; %% Fraction_left of leaves  and dead leaves 
fract_left_fr =Mpar.fract_left_fr; %% Fraction left of fruits 
fract_left_harAB =Mpar.fract_left_AB; %% Fraction left of harvested wood aboveground 
fract_left_harBG =Mpar.fract_left_BG; %% Fraction left of harvested wood belowground 
%%%
fract_log=Mpar.fract_log;
fire_eff = Mpar.fire_eff;  %% Effect of Fire 
funb_nit = Mpar.funb_nit; %% Fraction of unburned Nitrogen 
%%%%%%% 
FiS = St.FiS; %%[-] ; %%% Factor to increment nutrient reserve buffer 

%%ManI=; %%% Management Indicator (0) Nothing (1) Fire (-1) Logging 
%%%%%%%%%%%%%%%% Export Carbon %%%%%%%%%%%%%%%%5
%Sfr = %% Sfr Fruit maturation % [gC/m^2 d]
%Swm =  %% Wood mortality % [gC/m^2 d]
%Slf =  % Litterfall factor
%Sr= % Sr Mortality factor for fine roots
%Rexmy  %%% [Exudation , Mycorrhiza, Biol. Fix.] [gC / m^2 d]
%RB(1) Removed Live Leaves
%RB(2) Removed Sapwood
%RB(3) Removed Fine Roots 
%RB(4) Removed Carbohydrate Reserve
%RB(5) Removed Fruit and Flower
%RB(6) Removed Heartwood - Dead Sapwood
%RB(7) Removed Standing Dead Leaves 
%%%%% Exports  Nitrogen Plant [gN/m^2]
%(1-ftransL)*Slf./Nl + Sfr/Nf; % [gN/m^2 day]
%(1-ftransR)*Sr./Nr  ;  % [gN/m^2 day]
%fab*Swm/Nh; 
%fbe*Swm/Nh;
%%%%% Exports  Phosporus Plant [gP/m^2]
%(1-ftransL)*Slf./Phol + Sfr/Phof; % [gP/m^2 day]
%(1-ftransR)*Sr./Phor  ;  % [gP/m^2 day]
%fab*Swm/Phoh; 
%fbe*Swm/Phoh;
%%%%% Exports Potassium Plant [gK/m^2]
%(1-ftransL)*Slf./Kpotl + Sfr/Kpotf; % [gK/m^2 day]
%(1-ftransR)*Sr./Kpotr  ;  % [gK/m^2 day]
%fab*Swm/Kpoth; 
%fbe*Swm/Kpoth;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Total Nutrient Content -- basic stoichiometric 
TNIT =  B(1)/Nl + B(2)/Ns + B(3)/Nr + B(4)/Nc + B(5)/Nf  + B(6)./Nh  ;  %% Total Nitrogen Plant [gN/m^2]
TPHO =  B(1)/Phol + B(2)/Phos + B(3)/Phor + B(4)/Phoc + B(5)/Phof  + B(6)./Phoh  ;  %% Total Phosporus Plant [gP/m^2]
TPOT =  B(1)/Kpotl + B(2)/Kpots + B(3)/Kpotr + B(4)/Kpotc + B(5)/Kpotf  + B(6)./Kpoth  ;  %% Total Potassium Plant [gK/m^2]

TCAL =  B(1)/Call + B(2)/Cals + B(3)/Calr + B(4)/Calc + B(5)/Calf  + B(6)./Calh  ;  %% Total Calcium Plant [gCa/m^2]
TMAG =  B(1)/Magl + B(2)/Mags + B(3)/Magr + B(4)/Magc + B(5)/Magf  + B(6)./Magh  ;  %% Total Magnesium Plant [gMg/m^2]
TSIL =  B(1)/Sill + B(2)/Sils + B(3)/Silr + B(4)/Silc + B(5)/Silf  + B(6)./Silh  ;  %% Total Silicon Plant [gSi/m^2]


%%%%%%%%% Non-Structural Nutrient Cotent (Leaves, Fine Root, Reserve  Fruits)
TNITns =  B(1)/Nl +  B(3)/Nr + B(4)/Nc +  B(5)/Nf   ;  %% Total Nitrogen Plant [gN/m^2]
TPHOns =  B(1)/Phol +  B(3)/Phor + B(4)/Phoc + B(5)/Phof  ;  %% Total Phosporus Plant [gP/m^2]
TPOTns =  B(1)/Kpotl +  B(3)/Kpotr + B(4)/Kpotc  + B(5)/Kpotf ;  %% Total Potassium Plant [gK/m^2]

TCALns =  B(1)/Call +  B(3)/Calr + B(4)/Calc  + B(5)/Calf ;  %% Total Calcium Plant [gCa/m^2]
TMAGns =  B(1)/Magl +  B(3)/Magr + B(4)/Magc  + B(5)/Magf ;  %% Total Magnesium Plant [gMg/m^2]
TSILns =  B(1)/Sill +  B(3)/Silr + B(4)/Silc  + B(5)/Silf ;  %% Total Silicon Plant [gSi/m^2]

%%%% 
Nsto = FiS*B(4)*(1/Nl-1/Nc); %% Flexible storage Nitrogen Plant [gN/m^2]
Psto = FiS*B(4)*(1/Phol-1/Phoc); %%
Ksto = FiS*B(4)*(1/Kpotl-1/Kpotc); %%

Casto = FiS*B(4)*(1/Call-1/Calc); %%
Mgsto = FiS*B(4)*(1/Magl-1/Magc); %%
Sisto = FiS*B(4)*(1/Sill-1/Silc); %%

%%%%%
Nsto_tm1 = FiS*Btm1(4)*(1/Nl-1/Nc); %% Flexible storage Nitrogen Plant [gN/m^2]
Psto_tm1 = FiS*Btm1(4)*(1/Phol-1/Phoc); %%
Ksto_tm1 = FiS*Btm1(4)*(1/Kpotl-1/Kpotc); %%

Casto_tm1 = FiS*Btm1(4)*(1/Call-1/Calc); %%
Mgsto_tm1 = FiS*Btm1(4)*(1/Magl-1/Magc); %%
Sisto_tm1 = FiS*Btm1(4)*(1/Sill-1/Silc); %%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RB(RB<0) = 0; %% In case of sowing 
if sum(Btm1)==0
    Nsto_tm1 =1; Psto_tm1 =1; Ksto_tm1 =1; Casto_tm1 =1; Mgsto_tm1 =1; Sisto_tm1 =1; 
end 
if sum(B)==0
    Nsto =1; Psto =1; Ksto =1; Casto =1; Mgsto =1; Sisto =1;
end 
%%%% Real plant nutrient content 
%TNREAL = TNIT + Nreserve; 
%TPREAL = TPHO + Preserve; 
%TKREAL = TPOT + Kreserve; 
%%%%
AddN=0; AddP=0; AddK=0; AddCa=0; AddMg=0; AddSi=0; %%% Additional Nutrient (~= 0 only in case of some management)
%%%%% Nutrients on the standing lead biomass
if (B(7)+RB(7))==0 || (NuLITtm1(1)==0)
    Nlf = Nl;  Pholf = Phol; Kpotlf = Kpotl; Callf = Call; Maglf = Magl; Sillf = Sill;
else
    Nlf   = (B(7)+RB(7)*dtd)/NuLITtm1(1);  %% [gN/gC]
    Pholf = (B(7)+RB(7)*dtd)/NuLITtm1(2);   % [gP/gC]
    Kpotlf= (B(7)+RB(7)*dtd)/NuLITtm1(3);  %% [gK/gC]

    Callf= (B(7)+RB(7)*dtd)/NuLITtm1(4);  %% [gCa/gC]
    Maglf= (B(7)+RB(7)*dtd)/NuLITtm1(5);  %% [gMg/gC]
    Sillf= (B(7)+RB(7)*dtd)/NuLITtm1(6);  %% [gSi/gC]
end
%%%%%
LNIT = NuLITtm1(1) + (rNc*(1-ftransL*(rNc<=1))*Sll./Nl - Slf./Nlf -  RB(7)./Nlf)*dtd ; %% [gN/m^2]
LPHO = NuLITtm1(2) + (rPc*(1-ftransL*(rPc<=1))*Sll./Phol -  Slf./Pholf - RB(7)./Pholf)*dtd ;
LPOT = NuLITtm1(3) + (rKc*(1-ftransL*(rKc<=1))*Sll./Kpotl - Slf./Kpotlf -  RB(7)./Kpotlf)*dtd ;

LCAL = NuLITtm1(4) + (rCac*(1-ftransL*(rCac<=1))*Sll./Call - Slf./Callf -  RB(7)./Callf)*dtd ;
LMAG = NuLITtm1(5) + (rMgc*(1-ftransL*(rMgc<=1))*Sll./Magl - Slf./Maglf -  RB(7)./Maglf)*dtd ;
LSIL = NuLITtm1(6) + (rSic*(1-ftransL*(rSic<=1))*Sll./Sill - Slf./Sillf -  RB(7)./Sillf)*dtd ;

if (B(7)+RB(7))==0 %%  No Carbon 
    AddN=LNIT; AddP=LPHO; AddK=LPOT; AddCa=LCAL; AddMg=LMAG; AddSi=LSIL;
    LNIT=0; LPHO =0; LPOT=0; LCAL=0;  LMAG=0; LSIL=0;
end
NuLIT = [LNIT LPHO LPOT LCAL LMAG LSIL];  
%%%%%
%%%%%% Stop Uptake because too many reserves are available [Nsto to  Nsto+0.8*TNITns] Linear
%%%%%% 
%SupN = ((Nreserve - Nsto)./TNITns)/0.8 ; SupN(SupN>1)=1; SupN(SupN<0)=0;
%SupP=  ((Preserve - Psto)./TPHOns)/0.8 ; SupP(SupP>1)=1; SupP(SupP<0)=0;
%SupK = ((Kreserve - Ksto)./TPOTns)/0.8 ; SupK(SupK>1)=1; SupK(SupK<0)=0;
%%%%  Stop Uptake because too many reserves are available [0.8*Nsto to Nsto+0.6*TNITns] Linear
SupN = ((Nreserve - 0.8*Nsto)./(0.6*TNITns+0.2*Nsto)); SupN(SupN>1)=1; SupN(SupN<0)=0;
SupP=  ((Preserve - 0.8*Psto)./(0.6*TPHOns+0.2*Psto)); SupP(SupP>1)=1; SupP(SupP<0)=0;
SupK = ((Kreserve - 0.8*Ksto)./(0.6*TPOTns+0.2*Ksto)); SupK(SupK>1)=1; SupK(SupK<0)=0;
SupCa = ((Careserve - 0.8*Casto)./(0.6*TCALns+0.2*Casto)); SupCa(SupCa>1)=1; SupCa(SupCa<0)=0;

% if strcmp(id_location,'Kantola')
%     % kantola case
%     SupCa = ((Careserve - 0.8*Casto)./(0.1*TCALns+0.2*Casto)); SupCa(SupCa>1)=1; SupCa(SupCa<0)=0;
% else
%     %general
%     SupCa = ((Careserve - 0.8*Casto)./(0.6*TCALns+0.2*Casto)); SupCa(SupCa>1)=1; SupCa(SupCa<0)=0;
% end

SupMg = ((Mgreserve - 0.8*Mgsto)./(0.6*TMAGns+0.2*Mgsto)); SupMg(SupMg>1)=1; SupMg(SupMg<0)=0;
SupSi = ((Sireserve - 0.8*Sisto)./(0.6*TSILns+0.2*Sisto)); SupSi(SupSi>1)=1; SupSi(SupSi<0)=0;

%%%%
%%%%
TNIT = TNIT + LNIT; 
TPHO = TPHO + LPHO;
TPOT = TPOT + LPOT; 

TCAL = TCAL + LCAL; 
TMAG = TMAG + LMAG; 
TSIL = TSIL + LSIL; 

%%%%%
fll = Nreserve./Nsto_tm1; fll(fll>1)=1; fll(fll<0)=0; Ni = Nc/FiS + fll*(Nl-Nc)/FiS;
fll = Preserve./Psto_tm1; fll(fll>1)=1; fll(fll<0)=0; Phoi = Phoc/FiS + fll*(Phol-Phoc)/FiS;
fll = Kreserve./Ksto_tm1; fll(fll>1)=1; fll(fll<0)=0; Kpoti = Kpotc/FiS +  fll*(Kpotl-Kpotc)/FiS;

fll = Careserve./Casto_tm1; fll(fll>1)=1; fll(fll<0)=0; Cali = Calc/FiS +  fll*(Call-Calc)/FiS;
fll = Mgreserve./Mgsto_tm1; fll(fll>1)=1; fll(fll<0)=0; Magi = Magc/FiS +  fll*(Magl-Magc)/FiS;
fll = Sireserve./Sisto_tm1; fll(fll>1)=1; fll(fll<0)=0; Sili = Silc/FiS +  fll*(Sill-Silc)/FiS;

%%%%%%%

%%% Total Exported C-N-P-K
TexC = Sfr + Swm + Slf + Sr + sum(Rexmy) +RB(1)+RB(2)+RB(3)+RB(4)+RB(5)+RB(6)+RB(7); %% [gC/m^2 d]
%%%%
TexN = Slf./Nlf + rNc*Sfr/Nf + rNc*(1-ftransR*(rNc<=1))*Sr./Nr  + Swm/Nh + ...
rNc*RB(1)./Nl + RB(7)./Nlf  + RB(2)/Ns  + rNc*RB(3)/Nr + RB(4)/Ni + rNc*RB(5)/Nf  + RB(6)/Nh ; 
%%%
VarResN =  (RB(1)/Nl - rNc*RB(1)/Nl) +  (RB(3)/Nr -  rNc*RB(3)/Nr) + (RB(4)/Nc - RB(4)/Ni) + (RB(5)/Nf - rNc*RB(5)/Nf); %% % [gN/m^2 d] Variation on the reserve 
%%%%
TexP = Slf./Pholf + rPc*Sfr/Phof + rPc*(1-ftransR*(rPc<=1))*Sr./Phor  + Swm/Phoh + ...
rPc*RB(1)./Phol  +  RB(7)./Pholf  + RB(2)/Phos  + rPc*RB(3)/Phor + RB(4)/Phoi + rPc*RB(5)/Phof  +RB(6)/Phoh  ; 
%%%%
VarResP =  (RB(1)/Phol - rPc*RB(1)/Phol) +  (RB(3)/Phor -  rPc*RB(3)/Phor) + (RB(4)/Phoc - RB(4)/Phoi) + (RB(5)/Phof - rPc*RB(5)/Phof); %%% % [gP/m^2 d] Variation on the reserve 
%%%%
TexK = Slf./Kpotlf + rKc*Sfr/Kpotf + rKc*(1-ftransR*(rKc<=1))*Sr./Kpotr  + Swm/Kpoth + ...
rKc*RB(1)./Kpotl + RB(7)./Kpotlf + RB(2)/Kpots  + rKc*RB(3)/Kpotr + RB(4)/Kpoti + rKc*RB(5)/Kpotf  +RB(6)/Kpoth ; 
%%%%
VarResK =  (RB(1)/Kpotl - rKc*RB(1)/Kpotl) +  (RB(3)/Kpotr -  rKc*RB(3)/Kpotr) + (RB(4)/Kpotc - RB(4)/Kpoti) + (RB(5)/Kpotf - rKc*RB(5)/Kpotf); %%% % [gK/m^2 d] Variation on the reserve 
%%%%

TexCa = Slf./Callf + rCac*Sfr/Calf + rCac*(1-ftransR*(rCac<=1))*Sr./Calr  + Swm/Calh + ...
rCac*RB(1)./Call + RB(7)./Callf + RB(2)/Cals  + rCac*RB(3)/Calr + RB(4)/Cali + rCac*RB(5)/Calf  +RB(6)/Calh ; 
%%%%
VarResCa =  (RB(1)/Call - rCac*RB(1)/Call) +  (RB(3)/Calr -  rCac*RB(3)/Calr) + (RB(4)/Calc - RB(4)/Cali) + (RB(5)/Calf - rCac*RB(5)/Calf); %%% % [gCa/m^2 d] Variation on the reserve 
%%%%

TexMg = Slf./Maglf + rMgc*Sfr/Magf + rMgc*(1-ftransR*(rMgc<=1))*Sr./Magr  + Swm/Magh + ...
rMgc*RB(1)./Magl + RB(7)./Maglf + RB(2)/Mags  + rMgc*RB(3)/Magr + RB(4)/Magi + rMgc*RB(5)/Magf  +RB(6)/Magh ; 
%%%%
VarResMg =  (RB(1)/Magl - rMgc*RB(1)/Magl) +  (RB(3)/Magr -  rMgc*RB(3)/Magr) + (RB(4)/Magc - RB(4)/Magi) + (RB(5)/Magf - rMgc*RB(5)/Magf); %%% % [gMg/m^2 d] Variation on the reserve 
%%%%

TexSi = Slf./Sillf + rSic*Sfr/Silf + rSic*(1-ftransR*(rSic<=1))*Sr./Silr  + Swm/Silh + ...
rSic*RB(1)./Sill + RB(7)./Sillf + RB(2)/Sils  + rSic*RB(3)/Silr + RB(4)/Sili + rSic*RB(5)/Silf  +RB(6)/Silh ; 
%%%%
VarResSi =  (RB(1)/Sill - rSic*RB(1)/Sill) +  (RB(3)/Silr -  rSic*RB(3)/Silr) + (RB(4)/Silc - RB(4)/Sili) + (RB(5)/Silf - rSic*RB(5)/Silf); %%% % [gSi/m^2 d] Variation on the reserve 
%%%%

%%% Parton et all 1988 - Krinner et al 2005  Orwin et al 2011 
%Vmet_struct_CN = 5; %% Ratio between metabolic and structural C/N
%Lig_fr_x = ; % Lignin to Carbon fraction  [g Lig / gC]
%frac_to_metabolic = 0.85 - 0.018*(Nx*Lig_fr_x); 
%frac_to_struct = 1-  frac_to_metabolic; 
frac_to_metabolic_l = 0.85 - 0.018*(Nl*2*Lig_fr_l/rNc); 
frac_to_metabolic_l(frac_to_metabolic_l<0)=0;
frac_to_metabolic_fr = 0.85 - 0.018*(Nf*2*Lig_fr_fr/rNc); 
frac_to_metabolic_fr(frac_to_metabolic_fr<0)=0;
frac_to_metabolic_h = 0.85 - 0.018*(Nh*2*Lig_fr_h); 
frac_to_metabolic_h(frac_to_metabolic_h<0)=0;
frac_to_metabolic_r = 0.85 - 0.018*(Nr*2*Lig_fr_r/rNc); 
frac_to_metabolic_r(frac_to_metabolic_r<0)=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Managing the reserve issue
%%ManI=; %%% Management Indicator (0) Nothing (1) Fire (-1) Logging
%%% No Management
%%% Grazing, Grass cutting and Fruit harvesting are assumed to not modify
%%% nutrient reserves
NreserveM=Nreserve + VarResN;
PreserveM=Preserve + VarResP;
KreserveM=Kreserve + VarResK;

CareserveM=Careserve + VarResCa;
MgreserveM=Mgreserve + VarResMg;
SireserveM=Sireserve + VarResSi;

if OPT_SB == 1 %%% Only with soil-biogeochemistry on - otherwise reserves are depleted
    switch ManI
        case 1  %%% Fire
            AddN = AddN+ fire_eff*fab*NreserveM/dtd;
            AddP = AddP+ fire_eff*fab*PreserveM/dtd;
            AddK = AddK+ fire_eff*fab*KreserveM/dtd;
            
            AddCa = AddCa+ fire_eff*fab*CareserveM/dtd;
            AddMg = AddMg+ fire_eff*fab*MgreserveM/dtd;
            AddSi = AddSi+ fire_eff*fab*SireserveM/dtd;

            TexN = TexN  + AddN;
            TexP = TexP  + AddP;
            TexK = TexK  + AddK;

            TexCa = TexCa  + AddCa;
            TexMg = TexMg  + AddMg;
            TexSi = TexSi  + AddSi;

            NreserveM=(1-fire_eff*fab)*NreserveM;
            PreserveM=(1-fire_eff*fab)*PreserveM;
            KreserveM=(1-fire_eff*fab)*KreserveM;

            CareserveM=(1-fire_eff*fab)*CareserveM;
            MgreserveM=(1-fire_eff*fab)*MgreserveM;
            SireserveM=(1-fire_eff*fab)*SireserveM;

        case -1 %%% Logging
            AddN =  AddN+fract_log*NreserveM/dtd;
            AddP =  AddP+fract_log*PreserveM/dtd;
            AddK =  AddK+fract_log*KreserveM/dtd;

            AddCa =  AddCa+fract_log*CareserveM/dtd;
            AddMg =  AddMg+fract_log*MgreserveM/dtd;
            AddSi =  AddSi+fract_log*SireserveM/dtd;

            TexN = TexN  + AddN;
            TexP = TexP  + AddP;
            TexK = TexK  + AddK;

            TexCa = TexCa  + AddCa;
            TexMg = TexMg  + AddMg;
            TexSi = TexSi  + AddSi;

            NreserveM=(1-fract_log)*NreserveM;
            PreserveM=(1-fract_log)*PreserveM;
            KreserveM=(1-fract_log)*KreserveM;

            CareserveM=(1-fract_log)*CareserveM;
            MgreserveM=(1-fract_log)*MgreserveM;
            SireserveM=(1-fract_log)*SireserveM;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Fire_effect on C 
if ManI == 1 
    fract_left =0*fire_eff + (1-fire_eff)*fract_left; %%
    fract_left_fr =0*fire_eff + (1-fire_eff)*fract_left_fr; %% 
    fract_left_harAB =0*fire_eff +  (1-fire_eff)*fract_left_harAB; %% 
    fract_left_harBG =1*fire_eff + (1-fire_eff)*fract_left_harBG; %%
end
%%%% Carbon Inputs
%%%  Surface Litter 
IS_C_met_sur_lit = frac_to_metabolic_l*(Slf) +  frac_to_metabolic_fr*(Sfr); 
IS_C_str_sur_lit_lig = (1-frac_to_metabolic_l)*(Slf)*Lig_fr_l + (1-frac_to_metabolic_fr)*(Sfr)*Lig_fr_fr; 
IS_C_str_sur_lit_nlig = (1-frac_to_metabolic_l)*(Slf)*(1-Lig_fr_l) + (1-frac_to_metabolic_fr)*(Sfr)*(1-Lig_fr_fr); 
%%%% Surface Wood 
IS_C_wod_sur_lit_lig =  fab*Swm*Lig_fr_h ; 
IS_C_wod_sur_lit_nlig = fab*Swm*(1-Lig_fr_h); 
%%%% Suburface Litter 
IS_C_met_ssr_lit =  Rexmy(1) +  frac_to_metabolic_r*Sr + (frac_to_metabolic_h)*fbe*Swm;   
IS_C_str_ssr_lit_lig = (1-frac_to_metabolic_r)*Sr*Lig_fr_r + (1-frac_to_metabolic_h)*fbe*Swm*Lig_fr_h; 
IS_C_str_ssr_lit_nlig = Rexmy(3) + (1-frac_to_metabolic_r)*Sr*(1-Lig_fr_r) + (1-frac_to_metabolic_h)*fbe*Swm*(1-Lig_fr_h); 
%%%
IS_C_myc = Rexmy(2);  
%%%%%%% Litter due to management 
%%%  Surface Litter 
IS_C_met_sur_lit2 =   frac_to_metabolic_l*(fract_left*(RB(1)+RB(7))) +  frac_to_metabolic_fr*(fract_left_fr*RB(5)); 
IS_C_str_sur_lit_lig2 =  (1-frac_to_metabolic_l)*(fract_left*(RB(1)+RB(7)))*Lig_fr_l + (1-frac_to_metabolic_fr)*(fract_left_fr*RB(5))*Lig_fr_fr; 
IS_C_str_sur_lit_nlig2 = (1-frac_to_metabolic_l)*(fract_left*(RB(1)+RB(7)))*(1-Lig_fr_l) + (1-frac_to_metabolic_fr)*(fract_left_fr*RB(5))*(1-Lig_fr_fr); 
%%%% Surface Wood 
IS_C_wod_sur_lit_lig2 =   Lig_fr_h*fab*fract_left_harAB*(RB(2)+RB(4)+RB(6)); 
IS_C_wod_sur_lit_nlig2 =  (1-Lig_fr_h)*fab*fract_left_harAB*(RB(2)+RB(4)+RB(6)); 
%%%% Suburface Litter 
IS_C_met_ssr_lit2 =  frac_to_metabolic_r*fract_left_harBG*RB(3) + (frac_to_metabolic_h)*fbe*fract_left_harBG*(RB(2)+RB(4)+RB(6));
IS_C_str_ssr_lit_lig2 = (1-frac_to_metabolic_r)*fract_left_harBG*RB(3)*Lig_fr_r + (1-frac_to_metabolic_h)*fbe*fract_left_harBG*(RB(2)+RB(4)+RB(6))*Lig_fr_h;
IS_C_str_ssr_lit_nlig2 = (1-frac_to_metabolic_r)*fract_left_harBG*RB(3)*(1-Lig_fr_r) + (1-frac_to_metabolic_h)*fbe*fract_left_harBG*(RB(2)+RB(4)+RB(6))*(1-Lig_fr_h);
%%%%%% Fire_effect on N
if  ManI == 1  %% Hyp. %  a fraction of N left in the field
    fract_left =Mpar.fract_left; %% Fraction_left of leaves  and dead leaves
    fract_left_fr =Mpar.fract_left_fr; %% Fraction left of fruits
    fract_left_harAB =Mpar.fract_left_AB; %% Fraction left of harvested wood aboveground
    fract_left_harBG =Mpar.fract_left_BG; %% Fraction left of harvested wood belowground
    %%%
    fract_left =funb_nit*fire_eff + (1-fire_eff)*fract_left; %%
    fract_left_fr =funb_nit*fire_eff + (1-fire_eff)*fract_left_fr; %%
    fract_left_harAB =funb_nit*fire_eff +  (1-fire_eff)*fract_left_harAB; %%
    fract_left_harBG =1*fire_eff + (1-fire_eff)*fract_left_harBG; %%
end
%%%
%%%% Nitrogen Inputs 
IS_N_sur_lit =  Slf./Nlf  + rNc*Sfr/Nf ; % [gN/m^2 day]   
IS_N_wod_lit =  fab*Swm/Nh ; 
IS_N_ssr_lit =  rNc*(1-ftransR*(rNc<=1))*Sr./Nr + fbe*Swm/Nh ;  % [gN/m^2 day]
%%%% Managed 
IS_N_sur_lit2 =  rNc*(fract_left*RB(1))./Nl + fract_left*RB(7)/Nlf  + rNc*fract_left_fr*RB(5)/Nf + fab*AddN ; % [gN/m^2 day]   
IS_N_wod_lit2 =  fract_left_harAB*fab*(RB(2)/Ns + RB(4)/Ni + RB(6)/Nh) ; 
IS_N_ssr_lit2 =  rNc*fract_left_harBG*RB(3)/Nr + fract_left_harBG*fbe*(RB(2)/Ns+ RB(4)/Ni +RB(6)/Nh) +fbe*AddN;  % [gN/m^2 day]
%%%
%%%%%% Fire_effect on P and K 
if ManI == 1
    fract_left =Mpar.fract_left; %% Fraction_left of leaves  and dead leaves
    fract_left_fr =Mpar.fract_left_fr; %% Fraction left of fruits
    fract_left_harAB =Mpar.fract_left_AB; %% Fraction left of harvested wood aboveground
    fract_left_harBG =Mpar.fract_left_BG; %% Fraction left of harvested wood belowground
    %%%%
    fract_left =1*fire_eff + (1-fire_eff)*fract_left; %%
    fract_left_fr =1*fire_eff + (1-fire_eff)*fract_left_fr; %%
    fract_left_harAB =1*fire_eff +  (1-fire_eff)*fract_left_harAB; %%
    fract_left_harBG =1*fire_eff + (1-fire_eff)*fract_left_harBG; %%
end
%%%% Phosoporus Inputs 
IS_P_sur_lit = Slf./Pholf  + rPc*Sfr/Phof  ; % [gP/m^2 day]
IS_P_wod_lit = fab*Swm/Phoh ; 
IS_P_ssr_lit = rPc*(1-ftransR*(rPc<=1))*Sr./Phor + fbe*Swm/Phoh;   % [gP/m^2 day]
%%%% Managed
IS_P_sur_lit2 =  rPc*(fract_left*RB(1))./Phol + fract_left*RB(7)/Pholf +  rPc*fract_left_fr*RB(5)/Phof + fab*AddP ;  % [gP/m^2 day] 
IS_P_wod_lit2 =  fract_left_harAB*fab*(RB(2)/Phos+ RB(4)/Phoi + RB(6)/Phoh); 
IS_P_ssr_lit2 =  rPc*fract_left_harBG*RB(3)/Phor + fract_left_harBG*fbe*(RB(2)/Phos+ RB(4)/Phoi +RB(6)/Phoh) +fbe*AddP;   % [gP/m^2 day]
%%%% Potassium Inputs 
IS_K_sur_lit =  Slf./Kpotlf  + rKc*Sfr/Kpotf ; % [gK/m^2 day]
IS_K_wod_lit = fab*Swm/Kpoth ; 
IS_K_ssr_lit = rKc*(1-ftransR*(rKc<=1))*Sr./Kpotr + fbe*Swm/Kpoth;   % [gK/m^2 day]
%%%%% Managed 
IS_K_sur_lit2 =  rKc*(fract_left*RB(1))./Kpotl + fract_left*RB(7)/Kpotlf + rKc*fract_left_fr*RB(5)/Kpotf +  fab*AddK;  % [gK/m^2 day]
IS_K_wod_lit2 =  fract_left_harAB*fab*(RB(2)/Kpots + RB(4)/Kpoti + RB(6)/Kpoth); 
IS_K_ssr_lit2 =  rKc*fract_left_harBG*RB(3)/Kpotr + fract_left_harBG*fbe*(RB(2)/Kpots+ RB(4)/Kpoti +RB(6)/Kpoth) +fbe*AddK;   % [gK/m^2 day]

%%%% Calcium Inputs 
IS_Ca_sur_lit =  Slf./Callf  + rCac*Sfr/Calf ; % [gCa/m^2 day]
IS_Ca_wod_lit = fab*Swm/Calh ; 
IS_Ca_ssr_lit = rCac*(1-ftransR*(rCac<=1))*Sr./Calr + fbe*Swm/Calh;   % [gCa/m^2 day]
%%%%% Managed 
IS_Ca_sur_lit2 =  rCac*(fract_left*RB(1))./Call + fract_left*RB(7)/Callf + rCac*fract_left_fr*RB(5)/Calf +  fab*AddCa;  % [gCa/m^2 day]
IS_Ca_wod_lit2 =  fract_left_harAB*fab*(RB(2)/Cals + RB(4)/Cali + RB(6)/Calh); 
IS_Ca_ssr_lit2 =  rCac*fract_left_harBG*RB(3)/Calr + fract_left_harBG*fbe*(RB(2)/Cals+ RB(4)/Cali +RB(6)/Calh) +fbe*AddCa;   % [gCa/m^2 day]

%%%% Magnesium Inputs 
IS_Mg_sur_lit =  Slf./Maglf  + rMgc*Sfr/Magf ; % [gMg/m^2 day]
IS_Mg_wod_lit = fab*Swm/Magh ; 
IS_Mg_ssr_lit = rMgc*(1-ftransR*(rMgc<=1))*Sr./Magr + fbe*Swm/Magh;   % [gMg/m^2 day]
%%%%% Managed 
IS_Mg_sur_lit2 =  rMgc*(fract_left*RB(1))./Magl + fract_left*RB(7)/Maglf + rMgc*fract_left_fr*RB(5)/Magf +  fab*AddMg;  % [gMg/m^2 day]
IS_Mg_wod_lit2 =  fract_left_harAB*fab*(RB(2)/Mags + RB(4)/Magi + RB(6)/Magh); 
IS_Mg_ssr_lit2 =  rMgc*fract_left_harBG*RB(3)/Magr + fract_left_harBG*fbe*(RB(2)/Mags+ RB(4)/Magi +RB(6)/Magh) +fbe*AddMg;   % [gMg/m^2 day]

%%%% Silicon Inputs 
IS_Si_sur_lit =  Slf./Sillf  + rSic*Sfr/Silf ; % [gSi/m^2 day]
IS_Si_wod_lit = fab*Swm/Silh ; 
IS_Si_ssr_lit = rSic*(1-ftransR*(rSic<=1))*Sr./Silr + fbe*Swm/Silh;   % [gSi/m^2 day]
%%%%% Managed 
IS_Si_sur_lit2 =  rSic*(fract_left*RB(1))./Sill + fract_left*RB(7)/Sillf + rSic*fract_left_fr*RB(5)/Silf +  fab*AddSi;  % [gSi/m^2 day]
IS_Si_wod_lit2 =  fract_left_harAB*fab*(RB(2)/Sils + RB(4)/Sili + RB(6)/Silh); 
IS_Si_ssr_lit2 =  rSic*fract_left_harBG*RB(3)/Silr + fract_left_harBG*fbe*(RB(2)/Sils+ RB(4)/Sili +RB(6)/Silh) +fbe*AddSi;   % [gSi/m^2 day]



%%%%%%%%%%%%%%% Re-summing all litter inputs 
IS_C_met_sur_lit =  IS_C_met_sur_lit +IS_C_met_sur_lit2; 
IS_C_str_sur_lit_lig =  IS_C_str_sur_lit_lig  +IS_C_str_sur_lit_lig2; 
IS_C_str_sur_lit_nlig  = IS_C_str_sur_lit_nlig + IS_C_str_sur_lit_nlig2; 
IS_C_wod_sur_lit_lig =   IS_C_wod_sur_lit_lig + IS_C_wod_sur_lit_lig2; 
IS_C_wod_sur_lit_nlig =  IS_C_wod_sur_lit_nlig + IS_C_wod_sur_lit_nlig2; 
IS_C_met_ssr_lit =  IS_C_met_ssr_lit + IS_C_met_ssr_lit2; 
IS_C_str_ssr_lit_lig =  IS_C_str_ssr_lit_lig + IS_C_str_ssr_lit_lig2; 
IS_C_str_ssr_lit_nlig =  IS_C_str_ssr_lit_nlig + IS_C_str_ssr_lit_nlig2; 
IS_N_sur_lit =   IS_N_sur_lit +  IS_N_sur_lit2; 
IS_N_wod_lit =   IS_N_wod_lit +  IS_N_wod_lit2; 
IS_N_ssr_lit =   IS_N_ssr_lit + IS_N_ssr_lit2 ;
IS_P_sur_lit =  IS_P_sur_lit+ IS_P_sur_lit2; 
IS_P_wod_lit = IS_P_wod_lit + IS_P_wod_lit2; 
IS_P_ssr_lit =  IS_P_ssr_lit + IS_P_ssr_lit2; 
IS_K_sur_lit =   IS_K_sur_lit + IS_K_sur_lit2; 
IS_K_wod_lit =  IS_K_wod_lit + IS_K_wod_lit2; 
IS_K_ssr_lit =  IS_K_ssr_lit + IS_K_ssr_lit2; 

IS_Ca_sur_lit =   IS_Ca_sur_lit + IS_Ca_sur_lit2; 
IS_Ca_wod_lit =  IS_Ca_wod_lit + IS_Ca_wod_lit2; 
IS_Ca_ssr_lit =  IS_Ca_ssr_lit + IS_Ca_ssr_lit2; 

IS_Mg_sur_lit =   IS_Mg_sur_lit + IS_Mg_sur_lit2; 
IS_Mg_wod_lit =  IS_Mg_wod_lit + IS_Mg_wod_lit2; 
IS_Mg_ssr_lit =  IS_Mg_ssr_lit + IS_Mg_ssr_lit2; 

IS_Si_sur_lit =   IS_Si_sur_lit + IS_Si_sur_lit2; 
IS_Si_wod_lit =  IS_Si_wod_lit + IS_Si_wod_lit2; 
IS_Si_ssr_lit =  IS_Si_ssr_lit + IS_Si_ssr_lit2;

%%%%%%%%%%%%%%%
IS=[IS_C_met_sur_lit;IS_C_str_sur_lit_lig;IS_C_str_sur_lit_nlig;
    IS_C_wod_sur_lit_lig;IS_C_wod_sur_lit_nlig;
    IS_C_met_ssr_lit;IS_C_str_ssr_lit_lig;IS_C_str_ssr_lit_nlig;
    IS_C_myc;
    IS_N_sur_lit;IS_N_wod_lit;IS_N_ssr_lit;
    IS_P_sur_lit;IS_P_wod_lit;IS_P_ssr_lit;
    IS_K_sur_lit;IS_K_wod_lit;IS_K_ssr_lit;
    IS_Ca_sur_lit;IS_Ca_wod_lit;IS_Ca_ssr_lit;
    IS_Mg_sur_lit;IS_Mg_wod_lit;IS_Mg_ssr_lit;
    IS_Si_sur_lit;IS_Si_wod_lit;IS_Si_ssr_lit];
%%%%%%
%%%%%%%%
if isreal(sum(IS))==0 || isnan(sum(IS)) == 1
    disp('NaN in Litter export')
    return
end
return 
