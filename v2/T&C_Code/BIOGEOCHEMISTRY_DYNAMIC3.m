%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction BIOGEOCHEMISTRY_DYNAMIC   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% References Manzoni and Porporato, 2007; 2009;  Dickinson et al 2002;
%%% Xu-Ri and Prentice 2008; Li et al., 1992; 2000
%%% Porporato et al., 2003; Zahele and Friend 2010;
%%% Orwin et al 2011 ; Wang et al 2012
function[dB,R_litter,R_microbe,R_litter_sur,R_ew,VOL,N2flx,Min_N,Min_P,...
    R_bacteria,RmycAM,RmycEM,Prod_B,Prod_F,N2Oflx]= BIOGEOCHEMISTRY_DYNAMIC3(t,B,ZBIOG,rsd,...
    ISOIL,Ts,Ta,Psi_s,PH,Se,Se_fc,...
    FertN,DepN,BfixN,FertP,DepP,FertK,DepK,...
    FertCa,DepCa,FertMg,DepMg,FertSi,DepSi,...
    NH4_Uptake,NO3_Uptake,P_Uptake,K_Uptake,...
    Ca_Uptake,Mg_Uptake,Si_Uptake,...
    LEAK_DOC,LEAK_NH4,LEAK_NO3,LEAK_P,LEAK_K,LEAK_DON,LEAK_DOP,...
    LEAK_Ca,LEAK_Mg,LEAK_Si,...
    Tup_P,Tup_K,...
    Tup_Ca,Tup_Mg,Tup_Si, ...
    ExEM,Pcla,Psan,BiogeoPar,SC_par,IMAN,opt_cons_CUE,opt_ERW,L_pore)
%%% OUTPUT
%VOL ;
%N2flx ;
%R_microbe
%R_litter
%%%% INPUT
dtd=1; %% [day]
Psi_s = Psi_s*1000; %% [kPa] Soil Water Potential
%Lk = [mm/d]
%V=[mm]
%%%%%%%%%%%
%ZBIOG = 0.25; %%[m]
%rsd =1300; %% density dry soil [kg/m^3]
rsd=rsd*1000; %% [gsoil /m3]
%%%% Carbon Inputs %% [gC/gsoil day]
IS.C_met_sur_lit =  ISOIL(1)/(ZBIOG*rsd) + IMAN(1)/(ZBIOG*rsd);
IS.C_str_sur_lit_lig = ISOIL(2)/(ZBIOG*rsd) + IMAN(2)/(ZBIOG*rsd);
IS.C_str_sur_lit_nlig = ISOIL(3)/(ZBIOG*rsd)+ IMAN(3)/(ZBIOG*rsd);
%%%%
IS.C_wod_sur_lit_lig =  ISOIL(4)/(ZBIOG*rsd);
IS.C_wod_sur_lit_nlig =  ISOIL(5)/(ZBIOG*rsd);
%%%%
IS.C_met_ssr_lit =  ISOIL(6)/(ZBIOG*rsd);
IS.C_str_ssr_lit_lig = ISOIL(7)/(ZBIOG*rsd);
IS.C_str_ssr_lit_nlig = ISOIL(8)/(ZBIOG*rsd);
%%%
IS.C_myc = ISOIL(9)/(ZBIOG*rsd);
%%%% Nitrogen Inputs %  [gN/gsoil day]
IS.N_sur_lit = ISOIL(10)/(ZBIOG*rsd)+ IMAN(4)/(ZBIOG*rsd);
IS.N_wod_lit = ISOIL(11)/(ZBIOG*rsd);
IS.N_ssr_lit = ISOIL(12)/(ZBIOG*rsd);
%%%% Phosoporus Inputs   [gP/gsoil day]
IS.P_sur_lit = ISOIL(13)/(ZBIOG*rsd) + IMAN(5)/(ZBIOG*rsd);
IS.P_wod_lit = ISOIL(14)/(ZBIOG*rsd);
IS.P_ssr_lit = ISOIL(15)/(ZBIOG*rsd);
%%%% Potassium Inputs   [gK/gsoil day]
IS.K_sur_lit =  ISOIL(16)/(ZBIOG*rsd) + IMAN(6)/(ZBIOG*rsd) ;
IS.K_wod_lit =  ISOIL(17)/(ZBIOG*rsd);
IS.K_ssr_lit = ISOIL(18)/(ZBIOG*rsd);

%%%% Calcium Inputs   [gCa/gsoil day]
IS.Ca_sur_lit =  ISOIL(19)/(ZBIOG*rsd) + IMAN(7)/(ZBIOG*rsd) ;
IS.Ca_wod_lit =  ISOIL(20)/(ZBIOG*rsd);
IS.Ca_ssr_lit = ISOIL(21)/(ZBIOG*rsd);
%%%% Magnesium Inputs   [gMg/gsoil day]
IS.Mg_sur_lit =  ISOIL(22)/(ZBIOG*rsd) + IMAN(8)/(ZBIOG*rsd) ;
IS.Mg_wod_lit =  ISOIL(23)/(ZBIOG*rsd);
IS.Mg_ssr_lit = ISOIL(24)/(ZBIOG*rsd);
%%%% Silicate Inputs   [gSi/gsoil day]
IS.Si_sur_lit =  ISOIL(25)/(ZBIOG*rsd) + IMAN(9)/(ZBIOG*rsd) ;
IS.Si_wod_lit =  ISOIL(26)/(ZBIOG*rsd);
IS.Si_ssr_lit = ISOIL(27)/(ZBIOG*rsd);

%%%%%%%%%%%%%%%%%
%%%% Tectonic Uplift %% [gX/gsoil d]
Tup_P= Tup_P/(ZBIOG*rsd);
Tup_K= Tup_K/(ZBIOG*rsd);

Tup_Ca= Tup_Ca/(ZBIOG*rsd);
Tup_Mg= Tup_Mg/(ZBIOG*rsd);
Tup_Si= Tup_Si/(ZBIOG*rsd);

%%%%%%%% Plant and Mychorrizal uptake  %% [gX/gsoil d]
NH4_Uptake= NH4_Uptake/(ZBIOG*rsd);
NO3_Uptake=NO3_Uptake/(ZBIOG*rsd);
P_Uptake=P_Uptake/(ZBIOG*rsd);
K_Uptake=K_Uptake/(ZBIOG*rsd);

Ca_Uptake=Ca_Uptake/(ZBIOG*rsd);
Mg_Uptake=Mg_Uptake/(ZBIOG*rsd);
Si_Uptake=Si_Uptake/(ZBIOG*rsd);

%%%%%%%%%%%%%% Leakage of Nutrient Compounds
LEAK_NH4 = LEAK_NH4/(ZBIOG*rsd) ; %% [gN/gsoil d]
LEAK_NO3 = LEAK_NO3/(ZBIOG*rsd) ; %% [gN/gsoil d]
LEAK_P = LEAK_P/(ZBIOG*rsd);   %% [gP/gsoil d]
LEAK_K = LEAK_K/(ZBIOG*rsd);  %% [gK/gsoil d]
LEAK_DOC = LEAK_DOC/(ZBIOG*rsd);
LEAK_DOP = LEAK_DOP/(ZBIOG*rsd);
LEAK_DON = LEAK_DON/(ZBIOG*rsd);

LEAK_Ca = LEAK_Ca/(ZBIOG*rsd);  %% [gCa/gsoil d]
LEAK_Mg = LEAK_Mg/(ZBIOG*rsd);  %% [gMg/gsoil d]
LEAK_Si = LEAK_Si/(ZBIOG*rsd);  %% [gSi/gsoil d]

%%% METABOLIC (fast)  Proteins  Carbohydrate Reserve (Starch Fructants)
%%% Chlorophyll
%%% STRUCTURAL (Cellulose, hemicelllulose, pectin Lignin Tannin Polyphenolic
%%% Lipids Cutine Suberine
%%% OUTPUT
%%% dB [gX /m^2 d]
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% PARAMETERS
%%%%%%% CARBON POOL %%%%%%%%%%%
%%% B1 Above-ground Litter Metabolic
%%% B2 Above-ground Litter Structural - Cellulose/Hemicellulose
%%% B3 Above-ground Litter Structura - Lignin
%%% B4 Above-ground Woody  - Cellulose/Hemicellulose
%%% B5 Above-ground Woody - Lignin
%%% B6 Below-ground Litter Metabolic
%%% B7 Below-ground Litter Structural - Cellulose/Hemicellulose
%%% B8 Below-ground Litter Structura - Lignin
%%% B9  SOM-POC- lignin
%%% B10 SOM-POC -Cellulose/Hemicellulose
%%% B11 SOM-MOC
%%% B12 DOC - for bacteria
%%% B13 DOC - for fungi
%%% B14 Enzyme for decomposition of POC-Bact
%%% B15 Enzyme for decomposition of POC-Fung
%%% B16 Enzyme for decomposition of MOC-Bact
%%% B17 Enzyme for decomposition of MOC-Fung
%%% B18 Bacteria pool
%%% B19 Fungi saprotrophic
%%% B20 AM-Mycorrhizal - C
%%% B21 EM-Mycorrhizal - C
%%% B22 Earthworms - C
%%%%%% NITROGEN POOL
%%% B23 Nitrogen Above-ground Litter
%%% B24 Nitrogen Above-ground Woody
%%% B25 Nitrogen Below-ground Litter
%%% B26 Nitrogen SOM
%%% B27 Nitrogen Bacteria
%%% B28 Nitrogen Fungi
%%% B29 AM Mycorrhizal - N
%%% B30 EM Mycorrhizal - N
%%% B31 Nitrogen Ion Ammonium NH4+
%%% B32 Nitrogen Nitrate NO3-
%%% B33 DON
%%% B34 Earthworms - N
%%%%%% PHOSPHORUS POOL
%%% B35 phosphorus Above-ground Litter
%%% B36 phosphorus Above-ground Woody
%%% B37 phosphorus Below-ground Litter
%%% B38 phosphorus SOM
%%% B39 phosphorus Bacteria
%%% B40 phosphorus Fungi
%%% B41 AM - Mycorrhizal - P
%%% B42 EM - Mycorrhizal - P
%%% B43 phosphorus Mineral
%%% B44 phosphorus primary
%%% B45 phosphorus secondary
%%% B46 phosphorus occluded
%%% B47 DOP
%%%%%% POTASSIUM POOL
%%% B48 Potassium Above-ground Litter -
%%% B49 Potassium  Above-ground Woody -
%%% B50 Potassium  Below-ground Litter -
%%% B51 Potassium SOM
%%% B52 Potassium  Mineral  solution -
%%% B53 Potassium  exchangeable
%%% B54 Potassium fixed or non-exchangeable
%%% B55  Potassium in the lattice of certain primary minerals
%%%%%%%%
%%%%%%%%%%%%
%%%%%% CALCIUM POOL
%%% B56 Calcium Above-ground Litter -
%%% B57 Calcium  Above-ground Woody -
%%% B58 Calcium  Below-ground Litter -
%%% B59 Calcium SOM
%%% B60 Calcium  Mineral  solution -
%%% B61 Calcium  exchangeable
%%% B62 Calcium fixed or non-exchangeable
%%% B63  Calcium in the lattice of certain primary minerals
%%%%%% MAGNESIUM POOL
%%% B64 Magnesium Above-ground Litter -
%%% B65 Magnesium  Above-ground Woody -
%%% B66 Magnesium  Below-ground Litter -
%%% B67 Magnesium SOM
%%% B68 Magnesium  Mineral  solution -
%%% B69 Magnesium  exchangeable
%%% B70 Magnesium fixed or non-exchangeable
%%% B71  Magnesium in the lattice of certain primary minerals
%%%%%% SILICATE POOL
%%% B72 Silicate Above-ground Litter -
%%% B73 Silicate  Above-ground Woody -
%%% B74 Silicate  Below-ground Litter -
%%% B75 Silicate SOM
%%% B76 Silicate  Mineral  solution -
%%% B77 Silicate  exchangeable
%%% B78 Silicate fixed or non-exchangeable
%%% B79  Silicate in the lattice of certain primary minerals
%%%%%%%%
%%%%%%%%%%%%
B= B/(ZBIOG*rsd); %% [gX /g soil]
%%%%%%%% subtract quantities uptaken for internal computations
B(31)= B(31)- NH4_Uptake;%  Nitrogen Ione Ammonium NH4+
B(32)=B(32)-NO3_Uptake;%%   Nitrogen Nitrate NO3-
B(43)=B(43)-P_Uptake;%  phosphorus Mineral
B(52)=B(52)-K_Uptake;%% Potassium  Mineral  solution

B(60)=B(60)-Ca_Uptake;%% Calcium  Mineral  solution
B(68)=B(68)-Mg_Uptake;%% Magnesium  Mineral  solution
B(76)=B(76)-Si_Uptake;%% Silicate  Mineral  solution

%%% Retain physical inventory after plant uptake
%% Pools arrive as P = Ptm1 + dP from the previous step. The donor
%% limiters keep them nonnegative, but a pool driven to exactly zero can
%% land a rounding step below it. Tolerate that and clamp; a meaningful
%% negative still aborts, and reports how negative it was.
pool_tol = 1e-12*max(1,max(abs(B)));
if any(~isfinite(B)) || any(B < -pool_tol)
    error('BIOGEOCHEMISTRY_DYNAMIC3:InvalidPoolState', ...
        ['Biogeochemistry pools must be finite and nonnegative after ' ...
         'uptake; most negative is %.6g against a %.3g tolerance.'], ...
        min(B),pool_tol);
end
B = max(0,B);
B_physical = B;
%%%% Decay Constant on a 40�C basis (Kirschbaum 2002)
k_met_sur = BiogeoPar.k_met_sur; %%[1/day]
k_str_sur = BiogeoPar.k_str_sur; %% %[1/day]
k_met_ssr = BiogeoPar.k_met_ssr;
k_str_ssr = BiogeoPar.k_str_ssr;
k_wod_sur = BiogeoPar.k_wod_sur; %% [1/92 - 1/460]
%%%%%%%%%%%%%%%%%%%%
%%%%%% Respiration Coefficient for Litter [-]
%rr_met_sur=0.55;
%rr_str_sur=0.45;
%rr_wod_sur=0.55;
%rr_met_ssr=0.55;
%rr_str_ssr=0.55;
%%%%%%% Critical Ratios
R_CN_bac = BiogeoPar.R_CN_bac ; %%[C/N]  8.6
R_CN_fun = BiogeoPar.R_CN_fun ;
R_CN_mycAM = BiogeoPar.R_CN_myc;
R_CN_mycEM = BiogeoPar.R_CN_myc;
R_CN_ew = BiogeoPar.R_CN_ew;
R_CP_bac= BiogeoPar.R_CP_bac; %%[C/P] 60
R_CP_fun= BiogeoPar.R_CP_fun;
R_CP_mycAM = BiogeoPar.R_CP_myc;
R_CP_mycEM = BiogeoPar.R_CP_myc;
%%%%%%%%%%%%%%%%%%
%ExEM = [0-1] Fraction of EM vs AM
%%% Parameters
gd =BiogeoPar.gd; % [0.2-0.8] fraction of dead B allocated to D
fd =BiogeoPar.fd; % [0.2-0.8]  fraction of decomposed P allocated to D
%Ecb =0.47; % carbon use efficiency [-]
%Ecf =0.47;
Vdb = BiogeoPar.Vdb; %maximum specific uptake rate of D for growth of B [gC-D/gC-Bday]
Vdf = BiogeoPar.Vdf;
Kdb =BiogeoPar.Kdb ; % [0.00014 - 0.00038] half-saturation constant of uptake of D for growth of B [ g C/g soil]
Kdf =BiogeoPar.Kdf ;
mrb = SC_par(1)*BiogeoPar.mrb; % [ 1.20*1e-3 - 30*1e-3] specific maintenance factor or rate [ g C /g C day]
mrf = SC_par(2)*BiogeoPar.mrf;
mram = SC_par(2)*BiogeoPar.mram;
mrem = SC_par(2)*BiogeoPar.mrem;
Vpc =  BiogeoPar.Vpc; % maximum specific decomposition rate for P by EP [gC-P/gC-EP day]
Vpl =  BiogeoPar.Vpl;
Kpl =BiogeoPar.Kpl; % half-saturation constant for decomposition of P [gC / gsoil]
Kpc =BiogeoPar.Kpc ;
Vm = BiogeoPar.Vm;  %maximum specific decomposition rate for M by [gC-M /gC-EM day]
Km =BiogeoPar.Km; %half-saturation constant for decomposition of [ gC /gsoil];
rem = BiogeoPar.rem ; % [0.0075 - 0.075] turnover rate of EM [g C /g C day]
rep =BiogeoPar.rep ; %  [0.0075 - 0.075] turnover rate of EP [g C /g C day]
pepb = SC_par(3)*BiogeoPar.pepb; % [0.0031 - 0.031] fraction of mR for production of EP
pemb =SC_par(3)*BiogeoPar.pemb;% [0.0031 - 0.031] fraction of mR for production of EM
pepf = SC_par(3)*BiogeoPar.pepf; %
pemf =SC_par(3)*BiogeoPar.pemf;%
pepem = SC_par(3)*BiogeoPar.pepem; %
pemem =SC_par(3)*BiogeoPar.pemem;%
%%%%%%%%%% Leaching Coefficients
lambda_C = BiogeoPar.lambda_C;
lambda_N = BiogeoPar.lambda_N;
lambda_P = BiogeoPar.lambda_P;
lambda_K = BiogeoPar.lambda_K;

lambda_Ca = BiogeoPar.lambda_Ca;
lambda_Mg = BiogeoPar.lambda_Mg;
lambda_Si = BiogeoPar.lambda_Si;

f_org_lea = BiogeoPar.f_org_lea; %% Fraction of Organic leaching
%%%% Earthworms parameter
%kse = BiogeoPar.kse ; % half-saturation constant for consumption of litter
VmaxE=  SC_par(4)*BiogeoPar.VmaxE ;%%% [gC / gC-food d] Curry and Schmidt 2007
rmanEa = BiogeoPar.rmanEa ; % %% [gC / gC-EW d]
rmanEp = BiogeoPar.rmanEp ; % %% [gC / gC-EW d]
dew = BiogeoPar.dew ; % %% [gC / gC-EW d]
fabsE = BiogeoPar.fabsE ; %% [-]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Soil Moisture regulation %%%%%%%%%%%%%
%%%% Moyano et al 2013; Manzoni et al 2012
%%% Psi_s in [kPa]
if Psi_s == 0
    fSM_litter =0 ; fSM_microbe = 0; Aew=0;
else
    Psi_s(Psi_s>-10)=-10.1;
    Psi_opt = -3;  Psi_th = -15800;  alpha = 1.47;
    fSM_microbe  = 1 - ((log(Psi_s)-log(Psi_opt))./(log(Psi_th)-log(Psi_opt))).^alpha;
    Psi_opt = -10;  Psi_th = -28800;  alpha = 1.0;
    fSM_litter =  1 - ((log(Psi_s)-log(Psi_opt))./(log(Psi_th)-log(Psi_opt))).^alpha;
    %%%%
    fSM_microbe(fSM_microbe <0)=0;
    fSM_litter(fSM_litter<0)=0;
    %%%%
    Aew = 0.000008575*exp(11.67*Se./Se_fc)*(Se>0.2);
    Aew(Se>0.9)=0.002729*Se(Se>0.9).^(-56.03); Aew(Aew>1)=1;
end
%if Se <= Sfc
%    fSM = Se/Sfc;
%else
%    fSM =  Sfc/Se;  %% [-] Porporato et al 2003
%end

%%% Wang et al 2012b - varying sensitivity in pH=4.5-7.5
PHopt=6; Psen=2;
fPH= exp(-((PH-PHopt)/Psen)^2) ;

% %%% Orwin et al 2011  [min 2 4.5 optm 7.5 12 max]
if PH<=4.5
    fPH = 0.4*PH - 0.8;
else
    if PH >=7.5
        fPH = -0.2222*PH + 2.666;
    else
        fPH=1;
    end
end

fPH(fPH>1)=1; fPH(fPH<0)=0;
%%%%%%%%% Clay Dependence of half-saturation and fraction of dead B (Wieder et al. 2015, GMD, Six et al 2002 PLSO)
fClay =   0.6940+1.36*(Pcla);
%fClay = (2.582*exp(-2*sqrt(Pcla))).^(-1);
fClay2 = (1.2*exp(-0.8*(Pcla)));
%%%%%%%%%%%%%%%%%%%%%% Temperature Effects Kinetiks
R= 8.314*1e-3 ; %[kJ mol-1 K-1]
%V=Vref*exp(-Ea/R*(1/(Ts+273.15) - 1/(Tref+273.15)));
Tref=BiogeoPar.T1;
%%%%%%%
Vpc = Vpc*exp(-37/R*(1/(Ts+273.15) - 1/(Tref+273.15)))*fSM_microbe*fPH;
Vpl = Vpl*exp(-53/R*(1/(Ts+273.15) - 1/(Tref+273.15)))*fSM_microbe*fPH;
Vm = Vm*exp(-47/R*(1/(Ts+273.15) - 1/(Tref+273.15)))*fSM_microbe*fPH;
Vdb = Vdb*exp(-47/R*(1/(Ts+273.15) - 1/(Tref+273.15)))*fSM_microbe*fPH;
Vdf = Vdf*exp(-47/R*(1/(Ts+273.15) - 1/(Tref+273.15)))*fSM_microbe*fPH;
Kpl = Kpl*exp(-30/R*(1/(Ts+273.15) - 1/(Tref+273.15)))*fClay;
Kpc = Kpc*exp(-30/R*(1/(Ts+273.15) - 1/(Tref+273.15)))*fClay;
Km = Km*exp(-30/R*(1/(Ts+273.15) - 1/(Tref+273.15)))*fClay;
Kdb= Kdb*exp(-30/R*(1/(Ts+273.15) - 1/(Tref+273.15)));
Kdf= Kdf*exp(-30/R*(1/(Ts+273.15) - 1/(Tref+273.15)));
mrb=mrb*exp(-20/R*(1/(Ts+273.15) - 1/(Tref+273.15)));
mrf=mrf*exp(-20/R*(1/(Ts+273.15) - 1/(Tref+273.15)));
mram=mram*exp(-20/R*(1/(Ts+273.15) - 1/(Tref+273.15)));
mrem=mrem*exp(-20/R*(1/(Ts+273.15) - 1/(Tref+273.15)));
%%%%
gd = gd*fClay2; gd(gd>1)=1; 
%%%% Control of available space for accumulating MOC 
fd = ((1-fd)*((B(11)*1000)./(4.825*(100-Psan*100).^0.6287)).^3.322+fd); fd(fd>1)=1; 
%%%%%% Soil moisture accessibility of DOC
%%% Unresolved issue if to use or not
fsm= 1;% (Se./Se_fc); fsm(fsm>1)=1; %%[-]
%%%%
%rem
%pemb pemf
%rep
%pepb pepf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Scaling of enzyme production with Ninorg or Norg/Ninorg 
%ROI=(B(:,26)+B(:,27)+B(:,28)+B(:,29)+B(:,30)+B(:,33))./(B(:,31)+B(:,32)); %%[-] 
knp=1;  %knp =  1-exp(-ROI/500); 
%%%% Scaling enzyme production with harsh conditions ;; Fundamental!!
%%%%%
%Cpep = 5.1285e-004*(B(18))^(-0.78); Cpep(Cpep>10)=10;
%Cpep = 4.2433e-004*(B(18))^(-0.85); Cpep(Cpep>10)=10;
%Cpep = 2.736e-004*(B(18))^(-0.85); Cpep(Cpep>6)=6;
%Cpep = 5.3409e-004*(B(18))^(-0.75); Cpep(Cpep>6)=6;
Cpep = 0.0078*(B(18))^(-0.5); Cpep(Cpep>6)=6;
%Cpep = 3.5e-004*(B(18))^(-0.85); Cpep(Cpep>5)=5;
pepb = knp*Cpep*pepb; %
pemb = knp*Cpep*pemb;%
%Cpep=0.0070*(B(19))^(-0.6); Cpep(Cpep>10)=10;
%Cpep = 9.0331e-004*(B(19))^(-0.85); Cpep(Cpep>10)=10;
%Cpep = 6.419e-004*(B(19))^(-0.85); Cpep(Cpep>3.5)=3.5;
%Cpep = 0.0016*(B(19))^(-0.75); Cpep(Cpep>3.5)=3.5;
Cpep = 0.0134*(B(19))^(-0.5); Cpep(Cpep>3.0)=3.0;
%Cpep = 8e-004*(B(19))^(-0.85); Cpep(Cpep>3)=3;
pepf = knp*Cpep*pepf; %
pemf = knp*Cpep*pemf;%
%%%% enzyme production from EM fungi 
pepem = knp*pepem; %
pemem = knp*pemem; %
%%%%%
%%%%%%
%pepb = rep*B(14)./(mrb*B(18)); pepb(pepb<0.0031)=0.0031; pepb(pepb>0.031)=0.031;
%pemb = rem*B(16)./(mrb*B(18)); pemb(pemb<0.0031)=0.0031; pemb(pemb>0.031)=0.031;
%pepf = rep*B(15)./(mrf*B(19)); pepf(pepf<0.0031)=0.0031; pepf(pepf>0.031)=0.031;
%pemf = rem*B(17)./(mrf*B(19)); pemf(pemf<0.0031)=0.0031; pemf(pemf>0.031)=0.031;
%%%%%%%% Litter Decoomposition on a 40�C basis
Tref4=BiogeoPar.T4;
fT2 = exp(3.36*(Ts-Tref4)/(Ts+31.79)); %% Kirschbaum 1995 2000 2002
fT1 = exp(3.36*(Ta-Tref4)/(Ta+31.79)); %% Kirschbaum 1995 2000 2002
fT2(Ts <-31)=0; fT1(Ta <-31)=0;
%fT = 2^((Ts-30)/10); %%% [-] %% Foley et al 1995 ;   Krinner et al., 2005

%%%%%%%%%%%%% Concentration in the soil
%%%%%%% C/N Ratios
%%% Parton et all 1988 - Krinner et al 2005  Orwin et al 2011
r_met_str_CN = 5; %% Ratio between metabolic and structural [C/N]
f_str_sur = safe_fraction(B(2)+B(3), ...
    r_met_str_CN*B(1)+B(2)+B(3));
f_str_ssr = safe_fraction(B(7)+B(8), ...
    r_met_str_CN*B(6)+B(7)+B(8));
%%%%
CN_met_sur_lit= safe_ratio(B(1),(1-f_str_sur)*B_physical(23)); %% [gC/gN]
CN_str_sur_lit= safe_ratio(B(2)+B(3),f_str_sur*B_physical(23)); %% [gC/gN]
CN_wod_lit= safe_ratio(B(4)+B(5),B_physical(24)); %% [gC/gN]
CN_met_ssr_lit= safe_ratio(B(6),(1-f_str_ssr)*B_physical(25)); %% [gC/gN]
CN_str_ssr_lit= safe_ratio(B(7)+B(8),f_str_ssr*B_physical(25)); %% [gC/gN]
CN_som = safe_ratio(B(9)+B(10)+B(11),B_physical(26)); %% [gC/gN]
CN_bac = safe_ratio(B(18),B_physical(27)); %% [gC/gN]
CN_fun = safe_ratio(B(19),B_physical(28)); %% [gC/gN]
%CN_ew = B(22)/B(34); %% [gC/gN]
CN_mycAM = safe_ratio(B(20),B_physical(29)); %% [gC/gN]
CN_mycEM = safe_ratio(B(21),B_physical(30)); %% [gC/gN]
%%%  C/P Ratios
CP_bac = safe_ratio(B(18),B_physical(39)); %% [gC/gP]
CP_fun = safe_ratio(B(19),B_physical(40)); %% [gC/gP]
CP_mycAM = safe_ratio(B(20),B_physical(41)); %% [gC/gP]
CP_mycEM = safe_ratio(B(21),B_physical(42)); %% [gC/gP]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% CUE-Litter and CUE-Microbe/Fungi
%%% Sinsabaguh et al 2013 ; Manzoni et al 2008
%rLO=1./L_cn; rb=1/B_cn;
%E= 0.45*rLO.^0.76./rb;
CUEmax= BiogeoPar.CUEmax;
%E = CUEmax./(1 +0.015*L_cn);
%TER_cn = 2.33*CN.^0.78;
%TER_cp = 2.91*CP.^0.83;
Ae=BiogeoPar.Ae;
%E = min(CUEmax,min(Ae*B_cb./TER_cn,Ae*B_cp./TER_cp));
%%%%%%%%%%%%%%%%%% SC_par(4)
CUE_met_sur=CUEmax./(1 +0.015*CN_met_sur_lit);%1-0.55;
CUE_str_sur=CUEmax./(1 +0.015*CN_str_sur_lit);%1-0.45;
CUE_wod_sur=CUEmax./(1 +0.015*CN_wod_lit);%1-0.55;
CUE_met_ssr=CUEmax./(1 +0.015*CN_met_ssr_lit);%1-0.55;
CUE_str_ssr=CUEmax./(1 +0.015*CN_str_ssr_lit);%1-0.55;
%%%%%
%%% Microbe
Ecb = BiogeoPar.Ecb; % carbon use efficiency [-]
%Ecb = min(CUEmax,min(Ae*CN_bac./TER_cn,Ae*CP_bac./TER_cp));
%%% Fungi
Ecf =BiogeoPar.Ecf;
%Ecf =  min(CUEmax,min(Ae*CN_fun./TER_cn,Ae*CP_fun./TER_cp));
%%%%%%%%%%%%%%%%%%%%
%%% Temperature effects on CUE
Tref=BiogeoPar.T2;
mtem = BiogeoPar.mtem; %%[�C] [0 - 0.016] Allison et al 2010; Li et al 2014
rr_met_sur = 1 - min(CUEmax,(CUE_met_sur + mtem*(Ts-Tref)));
rr_str_sur= 1 - min(CUEmax,(CUE_str_sur + mtem*(Ts-Tref)));
rr_wod_sur=1 - min(CUEmax,(CUE_wod_sur + mtem*(Ts-Tref)));
rr_met_ssr=1 - min(CUEmax,(CUE_met_ssr + mtem*(Ts-Tref)));
rr_str_ssr=1 - min(CUEmax,(CUE_str_ssr + mtem*(Ts-Tref)));
Ecb = min(CUEmax,Ecb + mtem*(Ts-Tref));
Ecf = min(CUEmax,Ecf + mtem*(Ts-Tref));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% opt_cons_CUE=0;
if opt_cons_CUE == 1
    Ecf = BiogeoPar.Ecf;
    Ecb = BiogeoPar.Ecb;
    rr_met_sur = BiogeoPar.rr_met_sur;
    rr_str_sur=   BiogeoPar.rr_str_sur;
    rr_wod_sur=  BiogeoPar.rr_wod_sur;
    rr_met_ssr= BiogeoPar.rr_met_ssr;
    rr_str_ssr=  BiogeoPar.rr_str_ssr;
end
%%%%%%%% Decomposition rates  %%[1/day]
D_met_sur_lit = k_met_sur*fT1;
D_str_sur_lit = k_str_sur*fT1*exp(-5*safe_fraction(B(3),B(2)+B(3)));
D_met_ssr_lit = k_met_ssr*fT2*fSM_litter;
D_str_ssr_lit = k_str_ssr*fT2*fSM_litter* ...
    exp(-5*safe_fraction(B(8),B(7)+B(8)));
D_wod_sur_lit = k_wod_sur*fT1*exp(-5*safe_fraction(B(5),B(4)+B(5)));

%%% Decomposition Fluxes  [gC/gsoil  day]
D1=D_met_sur_lit*B(1);
D2=D_str_sur_lit*B(2);
D3=D_str_sur_lit*B(3);
D4=D_wod_sur_lit*B(4);
D5=D_wod_sur_lit*B(5);
D6=D_met_ssr_lit*B(6);
D7=D_str_ssr_lit*B(7);
D8=D_str_ssr_lit*B(8);
D1=limit_outgoing_fluxes(B_physical(1),0,D1,dtd);
D2=limit_outgoing_fluxes(B_physical(2),0,D2,dtd);
D3=limit_outgoing_fluxes(B_physical(3),0,D3,dtd);
D4=limit_outgoing_fluxes(B_physical(4),0,D4,dtd);
D5=limit_outgoing_fluxes(B_physical(5),0,D5,dtd);
D6=limit_outgoing_fluxes(B_physical(6),0,D6,dtd);
D7=limit_outgoing_fluxes(B_physical(7),0,D7,dtd);
D8=limit_outgoing_fluxes(B_physical(8),0,D8,dtd);
%%% Input to SOC
Ip_som_poc_lig = (1-rr_str_sur)*D3 + (1-rr_wod_sur)*D5 +  (1-rr_str_ssr)*D8 ; %%% [gC/gsoil d]
Ip_som_poc_cel = (1-rr_met_sur)*D1 + (1-rr_str_sur)*D2 +  (1-rr_met_ssr)*D6  + (1-rr_str_ssr)*D7  + (1-rr_wod_sur)*D4 ; %% [gC/gsoil d]
%%% leaching-C
Id = lambda_C*Ip_som_poc_cel;  %% [gC/gsoil d]
Ip_som_poc_cel = (1-lambda_C)*Ip_som_poc_cel;
%%%%
R_litter=  (rr_str_sur)*D3 + (rr_wod_sur)*D5 +  (rr_str_ssr)*D8 + (rr_met_sur)*D1 + (rr_str_sur)*D2 +  (rr_met_ssr)*D6  + (rr_str_ssr)*D7  + (rr_wod_sur)*D4; % [gC/gsoil d]
R_litter_sur =   (rr_str_sur)*D3 + (rr_wod_sur)*D5  + (rr_wod_sur)*D4 + (rr_met_sur)*D1 + (rr_str_sur)*D2 ; % [gC/gsoil d]
%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Soil Carbon Fluxes  %% [gC/gsoil d]
F1b = (1/Ecb)*(Vdb+mrb)*(fsm*B(12))*B(18)/(Kdb+(fsm*B(12))); %% Flux from DOC to Bacteria
F1f = (1/Ecf)*(Vdf+mrf)*(fsm*B(13))*B(19)/(Kdf+(fsm*B(13))); %% Flux from DOC to Fungi
F2c_b =  Vpc*B(14)*B(10)/(Kpc + B(10)); %% Flux from POC-Cellulose to DOC - Bacteria
F2c_f =  Vpc*B(15)*B(10)/(Kpc + B(10)); %% Flux from POC-Cellulose to DOC - Fungi
F2l =  Vpl*B(15)*B(9)/(Kpl + B(9));  %% Flux from POC-Lignin to DOC  (Fungi only)
F3_b =   Vm*B(16)*B(11)/(Km + B(11)); %%% Flux from MOC to DOC-Bacteria
F3_f =   Vm*B(17)*B(11)/(Km + B(11)); %%% Flux from MOC to DOC-Fungi
F4b = (1/Ecb - 1)*Vdb*B(18)*(fsm*B(12))/(Kdb+(fsm*B(12)));  %% Respiration growth Bacteria
F4f = (1/Ecf - 1)*Vdf*B(19)*(fsm*B(13))/(Kdf+(fsm*B(13)));  %% Respiration growth  Fungi
F5b = (1/Ecb - 1)*mrb*B(18)*(fsm*B(12))/(Kdb+(fsm*B(12)));  %% Respiration Maintenance Bacteria
F5f = (1/Ecf - 1)*mrf*B(19)*(fsm*B(13))/(Kdf+(fsm*B(13)));  %% Respiration Maintenance  Fungi
F8b = (1-pepb -pemb)*mrb*B(18);  %% Microbe Bacteria  mortality
F8f = (1-pepf -pemf)*mrf*B(19); %% Fungi  mortality
F8am = mram*B(20);  %% AM  mortality
F8em = (1- pepem - pemem)*mrem*B(21);  %% EM  mortality
F9ep_b= pepb*mrb*B(18);  %% Enzyme Production
F9em_b= pemb*mrb*B(18);
F9ep_f =pepf*mrf*B(19);
F9em_f= pemf*mrf*B(19);
F9ep_em =pepem*mrem*B(21);
F9em_em= pemem*mrem*B(21);
F10ep_b= rep*B(14); %% Enzyme Turnover
F10ep_f= rep*B(15); %% Enzyme Turnover
F10em_b= rem*B(16);
F10em_f= rem*B(17);
%%%%%%%%%%%%%%%%%%%%%
Resp_mycAM = mram*B(20); %%% Respiration maintenance AM-Mycorrhizal
Resp_mycEM = mrem*B(21); %%% Respiration maintenance EM-Mycorrhizal
F10ep_b=limit_outgoing_fluxes(B_physical(14),0,F10ep_b,dtd);
F10ep_f=limit_outgoing_fluxes(B_physical(15),0,F10ep_f,dtd);
F10em_b=limit_outgoing_fluxes(B_physical(16),0,F10em_b,dtd);
F10em_f=limit_outgoing_fluxes(B_physical(17),0,F10em_f,dtd);
turnover_b = limit_outgoing_fluxes(B_physical(18),0, ...
    [F8b F9ep_b F9em_b],dtd);
F8b = turnover_b(1); F9ep_b = turnover_b(2); F9em_b = turnover_b(3);
turnover_f = limit_outgoing_fluxes(B_physical(19),0, ...
    [F8f F9ep_f F9em_f],dtd);
F8f = turnover_f(1); F9ep_f = turnover_f(2); F9em_f = turnover_f(3);
turnover_am = limit_outgoing_fluxes(B_physical(20),0, ...
    [F8am Resp_mycAM],dtd);
F8am = turnover_am(1); Resp_mycAM = turnover_am(2);
turnover_em = limit_outgoing_fluxes(B_physical(21),0, ...
    [F8em Resp_mycEM F9ep_em F9em_em],dtd);
F8em = turnover_em(1); Resp_mycEM = turnover_em(2);
F9ep_em = turnover_em(3); F9em_em = turnover_em(4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Earthworms
fpal = 1.6 - (0.04*CN_som); fpal(fpal<0)=0;
Tref=BiogeoPar.T3; Tref_k=Tref+273.15;
Ts_k=Ts+273.15;
fT3 =exp(50*(Ts_k-Tref_k)./(Tref_k*R*Ts_k)).*(1+exp((Tref_k*0.649 - 200)/(Tref_k*R)))./(1+exp((Ts_k*0.649-200)./(Ts_k*R))); %% Assimilation
fT4 = exp(-20/R*(1/(Ts+273.15) - 1/(Tref+273.15))); %% Respiration
fClay3 =  2.938*exp(-4.82*Pcla) + -2.509*exp(-11.82*Pcla); fClay3(Pcla<0.15)=1;
%Few  = fabsE*B(22)*B(10)/(kse + B(10))*VmaxE*fpal*fT3*Aew; %%% [gC/ gsoil d] POC consumed by earthworms
Few  = fabsE*B(10)*VmaxE*fpal*fT3*fClay3*fPH*Aew; %%% % [gC/ gsoil d] POC consumed by earthworms
RespEw_m =  Aew*rmanEa*B(22)*fT4 + (1-Aew)*rmanEp*B(22)*fT4 ; %%% Maintenance Respiration
M_ew = dew*B(22)*fT4; %% % [gC/gsoil d]  Death of Earthworms
earthworm_turnover = limit_outgoing_fluxes(B_physical(22),0, ...
    [RespEw_m M_ew],dtd);
RespEw_m = earthworm_turnover(1); M_ew = earthworm_turnover(2);
DOC_mortality = F8b+F8f+F8am+F8em+M_ew;
F2l = limit_outgoing_fluxes(B_physical(9),Ip_som_poc_lig,F2l,dtd);
poc_flux = limit_outgoing_fluxes(B_physical(10), ...
    Ip_som_poc_cel+(1-gd)*DOC_mortality,[F2c_b F2c_f Few],dtd);
F2c_b = poc_flux(1); F2c_f = poc_flux(2); Few = poc_flux(3);
moc_source = (1-fd)*(F2l+F2c_b+F2c_f);
moc_flux = limit_outgoing_fluxes(B_physical(11),moc_source, ...
    [F3_b F3_f],dtd);
F3_b = moc_flux(1); F3_f = moc_flux(2);
RespEw_g = (1-CUEmax)*Few ; %% Growth Respiration
%CUEapp  = CUEmax*((Few + RespEw_m)./Few);
RespE = RespEw_g + RespEw_m ; %% [gC/gsoil d]  Total EW Respiration
FewN =  (Few - RespE - M_ew)/R_CN_ew; %% [gN/gsoil d]  N-requirement
FewN_loss = limit_outgoing_fluxes(B_physical(34),0,max(0,-FewN),dtd);
FewN = max(0,FewN)-FewN_loss;
DOC_total = B_physical(12)+B_physical(13);
if DOC_total > 0
    Rbf = B_physical(12)/DOC_total; %% Ratio of DOC - Bacteria to B+F
else
    Rbf = 0.5;
end
DOC_source_b = Rbf*Id + fd*F2c_b + Rbf*gd*DOC_mortality + ...
    F3_b + F10ep_b + F10em_b;
DOC_source_f = (1-Rbf)*Id + fd*(F2l+F2c_f) + ...
    (1-Rbf)*gd*DOC_mortality + F3_f + F10ep_f + F10em_f;
DOC_avail_b = max(0,B_physical(12)/dtd+DOC_source_b-Rbf*LEAK_DOC);
DOC_avail_f = max(0,B_physical(13)/dtd+DOC_source_f-(1-Rbf)*LEAK_DOC);
if F1b > DOC_avail_b
    DOC_scale_b = DOC_avail_b/F1b;
    F1b = F1b*DOC_scale_b;
    F4b = F4b*DOC_scale_b;
    F5b = F5b*DOC_scale_b;
end
if F1f > DOC_avail_f
    DOC_scale_f = DOC_avail_f/F1f;
    F1f = F1f*DOC_scale_f;
    F4f = F4f*DOC_scale_f;
    F5f = F5f*DOC_scale_f;
end
R_microbe =  F4b+F4f+F5f+F5b + Resp_mycAM + Resp_mycEM; % [gC/gsoil d]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Nitrogen Fluxes  %%% Kirschbaum et al 2002
N_mort_b = proportional_flux(F8b,B_physical(27),B_physical(18));
N_mort_f = proportional_flux(F8f,B_physical(28),B_physical(19));
N_mort_am = proportional_flux(F8am,B_physical(29),B_physical(20));
N_mort_em = proportional_flux(F8em,B_physical(30),B_physical(21));
if CN_bac < R_CN_bac
    %%% Minerelazation
    Nmin_b = B_physical(27) - B(18)/R_CN_bac;  %% [gN/gsoil]
    Nimm_b=0;
else
    %%% Immobilization
    Nimm_b = B(18)/R_CN_bac - B_physical(27);  %% [gN/gsoil]
    %if Nimm_b > 0.90*B(31)
    %    Nimm_b = 0.90*B(31);
    %end
    Nmin_b=0;
end
if CN_fun < R_CN_fun
    %%% Minerelazation
    Nmin_f = B_physical(28) - B(19)/R_CN_fun;  %% [gN/gsoil]
    Nimm_f=0;
else
    %%% Immobilization
    Nimm_f = B(19)/R_CN_fun - B_physical(28);  %% [gN/gsoil]
    Nmin_f=0;
    %if Nimm_f > 0.90*B(31)
    %    Nimm_f = 0.90*B(31);
    %end
end
if ExEM < 1
    if CN_mycAM < R_CN_mycAM
        %%% Minerelazation
        Nmin_am = B_physical(29) - B(20)/R_CN_mycAM;  %% [gN/gsoil]
        Nimm_am=0;
    else
        %%% Immobilization
        Nimm_am = B(20)/R_CN_mycAM - B_physical(29);  %% [gN/gsoil]
        Nmin_am=0;
        %if Nimm_am > 0.90*B(31)
        %    Nimm_am = 0.90*B(31);
        %end
    end
else
    Nimm_am=0; Nmin_am=0;
end
if ExEM > 0
    if CN_mycEM < R_CN_mycEM
        %%% Minerelazation
        Nmin_em = B_physical(30) - B(21)/R_CN_mycEM;  %% [gN/gsoil]
        Nimm_em=0;
    else
        %%% Immobilization
        Nimm_em = B(21)/R_CN_mycEM - B_physical(30);  %% [gN/gsoil]
        Nmin_em=0;
        %if Nimm_em > 0.90*B(31)
        %    Nimm_em = 0.90*B(31);
        %end
    end
else
    Nimm_em=0; Nmin_em=0;
end
Nmin_b = dtd*limit_outgoing_fluxes(B_physical(27), ...
    -N_mort_b,Nmin_b/dtd,dtd);
Nmin_f = dtd*limit_outgoing_fluxes(B_physical(28), ...
    -N_mort_f,Nmin_f/dtd,dtd);
Nmin_am = dtd*limit_outgoing_fluxes(B_physical(29), ...
    -N_mort_am,Nmin_am/dtd,dtd);
Nmin_em = dtd*limit_outgoing_fluxes(B_physical(30), ...
    -N_mort_em,Nmin_em/dtd,dtd);
%%% Input to SOC - Nitrogen
N_sur_flux = limit_outgoing_fluxes(B_physical(23),0, ...
    [proportional_flux(D1,(1-f_str_sur)*B_physical(23),B_physical(1)) ...
    proportional_flux(D2+D3,f_str_sur*B_physical(23), ...
    B_physical(2)+B_physical(3))],dtd);
N_wod_flux = limit_outgoing_fluxes(B_physical(24),0, ...
    proportional_flux(D4+D5,B_physical(24), ...
    B_physical(4)+B_physical(5)),dtd);
N_ssr_flux = limit_outgoing_fluxes(B_physical(25),0, ...
    [proportional_flux(D6,(1-f_str_ssr)*B_physical(25),B_physical(6)) ...
    proportional_flux(D7+D8,f_str_ssr*B_physical(25), ...
    B_physical(7)+B_physical(8))],dtd);
Ip_som_nit = sum(N_sur_flux)+N_wod_flux+sum(N_ssr_flux);
Ip_lea_nit = (lambda_N)*Ip_som_nit;  %% leaching of mineral-N
Ip_som_nit =  (1-lambda_N)*Ip_som_nit;
som_C = B_physical(9)+B_physical(10)+B_physical(11);
N_decomp_b = proportional_flux(fd*F2c_b+F3_b,B_physical(26),som_C);
N_decomp_f = proportional_flux(fd*(F2l+F2c_f)+F3_f, ...
    B_physical(26),som_C);
N_mortality = N_mort_b+N_mort_f+N_mort_am+N_mort_em;
N_som_flux = limit_outgoing_fluxes(B_physical(26), ...
    Ip_som_nit+N_mortality+max(0,-FewN), ...
    [N_decomp_b N_decomp_f max(0,FewN)],dtd);
N_decomp_b = N_som_flux(1); N_decomp_f = N_som_flux(2);
FewN = N_som_flux(3)+min(0,FewN);
%%%%%%
N_imm_fix_b = 0.011*(B(31))/dtd; %% [gN/gsoil day] %% Continous immobilization due to chemical stabilization
N_imm_fix_f = 0.011*(B(32))/dtd; %% [gN/gsoil day] %% Continous immobilization due to chemical stabilization
%%%%
NMIN = (Nmin_b + Nmin_f + Nmin_am + Nmin_em)/dtd; %% [gN/gsoil day]
NIMM = (Nimm_b + Nimm_f + Nimm_am + Nimm_em  + N_imm_fix_f*dtd + N_imm_fix_b*dtd)/dtd; %% [gN/gsoil day]
NMIN_IN = NMIN - NIMM; %% [gN/gsoil day]
%%%% External source, Biological fixation, deposition, fertilization
%%%  added as NH4 only to simplifiy -->  Dickinson et al., 2002
Incoming_N = (FertN + DepN + BfixN)/(ZBIOG*rsd);  %% [gN/ gsoil d]
%%%% Volatilization of Ammonia NH3
kv = BiogeoPar.kv; %% [1/d]
VOL = kv*B(31);  %% [gN/gsoil d]  Dickinson et al., 2002
%%%%%%% Nitrification (Aerobic)
KNmax = BiogeoPar.KNmaxn; %% [1/d] % Dickinson et al., 2002
ftemp=power((70-Ts)/(70-38),12)*exp(12*(Ts-38)/(70-38)); %%% Xu-Ri and Prentice 2008
fmoist = Se*(1-Se)/(0.25);  %% Dickinson et al., 2002
%%% pH response of nitrification rate - Parton et al 1996, Global Biogeochem Cycles 10:401
fpH_nit = 0.56 + atan(pi*0.45*(-5+PH))/pi;
fpH_nit(fpH_nit>1)=1; fpH_nit(fpH_nit<0)=0;
NO3flx = KNmax*ftemp*fmoist*fpH_nit*B(31); %% [gN/gsoil d]
NH4_reaction_avail = max(0,B(31)/dtd + Incoming_N + ...
    (1-f_org_lea)*Ip_lea_nit - LEAK_NH4);
NH4_reaction = NO3flx+VOL;
if NH4_reaction > NH4_reaction_avail
    reaction_scale = NH4_reaction_avail/NH4_reaction;
    NO3flx = NO3flx*reaction_scale;
    VOL = VOL*reaction_scale;
end
%%% N2O yield of nitrification - Goodroad and Keeney 1984; Khalil et al 2004
N2Oflx_nit = BiogeoPar.f_N2O_nit*NO3flx; %% [gN/gsoil d]
%%%%%%% Denitrification (Anaerobic)
Kc = BiogeoPar.Kc/(ZBIOG*rsd); %% [gC/ gsoil] %% Li et al., 1992
%Kn = BiogeoPar.Kn/(ZBIOG*rsd); %% [gN/ gsoil] % Li et al., 1992
ftemp=exp(308.56*(1/68.02-1/(Ts+46.02))); % Xu-Ri and Prentice 2008
%%% Anaerobic fraction. Under Brooks-Corey Se = (Pe/|Psi|)^L, so Se^(1/L)
%%% is Pe/|Psi|: the air-entry suction over the current suction, which is 1
%%% at saturation and falls as the profile drains. L is the pore size
%%% distribution index of the same retention curve that produces Se, so the
%%% two stay consistent. The literal 0.5 that stood here was a placeholder
%%% for L, which was never passed in; at this site L = 0.150, so it applied
%%% an exponent of 2 where the soil calls for 6.7. Se^2 returns 0.53 at
%%% field capacity, running denitrification at half rate in an aerobic
%%% soil on 99% of days; Se^(1/L) returns 0.12.
fmoist=Se^(1/L_pore); %% % Dickinson et al., 2002
%N2flx = (B(13)/(B(13) + Kc))*ftemp*fmoist*(B(28)/(B(28) + Kn));% [gN/gsoil d]
KNmax= BiogeoPar.KNmaxd; %% [1/d] % Dickinson et al., 2002
%%% pH response of denitrification rate - Heinen 2006, Geoderma 133:444 (ANIMO/SWAT-family)
fpH_den = (PH-3.5)/3.0;
fpH_den(fpH_den>1)=1; fpH_den(fpH_den<0)=0;
%%% Denitrifiers need NO3 as electron acceptor AND labile C as donor, so
%%% the rate carries a Michaelis-Menten term on dissolved organic carbon
%%% (Li et al. 1992). Without it denitrification and leaching are both
%%% strictly first order in B(32), their ratio is independent of pool
%%% size, and the split between them can only be set by rate constants.
%%% At this site DOC averages 91 gC/m2 against Kc = 17.5, so the term
%%% sits near saturation (0.84) and modulates weakly; it matters most
%%% when labile C is drawn down.
f_DOC_den = safe_fraction(B(12)+B(13),(B(12)+B(13))+Kc); %% [-]
DENITflx = KNmax*ftemp*fmoist*fpH_den*f_DOC_den*B(32); %% [gN/gsoil d] total denitrification (N2+N2O)
NO3_reaction_avail = max(0,B(32)/dtd + NO3flx - ...
    N2Oflx_nit - LEAK_NO3);
DENITflx = min(DENITflx,NO3_reaction_avail);
%%% N2O/(N2O+N2) product ratio: N2O reductase inhibited at low pH
%%% Simek and Cooper 2002; Cuhel et al 2010; Liu et al 2010
rN2O = 1/(1+exp(BiogeoPar.k_rN2O*(PH-BiogeoPar.pH0_rN2O))); %% [-]
N2flx = (1-rN2O)*DENITflx; %% [gN/gsoil d] N2 emission
N2Oflx = N2Oflx_nit + rN2O*DENITflx; %% [gN/gsoil d] total N2O emission
%%%%%
%%%%%%% Porporato et al., 2003
imm_cut = zeros(1,4); %%% [gN/gsoil d] unmet demand withheld from microbes
DON_source = f_org_lea*Ip_lea_nit + lambda_N* ...
    (N_decomp_b+N_decomp_f);
S_DON = min(0.01*B_physical(33), ...
    max(0,B_physical(33)/dtd+DON_source-LEAK_DON));
som_N_preimm = Ip_som_nit-(N_decomp_b+N_decomp_f) + ...
    N_mortality-FewN+S_DON;
if NMIN_IN >= 0
    NH4_in = 0.90*NMIN_IN; %%% [gN/gsoil d] Mineralization NH4
    NO3_in = 0.10*NMIN_IN; %%% [gN/gsoil d] Mineralization NO3
    Norg_imm = 0;
else
    %%% Immobilization is limited by mineral N remaining after other sinks
    NH4_avail = max(0,B(31) + dtd*(Incoming_N + ...
        (1-f_org_lea)*Ip_lea_nit - NO3flx - LEAK_NH4 - VOL));
    NO3_avail = max(0,B(32) + dtd*(NO3flx - N2Oflx_nit - ...
        LEAK_NO3 - DENITflx));
    mineral_avail = NH4_avail + NO3_avail;
    mineral_imm = min(-NMIN_IN*dtd,mineral_avail);
    if mineral_avail > 0
        if mineral_imm == mineral_avail
            NH4_imm = NH4_avail;
            NO3_imm = NO3_avail;
        else
            NH4_imm = mineral_imm*NH4_avail/mineral_avail;
            NO3_imm = mineral_imm-NH4_imm;
        end
    else
        NH4_imm = 0;
        NO3_imm = 0;
    end
    NH4_in = -NH4_imm/dtd;
    NO3_in = -NO3_imm/dtd;
    Norg_dem = (-NMIN_IN*dtd-mineral_imm)/dtd;
    Norg_avail = max(0,B_physical(26)/dtd+som_N_preimm);
    Norg_imm = min(Norg_dem,Norg_avail);
    N_imm_deficit = Norg_dem - Norg_imm;
    if N_imm_deficit > 0
        imm_gross = [Nimm_b + N_imm_fix_b*dtd, Nimm_f + N_imm_fix_f*dtd, ...
            Nimm_am, Nimm_em]/dtd;
        imm_gross_tot = sum(imm_gross);
        if N_imm_deficit >= imm_gross_tot
            imm_cut = imm_gross;
        elseif imm_gross_tot > 0
            imm_cut = N_imm_deficit*imm_gross/imm_gross_tot;
        end
    end
end
%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
%Phosphorus Fluxes
P_mort_b = proportional_flux(F8b,B_physical(39),B_physical(18));
P_mort_f = proportional_flux(F8f,B_physical(40),B_physical(19));
P_mort_am = proportional_flux(F8am,B_physical(41),B_physical(20));
P_mort_em = proportional_flux(F8em,B_physical(42),B_physical(21));
if CP_bac < R_CP_bac
    %%% Minerelazation
    Pmin_b = B_physical(39) - B(18)/R_CP_bac;  %% [gP/gsoil]
    Pimm_b=0;
else
    %%% Immobilization
    Pimm_b = B(18)/R_CP_bac - B_physical(39);  %% [gP/gsoil]
    Pmin_b=0;
    %if Pimm_b > B(43)
    %    Pimm_b = B(43);
    %end
end
if CP_fun < R_CP_fun
    %%% Minerelazation
    Pmin_f = B_physical(40) - B(19)/R_CP_fun;  %% [gP/gsoil]
    Pimm_f=0;
else
    %%% Immobilization
    Pimm_f = B(19)/R_CP_fun - B_physical(40);  %% [g/gsoil]
    Pmin_f =0;
    %if Pimm_f > B(43)
    %    Pimm_f = B(43);
    %end
end
if ExEM < 1
    if CP_mycAM < R_CP_mycAM
        %%% Minerelazation
        Pmin_am = B_physical(41) - B(20)/R_CP_mycAM;  %% [gP/gsoil]
        Pimm_am=0;
    else
        %%% Immobilization
        Pimm_am = B(20)/R_CP_mycAM - B_physical(41);  %% [gP/gsoil]
        Pmin_am =0;
        %if Pimm_am > B(43)
        %    Pimm_am = B(43);
        %end
    end
else
    Pimm_am=0;Pmin_am =0;
end
if ExEM > 0
    if CP_mycEM < R_CP_mycEM
        %%% Minerelazation
        Pmin_em = B_physical(42) - B(21)/R_CP_mycEM;  %% [gP/gsoil]
        Pimm_em=0;
    else
        %%% Immobilization
        Pimm_em = B(21)/R_CP_mycEM - B_physical(42);  %% [g/gsoil]
        Pmin_em =0;
        %if Pimm_em > B(43)
        %    Pimm_em = B(43);
        %end
    end
else
    Pimm_em=0;Pmin_em =0;
end
Pmin_b = dtd*limit_outgoing_fluxes(B_physical(39), ...
    -P_mort_b,Pmin_b/dtd,dtd);
Pmin_f = dtd*limit_outgoing_fluxes(B_physical(40), ...
    -P_mort_f,Pmin_f/dtd,dtd);
Pmin_am = dtd*limit_outgoing_fluxes(B_physical(41), ...
    -P_mort_am,Pmin_am/dtd,dtd);
Pmin_em = dtd*limit_outgoing_fluxes(B_physical(42), ...
    -P_mort_em,Pmin_em/dtd,dtd);
%%% Input to SOC - Phosporus
P_sur_flux = limit_outgoing_fluxes(B_physical(35),0, ...
    [proportional_flux(D1,(1-f_str_sur)*B_physical(35),B_physical(1)) ...
    proportional_flux(D2+D3,f_str_sur*B_physical(35), ...
    B_physical(2)+B_physical(3))],dtd);
P_wod_flux = limit_outgoing_fluxes(B_physical(36),0, ...
    proportional_flux(D4+D5,B_physical(36), ...
    B_physical(4)+B_physical(5)),dtd);
P_ssr_flux = limit_outgoing_fluxes(B_physical(37),0, ...
    [proportional_flux(D6,(1-f_str_ssr)*B_physical(37),B_physical(6)) ...
    proportional_flux(D7+D8,f_str_ssr*B_physical(37), ...
    B_physical(7)+B_physical(8))],dtd);
Ip_som_pho = sum(P_sur_flux)+P_wod_flux+sum(P_ssr_flux);
Ip_lea_pho = (lambda_P)*Ip_som_pho;
Ip_som_pho =  (1-lambda_P)*Ip_som_pho;
P_decomp_b = proportional_flux(fd*F2c_b+F3_b,B_physical(38),som_C);
P_decomp_f = proportional_flux(fd*(F2l+F2c_f)+F3_f, ...
    B_physical(38),som_C);
P_mortality = P_mort_b+P_mort_f+P_mort_am+P_mort_em;
P_som_flux = limit_outgoing_fluxes(B_physical(38), ...
    Ip_som_pho+P_mortality,[P_decomp_b P_decomp_f],dtd);
P_decomp_b = P_som_flux(1); P_decomp_f = P_som_flux(2);
%%%%%%
K1p= BiogeoPar.K1p; %%% [1/day]
K2p= BiogeoPar.K2p; %% [1/day]
K3p= BiogeoPar.K3p; %% [1/day]
K4p= BiogeoPar.K4p; %% [1/day]
PMIN = (Pmin_b + Pmin_f + Pmin_em + Pmin_am)/dtd ; PIMM =  (Pimm_b + Pimm_f + Pimm_em + Pimm_am)/dtd;  %% [gP/gsoil d]
PMIN_IN = (PMIN - PIMM );  %% [gP/gsoil d]
Incoming_P = (FertP + DepP)/(ZBIOG*rsd);  %% [gP/ gsoil d]
Psec_exch = K1p*B_physical(43) - K2p*B_physical(45);
Ppri_exch = K4p*B_physical(44); %% %% [gP/gsoil day]
Pocc_exch = K3p*B_physical(45);
P_external_base = Incoming_P + (1-f_org_lea)*Ip_lea_pho - LEAK_P;
P_flux = limit_internal_fluxes(B_physical([43 44 45 46]), ...
    [P_external_base+max(PMIN_IN,0) Tup_P 0 0],[1 2 3],[3 1 4], ...
    [Psec_exch Ppri_exch Pocc_exch],dtd);
Psec_exch = P_flux(1); Ppri_exch = P_flux(2); Pocc_exch = P_flux(3);
P_internal_solution = -Psec_exch+Ppri_exch;
P_decomp = P_decomp_b+P_decomp_f;
DOP_source = f_org_lea*Ip_lea_pho + lambda_P*P_decomp;
S_DOP = min(0.01*B_physical(47), ...
    max(0,B_physical(47)/dtd+DOP_source-LEAK_DOP));
P_som_preimm = Ip_som_pho-P_decomp+P_mortality+S_DOP;
Porg_imm = 0;
P_imm_deficit = 0;
P_imm_cut = zeros(1,4);
if PMIN_IN < 0
    Pmineral_avail = max(0,B_physical(43)/dtd + P_external_base + ...
        P_internal_solution);
    Pmineral_imm = min(-PMIN_IN,Pmineral_avail);
    Porg_dem = -PMIN_IN-Pmineral_imm;
    Porg_avail = max(0,B_physical(38)/dtd+P_som_preimm);
    Porg_imm = min(Porg_dem,Porg_avail);
    P_imm_deficit = Porg_dem-Porg_imm;
    if P_imm_deficit > 0 && PIMM > 0
        P_imm_gross = [Pimm_b Pimm_f Pimm_am Pimm_em]/dtd;
        if P_imm_deficit >= PIMM
            P_imm_cut = P_imm_gross;
        else
            P_imm_cut = P_imm_deficit*P_imm_gross/PIMM;
        end
    end
end
P_mineral_transfer = PMIN_IN+Porg_imm+P_imm_deficit;
%PIMM =  Pimm_b + Pimm_f + Pimm_em + Pimm_am;
%PMIN_IN = (PMIN - PIMM);  %% [gP/gsoil d]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%% Input to SOC - Potassium
K_sur_flux = limit_outgoing_fluxes(B_physical(48),0, ...
    [proportional_flux(D1,(1-f_str_sur)*B_physical(48),B_physical(1)) ...
    proportional_flux(D2+D3,f_str_sur*B_physical(48), ...
    B_physical(2)+B_physical(3))],dtd);
K_wod_flux = limit_outgoing_fluxes(B_physical(49),0, ...
    proportional_flux(D4+D5,B_physical(49), ...
    B_physical(4)+B_physical(5)),dtd);
K_ssr_flux = limit_outgoing_fluxes(B_physical(50),0, ...
    [proportional_flux(D6,(1-f_str_ssr)*B_physical(50),B_physical(6)) ...
    proportional_flux(D7+D8,f_str_ssr*B_physical(50), ...
    B_physical(7)+B_physical(8))],dtd);
Ip_som_pot = sum(K_sur_flux)+K_wod_flux+sum(K_ssr_flux);
Ip_lea_pot = (lambda_K)*Ip_som_pot;
Ip_som_pot =  (1-lambda_K)*Ip_som_pot;
K_decomp = proportional_flux(fd*(F2l+F2c_b+F2c_f)+F3_b+F3_f, ...
    B_physical(51),som_C);
K_decomp = limit_outgoing_fluxes(B_physical(51),Ip_som_pot, ...
    K_decomp,dtd);
%%%%%%
Kdk=BiogeoPar.Kdk;%%%% [1/day]
Kak= BiogeoPar.Kak; %% [1/day]
if opt_ERW == 1
    Kdk=0;%%%% [1/day]
    Kak= 0; %% [1/day]  
end

%%%%
K1k= BiogeoPar.K1k ; % []; %%% [1/day]
K2k= BiogeoPar.K2k; % []; %% [1/day]
K3k= BiogeoPar.K3k; %% [1/day]
Kex_sol =  Kdk*B_physical(53)- Kak*B_physical(52); %% Absrorptio desporption
Kfix_sol = K2k*B_physical(54) - K1k*B_physical(52);  %% From fixed to solution
Kmin_rel = K3k*B_physical(55) ; %% %% [gK/gsoil day] Conversion mineral


%%% Input to SOC - Calcium
Ca_sur_flux = limit_outgoing_fluxes(B_physical(56),0, ...
    [proportional_flux(D1,(1-f_str_sur)*B_physical(56),B_physical(1)) ...
    proportional_flux(D2+D3,f_str_sur*B_physical(56), ...
    B_physical(2)+B_physical(3))],dtd);
Ca_wod_flux = limit_outgoing_fluxes(B_physical(57),0, ...
    proportional_flux(D4+D5,B_physical(57), ...
    B_physical(4)+B_physical(5)),dtd);
Ca_ssr_flux = limit_outgoing_fluxes(B_physical(58),0, ...
    [proportional_flux(D6,(1-f_str_ssr)*B_physical(58),B_physical(6)) ...
    proportional_flux(D7+D8,f_str_ssr*B_physical(58), ...
    B_physical(7)+B_physical(8))],dtd);
Ip_som_cal = sum(Ca_sur_flux)+Ca_wod_flux+sum(Ca_ssr_flux);
Ip_lea_cal = (lambda_Ca)*Ip_som_cal;
Ip_som_cal =  (1-lambda_Ca)*Ip_som_cal;
Ca_decomp = proportional_flux(fd*(F2l+F2c_b+F2c_f)+F3_b+F3_f, ...
    B_physical(59),som_C);
Ca_decomp = limit_outgoing_fluxes(B_physical(59),Ip_som_cal, ...
    Ca_decomp,dtd);
%%%%%%
Kdca=BiogeoPar.Kdca;%%%% [1/day]
Kaca= BiogeoPar.Kaca; %% [1/day]

if opt_ERW == 1
    Kdca=0;%%%% [1/day]
    Kaca= 0; %% [1/day]  
end

%%%%
K1ca= BiogeoPar.K1ca ; % []; %%% [1/day]
K2ca= BiogeoPar.K2ca; % []; %% [1/day]
K3ca= BiogeoPar.K3ca; %% [1/day]
Caex_sol =  Kdca*B_physical(61)- Kaca*B_physical(60); %% Absrorptio desporption
Cafix_sol = K2ca*B_physical(62) - K1ca*B_physical(60);  %% From fixed to solution
Camin_rel = K3ca*B_physical(63) ; %% %% [gCa/gsoil day] Conversion mineral


%%% Input to SOC - Magnesium
Mg_sur_flux = limit_outgoing_fluxes(B_physical(64),0, ...
    [proportional_flux(D1,(1-f_str_sur)*B_physical(64),B_physical(1)) ...
    proportional_flux(D2+D3,f_str_sur*B_physical(64), ...
    B_physical(2)+B_physical(3))],dtd);
Mg_wod_flux = limit_outgoing_fluxes(B_physical(65),0, ...
    proportional_flux(D4+D5,B_physical(65), ...
    B_physical(4)+B_physical(5)),dtd);
Mg_ssr_flux = limit_outgoing_fluxes(B_physical(66),0, ...
    [proportional_flux(D6,(1-f_str_ssr)*B_physical(66),B_physical(6)) ...
    proportional_flux(D7+D8,f_str_ssr*B_physical(66), ...
    B_physical(7)+B_physical(8))],dtd);
Ip_som_mag = sum(Mg_sur_flux)+Mg_wod_flux+sum(Mg_ssr_flux);
Ip_lea_mag = (lambda_Mg)*Ip_som_mag;
Ip_som_mag =  (1-lambda_Mg)*Ip_som_mag;
Mg_decomp = proportional_flux(fd*(F2l+F2c_b+F2c_f)+F3_b+F3_f, ...
    B_physical(67),som_C);
Mg_decomp = limit_outgoing_fluxes(B_physical(67),Ip_som_mag, ...
    Mg_decomp,dtd);
%%%%%%
Kdmg=BiogeoPar.Kdmg;%%%% [1/day]
Kamg= BiogeoPar.Kamg; %% [1/day]
if opt_ERW == 1
    Kdmg=0;%%%% [1/day]
    Kamg= 0; %% [1/day]  
end
%%%%
K1mg= BiogeoPar.K1mg ; % []; %%% [1/day]
K2mg= BiogeoPar.K2mg; % []; %% [1/day]
K3mg= BiogeoPar.K3mg; %% [1/day]
Mgex_sol =  Kdmg*B_physical(69)- Kamg*B_physical(68); %% Absrorptio desporption
Mgfix_sol = K2mg*B_physical(70) - K1mg*B_physical(68);  %% From fixed to solution
Mgmin_rel = K3mg*B_physical(71) ; %% %% [gMg/gsoil day] Conversion mineral


%%% Input to SOC - Silicate
Si_sur_flux = limit_outgoing_fluxes(B_physical(72),0, ...
    [proportional_flux(D1,(1-f_str_sur)*B_physical(72),B_physical(1)) ...
    proportional_flux(D2+D3,f_str_sur*B_physical(72), ...
    B_physical(2)+B_physical(3))],dtd);
Si_wod_flux = limit_outgoing_fluxes(B_physical(73),0, ...
    proportional_flux(D4+D5,B_physical(73), ...
    B_physical(4)+B_physical(5)),dtd);
Si_ssr_flux = limit_outgoing_fluxes(B_physical(74),0, ...
    [proportional_flux(D6,(1-f_str_ssr)*B_physical(74),B_physical(6)) ...
    proportional_flux(D7+D8,f_str_ssr*B_physical(74), ...
    B_physical(7)+B_physical(8))],dtd);
Ip_som_sil = sum(Si_sur_flux)+Si_wod_flux+sum(Si_ssr_flux);
Ip_lea_sil = (lambda_Si)*Ip_som_sil;
Ip_som_sil =  (1-lambda_Si)*Ip_som_sil;
Si_decomp = proportional_flux(fd*(F2l+F2c_b+F2c_f)+F3_b+F3_f, ...
    B_physical(75),som_C);
Si_decomp = limit_outgoing_fluxes(B_physical(75),Ip_som_sil, ...
    Si_decomp,dtd);
%%%%%%
Kdsi=BiogeoPar.Kdsi;%%%% [1/day]
Kasi= BiogeoPar.Kasi; %% [1/day]
if opt_ERW == 1
    Kdsi=0;%%%% [1/day]
    Kasi= 0; %% [1/day]  
end
%%%%
K1si= BiogeoPar.K1si ; % []; %%% [1/day]
K2si= BiogeoPar.K2si; % []; %% [1/day]
K3si= BiogeoPar.K3si; %% [1/day]
Siex_sol =  0; %% Absrorptio desporption
Sifix_sol = K2si*B_physical(78) - K1si*B_physical(76);  %% From fixed to solution
Simin_rel = K3si*B_physical(79) ; %% %% [gSi/gsoil day] Conversion mineral

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Incoming_K = (FertK + DepK)/(ZBIOG*rsd);  %% [gK/ gsoil d]

Incoming_Ca = (FertCa + DepCa)/(ZBIOG*rsd);  %% [gCa/ gsoil d]
Incoming_Mg = (FertMg + DepMg)/(ZBIOG*rsd);  %% [gMg/ gsoil d]
Incoming_Si = (FertSi + DepSi)/(ZBIOG*rsd);  %% [gSi/ gsoil d]

K_base = Incoming_K+Ip_lea_pot+K_decomp-LEAK_K;
K_flux = limit_internal_fluxes(B_physical(52:55),[K_base 0 0 Tup_K], ...
    [2 3 4],[1 1 1],[Kex_sol Kfix_sol Kmin_rel],dtd);
Kex_sol = K_flux(1); Kfix_sol = K_flux(2); Kmin_rel = K_flux(3);

Ca_base = Incoming_Ca+Ip_lea_cal+Ca_decomp-LEAK_Ca;
Ca_flux = limit_internal_fluxes(B_physical(60:63),[Ca_base 0 0 Tup_Ca], ...
    [2 3 4],[1 1 1],[Caex_sol Cafix_sol Camin_rel],dtd);
Caex_sol = Ca_flux(1); Cafix_sol = Ca_flux(2); Camin_rel = Ca_flux(3);

Mg_base = Incoming_Mg+Ip_lea_mag+Mg_decomp-LEAK_Mg;
Mg_flux = limit_internal_fluxes(B_physical(68:71),[Mg_base 0 0 Tup_Mg], ...
    [2 3 4],[1 1 1],[Mgex_sol Mgfix_sol Mgmin_rel],dtd);
Mgex_sol = Mg_flux(1); Mgfix_sol = Mg_flux(2); Mgmin_rel = Mg_flux(3);

Si_base = Incoming_Si+Ip_lea_sil+Si_decomp-LEAK_Si;
Si_flux = limit_internal_fluxes(B_physical(76:79),[Si_base 0 0 Tup_Si], ...
    [2 3 4],[1 1 1],[Siex_sol Sifix_sol Simin_rel],dtd);
Siex_sol = Si_flux(1); Sifix_sol = Si_flux(2); Simin_rel = Si_flux(3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% MASS BALANCE %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Carbon [gC/g soil d]
dB(1)= IS.C_met_sur_lit -  D1; %%% B1  Above-ground Litter Metabolic
dB(2) = IS.C_str_sur_lit_nlig  - D2 ;  %%% B2  Above-ground Litter Structural - Cellulose/Hemicellulose
dB(3)= IS.C_str_sur_lit_lig - D3 ; %%% B3  Above-ground Litter Structura - Lignin
dB(4)= IS.C_wod_sur_lit_nlig - D4 ;  %%% B4  Above-ground Woody  - Cellulose/Hemicellulose
dB(5)= IS.C_wod_sur_lit_lig - D5 ;  %%% B5  Above-ground Woody - Lignin
dB(6)= IS.C_met_ssr_lit - D6; %%% B6 Below-ground Litter Metabolic
dB(7)= IS.C_str_ssr_lit_nlig - D7 ;  %%% B7 Below-ground Litter Structural - Cellulose/Hemicellulose
dB(8)= IS.C_str_ssr_lit_lig - D8 ; %%% B8 Below-ground Litter Structura - Lignin
%%%%
dB(9)=  Ip_som_poc_lig  -F2l ; %%% B9 SOM-POC- lignin
dB(10)= Ip_som_poc_cel - F2c_b - F2c_f  + (1-gd)*(F8b+F8f+F8am+F8em+M_ew) - Few ;  %%% B10 SOM-POC -  Cellulose/Hemicellulose
dB(11)= (1-fd)*(F2l + F2c_b + F2c_f) - F3_b - F3_f;  %%% B11 SOM-MOC
dB(12)= Rbf*Id + fd*(F2c_b) + Rbf*gd*(F8b+F8f+F8am+F8em+M_ew) + F3_b + (F10ep_b + F10em_b) - F1b - Rbf*LEAK_DOC ;  %%% B12 DOC - for bacteria
dB(13)= (1-Rbf)*Id + fd*(F2l + F2c_f) + (1-Rbf)*gd*(F8b+F8f+F8am+F8em+M_ew) + F3_f + (F10ep_f + F10em_f) - F1f  - (1-Rbf)*LEAK_DOC ;  %%% B13 DOC - Fungi
dB(14)= F9ep_b + Rbf*F9ep_em - F10ep_b;  %%% B14 Enzyme for decomposition of POC-Bact
dB(15)= F9ep_f + (1-Rbf)*F9ep_em  - F10ep_f;  %%% B15 Enzyme for decomposition of POC-Fung
dB(16)= F9em_b + Rbf*F9em_em - F10em_b; %%% B16 Enzyme for decomposition of MOC-Bact
dB(17)= F9em_f + (1-Rbf)*F9em_em  - F10em_f;  %%% B17 Enzyme for decomposition of MOC-Fung
dB(18)= F1b - (F4b+F5b) - F8b - (F9ep_b + F9em_b) ;  %%% B18 Bacteria pool
dB(19)= F1f - (F4f+F5f) - F8f - (F9ep_f + F9em_f) ;  %%% B19 Fungi saprotrophic
dB(20)= (1-ExEM)*IS.C_myc - F8am - Resp_mycAM; %%% B20 AM-Mycorrhizal
dB(21)= (ExEM)*IS.C_myc - F8em - Resp_mycEM - (F9ep_em + F9em_em)  ; %%% B21 EM-Mycorrhizal
dB(22)= Few - RespE - M_ew; %% B22 Earthworm -C
%%%%%%
%%%%%%%%%%%%%% Nitrogen [gN/ gsoil  d]
%%%%%% NITROGEN POOL
dB(23)=  IS.N_sur_lit-sum(N_sur_flux); %%% B23 Nitrogen Above-ground Litter
dB(24)=  IS.N_wod_lit-N_wod_flux;%%% B24 Nitrogen Above-ground Woody
dB(25)=  IS.N_ssr_lit-sum(N_ssr_flux);%%% B25 Nitrogen Below-ground Litter
dB(26)= som_N_preimm - Norg_imm; %%% B26 Nitrogen SOM
dB(27)=  (1-lambda_N)*N_decomp_b-N_mort_b + (N_imm_fix_b*dtd+Nimm_b-Nmin_b)/dtd - imm_cut(1);  %%% B27 Nitrogen Bacteria
dB(28)=  (1-lambda_N)*N_decomp_f-N_mort_f + (N_imm_fix_f*dtd+Nimm_f-Nmin_f)/dtd - imm_cut(2); %%% B28 Nitrogen Fungi S
dB(29)=  (Nimm_am-Nmin_am)/dtd-N_mort_am-imm_cut(3); %% B29 Nitrogen AM Mycorrhizal
dB(30)=  (Nimm_em-Nmin_em)/dtd-N_mort_em-imm_cut(4); %% B30 Nitrogen EM Mycorrhizal
dB(31)=  Incoming_N + (1-f_org_lea)*Ip_lea_nit  +  NH4_in  - NO3flx  - NH4_Uptake - LEAK_NH4 - VOL;  %%% B31 Nitrogen Ione Ammonium NH4+
dB(32)=  NO3_in + NO3flx - N2Oflx_nit - NO3_Uptake - LEAK_NO3 - DENITflx;%%% B32 Nitrogen Nitrate NO3-
dB(33) = f_org_lea*Ip_lea_nit + lambda_N*(N_decomp_b+N_decomp_f) - LEAK_DON - S_DON;  %%% B33 DON
dB(34) = FewN; %%  B34 Earthworm -N
%%%%%%
%%%%%%%%%%%%%% Phosphorus [gP/gsoil  d]
%%%%%% PHOSPHORUS POOL
dB(35)=  IS.P_sur_lit-sum(P_sur_flux); %%%B35 phosphorus Above-ground Litter
dB(36)=  IS.P_wod_lit-P_wod_flux;%%%B36 phosphorus Above-ground Woody
dB(37)=  IS.P_ssr_lit-sum(P_ssr_flux);%%% B37 phosphorus Below-ground Litter
dB(38)= P_som_preimm - Porg_imm; %%% B38 phosphorus SOM
dB(39)=  (1-lambda_P)*P_decomp_b-P_mort_b + (Pimm_b - Pmin_b)/dtd - P_imm_cut(1);  %%%B39 phosphorus Bacteria
dB(40)=  (1-lambda_P)*P_decomp_f-P_mort_f + (Pimm_f - Pmin_f)/dtd - P_imm_cut(2); %%% B40 phosphorus Fungi
dB(41)= (Pimm_am - Pmin_am)/dtd-P_mort_am-P_imm_cut(3);%% B41 phosphorus AM- Mycorrhizal
dB(42)= (Pimm_em - Pmin_em)/dtd-P_mort_em-P_imm_cut(4);%% B42 phosphorus EM- Mycorrhizal
dB(43)=  Incoming_P + (1-f_org_lea)*Ip_lea_pho + P_mineral_transfer - P_Uptake- Psec_exch + Ppri_exch  - LEAK_P;  %%% B43 phosphorus Mineral
dB(44)=  Tup_P-Ppri_exch;%%%  B44 phosphorus primary
dB(45)= Psec_exch - Pocc_exch  ;  %%% B45 phosphorus secondary
dB(46)=  Pocc_exch  ;%%%% B46  phosphorus occluded
dB(47) = f_org_lea*Ip_lea_pho + lambda_P*P_decomp -LEAK_DOP - S_DOP;  %%% B47 DOP
%%%%%%

% to do: disable 53 and 54 for 52
%%%%%%%%%%%%%% Potassium [gK/gsoil d]
%%%%%%%%% POTASSIUM POOL
dB(48)=  IS.K_sur_lit-sum(K_sur_flux); %% B48 Potassium Above-ground Litter
dB(49)=  IS.K_wod_lit-K_wod_flux;%%%B49 Potassium  Above-ground Woody
dB(50)=  IS.K_ssr_lit-sum(K_ssr_flux);%%% B50 Potassium  Below-ground Litter
dB(51)=  Ip_som_pot-K_decomp ; %%%B51 Potassium SOM
dB(52)=  Incoming_K+Ip_lea_pot+K_decomp-K_Uptake+Kmin_rel+Kfix_sol+Kex_sol-LEAK_K;  %%% B52 phosphorus Mineral solution
dB(53)=  -Kex_sol; %%%% B53 the exchangeable K
dB(54)=  -Kfix_sol ;  %%%% B54  the fixed or non-exchangeble K
dB(55)= Tup_K -Kmin_rel ;  %%%% B55 K in the lattice of certain primary minerals (Kl).

%%%%%%
%%%%%%%%%%%%%% Calcium [gCa/gsoil d]
%%%%%%%%% CALCIUM POOL
dB(56)=  IS.Ca_sur_lit-sum(Ca_sur_flux); %% B56 Calcium Above-ground Litter
dB(57)=  IS.Ca_wod_lit-Ca_wod_flux;%%%B57 Calcium  Above-ground Woody
dB(58)=  IS.Ca_ssr_lit-sum(Ca_ssr_flux);%%% B58 Calcium  Below-ground Litter
dB(59)=  Ip_som_cal-Ca_decomp ; %%%B59 Calcium SOM
dB(60)=  Incoming_Ca+Ip_lea_cal+Ca_decomp-Ca_Uptake+Camin_rel+Cafix_sol+Caex_sol-LEAK_Ca;  %%% B60 Calcium Mineral solution
dB(61)=  -Caex_sol; %%%% B61 the exchangeable Ca
dB(62)=  -Cafix_sol ;  %%%% B62  the fixed or non-exchangeble Ca
dB(63)= Tup_Ca -Camin_rel ;  %%%% B63 Ca in the lattice of certain primary minerals.

%%%%%%%%%%%%%% Magnesium [gMg/gsoil d]
%%%%%%%%% MAGNESIUM POOL
dB(64)=  IS.Mg_sur_lit-sum(Mg_sur_flux); %% B64 Magnesium Above-ground Litter
dB(65)=  IS.Mg_wod_lit-Mg_wod_flux;%%%B65 Magnesium  Above-ground Woody
dB(66)=  IS.Mg_ssr_lit-sum(Mg_ssr_flux);%%% B66 Magnesium  Below-ground Litter
dB(67)=  Ip_som_mag-Mg_decomp ; %%%B67 Magnesium SOM
dB(68)=  Incoming_Mg+Ip_lea_mag+Mg_decomp-Mg_Uptake+Mgmin_rel+Mgfix_sol+Mgex_sol-LEAK_Mg;  %%% B68 Magnesium Mineral solution
dB(69)=  -Mgex_sol; %%%% B69 the exchangeable Mg
dB(70)=  -Mgfix_sol ;  %%%% B70  the fixed or non-exchangeble Mg
dB(71)= Tup_Mg -Mgmin_rel ;  %%%% B71 Mg in the lattice of certain primary minerals.

%%%%%%%%%%%%%% Silicate [gSi/gsoil d]
%%%%%%%%% SILICATE POOL
dB(72)=  IS.Si_sur_lit-sum(Si_sur_flux); %% B72 Silicate Above-ground Litter
dB(73)=  IS.Si_wod_lit-Si_wod_flux;%%%B73 Silicate  Above-ground Woody
dB(74)=  IS.Si_ssr_lit-sum(Si_ssr_flux);%%% B74 Silicate  Below-ground Litter
dB(75)=  Ip_som_sil-Si_decomp ; %%%B75 Silicate SOM
dB(76)=  Incoming_Si+Ip_lea_sil+Si_decomp-Si_Uptake+Simin_rel+Sifix_sol+Siex_sol-LEAK_Si;  %%% B76 Silicate Mineral solution
dB(77)=  0; %%%% B77 the exchangeable Si
dB(78)=  -Sifix_sol ;  %%%% B78  the fixed or non-exchangeble Si
dB(79)= Tup_Si -Simin_rel ;  %%%% B79 Si in the lattice of certain primary minerals.

%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
dB= dB*(ZBIOG*rsd); %% [gX / m2 day]
R_litter_sur = R_litter_sur*(ZBIOG*rsd); %% %% [gC/ m2 day]
R_litter = R_litter*(ZBIOG*rsd); %% %% [gC/ m2 day]
R_microbe = R_microbe*(ZBIOG*rsd); %% %% [gC/ m2 day]
R_bacteria = (F4b+F5b)*(ZBIOG*rsd); %% %% [gC/ m2 day]
RmycAM = Resp_mycAM*(ZBIOG*rsd) ; %% [gC/ m2 day]
RmycEM = Resp_mycEM*(ZBIOG*rsd) ; %% [gC/ m2 day]
Prod_B=F1b*(ZBIOG*rsd);  %% [gC/ m2 day]
Prod_F=F1f*(ZBIOG*rsd);  %% [gC/ m2 day]
R_ew = RespE*(ZBIOG*rsd); %% %% [gC/ m2 day]
Min_N = (NH4_in+NO3_in)*(ZBIOG*rsd); %% %% [gN/ m2 day]
Min_P = P_mineral_transfer*(ZBIOG*rsd); %% %% [gP/ m2 day]
VOL=VOL*(ZBIOG*rsd);
N2flx=N2flx*(ZBIOG*rsd);
N2Oflx=N2Oflx*(ZBIOG*rsd);
return
end

function flux = limit_outgoing_fluxes(pool,source,flux,dtd)
available = pool/dtd+source;
if available < 0
    error('BIOGEOCHEMISTRY_DYNAMIC3:NegativeDonorBase', ...
        'Daily source terms leave a negative donor inventory.');
end
withdrawable = max(0,available-16*eps(max(1,abs(available))));
total_flux = sum(flux);
if total_flux > withdrawable
    flux = flux*(withdrawable/total_flux);
end
end

function flux = proportional_flux(carbon_flux,nutrient_pool,carbon_pool)
if carbon_flux <= 0 || nutrient_pool <= 0 || carbon_pool <= 0
    flux = 0;
else
    carbon_fraction = min(realmax('double'),carbon_flux/carbon_pool);
    flux = min(realmax('double'),carbon_fraction*nutrient_pool);
end
end

function flux = limit_internal_fluxes(pool,base,from,to,flux,dtd)
available = pool+dtd*base;
if any(available < 0)
    error('BIOGEOCHEMISTRY_DYNAMIC3:NegativeBasePool', ...
        'External daily sinks overdraw a mineral pool before internal exchange.');
end
withdrawable = max(0,available-16*eps(max(1,abs(available))));
outgoing = zeros(size(pool));
donor = from;
for k = 1:numel(flux)
    if flux(k) < 0
        donor(k) = to(k);
    end
    outgoing(donor(k)) = outgoing(donor(k))+abs(flux(k));
end
scale = ones(size(pool));
limited = dtd*outgoing > withdrawable;
scale(limited) = withdrawable(limited)./(dtd*outgoing(limited));
for k = 1:numel(flux)
    flux(k) = flux(k)*scale(donor(k));
end
end

function ratio = safe_ratio(numerator,denominator)
if denominator <= 0
    ratio = realmax('double');
elseif numerator <= 0
    ratio = realmin('double');
else
    ratio = min(realmax('double'), ...
        max(realmin('double'),numerator/denominator));
end
end

function fraction = safe_fraction(numerator,denominator)
if denominator <= 0
    fraction = 0;
else
    fraction = numerator/denominator;
end
end
