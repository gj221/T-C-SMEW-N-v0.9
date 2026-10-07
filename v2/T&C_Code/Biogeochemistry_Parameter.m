%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction Biogeochemistry_Parameter     %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[BiogeoPar]=Biogeochemistry_Parameter(opt_cons_CUE,KNmaxd_override,overrides)
if nargin < 2
    KNmaxd_override = [];
end
if nargin < 3
    overrides = [];
end
if ~isempty(KNmaxd_override) && ...
        ~(isnumeric(KNmaxd_override) && isreal(KNmaxd_override) && ...
          isscalar(KNmaxd_override) && isfinite(KNmaxd_override) && ...
          KNmaxd_override >= 0)
    error('Biogeochemistry_Parameter:InvalidKNmaxd', ...
        'KNmaxd must be a finite nonnegative scalar in d^-1.');
end
%%%%%%%%%%%%%%
%%%% Decay Constant on a 40�C basis (Kirschbaum 2002)
BiogeoPar.k_met_sur = 1/12.5; %%[1/day]
BiogeoPar.k_str_sur = 1/46; %% %[1/day]
BiogeoPar.k_met_ssr =  1/10;
BiogeoPar.k_str_ssr = 1/37;
BiogeoPar.k_wod_sur = 1/150; %% [1/92 - 1/460]
%%%%%%%%%%%%%%%%%%%%
%%%%%% Respiration Coefficient for Litter [-]  Parton et al 1987 
%rr_met_sur=0.55;
%rr_str_sur=0.45;
%rr_wod_sur=0.55;
%rr_met_ssr=0.55;
%rr_str_ssr=0.55;
%%%%%%% Critical Ratios
BiogeoPar.R_CN_bac = 5.2; % 7 ;%  %%[C/N]  8.6
BiogeoPar.R_CN_fun = 6.5; % 12 ;% 
BiogeoPar.R_CN_myc = 18;
BiogeoPar.R_CN_ew = 10; 
BiogeoPar.R_CP_bac = 16; % 35; %%[C/P] 60
BiogeoPar.R_CP_fun = 40;% 80; %% 120 
BiogeoPar.R_CP_myc = 120; % 120;
%%%%%%%%%%%%%%%%%%
%%% Parameters  %%% Wang et al 2012  Waring-Averill Ecol. Let 2013 
BiogeoPar.gd =0.2;% 0.3;% 0.3;% [0.2-0.8] fraction of dead B allocated to D
BiogeoPar.fd =0.4;% 0.5; % [0.2-0.8]  fraction of decomposed P allocated to D
%Ecm =0.47; % carbon use efficiency [-]
%Ecf =0.47;
BiogeoPar.Vdb = 0.04; %0.120; %0.010;  % [0.0024-0.060] maximum specific uptake rate of DOC for growth of B [gC-D/gC-B day]
BiogeoPar.Vdf = 0.02; %0.025;% 0.012;  % 
BiogeoPar.Kdb = 0.00026; % [0.00014 - 0.00038] half-saturation constant of uptake of DOC for growth of B [ g C/g soil]
BiogeoPar.Kdf = 0.00026;
%% Legacy T&C-SMEW values, under which Data_Kantola_BG.mat was spun up. Its
%% SI Table S2 lists 0.005 and 0.002 as the originals.
BiogeoPar.mrb =  10.0*1e-3; % 15.0*1e-3; % [ 1.20*1e-3 - 30*1e-3] specific maintenance factor or rate of B [ g C /g C day]
BiogeoPar.mrf =  5.0*1e-3; % 5.0*1e-3;
BiogeoPar.mram = 1.5*1e-3; % 1.5*1e-3; specific maintenance factor of AM 
BiogeoPar.mrem = 1.5*1e-3; % 1.5*1e-3; specific maintenance factor of EM 


BiogeoPar.Vpc =200;% 95; % 120; %% [4.8 - 792] maximum specific decomposition rate for POC-C by EP [gC-P/gC-EP day]
BiogeoPar.Vpl =23; % 8; % 23;% 
BiogeoPar.Vm = 100; % 48;% 24;  %[1.2 -528] maximum specific decomposition rate for MOC by [gC-M /gC-EM day]
BiogeoPar.Kpl =0.05; % [0.023 - 0.094] half-saturation constant for decomposition of POC-L [gC / gsoil]
BiogeoPar.Kpc =0.05; 
BiogeoPar.Km =0.25; % [0.057 - 0.95] half-saturation constant for decomposition of MOC [ gC /gsoil];
BiogeoPar.rem = 0.018; % 0.025; % [0.0075 - 0.075] turnover rate of EM [g C/g C day]
BiogeoPar.rep = 0.018; % 0.025; %  [0.0075 - 0.075] turnover rate of EP [g C/g C day]
BiogeoPar.pepb = 0.012; % [0.0031 - 0.031] fraction of mR for production of EP
BiogeoPar.pemb = 0.005;% [0.0031 - 0.031] fraction of mR for production of EM
BiogeoPar.pepf = 0.006; % [0.0031 - 0.031] fraction of mR for production of EP
BiogeoPar.pemf = 0.010;% [0.0031 - 0.031] fraction of mR for production of EM
BiogeoPar.pepem = 0.006;
BiogeoPar.pemem = 0.010; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BiogeoPar.Ecf = 0.27;% 0.47; %%% Sinsbaugh et al 2016 
BiogeoPar.Ecb = 0.27;% 0.47;
%%% Sinsabaguh et al 2013 ; Manzoni et al 2008
BiogeoPar.CUEmax= 0.6;
BiogeoPar.Ae=1; %% (not used)

%%% Temperature effects on CUE
%Tref=20;
BiogeoPar.T1= 12;  %%%%%% Temperature Effects Kinetiks
BiogeoPar.T2= 20;  %%% Temperature effects on CUE
BiogeoPar.T3= 15; % Temperature Earthworms
BiogeoPar.T4= 40; %%% Ref. Temp. Litter Decoomposition
%%%%%%%%
%%%%%
BiogeoPar.mtem = -0.008; %%[�C] [0 - 0.016] Allison et al 2010; Li et al 2014
if opt_cons_CUE == 1
    BiogeoPar.rr_met_sur = 0.55;
    BiogeoPar.rr_str_sur=  0.45;
    BiogeoPar.rr_wod_sur=  0.45;
    BiogeoPar.rr_met_ssr= 0.55;
    BiogeoPar.rr_str_ssr= 0.55;
end
%%%%%%%%%%%%%% Leaching Coefficients 
BiogeoPar.lambda_C=0.0015;  %%% 0 up to 2500 mm/yr 0.5 at 4000 mm/yr linear between 
BiogeoPar.lambda_N=0.0015;
BiogeoPar.lambda_P=0.0001;
BiogeoPar.lambda_K=0.90;
BiogeoPar.f_org_lea=1; %% Fraction of Organic leaching 

BiogeoPar.lambda_Ca=0.90;
BiogeoPar.lambda_Mg=0.90;
BiogeoPar.lambda_Si=0.90;

%%%%%%%%%%%%%%%%%% Earthworm parameters 
%BiogeoPar.kse =  5e-004 ; % [gC/ gsoil] half-saturation constant for consumption of food %
%BiogeoPar.VmaxE = 0.038 ; %[gC / gC-EW d] Ref earthworm growth 15�C
BiogeoPar.VmaxE = 0.0050 ; %%  [gC / gC-food d]
BiogeoPar.rmanEa= 0.0300 ; %% [gC / gC-EW d] Maint. respiration 15�C 
BiogeoPar.rmanEp= 0.0146 ; %% [gC / gC-EW d] Maint. respiration 15�C 
BiogeoPar.dew = 0.010; % %% [gC / gC-EW d] 
BiogeoPar.fabsE =0.15; % [-] fraction of absorbed food 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Volatilization of Ammonia NH3
BiogeoPar.kv = 8.6*1e-5; %% [1/d]
%%%%%%% Nitrification (Aerobic)
BiogeoPar.KNmaxn = 0.086; %% [1/d] % Dickinson et al., 2002
%%%%%%% Denitrification (Anaerobic)
BiogeoPar.Kc = 17.5; %% [gC/m3] DOC half-saturation, Li et al., 1992
BiogeoPar.Kn = 82.5; %% [gN/m3] inactive NOx half-saturation, Li et al., 1992
%%%%%% Historical denitrification rate constant, not the literature value.
%% Dickinson et al. (2002) give 2.5e-6 s^-1 = 0.217 d^-1. That value was never
%% fitted to this site and overestimates Energy Farm N2O by 5.4x against the
%% chamber record digitised from Val Martin et al. (2023), figure 3a.
%% Calibrated on the 2017 and 2018 maize seasons over the chamber sampling
%% windows (DOY 114-263 and 117-243), the two years the model rotation and the
%% field trial share as maize: 8601 gN/ha observed against 9023 modelled,
%% a ratio of 1.05.
%% That aggregate conceals two compensating per-year errors: 2017 modelled
%% 3774 against 1561 observed (2.4x over) and 2018 modelled 5249 against 7040
%% (0.75x under). The model damps interannual amplitude, giving a 2018/2017
%% ratio of 1.39 against 4.51 observed and 2.23 for CLM5. The calibration is
%% therefore sound for multi-year cumulative fluxes and unreliable for any
%% single year.
%% Cost of the calibration: NO3 leaching rises from 42 to 68 kg N/ha/yr, the
%% upper end for tile-drained Midwest maize-soy. That is the binding
%% constraint if drain nitrate data become available.
%% The basalt and calcite N2O reductions are insensitive to this constant
%% (44% at 0.217; 39% and 43% at 0.0271), so the treatment comparison does
%% not rest on the calibration.
%% Re-checked 2026-09-03 after the CEC_tot_in correction 8.6 -> 19.5 (USDA SSURGO,
%% Dana Ap CEC7): the same windows give a ratio of 1.01 to observed, so the
%% constant stands. The basalt reduction moves to 27% at CEC 19.5; that is the
%% CEC effect on soil pH, not a KNmaxd effect.
%% Those checks used the earlier Se^(1/L) moisture response. The default is
%% retained as the old reference only; it is not calibrated for Se^2. The
%% corrected-moisture sweep supplies KNmaxd explicitly for each run.
BiogeoPar.KNmaxd= 0.0048; %% [1/d] calibrated to N2O/obs=1 (knscore, 2026-09-08); Dickinson et al. 2002 default is 0.217
if ~isempty(KNmaxd_override)
    BiogeoPar.KNmaxd = KNmaxd_override;
end
%%%%%% N2O pathway - pH effects
%%%%%% Parton et al 1997
%%%%%%
BiogeoPar.K1p= 1/600; %%% [1/day]
BiogeoPar.K2p= 1/13500; %% [1/day]
BiogeoPar.K3p= 1/(30*1e+6); %% [1/day]
BiogeoPar.K4p= 1/(4.38*1e+6);% 1/300000; %% [1/day]
%%%%%
%%%%%% Myself 
BiogeoPar.Kdk=1/2;%%%% [1/day] desorption exchangeable
BiogeoPar.Kak= 1/2; %% [1/day] absorption exchangeable
%%%%
BiogeoPar.K1k= 1/600 ; % []; %%% [1/day] Fixation of  mineral K
BiogeoPar.K2k= 1/13500; % []; %% [1/day]  Release of fixed K 
BiogeoPar.K3k= 1/350000; %% [1/day] Release primary K 

%%%%%% Calcium - ref? 
BiogeoPar.Kdca=1/2;%%%% [1/day] desorption exchangeable
BiogeoPar.Kaca= 1/2; %% [1/day] absorption exchangeable
%%%%
BiogeoPar.K1ca= 1/600 ; % []; %%% [1/day] Fixation of  mineral Ca
BiogeoPar.K2ca= 1/13500; % []; %% [1/day]  Release of fixed Ca 
BiogeoPar.K3ca= 1/140000; %% [1/day] Release primary Ca
%%%%%% Magnesium - ref? 
BiogeoPar.Kdmg=1/2;%%%% [1/day] desorption exchangeable
BiogeoPar.Kamg= 1/2; %% [1/day] absorption exchangeable
%%%%
BiogeoPar.K1mg= 1/600 ; % []; %%% [1/day] Fixation of  mineral Mg
BiogeoPar.K2mg= 1/13500; % []; %% [1/day]  Release of fixed Mg
BiogeoPar.K3mg= 1/175000; %% [1/day] Release primary Mg
%%%%%% Silicate - ref? 
BiogeoPar.Kdsi=1/2;%%%% [1/day] desorption exchangeable
BiogeoPar.Kasi= 1/2; %% [1/day] absorption exchangeable
%%%%
BiogeoPar.K1si= 1/600 ; % []; %%% [1/day] Fixation of  mineral Si
BiogeoPar.K2si= 1/13500; % []; %% [1/day]  Release of fixed Si
BiogeoPar.K3si= 1/7000000; %% [1/day] Release primary Si


%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Solubility coefficients for leakage and passive uptake 
BiogeoPar.aNH4 =0.05;  %%0.1 CLM4
BiogeoPar.aNO3 = 1;  
BiogeoPar.aP = 0.005;
BiogeoPar.aK = 1; 
BiogeoPar.aDOC = 0.04; 
BiogeoPar.aDON= 1;
BiogeoPar.aDOP= 1;

%%% P(60) and P(68) hold only dissolved Ca and Mg since the ERW module took
%%% over the exchange partition (ERW_Daily_Step_2D writes them back from the
%%% pore water), so they leach and are taken up at their concentration like
%%% NO3 and K. The legacy 0.4 belonged to a lumped pool and kept 60 % of the
%%% weathered cations in the box (budget of package m, 1 October 2026).
BiogeoPar.aCa = 1;
BiogeoPar.aMg = 1;
BiogeoPar.aSi = 0.001;

%%%%%%%%%%%%%%% Uptake coefficients 
%%% Root EM AM  
%%%--< These are dependent on the exudation/mycorrhiza allocations costs 
% Cex_root = 5; 
% kc1 = 10; kc2 = 300; kc3 = 50 ;
% kn1 = 600; kn2 = 50; kn3 = 100;
%kp1 = kn1/8 ; kp2 = kn2/8;  kp3 = kn3/8;
%kk1 = kn1/2.5; kk2 = kn2/2.5; kk3 = kn3/2.5;
%%%%%%%
% BiogeoPar.vem  = BiogeoPar.mrm/kc2; %  5e-006 ;  %%[gN /gC^2 day]
% BiogeoPar.Kvem = (kn2/kc2/365); % 45.6e-005 ; %%[gN /gC ] 
% BiogeoPar.vam =  BiogeoPar.mrm/kc3;% 3e-005 ;  %%[gN /gC^2 day]
% BiogeoPar.Kvam = (kn3/kc3/365);% 0.0055 ; %%[gN /gC ] 
% BiogeoPar.vr = (1/kc1)*Cex_root/365;% 0.0014; %%[gN /gC day]
% BiogeoPar.Kvr = (kn1/kc1/365);% 0.1644 ; %%[gN /gC ] 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kn1 = 7000; kn2 = 3.5; kn3 = 2.0;
Zr = 0.25; % [m]
rroot =  0.5*1e-3 ; % 3.3*1e-4 ;%% [0.5-6 *10^-3] [m] root radius
rmyc_am =  5*1e-6 ; % [m] hyphae radius
rmyc_em =  5*1e-6 ; % [m] hyphae radius
rho_r = 122*1000; %% [gC / m^3] Root density  Jackson et al., 1997
rho_m = 200*1000; %% [gC / m^3]
BiogeoPar.g_r = (sqrt(2)*(rroot.^2)*(Zr.^1.5)*(rho_r.^1.5));
BiogeoPar.g_em = (sqrt(2)*(rmyc_em.^2)*(Zr.^1.5)*(rho_m.^1.5));
BiogeoPar.g_am = (sqrt(2)*(rmyc_am.^2)*(Zr.^1.5)*(rho_m.^1.5));
BiogeoPar.vr =   2.3e-7; %(Cex_root/365)*BiogeoPar.g_r/kn1 ;%% [m^2/day]
BiogeoPar.vem =  2.3e-7; %BiogeoPar.g_em*BiogeoPar.mrm/kn2 ;%% [m^2/day]
BiogeoPar.vam =  2.3e-7; %BiogeoPar.g_am*BiogeoPar.mrm/kn3 ; %% [m^2/day]

%%%%%% Uptake coefficients 
%BiogeoPar.vNH4 = 1.5*1e-4  ;  %%[gN /gC day]  2.33*1e-3  5.14*1e-6 
%BiogeoPar.vNO3 =  1.5*1e-4 ;  
%BiogeoPar.vP =  1.5*1e-4 ;
%BiogeoPar.vK =  1.5*1e-4 ;
%%%%%%%%%% half-constants  for uptakes 
%BiogeoPar.KnNH4 = 0.5; %% [0.83-1.0]  [gN/m2] 
%BiogeoPar.KnNO3 = 0.2; 
%BiogeoPar.KnP=0.05;
%BiogeoPar.KnK=0.1;

%%%%%% Shape of the nitrification pH response used in BIOGEOCHEMISTRY_DYNAMIC3.

%% Run-local replacement of the twelve decomposition literals the T&C-SMEW
%% uncertainty analysis varied (Table S2 of its supplement): applied last,
%% because nothing above derives another field from them. KNmaxd has its
%% own argument; any other name, or a value outside [0, Inf), is an error.
overridable = {'Vm','Vpc','Km','rem','rep','pepf','mrb','mrf','fd', ...
    'Vpl','Vdb','Vdf'};
if ~isempty(overrides)
    if ~(isstruct(overrides) && isscalar(overrides))
        error('Biogeochemistry_Parameter:InvalidOverride', ...
            'Overrides must be a scalar struct of parameter names.');
    end
    names = fieldnames(overrides);
    for k = 1:numel(names)
        if strcmp(names{k},'KNmaxd')
            error('Biogeochemistry_Parameter:InvalidOverride', ...
                'KNmaxd is set through its own argument, not the override struct.');
        end
        if ~ismember(names{k},overridable)
            error('Biogeochemistry_Parameter:UnknownOverride', ...
                '%s is not one of the overridable literals (%s).', ...
                names{k},strjoin(overridable,', '));
        end
        value = overrides.(names{k});
        if ~(isnumeric(value) && isreal(value) && isscalar(value) && ...
                isfinite(value) && value >= 0)
            error('Biogeochemistry_Parameter:InvalidOverride', ...
                '%s must be a finite nonnegative scalar.',names{k});
        end
        BiogeoPar.(names{k}) = value;
    end
end
return


