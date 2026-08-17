%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PFT_Prairie_Margin.m
%%
%% Tallgrass-prairie PFT pair for the 5 m field margin cells.  Run as a
%% script via Load_Cell_Parameters AFTER the base parameter file
%% (MOD_PARAM_Kantola_Maize_Soy.m), so it replaces only the vegetation
%% parameters and inherits soil hydraulics, snow and interception from
%% the base file -- the margin sits on the same soil as the field.
%%
%% PFT ASSIGNMENT (set by morphology, not by convenience)
%%
%%   HIGH  Big bluestem        Andropogon gerardii
%%         Warm-season C4 tallgrass, the dominant of the true tallgrass
%%         prairie.  Culms reach 2.0-2.5 m at flowering, overtopping
%%         every associated forb, and its dense fibrous root system
%%         reaches 2-3 m.  It is the taller, deeper-rooted, canopy-
%%         forming species, so it takes the HIGH (overstory) slot.
%%
%%   LOW   Purple coneflower   Echinacea purpurea
%%         C3 perennial forb, 0.6-1.2 m, with a deep fibrous root system.
%%         It grows beneath the grass canopy, so it takes the LOW
%%         (understory) slot.
%%
%%   Both species occupy crown 1 (Ccrown = [1 0]) as an overstory-grass /
%%   understory-forb pair sharing the same ground area.  Slot 2 is the
%%   base file's soybean slot: retained structurally so no caller needs
%%   re-indexing, but carries zero cover and zero rooting depth.
%%
%%   The species identity, pathway and high/low assignment are
%%   morphology based. The exact physiology and ZR95 values below are
%%   model priors, not measurements from the Kantola field.
%%
%% WHY A GRASS + FORB PAIR RATHER THAN A GRASS MONOCULTURE
%%   The margin treatments are being compared for their effect on N2O
%%   and NO3.  A C4 grass alone gives deep roots but a single phenology;
%%   adding a C3 forb spreads uptake across the season (cool-season forb
%%   activity before and after the C4 grass peak), which is the
%%   mechanism by which a diverse margin can intercept nitrogen the
%%   crop has left behind.  Both are non-leguminous, so neither fixes
%%   N and the margin adds no new N to the budget.
%%
%% References
%%   Weaver 1968, Prairie Plants and Their Environment (root depths)
%%   Jackson et al. 1996, Oecologia 108:389-411 (root profiles)
%%   Craine et al. 2002, Functional Ecology 16:563-574 (prairie traits)
%%   Tjoelker et al. 2005, New Phytologist 167:493-508 (C4 grass gas exchange)
%%   Swinnen et al. 1994, Plant and Soil 164:267-277 (grass exudation)
%%   Guo & Gifford 2002, Oecologia 131:341-350 (root C:N)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% --- Canopy cover: one prairie community in crown 1, slot 2 disabled
%% The margin is fully vegetated, so the bare fraction the base file
%% carried for a row crop goes to zero.  Land-cover fractions must still
%% sum to one or Check_Land_Cover_Fractions rejects the cell.
Ccrown = [1.0 0.0];
Cbare  = 1 - sum(Ccrown) - Cwat - Curb - Crock;

%% --- Root depth  [mm]
%%   Big bluestem roots can extend beyond 2 m. E. purpurea is fibrous and
%%   shallower. These ZR95 values preserve that ordering; neither is a
%%   site measurement.
ZR95_H  = [2500 0];
ZR50_H  = [NaN  NaN];
ZRmax_H = [NaN  NaN];
ZR95_L  = [1200 0];
ZR50_L  = [NaN  NaN];
ZRmax_L = [NaN  NaN];

%% Honour the active depth the caller injected, exactly as the crop
%% parameter file does. Run_One_Subfield sets MOD_PARAM_input.Zbio so every
%% column in the domain shares one depth; MAIN_FRAME_BG_ERW_2D requires
%% that, because lateral exchange is undefined across differing depths.
%% Deriving the depth from the deepest root here returned 2500 mm for big
%% bluestem while the crop columns stayed at 300 mm, which aborted runs 7-9
%% with ActiveDepthMismatch. 2500 mm is also the depth that diluted DOC per
%% gram of soil and starved the microbial pools, so it must not come back
%% in through the margin. Roots still reach 2500 mm; only the
%% biogeochemically active zone is held at the common depth.
%% MOD_PARAM_input is a table, so read the column the same way the crop
%% parameter file does. isfield() is always false on a table and would
%% silently send us down a fallback path.
Zbio = MOD_PARAM_input.Zbio;
if ~(isscalar(Zbio) && isfinite(Zbio) && Zbio > 0)
    error('PFT_Prairie_Margin:NoActiveDepth', ...
        ['The caller must supply a single positive MOD_PARAM_input.Zbio. ' ...
         'Deriving it from the deepest root here is what made the margin ' ...
         'columns disagree with the crop columns.']);
end
Bio_Zs = Evaporation_layers(Zs,Zbio);

%% Recompute root fraction profiles for the new depths
[RfH_Zs,RfL_Zs] = Root_Fraction_General(Zs,CASE_ROOT, ...
    ZR95_H,ZR50_H,ZR95_L,ZR50_L,ZRmax_H,ZRmax_L);

%% --- Photosynthetic pathway
CT_H = [4 4];        %% C4  big bluestem
CT_L = [3 3];        %% C3  purple coneflower

%% --- Quantum efficiency  [umolCO2/umolPhotons]
FI_H = [0.055 0.055];   %% C4
FI_L = [0.081 0.081];   %% C3

%% --- Stomatal / mesophyll kinetics
Do_H   = [1000 1000];  Do_L   = [800  800];
a1_H   = [7    7];     a1_L   = [7    7];
go_H   = [0.01 0.01];  go_L   = [0.01 0.01];
DSE_H  = [0.649 0.649];DSE_L  = [0.649 0.649];
Ha_H   = [72   72];    Ha_L   = [72   72];
gmes_H = [Inf  Inf];   gmes_L = [Inf  Inf];
rjv_H  = [1.97 1.97];  rjv_L  = [2.10 2.10];

%% --- Maximum carboxylation  [umol/m2/s]
%%   Below the maize 78 of the base file: prairie species run lower leaf
%%   N and lower Vmax than a fertilised annual crop.
Vmax_H = [50 0];
Vmax_L = [45 0];

%% --- Leaf traits
Sl_H = [0.020 0.020];   %% [m2/gC] narrow C4 grass blade
Sl_L = [0.030 0.030];   %% [m2/gC] broader forb lamina
Nl_H = [35 35];         %% [gC/gN] low-N C4 grass foliage
Nl_L = [22 22];         %% [gC/gN] higher-N forb foliage
[Stoich_H(1)] = Veg_Stoichiometric_Parameter(Nl_H(1));
[Stoich_H(2)] = Veg_Stoichiometric_Parameter(Nl_H(2));
[Stoich_L(1)] = Veg_Stoichiometric_Parameter(Nl_L(1));
[Stoich_L(2)] = Veg_Stoichiometric_Parameter(Nl_L(2));
KnitH = [0.20 0.20];
KnitL = [0.15 0.15];
mSl_H = [0.0 0.0];
mSl_L = [0.0 0.0];
d_leaf_H = [1.5 1.5];   %% [cm] narrow grass blade
d_leaf_L = [6.0 6.0];   %% [cm] broad forb leaf

%% --- Plant type: 2 = perennial herbaceous PFT, no annual sowing or harvest
%% The model has no separate perennial-forb state, so type 2 is also used
%% for coneflower.
aSE_H = [2 2];
aSE_L = [2 2];

%% --- Respiration
r_H  = [0.040 0.040];   %% [gC/gN d] perennial resource conservation
r_L  = [0.050 0.050];
gR_H = [0.25 0.25];
gR_L = [0.25 0.25];

%% --- Turnover  [1/d]
drn_H = [1/730 1/730];  %% 2-year fine-root turnover (perennial grass)
drn_L = [1/548 1/548];  %% 1.5-year (forb fine roots)
dsn_H = [0 0];
dsn_L = [0 0];
Klf_H = [1/30 1/30];    %% dead-leaf fall
Klf_L = [1/25 1/25];
Mf_H  = [0 0];   Mf_L = [0 0];
Wm_H  = [0 0];   Wm_L = [0 0];

%% --- Mortality and stress
dd_max_H = [1/50 1/50];
dd_max_L = [1/50 1/50];
dc_C_H   = [7/365 7/365];
dc_C_L   = [7/365 7/365];
Tcold_H  = [-5 -5];     %% both are frost-tolerant temperate perennials
Tcold_L  = [-5 -5];

%% --- Phenology
%%   The C3 forb leafs out earlier and senesces later than the C4 grass;
%%   that offset is the seasonal-complementarity mechanism above.
Tlo_H = [13 13];        Tlo_L = [8  8];
Tls_H = [NaN NaN];      Tls_L = [NaN NaN];
LDay_min_H = [12.0 12.0];  LDay_min_L = [11.0 11.0];
LDay_cr_H  = [11.5 11.5];  LDay_cr_L  = [10.5 10.5];
mjDay_H    = [250 250];    mjDay_L    = [280 280];
age_cr_H   = [140 140];    age_cr_L   = [120 120];
dmg_H      = [45 45];      dmg_L      = [40 40];
Bfac_lo_H  = [0.95 0.95];  Bfac_lo_L  = [0.95 0.95];
Bfac_ls_H  = [NaN NaN];    Bfac_ls_L  = [NaN NaN];
PAR_th_H   = [NaN NaN];    PAR_th_L   = [NaN NaN];
LAI_min_H  = [0.1 0.1];    LAI_min_L  = [0.05 0.05];
Trr_H      = [3.5 3.5];    Trr_L      = [3.0 3.0];
eps_ac_H   = [0.5 0.5];    eps_ac_L   = [0.5 0.5];

%% --- Allocation: perennial prairie invests heavily belowground
LtR_H = [1.8 1.8];      %% leaf-to-root ratio maximum, C4 bunchgrass
LtR_L = [1.5 1.5];      %% forb with belowground storage
fab_H = [1.0 1.0];  fbe_H = [0.0 0.0];
fab_L = [1.0 1.0];  fbe_L = [0.0 0.0];
ff_r_H = [0.05 0.05];   %% little allocation to reproduction
ff_r_L = [0.10 0.10];

%% --- Crop-only parameters: not applicable to a perennial stand
soCrop_H     = [NaN NaN];  soCrop_L     = [NaN NaN];
Sl_emecrop_H = [NaN NaN];  Sl_emecrop_L = [NaN NaN];
MHcrop_H     = [NaN NaN];  MHcrop_L     = [NaN NaN];

%% --- Water relations  [MPa]
%%   Big bluestem is markedly drought-tolerant; the forb closes earlier.
Psi_sto_00_H = [-1.0 -1.0];  Psi_sto_50_H = [-3.2 -3.2];
Psi_sto_00_L = [-0.8 -0.8];  Psi_sto_50_L = [-2.5 -2.5];
PsiL00_H = [-1.5 -1.5];      PsiL50_H = [-4.0 -4.0];
PsiL00_L = [-1.2 -1.2];      PsiL50_L = [-3.2 -3.2];
Kleaf_max_H = [5 5];         Kleaf_max_L = [5 5];
Cl_H = [1200 1200];          Cl_L = [1200 1200];
Axyl_H = [0.0 0.0];          Axyl_L = [0.0 0.0];
Kx_max_H = [80000 80000];    Kx_max_L = [80000 80000];
PsiX50_H = [-5.0 -5.0];      PsiX50_L = [-4.0 -4.0];
Cx_H = [150 150];            Cx_L = [150 150];

%% --- Growth water limitation  [MPa]
PsiG50_H = [-1.5 -1.5];  PsiG99_H = [-4.0 -4.0];
PsiG50_L = [-1.2 -1.2];  PsiG99_L = [-3.2 -3.2];
gcoef_H  = [3.5 3.5];    gcoef_L  = [3.5 3.5];

%% --- Optical properties
%%   14 = C4 grass (big bluestem); 13 = C3 grass/herbaceous (coneflower).
%%   The base file used 16 (Crop2) for maize and soybean.
[PFT_opt_H(1)] = Veg_Optical_Parameter(14);
[PFT_opt_H(2)] = Veg_Optical_Parameter(14);
[PFT_opt_L(1)] = Veg_Optical_Parameter(13);
[PFT_opt_L(2)] = Veg_Optical_Parameter(13);
OM_H = [1 1];
OM_L = [1 1];

%% --- Interception
Sp_LAI_H_In = [0.2 0.2];
Sp_LAI_L_In = [0.3 0.3];   %% broader forb leaf holds more water

%% --- Exudation and biological fixation
%% Exudation_Parameter's argument is the biological-fixation flag.
[ParEx_H(1)] = Exudation_Parameter(0);
[ParEx_H(2)] = Exudation_Parameter(0);
[ParEx_L(1)] = Exudation_Parameter(0);
[ParEx_L(2)] = Exudation_Parameter(0);
assert(all([ParEx_H.bfix] == 0) && all([ParEx_L.bfix] == 0), ...
    'PFT_Prairie_Margin:UnexpectedBNF', ...
    'Big bluestem and purple coneflower must both be non-fixing PFTs.');

%% --- Management: an established perennial margin is neither sown nor
%%     harvested during the simulation, and receives no fertiliser.
%% Overwrite fields on the management structs the base file already
%% built.  Re-calling Vegetation_Management_Parameter here would fail:
%% the base file adds crop-rotation fields (Date_sowing, Crop_B, ...) to
%% its structs, so a fresh one no longer matches field for field.
for islot = 1:2
    Mpar_H(islot).Date_sowing     = [];
    Mpar_H(islot).Date_harvesting = [];
    Mpar_H(islot).Crop_B          = [0 0];
    Mpar_H(islot).fract_left      = 1;
    Mpar_H(islot).fract_left_fr   = 1;
    Mpar_H(islot).fract_left_AB   = 1;
    Mpar_H(islot).fract_left_BG   = 1;
    Mpar_L(islot).Date_sowing     = [];
    Mpar_L(islot).Date_harvesting = [];
    Mpar_L(islot).Crop_B          = [0 0];
    Mpar_L(islot).fract_left      = 1;
    Mpar_L(islot).fract_left_fr   = 1;
    Mpar_L(islot).fract_left_AB   = 1;
    Mpar_L(islot).fract_left_BG   = 1;
end
Mpar_H(1).Crop_crown = ones(1,12);    %% present all year, every year
Mpar_H(2).Crop_crown = zeros(1,12);
Mpar_L(1).Crop_crown = ones(1,12);
Mpar_L(2).Crop_crown = zeros(1,12);
%% Litter returned to soil  [N P K Ca Mg Si] of residue
Mpar_H(1).NPK_res_ini = [10 1.5 12 3.0 2.0 8.0];   %% C4 grass litter, low N
Mpar_L(1).NPK_res_ini = [18 2.5 18 8.0 3.0 1.0];   %% forb litter, higher N and Ca
Mpar_H(2).NPK_res_ini = [0 0 0 0 0 0];
Mpar_L(2).NPK_res_ini = [0 0 0 0 0 0];

%% --- No fertiliser on the margin in any treatment
if exist('B_IO','var') && isstruct(B_IO)
    B_IO.FertN  = zeros(1,366);
    B_IO.FertP  = zeros(1,366);
    B_IO.FertK  = zeros(1,366);
    B_IO.FertCa = zeros(1,366);
    B_IO.FertMg = zeros(1,366);
    B_IO.FertSi = zeros(1,366);
end

%% --- Initial stand state: an established margin, not a bare seedbed
LAI_H(1,:)  = [0.3 0.0];
LAI_L(1,:)  = [0.2 0.0];
hc_H(1,:)   = [1.8 0.0];    %% [m] mature big bluestem culm height
hc_L(1,:)   = [0.8 0.0];    %% [m] coneflower flowering height
SAI_H(1,:)  = [0.05 0.0];
SAI_L(1,:)  = [0.01 0.0];
PHE_S_H(1,:)= [1 0];
PHE_S_L(1,:)= [1 0];
%% Live leaves / sapwood / fine roots / reserve / fruit / heartwood / dead leaves
B_H(1,1,:)  = [15 0 250 120 0 0 10 0];
B_L(1,1,:)  = [10 0 120  80 0 0  5 0];
B_H(1,2,:)  = zeros(1,8);
B_L(1,2,:)  = zeros(1,8);
AgeL_H(1,:) = [30 0];
AgeL_L(1,:) = [30 0];

fprintf(['[PFT_Prairie_Margin] big bluestem (C4, HIGH, ZR95=%d mm) + ' ...
    'purple coneflower (C3, LOW, ZR95=%d mm), perennial, unfertilised\n'], ...
    ZR95_H(1), ZR95_L(1));
