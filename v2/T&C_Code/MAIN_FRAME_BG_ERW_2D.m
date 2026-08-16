%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%  MULTI-COLUMN MAIN FRAME  --  T&C-SMEW-N sub-plot domains %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Drives one soil-vegetation column per sub-plot cell of a
%% Build_Cell_Topology partition, coupling the columns laterally every
%% hydrologic timestep.
%%
%% RELATION TO THE SINGLE-COLUMN MAIN FRAME
%%   The science is unchanged: this script calls exactly the same unit
%%   functions as MAIN_FRAME_BG_ERW_10mins (HYDROLOGIC_UNIT, BIOGEO_UNIT,
%%   VEGGIE_UNIT, Plant_Exports, CropHeightType, ERW_function) with the
%%   same arguments.  What changes is that every state array carries a
%%   trailing column dimension and the columns exchange water and
%%   dissolved solute and heat after each hourly step.  With N_col = 1 and a
%%   single-cell topology it reproduces the single-column result, which
%%   is how it is verified (test_2d_driver_equivalence).
%%
%% WHAT IS RESOLVED IN 2D, AND WHAT IS NOT
%%   Water     lateral Buckingham-Darcy exchange between adjacent cells,
%%             per soil layer, every hour (Darcy_Lateral_Topology).
%%   Mass      advective transport of the aqueous phase with that water,
%%             upwind (Solute_Lateral_Transport).  Sorbed and solid
%%             phases do not move laterally. ERW-only Na, Al, residual
%%             anions and DIC use the same water flux in the active layer.
%%   Energy    lateral Fourier conduction and liquid-water sensible heat
%%             advection couple the soil thermal profiles.
%%   Gas       soil-air CO2 diffuses across cell faces with the existing
%%             Millington-Quirk effective diffusivity.
%%
%% METEOROLOGY
%%   The plot is flat and small, so all columns share one meteorological
%%   forcing series.  Sub-domains diverge through vegetation, soil state
%%   and treatment, not through their driving weather.
%%
%% EXPECTED IN THE CALLER WORKSPACE
%%   topo                Build_Cell_Topology struct (N_col = topo.N_cells)
%%   ColCfg              [N_col x 1] struct: .PARAM_IC .PFT_OVERRIDE
%%                       .apply_rock (logical)
%%   Parameter_Directory path holding the parameter scripts
%%   MOD_PARAM_input, ERW_PARAM_input, NN, NNd, ms, Zs, cc, dt, dth,
%%   Datam, Date, and the meteorological series, as for the 1-D frame.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N_col = topo.N_cells;
fprintf('\n=== 2D driver: %d columns over %.0f m2 ===\n', ...
    N_col, sum(topo.cell_area_m2));

if ~exist('OPT_SoilBiogeochemistry','var'); OPT_SoilBiogeochemistry = 1; end
if ~exist('OPT_ERW','var');                 OPT_ERW = 1;                 end
if ~exist('IrD','var') || isempty(IrD);     IrD = zeros(NN,1);           end

NCP = 8;    %% number of carbon pools

%% ========================================================================
%% PER-COLUMN PARAMETERS
%%   Each column's parameter script runs in its own private workspace, so
%%   columns may differ in PFT, management, rooting depth and soil.
%% ========================================================================
%% Build field by field: struct('f',M) with a non-scalar M would create a
%% struct ARRAY the size of M, not a scalar struct holding M.
ctx = struct();
ctx.MOD_PARAM_input = MOD_PARAM_input;
ctx.ERW_PARAM_input = ERW_PARAM_input;
ctx.NN = NN;   ctx.NNd = NNd;   ctx.ms = ms;
ctx.Datam = Datam;
ctx.Lat = Lat; ctx.Lon = Lon;   ctx.DeltaGMT = DeltaGMT;
ctx.t_bef = t_bef; ctx.t_aft = t_aft;
ctx.Ta = Ta;   ctx.Ca = Ca;
ctx.OPT_SoilBiogeochemistry = OPT_SoilBiogeochemistry;
ctx.NN_case_study_ERW_MOD = NN_case_study_ERW_MOD;
if exist('B_out','var');      ctx.B_out = B_out;           end
if exist('O_out','var');      ctx.O_out = O_out;           end
if exist('NavlI_out','var');  ctx.NavlI_out = NavlI_out;   end
if exist('RexmyI_out','var'); ctx.RexmyI_out = RexmyI_out; end
if exist('NupI_L_out','var'); ctx.NupI_L_out = NupI_L_out; end

PAR = cell(N_col,1);
for c = 1:N_col
    ovr = '';
    if isfield(ColCfg(c),'PFT_OVERRIDE') && ~isempty(ColCfg(c).PFT_OVERRIDE)
        ovr = fullfile(Parameter_Directory, [char(ColCfg(c).PFT_OVERRIDE) '.m']);
    end
    PAR{c} = Load_Cell_Parameters(char(ColCfg(c).PARAM_IC), ovr, ctx);
    fprintf('  column %2d  %-14s  %s\n', c, topo.cell_name{c}, ...
        ternary(isempty(ovr),'base parameters','PFT override'));
end

%% Soil layering must agree across columns: the lateral solver exchanges
%% water layer by layer, which is only defined on a shared discretisation.
for c = 2:N_col
    if ~isequal(PAR{c}.Zs, PAR{1}.Zs)
        error('MAIN_FRAME_BG_ERW_2D:LayerMismatch', ...
            ['Column %d uses a different soil layer mesh from column 1. ' ...
             'Lateral exchange is computed layer by layer and is undefined ' ...
             'across differing meshes.'], c);
    end
end
if any(abs(cellfun(@(p) p.Zbio,PAR)-PAR{1}.Zbio) > 1e-9)
    error('MAIN_FRAME_BG_ERW_2D:ActiveDepthMismatch', ...
        'All columns must use the same active biogeochemical depth.');
end
Zlayer_m = diff(PAR{1}.Zs)/1000;      %% [m] layer thickness, shared

%% This retained sub-field topology has no elevation or roughness closure.
%% A sloped case requires T&C's raster routing machinery and cannot be
%% represented by the flat adjacency graph without inventing flow paths.
if N_col > 1
    for c = 1:N_col
        if abs(PAR{c}.Slo_top) > 1e-12 || any(abs(PAR{c}.Slo_pot) > 1e-12)
            error('MAIN_FRAME_BG_ERW_2D:SlopedDomainUnsupported', ...
                ['The sub-field graph requires zero Slo_top and Slo_pot. ' ...
                 'Use raster routing when a DTM and surface-flow parameters are available.']);
        end
    end
end

bio_active_fraction = zeros(N_col,ms);
for c = 1:N_col
    active_thickness_m = PAR{c}.Bio_Zs*(PAR{c}.Zbio/1000);
    bio_active_fraction(c,:) = active_thickness_m./Zlayer_m;
end
if any(~isfinite(bio_active_fraction(:))) || ...
        any(bio_active_fraction(:) < -1e-12) || ...
        any(bio_active_fraction(:) > 1+1e-12)
    error('MAIN_FRAME_BG_ERW_2D:BadActiveLayer', ...
        'Biogeochemical active-layer fractions must lie within each soil layer.');
end
bio_active_fraction = min(1,max(0,bio_active_fraction));

%% Per-column soil parameters for the lateral solver
soilpar_cell = Soil_Param_To_Routing(PAR{1}.Soil_Param, PAR{1}.SPAR);
soilpar_cell = repmat(soilpar_cell, N_col, 1);
for c = 1:N_col
    soilpar_cell(c) = Soil_Param_To_Routing(PAR{c}.Soil_Param, PAR{c}.SPAR);
end
A_cell = topo.cell_area_m2(:);             %% [m2]
A_tot  = sum(A_cell);

%% ========================================================================
%% STATE ALLOCATION   (trailing dimension = column)
%% ========================================================================
z2 = @(a,b) zeros(a,b,N_col);
z3 = @(a,b,d) zeros(a,b,d,N_col);

V=z2(NN,ms); O=z2(NN,ms); Qi_out=z2(NN,ms); WTR=z2(NN,ms); POT=z2(NN,ms);
Tdp=z2(NN,ms); Vice=z2(NN,ms); Oice=z2(NN,ms);
OF=z2(NN,1); OS=z2(NN,1); ZWT=z2(NN,1);
Pr_sno=z2(NN,1); Pr_liq=z2(NN,1); er=z2(NN,1);
Ts=z2(NN,1); ra=z2(NN,1); r_soil=z2(NN,1); b_soil=z2(NN,1); Tice=z2(NN,1);
Rn=z2(NN,1); H=z2(NN,1); QE=z2(NN,1); Qv=z2(NN,1); Lpho=z2(NN,1);
EG=z2(NN,1); ESN=z2(NN,1); ESN_In=z2(NN,1); EWAT=z2(NN,1);
EIn_urb=z2(NN,1); EIn_rock=z2(NN,1); dw_SNO=z2(NN,1);
G=z2(NN,1); SWE=z2(NN,1); SND=z2(NN,1); ros=z2(NN,1);
In_SWE=z2(NN,1); SP_wc=z2(NN,1); WR_SP=z2(NN,1); U_SWE=z2(NN,1); NIn_SWE=z2(NN,1);
dQ=z2(NN,1); Qfm=z2(NN,1); t_sls=z2(NN,1); DQ=z2(NN,1); DT=z2(NN,1);
In_urb=z2(NN,1); In_rock=z2(NN,1); SE_rock=z2(NN,1); SE_urb=z2(NN,1);
Csno=z2(NN,1); WIS=z2(NN,1); Lk=z2(NN,1); f=z2(NN,1); Rh=z2(NN,1); Rd=z2(NN,1);
NDVI=z2(NN,1); alp_soil=z2(NN,1); tau_sno=z2(NN,1); e_sno=z2(NN,1); ALB=z2(NN,1);
EK=z2(NN,1); dQVEG=z2(NN,1); TsV=z2(NN,1); HV=z2(NN,1); QEV=z2(NN,1);
Cice=z2(NN,1); Lk_wat=z2(NN,1); Lk_rock=z2(NN,1); EICE=z2(NN,1);
WAT=z2(NN,1); ICE=z2(NN,1); ICE_D=z2(NN,1); WR_IP=z2(NN,1); NIce=z2(NN,1);
Cicew=z2(NN,1); IP_wc=z2(NN,1); Csnow=z2(NN,1); FROCK=z2(NN,1);
Imelt=z2(NN,1); Smelt=z2(NN,1); Tdamp=z2(NN,1); Gfin=z2(NN,1);
In_Litter=z2(NN,1); ELitter=z2(NN,1); Ws_under=z2(NN,1);
r_litter=z2(NN,cc);
In_H=z2(NN,cc); In_L=z2(NN,cc); OH=z2(NN,cc); OL=z2(NN,cc);
Psi_s_H=z2(NN,cc); Psi_s_L=z2(NN,cc); Tdp_H=z2(NN,cc); Tdp_L=z2(NN,cc);
rap_H=z2(NN,cc); rap_L=z2(NN,cc); Dr_H=z2(NN,cc); Dr_L=z2(NN,cc);
T_H=z2(NN,cc); T_L=z2(NN,cc); EIn_H=z2(NN,cc); EIn_L=z2(NN,cc);
rb_H=z2(NN,cc); rb_L=z2(NN,cc);
An_H=z2(NN,cc); Rdark_H=z2(NN,cc); An_L=z2(NN,cc); Rdark_L=z2(NN,cc);
Ci_sunH=z2(NN,cc); Ci_sunL=z2(NN,cc); Ci_shdH=z2(NN,cc); Ci_shdL=z2(NN,cc);
rs_sunH=z2(NN,cc); rs_sunL=z2(NN,cc); rs_shdH=z2(NN,cc); rs_shdL=z2(NN,cc);
Vx_H=z2(NN,cc); Vl_H=z2(NN,cc); Vx_L=z2(NN,cc); Vl_L=z2(NN,cc);
Psi_x_H=z2(NN,cc); Psi_l_H=z2(NN,cc); Psi_x_L=z2(NN,cc); Psi_l_L=z2(NN,cc);
gsr_H=z2(NN,cc); Jsx_H=z2(NN,cc); Jxl_H=z2(NN,cc); Kleaf_H=z2(NN,cc); Kx_H=z2(NN,cc);
gsr_L=z2(NN,cc); Jsx_L=z2(NN,cc); Jxl_L=z2(NN,cc); Kleaf_L=z2(NN,cc); Kx_L=z2(NN,cc);
fapar_H=z2(NN,cc); fapar_L=z2(NN,cc); SIF_H=z2(NN,cc); SIF_L=z2(NN,cc);
Tdeb = zeros(NN, max(1,length(PAR{1}.Zs_deb)-1), N_col);

%% Daily vegetation and biogeochemistry state
LAI_L=z2(NNd,cc); B_L=z3(NNd,cc,NCP); NPP_L=z2(NNd,cc); Rg_L=z2(NNd,cc);
RA_L=z2(NNd,cc); Rmc_L=z2(NNd,cc); ANPP_L=z2(NNd,cc);
Rms_L=z2(NNd,cc); Rmr_L=z2(NNd,cc); PHE_S_L=z2(NNd,cc); dflo_L=z2(NNd,cc);
AgeL_L=z2(NNd,cc); e_rel_L=ones(NNd,cc,N_col); SAI_L=z2(NNd,cc); hc_L=z2(NNd,cc);
e_relN_L=ones(NNd,cc,N_col); BA_L=z2(NNd,cc);
LAI_H=z2(NNd,cc); B_H=z3(NNd,cc,NCP); NPP_H=z2(NNd,cc); Rg_H=z2(NNd,cc);
RA_H=z2(NNd,cc); Rms_H=z2(NNd,cc); Rmr_H=z2(NNd,cc); ANPP_H=z2(NNd,cc);
PHE_S_H=z2(NNd,cc); dflo_H=z2(NNd,cc); Rmc_H=z2(NNd,cc);
AgeL_H=z2(NNd,cc); e_rel_H=ones(NNd,cc,N_col); SAI_H=z2(NNd,cc); hc_H=z2(NNd,cc);
e_relN_H=ones(NNd,cc,N_col); BA_H=z2(NNd,cc);
Rrootl_H=z2(NNd,cc); Rrootl_L=z2(NNd,cc);
Bfac_dayH=ones(NNd,cc,N_col); Bfac_weekH=ones(NNd,cc,N_col);
NPPI_H=z2(NNd,cc); TdpI_H=z2(NNd,cc);
Bfac_dayL=ones(NNd,cc,N_col); Bfac_weekL=ones(NNd,cc,N_col);
NPPI_L=z2(NNd,cc); TdpI_L=z2(NNd,cc);
Sr_H=z2(NNd,cc); Slf_H=z2(NNd,cc); Sfr_H=z2(NNd,cc); Swm_H=z2(NNd,cc); Sll_H=z2(NNd,cc);
Sr_L=z2(NNd,cc); Slf_L=z2(NNd,cc); Sfr_L=z2(NNd,cc); Swm_L=z2(NNd,cc); Sll_L=z2(NNd,cc);
LAIdead_H=z2(NNd,cc); LAIdead_L=z2(NNd,cc);
Rexmy_H=z3(NNd,cc,3); AgeDL_H=z2(NNd,cc);
Nreserve_H=z2(NNd,cc); Preserve_H=z2(NNd,cc); Kreserve_H=z2(NNd,cc);
Careserve_H=z2(NNd,cc); Mgreserve_H=z2(NNd,cc); Sireserve_H=z2(NNd,cc);
rNc_H=z2(NNd,cc); rPc_H=z2(NNd,cc); rKc_H=z2(NNd,cc);
rCac_H=z2(NNd,cc); rMgc_H=z2(NNd,cc); rSic_H=z2(NNd,cc);
Rexmy_L=z3(NNd,cc,3); AgeDL_L=z2(NNd,cc);
Nreserve_L=z2(NNd,cc); Preserve_L=z2(NNd,cc); Kreserve_L=z2(NNd,cc);
Careserve_L=z2(NNd,cc); Mgreserve_L=z2(NNd,cc); Sireserve_L=z2(NNd,cc);
rNc_L=z2(NNd,cc); rPc_L=z2(NNd,cc); rKc_L=z2(NNd,cc);
rCac_L=z2(NNd,cc); rMgc_L=z2(NNd,cc); rSic_L=z2(NNd,cc);
NBLeaf_H=z2(NNd,cc); PARI_H=z3(NNd,cc,3); NBLI_H=z2(NNd,cc);
NBLeaf_L=z2(NNd,cc); PARI_L=z3(NNd,cc,3); NBLI_L=z2(NNd,cc);
NupI_H=z3(NNd,cc,6); NupI_L=z3(NNd,cc,6); NavlI=z2(NNd,6);
RB_L=z3(NNd,cc,7); RB_H=z3(NNd,cc,7);
NuLit_H=z3(NNd,cc,6); NuLit_L=z3(NNd,cc,6);
BLit=z2(NNd,cc);
Nuptake_H=z2(NNd,cc); NH4uptake_H=z2(NNd,cc); NO3uptake_H=z2(NNd,cc);
Puptake_H=z2(NNd,cc); Kuptake_H=z2(NNd,cc);
Cauptake_H=z2(NNd,cc); Mguptake_H=z2(NNd,cc); Siuptake_H=z2(NNd,cc);
FNC_H=z2(NNd,cc);
TexC_H=z2(NNd,cc); TexN_H=z2(NNd,cc); TexP_H=z2(NNd,cc); TexK_H=z2(NNd,cc);
TNIT_H=z2(NNd,cc); TPHO_H=z2(NNd,cc); TPOT_H=z2(NNd,cc);
SupN_H=z2(NNd,cc); SupP_H=z2(NNd,cc); SupK_H=z2(NNd,cc);
TexCa_H=z2(NNd,cc); TexMg_H=z2(NNd,cc); TexSi_H=z2(NNd,cc);
TCAL_H=z2(NNd,cc); TMAG_H=z2(NNd,cc); TSIL_H=z2(NNd,cc);
SupCa_H=z2(NNd,cc); SupMg_H=z2(NNd,cc); SupSi_H=z2(NNd,cc);
Nuptake_L=z2(NNd,cc); NH4uptake_L=z2(NNd,cc); NO3uptake_L=z2(NNd,cc);
Puptake_L=z2(NNd,cc); Kuptake_L=z2(NNd,cc);
Cauptake_L=z2(NNd,cc); Mguptake_L=z2(NNd,cc); Siuptake_L=z2(NNd,cc);
FNC_L=z2(NNd,cc);
TexC_L=z2(NNd,cc); TexN_L=z2(NNd,cc); TexP_L=z2(NNd,cc); TexK_L=z2(NNd,cc);
TNIT_L=z2(NNd,cc); TPHO_L=z2(NNd,cc); TPOT_L=z2(NNd,cc);
SupN_L=z2(NNd,cc); SupP_L=z2(NNd,cc); SupK_L=z2(NNd,cc);
TexCa_L=z2(NNd,cc); TexMg_L=z2(NNd,cc); TexSi_L=z2(NNd,cc);
TCAL_L=z2(NNd,cc); TMAG_L=z2(NNd,cc); TSIL_L=z2(NNd,cc);
SupCa_L=z2(NNd,cc); SupMg_L=z2(NNd,cc); SupSi_L=z2(NNd,cc);
ISOIL_H=z3(NNd,cc,27); ISOIL_L=z3(NNd,cc,27);
P=z2(NNd,79);
R_litter=z2(NNd,1); R_microbe=z2(NNd,1); R_litter_sur=z2(NNd,1); R_ew=z2(NNd,1);
VOL=z2(NNd,1); N2flx=z2(NNd,1); N2Oflx=z2(NNd,1); PH_d=z2(NNd,1);
Min_N=z2(NNd,1); Min_P=z2(NNd,1); R_bacteria=z2(NNd,1);
RmycAM=z2(NNd,1); RmycEM=z2(NNd,1); Prod_B=z2(NNd,1); Prod_F=z2(NNd,1);
BfixN=z2(NNd,1); LitFirEmi=z2(NNd,2);
LEAK_NH4=z2(NNd,1); LEAK_NO3=z2(NNd,1); LEAK_P=z2(NNd,1); LEAK_K=z2(NNd,1);
LEAK_DOC=z2(NNd,1); LEAK_DON=z2(NNd,1); LEAK_DOP=z2(NNd,1);
LEAK_Ca=z2(NNd,1); LEAK_Mg=z2(NNd,1); LEAK_Si=z2(NNd,1);
RexmyI=z2(NNd,3);
AgrHarNut=z2(NNd,6);
Ccrown_t=ones(NNd,cc,N_col);
AgePl_H=z2(NNd,cc); AgePl_L=z2(NNd,cc); Tden_H=z2(NNd,cc); Tden_L=z2(NNd,cc);
TBio_Ht=z2(NNd,cc); TBio_Lt=z2(NNd,cc);
q_runon=z2(NN,1); Qi_in=z2(NN,ms);
ManIH=zeros(cc,N_col); ManIL=zeros(cc,N_col);

%% Per-column scalars carried between steps
Ccrown_c = zeros(N_col,cc);  Cbare_c = zeros(N_col,1);
ZR95_H_c = zeros(N_col,cc);  ZR95_L_c = zeros(N_col,cc);
RfH_Zs_c = cell(N_col,1);    RfL_Zs_c = cell(N_col,1);
PH_c     = zeros(N_col,1);   Tstm0_c  = zeros(N_col,1);
Bam_c    = zeros(N_col,1);   Bem_c    = zeros(N_col,1);
snow_alb_c = cell(N_col,1);
B_IO_c   = cell(N_col,1);
TBio_H_c = zeros(N_col,cc);  TBio_L_c = zeros(N_col,cc);

for c = 1:N_col
    p = PAR{c};
    Ccrown_c(c,:) = p.Ccrown;   Cbare_c(c) = p.Cbare;
    ZR95_H_c(c,:) = p.ZR95_H;   ZR95_L_c(c,:) = p.ZR95_L;
    RfH_Zs_c{c} = p.RfH_Zs;     RfL_Zs_c{c} = p.RfL_Zs;
    PH_c(c) = p.PH(1);          snow_alb_c{c} = p.snow_alb;
    B_IO_c{c} = p.B_IO;
    TBio_H_c(c,:) = p.TBio_H;   TBio_L_c(c,:) = p.TBio_L;

    %% Seed row 1 of every state array from this column's parameter file
    O(1,:,c) = p.O(1,:);  V(1,:,c) = p.V(1,:);
    Ts(1,c) = p.Ts(1);  Tstm0_c(c) = p.Ts(1);
    Tdp(1,:,c) = p.Tdp(1,:);  Tdamp(1,c) = p.Tdamp(1);
    Tdp_H(1,:,c) = p.Tdp_H(1,:);  Tdp_L(1,:,c) = p.Tdp_L(1,:);
    TdpI_H(1,:,c) = p.TdpI_H(1,:); TdpI_L(1,:,c) = p.TdpI_L(1,:);
    SWE(1,c)=p.SWE(1); SND(1,c)=p.SND(1); ros(1,c)=p.ros(1);
    In_SWE(1,c)=p.In_SWE(1); SP_wc(1,c)=p.SP_wc(1); t_sls(1,c)=p.t_sls(1);
    e_sno(1,c)=p.e_sno(1); tau_sno(1,c)=p.tau_sno(1); EK(1,c)=p.EK(1);
    WAT(1,c)=p.WAT(1); ICE(1,c)=p.ICE(1); IP_wc(1,c)=p.IP_wc(1);
    ICE_D(1,c)=p.ICE_D(1); FROCK(1,c)=p.FROCK(1); Ws_under(1,c)=p.Ws_under(1);
    In_L(1,:,c)=p.In_L(1,:); In_H(1,:,c)=p.In_H(1,:);
    In_urb(1,c)=p.In_urb(1); In_rock(1,c)=p.In_rock(1); In_Litter(1,c)=p.In_Litter(1);
    Ci_sunL(1,:,c)=p.Ci_sunL(1,:); Ci_sunH(1,:,c)=p.Ci_sunH(1,:);
    Ci_shdL(1,:,c)=p.Ci_shdL(1,:); Ci_shdH(1,:,c)=p.Ci_shdH(1,:);
    LAI_H(1,:,c)=p.LAI_H(1,:); LAI_L(1,:,c)=p.LAI_L(1,:);
    B_H(1,:,:,c)=p.B_H(1,:,:);  B_L(1,:,:,c)=p.B_L(1,:,:);
    Rrootl_H(1,:,c)=p.Rrootl_H(1,:); Rrootl_L(1,:,c)=p.Rrootl_L(1,:);
    PHE_S_H(1,:,c)=p.PHE_S_H(1,:); PHE_S_L(1,:,c)=p.PHE_S_L(1,:);
    dflo_H(1,:,c)=p.dflo_H(1,:); dflo_L(1,:,c)=p.dflo_L(1,:);
    AgeL_H(1,:,c)=p.AgeL_H(1,:); AgeL_L(1,:,c)=p.AgeL_L(1,:);
    e_rel_H(1,:,c)=p.e_rel_H(1,:); e_rel_L(1,:,c)=p.e_rel_L(1,:);
    hc_H(1,:,c)=p.hc_H(1,:); hc_L(1,:,c)=p.hc_L(1,:);
    SAI_H(1,:,c)=p.SAI_H(1,:); SAI_L(1,:,c)=p.SAI_L(1,:);
    Nreserve_H(1,:,c)=p.Nreserve_H(1,:); Preserve_H(1,:,c)=p.Preserve_H(1,:);
    Kreserve_H(1,:,c)=p.Kreserve_H(1,:);
    Careserve_H(1,:,c)=p.Careserve_H(1,:); Mgreserve_H(1,:,c)=p.Mgreserve_H(1,:);
    Sireserve_H(1,:,c)=p.Sireserve_H(1,:);
    Nreserve_L(1,:,c)=p.Nreserve_L(1,:); Preserve_L(1,:,c)=p.Preserve_L(1,:);
    Kreserve_L(1,:,c)=p.Kreserve_L(1,:);
    Careserve_L(1,:,c)=p.Careserve_L(1,:); Mgreserve_L(1,:,c)=p.Mgreserve_L(1,:);
    Sireserve_L(1,:,c)=p.Sireserve_L(1,:);
    FNC_H(1,:,c)=p.FNC_H(1,:); FNC_L(1,:,c)=p.FNC_L(1,:);
    NupI_H(1,:,:,c)=p.NupI_H(1,:,:); NupI_L(1,:,:,c)=p.NupI_L(1,:,:);
    RexmyI(1,:,c)=p.RexmyI(1,:); NavlI(1,:,c)=p.NavlI(1,:);
    Vx_H(1,:,c)=p.Vx_H; Vl_H(1,:,c)=p.Vl_H;
    Vx_L(1,:,c)=p.Vx_L; Vl_L(1,:,c)=p.Vl_L;
    BLit(1,:,c)=p.BLit(1);
    if OPT_SoilBiogeochemistry == 1
        P(1,:,c) = p.P(1,:);
        PH_d(1,c) = p.PH(1);
    end
    Check_Land_Cover_Fractions(p.Crock,p.Curb,p.Cwat,p.Cbare,p.Ccrown);
end
CcrownFIX = Ccrown_c;
[~,~,lateral_species,lateral_pool_idx,lateral_fraction_fields] = Biogeo_Lateral_Pools( ...
    reshape(P(1,:,:),79,N_col).');
N_lateral_species = numel(lateral_species);
i_lateral_NO3 = find(strcmp(lateral_species,'NO3'),1);
lateral_fraction = zeros(N_col,N_lateral_species);
BiogeoPar_lateral = Biogeochemistry_Parameter(1);
for s = 1:N_lateral_species
    lateral_fraction(:,s) = BiogeoPar_lateral.(lateral_fraction_fields{s});
end
if any(~isfinite(lateral_fraction(:))) || any(lateral_fraction(:) < 0) || ...
        any(lateral_fraction(:) > 1)
    error('MAIN_FRAME_BG_ERW_2D:BadMobileFraction', ...
        'Biogeo_Leakage mobility fractions must be finite and within [0,1].');
end

%% ========================================================================
%% NUMERICAL OPTIONS  (identical to the single-column frame)
%% ========================================================================
Opt_CR = optimset('TolFun',1);
Opt_ST = optimset('TolFun',0.1);
OPT_SM  = odeset('AbsTol',0.05,'MaxStep',dth);
OPT_VD  = odeset('AbsTol',0.05);
OPT_STh = odeset('AbsTol',0.02);
OPT_PH  = odeset('AbsTol',0.01);
OPT_FR_SOIL = 1; OPT_VegSnow = 1; OPT_SoilTemp = 1;
if ~exist('OPT_VCA','var');  OPT_VCA = 0; OPT_ALLOME = 0; end
if ~exist('OPT_PlantHydr','var');      OPT_PlantHydr = 0;      end
if ~exist('OPT_EnvLimitGrowth','var'); OPT_EnvLimitGrowth = 1; end
OPT_WET = 0;

%% ERW per-column state
%% Cell array, not a struct array: each column's ERW state gains fields
%% after its first call, which a struct array cannot absorb elementwise.
ERW = cell(N_col,1);
for c = 1:N_col; ERW{c} = struct('init',false,'has_update',false); end
conv_mol = 1e6; conv_Al = 1e3;
Tdp_bio_c = zeros(N_col,1); VSUM_c = zeros(N_col,1); VTSUM_c = zeros(N_col,1);
lateral_diag = struct('water_buf_to_main', zeros(NNd,1), ...
                      'mass_buf_to_main',  zeros(NNd,N_lateral_species), ...
                      'NO3_buf_to_main',   zeros(NNd,1), ...
                      'heat_buf_to_main_J',zeros(NNd,1), ...
                      'heat_residual_J',   zeros(NNd,1), ...
                      'heat_residual_relative',zeros(NNd,1), ...
                      'heat_sub_steps_max',zeros(NNd,1), ...
                      'erw_mass_buf_to_main_umol',zeros(NNd,4), ...
                      'erw_CO2_buf_to_main_umol',zeros(NNd,1), ...
                      'erw_mass_residual_umol',zeros(NNd,4), ...
                      'sub_steps_max',     zeros(NNd,1), ...
                      'species_names',     {lateral_species}, ...
                      'erw_species_names', {{'Na','Al','An','DIC'}});

%% Domain-boundary accounting. The subsurface lateral perimeter is closed;
%% bottom leakage and unrouted surface runoff leave the modeled plot.
exit_water_vert_m3 = zeros(N_col,1);
exit_water_surface_m3 = zeros(N_col,1);

j = 1;
tic;

%% ========================================================================
%% TIME LOOP
%% ========================================================================
for i = 2:NN
    if (mod(i,1000) == 0) || (i == 2)
        fprintf('Iter %d / %d   (%.1f%%)\n', i, NN, 100*i/NN);
        %% Pool progress, printed to the run log only. Read-only: it
        %% touches no state and cannot alter the simulation. Column 1 is
        %% enough to spot a pool heading to zero or running away.
        if OPT_SoilBiogeochemistry == 1 && j > 1
            fprintf(['  POOLS  SOMC %.4g  POC %.4g  MOC %.4g  DOC %.4g  ' ...
                     'Bact %.4g  Fung %.4g  AM %.4g  SOMN %.4g  ' ...
                     'NH4 %.4g  NO3 %.4g  Pmin %.4g  Ca %.4g\n'], ...
                P(j-1,9,1), P(j-1,10,1), P(j-1,11,1), ...
                P(j-1,12,1)+P(j-1,13,1), P(j-1,18,1), P(j-1,19,1), ...
                P(j-1,20,1), P(j-1,26,1), P(j-1,31,1), P(j-1,32,1), ...
                P(j-1,43,1), P(j-1,60,1));
        end
    end
    pdind = max(1,i-24):i-1;
    h_S = NaN;

    %% ------------------------------------------------------------------
    %% DAILY BLOCK
    %% ------------------------------------------------------------------
    if (Datam(i,4) == 1)
        j = j+1;
        [jDay_i] = JULIAN_DAY(Datam(i,:));
        [~,~,~,~,~,L_day_j] = SetSunVariables(Datam(i,:),DeltaGMT,Lon,Lat,t_bef,t_aft);
        jDay(j) = jDay_i; L_day(j) = L_day_j; %#ok<SAGROW>

        for c = 1:N_col
            p = PAR{c};
            Ccrown = Ccrown_c(c,:); Cbare = Cbare_c(c);
            ZR95_H = ZR95_H_c(c,:); ZR95_L = ZR95_L_c(c,:);
            RfH_Zs = RfH_Zs_c{c};   RfL_Zs = RfL_Zs_c{c};
            B_IO_c_now = B_IO_c{c};

            %% Rotation fertiliser rule for cropped cells
            if isempty(char(ColCfg(c).PFT_OVERRIDE))
                if ismember(Datam(i,1),[2009 2011 2012 2014 2015 2017 2018 2020])
                    B_IO_c_now.FertN(81:100) = 0.84;
                else
                    B_IO_c_now.FertN = zeros(1,366);
                end
                B_IO_c{c} = B_IO_c_now;
            end

            [Se_bio,Se_fc,Psi_bio,Tdp_bio,VSUM,VTSUM] = Biogeo_environment( ...
                Tdp(pdind,:,c),O(pdind,:,c),V(pdind,:,c),p.Soil_Param,p.Phy,p.SPAR,p.Bio_Zs);
            Tdp_bio_c(c) = Tdp_bio; VSUM_c(c) = VSUM; VTSUM_c(c) = VTSUM;

            if OPT_SoilBiogeochemistry == 1
                IS    = Ccrown*squeeze(ISOIL_L(j-1,:,:,c)) + Ccrown*squeeze(ISOIL_H(j-1,:,:,c));
                REXMY = Ccrown*squeeze(Rexmy_L(j-1,:,:,c)) + Ccrown*squeeze(Rexmy_H(j-1,:,:,c));
                FireA = 1*((sum(ManIH(:,c)==-5) + sum(ManIL(:,c)==-5)) > 0);
                if sum(Ccrown) > 0
                    OPT_BfixN = p.ParEx_L(Ccrown>0).bfix;
                else
                    OPT_BfixN = 0;
                end
                PH_d(j,c) = PH_c(c);

                [P(j,:,c),LEAK_NH4(j,c),LEAK_NO3(j,c),LEAK_P(j,c),LEAK_K(j,c), ...
                 LEAK_DOC(j,c),LEAK_DON(j,c),LEAK_DOP(j,c), ...
                 LEAK_Ca(j,c),LEAK_Mg(j,c),LEAK_Si(j,c), ...
                 Nuptake_H(j,:,c),NH4uptake_H(j,:,c),NO3uptake_H(j,:,c), ...
                 Puptake_H(j,:,c),Kuptake_H(j,:,c), ...
                 Nuptake_L(j,:,c),NH4uptake_L(j,:,c),NO3uptake_L(j,:,c), ...
                 Puptake_L(j,:,c),Kuptake_L(j,:,c), ...
                 Cauptake_H(j,:,c),Mguptake_H(j,:,c),Siuptake_H(j,:,c), ...
                 Cauptake_L(j,:,c),Mguptake_L(j,:,c),Siuptake_L(j,:,c), ...
                 RexmyI(j,:,c), ...
                 R_litter(j,c),R_microbe(j,c),R_litter_sur(j,c),R_ew(j,c), ...
                 VOL(j,c),N2flx(j,c),Min_N(j,c),Min_P(j,c), ...
                 R_bacteria(j,c),RmycAM(j,c),RmycEM(j,c),Prod_B(j,c),Prod_F(j,c), ...
                 BfixN(j,c),NavlI(j,:,c),LitFirEmi(j,:,c),N2Oflx(j,c)] = ...
                 BIOGEO_UNIT(P(j-1,:,c),IS,p.Zbio,sum(p.Bio_Zs.*p.rsd), ...
                    PH_c(c),Tdp_bio,mean(Ta(pdind)),Psi_bio,Se_bio,Se_fc,VSUM,VTSUM, ...
                    Ccrown,p.Bio_Zs,RfH_Zs,RfL_Zs, ...
                    sum(Lk(pdind,1,c)),sum(T_H(pdind,:,c),1),sum(T_L(pdind,:,c),1), ...
                    B_H(j-1,:,3,c),B_L(j-1,:,3,c),LAI_H(j-1,:,c),LAI_L(j-1,:,c), ...
                    SupN_H(j-1,:,c),SupP_H(j-1,:,c),SupK_H(j-1,:,c), ...
                    SupN_L(j-1,:,c),SupP_L(j-1,:,c),SupK_L(j-1,:,c), ...
                    SupCa_H(j-1,:,c),SupMg_H(j-1,:,c),SupSi_H(j-1,:,c), ...
                    SupCa_L(j-1,:,c),SupMg_L(j-1,:,c),SupSi_L(j-1,:,c), ...
                    REXMY,RexmyI(j-1,:,c),p.ExEM,NavlI(j-1,:,c),p.Pcla,p.Psan, ...
                    B_IO_c_now,Datam(i,:),FireA,0, ...
                    OPT_ERW,OPT_BfixN);

                BLit(j,:,c) = 0.002*sum(P(j,1:5,c))*Ccrown';
                Bam_c(c) = P(j,20,c);
                Bem_c(c) = P(j,21,c);
            else
                BLit(j,:,c) = 0;
                NavlI(j,:,c) = 0; P(j,:,c) = 0;
                Bam_c(c) = 0; Bem_c(c) = 0;
            end

            %% ---- Vegetation, per PFT slot -----------------------------
            for k = 1:length(Ccrown)
                if ZR95_H(k) > 0
                    [LAI_H(j,k,c),B_H(j,k,:,c),NPP_H(j,k,c),ANPP_H(j,k,c),Rg_H(j,k,c), ...
                     RA_H(j,k,c),Rms_H(j,k,c),Rmr_H(j,k,c),Rmc_H(j,k,c),PHE_S_H(j,k,c), ...
                     dflo_H(j,k,c),AgeL_H(j,k,c),e_rel_H(j,k,c),e_relN_H(j,k,c), ...
                     LAIdead_H(j,k,c),NBLeaf_H(j,k,c), ...
                     Sr_H(j,k,c),Slf_H(j,k,c),Sfr_H(j,k,c),Swm_H(j,k,c),Sll_H(j,k,c), ...
                     Rexmy_H(j,k,:,c),Rrootl_H(j,k,c), ...
                     AgeDL_H(j,k,c),Bfac_dayH(j,k,c),Bfac_weekH(j,k,c),NPPI_H(j,k,c), ...
                     TdpI_H(j,k,c),NupI_H(j,k,:,c),PARI_H(j,k,:,c),NBLI_H(j,k,c), ...
                     RB_H(j,k,:,c),FNC_H(j,k,c), ...
                     Nreserve_H(j,k,c),Preserve_H(j,k,c),Kreserve_H(j,k,c), ...
                     Careserve_H(j,k,c),Mgreserve_H(j,k,c),Sireserve_H(j,k,c), ...
                     rNc_H(j,k,c),rPc_H(j,k,c),rKc_H(j,k,c), ...
                     rCac_H(j,k,c),rMgc_H(j,k,c),rSic_H(j,k,c),ManIH(k,c)] = ...
                     VEGGIE_UNIT( ...
                        B_H(j-1,k,:,c),PHE_S_H(j-1,k,c),dflo_H(j-1,k,c), ...
                        AgeL_H(j-1,k,c),AgeDL_H(j-1,k,c), ...
                        Ta(pdind),Tdp_H(pdind,k,c),PARB(pdind)+PARD(pdind), ...
                        Psi_x_H(pdind,k,c),Psi_l_H(pdind,k,c), ...
                        An_H(pdind,k,c),Rdark_H(pdind,k,c),NPP_H(j-1,k,c),jDay(j),Datam(i,:), ...
                        NPPI_H(j-1,k,c),TdpI_H(j-1,k,c),Bfac_weekH(j-1,k,c), ...
                        NupI_H(j-1,k,:,c),NavlI(j,:,c),PARI_H(j-1,k,:,c), ...
                        NBLI_H(j-1,k,c),NBLeaf_H(j-1,k,c), ...
                        L_day(j),p.Lmax_day,p.VegH_Param_Dyn,k, ...
                        Nreserve_H(j-1,k,c),Preserve_H(j-1,k,c),Kreserve_H(j-1,k,c), ...
                        Careserve_H(j-1,k,c),Mgreserve_H(j-1,k,c),Sireserve_H(j-1,k,c), ...
                        Nuptake_H(j,k,c),Puptake_H(j,k,c),Kuptake_H(j,k,c), ...
                        Cauptake_H(j,k,c),Mguptake_H(j,k,c),Siuptake_H(j,k,c), ...
                        FNC_H(j-1,k,c),Se_bio,Tdp_bio, ...
                        p.ParEx_H(k),p.ExEM,Bam_c(c),Bem_c(c),p.Mpar_H(k), ...
                        TBio_Ht(j-1,k,c),OPT_EnvLimitGrowth,OPT_VCA,OPT_VD,OPT_SoilBiogeochemistry);

                    [TexC_H(j,k,c),TexN_H(j,k,c),TexP_H(j,k,c),TexK_H(j,k,c), ...
                     TexCa_H(j,k,c),TexMg_H(j,k,c),TexSi_H(j,k,c), ...
                     TNIT_H(j,k,c),TPHO_H(j,k,c),TPOT_H(j,k,c), ...
                     TCAL_H(j,k,c),TMAG_H(j,k,c),TSIL_H(j,k,c), ...
                     NuLit_H(j,k,:,c),Nreserve_H(j,k,c),Preserve_H(j,k,c),Kreserve_H(j,k,c), ...
                     Careserve_H(j,k,c),Mgreserve_H(j,k,c),Sireserve_H(j,k,c), ...
                     SupN_H(j,k,c),SupP_H(j,k,c),SupK_H(j,k,c), ...
                     SupCa_H(j,k,c),SupMg_H(j,k,c),SupSi_H(j,k,c), ...
                     ISOIL_H(j,k,:,c)] = Plant_Exports(B_H(j,k,:,c),B_H(j-1,k,:,c), ...
                        NuLit_H(j-1,k,:,c),Slf_H(j,k,c),Sfr_H(j,k,c),Swm_H(j,k,c), ...
                        Sll_H(j,k,c),Sr_H(j,k,c), ...
                        Rexmy_H(j,k,:,c),p.Stoich_H(k),p.Mpar_H(k),p.fab_H(k),p.fbe_H(k), ...
                        RB_H(j,k,:,c),Nreserve_H(j,k,c),Preserve_H(j,k,c),Kreserve_H(j,k,c), ...
                        Careserve_H(j,k,c),Mgreserve_H(j,k,c),Sireserve_H(j,k,c), ...
                        rNc_H(j,k,c),rPc_H(j,k,c),rKc_H(j,k,c), ...
                        rCac_H(j,k,c),rMgc_H(j,k,c),rSic_H(j,k,c), ...
                        ManIH(k,c),OPT_SoilBiogeochemistry,id_location);

                    SAI_H(j,k,c) = SAI_H(j-1,k,c);
                    if p.aSE_H(k) == 2
                        hc_H(j,k,c) = GrassHeight(LAI_H(j,k,c),LAIdead_H(j,k,c));
                    elseif p.aSE_H(k) == 5
                        [hc_H(j,k,c),SAI_H(j,k,c),B_H(j,:,:,c),Ccrown,ZR95_H, ...
                         Nreserve_H(j-1:j,:,c),Preserve_H(j-1:j,:,c),Kreserve_H(j-1:j,:,c), ...
                         Careserve_H(j-1:j,:,c),Mgreserve_H(j-1:j,:,c),Sireserve_H(j-1:j,:,c), ...
                         AgrHarNut(j,:,c),RfH_Zs] = CropHeightType( ...
                            LAI_H(j,k,c),LAIdead_H(j,k,c),k,ZR95_H,B_H(j,:,:,c), ...
                            p.Zs,p.CASE_ROOT,Ccrown,Cbare+p.Cwat+p.Curb+p.Crock, ...
                            Nreserve_H(j,:,c),Preserve_H(j,:,c),Kreserve_H(j,:,c), ...
                            Careserve_H(j,:,c),Mgreserve_H(j,:,c),Sireserve_H(j,:,c), ...
                            ManIH(:,c),p.Mpar_H,p.VegH_Param_Dyn,OPT_SoilBiogeochemistry);
                    else
                        hc_H(j,k,c) = hc_H(j-1,k,c);
                    end
                else
                    LAI_H(j,k,c)=0; B_H(j,k,:,c)=0; NPP_H(j,k,c)=0; ANPP_H(j,k,c)=0;
                    Rg_H(j,k,c)=0; RA_H(j,k,c)=0; Rms_H(j,k,c)=0; Rmr_H(j,k,c)=0;
                    Rmc_H(j,k,c)=0; PHE_S_H(j,k,c)=0; dflo_H(j,k,c)=0; AgeL_H(j,k,c)=0;
                    e_rel_H(j,k,c)=0; e_relN_H(j,k,c)=0; SAI_H(j,k,c)=0; hc_H(j,k,c)=0;
                    LAIdead_H(j,k,c)=LAIdead_H(j-1,k,c);
                    Sr_H(j,k,c)=0; Slf_H(j,k,c)=0; Sfr_H(j,k,c)=0; Swm_H(j,k,c)=0;
                    Sll_H(j,k,c)=0; Rrootl_H(j,k,c)=0; Rexmy_H(j,k,:,c)=0;
                    AgeDL_H(j,k,c)=0; FNC_H(j,k,c)=1;
                    Nreserve_H(j,k,c)=0; Preserve_H(j,k,c)=0; Kreserve_H(j,k,c)=0;
                    Careserve_H(j,k,c)=0; Mgreserve_H(j,k,c)=0; Sireserve_H(j,k,c)=0;
                    rNc_H(j,k,c)=1; rPc_H(j,k,c)=1; rKc_H(j,k,c)=1;
                    rCac_H(j,k,c)=1; rMgc_H(j,k,c)=1; rSic_H(j,k,c)=1;
                    Bfac_dayH(j,k,c)=0; Bfac_weekH(j,k,c)=0; NPPI_H(j,k,c)=0;
                    TdpI_H(j,k,c)=0; NupI_H(j,k,:,c)=0; RB_H(j,k,:,c)=0;
                    TexC_H(j,k,c)=0; TexN_H(j,k,c)=0; TexP_H(j,k,c)=0; TexK_H(j,k,c)=0;
                    TexCa_H(j,k,c)=0; TexMg_H(j,k,c)=0; TexSi_H(j,k,c)=0;
                    TNIT_H(j,k,c)=0; TPHO_H(j,k,c)=0; TPOT_H(j,k,c)=0;
                    TCAL_H(j,k,c)=0; TMAG_H(j,k,c)=0; TSIL_H(j,k,c)=0;
                    ISOIL_H(j,k,:,c)=0; NuLit_H(j,k,:,c)=0;
                end

                if ZR95_L(k) > 0
                    [LAI_L(j,k,c),B_L(j,k,:,c),NPP_L(j,k,c),ANPP_L(j,k,c),Rg_L(j,k,c), ...
                     RA_L(j,k,c),Rms_L(j,k,c),Rmr_L(j,k,c),Rmc_L(j,k,c),PHE_S_L(j,k,c), ...
                     dflo_L(j,k,c),AgeL_L(j,k,c),e_rel_L(j,k,c),e_relN_L(j,k,c), ...
                     LAIdead_L(j,k,c),NBLeaf_L(j,k,c), ...
                     Sr_L(j,k,c),Slf_L(j,k,c),Sfr_L(j,k,c),Swm_L(j,k,c),Sll_L(j,k,c), ...
                     Rexmy_L(j,k,:,c),Rrootl_L(j,k,c), ...
                     AgeDL_L(j,k,c),Bfac_dayL(j,k,c),Bfac_weekL(j,k,c),NPPI_L(j,k,c), ...
                     TdpI_L(j,k,c),NupI_L(j,k,:,c),PARI_L(j,k,:,c),NBLI_L(j,k,c), ...
                     RB_L(j,k,:,c),FNC_L(j,k,c), ...
                     Nreserve_L(j,k,c),Preserve_L(j,k,c),Kreserve_L(j,k,c), ...
                     Careserve_L(j,k,c),Mgreserve_L(j,k,c),Sireserve_L(j,k,c), ...
                     rNc_L(j,k,c),rPc_L(j,k,c),rKc_L(j,k,c), ...
                     rCac_L(j,k,c),rMgc_L(j,k,c),rSic_L(j,k,c),ManIL(k,c)] = ...
                     VEGGIE_UNIT( ...
                        B_L(j-1,k,:,c),PHE_S_L(j-1,k,c),dflo_L(j-1,k,c), ...
                        AgeL_L(j-1,k,c),AgeDL_L(j-1,k,c), ...
                        Ta(pdind),Tdp_L(pdind,k,c),PARB(pdind)+PARD(pdind), ...
                        Psi_x_L(pdind,k,c),Psi_l_L(pdind,k,c), ...
                        An_L(pdind,k,c),Rdark_L(pdind,k,c),NPP_L(j-1,k,c),jDay(j),Datam(i,:), ...
                        NPPI_L(j-1,k,c),TdpI_L(j-1,k,c),Bfac_weekL(j-1,k,c), ...
                        NupI_L(j-1,k,:,c),NavlI(j,:,c),PARI_L(j-1,k,:,c), ...
                        NBLI_L(j-1,k,c),NBLeaf_L(j-1,k,c), ...
                        L_day(j),p.Lmax_day,p.VegL_Param_Dyn,k, ...
                        Nreserve_L(j-1,k,c),Preserve_L(j-1,k,c),Kreserve_L(j-1,k,c), ...
                        Careserve_L(j-1,k,c),Mgreserve_L(j-1,k,c),Sireserve_L(j-1,k,c), ...
                        Nuptake_L(j,k,c),Puptake_L(j,k,c),Kuptake_L(j,k,c), ...
                        Cauptake_L(j,k,c),Mguptake_L(j,k,c),Siuptake_L(j,k,c), ...
                        FNC_L(j-1,k,c),Se_bio,Tdp_bio, ...
                        p.ParEx_L(k),p.ExEM,Bam_c(c),Bem_c(c),p.Mpar_L(k), ...
                        TBio_Lt(j-1,k,c),OPT_EnvLimitGrowth,OPT_VCA,OPT_VD,OPT_SoilBiogeochemistry);

                    [TexC_L(j,k,c),TexN_L(j,k,c),TexP_L(j,k,c),TexK_L(j,k,c), ...
                     TexCa_L(j,k,c),TexMg_L(j,k,c),TexSi_L(j,k,c), ...
                     TNIT_L(j,k,c),TPHO_L(j,k,c),TPOT_L(j,k,c), ...
                     TCAL_L(j,k,c),TMAG_L(j,k,c),TSIL_L(j,k,c), ...
                     NuLit_L(j,k,:,c),Nreserve_L(j,k,c),Preserve_L(j,k,c),Kreserve_L(j,k,c), ...
                     Careserve_L(j,k,c),Mgreserve_L(j,k,c),Sireserve_L(j,k,c), ...
                     SupN_L(j,k,c),SupP_L(j,k,c),SupK_L(j,k,c), ...
                     SupCa_L(j,k,c),SupMg_L(j,k,c),SupSi_L(j,k,c), ...
                     ISOIL_L(j,k,:,c)] = Plant_Exports(B_L(j,k,:,c),B_L(j-1,k,:,c), ...
                        NuLit_L(j-1,k,:,c),Slf_L(j,k,c),Sfr_L(j,k,c),Swm_L(j,k,c), ...
                        Sll_L(j,k,c),Sr_L(j,k,c), ...
                        Rexmy_L(j,k,:,c),p.Stoich_L(k),p.Mpar_L(k),p.fab_L(k),p.fbe_L(k), ...
                        RB_L(j,k,:,c),Nreserve_L(j,k,c),Preserve_L(j,k,c),Kreserve_L(j,k,c), ...
                        Careserve_L(j,k,c),Mgreserve_L(j,k,c),Sireserve_L(j,k,c), ...
                        rNc_L(j,k,c),rPc_L(j,k,c),rKc_L(j,k,c), ...
                        rCac_L(j,k,c),rMgc_L(j,k,c),rSic_L(j,k,c), ...
                        ManIL(k,c),OPT_SoilBiogeochemistry,id_location);

                    SAI_L(j,k,c) = SAI_L(j-1,k,c);
                    if p.aSE_L(k) == 2
                        hc_L(j,k,c) = GrassHeight(LAI_L(j,k,c),LAIdead_L(j,k,c));
                    elseif p.aSE_L(k) == 5
                        [hc_L(j,k,c),SAI_L(j,k,c),B_L(j,:,:,c),Ccrown,ZR95_L, ...
                         Nreserve_L(j,:,c),Preserve_L(j,:,c),Kreserve_L(j,:,c), ...
                         Careserve_L(j,:,c),Mgreserve_L(j,:,c),Sireserve_L(j,:,c), ...
                         AgrHarNut(j,:,c),RfL_Zs] = CropHeightType( ...
                            LAI_L(j,k,c),LAIdead_L(j,k,c),k,ZR95_L,B_L(j,:,:,c), ...
                            p.Zs,p.CASE_ROOT,Ccrown,Cbare+p.Cwat+p.Curb+p.Crock, ...
                            Nreserve_L(j,:,c),Preserve_L(j,:,c),Kreserve_L(j,:,c), ...
                            Careserve_L(j,:,c),Mgreserve_L(j,:,c),Sireserve_L(j,:,c), ...
                            ManIL(:,c),p.Mpar_L,p.VegL_Param_Dyn,OPT_SoilBiogeochemistry);
                    else
                        hc_L(j,k,c) = hc_L(j-1,k,c);
                    end
                else
                    LAI_L(j,k,c)=0; B_L(j,k,:,c)=0; NPP_L(j,k,c)=0; ANPP_L(j,k,c)=0;
                    Rg_L(j,k,c)=0; RA_L(j,k,c)=0; Rms_L(j,k,c)=0; Rmr_L(j,k,c)=0;
                    PHE_S_L(j,k,c)=0; dflo_L(j,k,c)=0; AgeL_L(j,k,c)=0;
                    e_rel_L(j,k,c)=0; e_relN_L(j,k,c)=0; SAI_L(j,k,c)=0; hc_L(j,k,c)=0;
                    LAIdead_L(j,k,c)=LAIdead_L(j-1,k,c);
                    Sr_L(j,k,c)=0; Slf_L(j,k,c)=0; Sfr_L(j,k,c)=0; Swm_L(j,k,c)=0;
                    Sll_L(j,k,c)=0; Rrootl_L(j,k,c)=0; Rexmy_L(j,k,:,c)=0;
                    AgeDL_L(j,k,c)=0; FNC_L(j,k,c)=1;
                    Nreserve_L(j,k,c)=0; Preserve_L(j,k,c)=0; Kreserve_L(j,k,c)=0;
                    Careserve_L(j,k,c)=0; Mgreserve_L(j,k,c)=0; Sireserve_L(j,k,c)=0;
                    rNc_L(j,k,c)=1; rPc_L(j,k,c)=1; rKc_L(j,k,c)=1;
                    rCac_L(j,k,c)=1; rMgc_L(j,k,c)=1; rSic_L(j,k,c)=1;
                    Bfac_dayL(j,k,c)=0; Bfac_weekL(j,k,c)=0; NPPI_L(j,k,c)=0;
                    TdpI_L(j,k,c)=0; NupI_L(j,k,:,c)=0; RB_L(j,k,:,c)=0;
                    TexC_L(j,k,c)=0; TexN_L(j,k,c)=0; TexP_L(j,k,c)=0; TexK_L(j,k,c)=0;
                    TexCa_L(j,k,c)=0; TexMg_L(j,k,c)=0; TexSi_L(j,k,c)=0;
                    TNIT_L(j,k,c)=0; TPHO_L(j,k,c)=0; TPOT_L(j,k,c)=0;
                    TCAL_L(j,k,c)=0; TMAG_L(j,k,c)=0; TSIL_L(j,k,c)=0;
                    ISOIL_L(j,k,:,c)=0; NuLit_L(j,k,:,c)=0;
                end
            end

            %% Crops may resize Ccrown / rooting depth: persist per column
            Ccrown_c(c,:) = Ccrown;
            ZR95_H_c(c,:) = ZR95_H;  ZR95_L_c(c,:) = ZR95_L;
            RfH_Zs_c{c} = RfH_Zs;    RfL_Zs_c{c} = RfL_Zs;
            if OPT_VCA == 0
                Ccrown_t(j,:,c) = Ccrown;
                Cbare_c(c) = 1 - sum(Ccrown) - p.Cwat - p.Curb - p.Crock;
                Check_Land_Cover_Fractions(p.Crock,p.Curb,p.Cwat,Cbare_c(c),Ccrown);
            end
        end
    end

    %% ------------------------------------------------------------------
    %% HOURLY BLOCK  --  vertical physics, one column at a time
    %% ------------------------------------------------------------------
    for c = 1:N_col
        p = PAR{c};
        Ccrown = Ccrown_c(c,:);
        if N_col > 1
            Qi_in(i,:,c) = 0;
            Slo_top_cell = 0;
            Slo_pot_cell = zeros(1,ms);
        else
            Qi_in(i,:,c) = Qi_out(i-1,:,c);
            Slo_top_cell = p.Slo_top;
            Slo_pot_cell = p.Slo_pot;
        end
        q_runon(i,1,c) = 0;

        [V(i,:,c),Vice(i,:,c),O(i,:,c),Oice(i,:,c),ZWT(i,c),OF(i,c),OS(i,c), ...
         OH(i,:,c),OL(i,:,c),Psi_s_H(i,:,c),Psi_s_L(i,:,c),Rd(i,c),Qi_out(i,:,c),WTR(i,:,c), ...
         Rh(i,c),Lk(i,c),f(i,c),WIS(i,c),Ts(i,c),Pr_sno(i,c),Pr_liq(i,c),Csno(i,c), ...
         Cice(i,c),NDVI(i,c),rb_H(i,:,c),rb_L(i,:,c),rs_sunH(i,:,c), ...
         rs_sunL(i,:,c),rs_shdH(i,:,c),rs_shdL(i,:,c),r_litter(i,:,c), ...
         An_L(i,:,c),An_H(i,:,c),Rdark_L(i,:,c),Rdark_H(i,:,c), ...
         Ci_sunH(i,:,c),Ci_sunL(i,:,c),Ci_shdH(i,:,c),Ci_shdL(i,:,c), ...
         rap_H(i,:,c),rap_L(i,:,c),r_soil(i,c),b_soil(i,c),alp_soil(i,c),ra(i,c),Rn(i,c), ...
         H(i,c),QE(i,c),Qv(i,c),Lpho(i,c),T_H(i,:,c),T_L(i,:,c),EIn_H(i,:,c),EIn_L(i,:,c), ...
         EG(i,c),ESN(i,c),ESN_In(i,c),ELitter(i,c),EWAT(i,c),EICE(i,c),EIn_urb(i,c), ...
         EIn_rock(i,c),dw_SNO(i,c), ...
         G(i,c),Gfin(i,c),Tdp(i,:,c),Tdeb(i,:,c),Tdamp(i,c),Tice(i,c), ...
         Tdp_H(i,:,c),Tdp_L(i,:,c),SWE(i,c),SND(i,c),ros(i,c),In_SWE(i,c),SP_wc(i,c), ...
         WR_SP(i,c),U_SWE(i,c),NIn_SWE(i,c),dQ(i,c),Qfm(i,c),t_sls(i,c),DQ(i,c),DT(i,c), ...
         WAT(i,c),ICE(i,c),ICE_D(i,c),IP_wc(i,c),WR_IP(i,c),NIce(i,c),Cicew(i,c), ...
         Csnow(i,c),FROCK(i,c),Imelt(i,c),Smelt(i,c), ...
         In_H(i,:,c),In_L(i,:,c),In_Litter(i,c),In_urb(i,c),In_rock(i,c), ...
         Dr_H(i,:,c),Dr_L(i,:,c),SE_rock(i,c),SE_urb(i,c),Lk_wat(i,c),Lk_rock(i,c),er(i,c), ...
         gsr_H(i,:,c),Psi_x_H(i,:,c),Psi_l_H(i,:,c),Jsx_H(i,:,c),Jxl_H(i,:,c), ...
         Kleaf_H(i,:,c),Kx_H(i,:,c),Vx_H(i,:,c),Vl_H(i,:,c), ...
         gsr_L(i,:,c),Psi_x_L(i,:,c),Psi_l_L(i,:,c),Jsx_L(i,:,c),Jxl_L(i,:,c), ...
         Kleaf_L(i,:,c),Kx_L(i,:,c),Vx_L(i,:,c),Vl_L(i,:,c), ...
         fapar_H(i,:,c),fapar_L(i,:,c),SIF_H(i,:,c),SIF_L(i,:,c), ...
         snow_alb_new,tau_sno(i,c),e_sno(i,c),Ws_under(i,c),dQVEG(i,c),HV(i,c), ...
         QEV(i,c),TsV(i,c),EK(i,c),POT(i,:,c)] = ...
         HYDROLOGIC_UNIT(V(i-1,:,c),Oice(i-1,:,c),p.aR,p.Zs, ...
            p.EvL_Zs,p.Inf_Zs,p.Zinf,RfH_Zs_c{c},RfL_Zs_c{c},p.dz,p.Dz,ms,p.Kbot, ...
            Pr(i),Ta(i),Ds(i),Ws(i),p.zatm,Ts(i-1,c),IrD(i),dt,dth,ea(i),N(i),Pre(i),Tstm0_c(c), ...
            LAI_H(j,:,c),SAI_H(j,:,c),LAI_L(j,:,c),SAI_L(j,:,c),LAIdead_H(j,:,c),LAIdead_L(j,:,c), ...
            Rrootl_H(j,:,c),Rrootl_L(j,:,c),BLit(j,:,c),p.Sllit,p.Kct, ...
            Datam(i,:),DeltaGMT,Lon,Lat,t_bef,t_aft,h_S, ...
            Ccrown,Cbare_c(c),p.Crock,p.Curb,p.Cwat, ...
            p.Soil_Param,p.Interc_Param,p.SnowIce_Param,p.VegH_Param,p.VegL_Param, ...
            p.Zs_deb,p.Deb_Par, ...
            ZR95_H_c(c,:),ZR95_L_c(c,:), ...
            SAB1(i),SAB2(i),SAD1(i),SAD2(i),PARB(i),PARD(i),p.SvF, ...
            SND(i-1,c),snow_alb_c{c},p.Color_Class,p.OM_H,p.OM_L, ...
            p.PFT_opt_H,p.PFT_opt_L,hc_H(j,:,c),hc_L(j,:,c),p.d_leaf_H,p.d_leaf_L, ...
            Ca(i),Oa,Ci_sunH(i-1,:,c),Ci_shdH(i-1,:,c),Ci_sunL(i-1,:,c),Ci_shdL(i-1,:,c), ...
            e_rel_H(j,:,c),e_relN_H(j,:,c),e_rel_L(j,:,c),e_relN_L(j,:,c), ...
            e_sno(i-1,c),In_H(i-1,:,c),In_L(i-1,:,c),In_Litter(i-1,c),In_urb(i-1,c), ...
            In_rock(i-1,c),SWE(i-1,c),In_SWE(i-1,c), ...
            Tdeb(i-1,:,c),Tdp(i-1,:,c),Tdamp(i-1,c),Tice(i-1,c), ...
            WAT(i-1,c),ICE(i-1,c),IP_wc(i-1,c),ICE_D(i-1,c),Cicew(i-1,c), ...
            Vx_H(i-1,:,c),Vl_H(i-1,:,c),Vx_L(i-1,:,c),Vl_L(i-1,:,c), ...
            Psi_x_H(i-1,:,c),Psi_l_H(i-1,:,c),Psi_x_L(i-1,:,c),Psi_l_L(i-1,:,c), ...
            FROCK(i-1,c),p.Krock,Ws_under(i-1,c), ...
            Tdew(i),t_sls(i-1,c),ros(i-1,c),SP_wc(i-1,c),p.fpr,Pr_sno(pdind,1,c), ...
            p.In_max_urb,p.In_max_rock,p.K_usle,tau_sno(i-1,c),Ta(pdind), ...
            Slo_top_cell,Slo_pot_cell,p.Asur,p.Ared,p.aTop,EK(i-1,c),q_runon(i,1,c),Qi_in(i,:,c), ...
            p.pow_dis,p.a_dis, ...
            p.SPAR,p.SN,OPT_VegSnow,OPT_SoilTemp,OPT_PlantHydr,Opt_CR,Opt_ST, ...
            OPT_SM,OPT_STh,OPT_FR_SOIL,OPT_PH);

        snow_alb_c{c} = snow_alb_new;
        Tstm0_c(c) = 2*Ts(i,c) - Ts(i-1,c);
        STOT = SAB1(i)+SAB2(i)+SAD1(i)+SAD2(i);
        ALB(i,c) = SAB1(i)/STOT*snow_alb_new.dir_vis + SAD1(i)/STOT*snow_alb_new.dif_vis + ...
                   SAB2(i)/STOT*snow_alb_new.dir_nir + SAD2(i)/STOT*snow_alb_new.dif_nir;

        %% Vertical exit at the base of the column [m3]
        exit_water_vert_m3(c) = exit_water_vert_m3(c) + (Lk(i,c)/1000)*A_cell(c)*dth;
        exit_water_surface_m3(c) = exit_water_surface_m3(c) + ...
            ((Rh(i,c)+Rd(i,c))/1000)*A_cell(c);
    end

    %% ------------------------------------------------------------------
    %% LATERAL COUPLING  --  operator-split, after every column has
    %% taken its vertical step.
    %% ------------------------------------------------------------------
    if N_col > 1
        theta_cell = reshape(O(i,:,:), ms, N_col).';        %% [N_col x ms]
        Oice_cell = reshape(Oice(i,:,:), ms, N_col).';
        T_cell = reshape(Tdp(i,:,:), ms, N_col).';

        %% Biogeo_Leakage defines which T&C pools are dissolved and mobile.
        %% The pools are vertically integrated [gX/m2]; the T&C leakage
        %% formulation treats them as well mixed over the column water.
        P_cells = reshape(P(j,:,:),79,N_col).';
        [~,mobile] = Biogeo_Lateral_Pools(P_cells);
        mobile = mobile .* lateral_fraction;
        C_aq = Biogeo_Lateral_Concentration( ...
            mobile,A_cell,theta_cell,Zlayer_m,bio_active_fraction);

        lambda_cell = zeros(N_col,ms);
        capp_cell = zeros(N_col,ms);
        for c = 1:N_col
            p = PAR{c};
            [lambda_cell(c,:),capp_cell(c,:)] = Soil_Lateral_Thermal_State( ...
                T_cell(c,:),theta_cell(c,:)+Oice_cell(c,:), ...
                Pre(i),p,OPT_FR_SOIL);
        end
        thermal = struct('T_cell',T_cell,'lambda_cell',lambda_cell, ...
            'capp_cell',capp_cell,'Oice_cell',Oice_cell);
        thermal.parameter_cell = PAR;
        thermal.Pre = Pre(i);
        thermal.OPT_FR_SOIL = OPT_FR_SOIL;
        [theta_new, dM_cell, ldiag, T_new, Oice_new] = Lateral_Routing_Step( ...
            topo,theta_cell,C_aq,soilpar_cell,Zlayer_m,dth,thermal, ...
            bio_active_fraction);

        %% Write the redistributed moisture back, keeping V consistent
        %% with O through the column's own residual water content.
        for c = 1:N_col
            p = PAR{c};
            O(i,:,c) = theta_new(c,:);
            Oice(i,:,c) = Oice_new(c,:);
            V(i,:,c) = (theta_new(c,:) - p.Ohy) .* p.dz;
            Vice(i,:,c) = Oice_new(c,:) .* p.dz;
            Tdp(i,:,c) = T_new(c,:);
            [Tdp_H(i,:,c),Tdp_L(i,:,c)] = RootZone_Temp( ...
                T_new(c,:),RfH_Zs_c{c},RfL_Zs_c{c});
        end

        %% Return the advected mass to its T&C pools [gX/m2].
        dM_area = squeeze(sum(dM_cell,2)) ./ A_cell;
        P_cells = Biogeo_Lateral_Pools(P_cells,dM_area);
        P(j,:,:) = reshape(P_cells.',1,79,N_col);
        transported = unique([lateral_pool_idx{:}]);
        Pj = P(j,transported,:);
        if any(Pj(:) < -1e-10)
            error('MAIN_FRAME_BG_ERW_2D:NegativePool', ...
                'Conservative lateral transport produced a negative pool at i=%d.',i);
        end

        lateral_diag.water_buf_to_main(j) = lateral_diag.water_buf_to_main(j) + ldiag.buf_to_main.water;
        lateral_diag.mass_buf_to_main(j,:) = lateral_diag.mass_buf_to_main(j,:) + ldiag.buf_to_main.net.';
        lateral_diag.NO3_buf_to_main(j)   = lateral_diag.NO3_buf_to_main(j) + ldiag.buf_to_main.net(i_lateral_NO3);
        lateral_diag.sub_steps_max(j)     = max(lateral_diag.sub_steps_max(j), ldiag.sub_steps);
        lateral_diag.heat_buf_to_main_J(j)= lateral_diag.heat_buf_to_main_J(j) + ldiag.buf_to_main.heat_J;
        lateral_diag.heat_residual_J(j)   = lateral_diag.heat_residual_J(j) + ldiag.heat.energy_residual_J;
        lateral_diag.heat_residual_relative(j) = max( ...
            lateral_diag.heat_residual_relative(j), ...
            ldiag.heat.energy_residual_relative);
        lateral_diag.heat_sub_steps_max(j)= max(lateral_diag.heat_sub_steps_max(j),ldiag.heat.sub_steps);
    end

    %% ------------------------------------------------------------------
    %% ERW  --  once per day, per column
    %% ------------------------------------------------------------------
    if OPT_ERW == 1 && ((Datam(i,4) == 1) || i == NN)
        for c = 1:N_col
            ERW_PARAM_cell = ERW_PARAM_input;
            if isfield(ColCfg,'ERW_PARAM')
                ERW_PARAM_cell = ColCfg(c).ERW_PARAM;
            end
            ERW{c} = ERW_Daily_Step_2D(ERW{c}, ColCfg(c), PAR{c}, ...
                i, j, NN, NNd, Datam, Date, Ca, ...
                O(:,:,c), Lk(:,1,c), f(:,1,c), T_H(:,:,c), T_L(:,:,c), ...
                Tdp_bio_c(c), VSUM_c(c), VTSUM_c(c), ...
                R_microbe(:,1,c), R_litter(:,1,c), R_ew(:,1,c), R_litter_sur(:,1,c), ...
                Rmr_H(:,:,c), Rmr_L(:,:,c), Ccrown_c(c,:), ...
                P(:,:,c), MOD_PARAM_input, ERW_PARAM_cell, ...
                conv_mol, conv_Al);
            if ERW{c}.has_update
                PH_c(c)   = ERW{c}.PH;
                P(j,60,c) = ERW{c}.P60;  P(j,68,c) = ERW{c}.P68;
                P(j,52,c) = ERW{c}.P52;  P(j,76,c) = ERW{c}.P76;
                P(j,43,c) = ERW{c}.P43;
            end
        end
    end

    if N_col > 1 && OPT_ERW == 1 && ...
            all(cellfun(@(x) isfield(x,'init') && x.init,ERW))
        [ERW,erw_ldiag] = ERW_Lateral_Transport( ...
            ERW,topo,ldiag.Q_edge,theta_cell,Zlayer_m,bio_active_fraction, ...
            dth,Ca(i)/22.41);
        lateral_diag.erw_mass_buf_to_main_umol(j,:) = ...
            lateral_diag.erw_mass_buf_to_main_umol(j,:) + erw_ldiag.buf_to_main.aqueous_umol.';
        lateral_diag.erw_CO2_buf_to_main_umol(j) = ...
            lateral_diag.erw_CO2_buf_to_main_umol(j) + erw_ldiag.buf_to_main.CO2_umol;
        lateral_diag.erw_mass_residual_umol(j,:) = ...
            lateral_diag.erw_mass_residual_umol(j,:) + erw_ldiag.mass_residual_umol.';
    end
end

Computational_Time = toc;
fprintf('\nCOMPUTATIONAL TIME [h] %.3f   (%.1f ms per column-step)\n', ...
    Computational_Time/3600, 1000*Computational_Time/(NN*N_col));

%% ========================================================================
%% DOMAIN AGGREGATION  --  area-weighted over the sub-plot cells
%% ========================================================================
w_col = reshape(A_cell/A_tot, 1, N_col);
Domain.area_m2       = A_tot;
Domain.area_main_m2  = sum(A_cell(strcmp(topo.cell_region,'main')));
Domain.area_buffer_m2= sum(A_cell(~strcmp(topo.cell_region,'main')));
Domain.N2Oflx        = squeeze(N2Oflx)  * w_col.';
Domain.LEAK_NO3      = squeeze(LEAK_NO3)* w_col.';
Domain.LEAK_NH4      = squeeze(LEAK_NH4)* w_col.';
Domain.PH_d          = squeeze(PH_d)    * w_col.';
Domain.NO3_pool      = squeeze(P(:,32,:))* w_col.';
Domain.lateral       = lateral_diag;
Domain.exit_water_vert_m3 = exit_water_vert_m3;
Domain.exit_water_surface_m3 = exit_water_surface_m3;
Domain.exit_water_total_m3 = exit_water_vert_m3+exit_water_surface_m3;
Domain.cell_region   = topo.cell_region;
Domain.cell_area_m2  = A_cell;

fprintf('Domain: %.0f m2 total (%.0f main, %.0f margin)\n', ...
    Domain.area_m2, Domain.area_main_m2, Domain.area_buffer_m2);
fprintf('Cumulative N2O   : %.4f gN/m2\n', sum(Domain.N2Oflx,'omitnan'));
fprintf('Cumulative NO3   : %.4f gN/m2\n', sum(Domain.LEAK_NO3,'omitnan'));
fprintf('Net lateral water margin->field: %.4f m3\n', ...
    sum(lateral_diag.water_buf_to_main,'omitnan'));
fprintf('Net lateral NO3   margin->field: %.4f g\n', ...
    sum(lateral_diag.NO3_buf_to_main,'omitnan'));

%% ========================================================================
function out = ternary(cond, a, b)
if cond; out = a; else; out = b; end
end
