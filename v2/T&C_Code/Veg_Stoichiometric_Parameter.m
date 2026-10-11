%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction Stoichiometric_Parameter     %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[Stoich]=Veg_Stoichiometric_Parameter(Nl,crop_type)
%%%%%%
%Nl = [15 -42 ]; %[gC/gN ] Leaf Carbon-Nitrogen ratio 
%Ns  = [50 50];%%%[ 50-330] Sapwood Carbon Nitrogen  [gC/gN] Sapwood
%N = [58 58]; %% [30- 60] [gC/gN] Fine root  Carbon Nitrogen
%Nc = Ns 
%%% FRUIT (grain) : LEAF concentration ratios. -----------------------------
%%% Every other organ here carries a measured ratio to leaf with a DOI cited
%%% below - stem 0.195, root 0.25 for Ca. The fruit column alone sat at 1.0
%%% for N, P, K, Ca, Mg and Si, uncited: an unparameterised placeholder, not
%%% a measurement. It matters because harvest exports ~97% fruit
%%% (fract_left_fr = 0) while the other fractions return residue, so these
%%% ratios set the ENTIRE harvest nutrient export.
%%%
%%% At 1.0 the result is not merely inaccurate but physically impossible:
%%% the model removed ~330 kg N/ha/yr in maize grain, more nitrogen than the
%%% field receives as fertiliser. Maize grain Ca came out at 0.800% against
%%% a published 0.007-0.02%, and soybean seed Ca 2.586% against ~0.25-0.30%.
%%%
%%% Ratios are f = (leaf divisor) * [X]_grain / f_C_grain, from published
%%% maize grain and soybean seed concentrations on a dry-matter basis.
%%%
%%% CAVEAT applying to every value here: f_C_grain, the carbon fraction of
%%% grain dry matter, is NOT defined anywhere in the codebase. These use
%%% 0.45 (maize, starch-rich) and 0.53 (soybean, oil/protein-rich), carrying
%%% roughly +/-15% into each ratio. That is immaterial against corrections
%%% of 2-100x, but the ABSOLUTE export magnitudes should not be trusted.
%%% Before publication: dump simulated fruit carbon at each harvest, divide
%%% by those fractions, compare with measured Energy Farm yields 2016-2020,
%%% and rescale if it disagrees by more than ~10%.
%%%
%%% Sources: Bender, Haegele, Ruffo & Below 2013 Agron. J. 105:161-170;
%%% Ciampitti, Camberato, Murrell & Vyn 2013 Agron. J. 105:783-795;
%%% Nafziger 2017 Illinois survey (n=2,335); Ohio State Ext. ANR-74
%%% (n=2,338); Bruns 2016 Comm. Soil Sci. Plant Anal. 47:2008-2016.
%%%
%%% NOT changed: the stem/root/heartwood ratios (legacy), and the biological
%%% uptake depth weighting (legacy). Those are modelling choices; this is a
%%% placeholder correction.
if exist('crop_type','var') && strcmp(crop_type,'Maize')
    f_N_fruit  = 0.48;   %% grain N 1.44% DM
    f_P_fruit  = 1.54;   %% grain P 0.33% DM (phytate-rich, ABOVE leaf)
    f_K_fruit  = 0.28;   %% grain K 0.42% DM vs leaf 1.5%
    f_Ca_fruit = 0.019;  %% grain Ca 0.015% DM (range 0.010-0.029)
    f_Mg_fruit = 0.40;   %% grain Mg 0.141% vs published 0.126-0.156%
elseif exist('crop_type','var') && strcmp(crop_type,'Soybean')
    f_N_fruit  = 1.40;   %% seed N 6.18% DM (protein-rich, ABOVE leaf)
    f_P_fruit  = 1.84;   %% seed P 0.58% DM
    f_K_fruit  = 1.00;   %% seed K ~1.8% genuinely matches leaf K (three
                         %% datasets within 5%) - left at 1.0 deliberately
    f_Ca_fruit = 0.13;   %% seed Ca 0.40% DM (range 0.095-0.18)
    f_Mg_fruit = 0.48;   %% seed Mg 0.283% DM; differs from maize because
                         %% soybean seed carries ~0.53 gC/gDM, not 0.45
else
    f_N_fruit = 1; f_P_fruit = 1; f_K_fruit = 1;
    f_Ca_fruit = 1; f_Mg_fruit = 1;   %% unchanged without a grain analysis
end
%%% -----------------------------------------------------------------------
Stoich.Nl = Nl;  %% Leaf Carbon Nitrogen 
Stoich.Ns= Nl/0.145;%%%[ 50 -330] Sapwood Carbon Nitrogen  [gC/gN] Sapwood
Stoich.Nr= Nl/0.860; %%% [30- 60] [gC/gN] Fine root  Carbon Nitrogen
Stoich.Nf = Nl/f_N_fruit;  %%%  [gC/gN]  Fruit/Reproduction  
Stoich.Nh = Nl/0.145;  %%%  [gC/gN] Heartwood/Dead sapwood 
%%%% 
Pl = Nl*14;
Stoich.Phol = Pl; %%%  [gC/gP]
Stoich.Phos = Pl/0.145;
Stoich.Phor = Pl/0.860; 
Stoich.Phof = Pl/f_P_fruit; 
Stoich.Phoh= Pl/0.145; 
%%%%%
Kl = Nl*2;
Stoich.Kpotl = Kl;  %%  [gC/gK]
Stoich.Kpots  = Kl/0.145;
Stoich.Kpotr = Kl/0.20;
Stoich.Kpotf = Kl/f_K_fruit;
Stoich.Kpoth = Kl/0.145;

% Leaf C-Ca/Mg/Si ratio: doi: 10.3389/fpls.2021.674932
% Sapwood C-Ca/Mg ratio: doi/10.1111/nph.13904#
%  Sapwood C-Si ratio: doi.org/10.1038/s41598-022-20662-1
% root C-Si ratio: doi.org/10.1186/s12870-023-04236-5
% root C-Ca/Mg ratio: doi.org/10.1007/s11104-023-06208-y 

switch nargin
    case 1
        Cal = Nl*1.3;
        Mgl = Nl*7.5;
        Sil = Nl*2;
    case 2
        if strcmp(crop_type,'Maize')
            Cal = Nl*3.75;
            Mgl = Nl*8.5;
            Sil = Nl*2;
            
        elseif strcmp(crop_type,'Soybean')
            Cal = Nl*1.45;
            Mgl = Nl*7.5;
            Sil = Nl*2;
        else
            % Cal = Nl*1.3;
            % Mgl = Nl*7.5;
            % Sil = Nl*2;
            disp('no corresponding crop type')
            return
        end
end


%%%%% Cal
% Cal = Nl*2;
Stoich.Call = Cal;  %%  [gC/gCa]
Stoich.Cals  = Cal/0.195;
Stoich.Calr = Cal/0.25;
Stoich.Calf = Cal/f_Ca_fruit;
Stoich.Calh = Cal/0.195;
%%%%% Mag
% Mgl = Nl*2;
Stoich.Magl = Mgl;  %%  [gC/gMg]
Stoich.Mags  = Mgl/0.132;
Stoich.Magr = Mgl/0.80;
Stoich.Magf = Mgl/f_Mg_fruit;
Stoich.Magh = Mgl/0.132;
%%%%% Sil
% Sil = Nl*2;
Stoich.Sill = Sil;  %%  [gC/gSi]
Stoich.Sils  = Sil/0.145;
Stoich.Silr = Sil/0.286;
Stoich.Silf = Sil;
Stoich.Silh = Sil/0.145;

%%%%%%%%%%%%%%%% Max. Translocation Rates of nutrients 
Stoich.ftransR = 0.2; 
Stoich.ftransL = 0.6; 
%%%%%% 
Stoich.FiS=1; %%[-] ; %%% Factor to increment or decrement nutrient reserve buffer 
%%%%%%
%%%%% Lignin Fraction 
%% Poorter 1994; Chave et al 2009 ; Roumet et al 2015 ; Fortunel et al 2009
Stoich.Lig_fr_l= 0.15; %% Lignin fraction in leaves  [g Lignin / g DM] 
Stoich.Lig_fr_fr=0.10; %%% Lignin fraction in fruit [g Lignin / g DM] 
Stoich.Lig_fr_h= 0.25; %% Lignin fraction in wood [g Lignin / g DM] 
Stoich.Lig_fr_r= 0.10; %% Lignin fraction in fine roots [g Lignin / g DM] 
return
