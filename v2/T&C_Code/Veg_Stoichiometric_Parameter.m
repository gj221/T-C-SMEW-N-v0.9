%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction Stoichiometric_Parameter     %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[Stoich]=Veg_Stoichiometric_Parameter(Nl,crop_type)
%%%%%%
%Nl = [15 -42 ]; %[gC/gN ] Leaf Carbon-Nitrogen ratio 
%Ns  = [50 50];%%%[ 50-330] Sapwood Carbon Nitrogen  [gC/gN] Sapwood
%N = [58 58]; %% [30- 60] [gC/gN] Fine root  Carbon Nitrogen
%Nc = Ns 
Stoich.Nl = Nl;  %% Leaf Carbon Nitrogen 
Stoich.Ns= Nl/0.145;%%%[ 50 -330] Sapwood Carbon Nitrogen  [gC/gN] Sapwood
Stoich.Nr= Nl/0.860; %%% [30- 60] [gC/gN] Fine root  Carbon Nitrogen
Stoich.Nf = Nl;  %%%  [gC/gN]  Fruit/Reproduction  
Stoich.Nh = Nl/0.145;  %%%  [gC/gN] Heartwood/Dead sapwood 
%%%% 
Pl = Nl*14;
Stoich.Phol = Pl; %%%  [gC/gP]
Stoich.Phos = Pl/0.145;
Stoich.Phor = Pl/0.860; 
Stoich.Phof = Pl; 
Stoich.Phoh= Pl/0.145; 
%%%%%
Kl = Nl*2;
Stoich.Kpotl = Kl;  %%  [gC/gK]
Stoich.Kpots  = Kl/0.145;
Stoich.Kpotr = Kl/0.20;
Stoich.Kpotf = Kl;
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
        if strcmp(crop_type,'Kelland')
            Cal = Nl*2.3;
            Mgl = Nl*7.5;
            Sil = Nl*1.3;
        
        elseif strcmp(crop_type,'Maize')
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
Stoich.Calf = Cal;
Stoich.Calh = Cal/0.195;
%%%%% Mag
% Mgl = Nl*2;
Stoich.Magl = Mgl;  %%  [gC/gMg]
Stoich.Mags  = Mgl/0.132;
Stoich.Magr = Mgl/0.80;
Stoich.Magf = Mgl;
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
