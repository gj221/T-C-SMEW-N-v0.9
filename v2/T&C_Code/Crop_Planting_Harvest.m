%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  CROP HARVESTING or PLANTING   %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[B,Nreserve,Preserve,Kreserve,Careserve,Mgreserve,Sireserve, RB,LAI,LAIdead,ManI,PHE_S,dflo,AgeL,AgeDL,SeedInput] ...
    = Crop_Planting_Harvest(B,Nreserve,Preserve,Kreserve,Careserve,Mgreserve,Sireserve,RBtm1,dtd,LAI,LAIdead,PHE_S,dflo,AgeL,AgeDL,Datam,Mpar)
%%%%%%
DD = datenum(Datam(1),Datam(2),Datam(3),Datam(4),0,0);
Date_sowing=Mpar.Date_sowing+15; %%% Emergence as Date of sowing + 15 days
Date_harvesting=Mpar.Date_harvesting;
Crop_B= Mpar.Crop_B;
%%%%%
%%%% INPUT
ManI=0; %%% Management Indicator (>0) Sowing season (-2) Harvest
Btm1=B;
SeedInput=zeros(1,7); %% [C N P K Ca Mg Si], [g element m-2 crown]
%%%%%%%%%% Planting
isowing=find(abs(DD-Date_sowing)<=0.49,1,'first');
if ~isempty(isowing) && sum(abs(Btm1(1:6)))<=sqrt(eps)
    if numel(Crop_B)~=2 || any(~isfinite(Crop_B)) || any(Crop_B<0)
        error('Crop_Planting_Harvest:InvalidSeedCarbon', ...
            'Crop_B must contain two finite nonnegative carbon pools [g C m-2 crown].');
    end
    SeedNut=[];
    if isfield(Mpar,'SeedNut_input')
        SeedNut=Mpar.SeedNut_input;
        if numel(SeedNut)~=6 || any(~isfinite(SeedNut)) || any(SeedNut<0)
            error('Crop_Planting_Harvest:InvalidSeedNutrients', ...
                'SeedNut_input must contain finite nonnegative [N P K Ca Mg Si] in g element m-2 crown.');
        end
        SeedNut=reshape(SeedNut,1,6);
    end
    if ~isfield(Mpar,'EstablishmentStoich')
        error('Crop_Planting_Harvest:MissingEstablishmentStoich', ...
            'Crop establishment requires root and storage element stoichiometry.');
    end
    if isempty(SeedNut)
        Stoich=Mpar.EstablishmentStoich;
        leaf_ratio=[Stoich.Nl Stoich.Phol Stoich.Kpotl Stoich.Call Stoich.Magl Stoich.Sill];
        store_ratio=[Stoich.Ns Stoich.Phos Stoich.Kpots Stoich.Cals Stoich.Mags Stoich.Sils];
        if any(~isfinite([leaf_ratio store_ratio])) || any([leaf_ratio store_ratio]<=0)
            error('Crop_Planting_Harvest:InvalidEstablishmentStoich', ...
                'Leaf and storage C:element ratios must be finite and positive.');
        end
        SeedNut=max(0,Crop_B(2).*(1./leaf_ratio-1./store_ratio));
    end
    SeedInput=Plant_Establishment_Element_Boundary( ...
        Crop_B,SeedNut,Mpar.EstablishmentStoich);
    ManI=isowing; %%% identify which sowing season
    %%%%%
    B(3) = Crop_B(1);
    B(4) = Crop_B(2);
    %%%% The managed emergence event is the crop leaf-onset boundary.
    PHE_S=2; dflo=1; AgeL=0; AgeDL=0;
    %%%%
    B(B<0)=0;
    %%%%
    Nreserve = Nreserve+SeedNut(1);
    Preserve = Preserve+SeedNut(2);
    Kreserve = Kreserve+SeedNut(3);

    Careserve = Careserve+SeedNut(4);
    Mgreserve = Mgreserve+SeedNut(5);
    Sireserve = Sireserve+SeedNut(6);
end
%%%%%%%%%% Harvesting
if  sum(abs(DD-Date_harvesting)<=0.49)>=1
    ManI=-2;
    %%%%% Completely harvested pools // crop removed
    B(1:7)=0; %%; %% [gC/ m^2 ]
    %%%%%%%
    B(B<0)=0;
    %%%%
    LAI=0; LAIdead=0;
end
%%%%%%%%%%%
if ManI~=0 %% planting or harvest
    RB(1)=(Btm1(1)-B(1))/dtd; %%  [gC/ m^2 day] %%% Leaves - Grass
    RB(2)=(Btm1(2)-B(2))/dtd; %%  [gC/ m^2 day] %%% Sapwood
    RB(3)=(Btm1(3)-B(3))/dtd; %%  [gC/ m^2 day] %%% Fine roots
    RB(4)=(Btm1(4)-B(4))/dtd; %%  [gC/ m^2 day] %%% Carbohydrate Reserve
    RB(5)=(Btm1(5)-B(5))/dtd; %%  [gC/ m^2 day] %%% Fruit and Flower
    RB(6)=(Btm1(6)-B(6))/dtd; %%  [gC/ m^2 day] %%% Heartwood - Dead Sapwood
    RB(7)=(Btm1(7)-B(7))/dtd; %%  [gC/ m^2 day] %%% Leaves - Grass -- Standing Dead
else
    %%% No Management
    RB(1:7)=0; %%
end
%%%% If crop still not out - re-initialize phenology
if (dflo==2) && (LAI==0)
    PHE_S=1; dflo= 0; AgeL= 0; AgeDL=0;
end
%%%%
RB=RB+RBtm1;
end
