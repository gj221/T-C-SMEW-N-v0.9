function SeedInput=Plant_Establishment_Element_Boundary(Crop_B,SeedNut,Stoich)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Element boundary required by initialized crop carbon and seed reserves.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if numel(Crop_B)~=2 || any(~isfinite(Crop_B)) || any(Crop_B<0)
    error('Plant_Establishment_Element_Boundary:InvalidCarbon', ...
        'Crop_B must contain finite nonnegative [root store] carbon in g C m-2 crown.');
end
if numel(SeedNut)~=6 || any(~isfinite(SeedNut)) || any(SeedNut<0)
    error('Plant_Establishment_Element_Boundary:InvalidNutrients', ...
        'Seed nutrients must contain finite nonnegative [N P K Ca Mg Si] in g element m-2 crown.');
end

Crop_B=reshape(Crop_B,1,2);
SeedNut=reshape(SeedNut,1,6);
root_ratio=[Stoich.Nr Stoich.Phor Stoich.Kpotr Stoich.Calr Stoich.Magr Stoich.Silr];
store_ratio=[Stoich.Ns Stoich.Phos Stoich.Kpots Stoich.Cals Stoich.Mags Stoich.Sils];
if any(~isfinite(root_ratio)) || any(root_ratio<=0) || ...
        any(~isfinite(store_ratio)) || any(store_ratio<=0)
    error('Plant_Establishment_Element_Boundary:InvalidStoichiometry', ...
        'Root and storage stoichiometry must be finite and positive in g C per g element.');
end

tissue_input=Crop_B(1)./root_ratio+Crop_B(2)./store_ratio;
SeedInput=[sum(Crop_B) tissue_input+SeedNut];
end
