function test_crop_seed_nutrient_boundary()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Crop establishment cannot create nutrients outside a seed boundary.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

code_dir=fileparts(fileparts(mfilename('fullpath')));
addpath(code_dir);

Mpar.Date_sowing=datenum(2020,1,1)-15;
Mpar.Date_harvesting=datenum(2020,12,31);
Mpar.Crop_B=[0.52 0.78];
Mpar.EstablishmentStoich=Veg_Stoichiometric_Parameter(15,"Maize");
Mpar.NPK_res_ini=[20 3 20 0.9 1.5 7];
Datam=[2020 1 1 0];

[B,reserve,SeedInput,ManI,PHE_S,dflo]=plant_once( ...
    Mpar,Datam,zeros(1,7),zeros(1,6));
assert(ManI==1,'The emergence event was not detected.');
assert(PHE_S==2 && dflo==1, ...
    'A successful emergence event did not activate maximum-growth phenology.');
Stoich=Mpar.EstablishmentStoich;
leaf_ratio=[Stoich.Nl Stoich.Phol Stoich.Kpotl Stoich.Call Stoich.Magl Stoich.Sill];
store_ratio=[Stoich.Ns Stoich.Phos Stoich.Kpots Stoich.Cals Stoich.Mags Stoich.Sils];
expected_reserve=max(0,Mpar.Crop_B(2).*(1./leaf_ratio-1./store_ratio));
assert(max(abs(reserve-expected_reserve))<1e-14, ...
    'Default emergence reserve did not support storage-to-leaf stoichiometry.');
transfer_scale=Plant_Reserve_Growth_Limit(reserve,1,zeros(1,6), ...
    Mpar.Crop_B(2).*(1./leaf_ratio-1./store_ratio));
assert(abs(transfer_scale-1)<1e-14, ...
    'Default emergence reserve did not permit storage-to-leaf transfer.');
expected=Plant_Establishment_Element_Boundary( ...
    Mpar.Crop_B,expected_reserve,Mpar.EstablishmentStoich);
assert(max(abs(SeedInput-expected))<1e-14, ...
    'Default seed nutrients were not recorded on the explicit boundary.');

Mpar.SeedNut_input=[0.12 0.03 0.08 0.02 0.01 0.04];
[B,reserve,SeedInput]=plant_once(Mpar,Datam,zeros(1,7),zeros(1,6));
assert(max(abs(reserve-Mpar.SeedNut_input))<1e-14, ...
    'Explicit seed nutrients did not enter the reserve exactly once.');
expected=Plant_Establishment_Element_Boundary( ...
    Mpar.Crop_B,Mpar.SeedNut_input,Mpar.EstablishmentStoich);
assert(max(abs(SeedInput-expected))<1e-14, ...
    'Explicit seed nutrients were not recorded on the boundary.');

[~,reserve_repeat,SeedInput_repeat,ManI_repeat]=plant_once( ...
    Mpar,Datam,B,reserve);
assert(ManI_repeat==0 && all(SeedInput_repeat==0), ...
    'A repeated call on the emergence date reapplied the seed boundary.');
assert(max(abs(reserve_repeat-reserve))<1e-14, ...
    'A repeated emergence-date call changed reserve storage.');

standing=zeros(1,7);
standing(7)=0.25;
[~,reserve_standing,SeedInput_standing,ManI_standing]=plant_once( ...
    Mpar,Datam,standing,zeros(1,6));
assert(ManI_standing==1 && any(SeedInput_standing>0), ...
    'Standing-dead leaf carbon prevented the next managed crop emergence.');
assert(max(abs(reserve_standing-Mpar.SeedNut_input))<1e-14, ...
    'Standing-dead leaf carbon changed the managed seed reserve boundary.');

f_old=0.25;
f_new=0.75;
reserve_scaled=reserve*f_old/f_new;
assert(max(abs(f_new*reserve_scaled-f_old*reserve))<1e-14, ...
    'One-to-many crown rescaling did not conserve plot-area reserve.');

plot_input=Plant_Seed_Input_Accumulator(zeros(1,7),SeedInput,0.4);
assert(max(abs(plot_input-0.4*SeedInput))<1e-14, ...
    'Per-crown seed input was not converted to plot-area storage.');

bad=Mpar;
bad.SeedNut_input=[1 2 3];
assert_throws(@() plant_once(bad,Datam,zeros(1,7),zeros(1,6)), ...
    'Crop_Planting_Harvest:InvalidSeedNutrients');

fprintf('PASS: crop establishment has an explicit, area-consistent element boundary.\n');

end

function [B,reserve,SeedInput,ManI,PHE_S,dflo]=plant_once(Mpar,Datam,B,reserve)
[B,reserve(1),reserve(2),reserve(3),reserve(4),reserve(5),reserve(6), ...
    ~,~,~,ManI,PHE_S,dflo,~,~,SeedInput]=Crop_Planting_Harvest( ...
    B,reserve(1),reserve(2),reserve(3),reserve(4),reserve(5),reserve(6), ...
    zeros(1,7),1,0,0,1,0,0,0,Datam,Mpar);
end

function assert_throws(fcn,identifier)
threw=false;
try
    fcn();
catch ME
    threw=strcmp(ME.identifier,identifier);
end
assert(threw,'Expected error %s was not raised.',identifier);
end
