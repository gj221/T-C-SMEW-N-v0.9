function test_crop_establishment_element_boundary()
%% Crop establishment must carry every element implied by initialized tissue.

code_dir = fileparts(fileparts(mfilename('fullpath')));
param_file = fullfile(code_dir,'Parameters','MOD_PARAM_Kantola_Maize_Soy.m');
crop_file = fullfile(code_dir,'Crop_Planting_Harvest.m');

maize = Veg_Stoichiometric_Parameter(15,"Maize");
soy = Veg_Stoichiometric_Parameter(12,"Soybean");

maize_c = [20 30];   %% legacy 50 g C m-2 seed initializer
soy_c   = [20 30];   %% legacy 50 g C m-2 seed initializer
zero_seed = zeros(1,6);
extra_seed = [0.03 0.004 0.01 0.002 0.001 0];

maize_input = Plant_Establishment_Element_Boundary(maize_c,zero_seed,maize);
soy_input = Plant_Establishment_Element_Boundary(soy_c,extra_seed,soy);

assert(abs(sum(maize_c)-50)<1e-12 && abs(sum(soy_c)-50)<1e-12, ...
    'Crop establishment C must match the legacy 50 g C m-2 seed initializer.');
assert(max(abs(maize_input-[sum(maize_c) implied_elements(maize_c,maize)]))<1e-14, ...
    'Maize establishment does not conserve the initialized tissue elements.');
assert(max(abs(soy_input-[sum(soy_c) implied_elements(soy_c,soy)+extra_seed]))<1e-14, ...
    'Explicit seed reserves must add once to the soybean tissue boundary.');
assert(all(Plant_Establishment_Element_Boundary([0 0],zero_seed,maize)==0), ...
    'Zero establishment must not create an element boundary.');

expect_error(@() Plant_Establishment_Element_Boundary([-1 1],zero_seed,maize), ...
    'Plant_Establishment_Element_Boundary:InvalidCarbon');
expect_error(@() Plant_Establishment_Element_Boundary(maize_c,[-1 zeros(1,5)],maize), ...
    'Plant_Establishment_Element_Boundary:InvalidNutrients');

param_source = fileread(param_file);
crop_source = fileread(crop_file);
assert(contains(param_source,'Crop_B=[20 30]'), ...
    'The legacy 50 g C m-2 crop initializer is not active.');
seed_lai = [0.034*maize_c(2) 0.033*soy_c(2)];
lai_floor_token = regexp(param_source,'LAI_min_L\s*=\s*\[([^\]]+)\]', ...
    'tokens','once');
assert(~isempty(lai_floor_token), ...
    'The managed-crop minimum LAI parameter is missing.');
lai_floor = sscanf(lai_floor_token{1},'%f').';
assert(numel(lai_floor)==2 && all(lai_floor>=0) && ...
        all(seed_lai>lai_floor), ...
    'The managed-crop LAI floor suppresses the seed-scale emergence canopy.');
assert(contains(crop_source,'Plant_Establishment_Element_Boundary'), ...
    'Crop planting does not record the tissue element boundary.');

plot_input = Plant_Seed_Input_Accumulator(zeros(1,7),maize_input,0.4);
assert(max(abs(plot_input-0.4*maize_input))<1e-14, ...
    'Per-crown establishment input was not converted once to ground area.');
assert(max(abs(Plant_Seed_Input_Accumulator(plot_input,zeros(1,7),0.4)-plot_input))<1e-14, ...
    'A repeated no-event call changed the establishment boundary.');

fprintf('test_crop_establishment_element_boundary: PASS\n');
end

function value = implied_elements(crop_c,stoich)
root_ratio = [stoich.Nr stoich.Phor stoich.Kpotr stoich.Calr stoich.Magr stoich.Silr];
store_ratio = [stoich.Ns stoich.Phos stoich.Kpots stoich.Cals stoich.Mags stoich.Sils];
value = crop_c(1)./root_ratio + crop_c(2)./store_ratio;
end

function expect_error(f,id)
raised = false;
try
    f();
catch ME
    raised = strcmp(ME.identifier,id);
end
assert(raised,'Expected error %s.',id);
end
