function test_prairie_mow_management
%% Prairie is grass (aSE==2), so VEGGIE_UNIT routes it through
%% Grass_Cut_Grazing, not Crop_Planting_Harvest (which runs only for aSE==5).
%% Haying is therefore configured with jDay_cut and LAI_cut on Mpar.
test_dir = fileparts(mfilename('fullpath'));
code_dir = fileparts(test_dir);

pfm = fileread(fullfile(code_dir,'Parameters','PFT_Prairie_Margin.m'));
assert(contains(pfm,'aSE_H = [2 2]') && contains(pfm,'aSE_L = [2 2]'), ...
    'Prairie must be grass (aSE==2) for Grass_Cut_Grazing to be reachable.');
assert(contains(pfm,'Mpar_H(islot).jDay_cut') && contains(pfm,'Mpar_L(islot).jDay_cut'), ...
    'Prairie must set jDay_cut, or no cut ever fires (default NaN).');
assert(contains(pfm,'Mpar_H(islot).LAI_cut') && contains(pfm,'Mpar_L(islot).LAI_cut'), ...
    'Prairie must set LAI_cut, the residual canopy after the cut.');
assert(contains(pfm,'Nreserve_L(1,:)  = [prairie_seed_L(1) 0];'), ...
    ['Prairie must seed its own reserves from the sowing state, not keep ' ...
     'the crop seed reserve, which would persist all run and trip the e_relN guard.']);

vu = fileread(fullfile(code_dir,'VEGGIE_UNIT.m'));
assert(contains(vu,'jDay_cut = Mpar.jDay_cut;') && contains(vu,'LAI_cut = Mpar.LAI_cut;'), ...
    'VEGGIE_UNIT must read both cut parameters from Mpar.');
assert(contains(vu,'Grass_Cut_Grazing('), ...
    'The grass cut path must still exist.');

cph = fileread(fullfile(code_dir,'Crop_Planting_Harvest.m'));
assert(~contains(cph,'MowKeepFraction'), ...
    'Crop_Planting_Harvest must stay legacy: prairie never reaches it.');

fprintf('PASS: prairie hay cut wired through Grass_Cut_Grazing.\n');
end
