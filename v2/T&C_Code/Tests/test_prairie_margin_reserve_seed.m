function test_prairie_margin_reserve_seed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The prairie margin PFT must start from seed, with exactly the nutrient
%% reserve a sowing event provides -- neither the crop file's unresettable
%% nutrient-reserve placeholder nor an invented established stand.
%%
%% MOD_PARAM_Kantola_Maize_Soy.m seeds Careserve/Mgreserve/Sireserve at
%% 1000 g/m2 for both PFT slots, and, for the low PFT slot, Nreserve and
%% Kreserve at 1000 g/m2 and Preserve at 100 g/m2 ("legacy").  The crop
%% clears every one of these at its first harvest (CropHeightType.m, the
%% ManI==-2 branch) and re-seeds them at every sowing
%% (Crop_Planting_Harvest.m: Nreserve = Nreserve + SeedNut(1), ...).
%% PFT_Prairie_Margin.m sets Date_sowing = Date_harvesting = [] (the crop
%% machinery never sows or harvests the margin), so without an explicit
%% override nothing ever clears the placeholder for the margin.
%% Nutrients_Available.m's relative-concentration factors (rCac =
%% 1 + (Careserve - Casto)/TCALns, the same form for every element) are
%% uncapped, so a reserve stuck at 1000 g/m2 against an actual plant
%% nutrient content of a few g/m2 evaluates to hundreds, which both
%% inflates litter nutrient return without bound and saturates
%% Plant_Nutrient_Suppression at 1 (uptake fully suppressed): a one-way
%% nutrient pump into solution with no sink.  Traced against a Bouchet
%% trajectory (control, margin_prairie, run 7): this was the source of
%% the >0.4 eq/m2/d unbalanced base-cation influx that drove pore-water
%% pH past 11 by day 153.
%%
%% The override starts the margin from seed: the stand state is the
%% post-sowing state Crop_Planting_Harvest.m creates (fine root
%% B(3) = Crop_B(1), carbohydrate reserve B(4) = Crop_B(2), nothing else)
%% and every element's reserve is the SeedNut a sowing with that Crop_B
%% adds, Crop_B(2) * (1/leaf - 1/storage) at the PFT's own stoichiometry.
%% That is exactly the flexible storage capacity
%% Xsto = FiS*B(4)*(1/leaf - 1/storage) Nutrients_Available.m computes for
%% the same B(4) (FiS = 1, reserve-pool ratios equal storage ratios), so
%% every relative-concentration factor starts at its neutral value 1 and
%% Plant_Nutrient_Suppression starts partially, not fully, engaged: the
%% state a sown crop starts from.  The annual hay cut is kept.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n=== test_prairie_margin_reserve_seed ===\n');
test_dir = fileparts(mfilename('fullpath'));
code_dir = fileparts(test_dir);
addpath(code_dir,fullfile(code_dir,'Parameters'));

prairie_source = fileread(fullfile(code_dir,'Parameters','PFT_Prairie_Margin.m'));
%% One-line form: continuation ellipses removed and whitespace runs
%% collapsed, so multi-line statements can be matched literally.
flat_source = regexprep(regexprep(prairie_source,'\.\.\.\s*',''),'\s+',' ');

%% Source check: seed carbon, the seed-nutrient formula at each PFT's own
%% stoichiometry, and one reserve assignment per element and slot (slot 2
%% is empty and gets zero).
assert(contains(flat_source,'prairie_seed_B = [20 30];'), ...
    'test_prairie_margin_reserve_seed:MissingSeedB', ...
    'PFT_Prairie_Margin.m must define prairie_seed_B = [20 30] (sowing fine root and reserve, gC/m2 crown).');
seed_formula = ['max(0, prairie_seed_B(2) .* ([1/St_seed.Nl 1/St_seed.Phol ' ...
    '1/St_seed.Kpotl 1/St_seed.Call 1/St_seed.Magl 1/St_seed.Sill] - ' ...
    '[1/St_seed.Ns 1/St_seed.Phos 1/St_seed.Kpots 1/St_seed.Cals ' ...
    '1/St_seed.Mags 1/St_seed.Sils]));'];
reserve_vars = {'Nreserve','Preserve','Kreserve','Careserve','Mgreserve','Sireserve'};
pft_tags = {'H','L'};
for it = 1:numel(pft_tags)
    tag = pft_tags{it};
    assert(contains(flat_source,['St_seed = Stoich_' tag '(1); prairie_seed_' tag ' = ' seed_formula]), ...
        'test_prairie_margin_reserve_seed:MissingSeedFormula', ...
        'PFT_Prairie_Margin.m must compute prairie_seed_%s as Crop_B(2)*(1/leaf - 1/storage) at Stoich_%s(1).',tag,tag);
    for k = 1:numel(reserve_vars)
        pattern = sprintf('%s_%s\\(1,:\\)\\s*=\\s*\\[prairie_seed_%s\\(%d\\) 0\\];', ...
            reserve_vars{k},tag,tag,k);
        assert(~isempty(regexp(prairie_source,pattern,'once')), ...
            'test_prairie_margin_reserve_seed:MissingSeedReserve', ...
            'PFT_Prairie_Margin.m must set %s_%s(1,:) = [prairie_seed_%s(%d) 0].', ...
            reserve_vars{k},tag,tag,k);
    end
end

%% Source check: the stand starts from seed (the post-sowing state) and
%% the annual hay cut is kept.
stand_lines = { ...
    'LAI_H(1,:) = [0.0 0.0];','LAI_L(1,:) = [0.0 0.0];', ...
    'B_H(1,1,:) = [0 0 prairie_seed_B(1) prairie_seed_B(2) 0 0 0 0];', ...
    'B_L(1,1,:) = [0 0 prairie_seed_B(1) prairie_seed_B(2) 0 0 0 0];', ...
    'AgeL_H(1,:) = [0 0];','AgeL_L(1,:) = [0 0];', ...
    'Mpar_H(islot).jDay_cut = datenum(2009:1:2050,8,1);', ...
    'Mpar_L(islot).jDay_cut = datenum(2009:1:2050,8,1);'};
for k = 1:numel(stand_lines)
    assert(contains(flat_source,stand_lines{k}), ...
        'test_prairie_margin_reserve_seed:MissingSeedStand', ...
        'PFT_Prairie_Margin.m must contain "%s" (seed start, hay cut kept).',stand_lines{k});
end

%% Source check: the crop's own placeholder is untouched. The margin
%% override runs on top of it (Load_Cell_Parameters runs PARAM_IC then
%% PFT_OVERRIDE in the same workspace); if the crop file's seed ever
%% changed from 1000/100 the margin override above would need revisiting,
%% so pin the value this test's mechanism check assumes.
crop_source = fileread(fullfile(code_dir,'Parameters','MOD_PARAM_Kantola_Maize_Soy.m'));
assert(contains(crop_source,'Careserve_H(1,:)= [1000 1000]') && ...
       contains(crop_source,'Nreserve_L(1,:)= [1000 1000]') && ...
       contains(crop_source,'Preserve_L(1,:)= [100 100]'), ...
    'test_prairie_margin_reserve_seed:BaseFileChanged', ...
    'The crop base file''s reserve seed no longer matches what this fix assumes; re-check the override.');

%% Mechanism check, placeholder: replay Nutrients_Available.m and
%% Plant_Nutrient_Suppression.m directly on the coneflower's stoichiometry
%% (Veg_Stoichiometric_Parameter(Nl_L(1)), Nl_L(1) = 22 gC/gN) and the
%% established-stand biomass the margin used before the seed start
%% (B = [10 0 120 80 0 0 5 0], chosen so the replay stays below
%% Nutrients_Available's e_relN > 100 early return and the factors are
%% visible).  The unfixed 1000/100 g/m2 placeholder must produce
%% relative-concentration factors and a suppression saturated far outside
%% anything physical.
St = Veg_Stoichiometric_Parameter(22);
B = [10 0 120 80 0 0 5 0];
FNCtm1 = 1;
OPT_SoilBiogeochemistry = 1;
[~,e_relN_old,rNc_old,rPc_old,rKc_old,rCac_old,rMgc_old,rSic_old] = ...
    Nutrients_Available(B,FNCtm1,St,1000,100,1000,1000,1000,1000,OPT_SoilBiogeochemistry);
assert(rNc_old > 100 && rPc_old > 100 && rKc_old > 100 && ...
       rCac_old > 100 && rMgc_old > 100 && rSic_old > 100, ...
    'test_prairie_margin_reserve_seed:OldNotPathological', ...
    'The unfixed 1000/100 g/m2 seed no longer drives the relative-concentration factors past 100 (rNc=%.3g rPc=%.3g rKc=%.3g rCac=%.3g rMgc=%.3g rSic=%.3g); the reproduction case is stale.', ...
    rNc_old,rPc_old,rKc_old,rCac_old,rMgc_old,rSic_old);
assert(e_relN_old > 50, ...
    'test_prairie_margin_reserve_seed:OldNotPathological', ...
    'The unfixed seed no longer drives the photosynthesis N-efficiency factor (e_relN=%.3g) close to the 100 threshold that later raises "unrealistic Photosynthesis N efficiency".',e_relN_old);
Casto = St.FiS*B(4)*(1/St.Call-1/St.Cals);
Nsto  = St.FiS*B(4)*(1/St.Nl-1/St.Ns);
TCALns = B(1)/St.Call + B(3)/St.Calr + B(4)/St.Cals + B(5)/St.Calf;
TNITns = B(1)/St.Nl + B(3)/St.Nr + B(4)/St.Ns + B(5)/St.Nf;
SupCa_old = Plant_Nutrient_Suppression(1000,Casto,TCALns,'Ca');
SupN_old  = Plant_Nutrient_Suppression(1000,Nsto,TNITns,'N');
assert(SupCa_old == 1 && SupN_old == 1, ...
    'test_prairie_margin_reserve_seed:UptakeNotSuppressed', ...
    'The unfixed seed must fully suppress Ca and N uptake (SupCa=%.3g SupN=%.3g); the sink-removal half of the mechanism is stale.',SupCa_old,SupN_old);
fprintf(['PASS  unfixed 1000/100 g/m2 placeholder: e_relN=%.3g rCac=%.3g rMgc=%.3g rSic=%.3g ' ...
    'rNc=%.3g rPc=%.3g rKc=%.3g, SupCa=%.3g SupN=%.3g (uptake fully suppressed)\n'], ...
    e_relN_old,rCac_old,rMgc_old,rSic_old,rNc_old,rPc_old,rKc_old,SupCa_old,SupN_old);

%% Mechanism check, seed start: for each PFT evaluate the file's own seed
%% statement (extracted from the source, so the test exercises the text
%% the model runs) on the post-sowing biomass state and confirm that it
%% is the sowing SeedNut of Crop_Planting_Harvest.m for Crop_B =
%% prairie_seed_B, that it equals the flexible storage capacity
%% Nutrients_Available.m computes for the seed reserve pool, that every
%% relative-concentration factor is neutral, and that uptake suppression
%% starts partial (0 < Sup < 1), as for a sown crop.
prairie_seed_B = [20 30];
B_seed = [0 0 prairie_seed_B(1) prairie_seed_B(2) 0 0 0 0];
element = {'N','P','K','Ca','Mg','Si'};
cases = struct('tag',{'H','L'},'Nl',{35,22});   % PFT_Prairie_Margin.m: Nl_H(1), Nl_L(1)
for ic = 1:numel(cases)
    tag = cases(ic).tag;
    St_seed = Veg_Stoichiometric_Parameter(cases(ic).Nl);
    stmt = regexp(flat_source,['prairie_seed_' tag ' = max\(0,.*?\]\)\);'],'match','once');
    assert(~isempty(stmt), ...
        'test_prairie_margin_reserve_seed:MissingSeedFormula', ...
        'Could not extract the prairie_seed_%s statement from PFT_Prairie_Margin.m.',tag);
    eval(stmt);
    seed_pft = eval(['prairie_seed_' tag]);
    leaf_ratio  = [St_seed.Nl St_seed.Phol St_seed.Kpotl St_seed.Call St_seed.Magl St_seed.Sill];
    store_ratio = [St_seed.Ns St_seed.Phos St_seed.Kpots St_seed.Cals St_seed.Mags St_seed.Sils];
    root_ratio  = [St_seed.Nr St_seed.Phor St_seed.Kpotr St_seed.Calr St_seed.Magr St_seed.Silr];
    SeedNut = max(0,prairie_seed_B(2).*(1./leaf_ratio-1./store_ratio));   % Crop_Planting_Harvest.m sowing branch
    assert(isequal(size(seed_pft),[1 6]) && all(isfinite(seed_pft)) && ...
           all(abs(seed_pft-SeedNut) <= 1e-15), ...
        'test_prairie_margin_reserve_seed:SeedNotSowing', ...
        'prairie_seed_%s [%s] is not the sowing SeedNut [%s] Crop_Planting_Harvest.m would add for Crop_B = [%g %g].', ...
        tag,sprintf('%.4g ',seed_pft),sprintf('%.4g ',SeedNut),prairie_seed_B);
    Xsto = St_seed.FiS*B_seed(4)*(1./leaf_ratio-1./store_ratio);   % Nutrients_Available.m flexible storage (Nc = Ns, ...)
    TXns = B_seed(3)./root_ratio + B_seed(4)./store_ratio;           % non-structural inventory with B(1) = B(5) = 0
    assert(all(seed_pft > 0) && all(abs(seed_pft-Xsto) <= 1e-12*Xsto), ...
        'test_prairie_margin_reserve_seed:SeedNotCapacity', ...
        'prairie_seed_%s [%s] must equal the flexible storage capacity Xsto [%s] of the seed reserve pool.', ...
        tag,sprintf('%.4g ',seed_pft),sprintf('%.4g ',Xsto));
    [~,e_relN_new,rNc_new,rPc_new,rKc_new,rCac_new,rMgc_new,rSic_new] = ...
        Nutrients_Available(B_seed,FNCtm1,St_seed,seed_pft(1),seed_pft(2),seed_pft(3), ...
        seed_pft(4),seed_pft(5),seed_pft(6),OPT_SoilBiogeochemistry);
    r_new = [rNc_new rPc_new rKc_new rCac_new rMgc_new rSic_new];
    assert(all(abs(r_new-1) <= 1e-12) && abs(e_relN_new-1) <= 1e-12, ...
        'test_prairie_margin_reserve_seed:NewNotNeutral', ...
        'The %s seed reserve must leave every relative-concentration factor at its neutral value 1 (rNc=%.3g rPc=%.3g rKc=%.3g rCac=%.3g rMgc=%.3g rSic=%.3g e_relN=%.3g).', ...
        tag,r_new,e_relN_new);
    Sup = zeros(1,6);
    for k = 1:6
        Sup(k) = Plant_Nutrient_Suppression(seed_pft(k),Xsto(k),TXns(k),element{k});
        Sup_expected = 0.2*Xsto(k)/(0.6*TXns(k)+0.2*Xsto(k));   % reserve == capacity
        assert(Sup(k) > 0 && Sup(k) < 1 && abs(Sup(k)-Sup_expected) <= 1e-12, ...
            'test_prairie_margin_reserve_seed:SuppressionNotPartial', ...
            'The %s seed reserve must leave %s uptake partially suppressed (Sup=%.4g, expected %.4g, strictly between 0 and 1), as for a sown crop.', ...
            tag,element{k},Sup(k),Sup_expected);
    end
    fprintf(['PASS  %s seed (Nl=%g gC/gN): reserve [N P K Ca Mg Si] = [%s] g/m2 crown equals SeedNut and Xsto, ' ...
        'every rXc = 1, e_relN = 1, Sup = [%s]\n'],tag,cases(ic).Nl, ...
        strtrim(sprintf('%.4g ',seed_pft)),strtrim(sprintf('%.3f ',Sup)));
end
fprintf('PASS  prairie margin starts from seed (B = [%s]) with the sowing reserve; hay cut kept\n', ...
    strtrim(sprintf('%g ',B_seed)));
end
