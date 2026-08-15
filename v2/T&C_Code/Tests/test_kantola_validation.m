function test_kantola_validation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Kantola annual harvest-carbon observation comparison.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n=== test_kantola_validation ===\n');
n_pass = 0;
out_dir = tempname;
mkdir(out_dir);
figures_before = findall(groot,'Type','figure');
cleanup = onCleanup(@() tidy(out_dir,figures_before));

obs = Load_Kantola_Harvest_Observations();
check(height(obs) == 8 && all(obs.SampleSize == 5), ...
    'raw harvest observations retain eight treatment-year groups');
expected = [848.990 35.622692; 794.350 44.780172; ...
    775.166 42.480430; 770.664 30.716973; ...
    294.274 7.266192; 249.048 18.230858; ...
    830.070 22.465320; 732.096 45.031745];
check(max(abs([obs.Mean_gC_m2 obs.SEM_gC_m2]-expected),[],'all') < 5e-7, ...
    'harvest means and SEM are calculated from raw block data');
check(all(obs.Units == "g C m-2") && ...
    all(obs.SourceKind == "raw_block_data"), ...
    'harvest carbon uses measured units and raw-data provenance');
check(raises_error(@() Load_Kantola_Harvest_Observations( ...
    fullfile(out_dir,'missing.csv')), ...
    'Load_Kantola_Harvest_Observations:MissingFile'), ...
    'loader rejects a missing observation file');

obs_file = fullfile(out_dir,'harvest.csv');
write_fixture_harvest(obs_file);
write_fixture_results(out_dir);
output_file = fullfile(out_dir,'validation.png');
[validation,fig,plot_file] = Plot_Kantola_Validation( ...
    'ResultDirectory',out_dir,'Runs',[10 11], ...
    'ObservationFile',obs_file,'Visible','off','OutputFile',output_file);
check(height(validation) == 4, 'validation returns one row per year and treatment');
check(strcmp(fig.Visible,'off') && isfile(plot_file), ...
    'headless validation figure is written to disk');
check(isequal(sort(validation.ModelGrainC_gC_m2).',[30 60 60 120]), ...
    'grain C excludes the seed row and separates calendar years');
check(all(validation.ModelGrainC_gC_m2 < 1000), ...
    'prairie-margin cells are excluded from the cropped footprint');
check(isscalar(findall(fig,'Type','axes','Tag','KantolaGrainC')) && ...
    isscalar(findall(fig,'Type','axes','Tag','KantolaGrainCResidual')), ...
    'figure contains harvest-C comparison and residual subplots');
grain_axis = findall(fig,'Type','axes','Tag','KantolaGrainC');
check(contains(string(grain_axis.Title.String),'n = 2'), ...
    'figure reports the sample size from the selected observation file');
check(numel(findall(fig,'Tag','ModelGrainC')) == 2 && ...
    numel(findall(fig,'Tag','ObservedGrainC')) == 2 && ...
    isempty(findall(fig,'-regexp','Tag','GrainN')), ...
    'figure contains raw harvest observations and no grain-N points');
plot_axes = findall(fig,'Type','axes');
check(all(abs(vertcat(plot_axes.Color)-1) < 1e-12,'all') && ...
    all(abs(vertcat(plot_axes.XColor)) < 1e-12,'all') && ...
    all(abs(vertcat(plot_axes.YColor)) < 1e-12,'all'), ...
    'validation axes export with a white background and black labels');
close(fig);

empty_dir = tempname;
mkdir(empty_dir);
empty_cleanup = onCleanup(@() rmdir(empty_dir,'s'));
check(raises_error(@() Plot_Kantola_Validation( ...
    'ResultDirectory',empty_dir,'Runs',[10 11], ...
    'ObservationFile',obs_file,'Visible','off'), ...
    'Plot_Kantola_Validation:MissingResult'), ...
    'plot rejects a missing result file');
clear empty_cleanup

RM = Run_Matrix_12();
bad_file = fullfile(out_dir,sprintf('TC_2D_run%02d_%s.mat',11,RM(11).name));
S = load(bad_file);
S.Ccrown_t(3,1,1) = 1.1;
save(bad_file,'-struct','S');
open_before = findall(groot,'Type','figure');
check(raises_error(@() Plot_Kantola_Validation( ...
    'ResultDirectory',out_dir,'Runs',[10 11], ...
    'ObservationFile',obs_file,'Visible','off'), ...
    'Plot_Kantola_Validation:UnphysicalCarbonState'), ...
    'plot rejects an invalid crown fraction');
check(numel(findall(groot,'Type','figure')) == numel(open_before), ...
    'validation errors do not leak figures');

write_fixture_results(out_dir);
S = load(bad_file);
S.Ccrown_t(2,1,1) = 0.6;
S.Ccrown_t(2,2,1) = 0.6;
save(bad_file,'-struct','S');
check(raises_error(@() Plot_Kantola_Validation( ...
    'ResultDirectory',out_dir,'Runs',[10 11], ...
    'ObservationFile',obs_file,'Visible','off'), ...
    'Plot_Kantola_Validation:UnphysicalCarbonState'), ...
    'plot rejects overlapping crop crowns');

write_fixture_results(out_dir);
S = load(bad_file);
S.Ccrown_t = permute(S.Ccrown_t,[2 1 3]);
save(bad_file,'-struct','S');
check(raises_error(@() Plot_Kantola_Validation( ...
    'ResultDirectory',out_dir,'Runs',[10 11], ...
    'ObservationFile',obs_file,'Visible','off'), ...
    'Plot_Kantola_Validation:CarbonTimeLengthMismatch'), ...
    'plot rejects harvest arrays with transposed dimensions');

write_fixture_results(out_dir,datetime(2018,8,31));
check(raises_error(@() Plot_Kantola_Validation( ...
    'ResultDirectory',out_dir,'Runs',[10 11], ...
    'ObservationFile',obs_file,'Visible','off'), ...
    'Plot_Kantola_Validation:IncompleteModelYear'), ...
    'plot rejects a partial observation year');

write_fixture_results(out_dir);
outside_file = fullfile(out_dir,'outside_period.csv');
outside = readtable(obs_file,'TextType','string','VariableNamingRule','preserve');
outside.Year(outside.Year == 2018) = 2019;
writetable(outside,outside_file);
check(raises_error(@() Plot_Kantola_Validation( ...
    'ResultDirectory',out_dir,'Runs',[10 11], ...
    'ObservationFile',outside_file,'Visible','off'), ...
    'Plot_Kantola_Validation:ObservationOutsideModelPeriod'), ...
    'plot rejects observations outside the saved model period');

fprintf('  ---------------------------------------\n');
fprintf('  %d passed, 0 failed\n',n_pass);

    function check(ok,label)
        assert(ok,'test_kantola_validation:Failed','%s',label);
        n_pass = n_pass + 1;
        fprintf('  PASS  %s\n',label);
    end

end

function write_fixture_harvest(file)
Year = repelem([2017; 2017; 2018; 2018],2);
Crop = repelem(["Maize"; "Maize"; "Maize"; "Maize"],2);
Block = repmat([1; 2],4,1);
Treatment = repelem(["Control"; "Basalt"; "Control"; "Basalt"],2);
AbovegroundPeakBiomass = zeros(8,1);
BelowgroundPeakBiomass = zeros(8,1);
HarvestBiomassGrain = [10; 14; 16; 20; 20; 24; 36; 40];
T = table(Year,Crop,Block,Treatment,AbovegroundPeakBiomass, ...
    BelowgroundPeakBiomass,HarvestBiomassGrain,'VariableNames', ...
    {'Year','Crop','Block','Treatment','Aboveground Peak Biomass', ...
    'Belowground Peak Biomass','Harvest Biomass (Grain)'});
writetable(T,file);
end

function write_fixture_results(out_dir,last_day)
if nargin < 2
    last_day = datetime(2018,12,31);
end
RM = Run_Matrix_12();
times = (datetime(2017,1,1,1,0,0):hours(1):last_day+hours(23)).';
Datam = datevec(times);
daily_rows = find(Datam(:,4) == 1);
process_dates = datetime(Datam(daily_rows,1), ...
    Datam(daily_rows,2),Datam(daily_rows,3));
n_time = numel(daily_rows)+1;
cc = 2;
n_col = 2;
topo.cell_area_m2 = [3; 1];
ColCfg = repmat(struct('PFT_OVERRIDE',''),n_col,1);
ColCfg(2).PFT_OVERRIDE = 'PFT_Prairie_Margin';
for run_id = [10 11]
    R = RM(run_id);
    scale = run_id-9;
    RB_L = zeros(n_time,cc,7,n_col);
    Ccrown_t = zeros(n_time,cc,n_col);
    RB_L(1,1,5,1) = 3e6;
    Ccrown_t(1,1,1) = 1;
    harvest_rows = [find(process_dates == datetime(2017,12,31),1); ...
        find(process_dates == datetime(2018,1,1),1)]+1;
    RB_L(harvest_rows,1,5,1) = [30; 60]*scale;
    Ccrown_t(harvest_rows,1,1) = 1;
    RB_L(harvest_rows,1,5,2) = 3e6;
    Ccrown_t(harvest_rows,1,2) = 1;
    save(fullfile(out_dir,sprintf('TC_2D_run%02d_%s.mat',run_id,R.name)), ...
        'R','Datam','RB_L','Ccrown_t','topo','ColCfg');
end
end

function ok = raises_error(f,id)
ok = false;
try
    f();
catch ME
    ok = strcmp(ME.identifier,id);
end
end

function tidy(out_dir,figures_before)
figures = findall(groot,'Type','figure');
for k = 1:numel(figures)
    if ~any(arrayfun(@(x) isequal(x,figures(k)),figures_before))
        close(figures(k));
    end
end
if isfolder(out_dir)
    rmdir(out_dir,'s');
end
end
