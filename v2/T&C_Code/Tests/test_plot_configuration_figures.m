function test_plot_configuration_figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Configuration figures contain every treatment in every run group, the
%% parameter-uncertainty ensemble draws as bands and range bars under the
%% 1-D results, and scenario-tagged results are found by name.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

v2 = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(v2,'run'),fullfile(v2,'plotting'), ...
    fullfile(v2,'T&C_Code'),fullfile(v2,'T&C_Code','Parameters'));

out_dir = tempname;
mkdir(out_dir);
cleanup = onCleanup(@() tidy(out_dir));
RM = Run_Matrix_12();
n_time = 4;

wfps_test = Daily_WFPS_From_Hourly( ...
    repelem([0.1; 0.2; 0.3],24), ...
    struct('Bio_Zs',1,'Soil_Param',struct('Osat',0.5)),n_time);
%% Row j averages the hours BIOGEO used for row j: hour 1 for row 2, then
%% 01:00 of one day to 00:00 of the next.
assert(max(abs(wfps_test.'-[NaN 20 (23*20+40)/24 (23*40+60)/24]),[],'omitnan') < 1e-12 && ...
    isnan(wfps_test(1)), ...
    'WFPS must average the same hours as the BIOGEO row it is stored in.');
assert_error(@() Daily_WFPS_From_Hourly(zeros(71,1), ...
    struct('Bio_Zs',1,'Soil_Param',struct('Osat',0.5)),n_time), ...
    'Plot_1D_Then_All_12:WFPSLength');

for k = 1:numel(RM)
    write_fixture(fullfile(out_dir,sprintf('TC_2D_run%02d_%s__1d.mat', ...
        k,RM(k).name)),fixture_series(k,RM(k),n_time));
end

[~,~,~,~,fig,~,time_fig] = run_plot_script(10:12,out_dir,out_dir,1,'off','','');

bars = findall(fig(1),'Type','bar','Tag','CumulativeN2O');
assert(numel(bars) == 3 && isequal(sort([bars.UserData]),10:12), ...
    'A single configuration must retain three treatment bars.');
treatment_colors = [0 0 0; 0 0.25 0.60; 0.35 0.70 0.90];
for run_id = 10:12
    treatment = run_id-9;
    bar_handle = bars([bars.UserData] == run_id);
    assert(isscalar(bar_handle) && ...
        isequal(bar_handle.FaceColor,treatment_colors(treatment,:)), ...
        'Single-configuration bars must retain their treatment colours.');
    assert(isequal(bar_handle.YData,6*run_id), ...
        'Single-configuration bars must use the matching run result.');
end
close([fig(:); time_fig(:)]);

[fig1d,~,stats_1d,~,fig,plot_file,time_fig,time_plot_file] = ...
    run_plot_script([],out_dir,out_dir,1,'off','','');
check_group_figure(fig1d,stats_1d);
check_figure_set(fig,plot_file,time_fig,time_plot_file,v2);
assert(isempty(findall(fig1d,'Type','patch','Tag','UncertaintyBand')) && ...
    all(arrayfun(@(f) isempty(findall(f,'Type','patch','Tag','UncertaintyBand')),time_fig)) && ...
    all(arrayfun(@(f) isempty(findall(f,'Type','errorbar','Tag','UncertaintyRange')),fig)), ...
    'Without an uncertainty file no band or range bar is drawn.');
close([fig1d; fig(:); time_fig(:)]);

%% The 1-D parameter ensemble: three envelope members whose base is the
%% fixture itself, the other two symmetrically below and above it, plus a
%% sensitivity-only member far above that must not widen any band or
%% range, on the fixture's own four-day axis.
[uncertainty_file,U] = write_uncertainty_fixture(out_dir,RM,n_time,0);
assert_error(@() run_plot_script([],out_dir,out_dir,1,'off', ...
    fullfile(out_dir,'absent_summary.mat'),''), ...
    'Plot_1D_Then_All_12:MissingUncertainty');
shifted_file = write_uncertainty_fixture(out_dir,RM,n_time,1);
assert_error(@() run_plot_script([],out_dir,out_dir,1,'off',shifted_file,''), ...
    'Plot_1D_Then_All_12:UncertaintyAxis');
close(findall(groot,'Type','figure'));

[fig1d,~,stats_1d,~,fig,plot_file,time_fig,time_plot_file] = ...
    run_plot_script([],out_dir,out_dir,1,'off',uncertainty_file,'');
%% Everything the default run asserts still holds with the bands drawn.
check_group_figure(fig1d,stats_1d);
check_figure_set(fig,plot_file,time_fig,time_plot_file,v2);

group_tags = {'CaSolution','MgSolution','SiSolution','CumulativeGaseousNLoss','WFPS'};
for k = 1:numel(group_tags)
    check_bands(fig1d,group_tags{k},group_tags{k},U,treatment_colors);
end
layout = findall(fig1d,'Type','tiledlayout');
assert(isscalar(layout) && contains(string(layout.Subtitle.String), ...
    'parameter-uncertainty range (3 members)'), ...
    'The group figure must say what the shading is.');
time_tags = {'CumulativeN2O','CumulativeNO3Leaching','NO3Pool','BulkPH','NEE','DailyN2O'};
time_metrics = {'CumulativeN2O','CumulativeNO3Leaching','NO3Pool','BulkPH','NEE','DailyN2O_gN_ha_d'};
for k = 1:6
    check_bands(time_fig(k),time_tags{k},time_metrics{k},U,treatment_colors);
    notes = findall(time_fig(k),'Type','text');
    notes = notes(arrayfun(@(h) any(contains(string(h.String), ...
        'Shaded: parameter-uncertainty range (3 members)')),notes));
    metric_lines = findall(time_fig(k),'Type','line','Tag',time_tags{k});
    assert(isscalar(notes) && ~strcmp(notes.FontWeight,'bold') && ...
        isequal(notes.Parent,ancestor(metric_lines(1),'axes')), ...
        'Each time-series figure carries one plain note on the shading under the top panel title.');
end
check_bands(time_fig(4),'PoreWaterPH','PoreWaterPH',U,treatment_colors);
%% Bands are only ever drawn under the 1-D configuration.
for k = 1:6
    n_band = numel(findall(time_fig(k),'Type','patch','Tag','UncertaintyBand'));
    assert(n_band == 3+3*double(k == 4), ...
        'Only the three 1-D treatments (and the pore-water panel) carry bands.');
end

range_fields = {'CumulativeN2O_gN_m2','CumulativeNO3Leaching_gN_m2', ...
    'FinalNO3Pool_gN_m2','FinalBulkPH','CumulativeNEE_gC_m2'};
for k = 1:5
    ranges = findall(fig(k),'Type','errorbar','Tag','UncertaintyRange');
    assert(numel(ranges) == 3,'The 1-D bars carry one range bar per treatment.');
    [~,order] = sort([ranges.XData]);
    ranges = ranges(order);
    for treatment = 1:3
        r = find(U.run_ids == 9+treatment,1);
        expected_min = U.summary_range.(range_fields{k}).min(r);
        expected_max = U.summary_range.(range_fields{k}).max(r);
        assert(abs(ranges(treatment).YData-ranges(treatment).YNegativeDelta-expected_min) < 1e-9 && ...
            abs(ranges(treatment).YData+ranges(treatment).YPositiveDelta-expected_max) < 1e-9, ...
            'Range bars must span the ensemble min and max.');
        assert(isequal(ranges(treatment).Color,[0 0 0]) && ...
            strcmp(ranges(treatment).HandleVisibility,'off'), ...
            'Range bars are thin black lines outside the legend.');
        %% Each bar series holds one configuration per element, the 1-D
        %% configuration first (UserData lists its run ids in that order).
        bar_series = findall(fig(k),'Type','bar');
        one_d_bar = bar_series(arrayfun(@(h) h.UserData(1) == 9+treatment,bar_series));
        assert(isscalar(one_d_bar) && ...
            abs(one_d_bar.XEndPoints(1)-ranges(treatment).XData) < 1e-9, ...
            'Range bars sit on the 1-D bars.');
    end
    notes = findall(fig(k),'Type','text');
    assert(nnz(arrayfun(@(h) any(contains(string(h.String), ...
        'Range bars: parameter-uncertainty range (3 members)')),notes)) == 1, ...
        'Each summary figure carries one note on the range bars.');
end
%% The fixture's four November 2016 days hold no chamber observation day,
%% so the observation-day N2O mean is NaN for every member and gets no bar.
assert(isempty(findall(fig(6),'Type','errorbar','Tag','UncertaintyRange')), ...
    'A NaN observation-day mean must not draw a range bar.');

%% Observed pH points carry the SEM of the block means, in their own axes.
sem_file = fullfile(v2,'Experimental observations','Web plot digitizer', ...
    'Beerling2024_Fig1E_vector_extraction.csv');
opts = detectImportOptions(sem_file,'CommentStyle','#','Delimiter',',');
opts = setvartype(opts,{'treatment','depth_cm','period'},'string');
T = readtable(sem_file,opts);
periods = ["pre" "2017" "2018" "2019" "2020"];
series_list = {'control','basalt','control','basalt'; '0-10','0-10','10-30','10-30'; ...
    'control_pH_0_10','basalt_pH_0_10','control_pH_10_30','basalt_pH_10_30'};
for series = series_list
    point = findall(time_fig(4),'Type','scatter','Tag',['Observed_' series{3}]);
    bar_handle = findall(time_fig(4),'Type','errorbar','Tag',['Observed_' series{3} '_SEM']);
    assert(isscalar(bar_handle) && isequal(ancestor(bar_handle,'axes'),ancestor(point,'axes')) && ...
        isequal(bar_handle.XData(:),point.XData(:)) && strcmp(bar_handle.HandleVisibility,'off'), ...
        'Each observed pH series carries one SEM bar object in the observation panel.');
    rows = T.treatment == series{1} & T.depth_cm == series{2};
    expected_sem = nan(5,1);
    for j = 1:5
        expected_sem(j) = T.SEM(rows & T.period == periods(j));
    end
    assert(isequaln(bar_handle.YNegativeDelta(:),expected_sem) && ...
        isequaln(bar_handle.YPositiveDelta(:),expected_sem), ...
        'pH error bars must be the digitised SEM, absent where the figure draws none.');
end
close([fig1d; fig(:); time_fig(:)]);

%% Scenario results carry the scenario tag after the run name.
scenario_dir = fullfile(out_dir,'ssp245');
mkdir(scenario_dir);
for k = 10:12
    write_fixture(fullfile(scenario_dir,sprintf('TC_2D_run%02d_%s__ssp245__1d.mat', ...
        k,RM(k).name)),fixture_series(k,RM(k),n_time));
end
assert_error(@() run_plot_script([],scenario_dir,scenario_dir,1,'off','',''), ...
    'Plot_1D_Then_All_12:MissingResults');
assert_error(@() run_plot_script([],out_dir,out_dir,1,'off','','ssp245'), ...
    'Plot_1D_Then_All_12:MissingResults');
[fig1d,plot_file_1d,stats_1d,summary] = ...
    run_plot_script([],scenario_dir,scenario_dir,1,'off','','ssp245');
check_group_figure(fig1d,stats_1d);
assert(isempty(summary),'Stage 2 is skipped while the 2-D scenario runs are missing.');
assert(contains(plot_file_1d,'__ssp245__1d') && isfile(plot_file_1d), ...
    'Scenario figures carry the scenario tag.');
[~,~,~,summary,fig,plot_file,time_fig] = ...
    run_plot_script(10:12,scenario_dir,scenario_dir,1,'off','','ssp245');
assert(isequal(summary.RunID.',10:12) && all(cellfun(@isfile,plot_file)) && ...
    all(contains(plot_file,'__ssp245__1d')), ...
    'Scenario summary figures are found by the tagged result names.');
bars = findall(fig(1),'Type','bar','Tag','CumulativeN2O');
assert(numel(bars) == 3 && isequal(sort([bars.UserData]),10:12), ...
    'Scenario results plot like any other result set.');
close([fig1d; fig(:); time_fig(:)]);
fprintf('PASS: all treatments and configurations use one figure per metric.\n');
end

function [fig1d,plot_file_1d,stats_1d,summary,fig,plot_file,time_fig,time_plot_file] = ...
        run_plot_script(Runs,ResultDirectory,FigureDirectory,NDays,Visible,UncertaintyFile,Scenario)
%% Run the script with these workspace options and hand back its results.
fprintf(['Plot_1D_Then_All_12 with Runs %s, results %s, figures %s, NDays %d, ' ...
    'Visible %s, UncertaintyFile ''%s'', Scenario ''%s''\n'],mat2str(Runs), ...
    ResultDirectory,FigureDirectory,NDays,Visible,UncertaintyFile,Scenario);
summary = []; fig = []; plot_file = []; time_fig = []; time_plot_file = [];
fig1d = []; plot_file_1d = []; stats_1d = [];
Plot_1D_Then_All_12;
end

function F = fixture_series(k,R,n_time)
%% One run's saved state: every quantity the script reads, scaled by k.
F.R = R;
F.Domain.N2Oflx = k*[0; 1; 2; 3];
F.Domain.LEAK_NO3 = k*[0; 2; 4; 6];
F.Domain.NO3_pool = 10+k*[0; 1; 2; 3];
F.Domain.PH_d = 6.1*ones(n_time,1);
F.Domain.PH_bulk_d = F.Domain.PH_d;
F.Domain.PH_pore_d = F.Domain.PH_d+0.1;
F.PH_d = F.Domain.PH_d;
F.N2Oflx = F.Domain.N2Oflx;
F.N2Oflx_nit = 0.25*F.N2Oflx;
F.N2flx = k*[0; 2; 4; 6];
if k == 10
    F.N2Oflx_nit = F.N2Oflx;
    F.N2flx = zeros(n_time,1);
end
F.Datam = [2016 11 8 0];
F.NPP_H = reshape(k*[1; 2; 3; 4],n_time,1,1);
F.NPP_L = reshape(k*[0.5; 1; 1.5; 2],n_time,1,1);
F.Ccrown_t = ones(n_time,1,1);
F.R_litter = reshape(k*[0.1; 0.2; 0.3; 0.4],n_time,1,1);
F.R_microbe = reshape(k*[0.2; 0.4; 0.6; 0.8],n_time,1,1);
F.R_ew = zeros(n_time,1,1);
F.RB_L = zeros(n_time,1,7,1);
F.P = zeros(n_time,76);
F.P(:,60) = k+[0; 1; 2; 3];
F.P(:,68) = k+[0; 2; 4; 6];
F.P(:,76) = k+[0; 3; 6; 9];
F.O = repelem([0.1; 0.2; 0.3],24);
F.p.Bio_Zs = 1;
F.p.Soil_Param.Osat = 0.5;
F.topo.cell_area_m2 = 3600;
F.ColCfg.PFT_OVERRIDE = '';
end

function write_fixture(result_file,F)
save(result_file,'-struct','F');
end

function d = fixture_daily(F)
%% The daily quantities Plot_1D_Then_All_12 derives from one saved run,
%% in the summary file's names and units.
d.CaSolution = F.P(:,60); d.MgSolution = F.P(:,68); d.SiSolution = F.P(:,76);
d.CumulativeGaseousNLoss = cumsum(F.N2flx+F.N2Oflx);
d.WFPS = [NaN; 20; 40; 60];
d.CumulativeN2O = cumsum(F.Domain.N2Oflx);
d.CumulativeNO3Leaching = cumsum(F.Domain.LEAK_NO3);
d.NO3Pool = F.Domain.NO3_pool;
d.BulkPH = F.Domain.PH_bulk_d;
d.PoreWaterPH = F.Domain.PH_pore_d;
d.NEE = F.R_litter(:)+F.R_microbe(:)+F.R_ew(:)-(F.NPP_H(:)+F.NPP_L(:));
d.DailyN2O_gN_ha_d = 1e4*F.Domain.N2Oflx;
end

function [uncertainty_file,U] = write_uncertainty_fixture(out_dir,RM,n_time,day_shift)
%% Summarize_Uncertainty's contract on the fixture: member 1 is
%% the fixture itself, members 2 and 3 sit a spread below and above it,
%% member 4 is a sensitivity-only member three spreads above, outside the
%% envelope.
run_ids = [10 11 12];
n_members = 4;
in_envelope = [true true true false];
U = struct('Date',convertTo(datetime(2016,11,8)+caldays((0:n_time-1).'-2+day_shift),'datenum'), ...
    'run_ids',run_ids,'treatments',{{'control','basalt','calcite'}}, ...
    'in_envelope',in_envelope);
for r = 1:numel(run_ids)
    d = fixture_daily(fixture_series(run_ids(r),RM(run_ids(r)),n_time));
    for f = fieldnames(d).'
        base = d.(f{1});
        spread = 0.1*abs(base)+0.05;
        if r == 1; U.daily.(f{1}) = nan(n_time,n_members,numel(run_ids)); end
        U.daily.(f{1})(:,:,r) = [base base-spread base+spread base+3*spread];
    end
end
for f = fieldnames(U.daily).'
    x = U.daily.(f{1})(:,in_envelope,:);
    U.envelope.(f{1}) = struct('min',reshape(min(x,[],2),n_time,[]), ...
        'max',reshape(max(x,[],2),n_time,[]),'base',reshape(x(:,1,:),n_time,[]));
end
value = struct('CumulativeN2O_gN_m2',squeeze(U.daily.CumulativeN2O(end,:,:)), ...
    'CumulativeNO3Leaching_gN_m2',squeeze(U.daily.CumulativeNO3Leaching(end,:,:)), ...
    'FinalNO3Pool_gN_m2',squeeze(U.daily.NO3Pool(end,:,:)), ...
    'FinalBulkPH',squeeze(U.daily.BulkPH(end,:,:)), ...
    'CumulativeNEE_gC_m2',squeeze(sum(U.daily.NEE,1)));
for f = fieldnames(value).'
    x = value.(f{1});
    U.summary_range.(f{1}) = struct('min',min(x(in_envelope,:),[],1), ...
        'max',max(x(in_envelope,:),[],1),'base',x(1,:),'members',x);
end
U.summary_value = value;
U.obs = struct('years',2016:2020);
U.obs.pH.control_0_10 = struct('value',[6.2 6.28 6.31 6.13 5.98],'sem',0.06*ones(1,5));
U.obs.pH.basalt_0_10 = struct('value',[6.13 6.6 6.66 6.92 6.92],'sem',[0.07 0.06 0.04 NaN 0.03]);
U.obs.CDR_potential_tCO2_ha = struct('value',[0 3.8 6.8 11.5 10.5],'sem',[0 1.8 2.5 3.3 3.7]);
uncertainty_file = fullfile(out_dir,sprintf('uncertainty_summary_shift%d.mat',day_shift));
save(uncertainty_file,'-struct','U');
end

function check_bands(fig_handle,line_tag,metric,U,treatment_colors)
%% One translucent band per 1-D treatment beneath the lines of one panel,
%% its vertices the ensemble minimum and maximum of that treatment's run.
lines = findall(fig_handle,'Type','line','Tag',line_tag);
ax = ancestor(lines(1),'axes');
bands = findall(ax,'Type','patch','Tag','UncertaintyBand');
assert(numel(bands) == 3,'Expected one band per 1-D treatment under %s.',line_tag);
for treatment = 1:3
    band = bands(arrayfun(@(h) isequal(h.FaceColor,treatment_colors(treatment,:)),bands));
    assert(isscalar(band),'Bands take the treatment colour under %s.',line_tag);
    r = find(U.run_ids == 9+treatment,1);
    low = U.envelope.(metric).min(:,r); high = U.envelope.(metric).max(:,r);
    ok = isfinite(low) & isfinite(high);
    n = nnz(ok);
    y = band.YData(:);
    assert(numel(y) == 2*n && max(abs(y(1:n)-low(ok))) < 1e-9 && ...
        max(abs(flipud(y(n+1:end))-high(ok))) < 1e-9, ...
        'Band vertices must span the ensemble min and max under %s.',line_tag);
    assert(abs(band.FaceAlpha-0.18) < 1e-12 && strcmp(band.EdgeColor,'none') && ...
        strcmp(band.HandleVisibility,'off'), ...
        'Bands are translucent, edgeless and outside the legend.');
    %% allchild lists the most recently drawn object first, the hidden-handle
    %% bands included.
    stack = allchild(ax);
    band_index = find(stack == band,1);
    line_index = arrayfun(@(l) find(stack == l,1),lines);
    assert(band_index > max(line_index), ...
        'Bands are drawn beneath the lines under %s.',line_tag);
end
end

function check_group_figure(fig1d,stats_1d)
wfps_lines = findall(fig1d,'Type','line','Tag','WFPS');
assert(numel(wfps_lines) == 3, ...
    'The 1-D WFPS panel must contain three treatment traces.');
for k = 1:numel(wfps_lines)
    assert(isnan(wfps_lines(k).YData(1)) && max(abs(wfps_lines(k).YData- ...
        [NaN 20 (23*20+40)/24 (23*40+60)/24]),[],'omitnan') < 1e-12, ...
        'WFPS must average the same hours as the BIOGEO row it is stored in.');
end

expected_rN2O = 100*(0.75*3)/(6+0.75*3);
assert(isnan(stats_1d.control.rN2O) && ...
    abs(stats_1d.basalt.rN2O-expected_rN2O) < 1e-12 && ...
    abs(stats_1d.calcite.rN2O-expected_rN2O) < 1e-12, ...
    'The exact denitrification N2O share must be NaN when its denominator is zero.');

expected_cation_units = {'g Ca m^{-2}','g Mg m^{-2}','g Si m^{-2}'};
cation_tags = {'CaSolution','MgSolution','SiSolution'};
for measurement = 1:3
    plotted_line = findall(fig1d,'Type','line','Tag',cation_tags{measurement});
    assert(numel(plotted_line) == 3, ...
        'Each weathering metric must contain three treatment traces.');
    ax = ancestor(plotted_line(1),'axes');
    actual_unit = string(ax.YLabel.String);
    assert(actual_unit == expected_cation_units{measurement}, ...
        'Weathering-panel units must match the saved elemental pools. Expected %s, got %s.', ...
        expected_cation_units{measurement},actual_unit);
end
gas_line = findall(fig1d,'Type','line','Tag','CumulativeGaseousNLoss');
assert(numel(gas_line) == 3,'The gaseous-N panel must contain three traces.');
gas_axis = ancestor(gas_line(1),'axes');
gas_labels = string({findall(gas_axis,'Type','text').String});
assert(any(gas_labels == "Cumulative gaseous N loss (N_2 + N_2O)"), ...
    'The N2 + total N2O panel must not be labelled as pure denitrification.');
bottom_tags = {'CumulativeGaseousNLoss','WFPS'};
for measurement = 1:2
    bottom_line = findall(fig1d,'Type','line','Tag',bottom_tags{measurement});
    bottom_axis = ancestor(bottom_line(1),'axes');
    crop_labels = findall(bottom_axis,'Type','text','Tag','KantolaCropLabel');
    panel_title = setdiff(findall(bottom_axis,'Type','text'),crop_labels);
    panel_title = panel_title(arrayfun(@(h) strcmp(h.FontWeight,'bold'),panel_title));
    assert(isscalar(panel_title) && ~isempty(crop_labels) && ...
        panel_title.Position(2) < min(arrayfun(@(h) h.Position(2),crop_labels)), ...
        'Bottom-panel titles must sit below the crop-year labels.');
end
end

function check_figure_set(fig,plot_file,time_fig,time_plot_file,v2)
assert(numel(fig) == 6 && all(isgraphics(fig,'figure')), ...
    'Each summary metric must have its own figure.');
assert(numel(time_fig) == 6 && all(isgraphics(time_fig,'figure')), ...
    'Each time-series metric must have its own figure.');
assert(numel(plot_file) == 6 && all(cellfun(@isfile,plot_file)), ...
    'Each summary metric must have its own PNG.');
assert(numel(time_plot_file) == 6 && all(cellfun(@isfile,time_plot_file)), ...
    'Each time-series metric must have its own PNG.');
assert(all(cellfun(@companion_fig_exists,plot_file)) && ...
    all(cellfun(@companion_fig_exists,time_plot_file)), ...
    'Each metric PNG must have an editable FIG copy.');

expected_ids = [10 11 12 1 2 3 4 5 6 7 8 9];
expected_treatment_labels = {'Control','Basalt','Calcite'};
expected_configuration_labels = {'1-D','2-D whole field', ...
    '2-D buffer only','2-D native vegetation buffer'};
for k = 1:6
    assert(isscalar(findall(fig(k),'Type','axes')), ...
        'Summary metric figures must not contain subplots.');
    %% Metric 4 is pH: this fixture gives every run both bulk and pore
    %% series, so it draws as two stacked panels, plus a third panel
    %% comparing modelled years against the digitised Beerling pH
    %% observations (real data on disk; bulk_ph_obs.available does not
    %% depend on this fixture's tiny date range). Metric 6 would gain a
    %% second panel carrying the dated chamber N2O record, but
    %% Datam = [2016 11 8 0] with
    %% n_time = 4 gives a Date_d entirely inside 2016, which can never
    %% overlap the digitised chamber record's 2017/2018 dates -- the
    %% n2o_observations/clip_observations call inside plot_subfield is
    %% therefore always empty here, deterministically, independent of the
    %% fixture's actual flux values. That show_n2o_obs decision is local
    %% to plot_subfield (a local function in the file, not the top-level
    %% script Plot_1D_Then_All_12 itself), so this test cannot read the
    %% variable directly; the date-range argument above is what stands in
    %% its place. Every metric also carries one invisible overlay axes
    %% hosting the configuration (line-style) legend, alongside the
    %% treatment (colour) legend on the metric's own top data axes -- two
    %% short legends replaced the one that tried to spell out all twelve
    %% combinations and overlapped.
    expected_time_axes = 1+double(k == 4)*2+1;
    assert(numel(findall(time_fig(k),'Type','axes')) == expected_time_axes, ...
        'Time-series metric figures must not contain unexpected subplots.');
    summary_lgd = findall(fig(k),'Type','legend');
    assert(isscalar(summary_lgd) && ...
        isequal(summary_lgd.String,{'Control','Basalt','Calcite'}), ...
        'Summary bars must contain all three feedstock treatments.');
    lgd = findall(time_fig(k),'Type','legend');
    assert(numel(lgd) == 2, ...
        'Each time-series figure must carry a treatment legend and a configuration legend.');
    is_treatment = arrayfun(@(l) isequal(l.String,expected_treatment_labels),lgd);
    is_configuration = arrayfun(@(l) isequal(l.String,expected_configuration_labels),lgd);
    assert(nnz(is_treatment) == 1 && nnz(is_configuration) == 1, ...
        'The two legends must preserve the requested treatment and configuration order.');
end

bars = findall(fig(1),'Type','bar','Tag','CumulativeN2O');
assert(numel(bars) == 3 && ...
    isequal(sort([bars.UserData]),sort(expected_ids)), ...
    'Every summary figure must contain all twelve treatment configurations.');

lines = findall(time_fig(1),'Type','line','Tag','CumulativeN2O');
assert(numel(lines) == 12 && isequal(sort([lines.UserData]),sort(expected_ids)), ...
    'Every time-series figure must contain all twelve treatment configurations.');
for run_id = expected_ids
    line = lines([lines.UserData] == run_id);
    assert(isequal(line.YData,run_id*[0 1 3 6]), ...
        'Configuration traces must use the matching run result.');
end
assert(isequal(lines([lines.UserData] == 10).Color,[0 0 0]) && ...
    isequal(lines([lines.UserData] == 11).Color,[0 0.25 0.60]) && ...
    isequal(lines([lines.UserData] == 12).Color,[0.35 0.70 0.90]), ...
    'Control, basalt and calcite must retain their treatment colours.');
assert(strcmp(lines([lines.UserData] == 10).LineStyle,'-') && ...
    strcmp(lines([lines.UserData] == 1).LineStyle,'--') && ...
    strcmp(lines([lines.UserData] == 4).LineStyle,':') && ...
    strcmp(lines([lines.UserData] == 7).LineStyle,'-.'), ...
    'Line style must identify the 1-D and three 2-D configurations.');

pool_lines = findall(time_fig(3),'Type','line','Tag','NO3Pool');
assert(numel(pool_lines) == 12, ...
    'The soil nitrate pool figure must remain in the comparison set.');
pore_lines = findall(time_fig(4),'Type','line','Tag','PoreWaterPH');
assert(numel(pore_lines) == 12 && ...
    isequal(sort([pore_lines.UserData]),sort(expected_ids)), ...
    'The pH figure must retain pore-water pH for every treatment and configuration.');
bulk_lines = findall(time_fig(4),'Type','line','Tag','BulkPH');
bulk_axis = ancestor(bulk_lines(1),'axes');
pore_axis = ancestor(pore_lines(1),'axes');
assert(~isequal(bulk_axis,pore_axis), ...
    'Bulk and pore-water pH must be drawn in separate panels.');
assert(diff(bulk_axis.YLim) >= 0.1-1e-12 && diff(pore_axis.YLim) >= 0.1-1e-12, ...
    'A constant pH series must receive a readable axis range in both panels.');
%% bulk_pH_observations is a local function inside Plot_1D_Then_All_12.m;
%% local functions in a script file are not callable from another file, so
%% this test reads the same real data file directly instead, to check the
%% panel this test already exercised against real data.
beerling_file = fullfile(v2,'Experimental observations', ...
    'Web plot digitizer','Data_Beerling.mat');
S_beerling = load(beerling_file,'pH_Exp');
observed_pH = struct();
for tag = {'control_pH_0_10','control_pH_10_30','basalt_pH_0_10','basalt_pH_10_30'}
    h = findall(time_fig(4),'Type','scatter','Tag',['Observed_' tag{1}]);
    assert(isscalar(h) && numel(h.XData) == numel(S_beerling.pH_Exp.X), ...
        'The pH figure must carry one observed point per label for %s.',tag{1});
    observed_pH.(tag{1}) = h;
end
assert(isequal(observed_pH.control_pH_0_10.YData(:),double(S_beerling.pH_Exp.control_0(:))) && ...
       isequal(observed_pH.basalt_pH_10_30.YData(:),double(S_beerling.pH_Exp.basalt_10(:))), ...
    'The observed pH points must be the measured depth intervals, not a depth average.');
assert(isempty(findall(time_fig(4),'Type','scatter','Tag','Observed_control_pH_0_30')), ...
    'The depth-averaged pH observation series is retired.');
observed_control_pH = observed_pH.control_pH_0_10;
assert(~isequal(ancestor(observed_control_pH,'axes'),bulk_axis) && ...
    ~isequal(ancestor(observed_control_pH,'axes'),pore_axis), ...
    'The yearly pH observations must sit in their own panel, not the bulk or pore time series.');
end

function tf = companion_fig_exists(png_file)
[folder,name] = fileparts(png_file);
fig_file = fullfile(folder,[name '.fig']);
tf = isfile(fig_file);
if tf
    f = openfig(fig_file,'invisible');
    tf = isgraphics(f,'figure');
    close(f);
end
end

function tidy(out_dir)
close(findall(groot,'Type','figure'));
if isfolder(out_dir); rmdir(out_dir,'s'); end
end

function assert_error(action,identifier)
try
    action();
catch ME
    assert(strcmp(ME.identifier,identifier), ...
        'Expected %s but received %s.',identifier,ME.identifier);
    return
end
error('Expected error %s was not raised.',identifier);
end
