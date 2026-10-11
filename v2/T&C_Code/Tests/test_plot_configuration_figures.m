function test_plot_configuration_figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Configuration figures contain every treatment in every run group.
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
assert(isequaln(wfps_test.',[NaN 20 40 60]), ...
    'WFPS must align each 24-hour mean with the following retained state row.');
assert_error(@() Daily_WFPS_From_Hourly(zeros(71,1), ...
    struct('Bio_Zs',1,'Soil_Param',struct('Osat',0.5)),n_time), ...
    'Plot_1D_Then_All_12:WFPSLength');

for k = 1:numel(RM)
    R = RM(k);
    Domain.N2Oflx = k*[0; 1; 2; 3];
    Domain.LEAK_NO3 = k*[0; 2; 4; 6];
    Domain.NO3_pool = 10+k*[0; 1; 2; 3];
    Domain.PH_d = 6.1*ones(n_time,1);
    Domain.PH_bulk_d = Domain.PH_d;
    Domain.PH_pore_d = Domain.PH_d+0.1;
    PH_d = Domain.PH_d;
    N2Oflx = Domain.N2Oflx;
    N2Oflx_nit = 0.25*N2Oflx;
    N2flx = k*[0; 2; 4; 6];
    if k == 10
        N2Oflx_nit = N2Oflx;
        N2flx = zeros(n_time,1);
    end
    Datam = [2016 11 8 0];
    NPP_H = reshape(k*[1; 2; 3; 4],n_time,1,1);
    NPP_L = reshape(k*[0.5; 1; 1.5; 2],n_time,1,1);
    Ccrown_t = ones(n_time,1,1);
    R_litter = reshape(k*[0.1; 0.2; 0.3; 0.4],n_time,1,1);
    R_microbe = reshape(k*[0.2; 0.4; 0.6; 0.8],n_time,1,1);
    R_ew = zeros(n_time,1,1);
    RB_L = zeros(n_time,1,7,1);
    P = zeros(n_time,76);
    P(:,60) = k+[0; 1; 2; 3];
    P(:,68) = k+[0; 2; 4; 6];
    P(:,76) = k+[0; 3; 6; 9];
    O = repelem([0.1; 0.2; 0.3],24);
    p.Bio_Zs = 1;
    p.Soil_Param.Osat = 0.5;
    topo.cell_area_m2 = 3600;
    ColCfg.PFT_OVERRIDE = '';
    save(fullfile(out_dir,sprintf('TC_2D_run%02d_%s__1d.mat', ...
        k,R.name)),'R','Domain','Datam','NPP_H','NPP_L', ...
        'Ccrown_t','R_litter','R_microbe','R_ew','RB_L','P','O','p', ...
        'PH_d','N2Oflx','N2Oflx_nit','N2flx', ...
        'topo','ColCfg');
end

summary = []; fig = []; plot_file = []; time_fig = []; time_plot_file = []; %#ok<NASGU>
fig1d = []; plot_file_1d = []; stats_1d = []; %#ok<NASGU>
Runs = 10:12; %#ok<NASGU>
ResultDirectory = out_dir; %#ok<NASGU>
FigureDirectory = out_dir; %#ok<NASGU>
NDays = 1; %#ok<NASGU>
Visible = 'off'; %#ok<NASGU>
Plot_1D_Then_All_12;

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

Runs = []; %#ok<NASGU>
Plot_1D_Then_All_12;

wfps_lines = findall(fig1d,'Type','line','Tag','WFPS');
assert(numel(wfps_lines) == 3, ...
    'The 1-D WFPS panel must contain three treatment traces.');
for k = 1:numel(wfps_lines)
    assert(isequaln(wfps_lines(k).YData,[NaN 20 40 60]), ...
        'WFPS must align each 24-hour mean with the following retained state row.');
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
observed_control_pH = findall(time_fig(4),'Type','scatter', ...
    'Tag','Observed_control_pH_0_30');
observed_basalt_pH = findall(time_fig(4),'Type','scatter', ...
    'Tag','Observed_basalt_pH_0_30');
assert(isscalar(observed_control_pH) && isscalar(observed_basalt_pH) && ...
    numel(observed_control_pH.XData) == numel(S_beerling.pH_Exp.X), ...
    'The pH figure must carry one yearly observed point per treatment per label.');
assert(~isequal(ancestor(observed_control_pH,'axes'),bulk_axis) && ...
    ~isequal(ancestor(observed_control_pH,'axes'),pore_axis), ...
    'The yearly pH observations must sit in their own panel, not the bulk or pore time series.');

close([fig1d; fig(:); time_fig(:)]);
fprintf('PASS: all treatments and configurations use one figure per metric.\n');
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
