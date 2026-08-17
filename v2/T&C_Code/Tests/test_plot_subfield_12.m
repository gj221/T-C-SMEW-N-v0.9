function test_plot_subfield_12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Twelve-scenario summary and time-dynamics plot contract.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n=== test_plot_subfield_12 ===\n');
n_pass = 0;
%% Run_Matrix_12 and Plot_Subfield_12 live in v2, not in T&C_Code; without
%% this the test only passes when the caller already has v2 on the path.
v2 = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(v2,fullfile(v2,'T&C_Code'));
out_dir = tempname;
mkdir(out_dir);
subset_dir = tempname;
mkdir(subset_dir);
two_d_dir = tempname;
mkdir(two_d_dir);
cleanup = onCleanup(@() tidy(out_dir,subset_dir,two_d_dir));

RM = Run_Matrix_12();
expected_nee = zeros(numel(RM),4);
expected_crop_yield = zeros(numel(RM),4);
for k = 1:numel(RM)
    R = RM(k);
    Domain.N2Oflx = k*[0; 1; 2; 3];
    Domain.LEAK_NO3 = k*[0; 2; 4; 6];
    Domain.NO3_pool = 10+k*[0; 1; 2; 3];
    Domain.PH_d = 6+0.01*k*[0; 1; 2; 3];
    Datam = [2016 11 8 0];
    [NPP_H,NPP_L,Ccrown_t,R_litter,R_microbe,R_ew,RB_L,topo,ColCfg, ...
        expected_nee(k,:),expected_crop_yield(k,:)] = carbon_fixture(k,R);
    save(fullfile(out_dir,sprintf('TC_2D_run%02d_%s__1d.mat', ...
        k,R.name)),'R','Domain','Datam','NPP_H','NPP_L', ...
        'Ccrown_t','R_litter','R_microbe','R_ew','RB_L','topo','ColCfg');
end

[summary,fig,plot_file,time_fig,time_plot_file] = Plot_Subfield_12( ...
    'ResultDirectory',out_dir,'NDays',1,'Visible','off');
check(height(summary) == 12, 'plot summary contains all twelve scenarios');
check(isequal(summary.RunID.',1:12), 'plot summary preserves matrix order');
check(isequal(summary.CumulativeN2O_gN_m2.',6*(1:12)), ...
    'cumulative N2O values are read from each result');
check(isequal(summary.CumulativeNO3Leaching_gN_m2.',12*(1:12)), ...
    'cumulative nitrate leaching values are read from each result');
check(ismember('FinalPoreWaterPH',summary.Properties.VariableNames), ...
    'summary identifies PH_d as pore-water pH');
check(max(abs(summary.CumulativeNEE_gC_m2-sum(expected_nee,2))) < 1e-12, ...
    'summary area-weights cumulative NEE');
check(~ismember('CumulativeVegetationNPP_gC_m2', ...
    summary.Properties.VariableNames), ...
    'summary replaces vegetation NPP with NEE');
check(max(abs(summary.CumulativeCropYield_gC_m2- ...
    sum(expected_crop_yield,2))) < 1e-12, ...
    'summary area-weights crop yield');
check(all(summary.CumulativeNEE_gC_m2(7:9)./(7:9).' < ...
    summary.CumulativeNEE_gC_m2(4:6)./(4:6).'), ...
    'prairie PFT uptake contributes to whole-domain NEE');
check(all(summary.CumulativeCropYield_gC_m2(7:9)./(7:9).' < ...
    summary.CumulativeCropYield_gC_m2(4:6)./(4:6).'), ...
    'prairie PFTs contribute no crop yield');
check(isgraphics(fig,'figure') && isfile(plot_file), ...
    'plot is created and written to disk');
check(companion_fig_exists(plot_file), ...
    'summary plot is also written as a loadable FIG');
check(panel_has_values(fig,'Cumulative N_2O',6*(1:12)), ...
    'N2O panel plots all twelve scenarios');
check(panel_has_values(fig,'Cumulative NO_3 leaching',12*(1:12)), ...
    'nitrate-leaching panel plots all twelve scenarios');
check(panel_has_values(fig,'Final soil NO_3 pool',10+3*(1:12)), ...
    'final-NO3 panel plots all twelve scenarios');
check(panel_has_values(fig,'Final pore-water pH',6+0.03*(1:12)), ...
    'pore-water-pH panel plots all twelve scenarios');
check(panel_has_values(fig,'Cumulative NEE (+ to atmosphere)', ...
    sum(expected_nee,2)), ...
    'NEE panel plots all twelve scenarios');
check(panel_has_values(fig,'Cumulative crop yield', ...
    sum(expected_crop_yield,2)), ...
    'crop-yield panel plots all twelve scenarios');
expected_colors = [0 0 0; 0 0.25 0.60; 0.35 0.70 0.90];
check(object_colors_match(fig,'bar',expected_colors,6), ...
    'summary bars use black, dark blue and light blue treatments');
check(clean_axes_and_bottom_legend(fig,false,4), ...
    'summary plot has no grid and uses a bottom legend with black text');
check(isgraphics(time_fig,'figure') && isfile(time_plot_file), ...
    'time-dynamics plot is created and written to disk');
check(companion_fig_exists(time_plot_file), ...
    'time-dynamics plot is also written as a loadable FIG');
check(time_lines_match(time_fig,'CumulativeN2O', ...
    (1:12).'*[0 1 3 6]), ...
    'cumulative N2O dynamics contain all twelve scenarios');
check(time_lines_match(time_fig,'CumulativeNO3Leaching', ...
    (1:12).'*[0 2 6 12]), ...
    'cumulative nitrate-leaching dynamics contain all twelve scenarios');
check(time_lines_match(time_fig,'NO3Pool', ...
    10+(1:12).'*[0 1 2 3]), ...
    'soil-NO3 dynamics contain all twelve scenarios');
check(time_lines_match(time_fig,'PoreWaterPH', ...
    6+0.01*(1:12).'*[0 1 2 3]), ...
    'pore-water-pH dynamics contain all twelve scenarios');
check(time_lines_match(time_fig,'NEE',expected_nee), ...
    'daily NEE dynamics contain all twelve scenarios');
check(time_lines_match(time_fig,'CumulativeCropYield', ...
    cumsum(expected_crop_yield,2)), ...
    'cumulative crop-yield dynamics contain all twelve scenarios');
check(time_lines_have_daily_dates(time_fig,12), ...
    'time-dynamics traces use the saved daily date axis');
check(object_colors_match(time_fig,'line',expected_colors,24), ...
    'time traces use black, dark blue and light blue treatments');
check(feedstock_markers_match(time_fig,24), ...
    'each time panel marks the in-window feedstock application date');
check(clean_axes_and_bottom_legend(time_fig,true,4), ...
    'time plot uses at most three stacked rows with a bottom legend');
close([fig time_fig]);

for run_id = 1:9
    tag = sprintf('TC_2D_run%02d_%s__1d.mat',run_id,RM(run_id).name);
    copyfile(fullfile(out_dir,tag),fullfile(two_d_dir,tag));
end
[summary,fig,plot_file,time_fig,time_plot_file] = Plot_Subfield_12( ...
    'ResultDirectory',two_d_dir,'Runs',1:9,'NDays',1,'Visible','off');
check(isequal(summary.RunID.',1:9), ...
    'two-dimensional plot loads runs 1 to 9 without the 1-D results');
check(endsWith(plot_file,'Subfield_2D_summary__1d.png') && isfile(plot_file), ...
    'two-dimensional summary uses a distinct output filename');
check(endsWith(time_plot_file,'Subfield_2D_time_dynamics__1d.png') && ...
    isfile(time_plot_file), ...
    'two-dimensional dynamics use a distinct output filename');
check(companion_fig_exists(plot_file) && ...
    companion_fig_exists(time_plot_file), ...
    'two-dimensional summary and dynamics include FIG files');
check(time_lines_match(time_fig,'NEE',expected_nee(1:9,:)), ...
    'two-dimensional dynamics include area-weighted NEE');
check(clean_axes_and_bottom_legend(time_fig,true,3), ...
    'two-dimensional dynamics use no more than three stacked rows');
close([fig time_fig]);

for run_id = 10:12
    tag = sprintf('TC_2D_run%02d_%s__1d.mat',run_id,RM(run_id).name);
    copyfile(fullfile(out_dir,tag),fullfile(subset_dir,tag));
end
[summary,fig,plot_file,time_fig,time_plot_file] = Plot_Subfield_12( ...
    'ResultDirectory',subset_dir,'Runs',10:12,'NDays',1,'Visible','off');
check(isequal(summary.RunID.',10:12), ...
    'one-dimensional plot loads runs 10 to 12 without the 2-D results');
check(endsWith(plot_file,'Vertical_1D_summary__1d.png') && isfile(plot_file), ...
    'one-dimensional summary uses a distinct output filename');
check(endsWith(time_plot_file,'Vertical_1D_time_dynamics__1d.png') && ...
    isfile(time_plot_file), ...
    'one-dimensional dynamics use a distinct output filename');
check(companion_fig_exists(plot_file) && ...
    companion_fig_exists(time_plot_file), ...
    'one-dimensional summary and dynamics include FIG files');
check(panel_has_values(fig,'Cumulative NEE (+ to atmosphere)', ...
    sum(expected_nee(10:12,:),2)), ...
    'one-dimensional summary includes NEE');
check(panel_has_values(fig,'Cumulative crop yield', ...
    sum(expected_crop_yield(10:12,:),2)), ...
    'one-dimensional summary includes crop yield');
check(time_lines_match(time_fig,'NEE',expected_nee(10:12,:)), ...
    'one-dimensional dynamics include NEE');
check(time_lines_match(time_fig,'CumulativeCropYield', ...
    cumsum(expected_crop_yield(10:12,:),2)), ...
    'one-dimensional dynamics include cumulative crop yield');
check(object_colors_match(time_fig,'line',expected_colors,6), ...
    'one-dimensional traces retain the treatment colors');
check(feedstock_markers_match(time_fig,6), ...
    'one-dimensional subplots retain feedstock dates');
check(clean_axes_and_bottom_legend(time_fig,true,1), ...
    'one-dimensional dynamics use three rows and two columns');
close([fig time_fig]);

check(raises_error(@() Plot_Subfield_12( ...
    'ResultDirectory',subset_dir,'Runs',[10 11],'NDays',1,'Visible','off'), ...
    'Plot_Subfield_12:BadRuns'), ...
    'plot rejects an incomplete scenario treatment group');

bad_file = fullfile(out_dir,sprintf('TC_2D_run%02d_%s__1d.mat',12,RM(12).name));
S = load(bad_file);
S_missing = rmfield(S,'R_ew');
save(bad_file,'-struct','S_missing');
check(raises_error(@() Plot_Subfield_12( ...
    'ResultDirectory',out_dir,'NDays',1,'Visible','off'), ...
    'Plot_Subfield_12:IncompleteCarbonResult'), ...
    'plot rejects a result missing an NEE respiration component');
save(bad_file,'-struct','S');

S.Domain.PH_d(end) = NaN;
save(bad_file,'-struct','S');
check(raises_error(@() Plot_Subfield_12( ...
    'ResultDirectory',out_dir,'NDays',1,'Visible','off'), ...
    'Plot_Subfield_12:NonfiniteMeasurement'), ...
    'plot rejects a scenario with a nonfinite final measurement');
S.Domain.PH_d(end) = 6.36;
S.Domain.N2Oflx(:) = NaN;
save(bad_file,'-struct','S');
check(raises_error(@() Plot_Subfield_12( ...
    'ResultDirectory',out_dir,'NDays',1,'Visible','off'), ...
    'Plot_Subfield_12:NonfiniteMeasurement'), ...
    'plot rejects a scenario with no finite cumulative measurement');

S.Domain.N2Oflx = 12*[0; 1; 2; 3];
S.Domain.PH_d = S.Domain.PH_d(1:3);
save(bad_file,'-struct','S');
check(raises_error(@() Plot_Subfield_12( ...
    'ResultDirectory',out_dir,'NDays',1,'Visible','off'), ...
    'Plot_Subfield_12:TimeLengthMismatch'), ...
    'plot rejects measurements with mismatched daily lengths');

S.Domain.PH_d = 6+0.12*[0; 1; 2; 3];
S.RB_L = S.RB_L(1:3,:,:,:);
save(bad_file,'-struct','S');
check(raises_error(@() Plot_Subfield_12( ...
    'ResultDirectory',out_dir,'NDays',1,'Visible','off'), ...
    'Plot_Subfield_12:CropTimeLengthMismatch'), ...
    'plot rejects crop arrays with mismatched daily lengths');

fprintf('  ---------------------------------------\n');
fprintf('  %d passed, 0 failed\n',n_pass);

    function check(ok,label)
        assert(ok,'test_plot_subfield_12:Failed','%s',label);
        n_pass = n_pass + 1;
        fprintf('  PASS  %s\n',label);
    end

end

function ok = time_lines_match(fig,line_tag,expected)
ln = findall(fig,'Type','line','Tag',line_tag);
if numel(ln) ~= size(expected,1)
    ok = false;
    return
end
actual = zeros(numel(ln),size(expected,2));
for k = 1:numel(ln)
    actual(k,:) = ln(k).YData(:).';
end
ok = max(abs(sortrows(actual)-sortrows(expected)),[],'all') < 1e-12;
end

function ok = time_lines_have_daily_dates(fig,expected_count)
ln = findall(fig,'Type','line','Tag','CumulativeN2O');
ok = numel(ln) == expected_count;
for k = 1:numel(ln)
    x = ln(k).XData;
    ok = ok && isdatetime(x) && x(1) == datetime(2016,11,8) && ...
        isequal(days(diff(x(:))).',ones(1,3));
end
end

function [NPP_H,NPP_L,Ccrown_t,R_litter,R_microbe,R_ew,RB_L,topo, ...
        ColCfg,nee_daily,yield_daily] = carbon_fixture(k,R)
n_time = 4;
cc = 2;
if R.use_subdomains
    n_col = 2;
    topo.cell_area_m2 = [3; 1];
else
    n_col = 1;
    topo.cell_area_m2 = 4;
end
NPP_H = zeros(n_time,cc,n_col);
NPP_L = zeros(n_time,cc,n_col);
Ccrown_t = zeros(n_time,cc,n_col);
R_litter = zeros(n_time,1,n_col);
R_microbe = zeros(n_time,1,n_col);
R_ew = zeros(n_time,1,n_col);
RB_L = zeros(n_time,cc,7,n_col);
ColCfg = repmat(struct('PFT_OVERRIDE',''),n_col,1);

if n_col == 1
    NPP_L(:,1,1) = k*[0; 2; 4; 6];
    NPP_L(:,2,1) = 1e6;
    NPP_H(:,2,1) = 1e6;
    Ccrown_t(:,1,1) = 1;
    R_litter(:,1,1) = k;
    R_microbe(:,1,1) = 2*k;
    R_ew(:,1,1) = 3*k;
    RB_L(:,1,5,1) = k*[0; 0; 8; 0];
    nee_daily = k*[6 4 2 0];
    yield_daily = k*[0 0 8 0];
    return
end

NPP_L(:,1,1) = k*[0; 2; 4; 6];
NPP_L(:,2,1) = 1e6;
NPP_H(:,2,1) = 1e6;
Ccrown_t(:,1,1) = 1;
R_litter(:,1,1) = k;
R_microbe(:,1,1) = 2*k;
R_ew(:,1,1) = 3*k;
RB_L(:,1,5,1) = k*[0; 0; 8; 0];

NPP_L(:,1,2) = k*[0; 4; 8; 12];
Ccrown_t(:,1,2) = 1;
R_litter(:,1,2) = 2*k;
R_microbe(:,1,2) = 4*k;
R_ew(:,1,2) = 6*k;
RB_L(:,1,5,2) = k*[0; 0; 4; 0];
if strcmp(R.layout,'margin_prairie')
    NPP_H(:,1,2) = k*[0; 6; 12; 18];
    ColCfg(2).PFT_OVERRIDE = 'PFT_Prairie_Margin';
    nee_daily = k*[7.5 3.5 -0.5 -4.5];
    yield_daily = k*[0 0 6 0];
else
    nee_daily = k*[7.5 5 2.5 0];
    yield_daily = k*[0 0 7 0];
end
end

function ok = feedstock_markers_match(fig,expected_count)
marker = findall(fig,'Type','ConstantLine','Tag','FeedstockApplication');
ok = numel(marker) == expected_count;
for k = 1:numel(marker)
    ok = ok && marker(k).Value == datetime(2016,11,10) && ...
        strcmp(marker(k).LineStyle,':') && ...
        max(abs(marker(k).Color-[0.45 0.45 0.45])) < 1e-12;
end
end

function ok = clean_axes_and_bottom_legend(fig,check_time_ticks,n_group)
ax = findall(fig,'Type','axes');
ok = ~isempty(ax) && all(strcmp({ax.XGrid},'off')) && ...
    all(strcmp({ax.YGrid},'off'));
lgd = findall(fig,'Type','legend');
ok = ok && isscalar(lgd) && strcmp(lgd.Layout.Tile,'south') && ...
    isequal(lgd.TextColor,[0 0 0]);
if ~check_time_ticks
    upper = findall(fig,'Type','axes','-regexp','Tag','^SummaryRow[1-2]Column');
    bottom = findall(fig,'Type','axes','-regexp','Tag','^SummaryRow3Column');
    ok = ok && numel(upper) == 4 && numel(bottom) == 2 && ...
        all(arrayfun(@(a) isempty(a.XTickLabel),upper)) && ...
        all(arrayfun(@(a) ~isempty(a.XTickLabel),bottom));
    return
end
is_bottom = arrayfun(@(a) strcmp(a.XLabel.String,'Date'),ax);
upper = ax(~is_bottom);
bottom = ax(is_bottom);
contract = [numel(upper) == 4*n_group, numel(bottom) == 2*n_group, ...
    all(strcmp({upper.XTickLabelMode},'manual')), ...
    all(arrayfun(@(a) isempty(a.XTickLabel),upper)), ...
    all(strcmp({bottom.XTickLabelMode},'auto')), ...
    any(strcmp(string(lgd.String),'Feedstock application'))];
ok = ok && all(contract);
end

function ok = object_colors_match(fig,object_type,expected,n_repeat)
obj = findall(fig,'Type',object_type);
ok = numel(obj) == size(expected,1)*n_repeat;
actual = zeros(numel(obj),3);
for k = 1:numel(obj)
    if strcmp(object_type,'bar')
        color = obj(k).FaceColor;
    else
        color = obj(k).Color;
    end
    if ~isnumeric(color) || numel(color) ~= 3
        ok = false;
        return
    end
    actual(k,:) = color;
end
ok = ok && max(abs(sortrows(actual)-sortrows(repmat(expected,n_repeat,1))), ...
    [],'all') < 1e-12;
end

function ok = panel_has_values(fig,panel_title,expected)
ax = findall(fig,'Type','axes');
match = false(size(ax));
for k = 1:numel(ax)
    match(k) = strcmp(ax(k).Title.String,panel_title);
end
if sum(match) ~= 1
    ok = false;
    return
end
b = findall(ax(match),'Type','bar');
value = zeros(0,1);
for k = 1:numel(b)
    value = [value; b(k).YData(:)]; %#ok<AGROW>
end
ok = numel(value) == numel(expected) && ...
    max(abs(sort(value)-sort(expected(:)))) < 1e-12;
end

function ok = raises_error(f,expected_id)
ok = false;
try
    f();
catch ME
    ok = strcmp(ME.identifier,expected_id);
end
fig = findall(0,'Type','figure', ...
    '-regexp','Name','T&C-SMEW-N twelve-scenario');
if ~isempty(fig); close(fig); end
end

function ok = companion_fig_exists(png_file)
[folder,name] = fileparts(png_file);
fig_file = fullfile(folder,[name '.fig']);
if ~isfile(fig_file)
    ok = false;
    return
end
saved_fig = openfig(fig_file,'invisible');
ok = isgraphics(saved_fig,'figure');
close(saved_fig);
end

function tidy(varargin)
for k = 1:nargin
    if isfolder(varargin{k}); rmdir(varargin{k},'s'); end
end
end
