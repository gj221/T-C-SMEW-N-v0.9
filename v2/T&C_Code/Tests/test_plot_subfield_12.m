function test_plot_subfield_12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Twelve-scenario summary and time-dynamics plot contract.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n=== test_plot_subfield_12 ===\n');
n_pass = 0;
out_dir = tempname;
mkdir(out_dir);
cleanup = onCleanup(@() tidy(out_dir));

RM = Run_Matrix_12();
expected_vegetation_npp = zeros(numel(RM),4);
expected_crop_yield = zeros(numel(RM),4);
for k = 1:numel(RM)
    R = RM(k);
    Domain.N2Oflx = k*[0; 1; 2; 3];
    Domain.LEAK_NO3 = k*[0; 2; 4; 6];
    Domain.NO3_pool = 10+k*[0; 1; 2; 3];
    Domain.PH_d = 6+0.01*k*[0; 1; 2; 3];
    Datam = [2016 11 8 0];
    [NPP_H,NPP_L,Ccrown_t,RB_L,topo,ColCfg, ...
        expected_vegetation_npp(k,:),expected_crop_yield(k,:)] = ...
        vegetation_fixture(k,R);
    save(fullfile(out_dir,sprintf('TC_2D_run%02d_%s__1d.mat', ...
        k,R.name)),'R','Domain','Datam','NPP_H','NPP_L', ...
        'Ccrown_t','RB_L','topo','ColCfg');
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
check(max(abs(summary.CumulativeVegetationNPP_gC_m2- ...
    sum(expected_vegetation_npp,2))) < 1e-12, ...
    'summary area-weights cumulative vegetation NPP');
check(max(abs(summary.CumulativeCropYield_gC_m2- ...
    sum(expected_crop_yield,2))) < 1e-12, ...
    'summary area-weights crop yield');
check(all(summary.CumulativeVegetationNPP_gC_m2(7:9) > ...
    summary.CumulativeVegetationNPP_gC_m2(4:6)), ...
    'prairie PFTs contribute to whole-domain NPP');
check(all(summary.CumulativeCropYield_gC_m2(7:9) < ...
    summary.CumulativeCropYield_gC_m2(4:6)), ...
    'prairie PFTs contribute no crop yield');
check(isgraphics(fig,'figure') && isfile(plot_file), ...
    'plot is created and written to disk');
check(panel_has_values(fig,'Cumulative N_2O',6*(1:12)), ...
    'N2O panel plots all twelve scenarios');
check(panel_has_values(fig,'Cumulative NO_3 leaching',12*(1:12)), ...
    'nitrate-leaching panel plots all twelve scenarios');
check(panel_has_values(fig,'Final soil NO_3 pool',10+3*(1:12)), ...
    'final-NO3 panel plots all twelve scenarios');
check(panel_has_values(fig,'Final pore-water pH',6+0.03*(1:12)), ...
    'pore-water-pH panel plots all twelve scenarios');
check(panel_has_values(fig,'Cumulative vegetation NPP', ...
    sum(expected_vegetation_npp,2)), ...
    'vegetation-NPP panel plots all twelve scenarios');
check(panel_has_values(fig,'Cumulative crop yield', ...
    sum(expected_crop_yield,2)), ...
    'crop-yield panel plots all twelve scenarios');
expected_colors = [0 0 0; 0 0.25 0.60; 0.35 0.70 0.90];
check(object_colors_match(fig,'bar',expected_colors,6), ...
    'summary bars use black, dark blue and light blue treatments');
check(clean_axes_and_bottom_legend(fig,false), ...
    'summary plot has no grid and uses a bottom legend with black text');
check(isgraphics(time_fig,'figure') && isfile(time_plot_file), ...
    'time-dynamics plot is created and written to disk');
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
check(time_lines_match(time_fig,'VegetationNPP',expected_vegetation_npp), ...
    'daily vegetation-NPP dynamics contain all twelve scenarios');
check(time_lines_match(time_fig,'CumulativeCropYield', ...
    cumsum(expected_crop_yield,2)), ...
    'cumulative crop-yield dynamics contain all twelve scenarios');
check(time_lines_have_daily_dates(time_fig), ...
    'time-dynamics traces use the saved daily date axis');
check(object_colors_match(time_fig,'line',expected_colors,24), ...
    'time traces use black, dark blue and light blue treatments');
check(feedstock_markers_match(time_fig), ...
    'each time panel marks the in-window feedstock application date');
check(clean_axes_and_bottom_legend(time_fig,true), ...
    'time plot removes grids, repeated upper dates and uses a bottom legend');
close([fig time_fig]);

bad_file = fullfile(out_dir,sprintf('TC_2D_run%02d_%s__1d.mat',12,RM(12).name));
S = load(bad_file);
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
if numel(ln) ~= 12
    ok = false;
    return
end
actual = zeros(numel(ln),size(expected,2));
for k = 1:numel(ln)
    actual(k,:) = ln(k).YData(:).';
end
ok = max(abs(sortrows(actual)-sortrows(expected)),[],'all') < 1e-12;
end

function ok = time_lines_have_daily_dates(fig)
ln = findall(fig,'Type','line','Tag','CumulativeN2O');
ok = numel(ln) == 12;
for k = 1:numel(ln)
    x = ln(k).XData;
    ok = ok && isdatetime(x) && x(1) == datetime(2016,11,8) && ...
        isequal(days(diff(x(:))).',ones(1,3));
end
end

function [NPP_H,NPP_L,Ccrown_t,RB_L,topo,ColCfg,npp_daily,yield_daily] = ...
        vegetation_fixture(k,R)
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
RB_L = zeros(n_time,cc,7,n_col);
ColCfg = repmat(struct('PFT_OVERRIDE',''),n_col,1);

if n_col == 1
    NPP_L(:,1,1) = k*[0; 2; 4; 6];
    NPP_L(:,2,1) = 1e6;
    NPP_H(:,2,1) = 1e6;
    Ccrown_t(:,1,1) = 1;
    RB_L(:,1,5,1) = k*[0; 0; 8; 0];
    npp_daily = k*[0 2 4 6];
    yield_daily = k*[0 0 8 0];
    return
end

NPP_L(:,1,1) = k*[0; 2; 4; 6];
NPP_L(:,2,1) = 1e6;
NPP_H(:,2,1) = 1e6;
Ccrown_t(:,1,1) = 1;
RB_L(:,1,5,1) = k*[0; 0; 8; 0];

NPP_L(:,1,2) = k*[0; 4; 8; 12];
Ccrown_t(:,1,2) = 1;
RB_L(:,1,5,2) = k*[0; 0; 4; 0];
if strcmp(R.layout,'margin_prairie')
    NPP_H(:,1,2) = k*[0; 6; 12; 18];
    ColCfg(2).PFT_OVERRIDE = 'PFT_Prairie_Margin';
    npp_daily = k*[0 4 8 12];
    yield_daily = k*[0 0 6 0];
else
    npp_daily = k*[0 2.5 5 7.5];
    yield_daily = k*[0 0 7 0];
end
end

function ok = feedstock_markers_match(fig)
marker = findall(fig,'Type','ConstantLine','Tag','FeedstockApplication');
ok = numel(marker) == 24;
for k = 1:numel(marker)
    ok = ok && marker(k).Value == datetime(2016,11,10) && ...
        strcmp(marker(k).LineStyle,':') && ...
        max(abs(marker(k).Color-[0.45 0.45 0.45])) < 1e-12;
end
end

function ok = clean_axes_and_bottom_legend(fig,check_time_ticks)
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
upper = findall(fig,'Type','axes','-regexp','Tag','^TimeRow[1-5]Group');
bottom = findall(fig,'Type','axes','-regexp','Tag','^TimeRow6Group');
ok = ok && numel(upper) == 20 && numel(bottom) == 4 && ...
    all(arrayfun(@(a) isempty(a.XTickLabel),upper)) && ...
    all(arrayfun(@(a) ~isempty(a.XTickLabel),bottom)) && ...
    any(strcmp(string(lgd.String),'Feedstock application'));
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
ok = numel(value) == 12 && ...
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

function tidy(out_dir)
if isfolder(out_dir); rmdir(out_dir,'s'); end
end
