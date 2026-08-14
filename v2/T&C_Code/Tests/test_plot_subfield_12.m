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
for k = 1:numel(RM)
    R = RM(k);
    Domain.N2Oflx = k*[0; 1; 2; 3];
    Domain.LEAK_NO3 = k*[0; 2; 4; 6];
    Domain.NO3_pool = 10+k*[0; 1; 2; 3];
    Domain.PH_d = 6+0.01*k*[0; 1; 2; 3];
    Datam = [2016 1 1 0];
    save(fullfile(out_dir,sprintf('TC_2D_run%02d_%s__1d.mat', ...
        k,R.name)),'R','Domain','Datam');
end

[summary,fig,plot_file,time_fig,time_plot_file] = Plot_Subfield_12( ...
    'ResultDirectory',out_dir,'NDays',1,'Visible','off');
check(height(summary) == 12, 'plot summary contains all twelve scenarios');
check(isequal(summary.RunID.',1:12), 'plot summary preserves matrix order');
check(isequal(summary.CumulativeN2O_gN_m2.',6*(1:12)), ...
    'cumulative N2O values are read from each result');
check(isequal(summary.CumulativeNO3Leaching_gN_m2.',12*(1:12)), ...
    'cumulative nitrate leaching values are read from each result');
check(isgraphics(fig,'figure') && isfile(plot_file), ...
    'plot is created and written to disk');
check(panel_has_values(fig,'Cumulative N_2O',6*(1:12)), ...
    'N2O panel plots all twelve scenarios');
check(panel_has_values(fig,'Cumulative NO_3 leaching',12*(1:12)), ...
    'nitrate-leaching panel plots all twelve scenarios');
check(panel_has_values(fig,'Final soil NO_3 pool',10+3*(1:12)), ...
    'final-NO3 panel plots all twelve scenarios');
check(panel_has_values(fig,'Final soil pH',6+0.03*(1:12)), ...
    'pH panel plots all twelve scenarios');
expected_colors = [0 0 0; 0 0.25 0.60; 0.35 0.70 0.90];
check(object_colors_match(fig,'bar',expected_colors,4), ...
    'summary bars use black, dark blue and light blue treatments');
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
check(time_lines_match(time_fig,'SoilPH', ...
    6+0.01*(1:12).'*[0 1 2 3]), ...
    'soil-pH dynamics contain all twelve scenarios');
check(time_lines_have_daily_dates(time_fig), ...
    'time-dynamics traces use the saved daily date axis');
check(object_colors_match(time_fig,'line',expected_colors,16), ...
    'time traces use black, dark blue and light blue treatments');
close([fig time_fig]);

bad_file = fullfile(out_dir,sprintf('TC_2D_run%02d_%s__1d.mat',12,RM(12).name));
S = load(bad_file);
S.Domain.PH_d(end) = NaN;
R = S.R; Domain = S.Domain; Datam = S.Datam;
save(bad_file,'R','Domain','Datam');
check(raises_error(@() Plot_Subfield_12( ...
    'ResultDirectory',out_dir,'NDays',1,'Visible','off'), ...
    'Plot_Subfield_12:NonfiniteMeasurement'), ...
    'plot rejects a scenario with a nonfinite final measurement');
Domain.PH_d(end) = 6.36;
Domain.N2Oflx(:) = NaN;
save(bad_file,'R','Domain','Datam');
check(raises_error(@() Plot_Subfield_12( ...
    'ResultDirectory',out_dir,'NDays',1,'Visible','off'), ...
    'Plot_Subfield_12:NonfiniteMeasurement'), ...
    'plot rejects a scenario with no finite cumulative measurement');

Domain.N2Oflx = 12*[0; 1; 2; 3];
Domain.PH_d = Domain.PH_d(1:3);
save(bad_file,'R','Domain','Datam');
check(raises_error(@() Plot_Subfield_12( ...
    'ResultDirectory',out_dir,'NDays',1,'Visible','off'), ...
    'Plot_Subfield_12:TimeLengthMismatch'), ...
    'plot rejects measurements with mismatched daily lengths');

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
    ok = ok && isdatetime(x) && x(1) == datetime(2016,1,1) && ...
        isequal(days(diff(x(:))).',ones(1,3));
end
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
