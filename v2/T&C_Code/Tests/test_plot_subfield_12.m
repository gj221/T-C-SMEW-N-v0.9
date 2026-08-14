function test_plot_subfield_12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Twelve-scenario summary plot contract.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n=== test_plot_subfield_12 ===\n');
n_pass = 0;
out_dir = tempname;
mkdir(out_dir);
cleanup = onCleanup(@() tidy(out_dir));

RM = Run_Matrix_12();
for k = 1:numel(RM)
    R = RM(k);
    Domain.N2Oflx = [0; k];
    Domain.LEAK_NO3 = [0; 2*k];
    Domain.NO3_pool = [10; 10+k];
    Domain.PH_d = [6; 6+0.01*k];
    save(fullfile(out_dir,sprintf('TC_2D_run%02d_%s__1d.mat', ...
        k,R.name)),'R','Domain');
end

[summary,fig,plot_file] = Plot_Subfield_12( ...
    'ResultDirectory',out_dir,'NDays',1,'Visible','off');
check(height(summary) == 12, 'plot summary contains all twelve scenarios');
check(isequal(summary.RunID.',1:12), 'plot summary preserves matrix order');
check(isequal(summary.CumulativeN2O_gN_m2.',1:12), ...
    'cumulative N2O values are read from each result');
check(isequal(summary.CumulativeNO3Leaching_gN_m2.',2*(1:12)), ...
    'cumulative nitrate leaching values are read from each result');
check(isgraphics(fig,'figure') && isfile(plot_file), ...
    'plot is created and written to disk');
check(panel_has_values(fig,'Cumulative N_2O',1:12), ...
    'N2O panel plots all twelve scenarios');
check(panel_has_values(fig,'Cumulative NO_3 leaching',2*(1:12)), ...
    'nitrate-leaching panel plots all twelve scenarios');
check(panel_has_values(fig,'Final soil NO_3 pool',10+(1:12)), ...
    'final-NO3 panel plots all twelve scenarios');
check(panel_has_values(fig,'Final soil pH',6+0.01*(1:12)), ...
    'pH panel plots all twelve scenarios');
close(fig);

bad_file = fullfile(out_dir,sprintf('TC_2D_run%02d_%s__1d.mat',12,RM(12).name));
S = load(bad_file);
S.Domain.PH_d(end) = NaN;
R = S.R; Domain = S.Domain;
save(bad_file,'R','Domain');
check(raises_error(@() Plot_Subfield_12( ...
    'ResultDirectory',out_dir,'NDays',1,'Visible','off'), ...
    'Plot_Subfield_12:NonfiniteMeasurement'), ...
    'plot rejects a scenario with a nonfinite final measurement');
Domain.PH_d(end) = 6.12;
Domain.N2Oflx(:) = NaN;
save(bad_file,'R','Domain');
check(raises_error(@() Plot_Subfield_12( ...
    'ResultDirectory',out_dir,'NDays',1,'Visible','off'), ...
    'Plot_Subfield_12:NonfiniteMeasurement'), ...
    'plot rejects a scenario with no finite cumulative measurement');

fprintf('  ---------------------------------------\n');
fprintf('  %d passed, 0 failed\n',n_pass);

    function check(ok,label)
        assert(ok,'test_plot_subfield_12:Failed','%s',label);
        n_pass = n_pass + 1;
        fprintf('  PASS  %s\n',label);
    end

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
    'Name','T&C-SMEW-N twelve-scenario summary');
if ~isempty(fig); close(fig); end
end

function tidy(out_dir)
if isfolder(out_dir); rmdir(out_dir,'s'); end
end
