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
close(fig);

fprintf('  ---------------------------------------\n');
fprintf('  %d passed, 0 failed\n',n_pass);

    function check(ok,label)
        assert(ok,'test_plot_subfield_12:Failed','%s',label);
        n_pass = n_pass + 1;
        fprintf('  PASS  %s\n',label);
    end

end

function tidy(out_dir)
if isfolder(out_dir); rmdir(out_dir,'s'); end
end
