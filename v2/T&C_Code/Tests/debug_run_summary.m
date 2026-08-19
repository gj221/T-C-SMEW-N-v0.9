% Debug: summarize the 10-day EQUIV_1D_TMP run state (pore vs bulk pH)
clear; clc;
run_file = fullfile(pwd,'..','..','SIM','EQUIV_1D_TMP','TC_ERW_results_Kantola_crop_control.mat');
S = load(run_file);
names = fieldnames(S);
fprintf('workspace variables (%d):\n', numel(names));
fprintf('  %s\n', names{:});
nstep = 144;
if isfield(S,'ERW_biogeochem')
    bg = S.ERW_biogeochem;
    fprintf('\nERW_biogeochem fields: %s\n', strjoin(fieldnames(bg), ', '));
    pH = bg.pH(:);
    np = numel(pH);
    fprintf('pH(pore) length=%d\n', np);
    fprintf('  day1 start : %.4f\n', pH(1));
    fprintf('  day1 end   : %.4f\n', pH(min(nstep,np)));
    fprintf('  day2 end   : %.4f\n', pH(min(2*nstep,np)));
    fprintf('  day3 end   : %.4f\n', pH(min(3*nstep,np)));
    fprintf('  day5 end   : %.4f\n', pH(min(5*nstep,np)));
    fprintf('  day10 end  : %.4f\n', pH(min(10*nstep,np)));
    fprintf('  min pH     : %.4f at step %d\n', min(pH), find(pH==min(pH),1));
    fprintf('  max pH     : %.4f at step %d\n', max(pH), find(pH==max(pH),1));
    if isfield(bg,'pH_bulk')
        pb = bg.pH_bulk(:);
        fprintf('pH_bulk length=%d\n', numel(pb));
        fprintf('  day1 end : %.4f\n', pb(min(nstep,numel(pb))));
        fprintf('  day10 end: %.4f\n', pb(min(10*nstep,numel(pb))));
    end
    if isfield(bg,'An')
        An = bg.An(:);
        fprintf('An (umol_c/L): day1 start %.3g end %.3g ; day2 end %.3g ; day10 end %.3g ; min %.3g\n', ...
            An(1), An(min(nstep,np)), An(min(2*nstep,np)), An(min(10*nstep,np)), min(An));
    end
    if isfield(bg,'Ca')
        Ca = bg.Ca(:);
        fprintf('Ca: day1 start %.3g end %.3g ; day2 end %.3g ; day10 end %.3g\n', ...
            Ca(1), Ca(min(nstep,np)), Ca(min(2*nstep,np)), Ca(min(10*nstep,np)));
    end
    if isfield(bg,'f_CEC')
        fC = bg.f_CEC;
        fprintf('f_CEC row1 (seed): [%s]\n', sprintf('%.4f ', fC(1,:)));
        fprintf('f_CEC day1 end  : [%s]\n', sprintf('%.4f ', fC(min(nstep,np),:)));
        fprintf('f_CEC day2 end  : [%s]\n', sprintf('%.4f ', fC(min(2*nstep,np),:)));
        fprintf('f_CEC day10 end : [%s]\n', sprintf('%.4f ', fC(min(10*nstep,np),:)));
    end
end
if isfield(S,'PH')
    fprintf('\nPH (last-day bulk): [%s]\n', sprintf('%.4f ', S.PH));
end
if isfield(S,'PH_pore')
    fprintf('PH_pore (last-day pore): [%s]\n', sprintf('%.4f ', S.PH_pore));
end