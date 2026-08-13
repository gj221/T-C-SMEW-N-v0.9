%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Diagnostic plots: N2O emissions and NO3 leaching - CONTROL run only
%%%% Kantola (Urbana, US) 2016-2020, ERW_run case 5 (crop_control)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all;

col_N = [0 0.45 0.74]; % line/marker colour for the N plots

res_dir = fullfile(cd,'SIM','Kantola_pH_N2O');
res_file = fullfile(res_dir,'TC_ERW_results_Kantola_crop_control.mat');
if ~exist(res_file,'file')
    error('Missing result file: %s',res_file)
end
S = load(res_file,'N2Oflx','N2flx','VOL','LEAK_NO3','LEAK_NH4','P','PH_d','Datam','NNd');

% daily date axis (day 1 = simulation start; days are contiguous)
Date_d = datetime(S.Datam(1,1),S.Datam(1,2),S.Datam(1,3)) + caldays(0:S.NNd-1)';

figure('Position',[50 50 1250 1150]);
tiledlayout(4,2,'TileSpacing','compact','Padding','compact');

%%%% (1) daily N2O emission
nexttile; hold on; grid on;
plot(Date_d,S.N2Oflx,'Color',col_N,'LineWidth',1.1)
ylabel('N_2O flux [gN m^{-2} d^{-1}]')
title('Daily N_2O emission')

%%%% (2) cumulative N2O emission
nexttile; hold on; grid on;
plot(Date_d,cumsum(S.N2Oflx),'Color',col_N,'LineWidth',1.2)
ylabel('cumulative N_2O [gN m^{-2}]')
title('Cumulative N_2O emission')

%%%% (3) daily NO3 leaching
nexttile; hold on; grid on;
plot(Date_d,S.LEAK_NO3,'Color',col_N,'LineWidth',1.1)
ylabel('NO_3^- leaching [gN m^{-2} d^{-1}]')
title('Daily NO_3^- leaching')

%%%% (4) cumulative NO3 leaching
nexttile; hold on; grid on;
plot(Date_d,cumsum(S.LEAK_NO3),'Color',col_N,'LineWidth',1.2)
ylabel('cumulative NO_3^- leached [gN m^{-2}]')
title('Cumulative NO_3^- leaching')

%%%% (5) driving pore-water pH
nexttile; hold on; grid on;
plot(Date_d,S.PH_d,'Color',col_N,'LineWidth',1.1)
ylabel('pH [-]')
title('Pore-water pH driving the N cycle')

%%%% (6) soil NO3 pool
nexttile; hold on; grid on;
plot(Date_d,S.P(:,32),'Color',col_N,'LineWidth',1.1)
ylabel('NO_3^- pool [gN m^{-2}]')
title('Soil NO_3^- pool')

%%%% (7) daily N2O flux vs driving pH
nexttile; hold on; grid on;
scatter(S.PH_d,S.N2Oflx,6,col_N,'filled','MarkerFaceAlpha',0.35)
xlabel('pore-water pH [-]'); ylabel('N_2O flux [gN m^{-2} d^{-1}]')
title('N_2O emission vs pH')

%%%% (8) realized N2O/(N2O+N2) ratio vs pH, against the model partition curve
nexttile; hold on; grid on;
den_gas = S.N2Oflx+S.N2flx;
ok = den_gas > 0;
scatter(S.PH_d(ok),S.N2Oflx(ok)./den_gas(ok),6,col_N,'filled','MarkerFaceAlpha',0.35)
pH_th = linspace(4.5,8.5,200);
plot(pH_th,1./(1+exp(1.5*(pH_th-6.2))),'r--','LineWidth',1.2)
xlabel('pore-water pH [-]'); ylabel('N_2O/(N_2O+N_2) [-]')
title('Denitrification product ratio vs pH')
legend({'daily values','partition curve'},'Location','best','AutoUpdate','off')

exportgraphics(gcf,fullfile(res_dir,'Nflux_diagnostics_Kantola_control.png'),'Resolution',300)

%%%% run totals
fprintf('\nCONTROL totals 2016-2020: N2O=%.3f  N2=%.3f  NH3 vol.=%.3f  NO3 leach=%.3f  NH4 leach=%.3f [gN/m2]; mean pH=%.3f\n',...
    sum(S.N2Oflx),sum(S.N2flx),sum(S.VOL),sum(S.LEAK_NO3),sum(S.LEAK_NH4),mean(S.PH_d(S.PH_d>0)))
