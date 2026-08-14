%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Diagnostic plots: N2O emissions and NO3 leaching
%%%% Kantola (Urbana, US) 2016-2020: control vs basalt (Beerling) vs calcite
%%%% Run after the ERW_run cases (5, 3, 7) have been saved to
%%%% SIM\Kantola_pH_N2O ; scenarios not yet run are skipped.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all;

res_dir = fullfile(cd,'SIM','Kantola_pH_N2O');
all_files = {'TC_ERW_results_Kantola_crop_control.mat';...
    'TC_ERW_results_Kantola_crop_Beerling.mat';...
    'TC_ERW_results_Kantola_crop_calcite.mat'};
all_names = {'control','basalt','calcite'};
all_col = [0 0 0; 0 0.20 0.55; 0.35 0.70 0.90]; % control, basalt, calcite
col_appl = [0.45 0.45 0.45]; % rock/calcite application date markers
lw_appl = 1.8;               % their line width

vars = {'N2Oflx','N2flx','VOL','LEAK_NO3','LEAK_NH4','P','PH_d','Datam','NNd'};
scen_names = {}; scen_col = []; n_sc = 0;
for k = 1:3
    res_file = fullfile(res_dir,all_files{k});
    if exist(res_file,'file')
        n_sc = n_sc+1;
        S(n_sc) = load(res_file,vars{:});
        scen_names{n_sc} = all_names{k};
        scen_col(n_sc,:) = all_col(k,:);
    else
        fprintf('not yet available, skipping: %s\n',all_files{k});
    end
end
if n_sc == 0; error('No result files in %s',res_dir); end
i_ctr = find(strcmp(scen_names,'control'),1);

% daily date axis (day 1 = simulation start; days are contiguous)
Date_d = datetime(S(1).Datam(1,1),S(1).Datam(1,2),S(1).Datam(1,3)) + caldays(0:S(1).NNd-1)';
% rock/calcite application dates (10 Nov each year)
appl_dates = datetime([2016 11 10; 2017 11 10; 2018 11 10; 2019 11 10]);

figure('Position',[50 50 1250 950],'Color','w');
tiledlayout(3,2,'TileSpacing','compact','Padding','compact');
set(gcf,'DefaultAxesColor','w','DefaultAxesXColor','k','DefaultAxesYColor','k',...
    'DefaultTextColor','k','DefaultLegendColor','w','DefaultLegendTextColor','k',...
    'DefaultLegendEdgeColor','k')

%%%% (1) daily N2O emission
nexttile; hold on; grid on;
for sc = 1:n_sc
    plot(Date_d,S(sc).N2Oflx,'Color',scen_col(sc,:),'LineWidth',1.1)
end
xline(appl_dates,':','Color',col_appl,'LineWidth',lw_appl);
ylabel('N_2O flux [gN m^{-2} d^{-1}]')
title('Daily N_2O emission')
legend(scen_names,'Location','best','AutoUpdate','off')

%%%% (2) cumulative N2O emission
nexttile; hold on; grid on;
for sc = 1:n_sc
    plot(Date_d,cumsum(S(sc).N2Oflx),'Color',scen_col(sc,:),'LineWidth',1.2)
end
xline(appl_dates,':','Color',col_appl,'LineWidth',lw_appl);
ylabel('cumulative N_2O [gN m^{-2}]')
title('Cumulative N_2O emission')

%%%% (3) daily NO3 leaching
nexttile; hold on; grid on;
for sc = 1:n_sc
    plot(Date_d,S(sc).LEAK_NO3,'Color',scen_col(sc,:),'LineWidth',1.1)
end
xline(appl_dates,':','Color',col_appl,'LineWidth',lw_appl);
ylabel('NO_3^- leaching [gN m^{-2} d^{-1}]')
title('Daily NO_3^- leaching')

%%%% (4) cumulative NO3 leaching
nexttile; hold on; grid on;
for sc = 1:n_sc
    plot(Date_d,cumsum(S(sc).LEAK_NO3),'Color',scen_col(sc,:),'LineWidth',1.2)
end
xline(appl_dates,':','Color',col_appl,'LineWidth',lw_appl);
ylabel('cumulative NO_3^- leached [gN m^{-2}]')
title('Cumulative NO_3^- leaching')

%%%% (5) driving pore-water pH
nexttile; hold on; grid on;
for sc = 1:n_sc
    plot(Date_d,S(sc).PH_d,'Color',scen_col(sc,:),'LineWidth',1.1)
end
xline(appl_dates,':','Color',col_appl,'LineWidth',lw_appl);
ylabel('pH [-]')
title('Pore-water pH driving the N cycle')

%%%% (6) soil NO3 pool
nexttile; hold on; grid on;
for sc = 1:n_sc
    plot(Date_d,S(sc).P(:,32),'Color',scen_col(sc,:),'LineWidth',1.1)
end
xline(appl_dates,':','Color',col_appl,'LineWidth',lw_appl);
ylabel('NO_3^- pool [gN m^{-2}]')
title('Soil NO_3^- pool')

exportgraphics(gcf,fullfile(res_dir,'Nflux_diagnostics_Kantola.png'),'Resolution',300)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Figure 2: feedstock and pH influence on the N fluxes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Position',[80 80 1400 950],'Color','w');
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
set(gcf,'DefaultAxesColor','w','DefaultAxesXColor','k','DefaultAxesYColor','k',...
    'DefaultTextColor','k','DefaultLegendColor','w','DefaultLegendTextColor','k',...
    'DefaultLegendEdgeColor','k')

%%%% (1) total N loss budget by pathway and feedstock
nexttile; hold on; grid on;
loss_tot = zeros(5,n_sc);
for sc = 1:n_sc
    loss_tot(:,sc) = [sum(S(sc).N2Oflx); sum(S(sc).N2flx); sum(S(sc).VOL);...
        sum(S(sc).LEAK_NO3); sum(S(sc).LEAK_NH4)];
end
hb = bar(loss_tot);
for sc = 1:n_sc; hb(sc).FaceColor = scen_col(sc,:); end
set(gca,'XTick',1:5,'XTickLabel',{'N_2O','N_2','NH_3 vol.','NO_3^- leach','NH_4^+ leach'})
ylabel('total N loss 2016-2020 [gN m^{-2}]')
title('N loss budget by pathway')
legend(scen_names,'Location','best','AutoUpdate','off')

%%%% (2) pH change relative to control
nexttile; hold on; grid on;
if ~isempty(i_ctr) && n_sc > 1
    for sc = setdiff(1:n_sc,i_ctr)
        plot(Date_d,S(sc).PH_d-S(i_ctr).PH_d,'Color',scen_col(sc,:),'LineWidth',1.2)
    end
    xline(appl_dates,':','Color',col_appl,'LineWidth',lw_appl); yline(0,'k-');
    legend(scen_names(setdiff(1:n_sc,i_ctr)),'Location','best','AutoUpdate','off')
end
ylabel('\DeltapH vs control [-]')
title('Feedstock-driven pH change')

%%%% (3) daily N2O flux vs driving pH
nexttile; hold on; grid on;
for sc = 1:n_sc
    scatter(S(sc).PH_d,S(sc).N2Oflx,6,scen_col(sc,:),'filled','MarkerFaceAlpha',0.35)
end
xlabel('pore-water pH [-]'); ylabel('N_2O flux [gN m^{-2} d^{-1}]')
title('N_2O emission vs pH')

%%%% (4) cumulative N2O difference vs control
nexttile; hold on; grid on;
if ~isempty(i_ctr) && n_sc > 1
    for sc = setdiff(1:n_sc,i_ctr)
        plot(Date_d,cumsum(S(sc).N2Oflx)-cumsum(S(i_ctr).N2Oflx),'Color',scen_col(sc,:),'LineWidth',1.2)
    end
    xline(appl_dates,':','Color',col_appl,'LineWidth',lw_appl); yline(0,'k-');
end
ylabel('\Delta cumulative N_2O [gN m^{-2}]')
title('N_2O response to feedstock')

%%%% (5) cumulative NO3 leaching difference vs control
nexttile; hold on; grid on;
if ~isempty(i_ctr) && n_sc > 1
    for sc = setdiff(1:n_sc,i_ctr)
        plot(Date_d,cumsum(S(sc).LEAK_NO3)-cumsum(S(i_ctr).LEAK_NO3),'Color',scen_col(sc,:),'LineWidth',1.2)
    end
    xline(appl_dates,':','Color',col_appl,'LineWidth',lw_appl); yline(0,'k-');
end
ylabel('\Delta cumulative NO_3^- leached [gN m^{-2}]')
title('NO_3^- leaching response to feedstock')

%%%% (6) realized N2O/(N2O+N2) ratio vs pH, against the model partition curve
nexttile; hold on; grid on;
for sc = 1:n_sc
    den_gas = S(sc).N2Oflx+S(sc).N2flx;
    ok = den_gas > 0;
    scatter(S(sc).PH_d(ok),S(sc).N2Oflx(ok)./den_gas(ok),6,scen_col(sc,:),'filled','MarkerFaceAlpha',0.35)
end
pH_th = linspace(4.5,8.5,200);
plot(pH_th,1./(1+exp(1.5*(pH_th-6.2))),':','Color',col_appl,'LineWidth',lw_appl)
xlabel('pore-water pH [-]'); ylabel('N_2O/(N_2O+N_2) [-]')
title('Denitrification product ratio vs pH')
legend([scen_names {'partition curve'}],'Location','best','AutoUpdate','off')

exportgraphics(gcf,fullfile(res_dir,'Nflux_feedstock_pH_Kantola.png'),'Resolution',300)

%%%% run totals
fprintf('\n%-10s %13s %13s %13s %13s %13s %9s\n','scenario','N2O [gN/m2]','N2 [gN/m2]',...
    'NH3 vol.','NO3 leach','NH4 leach','mean pH')
for sc = 1:n_sc
    fprintf('%-10s %13.4f %13.4f %13.4f %13.4f %13.4f %9.3f\n',scen_names{sc},...
        sum(S(sc).N2Oflx),sum(S(sc).N2flx),sum(S(sc).VOL),sum(S(sc).LEAK_NO3),...
        sum(S(sc).LEAK_NH4),mean(S(sc).PH_d(S(sc).PH_d>0)))
end
