clear all
close all
clc

current_directory = cd;

% simulation folder dir
SCENARIO_TEST = '\FINAL_field';

%% choose site
id_location = 'Kantola'; 
NN_case_study_ERW_MOD = [3 4 5 6]; % No. of experiment setup in ERW_MOD sheet

%%%%%%%%%%%%%% read ERW module parameters fromexcel
ERW_PARAM_opts = detectImportOptions('ERW_MOD_PARAM_parameters.xlsx',"Sheet",'ERW_MOD');
ERW_PARAM_opts.DataRange ="A2:AH8";
ERW_PARAM_opts.VariableTypes ={'double','string','string','double','double',...
  'string','double','double','string','double',...
  'string','string','string','string','double',...
  'double','double','double','double','double',...
  'double','double','double','string','double',...
  'string','string','double','double','double',...
  'double','double','double','double'};
ERW_PARAM_opts = setvaropts(ERW_PARAM_opts,{'f_CEC_in','rock_filepath','mineral','rock_f_in'},'FillValue','NaN');

ERW_PARAM_input_t = readtable('ERW_MOD_PARAM_parameters.xlsx', ERW_PARAM_opts);

% find the corresponding case study
ERW_PARAM_input_B = ERW_PARAM_input_t(ERW_PARAM_input_t.Number == NN_case_study_ERW_MOD(1),:);
ERW_PARAM_input_L = ERW_PARAM_input_t(ERW_PARAM_input_t.Number == NN_case_study_ERW_MOD(2),:);
ERW_PARAM_input_untr = ERW_PARAM_input_t(ERW_PARAM_input_t.Number == NN_case_study_ERW_MOD(3),:);
ERW_PARAM_input_pre = ERW_PARAM_input_t(ERW_PARAM_input_t.Number == NN_case_study_ERW_MOD(4),:);

%% load data
% check which data you need
fl_path = dir(fullfile(current_directory,'SIM',SCENARIO_TEST));
aTable = struct2table(fl_path);
disp(aTable);

% Basalt_L
cell_row = contains(aTable.name,ERW_PARAM_input_L.ExperimentalDescription+".mat",'IgnoreCase',true);
file_name = fullfile(fl_path(cell_row).folder,fl_path(cell_row).name);
data_Lewis = load(file_name);

% Basalt_B
cell_row = contains(aTable.name,ERW_PARAM_input_B.ExperimentalDescription+".mat",'IgnoreCase',true);
file_name = fullfile(fl_path(cell_row).folder,fl_path(cell_row).name);
data_Beerling = load(file_name);

% Control
cell_row_untr = contains(aTable.name,ERW_PARAM_input_untr.ExperimentalDescription,'IgnoreCase',true);
file_name_untr = fullfile(fl_path(cell_row_untr).folder,fl_path(cell_row_untr).name);
data_untr = load(file_name_untr);

% pretreatment
cell_row_pre = contains(aTable.name,'ERW_off.mat','IgnoreCase',true);
file_name_pre = fullfile(fl_path(cell_row_pre).folder,fl_path(cell_row_pre).name);
data_pre = load(file_name_pre);


%% uncertainty results 

dir_SA = fullfile(current_directory,'SIM','Kantola_uncertainty');
addpath(dir_SA)
codedir = fullfile(dir_SA, 'frank-pk-DataViz-3.2.3.0','daboxplot');
addpath(codedir);
Data_SA = load('uc_post_results.mat');

Lewis_pH = [Data_SA.Lewis_uc.pH_uc]';
Lewis_CEC = [Data_SA.Lewis_uc.CEC_uc]';
CDR_uc_L =  [Data_SA.Lewis_uc.CDR_uc];
Lewis_cdr = [CDR_uc_L.CO2_pot_num]';
Lewis_Ca_diss = [CDR_uc_L.Ca_diss_TC]';
Lewis_Mg_diss = [CDR_uc_L.Mg_diss_TC]';


Beerling_pH = [Data_SA.Beerling_uc.pH_uc]';
Beerling_CEC = [Data_SA.Beerling_uc.CEC_uc]';
CDR_uc_B =  [Data_SA.Beerling_uc.CDR_uc];
Beerling_cdr = [CDR_uc_B.CO2_pot_num]';
Beerling_Ca_diss = [CDR_uc_B.Ca_diss_TC]';
Beerling_Mg_diss = [CDR_uc_B.Mg_diss_TC]';

condition_names = {'2016', '2017', '2018', '2019' ,'2020'};
group_names = {'Lewis', 'Beerling' };

Lewis_pH_std = std(Lewis_pH);
Beerling_pH_std = std(Beerling_pH);
Lewis_CEC_std = std(Lewis_CEC);
Beerling_CEC_std = std(Beerling_CEC);

data2{1} = Lewis_cdr(:,2:5);
data2{2} = Beerling_cdr(:,2:5);

data3{1} = Lewis_Ca_diss(:,2:5);
data3{2} = Beerling_Ca_diss(:,2:5);

data4{1} = Lewis_Mg_diss(:,2:5);
data4{2} = Beerling_Mg_diss(:,2:5);

%% load exp data
fl_path = dir(fullfile('Experimental observations','Web plot digitizer','*.mat'));

% check which data you need
aTable = struct2table(fl_path);
disp(aTable);

% Kantola et al. - Soybean
cell_row = contains(aTable.name,'Soybean','IgnoreCase',true);
file_name = fullfile(fl_path(cell_row).folder,fl_path(cell_row).name);
data_soy = load(file_name);

% Kantola et al. - Maize
cell_row = contains(aTable.name,'Maize','IgnoreCase',true);
file_name = fullfile(fl_path(cell_row).folder,fl_path(cell_row).name);
data_maize = load(file_name);

% Beerling et al.
cell_row = contains(aTable.name,'Beerling','IgnoreCase',true);
file_name = fullfile(fl_path(cell_row).folder,fl_path(cell_row).name);
data_exp = load(file_name);


%% pretreatment: 2009-2016

[data_pre_s_day, data_pre_m_day] = Cflux_post_processing(data_pre, [2009 2011 2012 2014 2015], [2010 2013 2016]);

data_pre_GPP_s_day = data_pre_s_day.data_GPP_s_day;
data_pre_NEE_s_day = data_pre_s_day.data_NEE_s_day;
data_pre_Re_s_day = data_pre_s_day.data_Re_s_day;

data_pre_Raut_s_day = data_pre_s_day.data_Raut_s_day;
data_pre_Rhet_s_day = data_pre_s_day.data_Rhet_s_day;

data_pre_GPP_m_day = data_pre_m_day.data_GPP_m_day;
data_pre_NEE_m_day = data_pre_m_day.data_NEE_m_day;
data_pre_Re_m_day = data_pre_m_day.data_Re_m_day;

data_pre_Raut_m_day = data_pre_m_day.data_Raut_m_day;
data_pre_Rhet_m_day = data_pre_m_day.data_Rhet_m_day;


%% Lewis basalt: 2017-2020

[data_Lewis_s_day, data_Lewis_m_day] = Cflux_post_processing(data_Lewis, [2017 2018 2020], [2019]);

data_Lewis_GPP_s_day = data_Lewis_s_day.data_GPP_s_day;
data_Lewis_NEE_s_day = data_Lewis_s_day.data_NEE_s_day;
data_Lewis_Re_s_day = data_Lewis_s_day.data_Re_s_day;

data_Lewis_Raut_s_day = data_Lewis_s_day.data_Raut_s_day;
data_Lewis_Rhet_s_day = data_Lewis_s_day.data_Rhet_s_day;

data_Lewis_GPP_m_day = data_Lewis_m_day.data_GPP_m_day;
data_Lewis_NEE_m_day = data_Lewis_m_day.data_NEE_m_day;
data_Lewis_Re_m_day = data_Lewis_m_day.data_Re_m_day;

data_Lewis_Raut_m_day = data_Lewis_m_day.data_Raut_m_day;
data_Lewis_Rhet_m_day = data_Lewis_m_day.data_Rhet_m_day;


%% Beerling basalt: 2017-2020

[data_Beerling_s_day, data_Beerling_m_day] = Cflux_post_processing(data_Beerling, [2017 2018 2020], [2019]);

data_Beerling_GPP_s_day = data_Beerling_s_day.data_GPP_s_day;
data_Beerling_NEE_s_day = data_Beerling_s_day.data_NEE_s_day;
data_Beerling_Re_s_day = data_Beerling_s_day.data_Re_s_day;

data_Beerling_Raut_s_day = data_Beerling_s_day.data_Raut_s_day;
data_Beerling_Rhet_s_day = data_Beerling_s_day.data_Rhet_s_day;

data_Beerling_GPP_m_day = data_Beerling_m_day.data_GPP_m_day;
data_Beerling_NEE_m_day = data_Beerling_m_day.data_NEE_m_day;
data_Beerling_Re_m_day = data_Beerling_m_day.data_Re_m_day;

data_Beerling_Raut_m_day = data_Beerling_m_day.data_Raut_m_day;
data_Beerling_Rhet_m_day = data_Beerling_m_day.data_Rhet_m_day;


%% concurrent: 2017-2020

[data_untr_s_day, data_untr_m_day] = Cflux_post_processing(data_untr, [2017 2018 2020], [2019]);

data_untr_GPP_s_day = data_untr_s_day.data_GPP_s_day;
data_untr_NEE_s_day = data_untr_s_day.data_NEE_s_day;
data_untr_Re_s_day = data_untr_s_day.data_Re_s_day;

data_untr_Raut_s_day = data_untr_s_day.data_Raut_s_day;
data_untr_Rhet_s_day = data_untr_s_day.data_Rhet_s_day;

data_untr_GPP_m_day = data_untr_m_day.data_GPP_m_day;
data_untr_NEE_m_day = data_untr_m_day.data_NEE_m_day;
data_untr_Re_m_day = data_untr_m_day.data_Re_m_day;

data_untr_Raut_m_day = data_untr_m_day.data_Raut_m_day;
data_untr_Rhet_m_day = data_untr_m_day.data_Rhet_m_day;


%% Figure : ecosystem C flux in pretreatment years

line_width = 3;

stat_func = {
    @(sim,obs) (rmse(sim,obs));...
    @(sim,obs) (corrcoef(sim,obs).^2);...
    @(sim,obs) (1-sum((obs-sim).^2)/sum((obs-mean(obs)).^2));...
    };


figure;
f1 = tiledlayout(3, 2);

color_obs = '#d95f02';
color_model = '#1b9e77';

axs(1) = nexttile;
% maize: plot daily values [gC/m2/day]
hold on
plot(data_pre_GPP_m_day,'Color',color_model,'LineWidth',line_width)

plot([1:366], data_maize.pre_GPP.y_daily,'--','Color',color_obs,'LineWidth',line_width,'MarkerSize',5,...
    'MarkerFaceColor','none','MarkerEdgeColor',color_obs)


title('Maize')
xlim([0 366])
ylim([0 30])
ylabel({'GPP' ;'[gC m^{2} d^{-1}]'})
set(gca,'box','on')
h=gca; h.XAxis.TickLength = [0 0];
% text(axs(1).XLim(1), axs(1).YLim(2)*0.95, '(a-1)')
text(0.01, 0.85, '(a-1)','FontWeight','bold','Units','normalized','VerticalAlignment','bottom')

obs=data_maize.pre_GPP.y_daily;
sim=data_pre_GPP_m_day;

R2 = stat_func{2}(sim,obs);
text(axs(1).XLim(2)*0.8, axs(1).YLim(2)*0.7, sprintf('R^2=%.2f\nNSE=%.2f',... 
    R2(1,2),...
    stat_func{3}(sim,obs)),'FontSize',8)


axs(2) = nexttile;
% soybean: plot daily values [gC/m2/day]

hold on
plot(data_pre_GPP_s_day,'Color',color_model,'LineWidth',line_width)

plot([1:366], data_soy.pre_GPP.y_daily,'--','Color',color_obs,'LineWidth',line_width,'MarkerSize',5,...
    'MarkerFaceColor','none','MarkerEdgeColor',color_obs)

title('Soybean')
xlim([0 366])
ylim([0 30])
set(gca,'box','on')
h=gca; h.XAxis.TickLength = [0 0];
% text(axs(2).XLim(1), axs(2).YLim(2)*0.95, '(a-2)')
text(0.01, 0.85, '(a-2)','FontWeight','bold','Units','normalized','VerticalAlignment','bottom')

obs=data_soy.pre_GPP.y_daily;
sim=data_pre_GPP_s_day;

R2 = stat_func{2}(sim,obs);
text(axs(2).XLim(2)*0.8, axs(2).YLim(2)*0.7, sprintf('R^2=%.2f\nNSE=%.2f',... 
    R2(1,2),...
    stat_func{3}(sim,obs)),'FontSize',8)

%%%% NEE daily
axs(3) = nexttile;
% maize: plot daily values [gC/m2/day]

hold on
plot(data_pre_NEE_m_day,'Color',color_model,'LineWidth',line_width)

plot([1:366], data_maize.pre_NEE.y_daily,'--','Color',color_obs,'LineWidth',line_width,'MarkerSize',5,...
    'MarkerFaceColor','none','MarkerEdgeColor',color_obs)
ylabel({'NEE'; '[gC m^{2} d^{-1}]'})
xlim([0 366])
ylim([-10 5])
set(gca,'box','on')
h=gca; h.XAxis.TickLength = [0 0];
% text(axs(3).XLim(1), axs(3).YLim(2)*0.9, '(b-1)')
text(0.01, 0.85, '(b-1)','FontWeight','bold','Units','normalized','VerticalAlignment','bottom')


obs=data_maize.pre_NEE.y_daily;
sim=data_pre_NEE_m_day;

R2 = stat_func{2}(sim,obs);
text(axs(3).XLim(2)*0.8, -4.5, sprintf('R^2=%.2f\nNSE=%.2f',... 
    R2(1,2),...
    stat_func{3}(sim,obs)),'FontSize',8)

axs(4) = nexttile;
% soybean: plot daily values [gC/m2/day]

hold on
plot(data_pre_NEE_s_day,'Color',color_model,'LineWidth',line_width)

plot([1:366], data_soy.pre_NEE.y_daily,'--','Color',color_obs,'LineWidth',line_width,'MarkerSize',5,...
    'MarkerFaceColor','none','MarkerEdgeColor',color_obs)
xlim([0 366])
ylim([-10 5])
set(gca,'box','on')
h=gca; h.XAxis.TickLength = [0 0];
% text(axs(4).XLim(1), axs(4).YLim(2)*0.9, '(b-2)')
text(0.01, 0.85, '(b-2)','FontWeight','bold','Units','normalized','VerticalAlignment','bottom')


obs=data_soy.pre_NEE.y_daily;
sim=data_pre_NEE_s_day;

R2 = stat_func{2}(sim,obs);
text(axs(4).XLim(2)*0.8, -4.5, sprintf('R^2=%.2f\nNSE=%.2f',... 
    R2(1,2),...
    stat_func{3}(sim,obs)),'FontSize',8)

%%% RE daily
axs(5) = nexttile;
% maize: plot daily values [gC/m2/day]

hold on
plot(data_pre_Re_m_day,'Color',color_model,'LineWidth',line_width)

plot([1:366], data_maize.pre_Re.y_daily,'--','Color',color_obs,'LineWidth',line_width,'MarkerSize',5,...
    'MarkerFaceColor','none','MarkerEdgeColor',color_obs)
ylabel({'Re' ;'[gC m^{2} d^{-1}]'})
xlim([0 366])
ylim([0 20])
set(gca,'box','on')
h=gca; h.XAxis.TickLength = [0 0];
% text(axs(5).XLim(1), axs(5).YLim(2)*0.95, '(c-1)')
text(0.01, 0.85, '(c-1)','FontWeight','bold','Units','normalized','VerticalAlignment','bottom')

obs=data_maize.pre_Re.y_daily;
sim=data_pre_Re_m_day;

R2 = stat_func{2}(sim,obs);
text(axs(5).XLim(2)*0.8, axs(5).YLim(2)*0.7, sprintf('R^2=%.2f\nNSE=%.2f',... 
    R2(1,2),...
    stat_func{3}(sim,obs)),'FontSize',8)

axs(6) = nexttile;
% soybean: plot daily values [gC/m2/day]

hold on
plot(data_pre_Re_s_day,'Color',color_model,'LineWidth',line_width)

plot([1:366], data_soy.pre_Re.y_daily,'--','Color',color_obs,'LineWidth',line_width,'MarkerSize',5,...
    'MarkerFaceColor','none','MarkerEdgeColor',color_obs)
xlim([0 366])
ylim([0 20])
set(gca,'box','on')
h=gca; h.XAxis.TickLength = [0 0];
% text(axs(6).XLim(1), axs(6).YLim(2)*0.95, '(c-2)')
text(0.01, 0.85, '(c-2)','FontWeight','bold','Units','normalized','VerticalAlignment','bottom')

obs=data_soy.pre_Re.y_daily;
sim=data_pre_Re_s_day;

R2 = stat_func{2}(sim,obs);
text(axs(6).XLim(2)*0.8, axs(6).YLim(2)*0.7, sprintf('R^2=%.2f\nNSE=%.2f',... 
    R2(1,2),...
    stat_func{3}(sim,obs)),'FontSize',8)


xlabel(f1, 'Day of Year')

set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.01 0.01 0.98 0.98]);
legend({'T&C-SMEW','Exp Flux'},'Location','north')
fontsize(gcf,scale=2.5)

%% SI figure - time series ecosystem C flux control vs basalt 

% maize
figure;
f8 = tiledlayout(3, 2);

axs(1) = nexttile;
% maize: plot daily values [gC/m2/day]
plot([1:366], data_maize.control_GPP.y_daily,'--','Color',color_obs,'LineWidth',line_width,'MarkerSize',5,...
    'MarkerFaceColor','none','MarkerEdgeColor',color_obs)
hold on
plot(data_untr_GPP_m_day,'Color',color_model,'LineWidth',line_width)
title('Maize - control')
xlim([0 366])
ylim([0 30])
ylabel('GPP [gC m^{2} d^{-1}]')
set(gca,'box','on')
h=gca; h.XAxis.TickLength = [0 0];
% text(axs(1).XLim(1), axs(1).YLim(2)*0.95, '(a-1)')
text(0.01, 0.85, '(a-1)','FontWeight','bold','Units','normalized','VerticalAlignment','bottom')

obs=data_maize.control_GPP.y_daily;
sim=data_untr_GPP_m_day;

R2 = stat_func{2}(sim,obs);
text(axs(1).XLim(2)*0.8, axs(1).YLim(2)*0.7, sprintf('R^2=%.2f\nNSE=%.2f',... 
    R2(1,2),...
    stat_func{3}(sim,obs)),'FontSize',8)


axs(2) = nexttile;
% soybean: plot daily values [gC/m2/day]
plot([1:366], data_maize.basalt_GPP.y_daily,'--','Color',color_obs,'LineWidth',line_width,'MarkerSize',5,...
    'MarkerFaceColor','none','MarkerEdgeColor',color_obs)
hold on
plot(data_Lewis_GPP_m_day,'Color',color_model,'LineWidth',line_width)
plot(data_Beerling_GPP_m_day,'Color',color_model,'LineWidth',line_width)

title('Maize - basalt')
xlim([0 366])
ylim([0 30])
set(gca,'box','on')
h=gca; h.XAxis.TickLength = [0 0];
% text(axs(2).XLim(1), axs(2).YLim(2)*0.95, '(a-2)')
text(0.01, 0.85, '(a-2)','FontWeight','bold','Units','normalized','VerticalAlignment','bottom')

obs=data_maize.basalt_GPP.y_daily;
sim=data_Lewis_GPP_m_day;

R2 = stat_func{2}(sim,obs);
text(axs(2).XLim(2)*0.8, axs(2).YLim(2)*0.7, sprintf('R^2=%.2f\nNSE=%.2f',... 
    R2(1,2),...
    stat_func{3}(sim,obs)),'FontSize',8)

%%%% NEE daily
axs(3) = nexttile;
% maize: plot daily values [gC/m2/day]
plot([1:366], data_maize.control_NEE.y_daily,'--','Color',color_obs,'LineWidth',line_width,'MarkerSize',5,...
    'MarkerFaceColor','none','MarkerEdgeColor',color_obs)
hold on
plot(data_untr_NEE_m_day,'Color',color_model,'LineWidth',line_width)
ylabel('NEE [gC m^{2} d^{-1}]')
xlim([0 366])
ylim([-15 5])
set(gca,'box','on')
h=gca; h.XAxis.TickLength = [0 0];
% text(axs(3).XLim(1), axs(3).YLim(2)*0.9, '(b-1)')
text(0.01, 0.85, '(b-1)','FontWeight','bold','Units','normalized','VerticalAlignment','bottom')


obs=data_maize.control_NEE.y_daily;
sim=data_untr_NEE_m_day;

R2 = stat_func{2}(sim,obs);
text(axs(3).XLim(2)*0.8, -4.5, sprintf('R^2=%.2f\nNSE=%.2f',... 
    R2(1,2),...
    stat_func{3}(sim,obs)),'FontSize',8)

axs(4) = nexttile;
% soybean: plot daily values [gC/m2/day]
plot([1:366], data_maize.basalt_NEE.y_daily,'--','Color',color_obs,'LineWidth',line_width,'MarkerSize',5,...
    'MarkerFaceColor','none','MarkerEdgeColor',color_obs)
hold on
plot(data_Lewis_NEE_m_day,'Color',color_model,'LineWidth',line_width)
plot(data_Beerling_NEE_m_day,'Color',color_model,'LineWidth',line_width)

xlim([0 366])
ylim([-15 5])
set(gca,'box','on')
h=gca; h.XAxis.TickLength = [0 0];
% text(axs(4).XLim(1), axs(4).YLim(2)*0.9, '(b-2)')
text(0.01, 0.85, '(b-2)','FontWeight','bold','Units','normalized','VerticalAlignment','bottom')

legend({'Exp Flux','T&C-SMEW'},'Location','southwest')
obs=data_maize.basalt_NEE.y_daily;
sim=data_Lewis_NEE_m_day;

R2 = stat_func{2}(sim,obs);
text(axs(4).XLim(2)*0.8, -4.5, sprintf('R^2=%.2f\nNSE=%.2f',... 
    R2(1,2),...
    stat_func{3}(sim,obs)),'FontSize',8)

%%% RE daily
axs(5) = nexttile;
% maize: plot daily values [gC/m2/day]
plot([1:366], data_maize.pre_Re.y_daily,'--','Color',color_obs,'LineWidth',line_width,'MarkerSize',5,...
    'MarkerFaceColor','none','MarkerEdgeColor',color_obs)
hold on
plot(data_pre_Re_m_day,'Color',color_model,'LineWidth',line_width)
ylabel('Re [gC m^{2} d^{-1}]')
xlim([0 366])
ylim([0 20])
set(gca,'box','on')
h=gca; h.XAxis.TickLength = [0 0];
% text(axs(5).XLim(1), axs(5).YLim(2)*0.95, '(c-1)')
text(0.01, 0.85, '(c-1)','FontWeight','bold','Units','normalized','VerticalAlignment','bottom')


obs=data_maize.pre_Re.y_daily;
sim=data_pre_Re_m_day;

R2 = stat_func{2}(sim,obs);
text(axs(5).XLim(2)*0.8, axs(5).YLim(2)*0.7, sprintf('R^2=%.2f\nNSE=%.2f',... 
    R2(1,2),...
    stat_func{3}(sim,obs)),'FontSize',8)

axs(6) = nexttile;
% soybean: plot daily values [gC/m2/day]
plot([1:366], data_maize.basalt_Re.y_daily,'--','Color',color_obs,'LineWidth',line_width,'MarkerSize',5,...
    'MarkerFaceColor','none','MarkerEdgeColor',color_obs)
hold on
plot(data_Lewis_Re_m_day,'Color',color_model,'LineWidth',line_width)
plot(data_Beerling_Re_m_day,'Color',color_model,'LineWidth',line_width)

xlim([0 366])
ylim([0 20])
set(gca,'box','on')
h=gca; h.XAxis.TickLength = [0 0];
% text(axs(6).XLim(1), axs(6).YLim(2)*0.95, '(c-2)')
text(0.01, 0.85, '(c-2)','FontWeight','bold','Units','normalized','VerticalAlignment','bottom')

obs=data_maize.basalt_Re.y_daily;
sim=data_Lewis_Re_m_day;

R2 = stat_func{2}(sim,obs);
text(axs(6).XLim(2)*0.8, axs(6).YLim(2)*0.7, sprintf('R^2=%.2f\nNSE=%.2f',... 
    R2(1,2),...
    stat_func{3}(sim,obs)),'FontSize',8)


xlabel(f8, 'Day of Year')

set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.01 0.01 0.98 0.98]);

fontsize(gcf,scale=2)

%% SI figure - time series ecosystem C flux control vs basalt 

% soy
figure;
f9 = tiledlayout(3, 2);

axs(1) = nexttile;
% maize: plot daily values [gC/m2/day]
plot([1:366], data_soy.control_GPP.y_daily,'--','Color',color_obs,'LineWidth',line_width,'MarkerSize',5,...
    'MarkerFaceColor','none','MarkerEdgeColor',color_obs)
hold on
plot(data_untr_GPP_s_day,'Color',color_model,'LineWidth',line_width)
title('Soybean - control')
xlim([0 366])
ylim([0 20])
set(gca,'box','on')
h=gca; h.XAxis.TickLength = [0 0];
ylabel('GPP [gC m^{2} d^{-1}]')
% text(axs(1).XLim(1), axs(1).YLim(2)*0.95, '(a-1)')
text(0.01, 0.85, '(a-1)','FontWeight','bold','Units','normalized','VerticalAlignment','bottom')

obs=data_soy.control_GPP.y_daily(1:365);
sim=data_untr_GPP_s_day;

R2 = stat_func{2}(sim,obs);
text(axs(1).XLim(2)*0.8, axs(1).YLim(2)*0.7, sprintf('R^2=%.2f\nNSE=%.2f',... 
    R2(1,2),...
    stat_func{3}(sim,obs)),'FontSize',8)


axs(2) = nexttile;
% soybean: plot daily values [gC/m2/day]
plot([1:366], data_soy.basalt_GPP.y_daily,'--','Color',color_obs,'LineWidth',line_width,'MarkerSize',5,...
    'MarkerFaceColor','none','MarkerEdgeColor',color_obs)
hold on
plot(data_Lewis_GPP_s_day,'Color',color_model,'LineWidth',line_width)
plot(data_Beerling_GPP_s_day,'Color',color_model,'LineWidth',line_width)

title('Soybean - basalt')
xlim([0 366])
ylim([0 20])
set(gca,'box','on')
h=gca; h.XAxis.TickLength = [0 0];
% text(axs(2).XLim(1), axs(2).YLim(2)*0.95, '(a-2)')
text(0.01, 0.85, '(a-2)','FontWeight','bold','Units','normalized','VerticalAlignment','bottom')

obs=data_soy.basalt_GPP.y_daily(1:365);
sim=data_Lewis_GPP_s_day;

R2 = stat_func{2}(sim,obs);
text(axs(2).XLim(2)*0.8, axs(2).YLim(2)*0.7, sprintf('R^2=%.2f\nNSE=%.2f',... 
    R2(1,2),...
    stat_func{3}(sim,obs)),'FontSize',8)

%%%% NEE daily
axs(3) = nexttile;
% soy: plot daily values [gC/m2/day]
plot([1:366], data_soy.control_NEE.y_daily,'--','Color',color_obs,'LineWidth',line_width,'MarkerSize',5,...
    'MarkerFaceColor','none','MarkerEdgeColor',color_obs)
hold on
plot(data_untr_NEE_s_day,'Color',color_model,'LineWidth',line_width)
ylabel('NEE [gC m^{2} d^{-1}]')
xlim([0 366])
ylim([-10 5])
set(gca,'box','on')
h=gca; h.XAxis.TickLength = [0 0];
% text(axs(3).XLim(1), axs(3).YLim(2)*0.9, '(b-1)')
text(0.01, 0.85, '(b-1)','FontWeight','bold','Units','normalized','VerticalAlignment','bottom')


obs=data_soy.control_NEE.y_daily(1:365);
sim=data_untr_NEE_s_day;

R2 = stat_func{2}(sim,obs);
text(axs(3).XLim(2)*0.8, -4.5, sprintf('R^2=%.2f\nNSE=%.2f',... 
    R2(1,2),...
    stat_func{3}(sim,obs)),'FontSize',8)

axs(4) = nexttile;
% soybean: plot daily values [gC/m2/day]
plot([1:366], data_soy.basalt_NEE.y_daily,'--','Color',color_obs,'LineWidth',line_width,'MarkerSize',5,...
    'MarkerFaceColor','none','MarkerEdgeColor',color_obs)
hold on
plot(data_Lewis_NEE_s_day,'Color',color_model,'LineWidth',line_width)
plot(data_Beerling_NEE_s_day,'Color',color_model,'LineWidth',line_width)

xlim([0 366])
ylim([-10 5])
set(gca,'box','on')
h=gca; h.XAxis.TickLength = [0 0];
% text(axs(4).XLim(1), axs(4).YLim(2)*0.9, '(b-2)')
text(0.01, 0.85, '(b-2)','FontWeight','bold','Units','normalized','VerticalAlignment','bottom')

obs=data_soy.basalt_NEE.y_daily(1:365);
sim=data_Lewis_NEE_s_day;

R2 = stat_func{2}(sim,obs);
text(axs(4).XLim(2)*0.8, -4.5, sprintf('R^2=%.2f\nNSE=%.2f',... 
    R2(1,2),...
    stat_func{3}(sim,obs)),'FontSize',8)

%%% RE daily
axs(5) = nexttile;
% soy: plot daily values [gC/m2/day]
plot([1:366], data_soy.pre_Re.y_daily,'--','Color',color_obs,'LineWidth',line_width,'MarkerSize',5,...
    'MarkerFaceColor','none','MarkerEdgeColor',color_obs)
hold on
plot(data_pre_Re_s_day,'Color',color_model,'LineWidth',line_width)
ylabel('Re [gC m^{2} d^{-1}]')
xlim([0 366])
ylim([0 20])
% text(axs(5).XLim(1), axs(5).YLim(2)*0.95, '(c-1)')
text(0.01, 0.85, '(c-1)','FontWeight','bold','Units','normalized','VerticalAlignment','bottom')
set(gca,'box','on')
h=gca; h.XAxis.TickLength = [0 0];
obs=data_soy.pre_Re.y_daily;
sim=data_pre_Re_s_day;

R2 = stat_func{2}(sim,obs);
text(axs(5).XLim(2)*0.8, axs(5).YLim(2)*0.7, sprintf('R^2=%.2f\nNSE=%.2f',... 
    R2(1,2),...
    stat_func{3}(sim,obs)),'FontSize',8)

axs(6) = nexttile;
% soybean: plot daily values [gC/m2/day]
plot([1:366], data_soy.basalt_Re.y_daily,'--','Color',color_obs,'LineWidth',line_width,'MarkerSize',5,...
    'MarkerFaceColor','none','MarkerEdgeColor',color_obs)
hold on
plot(data_Lewis_Re_s_day,'Color',color_model,'LineWidth',line_width)
plot(data_Beerling_Re_s_day,'Color',color_model,'LineWidth',line_width)

xlim([0 366])
ylim([0 20])
% text(axs(6).XLim(1), axs(6).YLim(2)*0.95, '(c-2)')
text(0.01, 0.85, '(c-2)','FontWeight','bold','Units','normalized','VerticalAlignment','bottom')
set(gca,'box','on')
h=gca; h.XAxis.TickLength = [0 0];

obs=data_soy.basalt_Re.y_daily(1:365);
sim=data_Lewis_Re_s_day;

R2 = stat_func{2}(sim,obs);
text(axs(6).XLim(2)*0.8, axs(6).YLim(2)*0.7, sprintf('R^2=%.2f\nNSE=%.2f',... 
    R2(1,2),...
    stat_func{3}(sim,obs)),'FontSize',8)


xlabel(f9, 'Day of Year')

set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.01 0.01 0.98 0.98]);
legend({'Exp Flux','T&C-SMEW'},'Location','north')
fontsize(gcf,scale=2)


%% Figure : ecosystem C flux in treatment years

line_width = 2;

stat_func = {
    @(sim,obs) (rmse(sim,obs));...
    @(sim,obs) (corrcoef(sim,obs).^2);...
    @(sim,obs) (1-sum((obs-sim).^2)/sum((obs-mean(obs)).^2));...
    };

maize_active = (data_maize.control_NEE.y_daily<0);
soy_active = (data_soy.control_NEE.y_daily<0);

figure;
f3 = tiledlayout(1, 2);

color_opt = {'#f6e8c3';'#d9d9d9';'#4292c6'};
categories = {'GPP','NEE','Re'};
xpos = [1 2 3];
bar_width = 0.4;

axs(1) = nexttile;
% maize: plot sum difference
bar(xpos-0.2, [mean(data_maize.control_GPP.y_daily(maize_active)),...
               mean(data_maize.control_NEE.y_daily(maize_active)),...
               mean(data_maize.control_Re.y_daily(maize_active))],...
              bar_width, 'FaceColor',color_opt{1})

hold on


scatter(xpos-0.2 , [mean(data_untr_GPP_m_day(maize_active)) ...
                mean(data_untr_NEE_m_day(maize_active)) ...
                mean(data_untr_Re_m_day(maize_active)) ],  120,'o', MarkerEdgeColor=color_opt{3},MarkerFaceColor='none', LineWidth=3)


bar(xpos+0.2, [mean(data_maize.basalt_GPP.y_daily(maize_active)),...
               mean(data_maize.basalt_NEE.y_daily(maize_active)),...
               mean(data_maize.basalt_Re.y_daily(maize_active))],...
              bar_width, 'FaceColor',color_opt{2})


scatter(xpos+0.15 , [mean(data_Lewis_GPP_m_day(maize_active)) ...
                mean(data_Lewis_NEE_m_day(maize_active)) ...
                mean(data_Lewis_Re_m_day(maize_active)) ],  120,'^', MarkerEdgeColor=color_opt{3},MarkerFaceColor=color_opt{3}, LineWidth=3)

scatter(xpos+0.25 , [mean(data_Beerling_GPP_m_day(maize_active)) ...
                mean(data_Beerling_NEE_m_day(maize_active)) ...
                mean(data_Beerling_Re_m_day(maize_active)) ],  120,'s', MarkerEdgeColor=color_opt{3},MarkerFaceColor=color_opt{3}, LineWidth=3)

ylim([-10 20])
xticks(xpos )
xticklabels(categories)
ylabel('[g m^{-2} d^{-1}]')
title('Maize (growing seasons)')
xlim([0 5])
text(0.01, 0.95, '(a)','FontWeight','bold','Units','normalized','VerticalAlignment','bottom')

axs(2) = nexttile;
bar(xpos-0.2, [mean(data_soy.control_GPP.y_daily(soy_active)),...
               mean(data_soy.control_NEE.y_daily(soy_active)),...
               mean(data_soy.control_Re.y_daily(soy_active))],...
               bar_width, 'FaceColor',color_opt{1})

hold on
scatter(xpos-0.2 , [mean(data_untr_GPP_s_day(soy_active)) ...
                mean(data_untr_NEE_s_day(soy_active)) ...
                mean(data_untr_Re_s_day(soy_active)) ],  120,'o', MarkerEdgeColor=color_opt{3},MarkerFaceColor='none', LineWidth=3)

bar(xpos+0.2, [mean(data_soy.basalt_GPP.y_daily(soy_active)),...
               mean(data_soy.basalt_NEE.y_daily(soy_active)),...
               mean(data_soy.basalt_Re.y_daily(soy_active))],...
               bar_width, 'FaceColor',color_opt{2})


scatter(xpos+0.15 , [mean(data_Lewis_GPP_s_day(soy_active)) ...
                mean(data_Lewis_NEE_s_day(soy_active)) ...
                mean(data_Lewis_Re_s_day(soy_active)) ],  120,'^', MarkerEdgeColor=color_opt{3},MarkerFaceColor=color_opt{3})

scatter(xpos+0.25 , [mean(data_Beerling_GPP_s_day(soy_active)) ...
                mean(data_Beerling_NEE_s_day(soy_active)) ...
                mean(data_Beerling_Re_s_day(soy_active)) ],  120,'s', MarkerEdgeColor=color_opt{3},MarkerFaceColor=color_opt{3})


ylim([-10 20])
xticks(xpos )
xticklabels(categories)
ylabel('[g m^{-2} d^{-1}]')
title('Soybean (growing seasons)')
xlim([0 5])
text(0.01, 0.95, '(b)','FontWeight','bold','Units','normalized','VerticalAlignment','bottom')
lgd=legend({'Exp Untr','Model Untr','Exp Basalt','Model Basalt_L','Model Basalt_B'},'Location','southeast',NumColumns=1);
% title(lgd,['no basalt   '   '      basalt'])

set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.01 0.01 0.98 0.98]);

fontsize(gcf,scale=2.8)



%% SI Figure: r_het + r_aut

% Kantola et al . Figure 3
figure;
f7 = tiledlayout(1, 2);

x_pos = [1 2];

axs(1) = nexttile;
Rhet_m_exp = [mean([1.16 0.78 1.51]) mean([1.24 1.40 1.35])];
Raut_m_exp = [mean([0.29 0.22 0.46]) mean([0.52 0.62 0.61])];

Rhet_m_exp_std = [sqrt([1.26-1.16 0.83-0.78 1.63-1.51].^2/3)...
                    sqrt([1.34-1.24 1.54-1.40 1.49-1.35].^2/3)];
Raut_m_exp_std = [sqrt([0.314-0.29 0.262-0.22 0.54-0.46].^2/3)...
    sqrt([0.579-0.52 0.73-0.62 0.718-0.61].^2/3)];

bar(x_pos, [Raut_m_exp(1) Rhet_m_exp(1)], bar_width, 'FaceColor',color_opt{1});
hold on
bar(x_pos + bar_width, [Raut_m_exp(2) Rhet_m_exp(2) ], bar_width, edgecolor='black',FaceColor=color_opt{2});

scatter(x_pos,[mean(data_untr_Raut_m_day)...
    mean(data_untr_Rhet_m_day)],  120,'o', MarkerEdgeColor=color_opt{3},MarkerFaceColor='none', LineWidth=3)

scatter(x_pos + bar_width-0.1,[mean(data_Lewis_Raut_m_day)...
    mean(data_Lewis_Rhet_m_day)],  120,'^', MarkerEdgeColor=color_opt{3},MarkerFaceColor=color_opt{3}, LineWidth=3)

scatter(x_pos + bar_width+0.1,[mean(data_Beerling_Raut_m_day)...
    mean(data_Beerling_Rhet_m_day)],  120,'s', MarkerEdgeColor=color_opt{3},MarkerFaceColor=color_opt{3}, LineWidth=3)

er = errorbar(x_pos ,[Raut_m_exp(1) Rhet_m_exp(1)],...
    [Raut_m_exp_std(1) Rhet_m_exp_std(1)],[Raut_m_exp_std(1) Rhet_m_exp_std(1)]);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';

er = errorbar(x_pos + bar_width,[Raut_m_exp(2) Rhet_m_exp(2) ],...
    [Raut_m_exp_std(2) Rhet_m_exp_std(2) ],[Raut_m_exp_std(2) Rhet_m_exp_std(2) ]);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';

xticks(x_pos+bar_width/2)
xticklabels({'Autotrophic (root)','Heterotrophic'})
ylabel('[gC m^{-2} d^{-1}]')
text(0.01, 0.95, '(a)','FontWeight','bold','Units','normalized','VerticalAlignment','bottom')
title('Maize')
set(gca,'box','on')
h=gca; h.XAxis.TickLength = [0 0];
legend('Exp Untr','Exp Basalt','Model Untr','Model Basalt_{L}','Model Basalt_{B}',...
    'Location','west')

axs(2) = nexttile;
Rhet_s_exp = [1.13 1.36];
Raut_s_exp = [0.34 0.61];

Rhet_s_exp_std = [1.225-1.13 1.49-1.36];
Raut_s_exp_std = [0.41-0.34 0.74-0.61];

hold on

bar(x_pos, [Raut_s_exp(1) Rhet_s_exp(1)], bar_width, 'FaceColor',color_opt{1});
hold on
bar(x_pos + bar_width, [Raut_s_exp(2) Rhet_s_exp(2) ], bar_width, edgecolor='black',FaceColor=color_opt{2});

scatter(x_pos,[mean(data_untr_Raut_s_day)...
    mean(data_untr_Rhet_s_day)],  120,'o', MarkerEdgeColor=color_opt{3},MarkerFaceColor='none', LineWidth=3)

scatter(x_pos + bar_width-0.1,[mean(data_Lewis_Raut_s_day)...
    mean(data_Lewis_Rhet_s_day)],  120,'^', MarkerEdgeColor=color_opt{3},MarkerFaceColor=color_opt{3}, LineWidth=3)
scatter(x_pos + bar_width+0.1,[mean(data_Beerling_Raut_s_day)...
    mean(data_Beerling_Rhet_s_day)],  120,'s', MarkerEdgeColor=color_opt{3},MarkerFaceColor=color_opt{3}, LineWidth=3)

er = errorbar(x_pos ,[Raut_s_exp(1) Rhet_s_exp(1)],...
    [Raut_s_exp_std(1) Rhet_s_exp_std(1)],[Raut_s_exp_std(1) Rhet_s_exp_std(1)]);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';

er = errorbar(x_pos + bar_width,[Raut_s_exp(2) Rhet_s_exp(2) ],...
    [Raut_s_exp_std(2) Rhet_s_exp_std(2) ],[Raut_s_exp_std(2) Rhet_s_exp_std(2) ]);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';

xticks(x_pos+bar_width/2)
xticklabels({'Autotrophic (root)','Heterotrophic'})
text(0.01, 0.95, '(b)','FontWeight','bold','Units','normalized','VerticalAlignment','bottom')
title('Soybean')
set(gca,'box','on')
h=gca; h.XAxis.TickLength = [0 0];

set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.01 0.01 0.98 0.98]);

fontsize(gcf,scale=2)


%% Figure : pH and base saturation comparison

cl = [
    216, 179, 101;...
    107,174,214;...
    49,130,189]./255;

color_TC_basalt_Lewis = '#6baed6';
color_TC_basalt_Beerling = '#08519c';
color_TC_untr = '#6baed6';
color_exp_0 = '#d8b365';%'#fee391';
color_exp_10 = '#8c510a';%'#cc4c02';

fig_colors = {color_TC_basalt_Lewis,...
            color_TC_basalt_Beerling,...
            color_TC_untr,...
            color_exp_0,...
            color_exp_10};

pH_exp_0 = data_exp.pH_Exp.basalt_0;
pH_exp_untr_0 = data_exp.pH_Exp.control_0;
pH_exp_10 = data_exp.pH_Exp.basalt_10;
pH_exp_untr_10 = data_exp.pH_Exp.control_10;

% SME
pH_exp_0_std = [0.0772	0.0551	0.04	0	0];
pH_exp_untr_0_std = [0.06	0.0735	0.0625	0.0772	0.0551];
pH_exp_10_std = [0.0699	0.0588	0.0515	0.05	0.05];
pH_exp_untr_10_std = [0.1066	0.0809	0.0772	0.0735	0.0735];

% Sim based on Lewis
pH_TC_Lewis = data_Lewis.ERW_biogeochem.pH;
pH_TC_untr = data_untr.ERW_biogeochem.pH;

CEC_TC_Lewis = data_Lewis.ERW_biogeochem.f_Ca+data_Lewis.ERW_biogeochem.f_Na+data_Lewis.ERW_biogeochem.f_K+data_Lewis.ERW_biogeochem.f_Mg;
CEC_TC_untr = data_untr.ERW_biogeochem.f_Ca+data_untr.ERW_biogeochem.f_Mg+data_untr.ERW_biogeochem.f_K+data_untr.ERW_biogeochem.f_Na;

% Sim based on Beerling
pH_TC_Beerling = data_Beerling.ERW_biogeochem.pH;
CEC_TC_Beerling = data_Beerling.ERW_biogeochem.f_Ca+data_Beerling.ERW_biogeochem.f_Na+data_Beerling.ERW_biogeochem.f_K+data_Beerling.ERW_biogeochem.f_Mg;


TT_TC_h = timetable(data_Lewis.TT_m.tt_h,pH_TC_Lewis,CEC_TC_Lewis*100,pH_TC_Beerling,CEC_TC_Beerling*100);
TT_TC = timetable(data_untr.TT_m.tt_h,pH_TC_untr,CEC_TC_untr*100);

TT_TC_h = retime(TT_TC_h,'daily','mean');
TT_TC = retime(TT_TC,'daily','mean');
TT_TC = synchronize(TT_TC_h,TT_TC,"daily","mean");
TT_TC = renamevars(TT_TC,1:width(TT_TC),["pH_TC_Lewis","CEC_TC_Lewis","pH_TC_Beerling","CEC_TC_Beerling","pH_TC_untr","CEC_TC_untr"]);

pH_TC_Lewis_std=[];
pH_TC_untr_std=[];
CEC_TC_Lewis_std=[];
CEC_TC_untr_std=[];
CEC_TC_Beerling_std=[];
pH_TC_Beerling_std=[];


TT_TC_year_pre = TT_TC((year(TT_TC.Time)==2016)&(month(TT_TC.Time)<11),:);
TT_TC_year_pre = retime(TT_TC_year_pre,'yearly','mean');
TT_TC_year = retime(TT_TC,'yearly','mean');

pH_TC_Lewis = [TT_TC_year_pre.pH_TC_Lewis;TT_TC_year.pH_TC_Lewis(2:end)];
pH_TC_untr = [TT_TC_year_pre.pH_TC_untr;TT_TC_year.pH_TC_untr(2:end)];

CEC_TC_Lewis = [TT_TC_year_pre.CEC_TC_Lewis; TT_TC_year.CEC_TC_Lewis(2:end)];
CEC_TC_untr = [TT_TC_year_pre.CEC_TC_untr; TT_TC_year.CEC_TC_untr(2:end)];

pH_TC_Beerling = [TT_TC_year_pre.pH_TC_Beerling; TT_TC_year.pH_TC_Beerling(2:end)];
CEC_TC_Beerling = [TT_TC_year_pre.CEC_TC_Beerling; TT_TC_year.CEC_TC_Beerling(2:end)];

ff=figure;
f4=tiledlayout(2, 2,'TileSpacing','compact');
% tiledlayout(2, 2,'TileSpacing','compact');

xpos = [1:5];
axs(1) = nexttile([1 1]);
title('No Basalt')
% legend('Model Untr','Exp Untr 0-10cm','Exp Untr 10-30cm','Location','best')
set(gca,'xcolor','w','ycolor','w','xtick',[],'ytick',[])
set(gca,'box','off')


axs(2) = nexttile;
title('With Basalt')
set(gca,'xcolor','w','ycolor','w','xtick',[],'ytick',[])
set(gca,'box','off')

CEC_exp_0 = data_exp.CEC_Exp.basalt_0;
CEC_exp_untr_0 = data_exp.CEC_Exp.control_0;
CEC_exp_10 = data_exp.CEC_Exp.basalt_10;
CEC_exp_untr_10 = data_exp.CEC_Exp.control_10;

% SME
CEC_exp_0_std = [1	0	1.4	0	0];
CEC_exp_untr_0_std = [1	1	1.2	1	0];
CEC_exp_10_std = [1.4 1 1 0 1.4];
CEC_exp_untr_10_std = [1.6	2.1 1.8 1 1.9];

axs(3) = nexttile;
% errorbar(xpos+0.1,CEC_TC_untr(1:5),CEC_TC_untr_std,'d','color',fig_colors{1},'MarkerSize',12,'MarkerFaceColor','w','LineWidth',line_width)
set(gca,'xcolor','w','ycolor','w','xtick',[],'ytick',[])
set(gca,'box','off')


axs(4) = nexttile;
set(gca,'xcolor','w','ycolor','w','xtick',[],'ytick',[])
set(gca,'box','off')

set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.01 0.01 0.98 0.98]);
fontsize(gcf,scale=2)

p = get(axs(1), 'Position');
h2 = axes('Parent', gcf, 'Position', [p(1) p(2)+0.11 p(3) .25]);
plot(xpos+0.1*xpos(1),pH_TC_untr(1:5),'d','color',fig_colors{3},'MarkerSize',12,'MarkerFaceColor','w','LineWidth',line_width)
hold on
errorbar(xpos-0.2*xpos(1),pH_exp_untr_0(1:5),pH_exp_untr_0_std(1:5),'o','color',fig_colors{4},'MarkerSize',12,'MarkerFaceColor','w','LineWidth',line_width)
errorbar(xpos-0.1*xpos(1),pH_exp_untr_10(1:5),pH_exp_untr_10_std(1:5),'o','color',fig_colors{5},'MarkerSize',12,'MarkerFaceColor','w','LineWidth',line_width)

xline(xpos(1:4)+0.5,'--')
xlim([0.5 5.5])
yticks([6:0.5:8])
ylim([5.5 8])
set(h2, 'XTick', []);
ylabel({'Soil pH';'(mean)'})
text(h2,0.01, 0.85, '(a-1)','FontWeight','bold','Units','normalized','VerticalAlignment','bottom')

fontsize(h2,scale=2)
legend('Model Untr','Exp Untr 0-10cm','Exp Untr 10-30cm','Location','best')

h1 = axes('Parent', gcf, 'Position', [p(1) p(2) p(3) .1]);
plot(TT_TC.Time,TT_TC.pH_TC_untr,'-','color',fig_colors{3},'LineWidth',line_width)
hold on
xline(TT_TC.Time(365:365:end-365),'--')

set(h1, 'Ylim', [5.5 6.4]);
yticks([5.5:0.5:6.4])
set(h1, 'XTick', []);
ylabel({'Soil pH';'(daily)'})
fontsize(h1,scale=2.5)
text(h1,0.01, 0.7, '(a-2)','FontSize',20,'FontWeight','bold','Units','normalized','VerticalAlignment','bottom')
legend('Model Untr','Location','best')


p = get(axs(2), 'Position');
h2 = axes('Parent', gcf, 'Position', [p(1) p(2)+0.11 p(3) .25]);
errorbar(xpos+0.1,pH_TC_Lewis(1:5),Lewis_pH_std,'^','color','#034e7b','MarkerSize',12,'MarkerFaceColor',fig_colors{1},'LineWidth',line_width)
hold on
errorbar(xpos+0.2,pH_TC_Beerling(1:5),Beerling_pH_std,'s','color',fig_colors{2},'MarkerSize',12,'MarkerFaceColor',fig_colors{2},'LineWidth',line_width)

errorbar(xpos-0.1,pH_exp_10(1:5),pH_exp_10_std(1:5),'o','color',fig_colors{5},'MarkerSize',12,'MarkerFaceColor',fig_colors{5},'LineWidth',line_width)
errorbar(xpos-0.2,pH_exp_0(1:5),pH_exp_0_std(1:5),'o','color',fig_colors{4},'MarkerSize',12,'MarkerFaceColor',fig_colors{4},'LineWidth',line_width)

xline(xpos(1:4)+0.5,'--')
xline(xpos(1)+0.35,':')
a = fill([xpos(1)+0.35 xpos(1)+0.5 xpos(1)+0.5 xpos(1)+0.35], [0 0 100 100], 'k','LineStyle','none');
a.FaceAlpha = 0.1;

xlim([0.5 5.5])
yticks([6:0.5:8])
ylim([5.5 8])
set(h2, 'XTick', []);
fontsize(h2,scale=2)
text(h2,0.01, 0.85, '(b-1)','FontSize',20,'FontWeight','bold','Units','normalized','VerticalAlignment','bottom')

h1 = axes('Parent', gcf, 'Position', [p(1) p(2) p(3) .1]);
plot(TT_TC.Time,TT_TC.pH_TC_Lewis,'-','color','#034e7b','LineWidth',line_width)
hold on
plot(TT_TC.Time,TT_TC.pH_TC_Beerling,'-','color',fig_colors{1},'LineWidth',line_width)
xline(TT_TC.Time(365:365:end-365),'--')
yticks([6:1:8])
set(h1, 'Ylim', [5.5 7.8]);
set(h1, 'XTick', []);
fontsize(h1,scale=2.5)
text(h1,0.01, 0.7, '(b-2)','FontSize',20,'FontWeight','bold','Units','normalized','VerticalAlignment','bottom')
xline(TT_TC.Time(315),':',{'initial','basalt'},'LabelHorizontalAlignment','center','LabelVerticalAlignment','middle')
legend({'Model Basalt_L','Model Basalt_B'},'Location','best','Orientation','horizontal')


p = get(axs(3), 'Position');
h2 = axes('Parent', gcf, 'Position', [p(1) p(2)+0.11 p(3) .25]);
plot(xpos+0.1,CEC_TC_untr(1:5),'d','color',fig_colors{3},'MarkerSize',12,'MarkerFaceColor','w','LineWidth',line_width)
hold on
errorbar(xpos-0.2,CEC_exp_untr_0(1:5),CEC_exp_untr_0_std(1:5),'o','color',fig_colors{4},'MarkerSize',12,'MarkerFaceColor','w','LineWidth',line_width)
errorbar(xpos-0.1,CEC_exp_untr_10(1:5),CEC_exp_untr_10_std(1:5),'o','color',fig_colors{5},'MarkerSize',12,'MarkerFaceColor','w','LineWidth',line_width)
xline(xpos(1:4)+0.5,'--')
xlim([0.5 5.5])
xticks([])
% xticklabels({[]})
ylabel({'Base Sat.', '[%] (mean)'})
ylim([65 100])
fontsize(h2,scale=2)
text(h2,0.01, 0.85, '(c-1)','FontSize',20,'FontWeight','bold','Units','normalized','VerticalAlignment','bottom')

h1 = axes('Parent', gcf, 'Position', [p(1) p(2) p(3) .1]);
plot(TT_TC.Time,TT_TC.CEC_TC_untr,'-','color',fig_colors{3},'LineWidth',line_width)
hold on
xline(TT_TC.Time(365:365:end-365),'--')

% set(h1, 'Ylim', [5.5 7.8]);
xticks(TT_TC.Time(180:365:end))
xticklabels({'2016','2017','2018','2019','2020'})
yticks([80 82])
h=gca; h.XAxis.TickLength = [0 0];
ylabel({'Base Sat.', '[%] (daily)'})
fontsize(h1,scale=2.5)
text(h1,0.01, 0.7, '(c-2)','FontSize',20,'FontWeight','bold','Units','normalized','VerticalAlignment','bottom')


p = get(axs(4), 'Position');

h2 = axes('Parent', gcf, 'Position', [p(1) p(2)+0.11 p(3) .25]);
errorbar(xpos+0.1,CEC_TC_Lewis(1:5),Lewis_CEC_std,'^','color','#034e7b','MarkerSize',12,'MarkerFaceColor',fig_colors{1},'LineWidth',line_width)
hold on
errorbar(xpos+0.2,CEC_TC_Beerling(1:5),Beerling_CEC_std,'s','color',fig_colors{2},'MarkerSize',12,'MarkerFaceColor',fig_colors{2},'LineWidth',line_width)
errorbar(xpos-0.2,CEC_exp_0(1:5),CEC_exp_0_std(1:5),'o','color',fig_colors{4},'MarkerSize',12,'MarkerFaceColor',fig_colors{4},'LineWidth',line_width)
errorbar(xpos-0.1,CEC_exp_10(1:5),CEC_exp_10_std(1:5),'o','color',fig_colors{5},'MarkerSize',12,'MarkerFaceColor',fig_colors{5},'LineWidth',line_width)
xline([xpos(1:4)+0.5],'--')
xline(xpos(1)+0.35 ,':')
a = fill([xpos(1)+0.35 xpos(1)+0.5 xpos(1)+0.5 xpos(1)+0.35], [0 0 100 100], 'k','LineStyle','none');
a.FaceAlpha = 0.1;
xlim([0.5 5.5])
xticks([1:5])
h=gca; h.XAxis.TickLength = [0 0];
ylim([65 100])
fontsize(h2,scale=2)
text(h2,0.01, 0.85, '(d-1)','FontSize',20,'FontWeight','bold','Units','normalized','VerticalAlignment','bottom')

legend({'Model Basalt_L','Model Basalt_B','Exp Basalt 0-10cm','Exp Basalt 10-30cm'},'Location','best')


h1 = axes('Parent', gcf, 'Position', [p(1) p(2) p(3) .1]);
plot(TT_TC.Time,TT_TC.CEC_TC_Lewis,'-','color','#034e7b','LineWidth',line_width)
hold on
plot(TT_TC.Time,TT_TC.CEC_TC_Beerling,'-','color',fig_colors{1},'LineWidth',line_width)
xline(TT_TC.Time(365:365:end-365),'--')

xticks(TT_TC.Time(180:365:end))
xticklabels({'2016','2017','2018','2019','2020'})
yticks([80 90])
h=gca; h.XAxis.TickLength = [0 0];
fontsize(h1,scale=2.5)
text(h1,0.01, 0.7, '(d-2)','FontSize',20,'FontWeight','bold','Units','normalized','VerticalAlignment','bottom')
xline(TT_TC.Time(315),':',{'initial','basalt'},'LabelHorizontalAlignment','center','LabelVerticalAlignment','middle')


%% Figure : Ca/Mg dissolution and CDR comparison

% Beerling et al.(2024) Supp excel 06
basalt_g_exp = [0:4]*5000; % g basalt/m2 per year in top layer
mix_g_exp = basalt_g_exp(2:end)./[0.023595168 0.046102539 0.067595642 0.088141529];
mix_g_exp = [0 mix_g_exp];

Ca_diss_TC_Lewis = (data_Lewis.ERW_weathering.EW*data_Lewis.ERW_const.min_st(:,1))*1e-3; % umol/m2/d -> mmol/m2/d
Mg_diss_TC_Lewis = (data_Lewis.ERW_weathering.EW*data_Lewis.ERW_const.min_st(:,2))*1e-3; % umol/m2/d -> mmol/m2/d
Ca_diss_TC_Beerling = (data_Beerling.ERW_weathering.EW*data_Beerling.ERW_const.min_st(:,1))*1e-3; % umol/m2/d -> mmol/m2/d
Mg_diss_TC_Beerling = (data_Beerling.ERW_weathering.EW*data_Beerling.ERW_const.min_st(:,2))*1e-3; % umol/m2/d -> mmol/m2/d

Ca_TC_Lewis_std=[];
Mg_TC_Lewis_std=[];
Ca_TC_Beerling_std=[];
Mg_TC_Beerling_std=[];


CaCO3_diss = cumsum(data_Beerling.ERW_TC_const.M_rock_add_CaCO3) - data_Beerling.ERW_biogeochem.CaCO3;
TT_TC_diss = timetable(data_Lewis.TT_m.tt_h+hours(1440),CaCO3_diss);
TT_TC_diss = retime(TT_TC_diss,'yearly','max');

TT_TC = timetable(data_Lewis.TT_m.tt_h,Ca_diss_TC_Lewis,Mg_diss_TC_Lewis,Ca_diss_TC_Beerling,Mg_diss_TC_Beerling);
TT_TC = retime(TT_TC,'yearly','mean'); % mmol/m2/d
TT_TC = renamevars(TT_TC,1:width(TT_TC),["Ca_diss_TC_Lewis","Mg_diss_TC_Lewis","Ca_diss_TC_Beerling","Mg_diss_TC_Beerling"]);

TT_TC.Ca_diss_TC_Lewis(2)=TT_TC.Ca_diss_TC_Lewis(2)+TT_TC.Ca_diss_TC_Lewis(1);
TT_TC.Ca_diss_TC_Lewis(1) = 0;
TT_TC.Mg_diss_TC_Lewis(2)=TT_TC.Mg_diss_TC_Lewis(2)+TT_TC.Mg_diss_TC_Lewis(1);
TT_TC.Mg_diss_TC_Lewis(1) = 0;
TT_TC.Ca_diss_TC_Beerling(2)=TT_TC.Ca_diss_TC_Beerling(2)+TT_TC.Ca_diss_TC_Beerling(1);
TT_TC.Ca_diss_TC_Beerling(1) = 0;
TT_TC.Mg_diss_TC_Beerling(2)=TT_TC.Mg_diss_TC_Beerling(2)+TT_TC.Mg_diss_TC_Beerling(1);
TT_TC.Mg_diss_TC_Beerling(1) = 0;

% Beerling et al.(2024) Supp excel 06
Ca_diss_exp = [0,429.93,1159.85,1791.11,1218.10]; % ug Ca/g mix
Mg_diss_exp = [0,237.25,156.23,343.32,540.11]; % ug Mg/g mix

Ca_diss_exp = Ca_diss_exp.*mix_g_exp*1e-3/40; % mmol/m2
Mg_diss_exp = Mg_diss_exp.*mix_g_exp*1e-3/24; % mmol/m2

% SEM
Ca_diss_exp_std = [0,1464.62/20, 1351.67/20,1334.53/21,1954.46/16].*mix_g_exp*1e-3/40;
Mg_diss_exp_std = [0,545.48/20,1241.83/20,1624.39/21,1391.67/16].*mix_g_exp*1e-3/24;


Ca_diss_TC_Lewis = cumsum(TT_TC.Ca_diss_TC_Lewis*365);
Mg_diss_TC_Lewis = cumsum(TT_TC.Mg_diss_TC_Lewis*365);
Ca_diss_TC_Beerling = cumsum((TT_TC.Ca_diss_TC_Beerling)*365)+TT_TC_diss.CaCO3_diss(1:5)*1e-3;
Mg_diss_TC_Beerling = cumsum(TT_TC.Mg_diss_TC_Beerling*365);


figure;
f5=tiledlayout(2, 2,'TileSpacing','compact');
% xpos = [2:5];
axs(1) = nexttile;
xpos = [0.7:1:3.7];
line_width = 2;
h3 = daboxplot(data3,'scatter',2,'whiskers',0,'boxalpha',0.7,...
    'xtlabels', condition_names(2:end),'mean',1,'color',cl(2:3,:),'scattersize',20,...
    'fill',0,'boxspacing',1.2,'scatteralpha',0.5,'whiskers',1,...
    'flipcolors',0,'outsymbol','kx'); 
for sc_n = 1:4
    h3.sc(sc_n,1).MarkerFaceColor = cl(2, :);
    h3.sc(sc_n,1).MarkerEdgeColor = cl(2, :);
    h3.bx(sc_n,1).LineWidth = 4;
    h3.ot(sc_n,1).LineWidth=5;
    h3.ot(sc_n,1).MarkerFaceColor = cl(2, :);
    h3.ot(sc_n,1).MarkerEdgeColor = cl(2, :);

    h3.sc(sc_n,2).MarkerFaceColor = cl(3, :);
    h3.sc(sc_n,2).MarkerEdgeColor = cl(3, :);
    h3.bx(sc_n,2).LineWidth = 4;
    h3.ot(sc_n,2).LineWidth=5;
    h3.ot(sc_n,2).MarkerFaceColor = cl(3, :);
    h3.ot(sc_n,2).MarkerEdgeColor = cl(3, :);

    h3.exp(sc_n,1) = errorbar(xpos(sc_n)-0.1,Ca_diss_exp(sc_n+1),Ca_diss_exp_std(sc_n+1),'o','color',cl(1,:),'MarkerSize',16,'MarkerFaceColor',cl(1,:),'LineWidth',line_width);
    % h3.b_L(sc_n,1) = plot(xpos(sc_n),Ca_diss_TC_Lewis(sc_n+1),'^','color','#034e7b','MarkerSize',16,'MarkerFaceColor',fig_colors{1},'LineWidth',line_width);
    h3.b_B(sc_n,1) = plot(xpos(sc_n)+0.3,Ca_diss_TC_Beerling(sc_n+1),'s','color',fig_colors{2},'MarkerSize',16,'MarkerFaceColor',fig_colors{2},'LineWidth',line_width);

end
xl = xlim; xlim([xl(1), xl(2)-0.4]);       % make space for the legend
xticks(xpos)
xticklabels({'2017','2018','2019','2020'})
title('Weathered Ca')
ylabel('[mmol m^{-2}]')
xline(xpos(1:3)+0.6,'--')
set(gca,'box','on')
h=gca; h.XAxis.TickLength = [0 0];
text(0.01, 0.9, '(a)','FontWeight','bold','Units','normalized','VerticalAlignment','bottom')

axs(2) = nexttile;
xpos = [0.7:1:3.7];
line_width = 2;
h3 = daboxplot(data4,'scatter',2,'whiskers',0,'boxalpha',0.7,...
    'xtlabels', condition_names(2:end),'mean',1,'color',cl(2:3,:),'scattersize',20,...
    'fill',0,'boxspacing',1.2,'scatteralpha',0.5,'whiskers',1,...
    'flipcolors',0,'outsymbol','kx'); 
for sc_n = 1:4
    h3.sc(sc_n,1).MarkerFaceColor = cl(2, :);
    h3.sc(sc_n,1).MarkerEdgeColor = cl(2, :);
    h3.bx(sc_n,1).LineWidth = 4;
    h3.ot(sc_n,1).LineWidth=5;
    h3.ot(sc_n,1).MarkerFaceColor = cl(2, :);
    h3.ot(sc_n,1).MarkerEdgeColor = cl(2, :);

    h3.sc(sc_n,2).MarkerFaceColor = cl(3, :);
    h3.sc(sc_n,2).MarkerEdgeColor = cl(3, :);
    h3.bx(sc_n,2).LineWidth = 4;
    h3.ot(sc_n,2).LineWidth=5;
    h3.ot(sc_n,2).MarkerFaceColor = cl(3, :);
    h3.ot(sc_n,2).MarkerEdgeColor = cl(3, :);

    h3.exp(sc_n,1) = errorbar(xpos(sc_n)-0.1,Mg_diss_exp(sc_n+1),Mg_diss_exp_std(sc_n+1),'o','color',cl(1,:),'MarkerSize',16,'MarkerFaceColor',cl(1,:),'LineWidth',line_width);
    % h3.b_L(sc_n,1) = plot(xpos(sc_n),Mg_diss_TC_Lewis(sc_n+1),'^','color','#034e7b','MarkerSize',16,'MarkerFaceColor',fig_colors{1},'LineWidth',line_width);
    h3.b_B(sc_n,1) = plot(xpos(sc_n)+0.3,Mg_diss_TC_Beerling(sc_n+1),'s','color',fig_colors{2},'MarkerSize',16,'MarkerFaceColor',fig_colors{2},'LineWidth',line_width);

end
xl = xlim; xlim([xl(1), xl(2)-0.4]); 
xticks(xpos)
xticklabels({'2017','2018','2019','2020'})
title('Weathered Mg')
xline(xpos(1:3)+0.6,'--')
% xlim([1.5 5.5])
% ylim([-30 45])
set(gca,'box','on')
h=gca; h.XAxis.TickLength = [0 0];
text(0.01, 0.9, '(b)','FontWeight','bold','Units','normalized','VerticalAlignment','bottom')
fontsize(gcf,scale=2)

% Beerling Supp Excel 06
CDR_pot_exp = [0 3.82 6.75 11.49 10.51]; % t/ha/yr

% % STD
% CDR_pot_exp_std = [0 8.00 11.69 14.59 15.03];

% SEM
CDR_pot_exp_std = [0 1.8261 2.5109 3.2609 3.7174];

num_avg_rel_Lewis = [cumsum(TT_TC.Ca_diss_TC_Lewis), cumsum(TT_TC.Mg_diss_TC_Lewis)]*365; % mmol/d -> mmol/yr
CO2_pot_num_Lewis = 1e-3*44*(2*(num_avg_rel_Lewis(:,1)+num_avg_rel_Lewis(:,2)))*1e-6*1e4; %(gCO2/m2) -> tCO2/ha

num_avg_rel_Beerling = [Ca_diss_TC_Beerling, Mg_diss_TC_Beerling];
CO2_pot_num_Beerling = 1e-3*44*(2*(num_avg_rel_Beerling(:,1)+num_avg_rel_Beerling(:,2)))*1e-6*1e4; %(gCO2/m2) -> tCO2/ha




xpos = [0.7:1:3.7];
line_width = 2;

axs(3) = nexttile([1 2]);
h2 = daboxplot(data2,'scatter',2,'whiskers',0,'boxalpha',0.7,...
    'xtlabels', condition_names(2:end),'mean',1,'color',cl(2:3,:),'scattersize',20,...
    'fill',0,'boxspacing',1.2,'scatteralpha',0.5,'whiskers',1,...
    'flipcolors',0,'outsymbol','kx'); 
for sc_n = 1:4
    h2.sc(sc_n,1).MarkerFaceColor = cl(2, :);
    h2.sc(sc_n,1).MarkerEdgeColor = cl(2, :);
    h2.bx(sc_n,1).LineWidth = 4;
    h2.ot(sc_n,1).LineWidth=5;
    h2.ot(sc_n,1).MarkerFaceColor = cl(2, :);
    h2.ot(sc_n,1).MarkerEdgeColor = cl(2, :);

    h2.sc(sc_n,2).MarkerFaceColor = cl(3, :);
    h2.sc(sc_n,2).MarkerEdgeColor = cl(3, :);
    h2.bx(sc_n,2).LineWidth = 4;
    h2.ot(sc_n,2).LineWidth=5;
    h2.ot(sc_n,2).MarkerFaceColor = cl(3, :);
    h2.ot(sc_n,2).MarkerEdgeColor = cl(3, :);

    h2.exp(sc_n,1) = errorbar(xpos(sc_n)-0.1,CDR_pot_exp(sc_n+1),CDR_pot_exp_std(sc_n+1),'o','color',cl(1,:),'MarkerSize',16,'MarkerFaceColor',cl(1,:),'LineWidth',line_width);
    % h2.b_L(sc_n,1) = plot(xpos(sc_n),CO2_pot_num_Lewis(sc_n+1),'^','color','#034e7b','MarkerSize',16,'MarkerFaceColor',fig_colors{1},'LineWidth',line_width);
    h2.b_B(sc_n,1) = plot(xpos(sc_n)+0.3,CO2_pot_num_Beerling(sc_n+1),'s','color',fig_colors{2},'MarkerSize',16,'MarkerFaceColor',fig_colors{2},'LineWidth',line_width);

end

ylabel('CDR');
xl = xlim; xlim([xl(1), xl(2)-0.2]); 

set(gca,'FontSize',20);
xline(xpos(1:3)+0.6,'--')
xticks(xpos)
xticklabels({'2017','2018','2019','2020'})
h=gca; h.XAxis.TickLength = [0 0];
title('Cumulative Potential CDR')
ylabel('[t CO_{2} ha^{-1}]')
text(0.01, 0.9, '(c)','FontWeight','bold','Units','normalized','VerticalAlignment','bottom','FontSize',18)

legend([h2.exp(1,:),h2.bx(1,1),h2.bx(1,2)],...
    [{'Exp','Basalt_{L,uc}','Basalt_{B,uc}'}],'Location','eastoutside'); 

set(gca,'box','on')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.01 0.01 0.98 0.98]);

fontsize(gcf,scale=1.5)

%% Figure : plant nutrient content comparison

St = data_untr.Stoich_L(1);
%%% Cal
Call = St.Call;  %%%  [gC/gCa]
Cals = St.Cals;  %%%  [gC/gCa]
Calr = St.Calr ; %%%  [gC/gCa]
Calc= Cals; %%%  [gC/gCa] Carbohydrate Reserve Carbon Calcium
Calf = St.Calf; %%%  [gC/gCa]
Calh = St.Calh;  %%%  [gC/gCa]
%%% Mag
Magl = St.Magl;  %%%  [gC/gMg]
Mags = St.Mags;  %%%  [gC/gMg]
Magr = St.Magr ; %%%  [gC/gMg]
Magc= Mags; %%%  [gC/gMg] Carbohydrate Reserve Carbon Magnesium
Magf = St.Magf; %%%  [gC/gMg]
Magh = St.Magh;  %%%  [gC/gMg]
%%% Sil
Sill = St.Sill;  %%%  [gC/gSi]
Sils = St.Sils;  %%%  [gC/gSi]
Silr = St.Silr ; %%%  [gC/gSi]
Silc= Sils; %%%  [gC/gSi] Carbohydrate Reserve Carbon Silicon 
Silf = St.Silf; %%%  [gC/gSi]
Silh = St.Silh;  %%%  [gC/gSi]

ind_growth = (data_untr.B_L(:,1,1)>0);
% Total Calcium Plant [gCa/m^2]
TCAL_tissue =  mean(data_untr.B_L(ind_growth,1,1)/Call + data_untr.B_L(ind_growth,1,3)/Calr + data_untr.B_L(ind_growth,1,4)/Calc )+mean(data_untr.Careserve_L(ind_growth,1));  
% Total Magnesium Plant [gMg/m^2]
TMAG_tissue =  mean(data_untr.B_L(ind_growth,1,1)/Magl + data_untr.B_L(ind_growth,1,3)/Magr + data_untr.B_L(ind_growth,1,4)/Magc)+mean(data_untr.Mgreserve_L(ind_growth,1)); 

values_p_untr = [TCAL_tissue, TMAG_tissue];

ind_growth = (data_Lewis.B_L(:,1,1)>0);
% Total Calcium Plant [gCa/m^2]
TCAL_tissue =  mean(data_Lewis.B_L(ind_growth,1,1)/Call + data_Lewis.B_L(ind_growth,1,3)/Calr + data_Lewis.B_L(ind_growth,1,4)/Calc )+mean(data_Lewis.Careserve_L(ind_growth,1));  
% Total Magnesium Plant [gMg/m^2]
TMAG_tissue =  mean(data_Lewis.B_L(ind_growth,1,1)/Magl + data_Lewis.B_L(ind_growth,1,3)/Magr + data_Lewis.B_L(ind_growth,1,4)/Magc )+mean(data_Lewis.Mgreserve_L(ind_growth,1)); 

values_p_Lewis = [TCAL_tissue, TMAG_tissue];

ind_growth = (data_Beerling.B_L(:,1,1)>0);
% Total Calcium Plant [gCa/m^2]
TCAL_tissue =  mean(data_Beerling.B_L(ind_growth,1,1)/Call + data_Beerling.B_L(ind_growth,1,3)/Calr + data_Beerling.B_L(ind_growth,1,4)/Calc )+mean(data_Beerling.Careserve_L(ind_growth,1));  
% Total Magnesium Plant [gMg/m^2]
TMAG_tissue =  mean(data_Beerling.B_L(ind_growth,1,1)/Magl + data_Beerling.B_L(ind_growth,1,3)/Magr + data_Beerling.B_L(ind_growth,1,4)/Magc )+mean(data_Beerling.Mgreserve_L(ind_growth,1)); 

values_p_Beerling = [TCAL_tissue, TMAG_tissue];


up_exp_untr = [3.72 2.86]; %g/m2 - maize Kelland Supp Excel
up_exp_untr_std = [1.08 0.99];

up_exp = [4.4 3.26]; %g/m2
up_exp_std = [1.39 1.11];

x_pos = [1:2];
bar_width = 0.4;

figure;
f6=tiledlayout(1, 2);
axs(1) = nexttile;
bar(x_pos, up_exp_untr, bar_width, 'FaceColor',color_opt{1});
hold on
bar(x_pos + bar_width, up_exp, bar_width, edgecolor='black',FaceColor=color_opt{2});

% TC-ERW
scatter(x_pos, values_p_untr,  120,'o', MarkerEdgeColor=fig_colors{3}, LineWidth=3)
scatter(x_pos + bar_width-0.05, values_p_Lewis,  120,'^', MarkerEdgeColor='#034e7b',MarkerFaceColor=fig_colors{1})
scatter(x_pos + bar_width+0.05, values_p_Beerling,  120,'s', MarkerEdgeColor=fig_colors{2},MarkerFaceColor=fig_colors{2})

er = errorbar(x_pos,up_exp_untr,up_exp_untr_std,up_exp_untr_std);    
er.Color = [0 0 0];                            
er.LineStyle = 'none'; 

er = errorbar(x_pos + bar_width,up_exp,up_exp_std,up_exp_std);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';

ylim([0 11])
xticks(x_pos + bar_width / 2)
xticklabels({'Ca','Mg'})
ylabel('[g m^{-2}]')
title('Maize')
text(0.01, 0.95, '(a)','FontWeight','bold','Units','normalized','VerticalAlignment','bottom')
set(gca,'box','on')
h=gca; h.XAxis.TickLength = [0 0];

legend('Exp Untr','Exp Basalt','Model Untr','Model Basalt_{L}','Model Basalt_{B}')

% text(axs(6).XLim(1), axs(6).YLim(2)*0.95, '(c)')

St = data_untr.Stoich_L(2);
%%% Cal
Call = St.Call;  %%%  [gC/gCa]
Cals = St.Cals;  %%%  [gC/gCa]
Calr = St.Calr ; %%%  [gC/gCa]
Calc= Cals; %%%  [gC/gCa] Carbohydrate Reserve Carbon Calcium
Calf = St.Calf; %%%  [gC/gCa]
Calh = St.Calh;  %%%  [gC/gCa]
%%% Mag
Magl = St.Magl;  %%%  [gC/gMg]
Mags = St.Mags;  %%%  [gC/gMg]
Magr = St.Magr ; %%%  [gC/gMg]
Magc= Mags; %%%  [gC/gMg] Carbohydrate Reserve Carbon Magnesium
Magf = St.Magf; %%%  [gC/gMg]
Magh = St.Magh;  %%%  [gC/gMg]
%%% Sil
Sill = St.Sill;  %%%  [gC/gSi]
Sils = St.Sils;  %%%  [gC/gSi]
Silr = St.Silr ; %%%  [gC/gSi]
Silc= Sils; %%%  [gC/gSi] Carbohydrate Reserve Carbon Silicon 
Silf = St.Silf; %%%  [gC/gSi]
Silh = St.Silh;  %%%  [gC/gSi]

ind_growth = (data_untr.B_L(:,2,1)>0);

% Total Calcium Plant [gCa/m^2]
TCAL_tissue =  mean(data_untr.B_L(ind_growth,2,1)/Call + data_untr.B_L(ind_growth,2,3)/Calr + data_untr.B_L(ind_growth,2,4)/Calc )+mean(data_untr.Careserve_L(ind_growth,2));  
% Total Magnesium Plant [gMg/m^2]
TMAG_tissue =  mean(data_untr.B_L(ind_growth,2,1)/Magl + data_untr.B_L(ind_growth,2,3)/Magr + data_untr.B_L(ind_growth,2,4)/Magc )+mean(data_untr.Mgreserve_L(ind_growth,2)); 

values_p_untr_s = [TCAL_tissue, TMAG_tissue];


ind_growth = (data_Lewis.B_L(:,2,1)>0);
% Total Calcium Plant [gCa/m^2]
TCAL_tissue =  mean(data_Lewis.B_L(ind_growth,2,1)/Call + data_Lewis.B_L(ind_growth,2,3)/Calr + data_Lewis.B_L(ind_growth,2,4)/Calc)+mean(data_Lewis.Careserve_L(ind_growth,2));  
% Total Magnesium Plant [gMg/m^2]
TMAG_tissue =  mean(data_Lewis.B_L(ind_growth,2,1)/Magl + data_Lewis.B_L(ind_growth,2,3)/Magr + data_Lewis.B_L(ind_growth,2,4)/Magc)+mean(data_Lewis.Mgreserve_L(ind_growth,2)); 

values_p_Lewis_s = [TCAL_tissue, TMAG_tissue];

ind_growth = (data_Beerling.B_L(:,2,1)>0);
% Total Calcium Plant [gCa/m^2]
TCAL_tissue =  mean(data_Beerling.B_L(ind_growth,2,1)/Call + data_Beerling.B_L(ind_growth,2,3)/Calr + data_Beerling.B_L(ind_growth,2,4)/Calc)+mean(data_Beerling.Careserve_L(ind_growth,2));  
% Total Magnesium Plant [gMg/m^2]
TMAG_tissue =  mean(data_Beerling.B_L(ind_growth,2,1)/Magl + data_Beerling.B_L(ind_growth,2,3)/Magr + data_Beerling.B_L(ind_growth,2,4)/Magc)+mean(data_Beerling.Mgreserve_L(ind_growth,2)); 

values_p_Beerling_s = [TCAL_tissue, TMAG_tissue];


up_exp_untr_s = [8.1 3.2]; %g/m2 - soybean Kelland Supp Excel
up_exp_untr_s_std = [1.19 0.65];

up_exp_s = [8.47 3.5]; %g/m2
up_exp_s_std = [1.53 1.32];

axs(2) = nexttile;
% TC-ERW
bar(x_pos, up_exp_untr_s, bar_width, 'FaceColor',color_opt{1});
hold on
scatter(x_pos, values_p_untr_s,  120,'o', MarkerEdgeColor=fig_colors{3}', LineWidth=3)

bar(x_pos + bar_width, up_exp_s, bar_width, edgecolor='black',FaceColor=color_opt{2});
scatter(x_pos + bar_width-0.05, values_p_Lewis_s,  120,'^', MarkerEdgeColor='#034e7b',MarkerFaceColor=fig_colors{1})
scatter(x_pos + bar_width+0.05, values_p_Beerling_s,  120,'s', MarkerEdgeColor=fig_colors{2},MarkerFaceColor=fig_colors{2})

er = errorbar(x_pos,up_exp_untr_s,up_exp_untr_s_std,up_exp_untr_s_std);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';

er = errorbar(x_pos + bar_width,up_exp_s,up_exp_s_std,up_exp_s_std);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';

ylim([0 11])
xticks(x_pos + bar_width / 2)
xticklabels({'Ca','Mg'})
text(0.01, 0.95, '(b)','FontWeight','bold','Units','normalized','VerticalAlignment','bottom')
title('Soybean')

set(gca,'box','on')
h=gca; h.XAxis.TickLength = [0 0];

set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.01 0.01 0.98 0.98]);
fontsize(gcf,scale=2.8)


