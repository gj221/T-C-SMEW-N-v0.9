clear all
close all
clc

current_directory = cd;

% simulation folder dir
SCENARIO_TEST = '\FINAL_mesocosm';

%% choose site
id_location = 'Kelland'; % choose case study: Dietzen, te_Pas, Kelland, Amann, Kantola 
NN_case_study_ERW_MOD = [1 2]; % No. of experiment setup in ERW_MOD sheet

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
ERW_PARAM_input = ERW_PARAM_input_t(ERW_PARAM_input_t.Number == NN_case_study_ERW_MOD(1),:);
ERW_PARAM_input_untr = ERW_PARAM_input_t(ERW_PARAM_input_t.Number == NN_case_study_ERW_MOD(2),:);

%% load data
fl_path = dir(fullfile(current_directory,'SIM',SCENARIO_TEST));

% check which data you need
aTable = struct2table(fl_path);
disp(aTable);

% basalt
cell_row = contains(aTable.name,ERW_PARAM_input.ExperimentalDescription,'IgnoreCase',true);
file_name = fullfile(fl_path(cell_row).folder,fl_path(cell_row).name);
data = load(file_name);

% control
cell_row_untr = contains(aTable.name,ERW_PARAM_input_untr.ExperimentalDescription,'IgnoreCase',true);
file_name_untr = fullfile(fl_path(cell_row_untr).folder,fl_path(cell_row_untr).name);
data_untr = load(file_name_untr);

%%
%--------------------------------------------------------------------
% ELEMENTAL RELEASES
%--------------------------------------------------------------------

t_end = data.MOD_PARAM_input.NumberOfDays;
t = data.t;
dt = data.dtw;

conv_mol = data.ERW_TC_const.conv_mol;

MM_Mg = 24/conv_mol; % magnesium
MM_Ca = 40/conv_mol; % calcium
MM_Na = 23/conv_mol; % sodium
MM_K = 39/conv_mol; % potassium
MM_Si = 28/conv_mol; % silicon 
MM_C = 12/conv_mol; % carbon
MM_Anions = 62/conv_mol; % average of associated anions
MM_Al = 27/conv_mol; % aluminium

min_st = data.ERW_const.min_st;
EW = data.ERW_weathering.EW;

%%%%%% TC %%%%%%%%%%%%%%%
%numerical elemental release
Ca_rel = EW*min_st(:, 1); %(mol-conv/ m2 d)
Mg_rel = EW*min_st(:, 2);
K_rel = EW*min_st(:, 3);
Na_rel = EW*min_st(:, 4);
Si_rel = EW*min_st(:, 6);

%Ca released by apatite (not numerically modeled)
M_apa = 0.03*10000; %(g/m2)
MM_apa = 422/conv_mol; %(g/mol-conv)
Ca_rel_apa = M_apa/MM_apa/t_end; %(mol-conv/ m2 d)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TC
%numerical averages and std (mmol/m2 d)
num_avg_rel = [mean(Ca_rel), mean(Ca_rel)+Ca_rel_apa, mean(Mg_rel), mean(K_rel), mean(Na_rel), mean(Si_rel)]*1e-3; %(mmol/ m2 d)
num_std_rel = [std(Ca_rel), std(Ca_rel), std(Mg_rel), std(K_rel), std(Na_rel), std(Si_rel)]*1e-3;

%experimental elemental release (mmol/m2 d)
exp_avg_rel = [18, 18, 4.2, 1.1, 0.33, 1.9];
exp_err_rel = [14, 14, 0.93, 1.4, 0.53, 0.4];


%--------------------------------------------------------------------
% CO2 SEQUESTRATION 
%--------------------------------------------------------------------
%EXPERIMENTS

%potential
CO2_pot_exp = 1e-3*44*(2*(exp_avg_rel(2)+exp_avg_rel(3))+(exp_avg_rel(4)+exp_avg_rel(5))); %(gCO2/m2d)
CO2_pot_exp_std = 1e-3*44*(2*(exp_err_rel(2)+exp_err_rel(3))+(exp_err_rel(4)+exp_err_rel(5)));

%effective (leaching)
Alk_L_exp = [-31, -0.35, 30, -0.62]; % (micromol/m2 d)
Alk_L_std_exp = [91, 5, 17, 11]; 
CO2_eff_exp_L = 44*(2*(Alk_L_exp(1)+Alk_L_exp(2))+(Alk_L_exp(3)+Alk_L_exp(4)))/conv_mol; %(gCO2/m2d)
CO2_eff_exp_L_std = 44*(2*(Alk_L_std_exp(1)+Alk_L_std_exp(2))+(Alk_L_std_exp(3)+Alk_L_std_exp(4)))/conv_mol; 
CO2_eff_exp = CO2_eff_exp_L;
CO2_eff_exp_std = CO2_eff_exp_L_std;

%NUMERICAL - TC

%potential
CO2_pot_num = 1e-3*44*(2*(num_avg_rel(2)+num_avg_rel(3))+(num_avg_rel(4)+num_avg_rel(5))); %(gCO2/m2d)
CO2_pot_num_std = 1e-3*44*(2*(num_std_rel(2)+num_std_rel(3))+(num_std_rel(4)+num_std_rel(5)));

%effective (leaching)
DIC_L_avg = mean((data.ERW_biogeochem.Leak_DIC-data_untr.ERW_biogeochem.Leak_DIC)); %(mol-conv/m2 d)
DIC_L_std = std((data.ERW_biogeochem.Leak_DIC-data_untr.ERW_biogeochem.Leak_DIC)); %(mol-conv/m2 d)
CO2_eff_num = 44*DIC_L_avg/conv_mol; %(gCO2/m2 d)
CO2_eff_num_std = 44*DIC_L_std/conv_mol; %(gCO2/m2d)


%figure for paper
figure;
f1 = tiledlayout(1, 3);
%------------------------------------------------------------------------------
% ELEMENT RELEASE
%------------------------------------------------------------------------------
axs(1) = nexttile;
categories = {'Ca','Mg', 'K', 'Na', 'Si'};
x_pos = [1:5];
ms = 16;
lw = 1.5;

errorbar(x_pos-0.1, exp_avg_rel(2:end), exp_err_rel(2:end), marker='s', linestyle='None', ...
                markersize=ms, markerfacecolor='white', markeredgecolor='black', color='black',linewidth=lw);
hold on
errorbar(x_pos+0.1, num_avg_rel(2:end), num_std_rel(2:end),  marker='v', linestyle='None', ...
                markersize=ms, markerfacecolor='#1f78b4', markeredgecolor='#1f78b4',color='#1f78b4',linewidth=lw);


xticks(axs(1),x_pos)
xticklabels(axs(1),categories)
ylabel('Element release [mmol m^{-2} d^{-1}]')
text(0.01, 0.95, '(a)','FontWeight','bold','Units','normalized','VerticalAlignment','bottom')
ylim([-0.5,50])

%-------------------------------------------------------------------
% pH
%-------------------------------------------------------------------
axs(2) = nexttile;
categories_pH = {'control','basalt'};

%experimental
pH_exp_avg = [6.56, 6.66];
pH_exp_std = [0.05, 0.05];
% TC
pH_end = [mean(data_untr.ERW_biogeochem.pH(100/dt:end)),mean(data.ERW_biogeochem.pH(100/dt:end))];
pH_std = [std(data_untr.ERW_biogeochem.pH(100/dt:end)),std(data.ERW_biogeochem.pH(100/dt:end))];

errorbar([1-0.1, 3-0.1], pH_exp_avg, pH_exp_std(1), marker='s', linestyle='None', ...
                markersize=ms, markerfacecolor='white', markeredgecolor='black', color='black',linewidth=lw)
hold on
errorbar([1+0.1, 3+0.1], pH_end, pH_std(2),  marker='v', linestyle='None', ...
                markersize=ms, markerfacecolor='#1f78b4', markeredgecolor='#1f78b4',color='#1f78b4',linewidth=lw)


xlim([0, 4])
xticks(axs(2),[1,3])
xticklabels(axs(2),categories_pH)
ylabel('final soil pH')
ylim([5.5, 8])
text(0.01, 0.95, '(b)','FontWeight','bold','Units','normalized','VerticalAlignment','bottom')

%---------------------------------------------
% CO2 sequestration
%---------------------------------------------
axs(3) = nexttile;
categories_CO2 = {'potential', 'effective'};

errorbar([1-0.1, 3-0.1], [CO2_pot_exp, CO2_eff_exp], [CO2_pot_exp_std, CO2_eff_exp_std], marker='s', linestyle='None', ...
                markersize=ms, markerfacecolor='white', markeredgecolor='black', color='black', linewidth=lw);
hold on
errorbar([1+0.1, 3+0.1], [CO2_pot_num, CO2_eff_num], [CO2_pot_num_std , CO2_eff_num_std],  marker='v', linestyle='None',... 
                markersize=ms, markerfacecolor='#1f78b4', markeredgecolor='#1f78b4', color='#1f78b4',linewidth=lw);


xlim([0, 4])
xticks(axs(3),[1,3])
xticklabels(axs(3),categories_CO2)
ylabel('CO_2 sequestration [g m^{-2} d^{-1}]')
legend({'Experiment','T&C-SMEW'},Location='northeast')
ylim([-0.25 3.5])
text(0.01, 0.95, '(c)','FontWeight','bold','Units','normalized','VerticalAlignment','bottom')

set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.01 0.01 0.98 0.98]);
fontsize(gcf,scale=2.8)

%%

%%%%%%%%%% TC
% untreated elemental 
%leachate 
K_L_untr_num = data_untr.LEAK_K;
Ca_L_untr_num = data_untr.LEAK_Ca;
Mg_L_untr_num = data_untr.LEAK_Mg;
Si_L_untr_num = data_untr.LEAK_Si;
values_L_untr = [sum(Ca_L_untr_num), sum(Mg_L_untr_num), sum(Si_L_untr_num)];
%in the column
Ca_untr_num = (data_untr.ERW_biogeochem.Ca_tot(end-1) + data_untr.ERW_biogeochem.CaCO3(end-1))*MM_Ca; %[g/m2]
Mg_untr_num = (data_untr.ERW_biogeochem.Mg_tot(end-1) + data_untr.ERW_biogeochem.MgCO3(end-1))*MM_Mg;
Si_untr_num = data_untr.ERW_biogeochem.Si_tot(end-1)*MM_Si; 
K_untr_num = data_untr.ERW_biogeochem.K_tot(end-1)*MM_K; 
values_untr = [ Ca_untr_num, Mg_untr_num, Si_untr_num];
%in the plant (and all sources of changes)
Ca_p_untr_num = sum(data_untr.ERW_biogeochem.dCa_tot)*MM_Ca; %[g/m2]
Mg_p_untr_num = sum(data_untr.ERW_biogeochem.dMg_tot)*MM_Mg;
Si_p_untr_num = sum(data_untr.ERW_biogeochem.dSi_tot)*MM_Si;
values_other_untr = -1*[Ca_p_untr_num, Mg_p_untr_num, Si_p_untr_num];


St = data.Stoich_L;
%%%
Kpotl = St.Kpotl;  %%%  [gC/gK]
Kpots = St.Kpots;  %%%  [gC/gK]
Kpotr = St.Kpotr ; %%%  [gC/gK]
Kpotc= Kpots; %%%  [gC/gK] Carbohydrate Reserve Carbon Potassium
Kpotf = St.Kpotf; %%%  [gC/gK]
Kpoth = St.Kpoth;  %%%  [gC/gK]
ftransL = St.ftransL;
ftransR = St.ftransR;

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

% Total Calcium Plant [gCa/m^2]
TCAL_tissue =  mean(data_untr.B_L(100:end,1,1)/Call + data_untr.B_L(100:end,1,3)/Calr + data_untr.B_L(100:end,1,4)/Calc + data_untr.B_L(100:end,1,5)/Calf);  
% Total Magnesium Plant [gMg/m^2]
TMAG_tissue =  mean(data_untr.B_L(100:end,1,1)/Magl + data_untr.B_L(100:end,1,3)/Magr + data_untr.B_L(100:end,1,4)/Magc + data_untr.B_L(100:end,1,5)/Magf); 
% Total Silicon Plant [gSi/m^2]
TSIL_tissue =  mean(data_untr.B_L(100:end,1,1)/Sill  + data_untr.B_L(100:end,1,3)/Silr + data_untr.B_L(100:end,1,4)/Silc + data_untr.B_L(100:end,1,5)/Silf); 

values_p_untr = [TCAL_tissue+mean(data_untr.Careserve_L(100:end,1)), TMAG_tissue+mean(data_untr.Mgreserve_L(100:end,1)), TSIL_tissue+mean(data_untr.Sireserve_L(100:end,1))];

%total
values_tot_untr = values_L_untr+values_untr+values_p_untr;

% basalt elemental
%leachate 
K_L_num = data.LEAK_K;
Ca_L_num = data.LEAK_Ca;
Mg_L_num = data.LEAK_Mg;
Si_L_num = data.LEAK_Si;

values_L = [sum(Ca_L_num), sum(Mg_L_num), sum(Si_L_num)];

%in the column
Ca_num = (data.ERW_biogeochem.Ca_tot(end-1) + data.ERW_biogeochem.CaCO3(end-1))*MM_Ca; %[g/m2]
Mg_num = (data.ERW_biogeochem.Mg_tot(end-1) + data.ERW_biogeochem.MgCO3(end-1))*MM_Mg; 
Si_num = data.ERW_biogeochem.Si_tot(end-1)*MM_Si ;
K_num = data.ERW_biogeochem.K_tot(end-1)*MM_K;
values = [ Ca_num, Mg_num, Si_num];
%in the plant
Ca_p_num = sum(data.ERW_biogeochem.dCa_tot)*MM_Ca; %[g/m2]
Mg_p_num = sum(data.ERW_biogeochem.dMg_tot)*MM_Mg;
Si_p_num = sum(data.ERW_biogeochem.dSi_tot)*MM_Si;
values_other = -1*([Ca_p_num, Mg_p_num, Si_p_num]);

% Total Calcium Plant [gCa/m^2]
TCAL_tissue =  mean(data.B_L(100:end,1,1)/Call + data.B_L(100:end,1,3)/Calr + data.B_L(100:end,1,4)/Calc + data.B_L(100:end,1,5)/Calf);  
% Total Magnesium Plant [gMg/m^2]
TMAG_tissue =  mean(data.B_L(100:end,1,1)/Magl + data.B_L(100:end,1,3)/Magr + data.B_L(100:end,1,4)/Magc + data.B_L(100:end,1,5)/Magf); 
% Total Silicon Plant [gSi/m^2]
TSIL_tissue =  mean(data.B_L(100:end,1,1)/Sill + data.B_L(100:end,1,3)/Silr + data.B_L(100:end,1,4)/Silc + data.B_L(100:end,1,5)/Silf); 

values_p = [TCAL_tissue+mean(data.Careserve_L(100:end,1)), TMAG_tissue+mean(data.Mgreserve_L(100:end,1)), TSIL_tissue+mean(data.Sireserve_L(100:end,1))];

%total
values_tot = values_L+values+values_p;%values+values_other;%


%experimental
conv = 1e-3/(80^2*pi*1e-6); % from mg/column to g/m2
data_exp = struct('element', ['Ca','Mg', 'Si'],...
        'L_untr', [150, 1.5, 15],... %leaching % [mg/column]
        'L_basalt', [150, 1.5, 15],...
        'untr', [54498, 1582, 95], ... %soil extraction
        'basalt', [56017, 1797, 117],...
        'p_untr', [270, 80, 400], ...  %plant
        'p_basalt', [300, 90, 500]);

Kelland_exp_std = struct('element', ['Ca','Mg', 'Si'],...
        'L_untr', [5.6, 0.1, 0.3],... %leaching % [mg/column]
        'L_basalt', [4.2, 0.05, 0.4],...
        'untr', [0.001, 38, 5], ... %soil extraction
        'basalt', [1300, 38, 1.8],...
        'p_untr', [5, 4.1, 19], ...  %plant
        'p_basalt', [0.1, 2.23, 17]);

Kelland_exp = data_exp;


%figure
figure;
f2 = tiledlayout(1, 4,'TileSpacing','compact');
categories = {'Ca','Mg', 'Si'};
% custom_colors = [0.4660 0.6740 0.1880;0.6350 0.0780 0.1840;0 0.4470 0.7410];
custom_colors = {'#f6e8c3','#d9d9d9','#1f78b4'};
bar_width = 0.4;
x_pos = [1:3];

%leachate panel
axs(4) = nexttile;

bar(x_pos, Kelland_exp.L_untr*conv, bar_width,'FaceColor',custom_colors{1});
hold on
er = errorbar(x_pos,Kelland_exp.L_untr*conv,Kelland_exp_std.L_untr*conv,Kelland_exp_std.L_untr*conv);    
er.Color = [0 0 0];                            
er.LineStyle = 'none'; 

scatter(x_pos, values_L_untr,  120,'o', MarkerEdgeColor=custom_colors{3}, LineWidth=3)

bar(x_pos + bar_width, Kelland_exp.L_basalt*conv, bar_width, EdgeColor='black',FaceColor=custom_colors{2});
er = errorbar(x_pos + bar_width,Kelland_exp.L_basalt*conv,Kelland_exp_std.L_basalt*conv,Kelland_exp_std.L_basalt*conv);    
er.Color = [0 0 0];                            
er.LineStyle = 'none'; 

% TC-ERW
scatter(x_pos + bar_width, values_L,  120,'o', MarkerEdgeColor='black',MarkerFaceColor=custom_colors{3}, LineWidth=1.2)

xlim([0.5 4])
xticks(x_pos + bar_width / 2)
xticklabels(categories)
ylabel('[g m^{-2}]')
% text(axs(4).XLim(1), axs(4).YLim(2)*0.95, '(a)')
text(0.01, 0.95, '(a)','FontWeight','bold','Units','normalized','VerticalAlignment','bottom')

title('cumulative leachate')

%column panel
axs(5) = nexttile;
bar(x_pos, Kelland_exp.untr*conv, bar_width, 'FaceColor',custom_colors{1});

hold on
er = errorbar(x_pos,Kelland_exp.untr*conv,Kelland_exp_std.untr*conv,Kelland_exp_std.untr*conv);    
er.Color = [0 0 0];                            
er.LineStyle = 'none'; 

bar(x_pos + bar_width, Kelland_exp.basalt*conv, bar_width, edgecolor='black',FaceColor=custom_colors{2});
er = errorbar(x_pos + bar_width,Kelland_exp.basalt*conv,Kelland_exp_std.basalt*conv,Kelland_exp_std.basalt*conv);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';

% TC-ERW
scatter(x_pos, values_untr,  120,'o', MarkerEdgeColor=custom_colors{3}, LineWidth=3)
scatter(x_pos + bar_width, values,  120,'o', MarkerEdgeColor='black',MarkerFaceColor=custom_colors{3}, LineWidth=1.2)

xlim([0.5 4])
xticks(x_pos + bar_width / 2) 
xticklabels(categories)
title('soil column')
text(0.01, 0.95, '(b)','FontWeight','bold','Units','normalized','VerticalAlignment','bottom')

%plant panel
axs(6) = nexttile;
bar(x_pos, Kelland_exp.p_untr*conv, bar_width, 'FaceColor',custom_colors{1});

hold on
er = errorbar(x_pos,Kelland_exp.p_untr*conv,Kelland_exp_std.p_untr*conv,Kelland_exp_std.p_untr*conv);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';

bar(x_pos + bar_width, Kelland_exp.p_basalt*conv, bar_width, edgecolor='black',FaceColor=custom_colors{2});
er = errorbar(x_pos + bar_width,Kelland_exp.p_basalt*conv,Kelland_exp_std.p_basalt*conv,Kelland_exp_std.p_basalt*conv);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';

% TC-ERW
scatter(x_pos, values_p_untr,  120,'o', MarkerEdgeColor=custom_colors{3}, LineWidth=3)
scatter(x_pos + bar_width, values_p,  120,'o', MarkerEdgeColor='black',MarkerFaceColor=custom_colors{3}', LineWidth=1.2)

xlim([0.5 4])
ylim([0 30])
xticks(x_pos + bar_width / 2)
xticklabels(categories)
title('plant')
text(0.01, 0.95, '(c)','FontWeight','bold','Units','normalized','VerticalAlignment','bottom')

%total
axs(7) = nexttile;
bar(x_pos, (Kelland_exp.L_untr+Kelland_exp.untr+Kelland_exp.p_untr)*conv, bar_width, 'FaceColor',custom_colors{1});

hold on


scatter(x_pos, values_tot_untr,...
               120,'o', MarkerEdgeColor=custom_colors{3}, LineWidth=3)

bar(x_pos + bar_width, (Kelland_exp.L_basalt+Kelland_exp.basalt+Kelland_exp.p_basalt)*conv, bar_width, edgecolor='black',FaceColor=custom_colors{2});


% TC-ERW
scatter(x_pos + bar_width, values_tot,...  
               120,'o', MarkerEdgeColor='black',MarkerFaceColor=custom_colors{3}, LineWidth=1.2)

er = errorbar(x_pos,(Kelland_exp.L_untr+Kelland_exp.untr+Kelland_exp.p_untr)*conv,sqrt(Kelland_exp_std.L_untr.^2+Kelland_exp_std.untr.^2+Kelland_exp_std.p_untr.^2)*conv,sqrt(Kelland_exp_std.L_untr.^2+Kelland_exp_std.untr.^2+Kelland_exp_std.p_untr.^2)*conv);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';

er = errorbar(x_pos + bar_width,(Kelland_exp.L_basalt+Kelland_exp.basalt+Kelland_exp.p_basalt)*conv,sqrt(Kelland_exp_std.L_basalt.^2+Kelland_exp_std.basalt.^2+Kelland_exp_std.p_basalt.^2)*conv,sqrt(Kelland_exp_std.L_basalt.^2+Kelland_exp_std.basalt.^2+Kelland_exp_std.p_basalt.^2)*conv);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';

xlim([0.5 4])
xticks(x_pos + bar_width / 2)
xticklabels(categories)
text(0.01, 0.95, '(d)','FontWeight','bold','Units','normalized','VerticalAlignment','bottom')

lgd = legend({'Experiment','T&C-SMEW','Experiment','T&C-SMEW'},'NumColumns',2,'Location','northoutside');
title(lgd,['no basalt   '   '      basalt'])

title('total')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.01 0.01 0.98 0.98]);
fontsize(gcf,scale=2.5)
