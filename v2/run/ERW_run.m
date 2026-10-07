%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% WORKING LAUNCH PAD HBM  %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% INPUT MANAGER 

clear all;
close all;
clc;

%%%%%% This launcher lives in v2/run; the repository root is one level up.
script_directory = fileparts(mfilename('fullpath'));
if isempty(script_directory); script_directory = cd; end
current_directory = fileparts(script_directory);
Directory = fullfile(current_directory,'T&C_Code');
Input_Directory = fullfile(Directory,'Inputs');
Parameter_Directory = fullfile(Directory,'Parameters');
addpath(Directory,Input_Directory,Parameter_Directory)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SCENARIO_TEST = 'FINAL_mesocosm';
%%%%%%%%%%%%%%%%%%%%%
id_location = 'Kantola';
NN_case_study_ERW_MOD = 5; % Kantola: 3 Beerling, 4 Lewis, 5 control, 6 pretreatment, 7 calcite
save_results = 1; % 1 - save results; 0 - no
save_fig = 1;
OPT_ERW = 1; % ERW module on or off
OPT_BG_spinup = 1; % p pool values based on spin-up
OPT_SoilBiogeochemistry = 1; % soil biogeochemistry module
%%%%%%%%%%%%%%%%%%%%%
%%% optional overrides for batch runs (environment variables survive the clear all above)
if ~isempty(getenv('TC_CASE')); NN_case_study_ERW_MOD = str2double(getenv('TC_CASE')); end
if ~isempty(getenv('TC_SCENARIO')); SCENARIO_TEST = getenv('TC_SCENARIO'); end
%%%%%%%%%%%%%%%%%%%%%
workbook = fullfile(Input_Directory,'ERW_MOD_PARAM_parameters.xlsx');
MOD_PARAM_opts = detectImportOptions(workbook,"Sheet",'MOD_PARAM');
MOD_PARAM_opts.DataRange = "A2:Y6";
MOD_PARAM_opts.VariableTypes = ["double","string","string","double","double",...
    "string","double","double","string","double",...
    "double","double","double","double","string", ...
    "double","double","string","double","double",...
    "double","double","double","string","string"];

MOD_PARAM_input = readtable(workbook,MOD_PARAM_opts);

%%%%%%%%%%%%%% read ERW module parameters fromexcel
ERW_PARAM_opts = detectImportOptions(workbook,"Sheet",'ERW_MOD');
ERW_PARAM_opts.DataRange ="A2:AJ7";
ERW_PARAM_opts.VariableTypes ={'double','string','string','double','double',...
  'string','double','double','string','double',...
  'string','string','string','string','double',...
  'double','double','double','double','double',...
  'double','double','double','string','double',...
  'string','string','double','double','double',...
  'double','double','double','double','double',...
  'double'};
ERW_PARAM_opts = setvaropts(ERW_PARAM_opts,{'f_CEC_in','rock_filepath','mineral','rock_f_in'},'FillValue','NaN');

ERW_PARAM_input = readtable(workbook, ERW_PARAM_opts);

% find the corresponding case study
ERW_PARAM_input = ERW_PARAM_input(ERW_PARAM_input.Number == NN_case_study_ERW_MOD,:);
MOD_PARAM_input=MOD_PARAM_input(MOD_PARAM_input.Number == NN_case_study_ERW_MOD,:);

% using spin-up data if any
if OPT_BG_spinup == 1
    load(fullfile(Input_Directory,['Data_' id_location '_BG.mat']))
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%
NN=MOD_PARAM_input.NumberOfDays*24;% %%% time Step [hrs]
%%%%%%%%% Time step 
dt=3600; %%[s] %%% 
dth=1; %%[h]
%%%%%%%%%%%%
eval(append('Zs','=',eval('MOD_PARAM_input.Zs'),';'))
ms = length(Zs)-1; %MOD_PARAM_input.ms; %%% Soil Layer 
cc = MOD_PARAM_input.cc; %% Crown area
%%%%%%%%% METEO INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(fullfile(Input_Directory,['Data_' id_location '.mat']))
% Date=D; clear D 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
if NN_case_study_ERW_MOD == 6 % field pretreatment
    x1 = 7; % starting from 01-Jan-2009
    NN = 4382*24;
    x2 = x1+NN-1;
    OPT_ERW = 0;
else
    x1 = 61351; % starting from 01-Jan-2016
    NN = 1826*24;
    x2 = x1+NN-1;
end

%%% optional short-run override for smoke tests
if ~isempty(getenv('TC_NDAYS')); NN = str2double(getenv('TC_NDAYS'))*24; x2 = x1+NN-1; end
%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%
Date=Date(x1:x2);
Pr=Pr(x1:x2);
Ta=Ta(x1:x2);
Ws=Ws(x1:x2); ea=ea(x1:x2);  SAD1=SAD1(x1:x2);
SAD2=SAD2(x1:x2); SAB1=SAB1(x1:x2); Pre=Pre(x1:x2);
SAB2=SAB2(x1:x2); N=N(x1:x2); Tdew=Tdew(x1:x2);esat=esat(x1:x2);
PARB=PARB(x1:x2); PARD = PARD(x1:x2);
%%%%%%%%%%%
t_bef= -0.67; t_aft= 1.67;
%%%%%%%%%%%%%%%%%%%%
Ds=esat-ea; %% [Pa] Vapor Pressure Deficit
Ds(Ds<0)=0;
%%%%%%%%%%%%%%%%%%%%%%%%%
load(fullfile(Input_Directory,'Ca_Data.mat'));
d1 = find(abs(Date_CO2-Date(1))<1/36);
d2 = find(abs(Date_CO2-Date(end))<1/36);
Ca=Ca(d1:d2); 

clear d1 d2 Date_CO2 
Oa= 210000;% Intercellular Partial Pressure Oxygen [umolO2/mol] -
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ws(Ws<=0)=0.01;
%%dt,Pr(i),Ta(i),Ws(i),ea(i),Pre(i),Rdir(i),Rdif(i),N(i),z,Tdew(i),esat(i),.
[YE,MO,DA,HO,MI,SE] = datevec(Date);
Datam(:,1) = YE; Datam(:,2)= MO; Datam(:,3)= DA; Datam(:,4)= HO;
clear YE MO DA HO MI SE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PARAM_IC = fullfile(Parameter_Directory,[char(MOD_PARAM_input.MOD_PARAM) '.m']);
ERW_PARAM_IC = fullfile(Parameter_Directory,'ERW_MOD_PARAM_excel.m');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd(Directory)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MAIN_FRAME_BG_ERW_10mins;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
cd(current_directory);

if save_results == 1
    % data dir
    fl_save = fullfile(current_directory,'SIM',SCENARIO_TEST);

    if ~exist(fl_save,'dir')
        mkdir(fl_save)
    end

    if strcmp(SCENARIO_TEST,'EQUIV_1D_TMP')
        if NN_case_study_ERW_MOD ~= 5
            error('ERW_run:EquivalenceCase', ...
                'The equivalence baseline requires Kantola control case 5.');
        end
        equivalence_fingerprint = ERW_Equivalence_Fingerprint( ...
            current_directory,NN/24);
    end

    if OPT_ERW == 1
        result_file = fullfile(fl_save, ...
            ['TC_ERW_results_' id_location '_' ...
            char(ERW_PARAM_input.ExperimentalDescription) '.mat']);
    else
        result_file = fullfile(fl_save, ...
            ['TC_ERW_results_' id_location '_ERW_off.mat']);
    end
    temporary_result = [tempname(fl_save) '.mat'];
    save(temporary_result)
    [move_ok,move_message] = movefile(temporary_result,result_file,'f');
    if ~move_ok
        if isfile(temporary_result); delete(temporary_result); end
        error('ERW_run:ResultMoveFailed','%s',move_message);
    end

end
