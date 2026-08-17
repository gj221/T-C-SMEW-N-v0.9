function Capture_2D(base_file,out)
%% Capture ERW_function inputs from the 2-D driver at N_col = 1.
%%
%% Replays the saved single-column workspace through MAIN_FRAME_BG_ERW_2D on
%% a one-cell topology, exactly as test_2d_driver_equivalence does, so the
%% two captures are directly comparable call for call.
%%
%% Arm with TC_ERW_CAPTURE=1 in the environment before calling.

here = fileparts(mfilename('fullpath'));
code = fullfile(here,'T&C_Code');
addpath(code,fullfile(code,'Inputs'),fullfile(code,'Parameters'));

R1 = load(base_file);
fprintf('[Capture_2D] baseline %s (%d h)\n',base_file,R1.NN);
run_driver_single(R1,code);

global ERW_CAP %#ok<GVMIS>
cap = ERW_CAP;
cap.vals = cap.vals(1:min(cap.k,cap.max),:);
cap.source = '2D MAIN_FRAME_BG_ERW_2D at N_col = 1';
cap.init = ERW_CAPI_get();
fprintf('[Capture_2D] %d ERW_function calls, %d scalar terms\n', ...
    cap.k,numel(cap.names));
save(out,'cap','-v7.3');
fprintf('[Capture_2D] saved %s\n',out);
end

%% =======================================================================
function run_driver_single(R1,code)
%% Rebuild the caller context the 2-D driver expects, then run it. Plain
%% local function so its workspace is not static and the driver script may
%% declare its own variables here.
NN = R1.NN; NNd = R1.NNd; ms = R1.ms; cc = R1.cc; Zs = R1.Zs; %#ok<NASGU>
dt = R1.dt; dth = R1.dth;                                      %#ok<NASGU>
Datam = R1.Datam; Date = R1.Date;                              %#ok<NASGU>
Pr=R1.Pr; Ta=R1.Ta; Ds=R1.Ds; Ws=R1.Ws; ea=R1.ea; N=R1.N; Pre=R1.Pre; %#ok<NASGU>
Tdew=R1.Tdew; SAB1=R1.SAB1; SAB2=R1.SAB2; SAD1=R1.SAD1; SAD2=R1.SAD2;  %#ok<NASGU>
PARB=R1.PARB; PARD=R1.PARD; Ca=R1.Ca; Oa=R1.Oa;                %#ok<NASGU>
Lat=R1.Lat; Lon=R1.Lon; DeltaGMT=R1.DeltaGMT;                  %#ok<NASGU>
t_bef=R1.t_bef; t_aft=R1.t_aft;                                %#ok<NASGU>
id_location = R1.id_location;                                  %#ok<NASGU>
NN_case_study_ERW_MOD = R1.NN_case_study_ERW_MOD;              %#ok<NASGU>
MOD_PARAM_input = R1.MOD_PARAM_input;
ERW_PARAM_input = R1.ERW_PARAM_input;
OPT_ERW = R1.OPT_ERW;                                          %#ok<NASGU>
OPT_SoilBiogeochemistry = R1.OPT_SoilBiogeochemistry;          %#ok<NASGU>
if isfield(R1,'B_out');      B_out = R1.B_out;           end   %#ok<NASGU>
if isfield(R1,'O_out');      O_out = R1.O_out;           end   %#ok<NASGU>
if isfield(R1,'NavlI_out');  NavlI_out = R1.NavlI_out;   end   %#ok<NASGU>
if isfield(R1,'RexmyI_out'); RexmyI_out = R1.RexmyI_out; end   %#ok<NASGU>
if isfield(R1,'NupI_L_out'); NupI_L_out = R1.NupI_L_out; end   %#ok<NASGU>
if isfield(R1,'IrD');        IrD = R1.IrD;               end   %#ok<NASGU>

Parameter_Directory = fullfile(code,'Parameters');
topo = Build_Single_Cell_Topology(3600);                       %#ok<NASGU>
ColCfg = struct( ...
    'PARAM_IC', fullfile(Parameter_Directory, [char(MOD_PARAM_input.MOD_PARAM) '.m']), ...
    'PFT_OVERRIDE', '', ...
    ... %% Mirror the single-column frame's own predicate. It applies rock
    ... %% when diss_f > 0; M_rock_in is the mesocosm's initial rock mass and
    ... %% is 0 for the Kantola field setup, which instead adds 5000 g/m2 on
    ... %% the time_application schedule. Testing M_rock_in > 0 therefore
    ... %% reports no rock for basalt and calcite, leaving cfg.apply_rock
    ... %% false while the chemistry row still carries a feedstock.
    'apply_rock', logical(R1.OPT_ERW && ERW_PARAM_input.diss_f > 0)); %#ok<NASGU>

old = cd(code);
cleanup = onCleanup(@() cd(old)); %#ok<NASGU>
MAIN_FRAME_BG_ERW_2D;
end
