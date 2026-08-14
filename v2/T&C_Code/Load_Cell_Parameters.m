function PAR = Load_Cell_Parameters(PARAM_IC, PFT_OVERRIDE, ctx)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load_Cell_Parameters
%%
%% Run one cell's parameter script in an isolated function workspace and
%% return every variable it defined as a single struct.
%%
%% The T&C parameter files are scripts that write ~200 loose variables
%% into whatever workspace calls them.  That is fine for a single-column
%% run but unusable when each sub-plot cell needs its own parameter set,
%% because the cells would overwrite one another.  Running the script
%% inside this function confines it to a private workspace, so the caller
%% gets N independent parameter structs and cells can differ in PFT,
%% management, rooting depth, soil texture and initial state.
%%
%% Inputs
%%   PARAM_IC      full path to the cell's MOD_PARAM_*.m parameter script
%%   PFT_OVERRIDE  optional full path to a PFT override script (e.g. the
%%                 prairie margin file) run AFTER PARAM_IC so it replaces
%%                 only the vegetation parameters and inherits soil,
%%                 snow and interception settings.  '' or [] to skip.
%%   ctx           struct of variables the parameter scripts expect to
%%                 already exist in the workspace (MOD_PARAM_input,
%%                 ERW_PARAM_input, NN, NNd, ms, Datam, Ca, Lat, Lon,
%%                 DeltaGMT, t_bef, t_aft, OPT_SoilBiogeochemistry,
%%                 NN_case_study_ERW_MOD, B_out and friends).  Every
%%                 field is unpacked into the private workspace first.
%%
%% Output
%%   PAR   struct with one field per variable the scripts defined.
%%
%% Note on Restating_parameters
%%   It is run here as well, so PAR carries the assembled Soil_Param,
%%   Interc_Param, SnowIce_Param, VegH_Param and VegL_Param structs that
%%   HYDROLOGIC_UNIT expects, already consistent with this cell's soil.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2; PFT_OVERRIDE = ''; end
if nargin < 3; ctx = struct(); end

%% Unpack the caller-supplied context into this private workspace.
ctx_names = fieldnames(ctx);
for ii = 1:numel(ctx_names)
    eval([ctx_names{ii} ' = ctx.(ctx_names{ii});']);
end
clear ii ctx_names

if ~isfile(PARAM_IC)
    error('Load_Cell_Parameters:MissingParamFile', ...
        'Parameter file not found: %s', PARAM_IC);
end
run(PARAM_IC);

%% Optional PFT override, applied on top of the base parameter set.
if ~isempty(PFT_OVERRIDE)
    if ~isfile(PFT_OVERRIDE)
        error('Load_Cell_Parameters:MissingPFTOverride', ...
            'PFT override file not found: %s', PFT_OVERRIDE);
    end
    run(PFT_OVERRIDE);
end

%% Assemble the parameter structs HYDROLOGIC_UNIT consumes.
Restating_parameters;

%% Package everything this workspace now holds.
vars = who;
drop = {'PARAM_IC','PFT_OVERRIDE','ctx','vars','drop','PAR','ii','ctx_names'};
PAR = struct();
for ii = 1:numel(vars)
    if any(strcmp(vars{ii}, drop)); continue; end
    PAR.(vars{ii}) = eval(vars{ii});
end

end
