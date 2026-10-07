function Run_Field_2D_Task
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run one Slurm array element of the five-year Kantola runs 1 to 9 on a
%% field grid: T&C-SMEW-N v1.0's polygon grid on the real field, its outer
%% band the margin, with or without topography from the DEM. The 2-D
%% counterpart of Run_Base_2D_Task for a grid saved by
%% Build_Field_Grid_Topology, named by TC_FIELD_TOPOLOGY.
%%
%% The grid file lives beside the deployed root, never inside it, and its
%% content enters the validation stamp, so a grid replaced after validation
%% stops production with StaleValidation. Its file name enters the campaign
%% tag, so a flat and a topographic grid, or two cell widths, each get a
%% results tree of their own.
%%
%% The chemistry is untouched: every calibrated number comes from
%% Campaign_Parameters in the five-year 1-D campaign directory, as for the
%% 12-column plot, and the provenance assertion is the same.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

one_d_tag = ...
    'Base_1D';
campaign_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(campaign_dir));
one_d_dir = fullfile(repo_root,'runtime',one_d_tag);
run_dir = fullfile(repo_root,'v2','run');
code_dir = fullfile(repo_root,'v2','T&C_Code');
addpath(run_dir,code_dir,fullfile(code_dir,'Inputs'), ...
    fullfile(code_dir,'Parameters'),campaign_dir,one_d_dir);

matlab_release = version('-release');
if ~strcmp(matlab_release,'2025a')
    error('Bouchet_Kantola_Field2D:WrongMatlabRelease', ...
        'This campaign requires MATLAB R2025a; loaded R%s.',matlab_release);
end
Bouchet_Single_Thread_Runtime('Bouchet_Kantola_Field2D');
P = Campaign_Parameters;

[topology_file,topology_name] = read_topology_file(repo_root);
campaign_tag = ['Field_2D__' topology_name];
topo = load_topology(topology_file);
n_columns = topo.N_cells;

required_files = { ...
    fullfile(code_dir,'Inputs','ERW_MOD_PARAM_parameters.xlsx'), ...
    fullfile(code_dir,'Inputs','Data_Kantola.mat'), ...
    fullfile(code_dir,'Inputs','Data_Kantola_BG.mat'), ...
    fullfile(code_dir,'Inputs','Ca_Data.mat'), ...
    fullfile(code_dir,'Inputs','All_deposition_data.mat'), ...
    fullfile(code_dir,'Inputs','blueridge_CDF_data.json'), ...
    fullfile(repo_root,'v2','SIM','EQUIV_1D_TMP', ...
        'TC_ERW_results_Kantola_crop_control.mat'), ...
    fullfile(one_d_dir,'Campaign_Parameters.m'), ...
    fullfile(one_d_dir,'Assert_Base_Provenance.m'), ...
    fullfile(one_d_dir,'Fingerprint_Digest.m'), ...
    topology_file};
for k = 1:numel(required_files)
    if ~isfile(required_files{k})
        error('Bouchet_Kantola_Field2D:MissingDeploymentFile', ...
            'Missing required deployment file: %s',required_files{k});
    end
end

run_id = read_positive_integer('SLURM_ARRAY_TASK_ID');
if run_id > 9
    error('Bouchet_Kantola_Field2D:BadRunId', ...
        'SLURM_ARRAY_TASK_ID must select run 1 through 9.');
end
allocated_cpus = read_positive_integer('SLURM_CPUS_PER_TASK');
%% The daily ERW columns run on a pool; the hourly loop is serial over
%% columns whatever the pool, so the worker count sets only the ERW share
%% of the run time. It is part of the stamp and of the provenance.
workers = Resolve_2D_Parfor_Workers(n_columns);
if workers ~= min(12,n_columns)
    error('Bouchet_Kantola_Field2D:BadWorkerCount', ...
        'TC_PARALLEL_WORKERS must be 12 for this campaign.');
end
if allocated_cpus < workers
    error('Bouchet_Kantola_Field2D:InsufficientCpus', ...
        'Reserve at least %d CPUs.',workers);
end
if ~strcmp(getenv('TC_CKPT_MINUTES'),'120')
    error('Bouchet_Kantola_Field2D:BadCheckpointInterval', ...
        'TC_CKPT_MINUTES must be 120.');
end
if ~isempty(getenv('TC_CKPT_FORCE_RESUME'))
    error('Bouchet_Kantola_Field2D:ForcedResumeForbidden', ...
        'TC_CKPT_FORCE_RESUME must remain unset for production.');
end
if ~strcmp(strtrim(getenv('TC_CAMPAIGN_TAG')),campaign_tag)
    error('Bouchet_Kantola_Field2D:BadCampaignTag', ...
        'TC_CAMPAIGN_TAG must be %s.',campaign_tag);
end
if ~license('test','Distrib_Computing_Toolbox')
    error('Bouchet_Kantola_Field2D:ParallelToolboxUnavailable', ...
        'Parallel Computing Toolbox is required for this campaign.');
end
validate_only = strtrim(getenv('TC_VALIDATE_ONLY'));
if ~isempty(validate_only) && ~strcmp(validate_only,'1')
    error('Bouchet_Kantola_Field2D:BadValidationMode', ...
        'TC_VALIDATE_ONLY must be unset or 1.');
end
%% The 1-D campaign directory supplies the parameter set, the provenance
%% assertion and the digest helper, and the grid file is this campaign's
%% own input; all of them are part of the gate.
stamp_dependencies = { ...
    fullfile(one_d_dir,'Campaign_Parameters.m'), ...
    fullfile(one_d_dir,'Assert_Base_Provenance.m'), ...
    fullfile(one_d_dir,'Fingerprint_Digest.m'), ...
    topology_file};
if strcmp(validate_only,'1')
    run_validation(code_dir,repo_root,campaign_dir,campaign_tag, ...
        one_d_dir,stamp_dependencies,workers,topology_file,n_columns);
    return
end

results_dir = fullfile(repo_root,'v2','SIM','Subfield_12',campaign_tag);
logs_dir = fullfile(repo_root,'v2','logs','Subfield_12',campaign_tag);
stamp_file = fullfile(logs_dir,'validation_stamp.txt');
if ~isfile(stamp_file)
    error('Bouchet_Kantola_Field2D:MissingValidation', ...
        'Run the R2025a validation job before production.');
end
expected_stamp = Bouchet_Campaign_Validation_Fingerprint( ...
    fullfile(repo_root,'v2'),campaign_dir,campaign_tag,n_columns,workers, ...
    stamp_dependencies{:});
if ~strcmp(strtrim(fileread(stamp_file)),expected_stamp)
    error('Bouchet_Kantola_Field2D:StaleValidation', ...
        'Code, inputs, grid, tests or numerical runtime changed after validation.');
end
stamp_digest = Fingerprint_Digest(expected_stamp);
setenv('TC_VALIDATION_FINGERPRINT',stamp_digest);
if ~exist(results_dir,'dir'); mkdir(results_dir); end
if ~exist(logs_dir,'dir'); mkdir(logs_dir); end

lock_cleanup = Bouchet_Acquire_Run_Lock( ...
    fullfile(logs_dir,'locks',sprintf('run_%02d.lock',run_id))); %#ok<NASGU>

if strcmp(Bouchet_Check_Existing_Result( ...
        results_dir,run_id,stamp_digest,n_columns,P.NDays),'complete')
    fprintf('Run %d needs no work under the current validation.\n',run_id);
    return
end

worker_log = fullfile(logs_dir,sprintf('worker_run%02d.log',run_id));
diary(worker_log);
diary_cleanup = onCleanup(@() diary('off')); %#ok<NASGU>
fprintf(['Campaign: %s\nGrid: %s (%d cells of %g m, topography %d)\n' ...
    'Slurm run: %d\nWorkers: %d\nAppelo alpha: %.17g\nKNmaxd: %.17g d^-1\n' ...
    'Kbot: %.17g mm h^-1\nCEC: %g cmolc kg^-1\n'], ...
    campaign_tag,topology_file,n_columns,topo.cell_size_m, ...
    isfield(topo,'topography_attached') && topo.topography_attached, ...
    run_id,workers,P.AppeloAlpha,P.KNmaxd,P.Kbot,P.CEC_cmolc_kg);

status_file = fullfile(logs_dir,sprintf('run_%02d.status',run_id));
Run_Subfield_Worker(run_id,status_file, ...
    'NDays',P.NDays,'StartDate',P.StartDate, ...
    'FieldTopologyFile',topology_file, ...
    'AppeloAlpha',P.AppeloAlpha,'AppeloLogKH',P.AppeloLogKH, ...
    'KNmaxd',P.KNmaxd,'Kbot',P.Kbot, ...
    'ParallelWorkers',workers,'SaveResults',true, ...
    'ResultDirectory',results_dir);

status = fileread(status_file);
if ~contains(status,'state=complete')
    error('Bouchet_Kantola_Field2D:IncompleteStatus', ...
        'Run %d did not write a complete status.',run_id);
end
%% Jobs execute from the shared source tree, which stays writable while
%% they run. Recompute the gate and reject the result if anything moved.
if ~strcmp(Bouchet_Campaign_Validation_Fingerprint( ...
        fullfile(repo_root,'v2'),campaign_dir,campaign_tag,n_columns,workers, ...
        stamp_dependencies{:}),expected_stamp)
    error('Bouchet_Kantola_Field2D:SourceChangedDuringRun', ...
        ['Code, inputs, grid, tests or numerical runtime changed while run %d ' ...
         'was executing. Discard its result.'],run_id);
end
if ~strcmp(Bouchet_Check_Existing_Result( ...
        results_dir,run_id,stamp_digest,n_columns,P.NDays),'complete')
    error('Bouchet_Kantola_Field2D:MissingResult', ...
        'Run %d reported completion but wrote no matching final MAT.',run_id);
end

existing = dir(fullfile(results_dir,sprintf('TC_2D_run%02d_*.mat',run_id)));
existing = existing(~[existing.isdir]);
result_file = fullfile(existing(1).folder,existing(1).name);
saved = load(result_file,'run_provenance');
Assert_Base_Provenance(saved.run_provenance,run_id,n_columns,workers, ...
    stamp_digest,campaign_tag,results_dir,P.NDays);
assert_field_grid_provenance(saved.run_provenance,topology_file,topo);
fprintf('Provenance: run %d carries the calibrated CEC 13 identity on %s.\n', ...
    run_id,topology_name);

fprintf('Lateral conservation over the whole run:\n');
Check_Lateral_Conservation(result_file,'Bouchet_Kantola_Field2D');
fprintf('Run %d complete under %s.\n',run_id,campaign_tag);
end

function run_validation(code_dir,repo_root,campaign_dir,campaign_tag, ...
    one_d_dir,stamp_dependencies,workers,topology_file,n_columns)
test_dir = fullfile(code_dir,'Tests');
addpath(test_dir,one_d_dir);
old_dir = pwd;
cleanup = onCleanup(@() cd(old_dir));
cd(test_dir);
%% The complete audit (the post-legacy suite, which carries the topography
%% and field-grid tests), the Linux deployment, and the grid itself for one
%% day across the first application, before any multi-day element starts.
run_post_legacy_equation_audit_tests(1);
cd(test_dir);
test_biogeochemistry_io_path_portability;
test_bouchet_linux_portability({campaign_dir,one_d_dir});
cd(old_dir);
smoke_dir = tempname;
mkdir(smoke_dir);
smoke_cleanup = onCleanup(@() rmdir(smoke_dir,'s'));
ckpt_minutes = getenv('TC_CKPT_MINUTES');
setenv('TC_CKPT_MINUTES','0');
Run_One_Subfield(5,'FieldTopologyFile',topology_file,'NDays',1, ...
    'StartDate','2016-11-10','ParallelWorkers',workers, ...
    'SaveResults',false,'ResultDirectory',smoke_dir);
setenv('TC_CKPT_MINUTES',ckpt_minutes);
clear smoke_cleanup
logs_dir = fullfile(repo_root,'v2','logs','Subfield_12',campaign_tag);
fingerprint = Bouchet_Campaign_Validation_Fingerprint( ...
    fullfile(repo_root,'v2'),campaign_dir,campaign_tag,n_columns,workers, ...
    stamp_dependencies{:});
Bouchet_Write_Validation_Stamp( ...
    fullfile(logs_dir,'validation_stamp.txt'),fingerprint);
fprintf('PASS: Bouchet R2025a CEC 13 field-grid 2-D validation (%s)\n', ...
    campaign_tag);
end

function [topology_file,topology_name] = read_topology_file(repo_root)
%% An absolute path to a grid saved by Build_Field_Grid_Topology, outside
%% the deployed root: a validated root never gains files.
topology_file = strtrim(getenv('TC_FIELD_TOPOLOGY'));
if isempty(topology_file) || ~isfile(topology_file)
    error('Bouchet_Kantola_Field2D:MissingTopology', ...
        'TC_FIELD_TOPOLOGY must name an existing grid file, not ''%s''.', ...
        topology_file);
end
[~,attributes] = fileattrib(topology_file);
topology_file = attributes.Name;
[~,root_attributes] = fileattrib(repo_root);
if startsWith(topology_file,[root_attributes.Name filesep])
    error('Bouchet_Kantola_Field2D:TopologyInsideRoot', ...
        'Keep the grid file beside the deployed root, not inside it.');
end
[~,topology_name,extension] = fileparts(topology_file);
if ~strcmp(extension,'.mat') || isempty(regexp(topology_name, ...
        '^[A-Za-z0-9_]{1,32}$','once'))
    error('Bouchet_Kantola_Field2D:BadTopologyName', ...
        ['The grid file must be <name>.mat with a name of at most 32 ' ...
         'letters, digits and underscores; it names the campaign.']);
end
end

function topo = load_topology(topology_file)
S = load(topology_file,'topo');
if ~isfield(S,'topo') || ~isstruct(S.topo) || ~isscalar(S.topo) || ...
        ~isfield(S.topo,'field_grid_build') || ~isfield(S.topo,'N_cells') || ...
        S.topo.N_cells < 2
    error('Bouchet_Kantola_Field2D:BadTopology', ...
        '%s is not a field grid saved by Build_Field_Grid_Topology.',topology_file);
end
topo = S.topo;
end

function assert_field_grid_provenance(p,topology_file,topo)
%% The result names the grid file it ran on and that file's content.
fid = fopen(topology_file,'rb');
bytes = fread(fid,Inf,'*uint8');
fclose(fid);
md = javaMethod('getInstance','java.security.MessageDigest','SHA-256');
md.update(typecast(bytes,'int8'));
hex = lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
attached = isfield(topo,'topography_attached') && topo.topography_attached;
ok = isfield(p,'field_topology_file') && strcmp(p.field_topology_file,topology_file) && ...
    isfield(p,'field_topology_file_sha256') && strcmp(p.field_topology_file_sha256,hex) && ...
    isfield(p,'opt_dem') && p.opt_dem == attached && ...
    isfield(p,'n_columns') && p.n_columns == topo.N_cells;
if ~ok
    error('Bouchet_Kantola_Field2D:WrongProvenance', ...
        'The result does not carry the grid it was validated on: %s.',topology_file);
end
end

function value = read_positive_integer(name)
raw = strtrim(getenv(name));
value = str2double(raw);
if isempty(raw) || ~(isfinite(value) && value == fix(value) && value >= 1)
    error('Bouchet_Kantola_Field2D:BadEnvironment', ...
        '%s must be a positive integer.',name);
end
end
