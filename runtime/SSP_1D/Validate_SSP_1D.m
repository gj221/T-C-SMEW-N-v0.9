function Validate_SSP_1D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Validate one Kantola SSP campaign in the Bouchet R2025a environment and
%% write its 1-D validation stamp.
%%
%% The audit this runs is the calibrated campaign's own, unchanged.  That is
%% deliberate: Validate_Preapplication_Identity and
%% Validate_All12_Smoke exercise the DEFAULT forcing path, with
%% no Scenario set, so passing them under the scenario-aware Run_One_Subfield
%% demonstrates that the scenario axis leaves the historical path untouched.
%% A regression there means the axis is not inert and the campaign must not
%% run, which is exactly what a gate is for.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[scenario,horizon] = read_scenario_and_horizon;
campaign_tag = sprintf( ...
    'SSP_1D_%s_%s%s', ...
    horizon,scenario,case_suffix(read_calcite_case));
two_d_tag = ...
    'SSP_2D';
one_d_tag = ...
    'Base_1D';
campaign_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(campaign_dir));
two_d_dir = fullfile(repo_root,'runtime',two_d_tag);
one_d_dir = fullfile(repo_root,'runtime',one_d_tag);
code_dir = fullfile(repo_root,'v2','T&C_Code');
test_dir = fullfile(code_dir,'Tests');
addpath(test_dir,code_dir,fullfile(code_dir,'Inputs'), ...
    fullfile(code_dir,'Parameters'),fullfile(repo_root,'v2','run'), ...
    campaign_dir,one_d_dir);

matlab_release = version('-release');
if ~strcmp(matlab_release,'2025a')
    error('Bouchet_Kantola_SSP:WrongMatlabRelease', ...
        'This validation requires MATLAB R2025a; loaded R%s.',matlab_release);
end
Bouchet_Single_Thread_Runtime('Bouchet_Kantola_SSP');
allocated_cpus = str2double(getenv('SLURM_CPUS_PER_TASK'));
if ~(isfinite(allocated_cpus) && allocated_cpus >= 6)
    error('Bouchet_Kantola_SSP:ValidationCpuCount', ...
        'The complete audit requires at least six allocated CPUs.');
end
if ~isempty(strtrim(getenv('TC_CKPT_FORCE_RESUME')))
    error('Bouchet_Kantola_SSP:ForcedResumeForbidden', ...
        'TC_CKPT_FORCE_RESUME must remain unset during validation.');
end
if ~strcmp(strtrim(getenv('TC_CAMPAIGN_TAG')),campaign_tag)
    error('Bouchet_Kantola_SSP:BadCampaignTag', ...
        'TC_CAMPAIGN_TAG must be %s.',campaign_tag);
end
if ~isfolder(two_d_dir)
    error('Bouchet_Kantola_SSP:MissingTwoDCampaign', ...
        'The companion 2-D campaign directory is missing: %s',two_d_dir);
end
if ~isfolder(one_d_dir)
    error('Bouchet_Kantola_SSP:MissingCalibratedCampaign', ...
        ['The five-year campaign directory supplies the calibrated set ' ...
        'and is missing: %s'],one_d_dir);
end
%% The scenario forcing is what distinguishes this campaign from the
%% calibrated one, so its absence is a validation failure, not a run-time
%% surprise three hours into production.
forcing_file = fullfile(code_dir,'Inputs', ...
    sprintf('Data_Kantola_%s.mat',scenario));
if ~isfile(forcing_file)
    error('Bouchet_Kantola_SSP:MissingScenarioForcing', ...
        'The scenario forcing is missing: %s',forcing_file);
end
equivalence_baseline = fullfile(repo_root,'v2','SIM','EQUIV_1D_TMP', ...
    'TC_ERW_results_Kantola_crop_control.mat');
if ~isfile(equivalence_baseline)
    error('Bouchet_Kantola_SSP:MissingEquivalenceBaseline', ...
        ['The equivalence baseline is missing and must not be rebuilt ' ...
         'here: %s'],equivalence_baseline);
end

%% Identical to the list both SSP tasks build, or production fails against
%% the stamp this writes.
stamp_dependencies = { ...
    fullfile(one_d_dir,'Campaign_Parameters.m'), ...
    fullfile(one_d_dir,'Assert_Base_Provenance.m'), ...
    fullfile(one_d_dir,'Fingerprint_Digest.m'), ...
    fullfile(one_d_dir,'test_all_12_configuration.m'), ...
    fullfile(campaign_dir,'submit_ssp.sh')};

old_dir = pwd;
cleanup = onCleanup(@() cd(old_dir));
cd(test_dir);
stage('post-legacy equation audit');
run_post_legacy_equation_audit_tests(1);
stage('Linux portability of the deployed campaign directories');
test_bouchet_linux_portability({campaign_dir,two_d_dir,one_d_dir});
stage('all-12 configuration');
test_all_12_configuration;
stage('pre-application identity, four 315-day runs on a 4-process pool');
Validate_Preapplication_Identity;
stage('three-day smoke on one and twelve columns');
Validate_All12_Smoke;
stage('validation stamp');
logs_dir = fullfile(repo_root,'v2','logs','Subfield_12',campaign_tag);
fingerprint = Bouchet_Campaign_Validation_Fingerprint( ...
    fullfile(repo_root,'v2'),campaign_dir,campaign_tag,1,0, ...
    stamp_dependencies{:});
Bouchet_Write_Validation_Stamp( ...
    fullfile(logs_dir,'validation_stamp.txt'),fingerprint);
fprintf('PASS: Bouchet R2025a CEC 13 SSP 1-D validation (%s)\n',campaign_tag);
end

function calcite_case = read_calcite_case
%% Optional campaign axis: an ERW workbook row for the calcite runs only,
%% such as 8, the agronomic 5 t/ha rate. Unset is the standard campaign.
raw = strtrim(getenv('TC_CALCITE_CASE'));
calcite_case = [];
if ~isempty(raw)
    calcite_case = str2double(raw);
    if ~(isfinite(calcite_case) && calcite_case == fix(calcite_case) && calcite_case >= 1)
        error('Bouchet_Kantola_SSP:BadCalciteCase', ...
            'TC_CALCITE_CASE must be unset or a positive integer.');
    end
end
end

function suffix = case_suffix(calcite_case)
suffix = '';
if ~isempty(calcite_case)
    suffix = sprintf('_case%d',calcite_case);
end
end

function [scenario,horizon] = read_scenario_and_horizon
scenario = strtrim(getenv('TC_SSP'));
if ~ismember(scenario,{'ssp126','ssp245','ssp585'})
    error('Bouchet_Kantola_SSP:BadScenario', ...
        'TC_SSP must be ssp126, ssp245 or ssp585, not ''%s''.',scenario);
end
horizon = strtrim(getenv('TC_HORIZON'));
if ~ismember(horizon,{'5yr','35yr'})
    error('Bouchet_Kantola_SSP:BadHorizon', ...
        'TC_HORIZON must be 5yr or 35yr, not ''%s''.',horizon);
end
end

function stage(name)
%% Timestamped banner so the Slurm output shows which stage a failed
%% validation stopped in.
fprintf('\n[%s] Validation stage: %s\n', ...
    char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')),name);
end
