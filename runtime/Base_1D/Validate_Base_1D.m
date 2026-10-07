function Validate_Base_1D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Validate the calibrated CEC 13 all-12 campaign in the Bouchet R2025a
%% environment and write the 1-D validation stamp.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

two_d_tag = ...
    'Base_2D';
campaign_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(campaign_dir));
two_d_dir = fullfile(repo_root,'runtime',two_d_tag);
code_dir = fullfile(repo_root,'v2','T&C_Code');
test_dir = fullfile(code_dir,'Tests');
addpath(test_dir,code_dir,fullfile(repo_root,'v2','run'),campaign_dir);
[~,case_suffix] = Calcite_Case;
campaign_tag = ['Base_1D' case_suffix];

matlab_release = version('-release');
if ~strcmp(matlab_release,'2025a')
    error('Bouchet_Kantola_CEC13:WrongMatlabRelease', ...
        'This validation requires MATLAB R2025a; loaded R%s.',matlab_release);
end
Bouchet_Single_Thread_Runtime('Bouchet_Kantola_CEC13');
allocated_cpus = str2double(getenv('SLURM_CPUS_PER_TASK'));
if ~(isfinite(allocated_cpus) && allocated_cpus >= 6)
    error('Bouchet_Kantola_CEC13:ValidationCpuCount', ...
        'The complete audit requires at least six allocated CPUs.');
end
if ~isempty(strtrim(getenv('TC_CKPT_FORCE_RESUME')))
    error('Bouchet_Kantola_CEC13:ForcedResumeForbidden', ...
        'TC_CKPT_FORCE_RESUME must remain unset during validation.');
end
if ~strcmp(strtrim(getenv('TC_CAMPAIGN_TAG')),campaign_tag)
    error('Bouchet_Kantola_CEC13:BadCampaignTag', ...
        'TC_CAMPAIGN_TAG must be %s.',campaign_tag);
end
if ~isfolder(two_d_dir)
    error('Bouchet_Kantola_CEC13:MissingTwoDCampaign', ...
        'The companion 2-D campaign directory is missing: %s',two_d_dir);
end
equivalence_baseline = fullfile(repo_root,'v2','SIM','EQUIV_1D_TMP', ...
    'TC_ERW_results_Kantola_crop_control.mat');
if ~isfile(equivalence_baseline)
    error('Bouchet_Kantola_CEC13:MissingEquivalenceBaseline', ...
        ['The equivalence baseline is missing and must not be rebuilt ' ...
         'here: %s'],equivalence_baseline);
end

old_dir = pwd;
cleanup = onCleanup(@() cd(old_dir));
cd(test_dir);
stage('post-legacy equation audit');
run_post_legacy_equation_audit_tests(1);
stage('Linux portability of the deployed campaign directories');
test_bouchet_linux_portability({campaign_dir,two_d_dir});
stage('all-12 configuration');
test_all_12_configuration;
stage('pre-application identity, four 315-day runs on a 4-process pool');
Validate_Preapplication_Identity;
stage('three-day smoke on one and twelve columns');
Validate_All12_Smoke;
stage('validation stamp');
logs_dir = fullfile(repo_root,'v2','logs','Subfield_12',campaign_tag);
fingerprint = Bouchet_Campaign_Validation_Fingerprint( ...
    fullfile(repo_root,'v2'),campaign_dir,campaign_tag,1,0);
Bouchet_Write_Validation_Stamp( ...
    fullfile(logs_dir,'validation_stamp.txt'),fingerprint);
fprintf('PASS: Bouchet R2025a calibrated CEC 13 1-D validation\n');
end

function stage(name)
%% Timestamped banner so the Slurm output shows which stage a failed
%% validation stopped in.
fprintf('\n[%s] Validation stage: %s\n', ...
    char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')),name);
end
