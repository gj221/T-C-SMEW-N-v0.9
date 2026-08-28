function Validate_Bouchet_1D_Campaign
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Validate the corrected 1-D source in the Bouchet R2025a environment.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

campaign_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(campaign_dir));
code_dir = fullfile(repo_root,'v2','T&C_Code');
test_dir = fullfile(code_dir,'Tests');
addpath(test_dir,code_dir,fullfile(repo_root,'v2','run'),campaign_dir);

matlab_release = version('-release');
if ~strcmp(matlab_release,'2025a')
    error('Bouchet_Kantola_1D:WrongMatlabRelease', ...
        'This validation requires MATLAB R2025a; loaded R%s.',matlab_release);
end
Bouchet_Single_Thread_Runtime('Bouchet_Kantola_1D');
allocated_cpus = str2double(getenv('SLURM_CPUS_PER_TASK'));
if ~(isfinite(allocated_cpus) && allocated_cpus >= 6)
    error('Bouchet_Kantola_1D:ValidationCpuCount', ...
        'The complete audit requires at least six allocated CPUs.');
end
if ~isempty(strtrim(getenv('TC_CKPT_FORCE_RESUME')))
    error('Bouchet_Kantola_1D:ForcedResumeForbidden', ...
        'TC_CKPT_FORCE_RESUME must remain unset during validation.');
end

old_dir = pwd;
cleanup = onCleanup(@() cd(old_dir));
cd(test_dir);
run_post_legacy_equation_audit_tests(1);
test_bouchet_linux_portability;
Validate_Bouchet_1D_Preapplication_Identity;
campaign_tag = 'Kantola_5yr_1D_Runs10to12_ScheduleFix_R2025a_Bouchet';
logs_dir = fullfile(repo_root,'v2','logs','Subfield_12',campaign_tag);
fingerprint = Bouchet_Campaign_Validation_Fingerprint( ...
    fullfile(repo_root,'v2'),campaign_dir,campaign_tag,1,0);
Bouchet_Write_Validation_Stamp( ...
    fullfile(logs_dir,'validation_stamp.txt'),fingerprint);
fprintf('PASS: Bouchet R2025a corrected 1-D validation\n');
end
