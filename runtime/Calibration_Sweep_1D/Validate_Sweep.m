function Validate_Sweep
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Validate the joint 1-D sweep under its Bouchet runtime.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

campaign_tag = ...
    'Calibration_Sweep_1D';
campaign_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(campaign_dir));
v2_dir = fullfile(repo_root,'v2');
code_dir = fullfile(v2_dir,'T&C_Code');
test_dir = fullfile(code_dir,'Tests');
run_dir = fullfile(v2_dir,'run');
addpath(test_dir,code_dir,run_dir,fullfile(code_dir,'Parameters'), ...
    campaign_dir);

matlab_release = version('-release');
if ~strcmp(matlab_release,'2025a')
    error('Bouchet_Joint_Sweep:WrongMatlabRelease', ...
        'This validation requires MATLAB R2025a; loaded R%s.',matlab_release);
end
Bouchet_Single_Thread_Runtime('Bouchet_Joint_Sweep');
allocated_cpus = str2double(getenv('SLURM_CPUS_PER_TASK'));
if ~(isfinite(allocated_cpus) && allocated_cpus >= 6)
    error('Bouchet_Joint_Sweep:ValidationCpuCount', ...
        'The complete audit requires at least six allocated CPUs.');
end
if ~isempty(strtrim(getenv('TC_CKPT_FORCE_RESUME')))
    error('Bouchet_Joint_Sweep:ForcedResumeForbidden', ...
        'TC_CKPT_FORCE_RESUME must remain unset during validation.');
end
if ~strcmp(strtrim(getenv('TC_CAMPAIGN_TAG')),campaign_tag)
    error('Bouchet_Joint_Sweep:CampaignIdentity', ...
        'TC_CAMPAIGN_TAG does not identify this sweep.');
end
validate_campaign_window(fullfile(code_dir,'Inputs','Data_Kantola.mat'));
Load_Sweep_Observations(v2_dir);
fingerprint_before = Sweep_Fingerprint( ...
    v2_dir,campaign_dir,campaign_tag);
if contains(fingerprint_before,':missing')
    error('Bouchet_Joint_Sweep:MissingFingerprintDependency', ...
        'A campaign fingerprint dependency is missing.');
end
setenv('TC_VALIDATION_FINGERPRINT', ...
    Fingerprint_Digest(fingerprint_before));

old_dir = pwd;
cleanup = onCleanup(@() cd(old_dir));
cd(test_dir);
run_post_legacy_equation_audit_tests(1);
test_bouchet_linux_portability(campaign_dir);
validate_joint_sweep_smoke;

logs_dir = fullfile(v2_dir,'logs','Subfield_12',campaign_tag);
fingerprint_after = Sweep_Fingerprint( ...
    v2_dir,campaign_dir,campaign_tag);
if ~strcmp(fingerprint_before,fingerprint_after)
    error('Bouchet_Joint_Sweep:SourceChangedDuringValidation', ...
        'Code, inputs, tests or numerical runtime changed during validation.');
end
Bouchet_Write_Validation_Stamp( ...
    fullfile(logs_dir,'validation_stamp.txt'),fingerprint_before);
% Every stamp is also kept under its own digest so that a later rescore can
% prove which files changed between the stamp a result carries and the
% stamp now on disk (Accepted_Stamp_Digests).
archive_file = fullfile(logs_dir,'validation_stamps', ...
    [Fingerprint_Digest(fingerprint_before) '.txt']);
if ~isfile(archive_file)
    Bouchet_Write_Validation_Stamp(archive_file,fingerprint_before);
end
fprintf('PASS: Bouchet R2025a joint 1-D parameter sweep validation\n');

end

function validate_campaign_window(forcing_file)
forcing = load(forcing_file,'Date');
if ~isfield(forcing,'Date') || ~isnumeric(forcing.Date) || ...
        ~isreal(forcing.Date) || any(~isfinite(forcing.Date(:)))
    error('Run_One_Subfield:BadWindow', ...
        'Data_Kantola.mat does not contain a finite real Date vector.');
end
dates = double(forcing.Date(:));
hourly_names = {'Pr','Ta','Ws','ea','SAD1','SAD2','SAB1','SAB2', ...
    'Pre','N','Tdew','esat','PARB','PARD'};
contents = whos('-file',forcing_file);
content_names = {contents.name};
for k = 1:numel(hourly_names)
    match = find(strcmp(content_names,hourly_names{k}));
    if numel(match) ~= 1 || ...
            ~ismember(contents(match).class,{'double','single'}) || ...
            prod(double(contents(match).size)) ~= numel(dates)
        error('Run_One_Subfield:BadWindow', ...
            'Hourly forcing %s must match Date.',hourly_names{k});
    end
end
start_date = datenum(2016,1,1); %#ok<DATNM>
[date_error,x1] = min(abs(dates-start_date));
n_hours = 1826*24;
x2 = x1+n_hours-1;
if isempty(x1) || date_error > 1/48 || x2 > numel(dates)
    error('Run_One_Subfield:BadWindow', ...
        'The 1,826-day production window is unavailable in forcing.');
end
date_tol = 64*eps(max(abs(dates(x1:x2))));
if any(abs(diff(dates(x1:x2))-1/24) > date_tol)
    error('Run_One_Subfield:BadWindow', ...
        'The production forcing window is not hourly and contiguous.');
end
expected_last_hour = datenum(2020,12,30,23,0,0); %#ok<DATNM>
if abs(dates(x2)-expected_last_hour) > 1/48
    error('Run_One_Subfield:BadWindow', ...
        'The production window does not end at 23:00 on 30 December 2020.');
end
end

function validate_joint_sweep_smoke
cases = Sweep_Cases;
smoke = [ ...
    select_case(cases,10,1.3,0.0048,NaN,1) ...
    select_case(cases,11,2.35,0.0032,0.03,1) ...
    select_case(cases,12,6,0.0064,0.01,1)];
root = tempname(tempdir);
mkdir(root);
old_ckpt = getenv('TC_CKPT_MINUTES');
old_force = getenv('TC_CKPT_FORCE_RESUME');
old_workers = getenv('TC_PARALLEL_WORKERS');
old_fingerprint = getenv('TC_VALIDATION_FINGERPRINT');
cleanup = onCleanup(@() restore_smoke_environment( ...
    root,old_ckpt,old_force,old_workers,old_fingerprint));
setenv('TC_CKPT_MINUTES','0');
setenv('TC_CKPT_FORCE_RESUME','');
setenv('TC_PARALLEL_WORKERS','0');

n_days = 3;
start_date = datenum(2016,11,9); %#ok<DATNM>
for k = 1:numel(smoke)
    case_cfg = smoke(k);
    result_dir = fullfile(root,case_cfg.tag);
    result_file = Run_One_Subfield(case_cfg.run_id,'NDays',n_days, ...
        'StartDate',start_date,'AppeloAlpha',case_cfg.AppeloAlpha, ...
        'AppeloLogKH',case_cfg.AppeloLogKH, ...
        'KNmaxd',case_cfg.KNmaxd,'Kbot',case_cfg.Kbot, ...
        'DissolutionFactor',case_cfg.DissolutionFactor, ...
        'ParallelWorkers',0,'SaveResults',true, ...
        'ResultDirectory',result_dir);
    S = load(result_file,'run_provenance','PAR','O');
    p = S.run_provenance;
    assert(isequal(p.KNmaxd,case_cfg.KNmaxd) && ...
           isequaln(p.Kbot,case_cfg.Kbot) && ...
           isequal(p.AppeloAlpha,case_cfg.AppeloAlpha) && ...
           isequal(reshape(p.AppeloLogKH,1,[]),case_cfg.AppeloLogKH) && ...
           isequal(p.dissolution_factor,case_cfg.DissolutionFactor) && ...
           p.n_columns == 1 && p.n_days == n_days && ...
           p.n_timesteps == 24*n_days && ...
           strcmp(p.numerical_identity,getenv('TC_NUMERICAL_IDENTITY')), ...
        'Bouchet_Joint_Sweep:SmokeProvenance', ...
        'The smoke result lost a requested sweep input.');
    assert(iscell(S.PAR) && isscalar(S.PAR) && ...
           isequaln(double(S.PAR{1}.Kbot),case_cfg.Kbot), ...
        'Bouchet_Joint_Sweep:SmokeKbot', ...
        'The bottom-boundary override did not reach the hydrology parameters.');
    [~,~,~,~,~,effective_alpha,effective_logKH] = ...
        Resolve_Appelo_Constants(S.PAR{1}.bulk_pH_params.appelo);
    assert(strcmpi(S.PAR{1}.bulk_pH_params.scheme,'appelo') && ...
           isequal(effective_alpha,case_cfg.AppeloAlpha) && ...
           isequal(effective_logKH,case_cfg.AppeloLogKH), ...
        'Bouchet_Joint_Sweep:SmokeAppelo', ...
        'The smoke run is not using the requested Appelo configuration.');
    assert(isnumeric(S.O) && ~isempty(S.O) && all(isfinite(S.O(:))), ...
        'Bouchet_Joint_Sweep:SmokeState', ...
        'The representative smoke run produced invalid soil water state.');
    expected_tags = {sprintf('__alpha%s',number_tag(case_cfg.AppeloAlpha)), ...
        sprintf('__knmaxd%s',number_tag(case_cfg.KNmaxd)), ...
        sprintf('__kbot%s',kbot_tag(case_cfg.Kbot)), ...
        sprintf('__fd%s',number_tag(case_cfg.DissolutionFactor))};
    assert(all(cellfun(@(x) contains(result_file,x),expected_tags)), ...
        'Bouchet_Joint_Sweep:SmokeFilename', ...
        'The smoke filename does not identify every sweep input.');
    remaining = dir(fullfile(result_dir,'checkpoints','ckpt_run*.mat'));
    assert(isempty(remaining), ...
        'Bouchet_Joint_Sweep:SmokeCheckpoint', ...
        'A completed smoke run retained a rolling checkpoint.');
end
clear cleanup
end

function case_cfg = select_case(cases,run_id,alpha,knmaxd,kbot,release_factor)
same_kbot = arrayfun(@(x) isequaln(x.Kbot,kbot),cases);
idx = find([cases.run_id] == run_id & ...
    abs([cases.AppeloAlpha]-alpha) <= 8*eps(max(1,abs(alpha))) & ...
    [cases.KNmaxd] == knmaxd & same_kbot & ...
    [cases.DissolutionFactor] == release_factor);
if numel(idx) ~= 1
    error('Bouchet_Joint_Sweep:SmokeCase', ...
        'The representative smoke case is missing or ambiguous.');
end
case_cfg = cases(idx);
end

function restore_smoke_environment(root,ckpt,force_resume,workers,fingerprint)
setenv('TC_CKPT_MINUTES',ckpt);
setenv('TC_CKPT_FORCE_RESUME',force_resume);
setenv('TC_PARALLEL_WORKERS',workers);
setenv('TC_VALIDATION_FINGERPRINT',fingerprint);
pool = gcp('nocreate');
if ~isempty(pool); delete(pool); end
if isfolder(root) && startsWith(root,tempdir)
    rmdir(root,'s');
end
end

function tag = number_tag(value)
tag = strrep(sprintf('%.12g',value),'.','p');
tag = strrep(tag,'-','m');
tag = strrep(tag,'+','');
end

function tag = kbot_tag(value)
if isnan(value)
    tag = 'free';
else
    tag = number_tag(value);
end
end
