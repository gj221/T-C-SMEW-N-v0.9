function Assert_Sweep_Provenance( ...
    p,case_cfg,stamp_digest,campaign_tag,n_days,start_date,results_dir)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Reject a result that does not have the exact joint-sweep identity.
%%
%% stamp_digest is the current validation-stamp digest, or a cell or
%% string array of digests when the summary scores results produced under
%% an accepted prior stamp (see Accepted_Stamp_Digests).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

accepted_digests = cellstr(string(stamp_digest));
accepted_digests = accepted_digests(:)';
if isempty(accepted_digests) || ...
        any(cellfun(@(d) isempty(regexp(d,'^[0-9a-f]{64}$','once')), ...
        accepted_digests))
    wrong_provenance('The accepted validation digests must be SHA-256 hex.');
end

if ~(isstruct(p) && isscalar(p))
    wrong_provenance('run_provenance must be a scalar structure.');
end
required = {'KNmaxd','Kbot','exchange_scheme','AppeloAlpha', ...
    'AppeloLogKH','CEC_cmolc_kg','CECpHSensitivity','dissolution_factor', ...
    'validation_fingerprint','matlab_release', ...
    'numerical_identity','max_threads','parallel_workers','n_columns', ...
    'n_days','n_timesteps','campaign', ...
    'date_start','date_end','result_directory','slurm_job_id', ...
    'slurm_restart_count','slurm_array_task_id','sweep_case_id','sweep_layout'};
for k = 1:numel(required)
    if ~isfield(p,required{k})
        wrong_provenance(sprintf('run_provenance.%s is missing.',required{k}));
    end
end

knmaxd = finite_scalar(p.KNmaxd,'KNmaxd');
kbot = kbot_scalar(p.Kbot);
alpha = finite_scalar(p.AppeloAlpha,'AppeloAlpha');
logKH = double(p.AppeloLogKH);
if ~isreal(logKH) || numel(logKH) ~= 5 || any(~isfinite(logKH))
    wrong_provenance('AppeloLogKH must contain five finite real values.');
end
logKH = reshape(logKH,1,[]);
cec = finite_scalar(p.CEC_cmolc_kg,'CEC_cmolc_kg');
cec_slope = finite_scalar(p.CECpHSensitivity,'CECpHSensitivity');
release_factor = finite_scalar(p.dissolution_factor,'dissolution_factor');
scheme = text_scalar(p.exchange_scheme,'exchange_scheme');
release = text_scalar(p.matlab_release,'matlab_release');
fingerprint = text_scalar(p.validation_fingerprint, ...
    'validation_fingerprint');
numerical_identity = text_scalar(p.numerical_identity,'numerical_identity');
campaign = text_scalar(p.campaign,'campaign');
result_directory = text_scalar(p.result_directory,'result_directory');
job_id = text_scalar(p.slurm_job_id,'slurm_job_id');
restart_count = str2double(text_scalar( ...
    p.slurm_restart_count,'slurm_restart_count'));
array_task = str2double(text_scalar( ...
    p.slurm_array_task_id,'slurm_array_task_id'));
sweep_case = str2double(text_scalar(p.sweep_case_id,'sweep_case_id'));
sweep_layout = text_scalar(p.sweep_layout,'sweep_layout');
if strcmp(sweep_layout,'independent')
    expected_array_task = case_cfg.id;
elseif strcmp(sweep_layout,'packed')
    expected_array_task = case_cfg.combination_id;
else
    expected_array_task = NaN;
end
date_start = finite_scalar(p.date_start,'date_start');
date_end = finite_scalar(p.date_end,'date_end');
expected_end = start_date+n_days-1/24;
date_tol = 16*eps(max(abs([start_date expected_end])));

valid = isequal(knmaxd,case_cfg.KNmaxd) && ...
    isequaln(kbot,case_cfg.Kbot) && ...
    strcmpi(scheme,'appelo') && ...
    isequal(alpha,case_cfg.AppeloAlpha) && ...
    isequal(logKH,case_cfg.AppeloLogKH) && ...
    isequal(cec,case_cfg.CEC_cmolc_kg) && ...
    isequal(release_factor,case_cfg.DissolutionFactor) && ...
    cec_slope == 0 && ...
    strcmp(release,'2025a') && ...
    strcmp(numerical_identity, ...
        'cpumodel=9575f;matlab_release=2025a;single_comp_thread=1;parallel_workers=0') && ...
    finite_scalar(p.max_threads,'max_threads') == 1 && ...
    finite_scalar(p.parallel_workers,'parallel_workers') == 0 && ...
    finite_scalar(p.n_columns,'n_columns') == 1 && ...
    finite_scalar(p.n_days,'n_days') == n_days && ...
    finite_scalar(p.n_timesteps,'n_timesteps') == 24*n_days && ...
    any(strcmp(fingerprint,accepted_digests)) && ...
    strcmp(campaign,char(campaign_tag)) && ...
    strcmp(result_directory,char(results_dir)) && ...
    abs(date_start-start_date) <= date_tol && ...
    abs(date_end-expected_end) <= date_tol && ...
    ~isempty(job_id) && isfinite(restart_count) && ...
    restart_count == fix(restart_count) && restart_count >= 0 && ...
    isfinite(array_task) && array_task == expected_array_task && ...
    isfinite(sweep_case) && sweep_case == case_cfg.id;
if ~valid
    wrong_provenance( ...
        sprintf('The result does not belong to validated case %s.',case_cfg.tag));
end

end

function value = finite_scalar(value,name)
if ~(isnumeric(value) && isreal(value) && isscalar(value) && isfinite(value))
    wrong_provenance(sprintf('%s must be a finite real scalar.',name));
end
value = double(value);
end

function value = kbot_scalar(value)
if ~(isnumeric(value) && isreal(value) && isscalar(value) && ...
        (isnan(value) || (isfinite(value) && value >= 0)))
    wrong_provenance( ...
        'Kbot must be NaN for free drainage or a finite nonnegative scalar.');
end
value = double(value);
end

function value = text_scalar(value,name)
if iscell(value) && isscalar(value)
    value = value{1};
end
if ~(ischar(value) || (isstring(value) && isscalar(value)))
    wrong_provenance(sprintf('%s must be scalar text.',name));
end
value = char(value);
end

function wrong_provenance(message)
error('Bouchet_Joint_Sweep:WrongProvenance','%s',message);
end
