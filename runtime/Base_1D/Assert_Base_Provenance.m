function Assert_Base_Provenance( ...
    p,run_id,n_columns,parallel_workers,stamp_digest,campaign_tag, ...
    results_dir,n_days)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Reject a result that does not carry the calibrated CEC 13 identity.
%%
%% Every saved MAT records run_provenance; this compares it with
%% Campaign_Parameters and with the campaign gate, so a result
%% produced with a default alpha, a different bottom boundary, another CEC
%% or another stamp cannot stand in the campaign results directory.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P = Campaign_Parameters;
if nargin < 8 || isempty(n_days); n_days = P.NDays; end
if ~(isstruct(p) && isscalar(p))
    wrong_provenance('run_provenance must be a scalar structure.');
end
required = {'KNmaxd','Kbot','exchange_scheme','AppeloAlpha', ...
    'AppeloLogKH','CEC_cmolc_kg','CECpHSensitivity', ...
    'validation_fingerprint','matlab_release','numerical_identity', ...
    'max_threads','parallel_workers','n_columns','n_days', ...
    'n_timesteps','campaign','date_start','date_end', ...
    'result_directory','slurm_job_id','slurm_restart_count'};
for k = 1:numel(required)
    if ~isfield(p,required{k})
        wrong_provenance(sprintf('run_provenance.%s is missing.',required{k}));
    end
end

knmaxd = finite_scalar(p.KNmaxd,'KNmaxd');
kbot = kbot_scalar(p.Kbot);
alpha = finite_vector(p.AppeloAlpha,'AppeloAlpha');
logKH = double(p.AppeloLogKH);
if ~isreal(logKH) || mod(numel(logKH),5) ~= 0 || any(~isfinite(logKH(:)))
    wrong_provenance('AppeloLogKH must hold five finite values per column.');
end
if isvector(logKH)
    logKH = reshape(logKH,1,5);
elseif size(logKH,2) ~= 5
    logKH = logKH.';
end
cec = finite_vector(p.CEC_cmolc_kg,'CEC_cmolc_kg');
cec_slope = finite_vector(p.CECpHSensitivity,'CECpHSensitivity');
schemes = text_list(p.exchange_scheme,'exchange_scheme');
release = text_scalar(p.matlab_release,'matlab_release');
fingerprint = text_scalar(p.validation_fingerprint,'validation_fingerprint');
numerical_identity = text_scalar(p.numerical_identity,'numerical_identity');
campaign = text_scalar(p.campaign,'campaign');
result_directory = text_scalar(p.result_directory,'result_directory');
job_id = text_scalar(p.slurm_job_id,'slurm_job_id');
restart_count = str2double(text_scalar( ...
    p.slurm_restart_count,'slurm_restart_count'));
expected_identity = sprintf(['cpumodel=9575f;matlab_release=2025a;' ...
    'single_comp_thread=1;parallel_workers=%d'],parallel_workers);
date_start = finite_scalar(p.date_start,'date_start');
date_end = finite_scalar(p.date_end,'date_end');
expected_end = P.StartDate+n_days-1/24;
date_tol = 16*eps(max(abs([P.StartDate expected_end])));

checks = { ...
    'KNmaxd', isequal(knmaxd,P.KNmaxd); ...
    'Kbot', isequaln(kbot,P.Kbot); ...
    'exchange_scheme', all(strcmpi(schemes,P.exchange_scheme)); ...
    'AppeloAlpha', all(alpha == P.AppeloAlpha); ...
    'AppeloLogKH', all(arrayfun(@(r) isequal(logKH(r,:),P.AppeloLogKH), ...
        1:size(logKH,1))); ...
    'CEC_cmolc_kg', all(cec == P.CEC_cmolc_kg); ...
    'CECpHSensitivity', all(cec_slope == P.CECpHSensitivity); ...
    'matlab_release', strcmp(release,'2025a'); ...
    'numerical_identity', strcmp(numerical_identity,expected_identity); ...
    'max_threads', finite_scalar(p.max_threads,'max_threads') == 1; ...
    'parallel_workers', finite_scalar(p.parallel_workers, ...
        'parallel_workers') == parallel_workers; ...
    'n_columns', finite_scalar(p.n_columns,'n_columns') == n_columns; ...
    'n_days', finite_scalar(p.n_days,'n_days') == n_days; ...
    'n_timesteps', finite_scalar(p.n_timesteps,'n_timesteps') == 24*n_days; ...
    'validation_fingerprint', strcmp(fingerprint,char(stamp_digest)); ...
    'campaign', strcmp(campaign,char(campaign_tag)); ...
    'result_directory', strcmp(result_directory,char(results_dir)); ...
    'date_start', abs(date_start-P.StartDate) <= date_tol; ...
    'date_end', abs(date_end-expected_end) <= date_tol; ...
    'slurm_job_id', ~isempty(job_id); ...
    'slurm_restart_count', isfinite(restart_count) && ...
        restart_count == fix(restart_count) && restart_count >= 0};
failed = checks(~[checks{:,2}],1);
if ~isempty(failed)
    wrong_provenance(sprintf( ...
        'Run %d does not carry the calibrated CEC 13 identity (%s).', ...
        run_id,strjoin(failed,', ')));
end
end

function value = kbot_scalar(value)
%% Kbot is NaN for free drainage (Leakage_Bottom: Lk = Ko), the legacy
%% boundary the spin-up used, so it cannot go through finite_scalar.
if ~(isnumeric(value) && isreal(value) && isscalar(value) && ...
        (isnan(value) || (isfinite(value) && value >= 0)))
    wrong_provenance( ...
        'Kbot must be NaN for free drainage or a finite nonnegative scalar.');
end
value = double(value);
end

function value = finite_scalar(value,name)
if ~(isnumeric(value) && isreal(value) && isscalar(value) && isfinite(value))
    wrong_provenance(sprintf('%s must be a finite real scalar.',name));
end
value = double(value);
end

function value = finite_vector(value,name)
if ~(isnumeric(value) && isreal(value) && ~isempty(value) && ...
        all(isfinite(value(:))))
    wrong_provenance(sprintf('%s must hold finite real values.',name));
end
value = double(value(:));
end

function value = text_list(value,name)
if ischar(value)
    value = cellstr(value);
elseif isstring(value)
    value = cellstr(value(:));
elseif iscell(value)
    value = cellfun(@char,value(:),'UniformOutput',false);
else
    wrong_provenance(sprintf('%s must be text.',name));
end
if isempty(value)
    wrong_provenance(sprintf('%s must not be empty.',name));
end
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
error('Bouchet_Kantola_CEC13:WrongProvenance','%s',message);
end
