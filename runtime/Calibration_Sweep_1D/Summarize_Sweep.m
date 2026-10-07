function [T,combination_table,summary_file] = ...
    Summarize_Sweep(summary_scope)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Verify and summarize the complete 1-D alpha, KNmaxd and Kbot sweep.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 1 || isempty(summary_scope)
    summary_scope = strtrim(getenv('TC_SWEEP_SUMMARY_SCOPE'));
end
if isempty(summary_scope); summary_scope = 'full'; end
summary_scope = char(summary_scope);
if ~ismember(summary_scope,{'full','pre_control'})
    error('Bouchet_Joint_Sweep:BadSummaryScope', ...
        'The summary scope must be full or pre_control.');
end

nd = 1826;
start_date = datenum(2016,1,1); %#ok<DATNM>
campaign_tag = ...
    'Calibration_Sweep_1D';
campaign_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(campaign_dir));
v2_dir = fullfile(repo_root,'v2');
code_dir = fullfile(v2_dir,'T&C_Code');
addpath(code_dir,fullfile(code_dir,'Tests'),fullfile(v2_dir,'run'), ...
    campaign_dir);

if ~strcmp(version('-release'),'2025a')
    error('Bouchet_Joint_Sweep:WrongMatlabRelease', ...
        'Summarization requires MATLAB R2025a.');
end
Bouchet_Single_Thread_Runtime('Bouchet_Joint_Sweep');

results_root = fullfile(v2_dir,'SIM','Subfield_12',campaign_tag);
logs_root = fullfile(v2_dir,'logs','Subfield_12',campaign_tag);
if ~exist(logs_root,'dir'); mkdir(logs_root); end
if strcmp(summary_scope,'pre_control')
    [T,summary_file] = summarize_pre_control( ...
        v2_dir,campaign_dir,campaign_tag,results_root,logs_root);
    combination_table = table();
    return
end
summary_file = fullfile(logs_root,'joint_sweep_trajectory_summary.csv');
combination_file = fullfile(logs_root, ...
    'joint_sweep_combination_summary.csv');
inventory_file = fullfile(logs_root, ...
    'joint_sweep_completion_inventory.csv');
mat_file = fullfile(logs_root,'joint_sweep_summary.mat');
status_file = fullfile(logs_root,'joint_sweep_summary.status');
complete_count = 0;
try
atomic_write_text(status_file,sprintf( ...
    'state=running\ncomplete=0\nexpected=336\nselection=none\n'));

stamp_file = fullfile(logs_root,'validation_stamp.txt');
if ~isfile(stamp_file)
    error('Bouchet_Joint_Sweep:MissingValidation', ...
        'The campaign validation stamp is missing: %s',stamp_file);
end
stored_stamp = strtrim(fileread(stamp_file));
current_stamp = Sweep_Fingerprint( ...
    v2_dir,campaign_dir,campaign_tag);
if ~strcmp(stored_stamp,current_stamp)
    error('Bouchet_Joint_Sweep:StaleValidation', ...
        'Code, tests, inputs or numerical runtime changed after validation.');
end
stamp_digest = Fingerprint_Digest(stored_stamp);
[accepted_digests,lineage,rejected_digests] = Accepted_Stamp_Digests( ...
    stored_stamp,logs_root,repo_root);
for k = 1:numel(lineage)
    fprintf(['Accepting results validated under prior stamp %s; ' ...
        'scoring-side files that differ from the current stamp:\n'], ...
        lineage(k).digest);
    fprintf('  %s\n',lineage(k).changed{:});
end
for k = 1:numel(rejected_digests)
    fprintf('Not accepting prior stamp %s: %s (%s)\n', ...
        rejected_digests(k).digest,rejected_digests(k).reason, ...
        rejected_digests(k).identifier);
end
accepted_prior = strjoin(accepted_digests(2:end),':');
if isempty(accepted_prior); accepted_prior = 'none'; end

cases = Sweep_Cases;
if numel(cases) ~= 336
    error('Bouchet_Joint_Sweep:CaseCount', ...
        'The frozen case table contains %d trajectories; expected 336.', ...
        numel(cases));
end
inventory = completion_inventory(cases,results_root,logs_root);
atomic_writetable(inventory,inventory_file);
complete_count = nnz(inventory.complete);
if ~all(inventory.complete)
    atomic_write_text(status_file,sprintf( ...
        'state=incomplete\ncomplete=%d\nexpected=%d\nselection=none\n', ...
        nnz(inventory.complete),height(inventory)));
    fprintf('Completed trajectory inventory: %d/%d\n', ...
        nnz(inventory.complete),height(inventory));
    disp(inventory(~inventory.complete,:));
    error('Bouchet_Joint_Sweep:IncompleteInventory', ...
        ['The sweep is incomplete or ambiguous. The inventory was saved, ' ...
         'and no parameter combination was selected.']);
end

observations = Load_Sweep_Observations(v2_dir);
observation_fingerprint = observation_identity(observations.files);

required = { ...
    'run_provenance','R','P','TNIT_H','TNIT_L','Nreserve_H', ...
    'Nreserve_L','FertN_d','DepN_d','BfixN','N2Oflx', ...
    'N2Oflx_nit','N2flx','VOL','LEAK_NO3','LEAK_NH4','LEAK_DON', ...
    'Lk','AgrHarNut','SeedInput','Nuptake_H','Nuptake_L', ...
    'DENITflx','NITflx','Min_N','Ccrown_t','ISOIL_H','ISOIL_L', ...
    'TexN_H','TexN_L','Date','PH_bulk_d','PH_pore_d','O','PAR', ...
    'ERW','RA_H','RA_L','NPP_H','NPP_L','R_litter','R_microbe', ...
    'R_ew','topo','ColCfg'};

ncase = numel(cases);
daily = empty_daily(nd,ncase,cases);
common_dates = [];
rows = struct([]);
for k = 1:ncase
    case_cfg = cases(k);
    result_file = char(inventory.result_file(k));
    contents = string({whos('-file',result_file).name});
    missing = setdiff(string(required),contents);
    if ~isempty(missing)
        error('Bouchet_Joint_Sweep:MissingField', ...
            'Result %s is missing: %s',result_file, ...
            char(strjoin(missing,', ')));
    end
    S = load(result_file,required{:});
    result_dir = fullfile(results_root,case_cfg.tag);
    assert_result_stamp_accepted(S.run_provenance,case_cfg.tag, ...
        accepted_digests,rejected_digests);
    Assert_Sweep_Provenance( ...
        S.run_provenance,case_cfg,accepted_digests,campaign_tag,nd, ...
        start_date,result_dir);
    treatment = assert_result_identity(S.R,case_cfg,result_file);
    dates = daily_dates(S.Date,nd,start_date,result_file);
    if isempty(common_dates)
        common_dates = dates;
        daily.date_datenum = dates;
    elseif ~isequal(dates,common_dates)
        error('Bouchet_Joint_Sweep:DateMismatch', ...
            'Case %s does not share the campaign daily time axis.', ...
            case_cfg.tag);
    end

    if ~isnumeric(S.P) || ~isreal(S.P) || ~ismatrix(S.P) || ...
            size(S.P,1) ~= nd+1 || size(S.P,2) < 34
        error('Bouchet_Joint_Sweep:PoolShape', ...
            'Case %s has an invalid retained pool array.',case_cfg.tag);
    end
    A = Kantola_N_Magnitude_Audit(S,treatment,nd);
    bulk_pH = daily_state(S.PH_bulk_d,nd,'PH_bulk_d',false);
    pore_pH = daily_state(S.PH_pore_d,nd,'PH_pore_d',false);
    assert_pH(bulk_pH,'bulk',case_cfg.tag);
    assert_pH(pore_pH,'pore-water',case_cfg.tag);
    pH_period = pH_period_means(dates,bulk_pH);
    if strcmp(treatment,'control')
        pre_score = Pre_Control_pH_Score( ...
            dates,bulk_pH,observations.pH);
        if pH_period(1) ~= pre_score.model_bulk_pH_pre
            error('Bouchet_Joint_Sweep:PreControlOperator', ...
                'The exploratory and full pre-period operators differ.');
        end
    end
    pH_skill = pH_observation_skill( ...
        pH_period,observations.pH,treatment);
    pH_season = pH_season_means(dates,bulk_pH);
    pH_season_skill = pH_season_observation_skill( ...
        pH_season,observations.pH,treatment);

    n2o = daily_state(S.N2Oflx,nd,'N2Oflx',true);
    n2o_nit = daily_state(S.N2Oflx_nit,nd,'N2Oflx_nit',true);
    n2 = daily_state(S.N2flx,nd,'N2flx',true);
    denit = daily_state(S.DENITflx,nd,'DENITflx',true);
    partition_residual = nitrogen_partition( ...
        n2o,n2o_nit,n2,denit,case_cfg.tag);
    n2o_skill = n2o_observation_skill( ...
        dates,n2o,observations.n2o,treatment,case_cfg.tag);
    n2o_cumulative = n2o_cumulative_skill( ...
        dates,n2o,observations.n2o,treatment,case_cfg.tag);

    no3_leak = daily_state(S.LEAK_NO3,nd,'LEAK_NO3',true);
    no3_pool = double(S.P(2:nd+1,32));
    if ~isreal(no3_pool) || any(~isfinite(no3_pool)) || ...
            any(no3_pool < -1e-12)
        error('Bouchet_Joint_Sweep:InvalidNO3Pool', ...
            'Case %s contains a nonfinite or negative NO3 pool.', ...
            case_cfg.tag);
    end

    drainage = daily_drainage(S.Lk,nd,case_cfg.tag);
    wfps = daily_wfps(S.O,S.PAR,nd,case_cfg.tag);
    charge = charge_diagnostics(S.ERW,nd,case_cfg.tag);
    base_sat_season = pH_season_means(dates,charge.base_sat_daily);
    base_sat_skill = pH_season_observation_skill( ...
        base_sat_season,observations.base_sat,treatment);
    carbon = carbon_diagnostics(S,nd,case_cfg.tag);
    annual = annual_carbon(dates,carbon);
    carbon_skill = carbon_observation_skill( ...
        dates,carbon,observations.carbon,treatment,case_cfg.tag);

    daily.bulk_pH(:,k) = bulk_pH;
    daily.pore_pH(:,k) = pore_pH;
    daily.base_sat_pct(:,k) = charge.base_sat_daily;
    daily.N2O_gN_m2_d(:,k) = n2o;
    daily.NO3_leaching_gN_m2_d(:,k) = no3_leak;
    daily.NO3_pool_gN_m2(:,k) = no3_pool;
    daily.drainage_mm_d(:,k) = drainage;
    daily.WFPS_pct(:,k) = wfps;
    daily.GPP_gC_m2_d(:,k) = carbon.GPP;
    daily.Re_gC_m2_d(:,k) = carbon.Re;
    daily.NEE_gC_m2_d(:,k) = carbon.NEE;

    row = trajectory_row(case_cfg,treatment,result_file,S.run_provenance, ...
        bulk_pH,pore_pH,pH_period,pH_skill,pH_season,pH_season_skill, ...
        base_sat_season,base_sat_skill, ...
        n2o,n2o_skill,n2o_cumulative,A,no3_pool, ...
        drainage,wfps,partition_residual,charge,annual,carbon, ...
        carbon_skill,observation_fingerprint);
    if isempty(rows)
        rows = repmat(row,ncase,1);
    else
        rows(k) = row;
    end
    fprintf('[%d/%d] accepted %s\n',k,ncase,case_cfg.tag);
    clear S
end

T = struct2table(rows);
T.scored_under_validation_fingerprint = repmat(string(stamp_digest),height(T),1);
combination_table = summarize_combinations(T);
% Objectives: the model does not resolve the wetting pulses that dominate
% the chamber record, and the soil cores were taken in August (peak biomass,
% Beerling et al. 2024 Methods) on a day not recorded, so N2O is scored as
% the cumulative emission over each chamber campaign and pH and base
% saturation as the August mean of the lumped 0-30 cm box set against the
% observed 0-10 and 10-30 cm intervals separately and pooled; the
% basalt-minus-control difference isolates the treatment signal from the
% common-mode variability of the annual soil samples. The day-matched N2O NRMSE, the annual pH RMSE and the 0-30 cm
% hydrogen-activity-average pH scores stay in the tables as diagnostics. Each
% treatment-year is scored over its own chamber window; the basalt campaigns
% cover May to June only, so the basalt and control totals of one row are
% not comparable with each other.
objectives = [combination_table.N2O_cumulative_log_ratio_RMS, ...
    combination_table.pH_season_discrete_RMSE, ...
    combination_table.pH_difference_discrete_RMSE, ...
    combination_table.base_sat_difference_discrete_RMSE, ...
    combination_table.NO3_2017_relative_RMSE, ...
    combination_table.carbon_GPP_NEE_NRMSE];
if any(~isfinite(objectives),'all')
    error('Bouchet_Joint_Sweep:IncompleteObjectives', ...
        'At least one observational objective is missing or nonfinite.');
end
combination_table.Pareto_nondominated = pareto_minimize(objectives);
combination_table.Front_N2O_pH_season = pareto_minimize(objectives(:,[1 2]));
combination_table.Front_N2O_pH_difference = ...
    pareto_minimize(objectives(:,[1 3]));
combination_table.Front_pH_difference_base_sat_difference = ...
    pareto_minimize(objectives(:,[3 4]));
combination_table.objective_set = repmat("cumulative_N2O_discrete_depth_pH_base_sat_v11", ...
    height(combination_table),1);
combination_table.selection = repmat("manual_review_required", ...
    height(combination_table),1);

if ~strcmp(observation_identity(observations.files), ...
        observation_fingerprint)
    error('Bouchet_Joint_Sweep:ObservationsChanged', ...
        'An observation file changed during summarization.');
end
if ~strcmp(Sweep_Fingerprint( ...
        v2_dir,campaign_dir,campaign_tag),stored_stamp)
    error('Bouchet_Joint_Sweep:SourceChangedDuringSummary', ...
        'Code, inputs, tests or numerical runtime changed during summary.');
end

selection_notice = ...
    "No winner is selected automatically. Pareto_nondominated minimizes " ...
    + "the cumulative N2O log-ratio RMS over the 2017 and 2018 chamber " ...
    + "campaigns of both treatments, the August (arithmetic mean, the " ...
    + "sampling month) lumped bulk pH RMSE pooled over the " ...
    + "observed 0 to 10 cm and 10 to 30 cm intervals, the basalt-minus-" ...
    + "control August pH difference RMSE pooled over the same two " ...
    + "intervals, the basalt-minus-control August base saturation " ...
    + "difference RMSE pooled over the same two intervals (Fig. 1F), the " ...
    + "2017 NO3 relative RMSE and the independent GPP plus " ...
    + "NEE carbon error. The 0 to 30 cm hydrogen-activity average of the " ...
    + "two intervals stays in the tables as pH_season_RMSE and " ...
    + "pH_difference_RMSE. Front_N2O_pH_season, " ...
    + "Front_N2O_pH_difference and Front_pH_difference_base_sat_difference " ...
    + "flag the two-objective fronts. The " ...
    + "day-matched N2O NRMSE and the annual pH RMSE (pre period 1 January " ...
    + "to 9 November 2016, then calendar years; model pH averaged on ideal " ...
    + "hydrogen concentration with gamma_H = 1) remain reported as " ...
    + "diagnostics. Each N2O treatment-year is scored over its own chamber " ...
    + "window, carried in the tables; the basalt campaigns cover May to " ...
    + "June only. Re is reported and checked through " ...
    + "NEE = Re - GPP, but is not counted again in the carbon objective. " ...
    + "The last complete 2020 model day is 30 December because " ...
    + "the available forcing ends at 17:00 on 31 December.";
atomic_writetable(T,summary_file);
atomic_writetable(combination_table,combination_file);
atomic_save(mat_file,T,combination_table,inventory,daily,observations, ...
    observation_fingerprint,selection_notice);
atomic_write_text(status_file,sprintf( ...
    ['state=complete\ncomplete=%d\nexpected=%d\n' ...
     'pareto_nondominated=%d\nselection=none\n' ...
     'objective_set=cumulative_N2O_discrete_depth_pH_base_sat_v11\n' ...
     'front_N2O_pH_season=%d\nfront_N2O_pH_difference=%d\n' ...
     'validation_fingerprint=%s\naccepted_prior_digests=%s\n'], ...
    height(T),ncase,nnz(combination_table.Pareto_nondominated), ...
    nnz(combination_table.Front_N2O_pH_season), ...
    nnz(combination_table.Front_N2O_pH_difference), ...
    stamp_digest,accepted_prior));

fprintf('PASS: summarized 336 validated Bouchet joint-sweep trajectories.\n');
fprintf('No winner was selected; inspect the Pareto diagnostics.\n');
fprintf(['Objectives: cumulative N2O log-ratio RMS, August pH ' ...
    'RMSE, basalt-minus-control pH RMSE, 2017 NO3 relative RMSE, ' ...
    'GPP+NEE NRMSE (%d Pareto-nondominated tuples).\n'], ...
    nnz(combination_table.Pareto_nondominated));
fprintf('Fronts: %d tuples on (N2O, pH season), %d on (N2O, pH difference).\n', ...
    nnz(combination_table.Front_N2O_pH_season), ...
    nnz(combination_table.Front_N2O_pH_difference));
fprintf('Saved: %s\n',summary_file);
fprintf('Saved: %s\n',combination_file);
fprintf('Saved: %s\n',mat_file);

catch ME
    if ~strcmp(ME.identifier,'Bouchet_Joint_Sweep:IncompleteInventory')
        try
            atomic_write_text(status_file,sprintf( ...
                ['state=failed\ncomplete=%d\nexpected=336\n' ...
                 'selection=none\nerror_identifier=%s\n'], ...
                complete_count,ME.identifier));
        catch
        end
    end
    rethrow(ME)
end
end

function [T,summary_file] = summarize_pre_control( ...
        v2_dir,campaign_dir,campaign_tag,results_root,logs_root)
nd = 366;
start_date = datenum(2016,1,1); %#ok<DATNM>
results_root = fullfile(results_root,'pre_control_gate');
logs_root = fullfile(logs_root,'pre_control_gate');
if ~exist(logs_root,'dir'); mkdir(logs_root); end
summary_file = fullfile(logs_root,'pre_control_pH_discrepancy.csv');
inventory_file = fullfile(logs_root, ...
    'pre_control_completion_inventory.csv');
mat_file = fullfile(logs_root,'pre_control_pH_discrepancy.mat');
status_file = fullfile(logs_root,'pre_control_pH_discrepancy.status');
complete_count = 0;
try
atomic_write_text(status_file,sprintf( ...
    'state=running\ncomplete=0\nexpected=1\ncomparison=none\n'));

stamp_file = fullfile(fileparts(logs_root),'validation_stamp.txt');
if ~isfile(stamp_file)
    error('Bouchet_Joint_Sweep:MissingValidation', ...
        'The campaign validation stamp is missing: %s',stamp_file);
end
stored_stamp = strtrim(fileread(stamp_file));
current_stamp = Sweep_Fingerprint( ...
    v2_dir,campaign_dir,campaign_tag);
if ~strcmp(stored_stamp,current_stamp)
    error('Bouchet_Joint_Sweep:StaleValidation', ...
        'Code, tests, inputs or numerical runtime changed after validation.');
end
stamp_digest = Fingerprint_Digest(stored_stamp);

cases = Sweep_Cases;
case_cfg = cases(1);
if numel(cases) ~= 336 || case_cfg.id ~= 1 || case_cfg.run_id ~= 10 || ...
        case_cfg.AppeloAlpha ~= 1.3 || case_cfg.KNmaxd ~= 0.0048 || ...
        ~isnan(case_cfg.Kbot) || case_cfg.CEC_cmolc_kg ~= 13 || ...
        case_cfg.DissolutionFactor ~= 1
    error('Bouchet_Joint_Sweep:PreControlCase', ...
        'Case 1 is not the accepted CEC 13 control tuple.');
end
inventory = completion_inventory(case_cfg,results_root,logs_root);
atomic_writetable(inventory,inventory_file);
complete_count = nnz(inventory.complete);
if ~inventory.complete
    atomic_write_text(status_file,sprintf( ...
        ['state=incomplete\ncomplete=0\nexpected=1\n' ...
         'comparison=none\nproblem=%s\n'],char(inventory.problem)));
    error('Bouchet_Joint_Sweep:IncompletePreControl', ...
        'The one-year control trajectory is incomplete or ambiguous.');
end

observations = Load_Sweep_Observations(v2_dir);
observation_fingerprint = observation_identity(observations.files);
result_file = char(inventory.result_file(1));
required = {'run_provenance','R','Date','PH_bulk_d','PH_pore_d'};
contents = string({whos('-file',result_file).name});
missing = setdiff(string(required),contents);
if ~isempty(missing)
    error('Bouchet_Joint_Sweep:MissingField', ...
        'Result %s is missing: %s',result_file,char(strjoin(missing,', ')));
end
S = load(result_file,required{:});
result_dir = fullfile(results_root,case_cfg.tag);
Assert_Sweep_Provenance( ...
    S.run_provenance,case_cfg,stamp_digest,campaign_tag,nd, ...
    start_date,result_dir);
treatment = assert_result_identity(S.R,case_cfg,result_file);
if ~strcmp(treatment,'control')
    error('Bouchet_Joint_Sweep:PreControlTreatment', ...
        'The exploratory result is not the control treatment.');
end
dates = daily_dates(S.Date,nd,start_date,result_file);
bulk_pH = daily_state(S.PH_bulk_d,nd,'PH_bulk_d',false);
pore_pH = daily_state(S.PH_pore_d,nd,'PH_pore_d',false);
assert_pH(bulk_pH,'bulk',case_cfg.tag);
assert_pH(pore_pH,'pore-water',case_cfg.tag);
score = Pre_Control_pH_Score( ...
    dates,bulk_pH,observations.pH);

if ~strcmp(observation_identity(observations.files), ...
        observation_fingerprint)
    error('Bouchet_Joint_Sweep:ObservationsChanged', ...
        'An observation file changed during summarization.');
end
if ~strcmp(Sweep_Fingerprint( ...
        v2_dir,campaign_dir,campaign_tag),stored_stamp)
    error('Bouchet_Joint_Sweep:SourceChangedDuringSummary', ...
        'Code, inputs, tests or numerical runtime changed during summary.');
end

p = S.run_provenance;
row.case_id = case_cfg.id;
row.run_id = case_cfg.run_id;
row.treatment = string(treatment);
row.campaign = string(p.campaign);
row.exchange_scheme = string(p.exchange_scheme);
row.CEC_model = "fixed_total";
row.CEC_cmolc_kg = p.CEC_cmolc_kg;
row.CECpHSensitivity = p.CECpHSensitivity;
row.AppeloAlpha = p.AppeloAlpha;
row.AppeloLogKH_Na = p.AppeloLogKH(1);
row.AppeloLogKH_K = p.AppeloLogKH(2);
row.AppeloLogKH_Ca = p.AppeloLogKH(3);
row.AppeloLogKH_Mg = p.AppeloLogKH(4);
row.AppeloLogKH_Al = p.AppeloLogKH(5);
row.KNmaxd_per_day = p.KNmaxd;
row.Kbot_mm_h = p.Kbot;
row.Kbot_boundary = "free_drainage";
row.DissolutionFactor = p.dissolution_factor;
row.n_days = score.n_days;
row.n_pre_days = score.n_pre_days;
row.date_start = string(datestr(score.date_start,'yyyy-mm-dd'));
row.date_end = string(datestr(score.date_end,'yyyy-mm-dd'));
row.pre_date_end = string(datestr(score.pre_date_end,'yyyy-mm-dd'));
row.initial_bulk_pH = bulk_pH(1);
row.final_bulk_pH = bulk_pH(end);
row.minimum_bulk_pH = min(bulk_pH);
row.maximum_bulk_pH = max(bulk_pH);
row.initial_pore_pH = pore_pH(1);
row.final_pore_pH = pore_pH(end);
row.modeled_bulk_pH_pre = score.model_bulk_pH_pre;
row.observed_bulk_pH_pre_0_10 = score.observed_bulk_pH_pre_0_10;
row.observed_bulk_pH_pre_10_30 = score.observed_bulk_pH_pre_10_30;
row.observed_bulk_pH_pre_0_30 = score.observed_bulk_pH_pre;
row.pH_residual_model_minus_observed = score.residual;
row.pH_residual_0_10 = score.residual_0_10;
row.pH_residual_10_30 = score.residual_10_30;
row.pH_absolute_error = score.rmse;
row.pH_RMSE_one_comparison = score.rmse;
row.pH_bias_one_comparison = score.bias;
row.pH_comparison_count = 1;
row.pH_temporal_aggregation = "ideal_hydrogen_concentration";
row.validation_fingerprint = string(p.validation_fingerprint);
row.observation_fingerprint = string(observation_fingerprint);
row.matlab_release = string(p.matlab_release);
row.numerical_identity = string(p.numerical_identity);
row.result_file = string(result_file);
T = struct2table(row);

atomic_writetable(T,summary_file);
atomic_save_pre_control(mat_file,T,inventory,score,dates,bulk_pH, ...
    pore_pH,observations,observation_fingerprint);
atomic_write_text(status_file,sprintf( ...
    ['state=complete\ncomplete=1\nexpected=1\ncomparison=pre_control\n' ...
     'pH_residual_model_minus_observed=%.17g\n' ...
     'pH_absolute_error=%.17g\nvalidation_fingerprint=%s\n'], ...
    score.residual,score.rmse,stamp_digest));

fprintf('Model pre bulk pH: %.6f\n',score.model_bulk_pH_pre);
fprintf('Observed pre 0-30 cm bulk pH: %.6f\n', ...
    score.observed_bulk_pH_pre);
fprintf('Model minus observation: %+.6f pH\n',score.residual);
fprintf('Absolute error: %.6f pH\n',score.rmse);
fprintf('PASS: one-year control trajectory and manuscript pre-pH comparison\n');
recheck = strtrim(getenv('TC_PRE_CONTROL_RECHECK'));
if isempty(recheck) || strcmp(recheck,'0')
    fprintf('STOP: no pilot or production sweep was submitted\n');
elseif strcmp(recheck,'1')
    fprintf('PASS: reviewed exploratory result matches the current validation stamp\n');
else
    error('Bouchet_Joint_Sweep:BadRecheckFlag', ...
        'TC_PRE_CONTROL_RECHECK must be 0 or 1.');
end
fprintf('Saved: %s\n',summary_file);

catch ME
    if ~strcmp(ME.identifier,'Bouchet_Joint_Sweep:IncompletePreControl')
        try
            atomic_write_text(status_file,sprintf( ...
                ['state=failed\ncomplete=%d\nexpected=1\n' ...
                 'comparison=none\nerror_identifier=%s\n'], ...
                complete_count,ME.identifier));
        catch
        end
    end
    rethrow(ME)
end
end

function inventory = completion_inventory(cases,results_root,logs_root)
n = numel(cases);
row = struct('case_id',NaN,'combination_id',NaN,'case_tag',"", ...
    'run_id',NaN,'status_file',"",'status_state',"", ...
    'status_detail',"",'result_count',NaN,'result_file',"", ...
    'result_bytes',NaN,'complete',false,'problem',"");
rows = repmat(row,n,1);
for k = 1:n
    c = cases(k);
    result_dir = fullfile(results_root,c.tag);
    log_dir = fullfile(logs_root,c.tag);
    status_file = fullfile(log_dir,sprintf('run_%02d.status',c.run_id));
    files = dir(fullfile(result_dir,sprintf( ...
        'TC_2D_run%02d_*.mat',c.run_id)));
    files = files(~[files.isdir]);
    rows(k).case_id = c.id;
    rows(k).combination_id = c.combination_id;
    rows(k).case_tag = string(c.tag);
    rows(k).run_id = c.run_id;
    rows(k).status_file = string(status_file);
    rows(k).result_count = numel(files);
    if isscalar(files)
        result_file = fullfile(files.folder,files.name);
        rows(k).result_file = string(result_file);
        rows(k).result_bytes = files.bytes;
    else
        result_file = '';
    end
    [state,detail,status_problem] = read_status(status_file);
    rows(k).status_state = string(state);
    rows(k).status_detail = string(detail);
    problem = strings(0,1);
    if numel(files) ~= 1
        problem(end+1) = sprintf('found %d result files',numel(files)); %#ok<AGROW>
    end
    if ~isempty(status_problem)
        problem(end+1) = string(status_problem); %#ok<AGROW>
    elseif ~strcmp(state,'complete')
        problem(end+1) = "status state is not complete"; %#ok<AGROW>
    elseif isscalar(files) && ~strcmp(detail,result_file)
        problem(end+1) = "status detail does not name the sole result"; %#ok<AGROW>
    end
    rows(k).complete = isempty(problem);
    rows(k).problem = strjoin(problem,'; ');
end
inventory = struct2table(rows);
end

function [state,detail,problem] = read_status(file)
state = '';
detail = '';
problem = '';
if ~isfile(file)
    problem = 'status file is missing';
    return
end
lines = splitlines(string(fileread(file)));
lines = strip(lines);
lines(lines == "") = [];
keys = strings(size(lines));
values = strings(size(lines));
for k = 1:numel(lines)
    p = strfind(char(lines(k)),'=');
    if isempty(p)
        problem = 'status file contains a line without =';
        return
    end
    keys(k) = extractBefore(lines(k),p(1));
    values(k) = extractAfter(lines(k),p(1));
end
if nnz(keys == "state") ~= 1 || nnz(keys == "detail") ~= 1
    problem = 'status file must contain one state and one detail';
    return
end
state = char(values(keys == "state"));
detail = char(values(keys == "detail"));
end

function assert_result_stamp_accepted(p,case_tag,accepted,rejected)
% Refuse a result whose validation stamp is neither the current one nor an
% accepted prior one, and say why the prior stamp was not accepted.
if ~(isstruct(p) && isscalar(p) && isfield(p,'validation_fingerprint'))
    error('Bouchet_Joint_Sweep:ResultStampNotAccepted', ...
        'Case %s carries no validation fingerprint.',case_tag);
end
digest = char(string(p.validation_fingerprint));
if any(strcmp(digest,accepted))
    return
end
hit = strcmp(digest,{rejected.digest});
if any(hit)
    reason = rejected(find(hit,1)).reason;
else
    reason = ['it is not the current stamp and no archived stamp of that ' ...
        'name was offered through TC_ACCEPTED_RESULT_DIGESTS'];
end
error('Bouchet_Joint_Sweep:ResultStampNotAccepted', ...
    ['Case %s was validated under stamp %s, which is not accepted: %s. ' ...
     'Run the rescore stage, which offers every archived stamp.'], ...
    case_tag,digest,reason);
end

function treatment = assert_result_identity(R,case_cfg,result_file)
if ~(isstruct(R) && isscalar(R) && isfield(R,'id') && ...
        isfield(R,'treatment') && isequal(double(R.id),case_cfg.run_id))
    error('Bouchet_Joint_Sweep:WrongResult', ...
        'Result identity does not match case %s: %s', ...
        case_cfg.tag,result_file);
end
treatment = lower(char(R.treatment));
expected = {'control','basalt','calcite'};
if ~strcmp(treatment,expected{case_cfg.run_id-9})
    error('Bouchet_Joint_Sweep:WrongTreatment', ...
        'Case %s contains treatment %s.',case_cfg.tag,treatment);
end
end

function dates = daily_dates(x,nd,start_date,result_file)
if ~isnumeric(x) || ~isreal(x)
    error('Bouchet_Joint_Sweep:InvalidDate', ...
        'Hourly dates are not real numeric values in %s.',result_file);
end
x = double(x(:));
if numel(x) ~= 24*nd
    error('Bouchet_Joint_Sweep:ShortDate', ...
        'Hourly dates have the wrong length in %s.',result_file);
end
date_tol = 64*eps(max(abs(x)));
if any(~isfinite(x)) || any(abs(diff(x)-1/24) > date_tol)
    error('Bouchet_Joint_Sweep:InvalidDate', ...
        'Hourly dates are nonfinite or do not advance hourly in %s.', ...
        result_file);
end
dates = floor(x(1:24:end));
expected = start_date+(0:nd-1).';
if numel(dates) ~= nd || ~isequal(dates,expected)
    error('Bouchet_Joint_Sweep:InvalidDate', ...
        'Hourly dates do not map to the expected %d campaign days in %s.', ...
        nd,result_file);
end
%% Retained row r holds the BIOGEO day that starts at 01:00 of
%% start_date+(r-3), and daily_state keeps rows 2..nd+1, so those rows are
%% dated one day earlier than the forcing days. Row 2 is the lone first hour.
dates = dates-1;
end

function x = daily_state(x,nd,name,nonnegative)
if ~isnumeric(x) || ~isreal(x) || ~ismatrix(x) || ...
        size(x,1) ~= nd+1 || size(x,2) ~= 1
    error('Bouchet_Joint_Sweep:ShortDailyState', ...
        '%s must contain one real column through retained day %d.',name,nd);
end
x = double(x(2:nd+1,1));
if any(~isfinite(x)) || (nonnegative && any(x < -1e-12))
    error('Bouchet_Joint_Sweep:InvalidDailyState', ...
        '%s is nonfinite or outside its physical range.',name);
end
end

function assert_pH(x,label,case_tag)
if any(~isfinite(x)) || any(x <= 0 | x >= 14)
    error('Bouchet_Joint_Sweep:InvalidPH', ...
        'Case %s contains invalid %s pH.',case_tag,label);
end
end

function value = pH_period_means(dates,pH)
dv = datevec(dates);
value = nan(1,5);
pick = dates >= datenum(2016,1,1) & dates < datenum(2016,11,10); %#ok<DATNM>
value(1) = -log10(mean(10.^(-pH(pick))));
for k = 1:4
    pick = dv(:,1) == 2016+k;
    value(k+1) = -log10(mean(10.^(-pH(pick))));
end
if any(~isfinite(value))
    error('Bouchet_Joint_Sweep:MissingPHPeriod', ...
        'The model does not contain every required pH comparison period.');
end
end

function skill = pH_observation_skill(model,obs,treatment)
skill = empty_skill;
if strcmp(treatment,'calcite'); return; end
target = obs.(treatment).profile_0_30;
residual = model(:)-target(:);
skill = residual_skill(residual,target);
skill.target = target(:).';
end

function residual = nitrogen_partition(n2o,n2o_nit,n2,denit,case_tag)
denit_n2o = n2o-n2o_nit;
tol = 1e-12;
if any(denit_n2o < -tol)
    error('Bouchet_Joint_Sweep:N2OComponents', ...
        'Case %s has nitrification N2O above total N2O.',case_tag);
end
residual = denit-n2-denit_n2o;
if max(abs(residual)) > 1e-10
    error('Bouchet_Joint_Sweep:NPartition', ...
        'Case %s has a %.12g g N m^-2 d^-1 partition residual.', ...
        case_tag,max(abs(residual)));
end
end

function skill = n2o_observation_skill(dates,model,obs,treatment,case_tag)
skill = empty_skill;
skill.model_mean = NaN;
skill.observed_mean = NaN;
if strcmp(treatment,'calcite'); return; end
pick = obs.treatment == string(treatment);
[found,index] = ismember(obs.date_datenum(pick),dates);
if ~all(found)
    error('Bouchet_Joint_Sweep:N2ODateMismatch', ...
        'Case %s does not cover every matched N2O observation date.', ...
        case_tag);
end
model_value = 1e4*model(index);
target = obs.flux_gN_ha_d(pick);
residual = model_value-target;
skill = residual_skill(residual,target);
skill.model_mean = mean(model_value);
skill.observed_mean = mean(target);
end

function season = pH_season_means(dates,pH)
% August mean of a daily series per comparison period. The soil cores were
% taken inside the peak-biomass quadrats, collected in mid August (Beerling
% et al. 2024, Methods; SI Fig. S7), with no day recorded, so the model is
% averaged arithmetically over August of 2016 (pre-application) and of 2017
% to 2020, and its daily range over the month is retained as the
% sampling-date envelope. Used for bulk pH and for base saturation.
dv = datevec(dates);
season.years = 2016:2020;
season.mean = nan(1,5);
season.min = nan(1,5);
season.max = nan(1,5);
season.count = zeros(1,5);
for k = 1:5
    pick = dv(:,1) == season.years(k) & dv(:,2) == 8;
    if nnz(pick) ~= 31
        error('Bouchet_Joint_Sweep:MissingPHSeason', ...
            'The model holds %d days of August %d, not 31.', ...
            nnz(pick),season.years(k));
    end
    season.mean(k) = mean(pH(pick));
    season.min(k) = min(pH(pick));
    season.max(k) = max(pH(pick));
    season.count(k) = nnz(pick);
end
end

function skill = pH_season_observation_skill(season,obs,treatment)
% The model is one lumped 0-30 cm box (Zbio); the field was sampled at 0-10
% and 10-30 cm and holds a gradient the box cannot. The lumped
% August mean is therefore scored against each observed interval on
% its own and pooled over both (ten residuals per treatment); the pooled
% score is the pH objective. The 0-30 cm hydrogen-activity average of the
% two intervals, the earlier single target, stays as the top-level
% diagnostic. The soil samples carry no date within their year, so every
% target also gets the timing-agnostic envelope miss: the distance of the
% observation outside the model's daily range over the same window, zero
% when a sampling date exists on which the model matches.
skill = empty_skill;
skill.envelope_miss_sse = NaN;
skill.envelope_miss_rms = NaN;
skill.shallow = empty_skill;
skill.shallow.envelope_miss_sse = NaN;
skill.deep = empty_skill;
skill.deep.envelope_miss_sse = NaN;
skill.discrete = empty_skill;
skill.discrete.envelope_miss_sse = NaN;
skill.discrete.envelope_miss_rms = NaN;
if strcmp(treatment,'calcite'); return; end
target = obs.(treatment).profile_0_30(:).';
shallow = obs.(treatment).shallow_0_10(:).';
deep = obs.(treatment).deep_10_30(:).';
skill = residual_skill(season.mean-target,target);
skill.target = target;
miss = max(0,target-season.max)+max(0,season.min-target);
skill.envelope_miss_sse = sum(miss.^2);
skill.envelope_miss_rms = sqrt(mean(miss.^2));
skill.shallow = residual_skill(season.mean-shallow,shallow);
skill.shallow.target = shallow;
skill.deep = residual_skill(season.mean-deep,deep);
skill.deep.target = deep;
skill.discrete = residual_skill( ...
    [season.mean-shallow, season.mean-deep],[shallow, deep]);
miss_shallow = max(0,shallow-season.max)+max(0,season.min-shallow);
miss_deep = max(0,deep-season.max)+max(0,season.min-deep);
skill.shallow.envelope_miss_sse = sum(miss_shallow.^2);
skill.deep.envelope_miss_sse = sum(miss_deep.^2);
skill.discrete.envelope_miss_sse = sum([miss_shallow miss_deep].^2);
skill.discrete.envelope_miss_rms = ...
    sqrt(mean([miss_shallow miss_deep].^2));
end

function skill = n2o_cumulative_skill(dates,model,obs,treatment,case_tag)
% Cumulative N2O per chamber campaign year. The observed total is the
% trapezoid over the chamber dates, the convention of the field studies,
% and the model is integrated by the same rule over the same window, so the
% comparison does not depend on the timing of individual wetting pulses.
skill.years = [2017 2018];
skill.observed_kgN_ha = nan(1,2);
skill.modeled_kgN_ha = nan(1,2);
skill.log_ratio = nan(1,2);
skill.observation_count = zeros(1,2);
skill.window_start = strings(1,2);
skill.window_end = strings(1,2);
skill.log_ratio_ss = NaN;
skill.count = 0;
skill.log_ratio_rms = NaN;
if strcmp(treatment,'calcite'); return; end
model = 1e4*double(model(:));           % g N m^-2 d^-1 -> g N ha^-1 d^-1
observed_years = datevec(obs.date_datenum);
for k = 1:2
    pick = obs.treatment == string(treatment) & ...
        observed_years(:,1) == skill.years(k);
    od = obs.date_datenum(pick); of = obs.flux_gN_ha_d(pick);
    [od,order] = sort(od(:)); of = double(of(order));
    if numel(od) < 2 || any(diff(od) <= 0)
        error('Bouchet_Joint_Sweep:N2OCumulativeWindow', ...
            'Case %s has fewer than two distinct %d %s chamber dates.', ...
            case_tag,skill.years(k),treatment);
    end
    window = dates >= od(1) & dates <= od(end);
    if nnz(window) ~= od(end)-od(1)+1
        error('Bouchet_Joint_Sweep:N2OCumulativeWindow', ...
            'Case %s does not cover the %d %s chamber window.', ...
            case_tag,skill.years(k),treatment);
    end
    observed = trapz(od,of)/1000;         % g N ha^-1 -> kg N ha^-1
    modeled = trapz(dates(window),model(window))/1000;
    if ~(isfinite(observed) && observed > 0 && ...
            isfinite(modeled) && modeled > 0)
        error('Bouchet_Joint_Sweep:N2OCumulativeNonpositive', ...
            'Case %s has a nonpositive cumulative N2O in %d (%s).', ...
            case_tag,skill.years(k),treatment);
    end
    skill.observed_kgN_ha(k) = observed;
    skill.modeled_kgN_ha(k) = modeled;
    skill.log_ratio(k) = log(modeled/observed);
    skill.observation_count(k) = numel(od);
    skill.window_start(k) = string(datestr(od(1),'yyyy-mm-dd'));
    skill.window_end(k) = string(datestr(od(end),'yyyy-mm-dd'));
end
skill.log_ratio_ss = sum(skill.log_ratio.^2);
skill.count = 2;
skill.log_ratio_rms = sqrt(skill.log_ratio_ss/skill.count);
end

function label = cumulative_window_label(skill,k)
if strlength(skill.window_start(k)) == 0
    label = "";
else
    label = skill.window_start(k)+"_to_"+skill.window_end(k);
end
end

function drainage = daily_drainage(x,nd,case_tag)
if ~isnumeric(x) || ~isreal(x) || ~ismatrix(x) || ...
        size(x,1) ~= 24*nd || size(x,2) ~= 1
    error('Bouchet_Joint_Sweep:DrainageShape', ...
        'Case %s does not contain one complete drainage column.',case_tag);
end
x = double(x);
x = x(:,1);
if any(~isfinite(x)) || any(x < -1e-12)
    error('Bouchet_Joint_Sweep:InvalidDrainage', ...
        'Case %s contains nonfinite or negative drainage.',case_tag);
end
drainage = [x(1); sum(reshape(x(2:24*(nd-1)+1),24,nd-1),1).'];
end

function wfps = daily_wfps(O,PAR,nd,case_tag)
if ~isnumeric(O) || ~isreal(O)
    error('Bouchet_Joint_Sweep:MoistureShape', ...
        'Case %s has a nonnumeric or complex moisture state.',case_tag);
end
O = double(O);
if ndims(O) == 3
    if size(O,3) ~= 1
        error('Bouchet_Joint_Sweep:MoistureColumns', ...
            'Case %s contains more than one moisture column.',case_tag);
    end
    O = O(:,:,1);
end
if ~(ismatrix(O) && size(O,1) == 24*nd && all(isfinite(O),'all')) || ...
        ~(iscell(PAR) && isscalar(PAR))
    error('Bouchet_Joint_Sweep:MoistureShape', ...
        'Case %s has an incomplete hourly moisture state.',case_tag);
end
soil = PAR{1};
Bio_Zs = double(soil.Bio_Zs(:));
Osat = double(soil.Soil_Param.Osat(:));
if size(O,2) ~= numel(Bio_Zs) || numel(Osat) ~= numel(Bio_Zs) || ...
        ~isreal(Bio_Zs) || ~isreal(Osat) || ...
        any(~isfinite(Bio_Zs)) || any(Bio_Zs < 0) || ...
        abs(sum(Bio_Zs)-1) > 1e-10 || any(~isfinite(Osat)) || ...
        any(Osat <= 0)
    error('Bouchet_Joint_Sweep:MoistureLayers', ...
        'Case %s has inconsistent active-layer water parameters.',case_tag);
end
sat_water = sum(Bio_Zs.*Osat);
hourly = 100*(O*Bio_Zs)/sat_water;
if any(hourly < -1e-10 | hourly > 100+1e-10)
    error('Bouchet_Joint_Sweep:UnphysicalWFPS', ...
        'Case %s has WFPS outside 0 to 100 percent.',case_tag);
end
wfps = [hourly(1); mean(reshape(hourly(2:24*(nd-1)+1),24,nd-1),1).'];
end

function out = charge_diagnostics(ERW,nd,case_tag)
if ~iscell(ERW) || numel(ERW) ~= 1 || ~isfield(ERW{1},'bg') || ...
        ~isfield(ERW{1},'const')
    error('Bouchet_Joint_Sweep:MissingERW', ...
        'Case %s lacks its one-column ERW state.',case_tag);
end
bg = ERW{1}.bg;
required = {'f_CEC','charge_residual','Al_inventory_residual'};
if ~all(isfield(bg,required)) || ~isfield(ERW{1}.const,'Al_CEC_scale')
    error('Bouchet_Joint_Sweep:MissingERWField', ...
        'Case %s lacks a required ERW closure field.',case_tag);
end
f = double(bg.f_CEC);
charge = double(bg.charge_residual(:));
al = double(bg.Al_inventory_residual(:));
scale = double(ERW{1}.const.Al_CEC_scale);
% The weathering module retains its closure state per weathering step: a
% daily layout holds nd+1 rows with the initial state first, a sub-daily
% layout holds an integer number of steps per retained day (the campaign
% results hold 144 ten-minute steps per day). Both are accepted and the
% tolerances are applied to every retained step.
nrow = size(f,1);
if nrow == nd+1
    rows = 2:nd+1;
    steps_per_day = 1;
elseif nrow > nd && mod(nrow,nd) == 0
    rows = 1:nrow;
    steps_per_day = nrow/nd;
else
    rows = [];
    steps_per_day = NaN;
end
if ~isreal(f) || ~isreal(charge) || ~isreal(al) || ~isreal(scale) || ...
        ~ismatrix(f) || isempty(rows) || size(f,2) ~= 6 || ...
        numel(charge) ~= nrow || numel(al) ~= nrow || ...
        any(~isfinite(f),'all') || ...
        any(f < -1e-12 | f > 1+1e-12,'all') || ...
        ~(isscalar(scale) && isfinite(scale) && scale > 0) || ...
        isempty(charge) || any(~isfinite(charge)) || ...
        isempty(al) || any(~isfinite(al))
    error('Bouchet_Joint_Sweep:InvalidERWClosure', ...
        ['Case %s contains a nonfinite or unphysical ERW closure state, ' ...
         'or its %d closure rows do not match %d retained days.'], ...
        case_tag,nrow,nd);
end
f = f(rows,:);
charge = charge(rows);
al = al(rows);
out.erw_steps_per_day = steps_per_day;
% Daily base saturation (Ca, Mg, K, Na sites) for the Fig. 1F score.
out.base_sat_daily = mean(reshape(100*sum(f(:,1:4),2),steps_per_day,[]),1).';
out.max_exchange_sum_residual = max(abs(sum(f,2)-1));
out.max_charge_residual_umolc_L = max(abs(charge));
out.max_Al_residual_umol_m2 = max(abs(al))/scale;
if out.max_exchange_sum_residual > 1e-8 || ...
        out.max_charge_residual_umolc_L > 1e-6 || ...
        out.max_Al_residual_umol_m2 > 1e-6
    error('Bouchet_Joint_Sweep:ERWClosure', ...
        ['Case %s fails exchange, charge or Al closure: ' ...
         '%.3g, %.3g, %.3g.'],case_tag, ...
        out.max_exchange_sum_residual, ...
        out.max_charge_residual_umolc_L, ...
        out.max_Al_residual_umol_m2);
end
end

function out = carbon_diagnostics(S,nd,case_tag)
area = double(S.topo.cell_area_m2(:));
ncol = numel(area);
nstate = nd+1;
cc = size(S.Ccrown_t,2);
plant_names = {'RA_H','RA_L','NPP_H','NPP_L','Ccrown_t'};
soil_names = {'R_litter','R_microbe','R_ew'};
if ~isreal(area) || ncol < 1 || any(~isfinite(area)) || ...
        any(area <= 0) || ...
        numel(S.ColCfg) ~= ncol || cc < 1
    error('Bouchet_Joint_Sweep:CarbonGeometry', ...
        'Case %s has invalid carbon geometry.',case_tag);
end
for k = 1:numel(plant_names)
    x = S.(plant_names{k});
    if ~isnumeric(x) || ~isreal(x) || size(x,1) ~= nstate || ...
            numel(x) ~= nstate*cc*ncol || any(~isfinite(x(:)))
        error('Bouchet_Joint_Sweep:CarbonShape', ...
            'Case %s has invalid %s dimensions.',case_tag,plant_names{k});
    end
end
for k = 1:numel(soil_names)
    x = S.(soil_names{k});
    if ~isnumeric(x) || ~isreal(x) || size(x,1) ~= nstate || ...
            numel(x) ~= nstate*ncol || any(~isfinite(x(:)))
        error('Bouchet_Joint_Sweep:CarbonShape', ...
            'Case %s has invalid %s dimensions.',case_tag,soil_names{k});
    end
end

crown = reshape(double(S.Ccrown_t),nstate,cc,ncol);
RA_H = reshape(double(S.RA_H),nstate,cc,ncol);
RA_L = reshape(double(S.RA_L),nstate,cc,ncol);
NPP_H = reshape(double(S.NPP_H),nstate,cc,ncol);
NPP_L = reshape(double(S.NPP_L),nstate,cc,ncol);
R_litter = reshape(double(S.R_litter),nstate,ncol);
R_microbe = reshape(double(S.R_microbe),nstate,ncol);
R_ew = reshape(double(S.R_ew),nstate,ncol);
if any(crown < -1e-12 | crown > 1+1e-12,'all') || ...
        any(RA_H < -1e-12,'all') || any(RA_L < -1e-12,'all') || ...
        any(R_litter < -1e-12,'all') || ...
        any(R_microbe < -1e-12,'all') || any(R_ew < -1e-12,'all')
    error('Bouchet_Joint_Sweep:UnphysicalCarbon', ...
        'Case %s contains an unphysical crown or respiration state.', ...
        case_tag);
end

weights = area/sum(area);
GPP = zeros(nd,1);
Re = zeros(nd,1);
NEE_direct = zeros(nd,1);
for c = 1:ncol
    w = crown(2:end,:,c);
    ra = RA_H(2:end,:,c)+RA_L(2:end,:,c);
    npp = NPP_H(2:end,:,c)+NPP_L(2:end,:,c);
    hetero = R_litter(2:end,c)+R_microbe(2:end,c)+R_ew(2:end,c);
    GPP = GPP + weights(c)*sum(w.*(npp+ra),2);
    Re = Re + weights(c)*(hetero+sum(w.*ra,2));
    NEE_direct = NEE_direct + weights(c)*(hetero-sum(w.*npp,2));
end
scale = max([1;abs(GPP);abs(Re);abs(NEE_direct)]);
tol = 100*eps(scale);
if any(~isfinite([GPP;Re;NEE_direct])) || ...
        any(GPP < -tol) || any(Re < -tol)
    error('Bouchet_Joint_Sweep:UnphysicalCarbon', ...
        'Case %s has nonfinite or negative GPP or Re.',case_tag);
end
GPP(GPP < 0) = 0;
Re(Re < 0) = 0;
NEE = Re-GPP;
identity_residual = NEE-NEE_direct;
if max(abs(identity_residual)) > tol
    error('Bouchet_Joint_Sweep:CarbonIdentity', ...
        'Case %s violates NEE = Re - GPP by %.12g g C m^-2 d^-1.', ...
        case_tag,max(abs(identity_residual)));
end
out.GPP = GPP;
out.Re = Re;
out.NEE = NEE;
out.max_identity_residual = max(abs(identity_residual));
end

function out = annual_carbon(dates,carbon)
years = [2017 2018 2019 2020];
dv = datevec(dates);
out.years = years;
out.n_days = zeros(size(years));
out.GPP = nan(size(years));
out.Re = nan(size(years));
out.NEE = nan(size(years));
for k = 1:numel(years)
    pick = dv(:,1) == years(k);
    out.n_days(k) = nnz(pick);
    out.GPP(k) = sum(carbon.GPP(pick));
    out.Re(k) = sum(carbon.Re(pick));
    out.NEE(k) = sum(carbon.NEE(pick));
end
%% 29 December is the last complete day of 2020 the forcing holds.
if ~isequal(out.n_days,[365 365 365 364]) || ...
        any(~isfinite([out.GPP out.Re out.NEE]))
    error('Bouchet_Joint_Sweep:AnnualCarbon', ...
        'A required carbon comparison window is absent or nonfinite.');
end
end

function skill = carbon_observation_skill( ...
        dates,carbon,obs,treatment,case_tag)
skill = struct('GPP_SSE',NaN,'NEE_SSE',NaN,'Re_SSE',NaN, ...
    'independent_SSE',NaN,'count',0,'independent_count',0, ...
    'GPP_NRMSE',NaN,'NEE_NRMSE',NaN,'Re_NRMSE',NaN, ...
    'independent_NRMSE',NaN,'target_GPP',nan(1,4), ...
    'target_NEE',nan(1,4),'target_Re',nan(1,4));
if strcmp(treatment,'calcite'); return; end
years = [2017 2018 2019 2020];
crops = {'maize','maize','soybean','maize'};
dv = datevec(dates);
eGPP = nan(4,1);
eNEE = nan(4,1);
eRe = nan(4,1);
for k = 1:numel(years)
    pick = dv(:,1) == years(k);
    index = observation_day_index(dates(pick),years(k));
    reference = obs.(crops{k}).(treatment);
    skill.target_GPP(k) = sum(reference.GPP(index));
    skill.target_NEE(k) = sum(reference.NEE(index));
    skill.target_Re(k) = sum(reference.Re(index));
    scale = abs(skill.target_GPP(k));
    if ~(isfinite(scale) && scale > 0)
        error('Bouchet_Joint_Sweep:CarbonObservationScale', ...
            'Case %s has an invalid observed GPP scale.',case_tag);
    end
    eGPP(k) = (sum(carbon.GPP(pick))-skill.target_GPP(k))/scale;
    eNEE(k) = (sum(carbon.NEE(pick))-skill.target_NEE(k))/scale;
    eRe(k) = (sum(carbon.Re(pick))-skill.target_Re(k))/scale;
end
skill.GPP_SSE = sum(eGPP.^2);
skill.NEE_SSE = sum(eNEE.^2);
skill.Re_SSE = sum(eRe.^2);
skill.independent_SSE = skill.GPP_SSE+skill.NEE_SSE;
skill.count = numel(eGPP);
skill.independent_count = 2*numel(eGPP);
skill.GPP_NRMSE = sqrt(skill.GPP_SSE/skill.count);
skill.NEE_NRMSE = sqrt(skill.NEE_SSE/skill.count);
skill.Re_NRMSE = sqrt(skill.Re_SSE/skill.count);
skill.independent_NRMSE = sqrt( ...
    skill.independent_SSE/skill.independent_count);
end

function index = observation_day_index(date_value,year_value)
index = round(date_value-datenum(year_value,1,1)+1); %#ok<DATNM>
if any(index < 1 | index > 366)
    error('Bouchet_Joint_Sweep:CarbonObservationDate', ...
        'A model date cannot be mapped to the 366-day observation curve.');
end
end

function row = trajectory_row(c,treatment,result_file,p,bulk_pH,pore_pH, ...
        pH_period,pH_skill,pH_season,pH_season_skill, ...
        base_sat_season,base_sat_skill,n2o,n2o_skill, ...
        n2o_cumulative,A,no3_pool,drainage,wfps, ...
        partition_residual,charge,annual,carbon,carbon_skill, ...
        observation_fingerprint)
row.case_id = c.id;
row.combination_id = c.combination_id;
row.case_tag = string(c.tag);
row.combination_tag = string(c.combination_tag);
row.run_id = c.run_id;
row.treatment = string(treatment);
row.campaign = string(p.campaign);
row.exchange_scheme = string(p.exchange_scheme);
row.CEC_model = "fixed_total";
row.matlab_release = string(p.matlab_release);
row.numerical_identity = string(p.numerical_identity);
row.sweep_layout = string(p.sweep_layout);
row.AppeloAlpha = c.AppeloAlpha;
row.AppeloLogKH_Na = p.AppeloLogKH(1);
row.AppeloLogKH_K = p.AppeloLogKH(2);
row.AppeloLogKH_Ca = p.AppeloLogKH(3);
row.AppeloLogKH_Mg = p.AppeloLogKH(4);
row.AppeloLogKH_Al = p.AppeloLogKH(5);
row.CEC_cmolc_kg = p.CEC_cmolc_kg;
row.CECpHSensitivity = p.CECpHSensitivity;
row.KNmaxd_per_day = c.KNmaxd;
row.Kbot_mm_h = c.Kbot;
if isnan(c.Kbot); row.Kbot_boundary = "free_drainage";
else; row.Kbot_boundary = "finite_conductivity"; end
row.DissolutionFactor = c.DissolutionFactor;
row.result_file = string(result_file);
row.validation_fingerprint = string(p.validation_fingerprint);
row.observation_fingerprint = string(observation_fingerprint);
row.slurm_job_id = string(p.slurm_job_id);
row.slurm_restart_count = string(p.slurm_restart_count);
row.initial_bulk_pH = bulk_pH(1);
row.final_bulk_pH = bulk_pH(end);
row.initial_pore_pH = pore_pH(1);
row.final_pore_pH = pore_pH(end);
row.bulk_pH_pre = pH_period(1);
row.bulk_pH_2017 = pH_period(2);
row.bulk_pH_2018 = pH_period(3);
row.bulk_pH_2019 = pH_period(4);
row.bulk_pH_2020 = pH_period(5);
row.bulk_pH_temporal_aggregation = "ideal_hydrogen_concentration";
row.observed_bulk_pH_pre = pH_skill.target(1);
row.observed_bulk_pH_2017 = pH_skill.target(2);
row.observed_bulk_pH_2018 = pH_skill.target(3);
row.observed_bulk_pH_2019 = pH_skill.target(4);
row.observed_bulk_pH_2020 = pH_skill.target(5);
row.pH_observation_RMSE = pH_skill.rmse;
row.pH_observation_bias = pH_skill.bias;
row.pH_observation_SSE = pH_skill.sse;
row.pH_observation_count = pH_skill.count;
row.bulk_pH_season_window = "Aug_01_to_Aug_31_arithmetic_mean";
row.bulk_pH_season_pre = pH_season.mean(1);
row.bulk_pH_season_2017 = pH_season.mean(2);
row.bulk_pH_season_2018 = pH_season.mean(3);
row.bulk_pH_season_2019 = pH_season.mean(4);
row.bulk_pH_season_2020 = pH_season.mean(5);
row.bulk_pH_season_min_pre = pH_season.min(1);
row.bulk_pH_season_min_2017 = pH_season.min(2);
row.bulk_pH_season_min_2018 = pH_season.min(3);
row.bulk_pH_season_min_2019 = pH_season.min(4);
row.bulk_pH_season_min_2020 = pH_season.min(5);
row.bulk_pH_season_max_pre = pH_season.max(1);
row.bulk_pH_season_max_2017 = pH_season.max(2);
row.bulk_pH_season_max_2018 = pH_season.max(3);
row.bulk_pH_season_max_2019 = pH_season.max(4);
row.bulk_pH_season_max_2020 = pH_season.max(5);
row.pH_season_RMSE = pH_season_skill.rmse;
row.pH_season_bias = pH_season_skill.bias;
row.pH_season_SSE = pH_season_skill.sse;
row.pH_season_count = pH_season_skill.count;
row.pH_envelope_miss_RMS = pH_season_skill.envelope_miss_rms;
row.pH_envelope_miss_SSE = pH_season_skill.envelope_miss_sse;
row.observed_bulk_pH_season_0_10_pre = pH_season_skill.shallow.target(1);
row.observed_bulk_pH_season_0_10_2017 = pH_season_skill.shallow.target(2);
row.observed_bulk_pH_season_0_10_2018 = pH_season_skill.shallow.target(3);
row.observed_bulk_pH_season_0_10_2019 = pH_season_skill.shallow.target(4);
row.observed_bulk_pH_season_0_10_2020 = pH_season_skill.shallow.target(5);
row.observed_bulk_pH_season_10_30_pre = pH_season_skill.deep.target(1);
row.observed_bulk_pH_season_10_30_2017 = pH_season_skill.deep.target(2);
row.observed_bulk_pH_season_10_30_2018 = pH_season_skill.deep.target(3);
row.observed_bulk_pH_season_10_30_2019 = pH_season_skill.deep.target(4);
row.observed_bulk_pH_season_10_30_2020 = pH_season_skill.deep.target(5);
row.pH_season_0_10_RMSE = pH_season_skill.shallow.rmse;
row.pH_season_0_10_bias = pH_season_skill.shallow.bias;
row.pH_season_0_10_SSE = pH_season_skill.shallow.sse;
row.pH_season_10_30_RMSE = pH_season_skill.deep.rmse;
row.pH_season_10_30_bias = pH_season_skill.deep.bias;
row.pH_season_10_30_SSE = pH_season_skill.deep.sse;
row.pH_season_discrete_RMSE = pH_season_skill.discrete.rmse;
row.pH_season_discrete_bias = pH_season_skill.discrete.bias;
row.pH_season_discrete_SSE = pH_season_skill.discrete.sse;
row.pH_season_discrete_count = pH_season_skill.discrete.count;
row.pH_envelope_miss_0_10_SSE = pH_season_skill.shallow.envelope_miss_sse;
row.pH_envelope_miss_10_30_SSE = pH_season_skill.deep.envelope_miss_sse;
row.pH_envelope_miss_discrete_SSE = pH_season_skill.discrete.envelope_miss_sse;
row.pH_envelope_miss_discrete_RMS = pH_season_skill.discrete.envelope_miss_rms;
%% Base saturation (Fig. 1F) scored as the pH is: the August mean of the
%% lumped box against each sampled interval and pooled.
period = {'pre','2017','2018','2019','2020'};
for y = 1:5
    row.(['base_sat_season_' period{y}]) = base_sat_season.mean(y);
    row.(['observed_base_sat_0_10_' period{y}]) = base_sat_skill.shallow.target(y);
    row.(['observed_base_sat_10_30_' period{y}]) = base_sat_skill.deep.target(y);
end
row.base_sat_season_RMSE = base_sat_skill.rmse;
row.base_sat_season_0_10_RMSE = base_sat_skill.shallow.rmse;
row.base_sat_season_10_30_RMSE = base_sat_skill.deep.rmse;
row.base_sat_season_discrete_RMSE = base_sat_skill.discrete.rmse;
row.base_sat_season_discrete_bias = base_sat_skill.discrete.bias;
row.base_sat_season_discrete_SSE = base_sat_skill.discrete.sse;
row.base_sat_season_discrete_count = base_sat_skill.discrete.count;
row.cumulative_N2O_gN_m2 = sum(n2o);
row.cumulative_N2O_kgN_ha = 10*sum(n2o);
row.N2O_observation_model_mean_gN_ha_d = n2o_skill.model_mean;
row.N2O_observation_mean_gN_ha_d = n2o_skill.observed_mean;
row.N2O_observation_RMSE_gN_ha_d = n2o_skill.rmse;
row.N2O_observation_NRMSE = n2o_skill.nrmse;
row.N2O_observation_bias_gN_ha_d = n2o_skill.bias;
row.N2O_observation_SSE = n2o_skill.sse;
row.N2O_observation_target_sumsq = n2o_skill.target_sumsq;
row.N2O_observation_count = n2o_skill.count;
row.N2O_cumulative_2017_window = cumulative_window_label(n2o_cumulative,1);
row.N2O_cumulative_2017_observation_count = ...
    n2o_cumulative.observation_count(1);
row.N2O_cumulative_2017_observed_kgN_ha = n2o_cumulative.observed_kgN_ha(1);
row.N2O_cumulative_2017_modeled_kgN_ha = n2o_cumulative.modeled_kgN_ha(1);
row.N2O_cumulative_2017_log_ratio = n2o_cumulative.log_ratio(1);
row.N2O_cumulative_2018_window = cumulative_window_label(n2o_cumulative,2);
row.N2O_cumulative_2018_observation_count = ...
    n2o_cumulative.observation_count(2);
row.N2O_cumulative_2018_observed_kgN_ha = n2o_cumulative.observed_kgN_ha(2);
row.N2O_cumulative_2018_modeled_kgN_ha = n2o_cumulative.modeled_kgN_ha(2);
row.N2O_cumulative_2018_log_ratio = n2o_cumulative.log_ratio(2);
row.N2O_cumulative_log_ratio_SS = n2o_cumulative.log_ratio_ss;
row.N2O_cumulative_count = n2o_cumulative.count;
row.N2O_cumulative_log_ratio_RMS = n2o_cumulative.log_ratio_rms;
row.cumulative_NO3_leaching_kgN_ha = A.cumulative_NO3_leak_kgNha;
row.modeled_2017_NO3_leaching_kgN_ha = ...
    A.modeled_2017_NO3_leak_kgNha;
row.observed_2017_NO3_leaching_kgN_ha = ...
    A.observed_2017_NO3_leak_kgNha;
row.NO3_2017_relative_error = (A.modeled_2017_NO3_leak_kgNha- ...
    A.observed_2017_NO3_leak_kgNha)/A.observed_2017_NO3_leak_kgNha;
row.final_NO3_pool_kgN_ha = 10*no3_pool(end);
row.maximum_NO3_mobile_mgN_L = A.maximum_NO3_mobile_mgNL;
row.cumulative_drainage_mm = sum(drainage);
row.mean_daily_drainage_mm = mean(drainage);
row.mean_WFPS_pct = mean(wfps);
row.minimum_WFPS_pct = min(wfps);
row.maximum_WFPS_pct = max(wfps);
row.N_closure_residual_gN_m2 = A.N_closure_residual_gNm2;
row.net_mineralization_gN_m2 = A.mineralization_gNm2;
row.net_immobilization_gN_m2 = A.net_immobilization_gNm2;
row.net_immobilization_days = A.net_immobilization_days;
row.row1_placeholder_plant_N_gN_m2 = A.row1_placeholder_plant_N_gNm2;
row.maximum_N_partition_residual_gN_m2_d = ...
    max(abs(partition_residual));
row.ERW_steps_per_day = charge.erw_steps_per_day;
row.max_exchange_sum_residual = charge.max_exchange_sum_residual;
row.max_charge_residual_umolc_L = charge.max_charge_residual_umolc_L;
row.max_Al_residual_umol_m2 = charge.max_Al_residual_umol_m2;
row.carbon_days_2017 = annual.n_days(1);
row.carbon_days_2018 = annual.n_days(2);
row.carbon_days_2019 = annual.n_days(3);
row.carbon_days_2020 = annual.n_days(4);
row.GPP_2017_gC_m2 = annual.GPP(1);
row.GPP_2018_gC_m2 = annual.GPP(2);
row.GPP_2019_gC_m2 = annual.GPP(3);
row.GPP_2020_Jan01_Dec30_gC_m2 = annual.GPP(4);
row.Re_2017_gC_m2 = annual.Re(1);
row.Re_2018_gC_m2 = annual.Re(2);
row.Re_2019_gC_m2 = annual.Re(3);
row.Re_2020_Jan01_Dec30_gC_m2 = annual.Re(4);
row.NEE_2017_gC_m2 = annual.NEE(1);
row.NEE_2018_gC_m2 = annual.NEE(2);
row.NEE_2019_gC_m2 = annual.NEE(3);
row.NEE_2020_Jan01_Dec30_gC_m2 = annual.NEE(4);
row.observed_GPP_2017_gC_m2 = carbon_skill.target_GPP(1);
row.observed_GPP_2018_gC_m2 = carbon_skill.target_GPP(2);
row.observed_GPP_2019_gC_m2 = carbon_skill.target_GPP(3);
row.observed_GPP_2020_Jan01_Dec30_gC_m2 = carbon_skill.target_GPP(4);
row.observed_Re_2017_gC_m2 = carbon_skill.target_Re(1);
row.observed_Re_2018_gC_m2 = carbon_skill.target_Re(2);
row.observed_Re_2019_gC_m2 = carbon_skill.target_Re(3);
row.observed_Re_2020_Jan01_Dec30_gC_m2 = carbon_skill.target_Re(4);
row.observed_NEE_2017_gC_m2 = carbon_skill.target_NEE(1);
row.observed_NEE_2018_gC_m2 = carbon_skill.target_NEE(2);
row.observed_NEE_2019_gC_m2 = carbon_skill.target_NEE(3);
row.observed_NEE_2020_Jan01_Dec30_gC_m2 = carbon_skill.target_NEE(4);
row.carbon_2020_window = "2020-01-01_to_2020-12-30";
row.max_NEE_identity_residual_gC_m2_d = carbon.max_identity_residual;
row.carbon_GPP_NRMSE = carbon_skill.GPP_NRMSE;
row.carbon_NEE_NRMSE = carbon_skill.NEE_NRMSE;
row.carbon_Re_NRMSE_diagnostic = carbon_skill.Re_NRMSE;
row.carbon_GPP_normalized_SSE = carbon_skill.GPP_SSE;
row.carbon_NEE_normalized_SSE = carbon_skill.NEE_SSE;
row.carbon_Re_normalized_SSE_diagnostic = carbon_skill.Re_SSE;
row.carbon_GPP_NEE_normalized_SSE = carbon_skill.independent_SSE;
row.carbon_comparison_count = carbon_skill.count;
row.carbon_GPP_NEE_comparison_count = carbon_skill.independent_count;
end

function C = summarize_combinations(T)
ids = unique(T.combination_id,'stable');
rows = struct([]);
for k = 1:numel(ids)
    X = T(T.combination_id == ids(k),:);
    [~,order] = sort(X.run_id);
    X = X(order,:);
    if height(X) ~= 3 || ~isequal(X.run_id.',10:12)
        error('Bouchet_Joint_Sweep:CombinationRuns', ...
            'Combination %d does not contain exactly runs 10 to 12.',ids(k));
    end
    r.combination_id = ids(k);
    r.combination_tag = X.combination_tag(1);
    r.campaign = unique_text(X.campaign,'campaign',ids(k));
    r.exchange_scheme = unique_text( ...
        X.exchange_scheme,'exchange_scheme',ids(k));
    r.CEC_model = unique_text(X.CEC_model,'CEC_model',ids(k));
    r.matlab_release = unique_text( ...
        X.matlab_release,'matlab_release',ids(k));
    r.numerical_identity = unique_text( ...
        X.numerical_identity,'numerical_identity',ids(k));
    r.sweep_layout = unique_text(X.sweep_layout,'sweep_layout',ids(k));
    r.AppeloAlpha = unique_scalar(X.AppeloAlpha,'AppeloAlpha',ids(k));
    r.DissolutionFactor = unique_scalar( ...
        X.DissolutionFactor,'DissolutionFactor',ids(k));
    r.AppeloLogKH_Na = unique_scalar( ...
        X.AppeloLogKH_Na,'AppeloLogKH_Na',ids(k));
    r.AppeloLogKH_K = unique_scalar( ...
        X.AppeloLogKH_K,'AppeloLogKH_K',ids(k));
    r.AppeloLogKH_Ca = unique_scalar( ...
        X.AppeloLogKH_Ca,'AppeloLogKH_Ca',ids(k));
    r.AppeloLogKH_Mg = unique_scalar( ...
        X.AppeloLogKH_Mg,'AppeloLogKH_Mg',ids(k));
    r.AppeloLogKH_Al = unique_scalar( ...
        X.AppeloLogKH_Al,'AppeloLogKH_Al',ids(k));
    r.CEC_cmolc_kg = unique_scalar( ...
        X.CEC_cmolc_kg,'CEC_cmolc_kg',ids(k));
    r.CECpHSensitivity = unique_scalar( ...
        X.CECpHSensitivity,'CECpHSensitivity',ids(k));
    r.validation_fingerprint = unique_text( ...
        X.validation_fingerprint,'validation_fingerprint',ids(k));
    r.observation_fingerprint = unique_text( ...
        X.observation_fingerprint,'observation_fingerprint',ids(k));
    r.KNmaxd_per_day = unique_scalar( ...
        X.KNmaxd_per_day,'KNmaxd',ids(k));
    r.Kbot_mm_h = unique_scalar_nan(X.Kbot_mm_h,'Kbot',ids(k));
    r.Kbot_boundary = X.Kbot_boundary(1);
    r.pH_control_RMSE = X.pH_observation_RMSE(1);
    r.pH_basalt_RMSE = X.pH_observation_RMSE(2);
    r.pH_RMSE = combined_rmse( ...
        X.pH_observation_SSE(1:2),X.pH_observation_count(1:2));
    r.pH_season_control_RMSE = X.pH_season_RMSE(1);
    r.pH_season_basalt_RMSE = X.pH_season_RMSE(2);
    r.pH_season_RMSE = combined_rmse( ...
        X.pH_season_SSE(1:2),X.pH_season_count(1:2));
    r.pH_envelope_miss_RMS = combined_rmse( ...
        X.pH_envelope_miss_SSE(1:2),X.pH_season_count(1:2));
    model_difference = [ ...
        X.bulk_pH_season_pre(2)-X.bulk_pH_season_pre(1), ...
        X.bulk_pH_season_2017(2)-X.bulk_pH_season_2017(1), ...
        X.bulk_pH_season_2018(2)-X.bulk_pH_season_2018(1), ...
        X.bulk_pH_season_2019(2)-X.bulk_pH_season_2019(1), ...
        X.bulk_pH_season_2020(2)-X.bulk_pH_season_2020(1)];
    observed_difference = [ ...
        X.observed_bulk_pH_pre(2)-X.observed_bulk_pH_pre(1), ...
        X.observed_bulk_pH_2017(2)-X.observed_bulk_pH_2017(1), ...
        X.observed_bulk_pH_2018(2)-X.observed_bulk_pH_2018(1), ...
        X.observed_bulk_pH_2019(2)-X.observed_bulk_pH_2019(1), ...
        X.observed_bulk_pH_2020(2)-X.observed_bulk_pH_2020(1)];
    r.pH_difference_pre = model_difference(1);
    r.pH_difference_2017 = model_difference(2);
    r.pH_difference_2018 = model_difference(3);
    r.pH_difference_2019 = model_difference(4);
    r.pH_difference_2020 = model_difference(5);
    r.pH_difference_RMSE = sqrt(mean( ...
        (model_difference-observed_difference).^2));
    r.pH_season_0_10_control_RMSE = X.pH_season_0_10_RMSE(1);
    r.pH_season_0_10_basalt_RMSE = X.pH_season_0_10_RMSE(2);
    r.pH_season_0_10_RMSE = combined_rmse( ...
        X.pH_season_0_10_SSE(1:2),X.pH_season_count(1:2));
    r.pH_season_10_30_control_RMSE = X.pH_season_10_30_RMSE(1);
    r.pH_season_10_30_basalt_RMSE = X.pH_season_10_30_RMSE(2);
    r.pH_season_10_30_RMSE = combined_rmse( ...
        X.pH_season_10_30_SSE(1:2),X.pH_season_count(1:2));
    r.pH_season_discrete_control_RMSE = X.pH_season_discrete_RMSE(1);
    r.pH_season_discrete_basalt_RMSE = X.pH_season_discrete_RMSE(2);
    r.pH_season_discrete_RMSE = combined_rmse( ...
        X.pH_season_discrete_SSE(1:2),X.pH_season_discrete_count(1:2));
    observed_difference_0_10 = [ ...
        X.observed_bulk_pH_season_0_10_pre(2)-X.observed_bulk_pH_season_0_10_pre(1), ...
        X.observed_bulk_pH_season_0_10_2017(2)-X.observed_bulk_pH_season_0_10_2017(1), ...
        X.observed_bulk_pH_season_0_10_2018(2)-X.observed_bulk_pH_season_0_10_2018(1), ...
        X.observed_bulk_pH_season_0_10_2019(2)-X.observed_bulk_pH_season_0_10_2019(1), ...
        X.observed_bulk_pH_season_0_10_2020(2)-X.observed_bulk_pH_season_0_10_2020(1)];
    observed_difference_10_30 = [ ...
        X.observed_bulk_pH_season_10_30_pre(2)-X.observed_bulk_pH_season_10_30_pre(1), ...
        X.observed_bulk_pH_season_10_30_2017(2)-X.observed_bulk_pH_season_10_30_2017(1), ...
        X.observed_bulk_pH_season_10_30_2018(2)-X.observed_bulk_pH_season_10_30_2018(1), ...
        X.observed_bulk_pH_season_10_30_2019(2)-X.observed_bulk_pH_season_10_30_2019(1), ...
        X.observed_bulk_pH_season_10_30_2020(2)-X.observed_bulk_pH_season_10_30_2020(1)];
    r.pH_difference_0_10_RMSE = sqrt(mean( ...
        (model_difference-observed_difference_0_10).^2));
    r.pH_difference_10_30_RMSE = sqrt(mean( ...
        (model_difference-observed_difference_10_30).^2));
    r.pH_difference_discrete_RMSE = sqrt(mean([ ...
        model_difference-observed_difference_0_10, ...
        model_difference-observed_difference_10_30].^2));
    r.pH_difference_2020_model_over_observed_0_10 = ...
        model_difference(5)/observed_difference_0_10(5);
    r.pH_envelope_miss_0_10_RMS = combined_rmse( ...
        X.pH_envelope_miss_0_10_SSE(1:2),X.pH_season_count(1:2));
    r.pH_envelope_miss_10_30_RMS = combined_rmse( ...
        X.pH_envelope_miss_10_30_SSE(1:2),X.pH_season_count(1:2));
    r.pH_envelope_miss_discrete_RMS = combined_rmse( ...
        X.pH_envelope_miss_discrete_SSE(1:2),X.pH_season_discrete_count(1:2));
    r.base_sat_season_control_RMSE = X.base_sat_season_discrete_RMSE(1);
    r.base_sat_season_basalt_RMSE = X.base_sat_season_discrete_RMSE(2);
    r.base_sat_season_discrete_RMSE = combined_rmse( ...
        X.base_sat_season_discrete_SSE(1:2),X.base_sat_season_discrete_count(1:2));
    period = {'pre','2017','2018','2019','2020'};
    bs_model = zeros(1,5); bs_obs_0_10 = zeros(1,5); bs_obs_10_30 = zeros(1,5);
    for y = 1:5
        bs_model(y) = diff(X.(['base_sat_season_' period{y}])(1:2));
        bs_obs_0_10(y) = diff(X.(['observed_base_sat_0_10_' period{y}])(1:2));
        bs_obs_10_30(y) = diff(X.(['observed_base_sat_10_30_' period{y}])(1:2));
    end
    r.base_sat_difference_2020 = bs_model(5);
    r.base_sat_difference_discrete_RMSE = sqrt(mean([ ...
        bs_model-bs_obs_0_10, bs_model-bs_obs_10_30].^2));
    r.N2O_control_RMSE_gN_ha_d = X.N2O_observation_RMSE_gN_ha_d(1);
    r.N2O_basalt_RMSE_gN_ha_d = X.N2O_observation_RMSE_gN_ha_d(2);
    r.N2O_RMSE_gN_ha_d = combined_rmse( ...
        X.N2O_observation_SSE(1:2),X.N2O_observation_count(1:2));
    r.N2O_NRMSE = sqrt(sum(X.N2O_observation_SSE(1:2))/ ...
        sum(X.N2O_observation_target_sumsq(1:2)));
    r.N2O_cumulative_control_2017_window = X.N2O_cumulative_2017_window(1);
    r.N2O_cumulative_control_2018_window = X.N2O_cumulative_2018_window(1);
    r.N2O_cumulative_basalt_2017_window = X.N2O_cumulative_2017_window(2);
    r.N2O_cumulative_basalt_2018_window = X.N2O_cumulative_2018_window(2);
    r.N2O_cumulative_control_2017_observed_kgN_ha = ...
        X.N2O_cumulative_2017_observed_kgN_ha(1);
    r.N2O_cumulative_control_2018_observed_kgN_ha = ...
        X.N2O_cumulative_2018_observed_kgN_ha(1);
    r.N2O_cumulative_basalt_2017_observed_kgN_ha = ...
        X.N2O_cumulative_2017_observed_kgN_ha(2);
    r.N2O_cumulative_basalt_2018_observed_kgN_ha = ...
        X.N2O_cumulative_2018_observed_kgN_ha(2);
    r.N2O_cumulative_control_2017_kgN_ha = ...
        X.N2O_cumulative_2017_modeled_kgN_ha(1);
    r.N2O_cumulative_control_2018_kgN_ha = ...
        X.N2O_cumulative_2018_modeled_kgN_ha(1);
    r.N2O_cumulative_basalt_2017_kgN_ha = ...
        X.N2O_cumulative_2017_modeled_kgN_ha(2);
    r.N2O_cumulative_basalt_2018_kgN_ha = ...
        X.N2O_cumulative_2018_modeled_kgN_ha(2);
    r.N2O_cumulative_control_2017_ratio = ...
        exp(X.N2O_cumulative_2017_log_ratio(1));
    r.N2O_cumulative_control_2018_ratio = ...
        exp(X.N2O_cumulative_2018_log_ratio(1));
    r.N2O_cumulative_basalt_2017_ratio = ...
        exp(X.N2O_cumulative_2017_log_ratio(2));
    r.N2O_cumulative_basalt_2018_ratio = ...
        exp(X.N2O_cumulative_2018_log_ratio(2));
    r.N2O_cumulative_log_ratio_RMS = combined_rmse( ...
        X.N2O_cumulative_log_ratio_SS(1:2),X.N2O_cumulative_count(1:2));
    r.control_NO3_2017_kgN_ha = ...
        X.modeled_2017_NO3_leaching_kgN_ha(1);
    r.basalt_NO3_2017_kgN_ha = ...
        X.modeled_2017_NO3_leaching_kgN_ha(2);
    r.NO3_2017_relative_RMSE = sqrt(mean( ...
        X.NO3_2017_relative_error(1:2).^2));
    r.carbon_control_GPP_NRMSE = X.carbon_GPP_NRMSE(1);
    r.carbon_basalt_GPP_NRMSE = X.carbon_GPP_NRMSE(2);
    r.carbon_control_NEE_NRMSE = X.carbon_NEE_NRMSE(1);
    r.carbon_basalt_NEE_NRMSE = X.carbon_NEE_NRMSE(2);
    r.carbon_Re_NRMSE_diagnostic = sqrt(sum( ...
        X.carbon_Re_normalized_SSE_diagnostic(1:2))/sum( ...
        X.carbon_comparison_count(1:2)));
    r.carbon_GPP_NEE_NRMSE = sqrt(sum( ...
        X.carbon_GPP_NEE_normalized_SSE(1:2))/sum( ...
        X.carbon_GPP_NEE_comparison_count(1:2)));
    r.control_drainage_mm = X.cumulative_drainage_mm(1);
    r.basalt_drainage_mm = X.cumulative_drainage_mm(2);
    r.calcite_drainage_mm = X.cumulative_drainage_mm(3);
    r.maximum_abs_N_closure_gN_m2 = max(abs(X.N_closure_residual_gN_m2));
    r.maximum_abs_charge_residual_umolc_L = ...
        max(X.max_charge_residual_umolc_L);
    r.maximum_abs_NEE_identity_residual_gC_m2_d = ...
        max(X.max_NEE_identity_residual_gC_m2_d);
    r.physical_gates_passed = true;
    if isempty(rows)
        rows = repmat(r,numel(ids),1);
    else
        rows(k) = r;
    end
end
C = struct2table(rows);
end

function value = combined_rmse(sse,count)
if any(~isfinite(sse)) || any(~isfinite(count)) || sum(count) <= 0
    value = NaN;
else
    value = sqrt(sum(sse)/sum(count));
end
end

function value = unique_scalar(x,name,id)
u = unique(x);
if numel(u) ~= 1 || ~isfinite(u)
    error('Bouchet_Joint_Sweep:CombinationParameter', ...
        'Combination %d has inconsistent %s values.',id,name);
end
value = u;
end

function value = unique_scalar_nan(x,name,id)
if ~all(arrayfun(@(v) isequaln(v,x(1)),x))
    error('Bouchet_Joint_Sweep:CombinationParameter', ...
        'Combination %d has inconsistent %s values.',id,name);
end
value = x(1);
end

function value = unique_text(x,name,id)
u = unique(string(x));
if numel(u) ~= 1 || ismissing(u)
    error('Bouchet_Joint_Sweep:CombinationParameter', ...
        'Combination %d has inconsistent %s values.',id,name);
end
value = u;
end

function mask = pareto_minimize(objectives)
n = size(objectives,1);
mask = true(n,1);
for i = 1:n
    for j = 1:n
        if i ~= j && all(objectives(j,:) <= objectives(i,:)) && ...
                any(objectives(j,:) < objectives(i,:))
            mask(i) = false;
            break
        end
    end
end
end

function fp = observation_identity(files)
fp = ERW_Checkpoint_Fingerprint(0,1,'joint_sweep_observations',files{:});
end

function skill = residual_skill(residual,target)
residual = double(residual(:));
target = double(target(:));
if isempty(residual) || numel(residual) ~= numel(target) || ...
        any(~isfinite([residual;target]))
    error('Bouchet_Joint_Sweep:InvalidSkillInput', ...
        'An observational comparison is empty, mismatched or nonfinite.');
end
skill.sse = sum(residual.^2);
skill.target_sumsq = sum(target.^2);
skill.count = numel(residual);
skill.rmse = sqrt(skill.sse/skill.count);
skill.bias = mean(residual);
skill.nrmse = sqrt(skill.sse/skill.target_sumsq);
skill.model_mean = NaN;
skill.observed_mean = NaN;
end

function skill = empty_skill
skill = struct('sse',NaN,'target_sumsq',NaN,'count',0, ...
    'rmse',NaN,'bias',NaN,'nrmse',NaN,'target',nan(1,5));
end

function daily = empty_daily(nd,ncase,cases)
daily.date_datenum = nan(nd,1);
daily.case_id = [cases.id].';
daily.case_tag = string({cases.tag}).';
fields = {'bulk_pH','pore_pH','N2O_gN_m2_d', ...
    'NO3_leaching_gN_m2_d','NO3_pool_gN_m2','drainage_mm_d', ...
    'WFPS_pct','GPP_gC_m2_d','Re_gC_m2_d','NEE_gC_m2_d'};
for k = 1:numel(fields)
    daily.(fields{k}) = nan(nd,ncase);
end
end

function atomic_writetable(T,file)
folder = fileparts(file);
tmp = [tempname(folder) '.csv'];
cleanup = onCleanup(@() delete_if_present(tmp));
writetable(T,tmp);
movefile(tmp,file,'f');
end

function atomic_save_pre_control(file,T,inventory,score,dates,bulk_pH, ...
        pore_pH,observations,observation_fingerprint)
folder = fileparts(file);
tmp = [tempname(folder) '.mat'];
cleanup = onCleanup(@() delete_if_present(tmp));
save(tmp,'T','inventory','score','dates','bulk_pH','pore_pH', ...
    'observations','observation_fingerprint','-v7.3');
movefile(tmp,file,'f');
end

function atomic_save(file,T,combination_table,inventory,daily, ...
        observations,observation_fingerprint,selection_notice)
folder = fileparts(file);
tmp = [tempname(folder) '.mat'];
cleanup = onCleanup(@() delete_if_present(tmp));
save(tmp,'T','combination_table','inventory','daily','observations', ...
    'observation_fingerprint','selection_notice','-v7.3');
movefile(tmp,file,'f');
end

function atomic_write_text(file,value)
folder = fileparts(file);
tmp = [tempname(folder) '.txt'];
cleanup = onCleanup(@() delete_if_present(tmp));
fid = fopen(tmp,'wt');
if fid < 0
    error('Bouchet_Joint_Sweep:StatusWrite', ...
        'Cannot open temporary summary status file.');
end
file_cleanup = onCleanup(@() fclose(fid));
fprintf(fid,'%s',value);
clear file_cleanup
movefile(tmp,file,'f');
end

function delete_if_present(file)
if isfile(file); delete(file); end
end
