function test_erw_lateral_cation_daily_delta(n_days)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% test_erw_lateral_cation_daily_delta
%%
%% ERW_Daily_Step_2D runs once per day (at hour 1), but the 2-D driver's
%% lateral coupling (Biogeo_Lateral_Pools) runs every hour and keeps
%% writing into P(j,60/61/68/69/52/53/76,:) for all 24 hours of that day.
%% Before the fix this file guards against, the daily Ca/Mg/K/Si delta fed
%% to ERW was read as P_c(j,:) - P_c(j-1,:): P_c(j-1,:) is exactly what
%% BIOGEO_UNIT copied into P_c(j,:) before hour 1's lateral step, so it
%% cancels out of that subtraction and hours 2-24 of lateral transport
%% -- already folded into that same P_c(j-1,:) snapshot from the previous
%% day -- vanish from the signal. NO3/NH4/PO4 avoid this by re-anchoring
%% their ERW-side totals to P directly, every day; Ca/Mg/K/Si must do the
%% same. This test fails loudly on the pre-fix code, where ERW's internal
%% Ca_tot/Mg_tot/K_tot/Si_tot drift away from T&C's own P-pool totals in
%% any run with active field-margin lateral transport, and passes once
%% they are anchored the same way PO4/NH4/NO3 already are.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 1 || isempty(n_days); n_days = 6; end
fprintf('\n=== test_erw_lateral_cation_daily_delta (%d days) ===\n',n_days);
test_dir = fileparts(mfilename('fullpath'));
code_dir = fileparts(test_dir);
v2 = fileparts(code_dir);
addpath(fullfile(v2,'run'),code_dir, ...
    fullfile(code_dir,'Inputs'),fullfile(code_dir,'Parameters'));

n_pass = 0; n_fail = 0;
    function check(name,cond,detail)
        if nargin < 3; detail = ''; end
        if cond
            n_pass = n_pass + 1; fprintf('  PASS  %s\n',name);
        else
            n_fail = n_fail + 1; fprintf('  FAIL  %s  %s\n',name,detail);
        end
    end

%% ---------------------------------------------------------------------
%% The daily step must anchor Ca/Mg/K/Si_tot the same way it already
%% anchors PO4/NH4/NO3_tot, not just carry a delta forward from P_c(j-1,:).
%% ---------------------------------------------------------------------
daily_code = fileread(fullfile(code_dir,'ERW_Daily_Step_2D.m'));
check('the daily step anchors PO4/NH4/NO3_tot to P (reference pattern)', ...
    contains(daily_code,'bg.PO4_tot(anchor_row,:)') && ...
    contains(daily_code,'bg.NH4_tot(anchor_row,:)') && ...
    contains(daily_code,'bg.NO3_tot(anchor_row,:)'));
check('the daily step anchors Ca/Mg/K/Si_tot to P the same way', ...
    contains(daily_code,'bg.Ca_tot(anchor_row,:)') && ...
    contains(daily_code,'bg.Mg_tot(anchor_row,:)') && ...
    contains(daily_code,'bg.K_tot(anchor_row,:)') && ...
    contains(daily_code,'bg.Si_tot(anchor_row,:)'));

%% ---------------------------------------------------------------------
%% End-to-end: a genuine two-column run with active field-margin lateral
%% transport (main vs margin receive different rock treatment, which
%% produces a real Ca/Mg/K/Si concentration gradient between cells).
%% ---------------------------------------------------------------------
result_dir = tempname; mkdir(result_dir);
old_ckpt = getenv('TC_CKPT_MINUTES');
old_workers = getenv('TC_PARALLEL_WORKERS');
cleanup = onCleanup(@() restore(result_dir,old_ckpt,old_workers)); %#ok<NASGU>
setenv('TC_CKPT_MINUTES','0');
%% These checks read the run's full hourly and ten-minute records, which
%% a multi-square run keeps only when asked (package y).
old_records = getenv('TC_FULL_RECORDS');
records_cleanup = onCleanup(@() setenv('TC_FULL_RECORDS',old_records)); %#ok<NASGU>
setenv('TC_FULL_RECORDS','1');
setenv('TC_PARALLEL_WORKERS','');

result_file = Run_One_Subfield(7,'StartDate',datenum(2016,1,1), ...
    'NDays',n_days,'GridNx',2,'GridNy',2,'AppeloAlpha',1.3, ...
    'ParallelWorkers',6,'SaveResults',true,'ResultDirectory',result_dir);
S = load(result_file,'ERW','P','Datam','Domain');

check('the run has more than one lateral cell (N_col > 1)', ...
    numel(S.ERW) > 1,sprintf('%d column(s)',numel(S.ERW)));
lateral_ca_moved = sum(abs(S.Domain.lateral.mass_buf_to_main(:, ...
    strcmp(S.Domain.lateral.species_names,'Ca'))),'omitnan');
check('the run actually moves Ca laterally between cells', ...
    isfinite(lateral_ca_moved) && lateral_ca_moved > 0, ...
    sprintf('domain lateral Ca throughput %.6g',lateral_ca_moved));

daily_i = find(S.Datam(:,4) == 1);
worst_daily = struct('Ca',0,'Mg',0,'K',0,'Si',0);
worst_final = struct('Ca',0,'Mg',0,'K',0,'Si',0);
missing = {};
for c = 1:numel(S.ERW)
    E = S.ERW{c};
    if ~isfield(E,'dm_end') || ~isfield(E,'bg') || ...
            ~all(isfield(E.bg,{'Ca_tot','Mg_tot','K_tot','Si_tot'}))
        missing{end+1} = sprintf('column %d cation handoff',c); %#ok<AGROW>
        continue
    end
    expect_final = struct( ...
        'Ca',(S.P(end,60,c)+S.P(end,61,c))/40  *1e6, ...
        'Mg',(S.P(end,68,c)+S.P(end,69,c))/24.3*1e6, ...
        'K', (S.P(end,52,c)+S.P(end,53,c))/39.1*1e6, ...
        'Si', S.P(end,76,c)/28.1*1e6);
    actual_final = struct('Ca',E.bg.Ca_tot(E.dm_end),'Mg',E.bg.Mg_tot(E.dm_end), ...
        'K',E.bg.K_tot(E.dm_end),'Si',E.bg.Si_tot(E.dm_end));
    for f = {'Ca','Mg','K','Si'}
        name = f{1};
        rel = abs(actual_final.(name)-expect_final.(name))/ ...
            max(1,abs(expect_final.(name)));
        worst_final.(name) = max(worst_final.(name),rel);
    end
    for d = 1:numel(daily_i)
        boundary_row = find(E.tm < daily_i(d)-1,1,'last');
        if isempty(boundary_row)
            missing{end+1} = sprintf('column %d daily boundary %d',c,d); %#ok<AGROW>
            continue
        end
        expect_daily = struct( ...
            'Ca',(S.P(d+1,60,c)+S.P(d+1,61,c))/40  *1e6, ...
            'Mg',(S.P(d+1,68,c)+S.P(d+1,69,c))/24.3*1e6, ...
            'K', (S.P(d+1,52,c)+S.P(d+1,53,c))/39.1*1e6, ...
            'Si', S.P(d+1,76,c)/28.1*1e6);
        %% The last boundary precedes the final flush, which revisits the
        %% same P row without its delta (ERW_Should_Apply_TC_Daily_Delta):
        %% its anchor is the pre-flush total, while the saved P row holds
        %% the post-flush state. The flush adds the boundary supply, the
        %% leached-salt return of ERW_Leached_Salt_Resupply for Ca and Mg,
        %% for every row after the boundary, and in a control run that is
        %% its only source, so the saved total exceeds the anchor by exactly
        %% that supply. Remove it to compare the anchor with what it read.
        if d == numel(daily_i)
            flush_days = (E.dm_end-boundary_row)*E.dtw;
            expect_daily.Ca = expect_daily.Ca-E.const.I_Ca*flush_days;
            expect_daily.Mg = expect_daily.Mg-E.const.I_Mg*flush_days;
            expect_daily.K  = expect_daily.K -E.const.I_K *flush_days;
            expect_daily.Si = expect_daily.Si-E.const.I_Si*flush_days;
        end
        actual_daily = struct('Ca',E.bg.Ca_tot(boundary_row), ...
            'Mg',E.bg.Mg_tot(boundary_row),'K',E.bg.K_tot(boundary_row), ...
            'Si',E.bg.Si_tot(boundary_row));
        for f = {'Ca','Mg','K','Si'}
            name = f{1};
            rel = abs(actual_daily.(name)-expect_daily.(name))/ ...
                max(1,abs(expect_daily.(name)));
            worst_daily.(name) = max(worst_daily.(name),rel);
        end
    end
end
check('every column exposes the Ca/Mg/K/Si handoff state',isempty(missing), ...
    strjoin(missing,'; '));

for f = {'Ca','Mg','K','Si'}
    name = f{1};
    fprintf('  worst daily %s_tot shadow error %.3e\n',name,worst_daily.(name));
    check(sprintf('%s_tot tracks the full 24-hour P pool at every daily boundary', ...
        name),worst_daily.(name) < 1e-10, ...
        sprintf('relative error %.3e (>= 1e-10 means lateral hours are being lost)', ...
        worst_daily.(name)));
    fprintf('  worst final %s_tot shadow error %.3e\n',name,worst_final.(name));
    check(sprintf('%s_tot tracks the final P pool',name),worst_final.(name) < 1e-10, ...
        sprintf('relative error %.3e',worst_final.(name)));
end

fprintf('\n  %d passed, %d failed\n',n_pass,n_fail);
if n_fail > 0
    error('test_erw_lateral_cation_daily_delta:Failed', ...
        '%d lateral cation delta check(s) failed',n_fail);
end
fprintf(['PASS: ERW Ca/Mg/K/Si inventories reflect all 24 hours of daily ' ...
    'lateral transport, not just hour 1.\n']);
end

function restore(result_dir,old_ckpt,old_workers)
setenv('TC_CKPT_MINUTES',old_ckpt);
setenv('TC_PARALLEL_WORKERS',old_workers);
pool = gcp('nocreate');
if ~isempty(pool); delete(pool); end
if isfolder(result_dir); rmdir(result_dir,'s'); end
end
