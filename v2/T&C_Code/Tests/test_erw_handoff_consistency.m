function test_erw_handoff_consistency(n_days)
%%% ERW dissolved P, NH4 and NO3 track their biogeochemical pools.

if nargin < 1 || isempty(n_days); n_days = 3; end
fprintf('\n=== test_erw_handoff_consistency (%d days) ===\n',n_days);
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

BiogeoPar = Biogeochemistry_Parameter(1);
aP = BiogeoPar.aP;
aNH4 = BiogeoPar.aNH4;
aNO3 = BiogeoPar.aNO3;
fprintf('  handoff fractions: aP %.4g, aNH4 %.4g, aNO3 %.4g\n', ...
    aP,aNH4,aNO3);

step_2d = fileread(fullfile(code_dir,'ERW_Daily_Step_2D.m'));
frame_1d = fileread(fullfile(code_dir,'MAIN_FRAME_BG_ERW_10mins.m'));
check('both drivers anchor P, NH4 and NO3', ...
    contains(step_2d,'bg.PO4_tot(anchor_row,:)') && ...
    contains(step_2d,'bg.NH4_tot(anchor_row,:)') && ...
    contains(step_2d,'bg.NO3_tot(anchor_row,:)') && ...
    contains(frame_1d,'ERW_biogeochem.PO4_tot(anchor_row,:)') && ...
    contains(frame_1d,'ERW_biogeochem.NH4_tot(anchor_row,:)') && ...
    contains(frame_1d,'ERW_biogeochem.NO3_tot(anchor_row,:)'));

result_dir = tempname; mkdir(result_dir);
old_ckpt = getenv('TC_CKPT_MINUTES');
old_workers = getenv('TC_PARALLEL_WORKERS');
cleanup = onCleanup(@() restore(result_dir,old_ckpt,old_workers)); %#ok<NASGU>
setenv('TC_CKPT_MINUTES','0');
setenv('TC_PARALLEL_WORKERS','');

result_file = Run_One_Subfield(7,'StartDate',datenum(2016,1,1), ...
    'NDays',n_days,'GridNx',2,'GridNy',2,'AppeloAlpha',1.3, ...
    'ParallelWorkers',6,'SaveResults',true,'ResultDirectory',result_dir);
S = load(result_file,'ERW','P','Datam');

worst_p = 0; worst_nh4 = 0; worst_no3 = 0;
worst_daily_p = 0; worst_daily_nh4 = 0; worst_daily_no3 = 0;
worst_p_release = 0;
missing = {};
daily_i = find(S.Datam(:,4) == 1);
for c = 1:numel(S.ERW)
    E = S.ERW{c};
    if ~isfield(E,'dm_end') || ~isfield(E,'bg') || ~isfield(E,'weath') || ...
            ~isfield(E.bg,'PO4_tot') || ~isfield(E.bg,'NH4_tot') || ...
            ~isfield(E.bg,'NO3_tot') || ...
            ~isfield(E.weath,'P_weathering_release')
        missing{end+1} = sprintf('column %d handoff',c); %#ok<AGROW>
        continue
    end
    expect_p = aP * S.P(end,43,c)/31 * 1e6;
    expect_nh4 = aNH4 * S.P(end,31,c)/14 * 1e6;
    expect_no3 = aNO3 * S.P(end,32,c)/14 * 1e6;
    actual_p = E.bg.PO4_tot(E.dm_end);
    actual_nh4 = E.bg.NH4_tot(E.dm_end);
    actual_no3 = E.bg.NO3_tot(E.dm_end);
    worst_p = max(worst_p,abs(actual_p-expect_p)/max(1,abs(expect_p)));
    worst_nh4 = max(worst_nh4, ...
        abs(actual_nh4-expect_nh4)/max(1,abs(expect_nh4)));
    worst_no3 = max(worst_no3, ...
        abs(actual_no3-expect_no3)/max(1,abs(expect_no3)));
    for d = 1:numel(daily_i)
        boundary_row = find(E.tm < daily_i(d)-1,1,'last');
        if isempty(boundary_row)
            missing{end+1} = sprintf('column %d daily boundary %d',c,d); %#ok<AGROW>
            continue
        end
        expect_p = aP * S.P(d+1,43,c)/31 * 1e6;
        expect_nh4 = aNH4 * S.P(d+1,31,c)/14 * 1e6;
        expect_no3 = aNO3 * S.P(d+1,32,c)/14 * 1e6;
        worst_daily_p = max(worst_daily_p, ...
            abs(E.bg.PO4_tot(boundary_row)-expect_p)/max(1,abs(expect_p)));
        worst_daily_nh4 = max(worst_daily_nh4, ...
            abs(E.bg.NH4_tot(boundary_row)-expect_nh4)/max(1,abs(expect_nh4)));
        worst_daily_no3 = max(worst_daily_no3, ...
            abs(E.bg.NO3_tot(boundary_row)-expect_no3)/max(1,abs(expect_no3)));
    end
    worst_p_release = max(worst_p_release, ...
        max(abs(E.weath.P_weathering_release(:))));
    if any(E.const.min_st(:,7) ~= 0)
        missing{end+1} = sprintf('column %d P-bearing campaign mineral',c); %#ok<AGROW>
    end
end

check('every column exposes the handoff state',isempty(missing), ...
    strjoin(missing,'; '));
fprintf('  worst final P shadow error %.3e\n',worst_p);
fprintf('  worst final NH4 shadow error %.3e\n',worst_nh4);
fprintf('  worst final NO3 shadow error %.3e\n',worst_no3);
check('PO4_tot tracks aP * P(:,43)',worst_p < 1e-10, ...
    sprintf('relative error %.3e',worst_p));
check('NH4_tot tracks aNH4 * P(:,31)',worst_nh4 < 1e-10, ...
    sprintf('relative error %.3e',worst_nh4));
check('NO3_tot tracks aNO3 * P(:,32)',worst_no3 < 1e-10, ...
    sprintf('relative error %.3e',worst_no3));
fprintf('  worst daily P shadow error %.3e\n',worst_daily_p);
fprintf('  worst daily NH4 shadow error %.3e\n',worst_daily_nh4);
fprintf('  worst daily NO3 shadow error %.3e\n',worst_daily_no3);
check('PO4_tot tracks every normal daily boundary',worst_daily_p < 1e-10, ...
    sprintf('relative error %.3e',worst_daily_p));
check('NH4_tot tracks every normal daily boundary',worst_daily_nh4 < 1e-10, ...
    sprintf('relative error %.3e',worst_daily_nh4));
check('NO3_tot tracks every normal daily boundary',worst_daily_no3 < 1e-10, ...
    sprintf('relative error %.3e',worst_daily_no3));
check('run 7 control mineralogy releases no P',worst_p_release == 0, ...
    sprintf('maximum substep release %.3e umol P m^-2',worst_p_release));

species = {'Ca_tot','Mg_tot','K_tot','Na_tot','An_tot','Al_tot','Si_tot', ...
    'IC_tot','PO4_tot','NO3_tot','NH4_tot'};
bad = {};
for c = 1:numel(S.ERW)
    for s = 1:numel(species)
        if ~isfield(S.ERW{c}.bg,species{s})
            bad{end+1} = sprintf('%s missing in column %d',species{s},c); %#ok<AGROW>
            continue
        end
        v = S.ERW{c}.bg.(species{s});
        if any(~isfinite(v(:))) || any(v(:) < 0)
            bad{end+1} = sprintf('%s col %d',species{s},c); %#ok<AGROW>
        end
    end
end
check('every ERW inventory stays finite and nonnegative',isempty(bad), ...
    strjoin(bad,'; '));

fprintf('\n  %d passed, %d failed\n',n_pass,n_fail);
if n_fail > 0
    error('test_erw_handoff_consistency:Failed', ...
        '%d handoff check(s) failed',n_fail);
end
fprintf('PASS: ERW dissolved P, NH4 and NO3 inventories track their pools.\n');
end

function restore(result_dir,old_ckpt,old_workers)
setenv('TC_CKPT_MINUTES',old_ckpt);
setenv('TC_PARALLEL_WORKERS',old_workers);
pool = gcp('nocreate');
if ~isempty(pool); delete(pool); end
if isfolder(result_dir); rmdir(result_dir,'s'); end
end
