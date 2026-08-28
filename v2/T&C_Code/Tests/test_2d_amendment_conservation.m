function test_2d_amendment_conservation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The 2-D invariants must still hold while feedstock is being applied.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Everything else that checks the lateral scheme checks it in the quiet
%% period. test_lateral_routing drives it with synthetic states.
%% test_2d_driver_equivalence runs the control case for ten days from
%% 1 January, and the amendment does not arrive until day 314. So the code
%% is proven conservative and bounded only when nothing is dissolving.
%%
%% The amendment is when the scheme is most stressed: a 5000 g/m2 pulse
%% releases Ca, Mg, Na, Al, Si and DIC into solution, moves pH, and drives
%% the largest inter-cell gradients of the run. If lateral transport is
%% going to lose mass or produce an unphysical state, that is when.
%%
%% The application dates are fixed at 10 November, so starting on the 5th
%% puts the amendment at day 5 and the whole check costs a few minutes
%% rather than the four hours a run from January would take.

fprintf('\n=== test_2d_amendment_conservation ===\n');
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

result_dir = tempname; mkdir(result_dir);
old_ckpt = getenv('TC_CKPT_MINUTES');
old_workers = getenv('TC_PARALLEL_WORKERS');
cleanup = onCleanup(@() restore(result_dir,old_ckpt,old_workers)); %#ok<NASGU>
setenv('TC_CKPT_MINUTES','0');
setenv('TC_PARALLEL_WORKERS','');

%% Run 2 is basalt on the whole-field layout: rock on every cell, so the
%% amendment perturbs all twelve columns rather than only the margin.
result_file = Run_One_Subfield(2,'StartDate',datenum(2016,11,5), ...
    'NDays',8,'GridNx',2,'GridNy',2,'AppeloAlpha',1.3, ...
    'ParallelWorkers',6,'SaveResults',true,'ResultDirectory',result_dir);
S = load(result_file,'Domain','ERW','O','V','PH_d','P');

%% ---------------------------------------------------------------------
%% The amendment actually happened inside the window
%% ---------------------------------------------------------------------
add = S.ERW{1}.TC_const.M_rock_add;
applied = sum(add);
n_events = sum(add > 0);
check('exactly one application inside the window',n_events == 1, ...
    sprintf('got %d events',n_events));
check('the application delivers 5000 g/m2',abs(applied-5000) < 1e-9, ...
    sprintf('got %.6g',applied));
expected_row = 5*144 + 1;
check('the application lands on the 10 November row', ...
    ~isempty(find(add > 0,1)) && find(add > 0,1) == expected_row, ...
    sprintf('got row %d, expected %d',find(add > 0,1),expected_row));

%% Every column carries primary mineral once the rock is on.
has_mineral = arrayfun(@(c) any(S.ERW{c}.bg.M_min(end,:) > 0), ...
    (1:numel(S.ERW)).');
check('every column holds primary mineral after application', ...
    all(has_mineral),sprintf('%d of %d columns', ...
    sum(has_mineral),numel(has_mineral)));

%% ---------------------------------------------------------------------
%% Lateral conservation across the amendment
%% ---------------------------------------------------------------------
L = S.Domain.lateral;
species = {'Na','Al','An','DIC'};
for s = 1:min(numel(species),size(L.erw_mass_residual_umol,2))
    r = sum(L.erw_mass_residual_umol(:,s),'omitnan');
    moved = sum(abs(L.erw_mass_buf_to_main_umol(:, ...
        min(s,size(L.erw_mass_buf_to_main_umol,2)))),'omitnan');
    rel = abs(r)/max(1,moved);
    check(sprintf('lateral %s conserved across the amendment',species{s}), ...
        isfinite(rel) && rel < 1e-6, ...
        sprintf('residual %.6g of %.6g moved (relative %.3e)',r,moved,rel));
end
heat_rel = max(abs(L.heat_residual_relative),[],'omitnan');
check('lateral energy conserved across the amendment', ...
    isfinite(heat_rel) && heat_rel < 1e-6, ...
    sprintf('max relative %.3e',heat_rel));

%% Adaptive sub-stepping must not have saturated: a capped solver limits
%% fluxes instead of resolving them, which would hide a real gradient.
check('water sub-stepping did not run away', ...
    isfinite(max(L.sub_steps_max)) && max(L.sub_steps_max) < 1e4, ...
    sprintf('peak %g',max(L.sub_steps_max)));
check('heat sub-stepping did not run away', ...
    isfinite(max(L.heat_sub_steps_max)) && max(L.heat_sub_steps_max) < 1e4, ...
    sprintf('peak %g',max(L.heat_sub_steps_max)));

%% ---------------------------------------------------------------------
%% The state stays physical while the pulse dissolves
%% ---------------------------------------------------------------------
check('soil moisture is finite and nonnegative', ...
    all(isfinite(S.O(:))) && all(S.O(:) >= 0), 'O left the physical range');
check('soil water is finite and nonnegative', ...
    all(isfinite(S.V(:))) && all(S.V(:) >= 0), 'V left the physical range');
ph = S.PH_d(:); ph = ph(~isnan(ph));
check('pH stays inside 0 to 14', ...
    ~isempty(ph) && all(ph > 0 & ph < 14), ...
    sprintf('range %.4g to %.4g',min(ph),max(ph)));
check('biogeochemical pools are finite', ...
    all(isfinite(S.P(~isnan(S.P)))), 'P contains Inf');

for c = 1:numel(S.ERW)
    bg = S.ERW{c}.bg;
    bad = false;
    for name = {'Ca','Mg','Na','K','Al','Si','DIC','NO3','NH4'}
        v = bg.(name{1});
        if any(~isfinite(v(:))) || any(v(:) < 0); bad = true; end
    end
    if bad
        check(sprintf('column %d chemistry is finite and nonnegative',c), ...
            false,'negative or non-finite species');
        break
    end
end
if ~bad
    check('all columns keep finite nonnegative chemistry',true);
end

fprintf('\n  %d passed, %d failed\n',n_pass,n_fail);
if n_fail > 0
    error('test_2d_amendment_conservation:Failed', ...
        '%d amendment-period check(s) failed',n_fail);
end
fprintf(['PASS: the 2-D invariants hold while feedstock is applied.' ...
    '\n']);
end

function restore(result_dir,old_ckpt,old_workers)
setenv('TC_CKPT_MINUTES',old_ckpt);
setenv('TC_PARALLEL_WORKERS',old_workers);
pool = gcp('nocreate');
if ~isempty(pool); delete(pool); end
if isfolder(result_dir); rmdir(result_dir,'s'); end
end
