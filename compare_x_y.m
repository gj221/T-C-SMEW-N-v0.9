function ok = compare_x_y(path_x,path_y,x_root,y_root)
%% A run of package y against the same run of package x (compare_x_y.sbatch).
%%
%% Package y computes every number as x does and keeps less: a multi-square
%% y result holds daily stores, hourly field means and the final state
%% where x holds every hourly and ten-minute record. Daily_Store_Check
%% (v2/T&C_Code) checks each y output against the same reduction of x's
%% records and every variable both hold, all with isequaln; a full-record
%% y result (TC_FULL_RECORDS=1) is compared with x's arrays whole. Run
%% times, the results' locations, provenance and code identities differ by
%% design and are listed apart; paths into the two trees compare by what
%% they point to (x's root read as y's). It only reads.
%%
%%   ok = compare_x_y('<scratch>/x','<scratch>/y')    folders or result files
%%   ok = compare_x_y(...,x_root,y_root)              the roots, if not the
%%                                                    ones the runs recorded
%%
%% Package y's v2/T&C_Code must be on the path.

if isempty(which('Daily_Store_Check'))
    error('compare_x_y:Path','Put package y''s v2/T&C_Code on the path first.');
end
fx = one_result(path_x); fy = one_result(path_y);
fprintf('x: %s\ny: %s\n',fx,fy);
X = load(fx); Y = load(fy);
if isfield(X,'FULL_RECORDS')
    error('compare_x_y:NotX','%s is a package-y result, not package x''s.',fx);
end
if ~isfield(Y,'FULL_RECORDS')
    error('compare_x_y:NotY','%s is not a package-y result (no FULL_RECORDS).',fy);
end
fprintf('y keeps %s records.\n',ternary(Y.FULL_RECORDS,'full','daily-store'));

%% What differs between the two trees by design: run times, where each
%% result went (the options and the result and checkpoint paths), the
%% provenance (code digests, times), the code identities the checkpoint
%% fingerprints carry, and the loop counter the two drivers leave last.
by_design = {'run_provenance','Computational_Time','ERW_ckpt_last_toc', ...
    'ERW_ckpt_fingerprint','ERW_ckpt_runtime_identity','ERW_execution_identity', ...
    'ERW_ckpt_file','ERW_ckpt_tmp','result_file','tmp_result','ckpt_dir', ...
    'fl_save','opt','varargin','k'};
%% The roots: v2 is current_directory in Run_One_Subfield.
if nargin < 4
    x_root = ''; y_root = '';
    if isfield(X,'current_directory') && isfield(Y,'current_directory')
        x_root = fileparts(X.current_directory);
        y_root = fileparts(Y.current_directory);
    end
end
rename = {};
if ~isempty(x_root) && ~isempty(y_root)
    rename = {char(x_root),char(y_root)};
    fprintf('Paths compared with %s read as %s\n',rename{1},rename{2});
end
[bad,r] = Daily_Store_Check(X,Y,by_design,rename);

fprintf('Variables both hold: %d, identical: %d\n',r.n_both,r.n_same);
fprintf('Differ by design (times, paths, provenance, code identity): %s\n', ...
    list(r.ignored_differing));
fprintf('In x only: %s\n',list(r.x_only));
fprintf('In y only: %s\n',list(r.y_only));
fprintf('Hourly-derived daily stores: 7 series over %d daily steps and the tail\n', ...
    r.n_daily_steps);
fprintf('Hourly field means: %d hours%s\n',r.n_hours, ...
    ternary(isempty(r.field_hourly_unchecked),'', ...
    ' (GPP not checked: x keeps no crown cover per hour)'));
if isfield(r,'final_hours')
    fprintf('Final state: hours %d-%d, daily rows %d-%d\n',r.final_hours,r.final_days);
end
fprintf('ERW: %d squares, %d square-days of summaries\n',r.n_erw_squares,r.n_erw_days);

ok = isempty(bad);
if ok
    fprintf('PASS: y reproduces x; every daily output and the final state are identical.\n');
else
    fprintf('CHECK: %d difference(s):\n',numel(bad));
    fprintf('  %s\n',bad{:});
end
end

function f = one_result(d)
if isfile(d); f = d; return; end
L = dir(fullfile(d,'TC_2D_run*.mat'));
if numel(L) ~= 1
    error('compare_x_y:Result','Expected one result in %s, found %d.',d,numel(L));
end
f = fullfile(L.folder,L.name);
end

function s = list(c)
if isempty(c); s = 'none'; else; s = strjoin(c,', '); end
end

function out = ternary(cond,a,b)
if cond; out = a; else; out = b; end
end
