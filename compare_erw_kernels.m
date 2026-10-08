function ok = compare_erw_kernels(path_a,path_b,mode_a,mode_b)
%% Two runs of one package that differ only in TC_ERW_KERNEL
%% (compare_erw_kernels.sbatch): every saved variable must be identical,
%% down to every field of every struct and cell (a daily-store result keeps
%% most of its numbers in ERW_daily, Final_state and Field_hourly). Only the
%% run times and the paths each run wrote to differ by design, and are
%% listed apart; run_provenance must agree but for the kernel it names and
%% the result directory. Function handles compare by their text and the
%% values they captured, and solver option objects by their public
%% properties, since two loaded copies of either never compare isequaln.
%% It only reads.
%%
%%   ok = compare_erw_kernels('<scratch>/r2_matlab','<scratch>/r2_fast')
%%   ok = compare_erw_kernels(a,b,'matlab','fast')   also checks that each
%%                                                   provenance names its kernel
%% Folders (one TC_2D_run*.mat each) or result files.

if nargin < 3; mode_a = ''; end
if nargin < 4; mode_b = ''; end
fa = one_result(path_a); fb = one_result(path_b);
fprintf('a: %s\nb: %s\n',fa,fb);
A = load(fa); B = load(fb);

%% The roots: v2 is current_directory in Run_One_Subfield, so paths into
%% two copies of the package compare by what they point to.
rename = {};
if isfield(A,'current_directory') && isfield(B,'current_directory') && ...
        ~strcmp(A.current_directory,B.current_directory)
    rename = {fileparts(A.current_directory),fileparts(B.current_directory)};
    fprintf('Paths compared with %s read as %s\n',rename{1},rename{2});
end

by_design = {'Computational_Time','ERW_ckpt_last_toc','ERW_ckpt_file', ...
    'ERW_ckpt_tmp','opt','varargin','fl_save','result_file','tmp_result', ...
    'ckpt_dir'};
names = union(fieldnames(A),fieldnames(B));
bad = {}; differing_by_design = {}; n_same = 0;
for k = 1:numel(names)
    n = names{k};
    if ~isfield(A,n) || ~isfield(B,n)
        bad{end+1} = sprintf('%s (in one result only)',n); %#ok<AGROW>
    elseif strcmp(n,'run_provenance')
        continue
    elseif same_value(A.(n),B.(n),rename)
        n_same = n_same + 1;
    elseif any(strcmp(n,by_design))
        differing_by_design{end+1} = n; %#ok<AGROW>
    else
        bad{end+1} = sprintf('%s [%s]',n,what_differs(A.(n),B.(n),rename)); %#ok<AGROW>
    end
end

%% The provenance: the same run but for its kernel and where it went.
kernels = {'',''};
if isfield(A,'run_provenance') && isfield(B,'run_provenance')
    pa = A.run_provenance; pb = B.run_provenance;
    if isfield(pa,'erw_kernel'); kernels{1} = char(pa.erw_kernel); end
    if isfield(pb,'erw_kernel'); kernels{2} = char(pb.erw_kernel); end
    own = {'erw_kernel','result_directory'};
    qa = rmfield(pa,intersect(fieldnames(pa),own));
    qb = rmfield(pb,intersect(fieldnames(pb),own));
    if ~same_value(qa,qb,rename)
        bad{end+1} = sprintf('run_provenance [%s]',what_differs(qa,qb,rename));
    end
end
expected = {mode_a,mode_b};
for s = 1:2
    if ~isempty(expected{s}) && ~strcmp(kernels{s},expected{s})
        bad{end+1} = sprintf('run_provenance.erw_kernel of %s is ''%s'', not ''%s''', ...
            char('a'+s-1),kernels{s},expected{s}); %#ok<AGROW>
    end
end

fprintf('Kernels named by the provenance: a %s, b %s\n', ...
    quoted(kernels{1}),quoted(kernels{2}));
fprintf('Variables: %d, identical: %d\n',numel(names),n_same);
fprintf('Differ by design (run times, result paths): %s\n',list(differing_by_design));
if isfield(B,'FULL_RECORDS')
    fprintf('Records: %s\n',ternary(B.FULL_RECORDS,'full', ...
        'daily stores (ERW_daily, Final_state, Field_hourly compared field by field)'));
end
if isfield(A,'ERW') && iscell(A.ERW)
    fprintf('ERW squares compared: %d\n',numel(A.ERW));
end
ok = isempty(bad);
if ok
    fprintf('PASS: every saved value is identical under both kernels.\n');
else
    fprintf('CHECK: %d difference(s):\n',numel(bad));
    fprintf('  %s\n',bad{:});
end
end

function f = one_result(d)
if isfile(d); f = d; return; end
L = dir(fullfile(d,'TC_2D_run*.mat'));
if numel(L) ~= 1
    error('compare_erw_kernels:Result','Expected one result in %s, found %d.',d,numel(L));
end
f = fullfile(L.folder,L.name);
end

function tf = same_value(a,b,rename)
%% isequaln, except that function handles compare by their text and the
%% values they captured, objects by their public properties, and text
%% after rename{1} becomes rename{2}.
if isequaln(a,b); tf = true; return; end
tf = false;
if ~isempty(rename) && (ischar(a) || isstring(a)) && strcmp(class(a),class(b))
    tf = isequaln(strrep(a,rename{1},rename{2}),b);
    return
end
if isa(a,'function_handle') && isa(b,'function_handle')
    if ~strcmp(func2str(a),func2str(b)); return; end
    fa = functions(a); fb = functions(b);
    if isfield(fa,'workspace') ~= isfield(fb,'workspace'); return; end
    tf = ~isfield(fa,'workspace') || same_value(fa.workspace,fb.workspace,rename);
    return
end
if ~strcmp(class(a),class(b)) || ~isequal(size(a),size(b)); return; end
if isstruct(a)
    names = fieldnames(a);
    if ~isempty(setxor(names,fieldnames(b))); return; end
    for k = 1:numel(a)
        for q = 1:numel(names)
            if ~same_value(a(k).(names{q}),b(k).(names{q}),rename); return; end
        end
    end
    tf = true;
elseif iscell(a)
    for k = 1:numel(a)
        if ~same_value(a{k},b{k},rename); return; end
    end
    tf = true;
elseif isobject(a) && isscalar(a)
    props = properties(a);
    for q = 1:numel(props)
        if ~same_value(a.(props{q}),b.(props{q}),rename); return; end
    end
    tf = true;
end
end

function s = what_differs(a,b,rename)
if isstruct(a) && isstruct(b) && isscalar(a) && isscalar(b)
    f = union(fieldnames(a),fieldnames(b));
    d = f(cellfun(@(n) ~isfield(a,n) || ~isfield(b,n) || ...
        ~same_value(a.(n),b.(n),rename),f));
    s = ['fields ' strjoin(d',' ')];
elseif iscell(a) && iscell(b) && isequal(size(a),size(b))
    d = find(~cellfun(@(x,y) same_value(x,y,rename),a,b));
    s = sprintf('cells %s',mat2str(d(:).'));
elseif isnumeric(a) && isnumeric(b) && isequal(size(a),size(b))
    dd = abs(double(a(:))-double(b(:)));
    dd(isnan(dd) & ~(isnan(a(:)) & isnan(b(:)))) = Inf;
    dd(isnan(dd)) = 0;
    [m,at] = max(dd);
    s = sprintf('max |d| %.3g at element %d of %d',m,at,numel(dd));
else
    s = sprintf('%s %s against %s %s',class(a),mat2str(size(a)),class(b),mat2str(size(b)));
end
end

function s = list(c)
if isempty(c); s = 'none'; else; s = strjoin(c,', '); end
end

function s = quoted(t)
if isempty(t); s = '(none recorded)'; else; s = ['''' t '''']; end
end

function out = ternary(cond,a,b)
if cond; out = a; else; out = b; end
end
