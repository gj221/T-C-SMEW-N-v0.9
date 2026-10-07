function compare_u_v(dir_u,dir_v)
%% Every saved variable of a run of package u against the same run of
%% package v (compare_u_v.sbatch). v keeps the particle-class histories
%% (psd, d, lamb, SSA, delta_d) at its initial rows and at the last row of
%% every ERW day, listed in weath.particle_rows; they must equal u's full
%% arrays at those rows. When both results keep only those rows (v against
%% w), the rows and their values must be identical. Every other model value
%% must be identical; the fsolve options object is compared by class and
%% settings, since two loaded copies of one optimoptions object never
%% compare equal (as in test_parallel_erw_equivalence). Variables holding
%% run times, paths, provenance or code hashes differ by design and are
%% listed apart. It only reads.
%%
%%   compare_u_v('<scratch>/u','<scratch>/v')    folders or result files
fu = one_result(dir_u); fv = one_result(dir_v);
fprintf('u: %s\nv: %s\n',fu,fv);
U = load(fu); V = load(fv);
particle = {'psd','d','lamb','SSA','delta_d'};
names = union(fieldnames(U),fieldnames(V));
one_side = {}; numeric_bad = {}; timing = {}; other = {};
for k = 1:numel(names)
    n = names{k};
    if ~isfield(U,n) || ~isfield(V,n)
        one_side{end+1} = n; %#ok<AGROW>
        continue
    end
    if strcmp(n,'ERW') || isequaln(U.(n),V.(n)); continue; end
    if (isnumeric(U.(n)) || islogical(U.(n))) && ...
            isempty(regexp(n,'(?i)(time|toc)','once'))
        numeric_bad{end+1} = n; %#ok<AGROW>
    elseif isnumeric(U.(n))
        timing{end+1} = n; %#ok<AGROW>
    else
        other{end+1} = sprintf('%s [%s]',n,what_differs(U.(n),V.(n))); %#ok<AGROW>
    end
end
fprintf('In one file only: %s\n',list(one_side));
fprintf('Run times (differ by design): %s\n',list(timing));
fprintf('Text, structs and cells that differ (paths, provenance and code hashes are expected): %s\n',list(other));
fprintf('Numeric model variables that differ (must be none): %s\n',list(numeric_bad));

erw_ok = true;
if numel(U.ERW) ~= numel(V.ERW)
    fprintf('ERW: %d columns in u, %d in v\n',numel(U.ERW),numel(V.ERW));
    erw_ok = false;
end
for c = 1:min(numel(U.ERW),numel(V.ERW))
    Eu = U.ERW{c}; Ev = V.ERW{c};
    bad = {};
    for f = union(fieldnames(Eu),fieldnames(Ev))'
        fn = f{1};
        if ~isfield(Eu,fn) || ~isfield(Ev,fn)
            bad{end+1} = [fn ' (one side only)']; %#ok<AGROW>
        elseif strcmp(fn,'weath')
            for g = union(fieldnames(Eu.weath),fieldnames(Ev.weath))'
                gn = g{1};
                if any(strcmp(gn,[particle {'particle_rows'}])); continue; end
                if ~isfield(Eu.weath,gn) || ~isfield(Ev.weath,gn) || ...
                        ~isequaln(Eu.weath.(gn),Ev.weath.(gn))
                    bad{end+1} = ['weath.' gn]; %#ok<AGROW>
                end
            end
            thin_u = isfield(Eu.weath,'particle_rows');
            if thin_u && isfield(Ev.weath,'particle_rows')
                if ~isequaln(Eu.weath.particle_rows,Ev.weath.particle_rows)
                    bad{end+1} = 'weath.particle_rows'; %#ok<AGROW>
                end
                for p = particle
                    if ~isfield(Ev.weath,p{1}) || ~isfield(Eu.weath,p{1}) || ...
                            ~isequaln(Ev.weath.(p{1}),Eu.weath.(p{1}))
                        bad{end+1} = ['weath.' p{1}]; %#ok<AGROW>
                    end
                end
                fprintf('ERW column %d: %d particle rows kept in both\n', ...
                    c,numel(Ev.weath.particle_rows));
            elseif thin_u
                bad{end+1} = 'weath.particle_rows (u only)'; %#ok<AGROW>
            elseif isfield(Ev.weath,'particle_rows')
                rows = Ev.weath.particle_rows;
                for p = particle
                    if ~isfield(Ev.weath,p{1}) || ~isfield(Eu.weath,p{1}) || ...
                            ~isequaln(Ev.weath.(p{1}),Eu.weath.(p{1})(rows,:,:))
                        bad{end+1} = ['weath.' p{1} ' at the kept rows']; %#ok<AGROW>
                    end
                end
                fprintf('ERW column %d: %d particle rows kept of %d\n', ...
                    c,numel(rows),size(Eu.weath.psd,1));
            else
                for p = particle
                    if isfield(Eu.weath,p{1}) && (~isfield(Ev.weath,p{1}) || ...
                            ~isequaln(Eu.weath.(p{1}),Ev.weath.(p{1})))
                        bad{end+1} = ['weath.' p{1}]; %#ok<AGROW>
                    end
                end
                fprintf('ERW column %d: no particle histories kept as blocks\n',c);
            end
        elseif strcmp(fn,'var') && isstruct(Eu.var) && isstruct(Ev.var) && ...
                isfield(Eu.var,'fsolve_options') && isfield(Ev.var,'fsolve_options')
            if ~isequaln(rmfield(Eu.var,'fsolve_options'),rmfield(Ev.var,'fsolve_options')) || ...
                    ~same_solver_options(Eu.var.fsolve_options,Ev.var.fsolve_options)
                bad{end+1} = 'var'; %#ok<AGROW>
            end
        elseif ~isequaln(Eu.(fn),Ev.(fn))
            bad{end+1} = fn; %#ok<AGROW>
        end
    end
    if isempty(bad)
        fprintf('ERW column %d: identical\n',c);
    else
        fprintf('ERW column %d differs in: %s\n',c,strjoin(bad,', '));
        erw_ok = false;
    end
end
if erw_ok && isempty(numeric_bad) && isempty(one_side)
    fprintf('PASS: v reproduces u; only the particle histories are thinned.\n');
else
    fprintf('CHECK: differences above.\n');
end
end

function f = one_result(d)
if isfile(d); f = d; return; end
L = dir(fullfile(d,'TC_2D_run*.mat'));
if numel(L) ~= 1
    error('compare_u_v:Result','Expected one result in %s, found %d.',d,numel(L));
end
f = fullfile(L.folder,L.name);
end

function s = what_differs(a,b)
if isstruct(a) && isstruct(b) && isscalar(a) && isscalar(b)
    f = union(fieldnames(a),fieldnames(b));
    d = f(cellfun(@(n) ~isfield(a,n) || ~isfield(b,n) || ~isequaln(a.(n),b.(n)),f));
    s = ['fields ' strjoin(d',' ')];
else
    s = class(a);
end
end

function s = list(c)
if isempty(c); s = 'none'; else; s = strjoin(c,', '); end
end

function same = same_solver_options(a,b)
same = strcmp(class(a),class(b));
for prop = {'Algorithm','Display','FunctionTolerance','StepTolerance', ...
        'MaxIterations','MaxFunctionEvaluations'}
    if same && isprop(a,prop{1}) && isprop(b,prop{1})
        same = isequaln(a.(prop{1}),b.(prop{1}));
    end
end
end
