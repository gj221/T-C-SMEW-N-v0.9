function [bad,report] = Daily_Store_Check(X,Y,ignore)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The saved outputs of a 2-D run of package y (Y) against the same run
%% kept with full records (X: package x, or package y with
%% TC_FULL_RECORDS=1), both as loaded from their result files.
%%
%% A daily-store result keeps daily stores, hourly field means and the
%% final state where X keeps every hourly and ten-minute record, so each
%% of its outputs is checked against the same reduction of X's records,
%% made with the code the 2-D frame uses (Field_Daily_From_Hourly,
%% Field_Hour_Values with Field_Hourly_Means, ERW_Daily_Record), and must
%% be identical (isequaln):
%%   - every variable both results hold (the daily stores, Domain, the
%%     lateral ledgers, the carried per-square state, the run's inputs),
%%     except those named in ignore (identities, timing, paths);
%%   - the hourly-derived daily stores O_d, ZWT_d, Lk_d, f_d, q_runon_d,
%%     Rh_d, Rd_d and the hours after the last daily step (Hourly_tail);
%%   - the hourly field means (Field_hourly);
%%   - each square's ERW day summaries (ERW_daily) and its t = 0 row;
%%   - the final state: the last hours and days of every hourly and daily
%%     array (Final_state), and each square's ERW state, its last block
%%     against the same rows of X.
%% A full-record Y holds X's arrays themselves, compared whole.
%%
%%   [bad,report] = Daily_Store_Check(X,Y,ignore)
%%
%% bad lists the differences (empty when Y reproduces X); report holds
%% the counts and the variables found on one side only.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3; ignore = {}; end
if ~isfield(Y,'FULL_RECORDS')
    error('Daily_Store_Check:NotY','Y is not a result of package y (no FULL_RECORDS).');
end
if isfield(X,'FULL_RECORDS') && ~X.FULL_RECORDS
    error('Daily_Store_Check:NotFull','X keeps daily stores; it must hold the full records.');
end
bad = {};
report = struct();

%% ------------------------------------------------------------------------
%% Variables both results hold
%% ------------------------------------------------------------------------
names = setdiff(intersect(fieldnames(X),fieldnames(Y)),{'ERW'});
n_same = 0; ignored = {};
for k = 1:numel(names)
    n = names{k};
    if same_value(X.(n),Y.(n))
        n_same = n_same + 1;
    elseif any(strcmp(n,ignore))
        ignored{end+1} = n; %#ok<AGROW>
    else
        bad{end+1} = sprintf('%s [%s]',n,what_differs(X.(n),Y.(n))); %#ok<AGROW>
    end
end
x_only = setdiff(fieldnames(X),fieldnames(Y));
y_only = setdiff(fieldnames(Y),fieldnames(X));
if ~Y.FULL_RECORDS
    %% Their last rows are in Final_state (checked below).
    x_only = setdiff(x_only,[Y.HOURLY_NAMES(:); Y.DAILY_STATE_ONLY(:)]);
    missing = setdiff(Y.DAILY_STORE_NAMES,fieldnames(Y));
    if ~isempty(missing)
        bad{end+1} = ['daily stores missing: ' strjoin(missing,', ')];
    end
end
report.n_both = numel(names);
report.n_same = n_same;
report.ignored_differing = ignored;
report.x_only = x_only;
report.y_only = y_only;

%% ------------------------------------------------------------------------
%% Daily rows of the hourly records
%% ------------------------------------------------------------------------
NN = X.NN; dth = X.dth; N_col = X.N_col; nl = Y.n_bio_layers;
%% The daily steps: the hours at 01:00 from hour 2, where the loop starts.
day_i = find(X.Datam(2:NN,4) == 1) + 1;
stores = {'O_d','O','mean'; 'ZWT_d','ZWT','mean'; 'Lk_d','Lk','rate_total'; ...
    'f_d','f','rate_total'; 'q_runon_d','q_runon','rate_total'; ...
    'Rh_d','Rh','total'; 'Rd_d','Rd','total'};
if isempty(day_i) || ~isequal(Y.Hourly_tail.hours,[day_i(end) NN])
    bad{end+1} = 'Hourly_tail.hours';
    tail_hours = [];
else
    tail_hours = day_i(end):NN;
end
for s = 1:size(stores,1)
    src = X.(stores{s,2});
    if strcmp(stores{s,2},'O'); src = src(:,1:nl,:); end
    expect = nan(size(Y.(stores{s,1})));
    for d = 1:numel(day_i)
        hours = max(1,day_i(d)-24):day_i(d)-1;
        expect(d+1,:,:) = Field_Daily_From_Hourly(src(hours,:,:),stores{s,3},dth);
    end
    if ~isequaln(expect,Y.(stores{s,1}))
        bad{end+1} = sprintf('%s (from the hourly %s)',stores{s,1},stores{s,2}); %#ok<AGROW>
    end
    if ~isempty(tail_hours) && ~isequaln(Field_Daily_From_Hourly( ...
            src(tail_hours,:,:),stores{s,3},dth),Y.Hourly_tail.(stores{s,2}))
        bad{end+1} = sprintf('Hourly_tail.%s',stores{s,2}); %#ok<AGROW>
    end
end
report.n_daily_steps = numel(day_i);

%% ------------------------------------------------------------------------
%% Hourly field means
%% ------------------------------------------------------------------------
%% The crown cover in force at hour i is the initial one before the first
%% daily step, then the row the latest daily step wrote; X records it per
%% day only with a fixed vegetation cover (OPT_VCA = 0).
names_fh = Y.FIELD_HOURLY_NAMES;
crown_known = ~isfield(X,'OPT_VCA') || X.OPT_VCA == 0;
expect = struct();
for k = 1:numel(names_fh)
    expect.(names_fh{k}) = nan(size(Y.Field_hourly.(names_fh{k})));
end
for i = 1:NN
    jd = sum(day_i <= i) + 1;
    if jd >= 2
        Ccrown_hour = reshape(X.Ccrown_t(jd,:,:),X.cc,N_col).';
    else
        Ccrown_hour = X.CcrownFIX;
    end
    fh = Field_Hour_Values(i,N_col,X.cc,nl,Ccrown_hour, ...
        X.Rn,X.H,X.QE,X.G,X.Ts,X.EG,X.ELitter,X.ESN,X.ESN_In,X.EWAT,X.EICE, ...
        X.EIn_urb,X.EIn_rock,X.Lk,X.Rh,X.Rd,X.T_H,X.T_L,X.EIn_H,X.EIn_L, ...
        X.An_H,X.An_L,X.Rdark_H,X.Rdark_L,X.O,X.Tdp);
    fhm = Field_Hourly_Means(fh,Y.Field_hourly.weights);
    for k = 1:numel(names_fh)
        expect.(names_fh{k})(i,:) = fhm.(names_fh{k});
    end
end
report.field_hourly_unchecked = {};
for k = 1:numel(names_fh)
    if strcmp(names_fh{k},'GPP') && ~crown_known
        report.field_hourly_unchecked = {'GPP'};
        continue
    end
    if ~isequaln(expect.(names_fh{k}),Y.Field_hourly.(names_fh{k}))
        bad{end+1} = ['Field_hourly.' names_fh{k}]; %#ok<AGROW>
    end
end
report.n_hours = NN;

%% ------------------------------------------------------------------------
%% Final state of a daily-store result
%% ------------------------------------------------------------------------
if ~Y.FULL_RECORDS
    F = Y.Final_state;
    for k = 1:numel(Y.HOURLY_NAMES)
        n = Y.HOURLY_NAMES{k};
        if ~isequaln(F.hourly.(n),X.(n)(F.hour_first:NN,:,:))
            bad{end+1} = ['Final_state.hourly.' n]; %#ok<AGROW>
        end
    end
    for k = 1:numel(Y.DAILY_NAMES)
        n = Y.DAILY_NAMES{k};
        if ~isequaln(F.daily.(n),X.(n)(F.daily_first:X.j,:,:,:))
            bad{end+1} = ['Final_state.daily.' n]; %#ok<AGROW>
        end
    end
    report.final_hours = [F.hour_first NN];
    report.final_days = [F.daily_first X.j];
end

%% ------------------------------------------------------------------------
%% ERW, square by square
%% ------------------------------------------------------------------------
report.n_erw_squares = 0;
report.n_erw_days = 0;
if isfield(X,'ERW') || isfield(Y,'ERW')
    if ~isfield(X,'ERW') || ~isfield(Y,'ERW') || numel(X.ERW) ~= numel(Y.ERW)
        bad{end+1} = 'ERW: the squares differ in number';
    else
        for c = 1:numel(X.ERW)
            [bad_c,n_c] = erw_square(X.ERW{c},Y.ERW{c},Y.ERW_daily{c},Y.FULL_RECORDS,c);
            bad = [bad bad_c]; %#ok<AGROW>
            report.n_erw_days = report.n_erw_days + n_c;
        end
        report.n_erw_squares = numel(X.ERW);
    end
end
end

function [bad,n_days] = erw_square(Ex,Ey,D,full_records,c)
bad = {};
tag = sprintf('ERW{%d}',c);
n_days = 0;
if ~isfield(Ex,'bg')
    %% A square the run never reached with ERW: both must be untouched.
    if ~isequaln(Ex,Ey); bad{end+1} = [tag ' (no ERW state in X)']; end
    return
end
%% The day summaries, against the same summary of X's ten-minute rows.
for k = 1:D.n
    rows = D.rows(k,1):D.rows(k,2);
    rock = [sum(Ex.TC_const.M_rock_add(rows)) sum(Ex.TC_const.M_rock_add_CaCO3(rows))];
    r = ERW_Daily_Record('summarise',Ex.bg,Ex.weath,rows,Ex.dtw,rock);
    for kind = {'end','mean','total','max'}
        for q = fieldnames(r.(kind{1}))'
            if ~isequaln(r.(kind{1}).(q{1}),D.(kind{1}).(q{1})(k,:))
                bad{end+1} = sprintf('%s ERW day %d %s.%s',tag,k,kind{1},q{1}); %#ok<AGROW>
            end
        end
    end
end
n_days = D.n;
%% The days tile the grid from row 2 to the last row, one after another.
if D.n == 0 || D.rows(1,1) ~= 2 || D.rows(end,2) ~= Ex.n_t || ...
        any(D.rows(2:end,1) ~= D.rows(1:end-1,2)+1)
    bad{end+1} = [tag ' ERW days do not tile the grid'];
end
%% The t = 0 row, as the run finally holds it.
particle = {'psd','d','lamb','SSA','delta_d','particle_rows'};
for part = {'bg','weath'}
    R0 = Ey.daily_initial.(part{1});
    for f = fieldnames(R0)'
        if any(strcmp(f{1},particle)); continue; end
        if ~isequaln(R0.(f{1}),Ex.(part{1}).(f{1})(1,:,:))
            bad{end+1} = sprintf('%s t = 0 row %s.%s',tag,part{1},f{1}); %#ok<AGROW>
        end
    end
end

%% The state the square ends with. Package x's grids are rebuilt from
%% their definitions in y and not kept; y's own fields hold its records.
dropped = {'t','th','td','CO2_h'};
y_own = {'full_records','rock','daily_initial','daily_pending','daily_row', ...
    'bg_fill','weath_fill','T_ref_row'};
if full_records
    first = 1; last = Ex.n_t;
else
    first = D.rows(end,1) - 1; last = D.rows(end,2);
    dropped{end+1} = 'tm';
end
for f = union(fieldnames(Ex),fieldnames(Ey))'
    fn = f{1};
    if any(strcmp(fn,[dropped y_own])); continue; end
    if ~isfield(Ex,fn) || ~isfield(Ey,fn)
        bad{end+1} = sprintf('%s.%s (one side only)',tag,fn); %#ok<AGROW>
        continue
    end
    switch fn
        case {'bg','weath'}
            bad = [bad block_fields(Ex.(fn),Ey.(fn),first,last,full_records, ...
                [tag '.' fn],particle)]; %#ok<AGROW>
        case 'plant'
            %% Never written by the run: a daily-block run keeps row 1.
            for q = fieldnames(Ex.plant)'
                xq = Ex.plant.(q{1});
                if ~full_records; xq = xq(1,:,:); end
                if ~isfield(Ey.plant,q{1}) || ~isequaln(xq,Ey.plant.(q{1}))
                    bad{end+1} = sprintf('%s.plant.%s',tag,q{1}); %#ok<AGROW>
                end
            end
        case 'const'
            if ~isequaln(Ex.const.T_ref(first:last,:),Ey.const.T_ref)
                bad{end+1} = [tag '.const.T_ref']; %#ok<AGROW>
            end
            cx = rmfield(Ex.const,'T_ref'); cy = rmfield(Ey.const,'T_ref');
            if ~same_value(cx,cy)
                bad{end+1} = [tag '.const [' what_differs(cx,cy) ']']; %#ok<AGROW>
            end
        case 'TC'
            %% x kept the hourly water and CO2 series of the whole run here.
            tx = rmfield(Ex.TC,intersect(fieldnames(Ex.TC),{'Tw','Iw','Lw','CO2_atm'}));
            if ~same_value(tx,Ey.TC)
                bad{end+1} = [tag '.TC [' what_differs(tx,Ey.TC) ']']; %#ok<AGROW>
            end
        case 'TC_const'
            tx = Ex.TC_const; ty = Ey.TC_const;
            if ~full_records
                %% The rock schedule is kept as its rows (E.rock).
                tx = rmfield(tx,{'M_rock_add','M_rock_add_CaCO3'});
                ty = rmfield(ty,{'M_rock_add','M_rock_add_CaCO3'});
                add = Ex.TC_const.M_rock_add; add_c = Ex.TC_const.M_rock_add_CaCO3;
                rows = find(add ~= 0 | add_c ~= 0);
                if ~isequal(Ey.rock.rows,rows) || ~isequaln(Ey.rock.add,add(rows)) || ...
                        ~isequaln(Ey.rock.add_CaCO3,add_c(rows)) || ...
                        ~isequal(Ey.rock.positive_rows,find(add > 0))
                    bad{end+1} = [tag '.rock']; %#ok<AGROW>
                end
            end
            if ~same_value(tx,ty)
                bad{end+1} = [tag '.TC_const [' what_differs(tx,ty) ']']; %#ok<AGROW>
            end
        case 'var'
            if ~same_value(rmfield(Ex.var,'fsolve_options'),rmfield(Ey.var,'fsolve_options')) || ...
                    ~same_solver_options(Ex.var.fsolve_options,Ey.var.fsolve_options)
                bad{end+1} = [tag '.var']; %#ok<AGROW>
            end
        case 'dm_end'
            if Ey.dm_end + first - 1 ~= Ex.dm_end
                bad{end+1} = [tag '.dm_end']; %#ok<AGROW>
            end
        otherwise
            if ~same_value(Ex.(fn),Ey.(fn))
                bad{end+1} = sprintf('%s.%s [%s]',tag,fn,what_differs(Ex.(fn),Ey.(fn))); %#ok<AGROW>
            end
    end
end
end

function bad = block_fields(Sx,Sy,first,last,full_records,tag,particle)
%% Rows first:last of X against Y's block (the whole arrays with full
%% records). Particle-class histories are kept as rows on the grid
%% (particle_rows): each Y row must equal X's row of the same grid index.
bad = {};
by_rows = isfield(Sx,'particle_rows') || isfield(Sy,'particle_rows');
for f = union(fieldnames(Sx),fieldnames(Sy))'
    fn = f{1};
    if by_rows && any(strcmp(fn,particle)); continue; end
    if strcmp(fn,'CO2_atm') && ~full_records; continue; end    %% forcing, read from W
    if ~isfield(Sx,fn) || ~isfield(Sy,fn)
        bad{end+1} = sprintf('%s.%s (one side only)',tag,fn); %#ok<AGROW>
        continue
    end
    xv = Sx.(fn);
    if ~full_records
        xv = xv(first:last,:,:);
    end
    if ~isequaln(xv,Sy.(fn))
        bad{end+1} = sprintf('%s.%s [%s]',tag,fn,what_differs(xv,Sy.(fn))); %#ok<AGROW>
    end
end
if by_rows
    if ~isfield(Sx,'particle_rows') || ~isfield(Sy,'particle_rows')
        bad{end+1} = [tag '.particle_rows (one side only)'];
        return
    end
    [found,where] = ismember(Sy.particle_rows,Sx.particle_rows);
    if ~all(found) || (full_records && ~isequal(Sy.particle_rows,Sx.particle_rows))
        bad{end+1} = [tag '.particle_rows'];
        return
    end
    for p = particle(1:end-1)
        if isfield(Sx,p{1}) ~= isfield(Sy,p{1})
            bad{end+1} = [tag '.' p{1} ' (one side only)']; %#ok<AGROW>
        elseif isfield(Sy,p{1}) && ~isequaln(Sy.(p{1}),Sx.(p{1})(where,:,:))
            bad{end+1} = [tag '.' p{1}]; %#ok<AGROW>
        end
    end
end
end

function s = what_differs(a,b)
if isstruct(a) && isstruct(b) && isscalar(a) && isscalar(b)
    f = union(fieldnames(a),fieldnames(b));
    d = f(cellfun(@(n) ~isfield(a,n) || ~isfield(b,n) || ~same_value(a.(n),b.(n)),f));
    s = ['fields ' strjoin(d',' ')];
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

function tf = same_value(a,b)
%% isequaln, except that function handles compare by their text and the
%% values they captured: two loaded copies of one anonymous function are
%% never isequaln.
if isequaln(a,b); tf = true; return; end
tf = false;
if isa(a,'function_handle') && isa(b,'function_handle')
    if ~strcmp(func2str(a),func2str(b)); return; end
    fa = functions(a); fb = functions(b);
    if isfield(fa,'workspace') ~= isfield(fb,'workspace'); return; end
    tf = ~isfield(fa,'workspace') || same_value(fa.workspace,fb.workspace);
    return
end
if ~strcmp(class(a),class(b)) || ~isequal(size(a),size(b)); return; end
if isstruct(a)
    names = fieldnames(a);
    if ~isempty(setxor(names,fieldnames(b))); return; end
    for k = 1:numel(a)
        for q = 1:numel(names)
            if ~same_value(a(k).(names{q}),b(k).(names{q})); return; end
        end
    end
    tf = true;
elseif iscell(a)
    for k = 1:numel(a)
        if ~same_value(a{k},b{k}); return; end
    end
    tf = true;
end
end

function same = same_solver_options(a,b)
%% Solver option objects are compared by their settings.
same = strcmp(class(a),class(b));
for prop = {'Algorithm','Display','FunctionTolerance','StepTolerance', ...
        'OptimalityTolerance','MaxIterations','MaxFunctionEvaluations'}
    if same && isprop(a,prop{1}) && isprop(b,prop{1})
        same = isequaln(a.(prop{1}),b.(prop{1}));
    end
end
end
