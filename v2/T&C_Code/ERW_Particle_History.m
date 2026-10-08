function varargout = ERW_Particle_History(action,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Particle-class histories of an amended column: the last few ERW days
%% whole, and the last row of every earlier day.
%%
%% psd, d, lamb, SSA and delta_d are rows x size classes x applications:
%% one row per ERW step of the whole run, one slice per scheduled rock
%% application. Over 35 years (1840896 ten-minute rows, 35 applications)
%% each is about 3.6 GB, so the five held whole took 17 GiB of memory in
%% every amended column, and every checkpoint wrote them. ERW_function
%% reads only two things from them: the previous row, and for each
%% application the first row whose last size class is positive (its
%% initial distribution, psd_0(find(psd_0(:,end,k)>0,1),:,k)). Nothing
%% else in the model, the campaign checks or the read-outs reads them.
%%
%% The rows a day writes collect in a block, so ERW_Daily_Step_2D never
%% writes into an array the frame still holds. The blocks of the last
%% keep_days days (three unless 'start' is told otherwise) stay whole.
%% When an older block is dropped, its last written row is kept at its
%% latest value. The daily windows follow one another, so every row the
%% daily step reads comes from today's block or the day before it, exactly
%% as the full arrays gave it; a window repeating the day before would be
%% served exactly too. A read of any other row, which the daily step does
%% not make, stops the run with RowNotKept instead of returning a wrong
%% value.
%%
%% After the run the frame assembles what was kept: weath.psd and the
%% others hold the rows set at initialisation and the last row of every
%% ERW day, at their latest values, and weath.particle_rows gives each
%% one's row on the ten-minute grid.
%%
%% A multi-square run keeps no row of a dropped block (keep_rows false):
%% no reader uses them, the model reads only the previous row and each
%% application's first row, and at ~125 MB per basalt square over 35 years
%% they were most of a field run's particle record. The run then assembles
%% the initial rows and the last row of each block still held, its final
%% state. Its ten-minute state is itself one day's block (ERW_Daily_Step_2D),
%% so it keeps a column with a single application here too, where a
%% full-record run leaves that column's whole arrays in weath.
%% n_rows gives the length of the ten-minute grid when weath was
%% allocated with fewer rows than that (the 2-D daily-block layout).
%%   [H,weath] = ERW_Particle_History('start',weath,weath_names,keep_days,keep_rows,n_rows)
%%   day = ERW_Particle_History('open',H,first_row,n_rows)
%%   v   = ERW_Particle_History('row',H,day,weath,name,r)    row r, 1 x classes x apps
%%   F   = ERW_Particle_History('first',H,day,name)          psd or lamb
%%   day = ERW_Particle_History('write',H,day,out,r)
%%   H   = ERW_Particle_History('track',H,day,out,r,weath)
%%   H   = ERW_Particle_History('close',H,day)
%%   E   = ERW_Particle_History('assemble',E)
%% 'first' returns, for each application, the first row whose last size
%% class is positive (zeros where none is): ERW_function's
%% psd_0(find(psd_0(:,end,k)>0,1),:,k) gives the same row from it as from
%% the whole history, which is that application's initial distribution.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch action
    case 'start'
        [varargout{1},varargout{2}] = start(varargin{:});
    case 'open'
        [H,first_row,n_rows] = varargin{:};
        day = struct('first',first_row,'n',n_rows,'written',false(n_rows,1),'data',struct());
        for k = 1:numel(H.names)
            name = H.names{k};
            day.data.(name) = zeros(n_rows,H.size.(name)(2),H.size.(name)(3), ...
                'like',H.initial.(name));
        end
        varargout{1} = day;
    case 'row'
        varargout{1} = row(varargin{:});
    case 'first'
        [H,~,name] = varargin{:};
        varargout{1} = H.first.(name);
    case 'write'
        [H,day,out,r] = varargin{:};
        idx = r - day.first + 1;
        if idx < 1 || idx > day.n
            error('ERW_Particle_History:RowOutsideDay', ...
                'Row %d is outside the day''s rows %d-%d.',r,day.first,day.first+day.n-1);
        end
        for k = 1:numel(H.names)
            name = H.names{k};
            if isfield(out,name)
                day.data.(name)(idx,:,:) = out.(name);
            end
        end
        day.written(idx) = true;
        varargout{1} = day;
    case 'track'
        varargout{1} = track(varargin{:});
    case 'close'
        [H,day] = varargin{:};
        if any(day.written)
            H = refresh_kept(H,day);
            H.blocks{end+1} = day;
            while numel(H.blocks) > H.keep_days
                H = drop_oldest(H);
            end
        end
        varargout{1} = H;
    case 'assemble'
        E = varargin{1};
        if isstruct(E) && isfield(E,'hist')
            E = assemble(E);
        end
        varargout{1} = E;
    otherwise
        error('ERW_Particle_History:UnknownAction','Unknown action %s.',action);
end
end

function [H,weath] = start(weath,weath_names,keep_days,keep_rows,n_rows)
if nargin < 2; weath_names = fieldnames(weath); end
if nargin < 3 || isempty(keep_days); keep_days = 3; end
if nargin < 4 || isempty(keep_rows); keep_rows = true; end
if nargin < 5; n_rows = []; end
if ~(isnumeric(keep_days) && isscalar(keep_days) && isfinite(keep_days) && ...
        keep_days >= 1 && keep_days == fix(keep_days))
    error('ERW_Particle_History:BadKeepDays', ...
        'keep_days must be a whole number of days, at least one.');
end
candidates = {'psd','d','lamb','SSA','delta_d'};
keep = cellfun(@(n) isfield(weath,n) && any(strcmp(weath_names,n)) && ...
    (size(weath.(n),3) > 1 || (~keep_rows && size(weath.(n),3) == 1)),candidates);
H.names = candidates(keep);
H.blocks = {};
H.keep_days = keep_days;
H.keep_rows = logical(keep_rows);
H.size = struct();
H.initial = struct();
H.initial_rows = struct();
H.first = struct();
H.first_row = struct();
%% Rows kept from dropped blocks, one per day. The rows those days wrote
%% are recorded as ranges, so that a read of any other of them is refused.
H.kept = struct();
H.kept_rows = zeros(0,1);
H.n_kept = 0;
H.dropped = zeros(0,2);
for k = 1:numel(H.names)
    name = H.names{k};
    H.size.(name) = [size(weath.(name),1) size(weath.(name),2) size(weath.(name),3)];
    if ~isempty(n_rows)
        H.size.(name)(1) = n_rows;
    end
    rows = find(any(reshape(weath.(name),size(weath.(name),1),[]),2));
    H.initial_rows.(name) = rows;
    H.initial.(name) = weath.(name)(rows,:,:);
    H.kept.(name) = zeros(0,H.size.(name)(2),H.size.(name)(3),'like',H.initial.(name));
end
for name = intersect({'psd','lamb'},H.names)
    n = name{1};
    napp = size(weath.(n),3);
    H.first.(n) = zeros(1,size(weath.(n),2),napp);
    H.first_row.(n) = inf(1,napp);
    for a = 1:napp
        r = find(weath.(n)(:,end,a) > 0,1);
        if ~isempty(r)
            H.first.(n)(1,:,a) = weath.(n)(r,:,a);
            H.first_row.(n)(a) = r;
        end
    end
end
weath = rmfield(weath,H.names);
end

function v = row(H,day,~,name,r)
%% The latest value written to row r.
if r < 1 || r > H.size.(name)(1)
    error('ERW_Particle_History:RowOutsideHistory', ...
        'Row %d is outside the history''s rows 1-%d.',r,H.size.(name)(1));
end
[v,ok] = latest(H,name,r,day);
if ~ok
    error('ERW_Particle_History:RowNotKept', ...
        ['Row %d of %s was written on a day whose ten-minute rows are no ' ...
         'longer kept; only that day''s last row is.'],r,name);
end
end

function [v,ok] = latest(H,name,r,day)
%% The latest value written to row r: today's block, then the whole blocks
%% from the most recent, then the rows kept from dropped blocks, then the
%% array as initialised. ok is false for a row a dropped block wrote and
%% did not keep.
ok = true;
if nargin > 3
    idx = r - day.first + 1;
    if idx >= 1 && idx <= day.n && day.written(idx)
        v = day.data.(name)(idx,:,:);
        return
    end
end
for b = numel(H.blocks):-1:1
    block = H.blocks{b};
    idx = r - block.first + 1;
    if idx >= 1 && idx <= block.n && block.written(idx)
        v = block.data.(name)(idx,:,:);
        return
    end
end
idx = find(H.kept_rows(1:H.n_kept) == r,1);
if ~isempty(idx)
    v = H.kept.(name)(idx,:,:);
    return
end
if any(r >= H.dropped(:,1) & r <= H.dropped(:,2))
    ok = false;
    v = [];
    return
end
idx = find(H.initial_rows.(name) == r,1);
if isempty(idx)
    v = zeros(1,H.size.(name)(2),H.size.(name)(3),'like',H.initial.(name));
else
    v = H.initial.(name)(idx,:,:);
end
end

function H = refresh_kept(H,day)
%% A day that writes a row kept from a dropped block leaves the kept copy
%% at the value written last, as the array written in place would hold.
if H.n_kept == 0; return; end
written = find(day.written);
[hit,where] = ismember(H.kept_rows(1:H.n_kept),day.first + written - 1);
if ~any(hit); return; end
k = find(hit);
in_day = written(where(hit));
for q = 1:numel(H.names)
    name = H.names{q};
    H.kept.(name)(k,:,:) = day.data.(name)(in_day,:,:);
end
end

function H = drop_oldest(H)
%% Drop the oldest whole block: keep its last written row at its latest
%% value (unless the run keeps no rows), and record the rows it wrote as
%% no longer kept whole.
block = H.blocks{1};
rows = block.first + find(block.written) - 1;
if isfield(H,'keep_rows') && ~H.keep_rows
    H.dropped = add_range(H.dropped,rows);
    H.blocks(1) = [];
    return
end
r = rows(end);
values = cell(1,numel(H.names));
for q = 1:numel(H.names)
    values{q} = latest(H,H.names{q},r);
end
k = find(H.kept_rows(1:H.n_kept) == r,1);
if isempty(k)
    if H.n_kept == numel(H.kept_rows)
        capacity = max(366,2*numel(H.kept_rows));
        H.kept_rows(capacity,1) = 0;
        for q = 1:numel(H.names)
            name = H.names{q};
            H.kept.(name)(capacity,H.size.(name)(2),H.size.(name)(3)) = 0;
        end
    end
    H.n_kept = H.n_kept + 1;
    k = H.n_kept;
    H.kept_rows(k) = r;
end
for q = 1:numel(H.names)
    H.kept.(H.names{q})(k,:,:) = values{q};
end
H.dropped = add_range(H.dropped,rows);
H.blocks(1) = [];
end

function ranges = add_range(ranges,rows)
%% The union of row ranges and the rows given, merged where they touch.
rows = rows(:);
breaks = [0; find(diff(rows) > 1); numel(rows)];
for q = 1:numel(breaks)-1
    ranges(end+1,:) = [rows(breaks(q)+1) rows(breaks(q+1))]; %#ok<AGROW>
end
ranges = sortrows(ranges,1);
merged = ranges(1,:);
for q = 2:size(ranges,1)
    if ranges(q,1) <= merged(end,2) + 1
        merged(end,2) = max(merged(end,2),ranges(q,2));
    else
        merged(end+1,:) = ranges(q,:); %#ok<AGROW>
    end
end
ranges = merged;
end

function E = assemble(E)
%% The rows the column kept, at their latest values: those set at
%% initialisation and the last written row of every ERW day.
H = E.hist;
if ~isempty(H.names)
    candidates = reshape(H.kept_rows(1:H.n_kept),[],1);
    for b = 1:numel(H.blocks)
        candidates(end+1,1) = H.blocks{b}.first + ...
            find(H.blocks{b}.written,1,'last') - 1; %#ok<AGROW>
    end
    for q = 1:numel(H.names)
        candidates = [candidates; H.initial_rows.(H.names{q})(:)]; %#ok<AGROW>
    end
    candidates = unique(candidates);
    servable = true(numel(candidates),1);
    for c = 1:numel(candidates)
        [~,servable(c)] = latest(H,H.names{1},candidates(c));
    end
    rows = candidates(servable);
    for q = 1:numel(H.names)
        name = H.names{q};
        values = zeros(numel(rows),H.size.(name)(2),H.size.(name)(3), ...
            'like',H.initial.(name));
        for c = 1:numel(rows)
            values(c,:,:) = latest(H,name,rows(c));
        end
        E.weath.(name) = values;
    end
    E.weath.particle_rows = rows;
end
E = rmfield(E,'hist');
end

function H = track(H,day,out,r,weath)
%% Keep each application's first row with a positive last class.
for name = intersect({'psd','lamb'},H.names)
    n = name{1};
    if ~isfield(out,n); continue; end
    value = day.data.(n)(r-day.first+1,:,:);
    for a = 1:numel(H.first_row.(n))
        positive = value(1,end,a) > 0;
        if positive && r < H.first_row.(n)(a)
            H.first.(n)(1,:,a) = value(1,:,a);
            H.first_row.(n)(a) = r;
        elseif r == H.first_row.(n)(a)
            if positive
                H.first.(n)(1,:,a) = value(1,:,a);
            else
                %% The first positive row was overwritten with a value that
                %% is not: search the whole history again.
                [H.first.(n)(1,:,a),H.first_row.(n)(a)] = rescan(H,day,weath,n,a);
            end
        end
    end
end
end

function [value,first_row] = rescan(H,day,weath,name,a)
%% A search of the whole history needs every row, so it is possible only
%% while no block has been dropped. A run never needs it: the daily step
%% writes each row once, so an application's first positive row is never
%% overwritten.
if ~isempty(H.dropped)
    error('ERW_Particle_History:RescanNotKept', ...
        ['The first positive row of application %d of %s was overwritten ' ...
         'with a value that is not positive, and the rows to search again ' ...
         'are no longer all kept.'],a,name);
end
column = zeros(H.size.(name)(1),1,'like',H.initial.(name));
column(H.initial_rows.(name)) = H.initial.(name)(:,end,a);
blocks = [H.blocks {day}];
for b = 1:numel(blocks)
    block = blocks{b};
    rows = block.first + find(block.written) - 1;
    column(rows) = block.data.(name)(block.written,end,a);
end
first_row = find(column > 0,1);
value = zeros(1,H.size.(name)(2),'like',H.initial.(name));
if isempty(first_row)
    first_row = inf;
else
    v = row(H,day,weath,name,first_row);
    value = v(1,:,a);
end
end
