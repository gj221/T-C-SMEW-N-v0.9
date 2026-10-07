function test_erw_particle_history
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The particle-class histories of an amended column (rows x classes x
%% applications) are kept in daily blocks. Only the last keep_days blocks
%% stay whole; of every older day only its last row is kept. Everything
%% ERW_function reads must stay exactly what the full arrays would give:
%% the previous row, and the first row of each application whose last
%% class is positive (its initial size distribution, which ERW_function
%% finds with psd_0(find(psd_0(:,end,k)>0,1),:,k)), a window that repeats
%% the day before it included. A read the kept rows cannot answer must
%% fail rather than return a wrong value, and the assembled result must
%% hold the kept rows at the values the full arrays end with.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

code_dir = fileparts(fileparts(mfilename('fullpath')));
addpath(code_dir);

nn = 60; ncl = 3; napp = 4; keep_days = 2;
names = {'psd','d','lamb','SSA','delta_d'};
full = struct();
for k = 1:numel(names)
    full.(names{k}) = zeros(nn,ncl,napp);
end
%% Row 1 is the initial state; application 1 is present from the start.
for k = 1:numel(names)
    full.(names{k})(1,:,1) = k + [0.1 0.2 0.3];
end
full.EW = zeros(nn,2);            % two-dimensional fields stay in place
full.psd(40,:,4) = [0.5 0.6 0.7];
full.lamb(40,:,4) = [1.5 1.6 1.7];
[H,weath] = ERW_Particle_History('start',full,fieldnames(full),keep_days);
assert(~any(isfield(weath,names)), ...
    'The live weathering state must not retain the full particle histories.');
assert(isequal(sort(H.names),sort(names)) && H.keep_days == keep_days, ...
    'Only the three-dimensional histories are kept in blocks.');
[H3,~] = ERW_Particle_History('start',full);
assert(H3.keep_days == 3,'The daily step keeps the last three days whole.');
try
    ERW_Particle_History('start',full,fieldnames(full),0);
    error('test_erw_particle_history:MissingError','A history keeping no day was built.');
catch ME
    assert(strcmp(ME.identifier,'ERW_Particle_History:BadKeepDays'), ...
        'At least one day must be kept whole.');
end

%% Before the first application ERW_function returns scalar zeros; the
%% original array assignment expands them across all scheduled cohorts.
scalar_out = struct();
for k = 1:numel(names); scalar_out.(names{k}) = 0; end
before_day = ERW_Particle_History('open',H,2,1);
before_day = ERW_Particle_History('write',H,before_day,scalar_out,2);
before_first = ERW_Particle_History('track',H,before_day,scalar_out,2,weath);
assert(isequal(before_first.first,H.first) && ...
    isequal(before_first.first_row,H.first_row), ...
    'Scalar zeros before an application must preserve the first-row cache.');
assert(isequal(ERW_Particle_History('row',H,before_day,weath,'psd',2), ...
    zeros(1,ncl,napp)), ...
    'A scalar solver output must expand as in the full history.');

rng(7);
%% Days of five rows from row 2, two of them kept whole at a time. The
%% third window repeats rows 5-7 and overwrites row 5, the first positive
%% row of application 2, with a last class of zero, so the history is
%% searched again while every block is still whole. The seventh repeats
%% the day before it, as a final window may. Applications 2, 3 and 4 start
%% at rows 5, 13 and 29, and application 2 goes non-positive in the last
%% class at row 33.
windows = {2:6, 7:11, 5:7, 12:16, 17:21, 22:26, 22:26, 27:31, 32:36, 37:41};
for w = 1:numel(windows)
    rows = windows{w};
    day = ERW_Particle_History('open',H,rows(1),numel(rows));
    for r = rows
        %% What ERW_function would receive at this step.
        for k = {'psd','d','lamb','SSA'}
            expected = full.(k{1})(r-1,:,:);
            got = ERW_Particle_History('row',H,day,weath,k{1},r-1);
            assert(isequal(got,expected), ...
                'Previous row of %s at row %d differs from the full history.',k{1},r);
        end
        for k = {'psd','lamb'}
            first = ERW_Particle_History('first',H,day,k{1});
            for a = 1:napp
                expected = full.(k{1})(find(full.(k{1})(:,end,a)>0,1),:,a);
                got = first(find(first(:,end,a)>0,1),:,a);
                assert(isequal(got,expected), ...
                    'The initial %s of application %d at row %d differs.',k{1},a,r);
            end
        end
        %% This step's output.
        out = struct();
        for k = 1:numel(names)
            v = zeros(1,ncl,napp);
            v(1,:,1) = rand(1,ncl);
            if r >= 5;  v(1,:,2) = rand(1,ncl); end
            if r >= 13; v(1,:,3) = rand(1,ncl); end
            if r >= 29; v(1,:,4) = rand(1,ncl); end
            if r == 33; v(1,end,2) = 0; end
            out.(names{k}) = v;
        end
        if w == 3 && r == 5
            out.psd(1,end,2) = 0;     % overwrites application 2's first positive row
        end
        for k = 1:numel(names)
            full.(names{k})(r,:,:) = out.(names{k});
        end
        day = ERW_Particle_History('write',H,day,out,r);
        H = ERW_Particle_History('track',H,day,out,r,weath);
    end
    H = ERW_Particle_History('close',H,day);
    assert(numel(H.blocks) <= keep_days, ...
        'No more than %d days may stay whole.',keep_days);
    if w == 6
        saved = struct('weath',weath,'hist',H);
        checkpoint = [tempname '.mat'];
        cleanup_file = onCleanup(@() delete(checkpoint));
        save(checkpoint,'saved','-v7.3');
        loaded = load(checkpoint,'saved');
        assert(isequaln(saved,loaded.saved), ...
            'Saving and reloading the history must preserve the live state.');
        H = loaded.saved.hist; weath = loaded.saved.weath;
    end
end

%% A row of a dropped day is refused unless it is that day's last row; a
%% row no day wrote keeps its initial value.
empty_day = ERW_Particle_History('open',H,nn+1,0);
try
    ERW_Particle_History('row',H,empty_day,weath,'psd',10);
    error('test_erw_particle_history:MissingError','A row no longer kept was read.');
catch ME
    assert(strcmp(ME.identifier,'ERW_Particle_History:RowNotKept'), ...
        'A row no longer kept must be refused.');
end
for r = [11 26 35 40]
    assert(isequal(ERW_Particle_History('row',H,empty_day,weath,'psd',r), ...
        full.psd(r,:,:)), ...
        'Kept row %d must hold the value written to it last.',r);
end
assert(isequal(ERW_Particle_History('row',H,empty_day,weath,'psd',50), ...
    zeros(1,ncl,napp)),'A row no day wrote must read as initialised.');

%% A day that writes a row kept from a dropped day updates the kept copy.
refresh_day = ERW_Particle_History('open',H,16,1);
new16 = struct();
for k = 1:numel(names); new16.(names{k}) = rand(1,ncl,napp); end
refresh_day = ERW_Particle_History('write',H,refresh_day,new16,16);
H = ERW_Particle_History('close',H,refresh_day);
for k = 1:numel(names); full.(names{k})(16,:,:) = new16.(names{k}); end
k16 = find(H.kept_rows(1:H.n_kept) == 16);
assert(isscalar(k16) && isequal(H.kept.psd(k16,:,:),full.psd(16,:,:)), ...
    'A rewritten row kept from a dropped day must hold its latest value.');

%% Searching again for a first row needs rows that are no longer kept.
bad_day = ERW_Particle_History('open',H,6,1);
bad = struct();
for k = 1:numel(names); bad.(names{k}) = zeros(1,ncl,napp); end
bad.psd(1,end,2) = -1;
bad_day = ERW_Particle_History('write',H,bad_day,bad,6);
try
    ERW_Particle_History('track',H,bad_day,bad,6,weath);
    error('test_erw_particle_history:MissingError','A search of dropped rows succeeded.');
catch ME
    assert(strcmp(ME.identifier,'ERW_Particle_History:RescanNotKept'), ...
        'Searching again for a first row must fail once rows were dropped.');
end

E = struct('weath',weath,'hist',H);
E = ERW_Particle_History('assemble',E);
assert(~isfield(E,'hist'),'Assembly must leave no block store in the column state.');
%% The rows set at initialisation (1, and 40 for psd and lamb) and the last
%% row of every window, at the values the full arrays end with.
kept_rows = [1; 6; 7; 11; 16; 21; 26; 31; 36; 40; 41];
assert(isequal(E.weath.particle_rows,kept_rows), ...
    'The assembled rows must be the initial rows and each window''s last row.');
for k = 1:numel(names)
    assert(isequal(E.weath.(names{k}),full.(names{k})(kept_rows,:,:)), ...
        'Assembled %s differs from the full history at the kept rows.',names{k});
end
assert(isequal(E.weath.EW,full.EW),'Fields kept in place must be untouched.');
E2 = ERW_Particle_History('assemble',struct('weath',full));
assert(isequal(E2.weath,full),'A column without blocks is left as it is.');
assert(isequal(ERW_Particle_History('close',H,empty_day),H), ...
    'An empty final window must not add a history block.');
try
    ERW_Particle_History('write',H,empty_day,out,nn+1);
    error('test_erw_particle_history:MissingError','An out-of-window write succeeded.');
catch ME
    assert(strcmp(ME.identifier,'ERW_Particle_History:RowOutsideDay'), ...
        'An out-of-window write must report its row bounds.');
end
for r = [0 nn+1]
    try
        ERW_Particle_History('row',H,empty_day,weath,'psd',r);
        error('test_erw_particle_history:MissingError','An invalid history row was read.');
    catch ME
        assert(strcmp(ME.identifier,'ERW_Particle_History:RowOutsideHistory'), ...
            'An invalid history row must report the history bounds.');
    end
end
single_application = struct('psd',full.psd(:,:,1),'lamb',full.lamb(:,:,1));
[no_blocks,unchanged] = ERW_Particle_History('start',single_application);
assert(isempty(no_blocks.names) && isequal(unchanged,single_application), ...
    'A single application must retain the existing two-dimensional storage.');
E3 = ERW_Particle_History('assemble',struct('weath',unchanged,'hist',no_blocks));
assert(~isfield(E3,'hist') && isequal(E3.weath,single_application), ...
    'A column without particle histories must assemble unchanged.');

%% The daily step and the frame use it; ERW_function is unchanged.
step = fileread(fullfile(code_dir,'ERW_Daily_Step_2D.m'));
frame = fileread(fullfile(code_dir,'MAIN_FRAME_BG_ERW_2D.m'));
assert(contains(step,"ERW_Particle_History('start',") && ...
    contains(step,"ERW_Particle_History('first',") && ...
    contains(step,"ERW_Particle_History('write',") && ...
    contains(step,"ERW_Particle_History('close',") && ...
    contains(frame,"ERW{c} = ERW_Particle_History('assemble',ERW{c});"), ...
    'The daily step must keep the histories in blocks and the frame assemble them.');
fprintf('PASS: particle-class histories keep the rows ERW_function reads and each day''s last row.\n');
end
