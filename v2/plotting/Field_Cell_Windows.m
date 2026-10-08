function [win,tag,words] = Field_Cell_Windows(day_axis,interval,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Field_Cell_Windows
%%
%% Time windows over a run's daily rows (res.day of Field_Cell_Result).
%%
%%   [win,tag,words] = Field_Cell_Windows(res.day,'year')
%%
%% interval
%%   'all'               the whole run
%%   'year'              calendar years (YearStartMonth moves the start,
%%                       e.g. 11 for November-to-October years)
%%   'growing_season'    GrowingSeasonMonths of each year ([5 10]: May-Oct)
%%   'season'            DJF, MAM, JJA, SON
%%   'month', 'week', 'day'
%%   N                   consecutive N-day blocks from the first day
%%   a calendarDuration or duration, e.g. calmonths(3): consecutive blocks
%% Options
%%   'Windows'    [start stop] rows (datetime, or text datetime reads);
%%                replaces interval
%%   'Snapshots'  dates, each its own one-day window; replaces interval
%%   'Period'     [t0 t1]: clip every window to it
%%   'YearStartMonth' (1), 'GrowingSeasonMonths' ([5 10])
%%
%% win(k) has label, start, stop, range_label, short_range, rows (logical
%% over the daily rows), first_row, last_row, prev_row (the row before
%% the window: its state is the state at the window's start) and n_days.
%% Row 1, the initial state, is in no window. Windows the run does not
%% reach are dropped. tag names the windows in file names, words in
%% titles ('per year').
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
p.addParameter('Windows',[]);
p.addParameter('Snapshots',[]);
p.addParameter('Period',[]);
p.addParameter('YearStartMonth',1,@(x) isscalar(x) && any(x == 1:12));
p.addParameter('GrowingSeasonMonths',[5 10], ...
    @(x) numel(x) == 2 && all(ismember(x,1:12)));
p.parse(varargin{:});
o = p.Results;

if ~isdatetime(day_axis) || numel(day_axis) < 2
    error('Field_Cell_Windows:BadDays','The daily axis must be a datetime vector.');
end
day_axis = dateshift(day_axis(:),'start','day');
simulated = true(size(day_axis));
simulated(1) = false;
p0 = min(day_axis(simulated));
p1 = max(day_axis(simulated));
if ~isempty(o.Period)
    period = dateshift(to_datetime(o.Period),'start','day');
    if numel(period) ~= 2 || period(2) < period(1)
        error('Field_Cell_Windows:BadPeriod','Period must be [start stop].');
    end
    p0 = max(p0,period(1));
    p1 = min(p1,period(2));
    if p1 < p0
        error('Field_Cell_Windows:EmptyPeriod', ...
            'The period does not overlap the run (%s to %s).', ...
            day_label(min(day_axis(simulated))),day_label(max(day_axis(simulated))));
    end
end

if ~isempty(o.Snapshots)
    dates = dateshift(to_datetime(o.Snapshots(:)),'start','day');
    outside = dates < p0 | dates > p1;
    if any(outside)
        warning('Field_Cell_Windows:SnapshotOutside', ...
            '%d snapshot date(s) lie outside the run and are left out.',sum(outside));
    end
    starts = dates(~outside);
    stops = starts;
    labels = arrayfun(@day_label,starts,'UniformOutput',false);
    tag = 'snapshots';
    words = 'on the day';
elseif ~isempty(o.Windows)
    W = to_datetime(o.Windows);
    if size(W,2) ~= 2
        error('Field_Cell_Windows:BadWindows','Windows must be [start stop] rows.');
    end
    W = dateshift(W,'start','day');
    starts = max(W(:,1),p0);
    stops = min(W(:,2),p1);
    labels = arrayfun(@range_label,starts,stops,'UniformOutput',false);
    tag = 'custom';
    words = 'per window';
else
    [starts,stops,labels,tag,words] = regular_windows(interval,p0,p1,o);
end

win = struct('label',{},'start',{},'stop',{},'range_label',{}, ...
    'short_range',{},'rows',{},'first_row',{},'last_row',{},'prev_row',{},'n_days',{});
for k = 1:numel(starts)
    if stops(k) < starts(k); continue; end
    rows = simulated & day_axis >= starts(k) & day_axis <= stops(k);
    if ~any(rows); continue; end
    w.label = labels{k};
    w.start = starts(k);
    w.stop = stops(k);
    w.range_label = range_label(starts(k),stops(k));
    w.short_range = short_range(starts(k),stops(k));
    w.rows = rows;
    w.first_row = find(rows,1);
    w.last_row = find(rows,1,'last');
    w.prev_row = max(1,w.first_row-1);
    w.n_days = numel(unique(day_axis(rows)));
    win(end+1) = w; %#ok<AGROW>
end
if isempty(win)
    error('Field_Cell_Windows:NoWindows','No window falls inside the run.');
end
end

%% =======================================================================
function [starts,stops,labels,tag,words] = regular_windows(interval,p0,p1,o)
if isnumeric(interval)
    if ~(isscalar(interval) && interval >= 1 && interval == round(interval))
        error('Field_Cell_Windows:BadInterval','A numeric interval is a whole number of days.');
    end
    starts = (p0:caldays(interval):p1).';
    stops = min(starts+caldays(interval-1),p1);
    labels = arrayfun(@range_label,starts,stops,'UniformOutput',false);
    tag = sprintf('%dd',interval);
    words = sprintf('per %d days',interval);
    return
end
if isa(interval,'calendarDuration') || isa(interval,'duration')
    if ~isscalar(interval)
        error('Field_Cell_Windows:BadInterval','A duration interval must be scalar.');
    end
    starts = (p0:interval:p1).';
    stops = [starts(2:end)-caldays(1); p1];
    stops = dateshift(stops,'start','day');
    labels = arrayfun(@range_label,starts,stops,'UniformOutput',false);
    tag = 'blocks';
    words = 'per block';
    return
end
if ~(ischar(interval) || (isstring(interval) && isscalar(interval)))
    error('Field_Cell_Windows:BadInterval', ...
        'Interval is a name, a number of days or a duration.');
end
names = month_names();
switch lower(char(interval))
    case 'all'
        starts = p0; stops = p1;
        labels = {'Whole run'};
        tag = 'all'; words = 'over the run';
    case 'year'
        ys = o.YearStartMonth;
        b0 = datetime(year(p0),ys,1);
        if b0 > p0; b0 = b0-calyears(1); end
        bounds = (b0:calyears(1):p1).';
        starts = max(bounds,p0);
        stops = min(bounds+calyears(1)-caldays(1),p1);
        if ys == 1
            labels = arrayfun(@(b) sprintf('%d',year(b)),bounds,'UniformOutput',false);
        else
            labels = arrayfun(@(b) sprintf('%d/%02d',year(b),mod(year(b)+1,100)), ...
                bounds,'UniformOutput',false);
        end
        tag = 'year'; words = 'per year';
    case 'growing_season'
        m = o.GrowingSeasonMonths;
        years = (year(p0)-1:year(p1)).';
        first = datetime(years,m(1),1);
        last_year = years + (m(2) < m(1));
        last = datetime(last_year,m(2),1)+calmonths(1)-caldays(1);
        starts = max(first,p0);
        stops = min(last,p1);
        labels = arrayfun(@(y) sprintf('%s-%s %d',names{m(1)},names{m(2)},y), ...
            years,'UniformOutput',false);
        tag = 'growing_season'; words = 'per growing season';
    case 'season'
        season_start = [12 3 6 9];
        season_name = {'DJF','MAM','JJA','SON'};
        mo = month(p0);
        s = floor(mod(mo,12)/3)+1;
        y0 = year(p0) - double(mo <= 2);
        b0 = datetime(y0,season_start(s),1);
        bounds = (b0:calmonths(3):p1).';
        starts = max(bounds,p0);
        stops = min(bounds+calmonths(3)-caldays(1),p1);
        labels = cell(numel(bounds),1);
        for k = 1:numel(bounds)
            sk = floor(mod(month(bounds(k)),12)/3)+1;
            if sk == 1
                labels{k} = sprintf('DJF %d/%02d',year(bounds(k)),mod(year(bounds(k))+1,100));
            else
                labels{k} = sprintf('%s %d',season_name{sk},year(bounds(k)));
            end
        end
        tag = 'season'; words = 'per season';
    case 'month'
        bounds = (datetime(year(p0),month(p0),1):calmonths(1):p1).';
        starts = max(bounds,p0);
        stops = min(bounds+calmonths(1)-caldays(1),p1);
        labels = arrayfun(@(b) sprintf('%s %d',names{month(b)},year(b)),bounds, ...
            'UniformOutput',false);
        tag = 'month'; words = 'per month';
    case 'week'
        starts = (p0:caldays(7):p1).';
        stops = min(starts+caldays(6),p1);
        labels = arrayfun(@range_label,starts,stops,'UniformOutput',false);
        tag = 'week'; words = 'per week';
    case 'day'
        starts = (p0:caldays(1):p1).';
        stops = starts;
        labels = arrayfun(@day_label,starts,'UniformOutput',false);
        tag = 'day'; words = 'per day';
    otherwise
        error('Field_Cell_Windows:BadInterval', ...
            ['Unknown interval "%s": use all, year, growing_season, season, ' ...
            'month, week, day, a number of days or a duration.'],char(interval));
end
end

%% =======================================================================
function t = to_datetime(x)
if isdatetime(x)
    t = x;
elseif isnumeric(x)
    t = datetime(x,'ConvertFrom','datenum');
else
    t = datetime(x);
end
end

%% =======================================================================
function names = month_names()
names = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
end

function text = day_label(d)
names = month_names();
text = sprintf('%d %s %d',day(d),names{month(d)},year(d));
end

function text = range_label(a,b)
names = month_names();
if a == b
    text = day_label(a);
elseif year(a) == year(b)
    text = sprintf('%d %s - %d %s %d',day(a),names{month(a)},day(b),names{month(b)},year(b));
else
    text = [day_label(a) ' - ' day_label(b)];
end
end

function text = short_range(a,b)
names = month_names();
if a == b
    text = sprintf('%d %s',day(a),names{month(a)});
elseif year(a) == year(b)
    text = sprintf('%d %s - %d %s',day(a),names{month(a)},day(b),names{month(b)});
else
    text = range_label(a,b);
end
end
