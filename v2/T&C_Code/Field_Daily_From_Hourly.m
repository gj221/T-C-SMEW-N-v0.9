function d = Field_Daily_From_Hourly(X,kind,dth)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% One daily value per square from a day of hourly records (package y).
%%
%% X holds the day's hours, in order: hours x k x squares. The day is the
%% one a daily row covers: the hours max(1,i-24):i-1 its daily step
%% averaged (row 2 is the first hour alone). The 2-D frame passes the rows
%% of its hourly window; the x-y comparison passes the same hours of a
%% full record, so both reduce the identical array.
%%
%%   'mean'        mean over the hours                   O, ZWT
%%   'rate_total'  sum of rate x dth (mm h-1 -> mm)      Lk, f, q_runon
%%   'total'       sum of per-step amounts (mm)          Rh, Rd
%%
%% d is 1 x k x squares.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch kind
    case 'mean'
        d = mean(X,1);
    case 'rate_total'
        d = sum(X*dth,1);
    case 'total'
        d = sum(X,1);
    otherwise
        error('Field_Daily_From_Hourly:UnknownKind','Unknown kind %s.',kind);
end
end
