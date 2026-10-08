function varargout = Plant_Carbon_Deficit(action,B)
% Record, without altering, carbon pools that the integrator drove negative.
%
% Legacy TC-SMEW let B(4) go negative and carried on: the clamped ratio in
% Plant_Exports absorbed it downstream, so carbon spent beyond the pool was
% never counted. Nothing here changes that. This function does not modify B
% and does not affect any computed value; it only records what legacy
% discarded, in the manner of imm_cut in BIOGEOCHEMISTRY_DYNAMIC3.
%
% The distinction it exists to settle: a single overdraw of order 1e-3 gC/m2
% at an exhausted reserve is a boundary artifact of the adaptive step, while
% an overdraw recurring across a season is a systematic carbon sink that the
% campaign's conservation checks cannot see, because they audit ERW mass and
% soil-atmosphere exchange rather than the plant carbohydrate pool.
%
%   Plant_Carbon_Deficit('record',B)    observe one accepted state
%   s = Plant_Carbon_Deficit('report')  accumulated totals
%   Plant_Carbon_Deficit('reset')       clear the counters
%   e = Plant_Carbon_Deficit('take')    the events recorded since the last take
%   Plant_Carbon_Deficit('replay',e)    add taken events to the totals
%
% Call 'record' on the state returned by deval, not inside VEGETATION_DYNAMIC:
% ode45 evaluates the derivative at trial points it may reject, and a negative
% pool at a rejected trial state is not a physical overdraw.
%
% The counters are persistent, so they are per MATLAB process. Each campaign
% run is its own process, which gives the intended per-run isolation, but call
% 'reset' at the start of a run if several share a process.
%
% A thread-pool worker has counters of its own, so the 2-D frame takes each
% column's events inside its daily parfor and replays them on the client in
% column order, the order a serial loop recorded them in. Events wait in a
% pending list until they are taken, replayed or reported, and the totals
% add them one at a time in that order, so 'report' gives the same numbers,
% bit for bit, whichever way they arrived.

persistent n_reserve total_reserve worst_reserve n_other pending_reserve pending_other

if isempty(n_reserve)
    n_reserve = 0; total_reserve = 0; worst_reserve = 0; n_other = 0;
    pending_reserve = zeros(0,1); pending_other = 0;
end

switch lower(char(action))
    case 'record'
        if nargin < 2 || ~isnumeric(B) || numel(B) < 4
            error('Plant_Carbon_Deficit:InvalidState', ...
                'record requires the integrated carbon pool vector B.');
        end
        %% B(4) is the carbohydrate reserve: the pool the clamped ratio hides.
        if isfinite(B(4)) && B(4) < 0
            pending_reserve(end+1,1) = B(4);
        end
        %% B(8) is a working accumulator that is reset by the caller, so the
        %% structural pools are the only other ones worth watching.
        other = B([1:3, 5:min(7,numel(B))]);
        pending_other = pending_other + sum(isfinite(other) & other < 0);
    case 'take'
        varargout{1} = struct('reserve',pending_reserve,'other',pending_other);
        pending_reserve = zeros(0,1); pending_other = 0;
    case 'replay'
        if nargin < 2 || ~isstruct(B) || ~all(isfield(B,{'reserve','other'}))
            error('Plant_Carbon_Deficit:InvalidEvents', ...
                'replay requires the events returned by take.');
        end
        [n_reserve,total_reserve,worst_reserve,n_other] = add_events( ...
            n_reserve,total_reserve,worst_reserve,n_other,B.reserve,B.other);
    case 'report'
        [n_reserve,total_reserve,worst_reserve,n_other] = add_events( ...
            n_reserve,total_reserve,worst_reserve,n_other, ...
            pending_reserve,pending_other);
        pending_reserve = zeros(0,1); pending_other = 0;
        s.reserve_events = n_reserve;
        s.reserve_total_gC_m2 = total_reserve;
        s.reserve_worst_gC_m2 = worst_reserve;
        s.other_pool_events = n_other;
        varargout{1} = s;
    case 'reset'
        n_reserve = 0; total_reserve = 0; worst_reserve = 0; n_other = 0;
        pending_reserve = zeros(0,1); pending_other = 0;
    otherwise
        error('Plant_Carbon_Deficit:UnknownAction', ...
            'action must be record, report, reset, take or replay.');
end
end

function [n_reserve,total_reserve,worst_reserve,n_other] = add_events( ...
    n_reserve,total_reserve,worst_reserve,n_other,reserve,other)
for k = 1:numel(reserve)
    n_reserve = n_reserve + 1;
    total_reserve = total_reserve + reserve(k);
    worst_reserve = min(worst_reserve,reserve(k));
end
n_other = n_other + other;
end
