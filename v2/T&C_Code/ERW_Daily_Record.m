function varargout = ERW_Daily_Record(action,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Daily summaries of the ten-minute ERW chemistry, one per ERW day.
%%
%% An ERW day is the rows one daily call of ERW_Daily_Step_2D computes:
%% 01:00 to 00:50 (144 rows at six steps an hour), five rows on the first
%% call and its own entry for the final flush, so a run of NDays has
%% NDays+1 of them. Each is summarised from the values the rows finally
%% hold: the last row of a call is still moved by the lateral and overland
%% transport of the hours that follow and by the next call's anchors, so a
%% call is summarised when the next one begins (and the last after the
%% time loop), exactly as the full ten-minute arrays of a one-column run
%% hold them at the end. The fields and kinds are the list approved for
%% package y (8 Oct 2026):
%%
%%   end    value at the call's last row
%%   mean   mean over the call's rows
%%   total  sum over the call's rows; rates (umol per day) times dtw
%%   max    largest absolute closure residual over the call's rows
%%
%%   S = ERW_Daily_Record('spec')
%%   r = ERW_Daily_Record('summarise',bg,weath,rows,dtw,rock)   one call
%%   D = ERW_Daily_Record('store',D,r,capacity)                  append r
%%   D = ERW_Daily_Record('finish',D)                            trim
%%   [E,r] = ERW_Daily_Record('close',E)    the call E still holds open
%%
%% rows are the call's rows in bg and weath, rock its applied masses
%% [M_rock_add M_rock_add_CaCO3] summed over those rows. A store D holds,
%% per field, one row per ERW day (D.n of them), with D.rows the call's
%% first and last ten-minute rows on the run's grid and D.j its daily row.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch action
    case 'spec'
        varargout{1} = spec();
    case 'summarise'
        varargout{1} = summarise(varargin{:});
    case 'store'
        varargout{1} = store(varargin{:});
    case 'finish'
        varargout{1} = finish(varargin{:});
    case 'close'
        [varargout{1},varargout{2}] = close_pending(varargin{:});
    otherwise
        error('ERW_Daily_Record:UnknownAction','Unknown action %s.',action);
end
end

function S = spec()
%% End-of-day values: the read-outs and the handoff anchors, then the pore
%% water state for maps.
S.end = {'pH_bulk','f_Ca','f_Mg','f_K','f_Na','f_Al','f_H', ...
    'CaCO3','MgCO3','M_min', ...
    'Ca_tot','Mg_tot','K_tot','Si_tot','PO4_tot','NO3_tot','NH4_tot','IC_tot', ...
    'pH','H','Ca','Mg','K','Na','Al_w','An','Si','PO4','NO3','NH4', ...
    'Alk','DIC','CO2_air','CEC_tot'};
S.mean = {'pH','Ca','Mg','Alk','DIC','pCO2_atm','ionic_str','Omega_CaCO3', ...
    'f_Ca','f_Mg','f_K','f_Na','f_Al','f_H'};
%% Rates in umol per day, summed over the day's steps times dtw [umol].
S.rate_total = {'EW','Leak_Na','Leak_Al','Leak_An','Leak_DIC'};
%% Carbonate reaction rates [umol per day], dissolution (positive) and
%% precipitation (negative) apart, each times dtw [umol].
S.split_total = {'W_CaCO3','W_MgCO3'};
%% Per-step amounts, summed as they stand.
S.step_total = {'dAlk_sil','dAlk_carb','dAlk_total','P_weathering_release'};
S.max = {'charge_residual','Al_inventory_residual','f_sum'};
end

function r = summarise(bg,weath,rows,dtw,rock)
S = spec();
last = rows(end);
r = struct('end',struct(),'mean',struct(),'total',struct(),'max',struct());
for k = 1:numel(S.end)
    X = field_of(bg,weath,S.end{k});
    r.end.(S.end{k}) = X(last,:);
end
for k = 1:numel(S.mean)
    X = field_of(bg,weath,S.mean{k});
    r.mean.(S.mean{k}) = mean(X(rows,:),1);
end
for k = 1:numel(S.rate_total)
    X = field_of(bg,weath,S.rate_total{k});
    r.total.(S.rate_total{k}) = sum(X(rows,:),1)*dtw;
end
for k = 1:numel(S.split_total)
    name = S.split_total{k};
    X = field_of(bg,weath,name);
    r.total.([name '_dissolution']) = sum(max(X(rows,:),0),1)*dtw;
    r.total.([name '_precipitation']) = sum(min(X(rows,:),0),1)*dtw;
end
for k = 1:numel(S.step_total)
    X = field_of(bg,weath,S.step_total{k});
    r.total.(S.step_total{k}) = sum(X(rows,:),1);
end
r.total.M_rock_add = rock(1);
r.total.M_rock_add_CaCO3 = rock(2);
r.max.charge_residual = max(abs(bg.charge_residual(rows,:)),[],1);
r.max.Al_inventory_residual = max(abs(bg.Al_inventory_residual(rows,:)),[],1);
f_sum = bg.f_Ca(rows,:)+bg.f_Mg(rows,:)+bg.f_K(rows,:)+bg.f_Na(rows,:)+ ...
    bg.f_Al(rows,:)+bg.f_H(rows,:);
r.max.f_sum = max(abs(f_sum-1),[],1);
end

function X = field_of(bg,weath,name)
if isfield(bg,name)
    X = bg.(name);
elseif isfield(weath,name)
    X = weath.(name);
else
    error('ERW_Daily_Record:MissingField', ...
        'The ERW state has no field %s to summarise.',name);
end
end

function D = store(D,r,capacity)
if isempty(D)
    D = struct('n',0,'rows',nan(capacity,2),'j',nan(capacity,1));
    for kind = {'end','mean','total','max'}
        names = fieldnames(r.(kind{1}));
        for k = 1:numel(names)
            D.(kind{1}).(names{k}) = nan(capacity,size(r.(kind{1}).(names{k}),2));
        end
    end
end
n = D.n + 1;
if n > size(D.rows,1)
    error('ERW_Daily_Record:Capacity', ...
        'More ERW days than the %d the store was sized for.',size(D.rows,1));
end
D.n = n;
D.rows(n,:) = r.rows;
D.j(n) = r.j;
for kind = {'end','mean','total','max'}
    names = fieldnames(r.(kind{1}));
    for k = 1:numel(names)
        D.(kind{1}).(names{k})(n,:) = r.(kind{1}).(names{k});
    end
end
end

function [E,r] = close_pending(E)
%% The last call of a run, summarised from the state it leaves.
r = [];
if isfield(E,'daily_pending') && ~isempty(E.daily_pending)
    pending = E.daily_pending;
    r = summarise(E.bg,E.weath,pending.rows_local,E.dtw,pending.rock);
    r.rows = pending.rows_global;
    r.j = pending.j;
    E.daily_pending = [];
end
end

function D = finish(D)
if isempty(D); return; end
n = D.n;
D.rows = D.rows(1:n,:);
D.j = D.j(1:n);
for kind = {'end','mean','total','max'}
    names = fieldnames(D.(kind{1}));
    for k = 1:numel(names)
        D.(kind{1}).(names{k}) = D.(kind{1}).(names{k})(1:n,:);
    end
end
end
