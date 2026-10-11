function B = Kantola_Base_Cation_Budget(S,c,nd)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Kantola_Base_Cation_Budget
%%
%% Splits the change in the solution plus exchangeable Ca and Mg inventory
%% of one column over the retained days 1..nd into the terms that move it.
%% Everything is returned in mol charge m^-2 (2 charges per Ca or Mg).
%% One mol charge m^-2 equals 10 kmol charge ha^-1.
%%
%% Retained states P are g element m^-2.  Column indices follow
%% Biogeochemistry_Parameter: 56:59 Ca litter and SOM, 60 Ca solution,
%% 61 Ca exchangeable, 62 Ca fixed, 63 Ca primary; 64:71 the same for Mg.
%% Retained daily fluxes are g element m^-2 d^-1 on ground area, except the
%% per-crown uptake series which Biogeo_Per_Crown_Flux divides by Ccrown,
%% so they are multiplied back by Ccrown_t here as in
%% Kantola_N_Magnitude_Audit.
%%
%% The litter and decomposition return to solution, deposition and the
%% fixed-pool release are not retained as daily series, so they appear
%% together as `return_to_solution`, closed from the other terms.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

required = {'P','LEAK_Ca','LEAK_Mg','Cauptake_H','Cauptake_L', ...
    'Mguptake_H','Mguptake_L','AgrHarNut','Ccrown_t'};
for k = 1:numel(required)
    assert(isfield(S,required{k}), ...
        'Kantola_Base_Cation_Budget:MissingField', ...
        'Missing retained field %s.',required{k});
end
if nargin < 2 || isempty(c)
    c = 1;
end
P = column_matrix(S.P,c);
if nargin < 3 || isempty(nd)
    assert(isfield(S,'i'),'Kantola_Base_Cation_Budget:MissingDay', ...
        'A completed day or retained iteration is required.');
    nd = min(size(P,1)-1,floor(double(S.i)/24));
end
assert(isscalar(nd) && isfinite(nd) && nd == floor(nd) && nd >= 1, ...
    'Kantola_Base_Cation_Budget:BadDay','Completed day must be a positive integer.');
assert(size(P,1) >= nd+1 && size(P,2) >= 71, ...
    'Kantola_Base_Cation_Budget:ShortState', ...
    'The retained pool state does not contain day %d.',nd);
P = double(P(1:nd+1,:));
assert(all(isfinite(P(:,56:71)),'all'), ...
    'Kantola_Base_Cation_Budget:InvalidPool','Retained Ca or Mg pools are nonfinite.');

%% g element -> mol charge
qCa = 2/40;
qMg = 2/24.3;
dP = P(end,:) - P(1,:);
charge = @(ca_idx,mg_idx) qCa*sum(dP(ca_idx)) + qMg*sum(dP(mg_idx));

w = column_matrix(S.Ccrown_t,c);
B = struct();
B.column = c;
B.day = nd;
B.delta_solution_exchange = charge(60:61,68:69);
B.delta_fixed = charge(62,70);
B.delta_primary = charge(63,71);
B.delta_litter_som = charge(56:59,64:67);
B.leaching = qCa*sum_daily(column_series(S.LEAK_Ca,c),nd,'LEAK_Ca') + ...
    qMg*sum_daily(column_series(S.LEAK_Mg,c),nd,'LEAK_Mg');
B.plant_uptake = ...
    qCa*(sum_weighted_daily(column_matrix(S.Cauptake_H,c),w,nd,'Cauptake_H') + ...
         sum_weighted_daily(column_matrix(S.Cauptake_L,c),w,nd,'Cauptake_L')) + ...
    qMg*(sum_weighted_daily(column_matrix(S.Mguptake_H,c),w,nd,'Mguptake_H') + ...
         sum_weighted_daily(column_matrix(S.Mguptake_L,c),w,nd,'Mguptake_L'));
har = column_matrix(S.AgrHarNut,c);
B.harvest_reserve_export = qCa*sum_daily(har(:,4),nd,'AgrHarNut(:,4)') + ...
    qMg*sum_daily(har(:,5),nd,'AgrHarNut(:,5)');
%% Primary release is the primary pool drawdown (uplift supply is of order
%% 1e-5 g m^-2 d^-1 and is ignored).
B.primary_release = -B.delta_primary;
B.return_to_solution = B.delta_solution_exchange - B.primary_release + ...
    B.leaching + B.plant_uptake;
B.net_plant_removal = B.plant_uptake - B.return_to_solution;
B.units = 'mol charge m^-2 over days 1..nd; multiply by 10 for kmol charge ha^-1';
end

function x = column_matrix(x,c)
%% (day, k, column) series -> (day, k) for one column; 1-D records are
%% already (day, k).
if ndims(x) == 3
    x = x(:,:,c);
else
    assert(c == 1,'Kantola_Base_Cation_Budget:NoColumn', ...
        'The record holds one column but column %d was requested.',c);
end
end

function x = column_series(x,c)
%% (day, column) scalar series -> (day, 1) for one column.
assert(size(x,2) >= c,'Kantola_Base_Cation_Budget:NoColumn', ...
    'The record holds %d columns but column %d was requested.',size(x,2),c);
x = x(:,c);
end

function value = sum_daily(x,nd,name)
assert(size(x,1) >= nd+1, ...
    'Kantola_Base_Cation_Budget:ShortFlux','%s lacks the state row for day %d.',name,nd);
x = double(x(2:nd+1,:));
assert(all(isfinite(x(:))) && min(x(:)) >= -1e-12, ...
    'Kantola_Base_Cation_Budget:InvalidFlux','%s is nonfinite or negative.',name);
value = sum(x(:));
end

function value = sum_weighted_daily(x,w,nd,name)
assert(size(x,1) >= nd+1 && size(w,1) >= nd+1, ...
    'Kantola_Base_Cation_Budget:ShortFlux', ...
    '%s or Ccrown_t lacks the state row for day %d.',name,nd);
x = double(x(2:nd+1,:));
w = double(w(2:nd+1,:));
assert(isequal(size(x),size(w)), ...
    'Kantola_Base_Cation_Budget:CrownShape', ...
    '%s and Ccrown_t have different retained dimensions.',name);
assert(all(isfinite(x(:))) && min(x(:)) >= -1e-12, ...
    'Kantola_Base_Cation_Budget:InvalidFlux','%s is nonfinite or negative.',name);
assert(all(isfinite(w(:))) && min(w(:)) >= -1e-12 && max(w(:)) <= 1+1e-12, ...
    'Kantola_Base_Cation_Budget:InvalidCrown','Ccrown_t is nonfinite or outside [0,1].');
value = sum(x.*w,'all');
end
