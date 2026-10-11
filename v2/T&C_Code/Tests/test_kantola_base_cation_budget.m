%% Synthetic record: known fluxes must come back in mol charge m^-2 and the
%% closed return term must equal the litter return that built the pools.
here = fileparts(mfilename('fullpath'));
addpath(here);

nd = 3; cc = 2;
qCa = 2/40; qMg = 2/24.3;
leak_ca = [0 0.10 0.20 0.30]';          % g Ca m^-2 d^-1, row 1 unused
leak_mg = [0 0.05 0.05 0.05]';
up_ca   = [0 0; 0.2 0.1; 0.2 0.1; 0.2 0.1];   % per crown
up_mg   = [0 0; 0.1 0.0; 0.1 0.0; 0.1 0.0];
crown   = repmat([0.5 0.5],nd+1,1);
ret_ca  = 0.4; ret_mg = 0.2;            % litter return to solution per day
prim_ca = 0.01; prim_mg = 0.02;          % primary release per day

P = zeros(nd+1,79);
P(1,60) = 1; P(1,61) = 5; P(1,63) = 500;
P(1,68) = 1; P(1,69) = 3; P(1,71) = 500;
for j = 2:nd+1
    P(j,:) = P(j-1,:);
    P(j,60) = P(j-1,60) + ret_ca + prim_ca - leak_ca(j) - sum(up_ca(j,:).*crown(j,:));
    P(j,68) = P(j-1,68) + ret_mg + prim_mg - leak_mg(j) - sum(up_mg(j,:).*crown(j,:));
    P(j,63) = P(j-1,63) - prim_ca;
    P(j,71) = P(j-1,71) - prim_mg;
end

S = struct('P',P,'LEAK_Ca',leak_ca,'LEAK_Mg',leak_mg, ...
    'Cauptake_H',zeros(nd+1,cc),'Cauptake_L',up_ca, ...
    'Mguptake_H',zeros(nd+1,cc),'Mguptake_L',up_mg, ...
    'AgrHarNut',zeros(nd+1,6),'Ccrown_t',crown,'i',24*nd);
S.AgrHarNut(4,4) = 0.5; S.AgrHarNut(4,5) = 0.25;

B = Kantola_Base_Cation_Budget(S,1,nd);
tol = 1e-12;
assert(abs(B.leaching - (qCa*0.6 + qMg*0.15)) < tol,'leaching');
assert(abs(B.plant_uptake - (qCa*3*0.15 + qMg*3*0.05)) < tol,'uptake');
assert(abs(B.primary_release - (qCa*3*prim_ca + qMg*3*prim_mg)) < tol,'primary');
assert(abs(B.harvest_reserve_export - (qCa*0.5 + qMg*0.25)) < tol,'harvest');
assert(abs(B.return_to_solution - (qCa*3*ret_ca + qMg*3*ret_mg)) < tol,'closure');
assert(abs(B.delta_fixed) < tol && abs(B.delta_litter_som) < tol,'inert pools');

%% Day inferred from S.i when nd is omitted
B2 = Kantola_Base_Cation_Budget(S,1);
assert(B2.day == nd,'day from S.i');

%% A 2-D record with two columns reads the requested column
S2 = S;
S2.P = cat(3,zeros(nd+1,79),P);
S2.LEAK_Ca = [zeros(nd+1,1) leak_ca]; S2.LEAK_Mg = [zeros(nd+1,1) leak_mg];
S2.Cauptake_L = cat(3,zeros(nd+1,cc),up_ca); S2.Mguptake_L = cat(3,zeros(nd+1,cc),up_mg);
S2.Cauptake_H = zeros(nd+1,cc,2); S2.Mguptake_H = zeros(nd+1,cc,2);
S2.AgrHarNut = cat(3,zeros(nd+1,6),S.AgrHarNut);
S2.Ccrown_t = cat(3,crown,crown);
B3 = Kantola_Base_Cation_Budget(S2,2,nd);
assert(abs(B3.leaching - B.leaching) < tol && ...
    abs(B3.return_to_solution - B.return_to_solution) < tol,'column select');

%% Missing field and bad column are errors
bad = rmfield(S,'LEAK_Mg');
try
    Kantola_Base_Cation_Budget(bad,1,nd); error('no error');
catch ME
    assert(strcmp(ME.identifier,'Kantola_Base_Cation_Budget:MissingField'));
end
try
    Kantola_Base_Cation_Budget(S,2,nd); error('no error');
catch ME
    assert(strcmp(ME.identifier,'Kantola_Base_Cation_Budget:NoColumn'));
end
fprintf('PASS: Kantola_Base_Cation_Budget closes a synthetic record.\n');
