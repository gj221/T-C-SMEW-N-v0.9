function calcite_check(mat_file)
%% Carbonate state of one Kantola result (column 1): monthly pore and bulk
%% pH, calcite saturation, Ca, soil pCO2, exchanger acidity, applied and
%% dissolved primary calcite and the secondary (pedogenic) CaCO3 pool.
S = load(mat_file,'Datam','ERW','PH_bulk_d');
bg = S.ERW{1}.bg; weath = S.ERW{1}.weath; cst = S.ERW{1}.const; var_ = S.ERW{1}.var;
TC = S.ERW{1}.TC_const;
n_t = numel(bg.pH); n_d = numel(S.PH_bulk_d);
t0 = datenum(S.Datam(1,1),S.Datam(1,2),S.Datam(1,3));
ts = t0 + (1:n_t)'/(n_t/n_d); vs = datevec(ts);
td = t0 + (0:n_d-1)'; vd = datevec(td);
dt_d = mean(diff(ts));                         % step length, days
names = string(var_.mineral(:)');
kc = find(strcmpi(names,'calcite'));
fprintf('minerals: %s; calcite column %s\n', strjoin(names,' '), mat2str(kc));
if isfield(TC,'M_rock_add_CaCO3'), add2 = sum(TC.M_rock_add_CaCO3)/1e6; else, add2 = NaN; end
if ~isempty(kc)
    EWc = cumsum(weath.EW(:,kc))*dt_d/1e6;      % mol m-2 dissolved (primary)
    Mc = bg.M_min(:,kc);                        % g m-2 remaining (primary)
else
    EWc = zeros(n_t,1); Mc = zeros(n_t,1);
end
fprintf('direct-to-secondary carbonate added: %.3g mol m-2\n', add2);
fprintf('%7s %6s %6s %7s %7s %7s %6s %7s %9s %9s %9s\n','month','pHpw','pHblk','OmCaCO3','Ca mM','Alk meq','pCO2%','f_H%','calc_g_m2','diss_mol','sec_mol');
yrs = unique(vs(:,1))';
for y = yrs
    for m = 1:12
        w = vs(:,1) == y & vs(:,2) == m; if ~any(w), continue; end
        wd = vd(:,1) == y & vd(:,2) == m;
        k = find(w,1,'last');
        fprintf('%4d-%02d %6.2f %6.2f %7.3f %7.2f %7.3f %6.2f %7.2f %9.1f %9.3f %9.3f\n', y, m, ...
            mean(bg.pH(w)), mean(S.PH_bulk_d(wd)), mean(bg.Omega_CaCO3(w)), mean(bg.Ca(w))/1e3, ...
            mean(bg.Alk(w))/1e3, 100*mean(bg.pCO2_atm(w)), 100*mean(bg.f_H(w)), Mc(k), EWc(k), bg.CaCO3(k)/1e6);
    end
end
fprintf('(calc_g_m2: primary calcite remaining; diss_mol: cumulative primary calcite dissolved, mol m-2;\n');
fprintf(' sec_mol: secondary CaCO3 pool, mol m-2. Omega on the activity basis as the solver uses it.)\n');
if isfield(bg,'W_CaCO3') || isfield(weath,'W_CaCO3')
    if isfield(weath,'W_CaCO3'), Wc = weath.W_CaCO3; else, Wc = bg.W_CaCO3; end
    fprintf('secondary CaCO3: net precipitated (negative W) %.3g mol m-2, net dissolved %.3g mol m-2 over the run\n', ...
        -sum(Wc(Wc<0))*dt_d/1e6, sum(Wc(Wc>0))*dt_d/1e6);
end
end
