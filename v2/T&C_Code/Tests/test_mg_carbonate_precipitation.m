function test_mg_carbonate_precipitation
%%% Nesquehonite (MgCO3.3H2O) is the Mg carbonate that precipitates from
%%% soil water at field temperature; magnesite does not form below about
%%% 60 C (Hanchen et al., 2008; Harrison et al., 2019). Without that sink
%%% the 35-year basalt runs carried pore-water Mg to 0.12 mol/L and pH 9.8
%%% by 2036, the Davies ionic strength passed 1.9 mol/L and the activity
%%% guard in ERW_Calcite_Kinetics stopped the lower-set run 11.

test_dir = fileparts(mfilename('fullpath'));
code_dir = fullfile(test_dir,'..');
addpath(code_dir,fullfile(code_dir,'Inputs'),fullfile(code_dir,'Parameters'));
conv_mol = 1e6;

%% Closed-water bound: the calcite limiter applies unchanged to a 1:1
%% Mg carbonate. Start from the lower-set run 11 median of 2036
%% (Mg 0.12 mol/L, CO3 0.056 mol/L, pH 9.8, I 0.4 mol/L) at 25 C.
k1 = 10^(-6.35)*conv_mol;
k2 = 10^(-10.33)*conv_mol;
k_w = 10^(-14)*conv_mol^2;
K_nsq = 10^(-5.27)*conv_mol^2;
gamma_IAP = ERW_Activity_Coefficient(2,0.4)^2;
H = 10^(-9.8)*conv_mol;
Mg = 1.2e5;
CO3 = 5.6e4;
DIC = CO3*(H^2+k1*H+k1*k2)/(k1*k2);
Alk = DIC*(k1*H+2*k1*k2)/(H^2+k1*H+k1*k2)-H+k_w/H;
Omega0 = gamma_IAP*Mg*CO3/K_nsq;
water_l = 100;
assert(Omega0 > 50,'The 2036 pore water must be far above nesquehonite saturation.');
[primary,W_nsq,lim] = ERW_Calcite_Reaction_Limit(0,-1e9,1,Omega0,Mg,DIC,Alk, ...
    water_l,k1,k2,k_w,K_nsq,conv_mol,gamma_IAP);
dMg = lim.Ca_final-Mg;
assert(primary == 0 && W_nsq < 0 && lim.limited && abs(lim.Omega_final-1) <= 1e-8, ...
    'Nesquehonite precipitation must stop at saturation.');
assert(abs(lim.DIC_final-DIC-dMg) <= 1e-9*DIC && abs(lim.Alk_final-Alk-2*dMg) <= 1e-9*Alk && ...
    abs(dMg*water_l-W_nsq) <= 1e-9*abs(W_nsq), ...
    'Nesquehonite removes Mg, DIC and alkalinity in the ratio 1:1:2.');
%% The kernel forms Omega before calcite takes its carbonate, so the value
%% can be stale: above 1 on water that is by then undersaturated.
[~,W_stale] = ERW_Calcite_Reaction_Limit(0,-1e9,1,Omega0,1e3,DIC,Alk, ...
    water_l,k1,k2,k_w,K_nsq,conv_mol,gamma_IAP);
assert(gamma_IAP*1e3*CO3/K_nsq < 1 && W_stale == 0, ...
    'A stale supersaturated Omega must not precipitate from undersaturated water.');

%% Kernel replay from the stored ERW_function inputs of step 22053
%% (test_erw_appelo_exchange_floor), with the precipitation coefficient
%% the parameter file now sets.
erw_source = fileread(fullfile(code_dir,'ERW_function.m'));
param_source = fileread(fullfile(code_dir,'Parameters','ERW_MOD_PARAM_excel.m'));
F = load(fullfile(test_dir,'Fixtures','erw_appelo_al_floor_step22053.mat'),'fixture');
F = F.fixture;
header_start = regexp(erw_source,'ERW_function\(','end','once') + 1;
header_end = regexp(erw_source,'\)\r?\n\r?\n%%%% INPUT','start','once') - 1;
arg_names = strtrim(strsplit(regexprep( ...
    erw_source(header_start:header_end),'\.\.\.',''),','));
args = cell(1,numel(arg_names));
for k = 1:numel(arg_names)
    if ~strcmp(arg_names{k},'~')
        args{k} = F.inputs.(arg_names{k});
    end
end
arg = @(name) find(strcmp(arg_names,name),1);
args{arg('An_tot')} = F.inputs.An_tot*0.98;
args{arg('r_MgCO3')} = F.inputs.r_CaCO3;
warn_state = warning('off','ERW:NoConverge');
cleanup = onCleanup(@() warning(warn_state));

%% Solubility follows Harrison et al. (2019), Chem. Geol. 504, 123-135:
%% log Ksp -5.03, -5.27 and -5.34 (3 sigma 0.13, 0.15, 0.04) at 5, 25 and
%% 35 C, whatever K_MgCO3 the caller passes.
for case_T = [5 -5.03 0.13; 25 -5.27 0.15; 35 -5.34 0.04]'
    T_args = args;
    T_args{arg('temp_soil')} = case_T(1);
    T_args{arg('K_MgCO3')} = 10^(-7.46)*conv_mol^2;
    D = ERW_function(T_args{:});
    assert(abs(log10(D.K_MgCO3/conv_mol^2)-case_T(2)) <= case_T(3), ...
        'K_MgCO3 at %g C is 10^%.3f, outside Harrison et al. (2019).', ...
        case_T(1),log10(D.K_MgCO3/conv_mol^2));
end

%% Undersaturated water with no solid: no Mg carbonate flux.
[D,~,~,W] = ERW_function(args{:});
assert(D.Omega_MgCO3 < 1 && W.W_MgCO3 == 0, ...
    'Undersaturated water without nesquehonite must carry no Mg carbonate flux.');

%% Drifted state of the basalt runs: calcite has taken the Ca, and
%% weathering has added Mg and inorganic carbon 1:1 until nesquehonite is
%% about 150 times supersaturated at I = 0.46 mol/L.
args{arg('Ca_tot')} = 3e5;
args{arg('Mg_tot')} = F.inputs.Mg_tot+F.inputs.Ca_tot-3e5+1e7;
args{arg('IC_tot')} = F.inputs.IC_tot+1e7;
for name = {'I_Ca','I_Mg','I_Na','I_K','I_An','dCa_tot','dMg_tot','dK_tot', ...
        'dSi_tot','dPO4_tot','dNO3_tot','dNH4_tot'}
    args{arg(name{1})} = 0;
end
state_names = {'Na','An','CO2_air','Ca_tot','Mg_tot','K_tot', ...
    'Na_tot','An_tot','Al_tot','Si_tot','Al','AlOH4','R_alk','H','EW', ...
    'IC_tot','Wr','Fs','DIC','W_CaCO3','W_MgCO3','CaCO3','MgCO3','f_CEC', ...
    'f_Ca','f_Mg','f_K','f_Na','f_Al','f_H','M_min','PO4_tot','NO3_tot', ...
    'NH4_tot','CEC_tot'};
n_steps = 12;
I = zeros(1,n_steps);
Omega = zeros(1,n_steps);
Mg_total = zeros(1,n_steps);
for step = 1:n_steps
    [D,~,B,W] = ERW_function(args{:});
    I(step) = B.ionic_str;
    Omega(step) = B.Omega_MgCO3;
    Mg_total(step) = B.Mg_tot+B.MgCO3;
    assert(W.W_MgCO3 < 0 && D.gamma_IAP <= 1, ...
        'Step %d: supersaturated water must precipitate nesquehonite.',step);
    [~,W_again,again] = ERW_Calcite_Reaction_Limit(0,W.W_MgCO3,D.dtw, ...
        B.Omega_MgCO3,B.Mg,D.calcite_limit.DIC_final,D.calcite_limit.Alk_final, ...
        D.n*D.Zrw*D.s*1000,D.k1,D.k2,D.k_w,D.K_MgCO3,conv_mol,D.gamma_IAP);
    assert(W_again == W.W_MgCO3 && again.Omega_final >= 1-1e-8, ...
        'Step %d: the kernel flux must already lie inside the closed-water bound.',step);
    for name = state_names
        if isfield(B,name{1})
            args{arg(name{1})} = B.(name{1});
        elseif isfield(W,name{1})
            args{arg(name{1})} = W.(name{1});
        end
    end
    args{arg('time_step')} = args{arg('time_step')}+1;
end
assert(Omega(1) > 100 && all(diff(Omega) < 0) && all(diff(I) < 0), ...
    'Nesquehonite precipitation must lower the saturation state and ionic strength.');
assert(max(abs(Mg_total-Mg_total(1))) <= 1e-9*Mg_total(1), ...
    'Dissolved, exchangeable and carbonate Mg must be conserved.');

%% Production sources: the parameter file sets the same solubility at the
%% seed temperature and the calcite precipitation coefficient, and both
%% callers bound nesquehonite after calcite.
for needed = {'K_MgCO3 = 10^(-5.27)*exp(17170/R*(1/T_K-1/298.15))*conv_mol^2;', ...
        '0,W_MgCO3,dtw,Omega_MgCO3,Mg,calcite_limit.DIC_final,calcite_limit.Alk_final, ...', ...
        'n*Zrw*s*1000,k1,k2,k_w,K_MgCO3,conv_mol,gamma_IAP);'}
    assert(contains(erw_source,needed{1}),'ERW_function.m must carry: %s',needed{1});
end
for needed = {'K_MgCO3 = 10^(-5.27)*exp(17170/R*(1/T_K(1)-1/298.15))*conv_mol^2;', ...
        'r_MgCO3 = r_CaCO3;', ...
        '0,W_MgCO3(1),dtw_init,Omega_MgCO3(1),Mg(1),calcite_limit.DIC_final,calcite_limit.Alk_final, ...', ...
        'n(1)*Zrw*s(1)*1000,k1_s,k2_s,k_w_s,K_MgCO3,conv_mol,gamma_IAP_0);'}
    assert(contains(param_source,needed{1}),'ERW_MOD_PARAM_excel.m must carry: %s',needed{1});
end

fprintf('PASS: nesquehonite caps Mg; Omega %.3g -> %.3g, I %.3f -> %.3f mol/L over %d steps.\n', ...
    Omega(1),Omega(end),I(1),I(end),n_steps);
end
