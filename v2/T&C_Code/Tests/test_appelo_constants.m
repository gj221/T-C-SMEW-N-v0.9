function test_appelo_constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SCEPTER form of the Appelo H-referenced selectivity constants.
%%
%% Ground truth is scepter.f90: keqiex_all = 10^(-logkhaq_all) with
%% defaults logKH = [5.9 4.8 10.47 10.786 16.47] (Na K Ca Mg Al), alpha
%% = 3.4, and the
%% mass action f[X-Na] = K*[Na+]*f[X-H]/[H+] (get_maqads_all_v3); the
%% Appelo (1994) proton-loading term 10^(alpha*Z*f_H) is applied at use
%% and no part of alpha is folded into the constant.  The tabulated
%% logKH already carries the Z multiplicity (log KCa\H = 2*logKCa\Na -
%% 2*logKH\Na = -10.47 for the full 2-site reaction).
%%
%% An inverted form (10^(+Z*(logKH-alpha))) overstates Al selectivity by
%% ~56 orders of magnitude, Al-saturates the exchanger at pH 6, and makes
%% every base-cation pool unhostable - the failure cascade behind the
%% step-2 "no mass-conserving exchange solution" aborts.  This test pins
%% the constants and the physical equilibrium they imply.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n=== test_appelo_constants ===\n');
n_pass = 0;

v2 = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(v2, 'T&C_Code'));

[K_Na_H, K_K_H, K_Ca_H2, K_Mg_H2, K_Al_H3, alpha, logKH] = ...
    Resolve_Appelo_Constants([]);

%% ---- exact scepter.f90 values (umol/L basis carries 1e6 per extra charge)
tol = 1e-12;
check(abs(K_Na_H /10^(-5.9)        - 1) < tol, 'K_Na_H  = 10^(-5.9)');
check(abs(K_K_H  /10^(-4.8)        - 1) < tol, 'K_K_H   = 10^(-4.8)');
check(abs(K_Ca_H2/(10^(-10.47) *1e6)  - 1) < tol, 'K_Ca_H2 = 10^(-10.47)*1e6');
check(abs(K_Mg_H2/(10^(-10.786)*1e6)  - 1) < tol, 'K_Mg_H2 = 10^(-10.786)*1e6');
check(abs(K_Al_H3/(10^(-16.47) *1e12) - 1) < tol, 'K_Al_H3 = 10^(-16.47)*1e12');
check(alpha == 3.4, 'alpha = 3.4 (SCEPTER default exchanger)');
check(isequal(logKH,[5.9 4.8 10.47 10.786 16.47]), ...
    'effective SCEPTER default logKH is returned');

%% ---- SCEPTER field fit uses the reported static constants
% Table 9 reports Na, K, Ca and Mg.  Al was not fitted there, so retain
% the SCEPTER default Al constant rather than extending the fit.
table9_logKH = [4.027 4.474 9.032 8.704 16.47];
[K_Na_K,K_K_K,K_Ca_K,K_Mg_K,K_Al_K,alpha_K,logKH_K] = ...
    Resolve_Appelo_Constants(struct('logKH',table9_logKH,'alpha',1.3));
check(alpha_K == 1.3, 'SCEPTER Table 9 alpha is the explicit fitted value');
check(isequal(logKH_K,table9_logKH), ...
    'effective SCEPTER Table 9 logKH is returned');
check(abs(K_Na_K/10^(-table9_logKH(1))-1) < tol, 'SCEPTER Table 9 Na logKH is explicit');
check(abs(K_K_K/10^(-table9_logKH(2))-1) < tol, 'SCEPTER Table 9 K logKH is explicit');
check(abs(K_Ca_K/(10^(-table9_logKH(3))*1e6)-1) < tol, 'SCEPTER Table 9 Ca logKH is explicit');
check(abs(K_Mg_K/(10^(-table9_logKH(4))*1e6)-1) < tol, 'SCEPTER Table 9 Mg logKH is explicit');
check(abs(K_Al_K/(10^(-16.47)*1e12)-1) < tol, ...
    'SCEPTER Table 9 retains the unfitted default Al logKH');
check(max(abs(table9_logKH(1:4)-[4.027 4.474 9.032 8.704])) < tol, ...
    'SCEPTER Table 9 uses the reported field-fit cation constants');

%% ---- H selectivity direction: Appelo organic exchanger is H-loaded at pH 4
% Site sum 1 = f_H*(1 + K_Na*[Na]/[H]*gamma) for the Fig-2a-like Na-only
% system ([Na] = 1000 uM, pH 4): solve by fixed point.
H = 10^(-4)*1e6;   % umol/L
Na = 1000;
f_H = 0.9;
for it = 1:200
    gamma = 10^(alpha*f_H);
    f_H = 1/(1 + K_Na_H*Na/H*gamma);
end
check(f_H > 0.5, sprintf('Na-only exchanger at pH 4 is H-loaded (f_H = %.3f)', f_H));

%% ---- Kantola porewater equilibrium is base-saturated with negligible Al
% Porewater-like: Ca 1450, Mg 2270, K 7675, Na 500 uM at pH 6, with Al3+
% at amorphous Al(OH)3 saturation (10^10.8*[H+]^3 M = 0.039 uM).
H = 10^(-6)*1e6;
Ca = 1450; Mg = 2270; K = 7675; Na = 500;
Al3 = 10^(10.8) * (10^-6)^3 * 1e6;
% Site sum g(f_H) = f_H + a1*f_H + a2*f_H^2 + a3*f_H^3 - 1 with
% gamma = 10^(alpha*f_H) is monotone increasing on (0,1]: bisect.
site_gap = @(fh) fh ...
    + (K_Na_H*Na/H + K_K_H*K/H)*10^(alpha*fh)*fh ...
    + (K_Ca_H2*Ca/H^2 + K_Mg_H2*Mg/H^2)*10^(2*alpha*fh)*fh^2 ...
    + K_Al_H3*Al3/H^3*10^(3*alpha*fh)*fh^3 - 1;
lo = 0; hi = 1;
for it = 1:200
    f_H = (lo + hi)/2;
    if site_gap(f_H) > 0, hi = f_H; else, lo = f_H; end
end
gamma = 10^(alpha*f_H);
f_Na = K_Na_H*Na/H*gamma*f_H;
f_K  = K_K_H*K/H*gamma*f_H;
f_Ca = K_Ca_H2*Ca/H^2*gamma^2*f_H^2;
f_Mg = K_Mg_H2*Mg/H^2*gamma^2*f_H^2;
f_Al = K_Al_H3*Al3/H^3*gamma^3*f_H^3;
site_sum = f_H + f_Na + f_K + f_Ca + f_Mg + f_Al;
check(abs(site_sum - 1) < 1e-6, 'site sum closes at the fixed point');
base = f_Ca + f_Mg + f_K + f_Na;
check(base > 0.5, sprintf('Kantola pH-6 exchanger is base-saturated (BS = %.2f)', base));
check(f_Al < 1e-3, sprintf('Al occupies a negligible site fraction (f_Al = %.2e)', f_Al));
check(f_H > 0.05 && f_H < 0.6, sprintf('proton loading is moderate (f_H = %.3f)', f_H));

fprintf('  f_H=%.3f f_Ca=%.3f f_Mg=%.3f f_K=%.3f f_Na=%.4f f_Al=%.2e\n', ...
    f_H, f_Ca, f_Mg, f_K, f_Na, f_Al);
fprintf('  ---------------------------------------\n');
fprintf('  %d passed, 0 failed\n', n_pass);

    function check(ok, label)
        assert(ok, 'test_appelo_constants:Failed', '%s', label);
        n_pass = n_pass + 1;
        fprintf('  PASS  %s\n', label);
    end

end
