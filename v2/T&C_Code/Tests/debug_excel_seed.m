% Debug: run the REAL excel seed path standalone (fresh state, like a
% first-ever day-1 call) and print the seed-Alw decisions.
clear; clc;
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..')); addpath(fullfile(here,'..','Parameters')); cd(here);

S = load(fullfile('..','..','SIM','EQUIV_1D_TMP','TC_ERW_results_Kantola_crop_control.mat'));

t = 1;
TT = S.TT_m(1,:);
ERW_TC_const = struct();
ERW_TC_const.n = S.ERW_const.n;
ERW_TC_const.rsd_bio = S.ERW_const.rho_bulk/1e3;
ERW_TC_const.conv_mol = S.conv_mol;
ERW_TC_const.conv_Al = S.conv_Al;
ERW_TC_const.PH_in = 6.1;
ERW_TC_const.Zbio = S.ERW_const.Zrw*1000;
ERW_TC_const.M_rock_add = 0;
ERW_TC_const.M_rock_add_CaCO3 = 0;
ERW_TC_const.cov_moll = 1;

ERW_biogeochem = struct();
ERW_biogeochem.Mg(1) = 2268.73; ERW_biogeochem.Si(1) = 1;
ERW_biogeochem.Ca(1) = 530.237; ERW_biogeochem.K(1) = 7675.04;
ERW_biogeochem.PO4(1) = 54053.8; ERW_biogeochem.NO3(1) = 481.263;
ERW_biogeochem.NH4(1) = 3088.98;
% exchangeable pools feeding the Kantola preprocessor (per ERW_Daily_Step_2D:164-166)
ERW_biogeochem.excCa = 530.237; ERW_biogeochem.excMg = 2268.73; ERW_biogeochem.excK = 7675.04;

ERW_PARAM_input = struct();
xl = fullfile('..','..','T&C_Code','Inputs','ERW_MOD_PARAM_parameters.xlsx');
ERW_PARAM_input = readtable(xl, 'Sheet', 'ERW_MOD');
ERW_PARAM_input = ERW_PARAM_input(4,:);          % case 5 (Kantola control)
ERW_PARAM_input.f_CEC_in = string("NaN");
for rc = {'rock_filepath','rock_type','mineral','rock_f_in'}
    if ismember(rc{1}, ERW_PARAM_input.Properties.VariableNames)
        ERW_PARAM_input.(rc{1}) = string("NaN");
    end
end
erw_cols = {'Al_w_in','Si_in','Ca_in','Mg_in','K_in','Na_in','PO4_in','NO3_in','NH4_in'};
for c = 1:numel(erw_cols)
    if ~any(strcmp(ERW_PARAM_input.Properties.VariableNames, erw_cols{c}))
        ERW_PARAM_input.(erw_cols{c}) = NaN;
    end
end
ERW_PARAM_input.opt_TC_biochem = 1;
ERW_PARAM_input.num_rock_app = sum(ERW_TC_const.M_rock_add > 0);

fprintf('=== standalone day-1 excel seed ===\n');
try
    [B, P, W, EC, EV, bn, wn, pn, cn, vn, A] = ERW_MOD_PARAM_excel( ...
        t, TT, ERW_TC_const, ERW_biogeochem, ERW_PARAM_input, 1);
    bpp = EC.bulk_pH_params;
    fprintf('initial_Al_w      = %.6e nmol/L\n', bpp.initial_Al_w);
    fprintf('scheme            = %s  alpha=%g\n', bpp.scheme, bpp.appelo.alpha);
    fprintf('n_extract         = %g\n', bpp.n_extract);
    fprintf('seed: pH(pore)=%.4f pH_bulk=%.4f An=%.6g Al_w(1)=%.6e Al(1)=%.6g\n', ...
        B.pH(1), B.pH_bulk(1), B.An(1), B.Al_w(1), B.Al(1));
    fprintf('f_CEC seed = [%s]\n', sprintf('%.4f ', B.f_CEC(1,:)));
    fprintf('Ca=%.4g Mg=%.4g K=%.4g Na=%.4g\n', B.Ca(1), B.Mg(1), B.K(1), B.Na(1));
catch ME
    fprintf('EXCEL FAILED: %s\n', ME.message);
    for s = 1:numel(ME.stack)
        fprintf('  at %s line %d\n', ME.stack(s).name, ME.stack(s).line);
    end
end