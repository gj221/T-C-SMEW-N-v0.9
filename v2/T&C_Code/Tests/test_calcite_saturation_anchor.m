function test_calcite_saturation_anchor
%%% External anchor on the calcite-CO2-water system.
%%%
%%% Every other calcite test in this suite checks the code against itself:
%%% the source pins, the Davies identity and the linearity of Omega all hold
%%% whatever the constants and the units happen to be. None of them would
%%% have fired on the 2026-09-18 campaign, which applied 200 t/ha of calcite,
%%% dissolved 17 and reported a pore-water pH of 6.1. The Energy Farm has no
%%% calcite treatment, so no observation constrains that run either.
%%%
%%% This test supplies the missing constraint from thermodynamics rather than
%%% from the field. Calcite in water open to a fixed CO2 partial pressure is a
%%% closed-form system with one degree of freedom, solved here independently
%%% of the model and compared with the number any geochemistry text or
%%% speciation code gives. Two things follow from it and are asserted:
%%%
%%%   1. the model's constant set and unit conventions reproduce the textbook
%%%      state at atmospheric CO2 (pH 8.3, Ca 0.5 mM), which a unit slip in
%%%      Henry's law, a concentration/activity confusion or a wrong gamma_IAP
%%%      all break by 0.3 pH units or more;
%%%
%%%   2. soil CO2, not the applied dose, sets the pH a calcite-amended soil
%%%      can reach. The ceiling ladder below is the diagnostic the campaign
%%%      lacked: a calcite run sitting at pore pH 6.1 is reporting a soil
%%%      atmosphere near 1 atm CO2, not an exhausted amendment.
%%%
%%% Conventions follow test_erw_omega_activity_basis and ERW_function: the
%%% solver's species are CONCENTRATIONS formed with ionic-strength-
%%% conditional constants (ERW_Conditional_Constants), so its CO3 is a
%%% concentration; Omega = gamma_IAP*Ca*CO3/K_CaCO3 with gamma_IAP =
%%% gamma_div^2 and the thermodynamic K_CaCO3, and the rate law receives
%%% the proton activity gamma_mono*H.

test_dir = fileparts(mfilename('fullpath'));
code_dir = fullfile(test_dir,'..');
addpath(code_dir);

conv_mol = 1e6;
V_m      = 22.41;            % L/mol, the molar volume ERW_Calcite_Kinetics uses
T        = 298.15;
[K1,K2,Kw] = carbonate_constants(T);
Ksp      = 10^-8.35;         % ERW_MOD_PARAM_excel.m
k_H      = 0.83;             % dimensionless CO2 solubility at the reference T

%% --- 1. the constant set itself, against the published values at 25 degC
assert(abs(-log10(K1) - 6.35) < 0.01, 'pK1 must be 6.35 at 25 degC; got %.4f.',-log10(K1));
assert(abs(-log10(K2) - 10.33) < 0.01, 'pK2 must be 10.33 at 25 degC; got %.4f.',-log10(K2));
assert(abs(-log10(Kw) - 14.00) < 0.01, 'pKw must be 14.00 at 25 degC; got %.4f.',-log10(Kw));
% k_H is dimensionless on a GAS CONCENTRATION; the Henry constant it implies
% is k_H/V_m [mol L-1 atm-1] and must sit near the published 0.0339.
K_H_implied = k_H/V_m;
assert(abs(K_H_implied/0.0339 - 1) < 0.12, ...
    'Implied Henry constant %.5f mol/L/atm is not within 12 percent of 0.0339.',K_H_implied);

%% --- 2. atmospheric anchor: calcite + water + air
[pH_atm, Ca_atm] = calcite_open_equilibrium(400e-6);
assert(abs(pH_atm - 8.30) < 0.25, ...
    ['Calcite-saturated water at 400 ppm CO2 must give pH 8.3 +/- 0.25; ' ...
     'got %.2f. A value near 7.3 is the Henry unit slip (k_H used as ' ...
     'mol/L/atm); a value near 6.5 is a gamma or Ksp error.'], pH_atm);
assert(Ca_atm > 0.40e-3 && Ca_atm < 0.75e-3, ...
    'Calcite-saturated water at 400 ppm CO2 must give Ca 0.4-0.75 mM; got %.3f mM.',Ca_atm*1e3);

%% --- 3. the production Omega path agrees with the anchor
% Feed the anchored state through ERW_Calcite_Kinetics in the solver's own
% units and require Omega == 1 and a rate of exactly zero, then require a
% strictly positive rate once the state is undersaturated.
pCO2 = 0.03;                                   % 3 % of soil air
[pH_eq, Ca_eq, I_eq] = calcite_open_equilibrium(pCO2);
[Ca_u, CO3_u, H_u, CO2_air_u, Ksp_u, gamma_IAP] = ...
    model_units(pH_eq, Ca_eq, I_eq, pCO2, K1, K2, k_H, Ksp, conv_mol, V_m);

rate_args = {0.7, 1, 1e-4, 1e-6, 1e-5, 0.5, 0.3, 2.0, 0.1, 5000, 100.087, 1};
[Omega_eq, Wr_eq] = ERW_Calcite_Kinetics(Ca_u, CO3_u, Ksp_u, H_u, CO2_air_u, ...
    conv_mol, rate_args{:}, gamma_IAP);
assert(abs(Omega_eq - 1) < 1e-6, ...
    'The production Omega must be 1 at calcite equilibrium; got %.9f.',Omega_eq);

[Omega_lo, Wr_lo] = ERW_Calcite_Kinetics(0.5*Ca_u, CO3_u, Ksp_u, H_u, CO2_air_u, ...
    conv_mol, rate_args{:}, gamma_IAP);
assert(abs(Omega_lo - 0.5) < 1e-6 && Wr_lo > 0, ...
    'Halving Ca must halve Omega and open a positive dissolution rate.');

% The equilibrium solve lands on Omega = 1 to the tolerance above, not bit
% exactly, so the affinity term is a rounding residue rather than a hard
% zero. It must be negligible against the rate a genuinely undersaturated
% state carries; 1e-5 of it is four orders below the tolerance on Omega.
assert(Wr_eq >= 0 && Wr_eq < 1e-5*Wr_lo, ...
    'The affinity term must vanish at Omega = 1; got %g against %g.',Wr_eq,Wr_lo);

%% --- 4. the ceiling ladder: soil CO2 sets the attainable pH
% Reference values from the independent solve above; they are a property of
% the constants, not a fit. Monotone decreasing, and by 70 % CO2 the ceiling
% has fallen BELOW the untreated Energy Farm soil's own pH of about 6.1 --
% at which point no amount of calcite can raise the modelled pH at all.
ladder_pct = [0.04 1 3 10 30 70 100];
ladder_ref = [8.23 7.31 7.00 6.66 6.35 6.12 6.02];
ceiling = arrayfun(@(p) calcite_open_equilibrium(p/100), ladder_pct);
for i = 1:numel(ladder_pct)
    assert(abs(ceiling(i) - ladder_ref(i)) < 0.08, ...
        'Calcite-equilibrium pH at %.2f percent CO2 must be %.2f +/- 0.08; got %.2f.', ...
        ladder_pct(i), ladder_ref(i), ceiling(i));
end
assert(all(diff(ceiling) < 0), 'The pH ceiling must fall monotonically with soil CO2.');
assert(ceiling(end-1) < 6.2, ...
    ['At 70 percent soil CO2 the calcite-equilibrium pH must fall below 6.2, ' ...
     'i.e. below the untreated soil. The 2026-09-18 campaign ran at 70-109 percent.']);

%% --- 5. the ceiling is a CO2 property, not a dose property
% Omega, and therefore the affinity term, does not see the remaining rock at
% all: once the solution is at saturation the fifth dose behaves like the
% first. This is why 200 t/ha and 50 t/ha give the same equilibrium pH and
% why a dose increase is not a test of the amendment.
for M_calcite = [5000 20000 200000]
    args = rate_args; args{10} = M_calcite;
    [Omega_M, Wr_M] = ERW_Calcite_Kinetics(Ca_u, CO3_u, Ksp_u, H_u, CO2_air_u, ...
        conv_mol, args{:}, gamma_IAP);
    assert(Omega_M == Omega_eq && Wr_M == Wr_eq, ...
        'Omega and the affinity term must not depend on the calcite inventory.');
end

fprintf('PASS: calcite saturation anchored to the open-system equilibrium.\n');

%% ------------------------------------------------------------------ helpers
    function [K1_,K2_,Kw_] = carbonate_constants(T_K)
        % The temperature functions of ERW_MOD_PARAM_excel.m, restated here so
        % the anchor does not depend on the file it is checking.
        pk1_ = -(-356.309-0.0609*T_K+21834.37/T_K+126.8339*log10(T_K)-1684915/T_K^2);
        pk2_ = -(-107.887-0.032528*T_K+5151.79/T_K+38.92561*log10(T_K)-563713.9/T_K^2);
        pkw_ = -(-283.971+13323/T_K-0.0507*T_K+102.24447*log10(T_K)-1119669/T_K^2);
        K1_ = 10^-pk1_; K2_ = 10^-pk2_; Kw_ = 10^-pkw_;
    end

    function [pH_out, Ca_out, I_out] = calcite_open_equilibrium(pCO2_atm)
        % Calcite-saturated water open to a fixed pCO2, on the activity basis,
        % with Davies iterated to a self-consistent ionic strength. Independent
        % of ERW_function; uses only ERW_Activity_Coefficient.
        I_out = 1e-3; converged = false;
        for it = 1:300
            g1 = ERW_Activity_Coefficient(1,I_out);
            g2 = ERW_Activity_Coefficient(2,I_out);
            pH_out = fzero(@(pH) charge_residual(pH,pCO2_atm,g1,g2), [3 12]);
            [Ca_out, I_new] = state_at(pH_out, pCO2_atm, g1, g2);
            if abs(I_new - I_out) < 1e-12
                I_out = I_new; converged = true; break
            end
            I_out = 0.5*(I_out + I_new);
        end
        assert(converged, 'Calcite equilibrium did not converge at pCO2 = %g atm.',pCO2_atm);
    end

    function r = charge_residual(pH_t, pCO2_atm, g1, g2)
        [aH,aHCO3,aCO3,aOH] = speciate(pH_t, pCO2_atm);
        cCa = (Ksp/aCO3)/g2;
        r = 2*cCa + aH/g1 - (aHCO3/g1 + 2*aCO3/g2 + aOH/g1);
    end

    function [cCa, I_new] = state_at(pH_t, pCO2_atm, g1, g2)
        [aH,aHCO3,aCO3,aOH] = speciate(pH_t, pCO2_atm);
        cCa = (Ksp/aCO3)/g2;
        I_new = 0.5*(4*cCa + aH/g1 + aHCO3/g1 + 4*aCO3/g2 + aOH/g1);
    end

    function [aH,aHCO3,aCO3,aOH] = speciate(pH_t, pCO2_atm)
        aH = 10^-pH_t;
        aCO2 = k_H*pCO2_atm/V_m;      % neutral: activity == concentration
        aHCO3 = K1*aCO2/aH;
        aCO3 = K2*aHCO3/aH;
        aOH = Kw/aH;
    end

    function [Ca_u_,CO3_u_,H_u_,CO2_air_u_,Ksp_u_,gIAP_] = model_units( ...
            pH_e, Ca_e, I_e, pCO2_atm, K1_, K2_, kH_, Ksp_, conv_, Vm_)
        g1 = ERW_Activity_Coefficient(1,I_e);
        g2 = ERW_Activity_Coefficient(2,I_e);
        aH = 10^-pH_e;
        aCO2 = kH_*pCO2_atm/Vm_;
        aCO3 = K2_*K1_*aCO2/aH^2;
        % The solver's CO3, formed with the conditional constants, is the
        % concentration a_CO3/gamma_div, which the production helper must
        % reproduce; ERW_function hands the rate law gamma_mono*H, i.e.
        % a_H*conv_mol.
        kc = ERW_Conditional_Constants( ...
            [K1_*conv_,K2_*conv_,1e-14*conv_^2,1e-5*conv_,1e-5*conv_,1e-7*conv_,1e-6*conv_],[],I_e);
        CO3_c = kc(1)*kc(2)*(aCO2*conv_)/((aH/g1)*conv_)^2;
        assert(abs(CO3_c/((aCO3/g2)*conv_) - 1) < 1e-12, ...
            'Conditional speciation must give the carbonate concentration a_CO3/gamma_div.');
        H_u_       = aH*conv_;
        Ca_u_      = Ca_e*conv_;
        CO3_u_     = CO3_c;
        CO2_air_u_ = pCO2_atm/Vm_*conv_;
        Ksp_u_     = Ksp_*conv_^2;
        gIAP_      = g2^2;
    end
end
