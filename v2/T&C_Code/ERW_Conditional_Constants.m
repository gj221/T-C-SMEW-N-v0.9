function [k_cond,logKso_cond,gamma_z] = ERW_Conditional_Constants(k,logKso,I)
%%% Ionic-strength-conditional (concentration-quotient) forms of the ERW
%%% mass-action constants, from Davies single-ion activity coefficients
%%% gamma_z = [g1 g2 g3] at ionic strength I [mol/L]
%%% (ERW_Activity_Coefficient).
%%%
%%%   k       [k1 k2 k_w K1 K2 K3 K4], thermodynamic, on the model's umol/L
%%%           scale (k_w in (umol/L)^2)
%%%   logKso  log10 *Kso of Al(OH)3 + 3H+ = Al+++ + 3H2O, or [] when the
%%%           caller has no Al(OH)3 phase (logKso_cond is then [])
%%%
%%% With a_i = g(z_i)*c_i and neutral species at unit activity, each
%%% constant applied to concentrations becomes
%%%     CO2 + H2O = H+ + HCO3-           k1c   = k1/g1^2
%%%     HCO3- = H+ + CO3--               k2c   = k2/g2
%%%     H2O = H+ + OH-                   k_wc  = k_w/g1^2
%%%     Al+++ + H2O = AlOH++ + H+        K1c   = K1*g3/(g1*g2)
%%%     AlOH++ + H2O = Al(OH)2+ + H+     K2c   = K2*g2/g1^2
%%%     Al(OH)2+ + H2O = Al(OH)3 + H+    K3c   = K3
%%%     Al(OH)3 + H2O = Al(OH)4- + H+    K4c   = K4/g1^2
%%%     [Al+++]/[H+]^3 = *Kso            logKso_c = logKso + log10(g1^3/g3)
%%% Every factor is a ratio of dimensionless coefficients, so the conv_mol
%%% powers k carries are unchanged, and at I = 0 every factor is exactly 1.
%%% Molarity stands in for molality (dilute soil water).  NH4+ = NH3 + H+ is
%%% charge-symmetric and CO2 is neutral, so k_nh3 and k_H need no factor.

if ~(isnumeric(k) && numel(k) == 7 && all(isfinite(k(:))) && all(k(:) > 0))
    error('ERW:invalidParam', ...
        'k must hold seven positive finite constants [k1 k2 k_w K1 K2 K3 K4].');
end
if ~(isnumeric(logKso) && (isempty(logKso) || (isscalar(logKso) && isfinite(logKso))))
    error('ERW:invalidParam','logKso must be a finite scalar or empty.');
end
if ~(isnumeric(I) && isscalar(I))
    error('ERW:invalidIonicStrength','Ionic strength must be a numeric scalar.');
end

gamma_z = ERW_Activity_Coefficient(1:3,I);
g1 = gamma_z(1);
g2 = gamma_z(2);
g3 = gamma_z(3);
k_cond = reshape(k,1,7).*[1/g1^2, 1/g2, 1/g1^2, g3/(g1*g2), g2/g1^2, 1, 1/g1^2];
logKso_cond = logKso + log10(g1^3/g3);   % stays [] for an empty logKso
end
