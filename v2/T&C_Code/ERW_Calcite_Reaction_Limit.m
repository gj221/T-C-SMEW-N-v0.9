function [primary_rate,secondary_rate,diag] = ERW_Calcite_Reaction_Limit( ...
    primary_rate,secondary_rate,dtw,Omega0,Ca,DIC,Alk,water_l, ...
    k1,k2,k_w,K_CaCO3,conv_mol,gamma_IAP)
%%% Bound the combined calcite extent at closed-water equilibrium on the
%%% activity basis: the trial saturation state is gamma_IAP*Ca*CO3/K_CaCO3
%%% with gamma_IAP = gamma_div^2 the factor the caller formed
%%% Omega0 with, held at the step's incoming value through the bisection.
%%% It is required: the concentration basis it would default to is the
%%% defect this argument exists to remove.

if nargin < 14
    error('ERW:invalidCalciteLimit', ...
        'gamma_IAP is required; pass gamma_div^2 at the step ionic strength.');
end
if ~(isscalar(gamma_IAP) && isfinite(gamma_IAP) && gamma_IAP > 0 && gamma_IAP <= 1)
    error('ERW:invalidCalciteLimit', ...
        'Calcite activity-product factor must lie in (0,1]; got %g.',gamma_IAP);
end
values = [primary_rate,secondary_rate,dtw,Omega0,Ca,DIC,Alk,water_l, ...
    k1,k2,k_w,K_CaCO3,conv_mol];
if any(~isfinite(values)) || dtw <= 0 || water_l <= 0 || ...
        K_CaCO3 <= 0 || conv_mol <= 0 || Ca < 0 || DIC < 0
    error('ERW:invalidCalciteLimit','Calcite reaction limiter received an invalid state.');
end
if primary_rate < 0
    error('ERW:calcitePathConflict', ...
        'Applied calcite cannot receive secondary precipitation.');
end
if primary_rate > 0 && secondary_rate < 0
    error('ERW:calcitePathConflict', ...
        'Primary dissolution and secondary precipitation cannot occur together.');
end

raw_rate = primary_rate + secondary_rate;
raw_extent = raw_rate*dtw;
extent = raw_extent;
limited = false;

if raw_extent < 0
    inventory = min(Ca,DIC)*water_l;
    if -extent > inventory
        extent = -inventory;
        limited = true;
    end
end

Omega_trial = trial_omega(extent);
g0 = Omega0-1;
g1 = Omega_trial-1;
if raw_extent == 0 || abs(g0) <= 1e-12
    extent = 0;
    limited = raw_extent ~= 0;
elseif ~isfinite(Omega_trial)
    error('ERW:invalidCalciteLimit','Calcite trial state has no carbonate root.');
elseif (raw_extent > 0 && g0 > 0) || (raw_extent < 0 && g0 < 0)
    error('ERW:calciteDirection','Calcite reaction points away from saturation.');
elseif sign(g0) ~= sign(g1) && g1 ~= 0
    lo = 0;
    hi = extent/raw_extent;
    for it = 1:60
        mid = 0.5*(lo+hi);
        g_mid = trial_omega(mid*raw_extent)-1;
        if ~isfinite(g_mid)
            hi = mid;
        elseif sign(g_mid) == sign(g0)
            lo = mid;
        else
            hi = mid;
        end
    end
    extent = lo*raw_extent;
    limited = true;
end

if raw_extent == 0
    scale = 0;
else
    scale = extent/raw_extent;
end
primary_rate = primary_rate*scale;
secondary_rate = secondary_rate*scale;
[Omega_final,Ca_final,DIC_final,Alk_final] = trial_omega(extent);
diag = struct('limited',limited,'extent',extent,'Omega_final',Omega_final, ...
    'Ca_final',Ca_final,'DIC_final',DIC_final,'Alk_final',Alk_final);

    function [Omega,Ca_out,DIC_out,Alk_out] = trial_omega(xi)
        Ca_out = Ca+xi/water_l;
        DIC_out = DIC+xi/water_l;
        Alk_out = Alk+2*xi/water_l;
        %% The inventory limit sets the extent to -min(Ca,DIC)*water_l, and
        %% dividing it back by water_l can leave the exhausted pool a few
        %% ulps below zero. That is an empty pool, not a negative one: clamp
        %% it so the carbonate root of the exhausted state (Omega = 0) is
        %% found instead of reporting that no root exists. The tolerance is
        %% the larger of the historical absolute bound and eight ulps of the
        %% pool, so it stays a rounding allowance at any concentration.
        rounding = max(1e-10,8*eps(max(Ca,DIC)));
        if Ca_out < -rounding || DIC_out < -rounding
            Omega = NaN;
            return
        end
        Ca_out = max(Ca_out,0);
        DIC_out = max(DIC_out,0);
        [H_trial,CO3_trial] = carbonate_state(DIC_out,Alk_out);
        Omega = gamma_IAP*Ca_out*CO3_trial/K_CaCO3;
        if ~isfinite(H_trial)
            Omega = NaN;
        end
    end

    function [H_trial,CO3_trial] = carbonate_state(DIC_in,Alk_in)
        pH_lo = 1e-8;
        pH_hi = 14-1e-8;
        r_lo = carbonate_residual(pH_lo,DIC_in,Alk_in);
        r_hi = carbonate_residual(pH_hi,DIC_in,Alk_in);
        if ~isfinite(r_lo) || ~isfinite(r_hi) || sign(r_lo) == sign(r_hi)
            H_trial = NaN;
            CO3_trial = NaN;
            return
        end
        for jt = 1:60
            pH_mid = 0.5*(pH_lo+pH_hi);
            r_mid = carbonate_residual(pH_mid,DIC_in,Alk_in);
            if sign(r_mid) == sign(r_lo)
                pH_lo = pH_mid;
                r_lo = r_mid;
            else
                pH_hi = pH_mid;
            end
        end
        H_trial = 10^(-0.5*(pH_lo+pH_hi))*conv_mol;
        den = H_trial^2+k1*H_trial+k1*k2;
        CO3_trial = DIC_in*k1*k2/den;
    end

    function residual = carbonate_residual(pH_trial,DIC_in,Alk_in)
        H_trial = 10^(-pH_trial)*conv_mol;
        den = H_trial^2+k1*H_trial+k1*k2;
        HCO3_trial = DIC_in*k1*H_trial/den;
        CO3_trial = DIC_in*k1*k2/den;
        residual = HCO3_trial+2*CO3_trial-H_trial+k_w/H_trial-Alk_in;
    end
end
