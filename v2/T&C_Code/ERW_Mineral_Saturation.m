function [Omega, thermodynamic_defined] = ERW_Mineral_Saturation( ...
    mineral, K, Ca, Mg, Na, Si, Al, PO4, H, K_sp, conv_mol, Al_CEC_scale, I_str)
%%% Primary-mineral saturation ratio from the complete modeled reaction.
%%%
%%% Dissolved concentrations are in umol/L except Al, which arrives on the
%%% ERW solver's Al_CEC_scale basis (conv_Al, nmol/L) and is divided by
%%% Al_CEC_scale here; a caller that already holds Al in umol/L passes 1.
%%% Activities are Davies single-ion activity coefficients at the pore-water
%%% ionic strength I_str [mol/L] (ERW_Activity_Coefficient). Both are
%%% required: the values they would default to, 1 and 0, are the two
%%% defects these arguments exist to remove -- an Al scale read as umol/L
%%% and a concentration-basis saturation ratio.
%%% Every reaction quotient below is already the full quotient raised to
%%% 1/n_H, so substituting gamma*c for each ion is the exact activity
%%% form; Si (H4SiO4, neutral) takes no correction. K_sp is the
%%% proton-normalized equilibrium constant for forsterite, wollastonite,
%%% diopside, albite and anorthite. Thermoddem reports the zoisite constant
%%% for the complete 13-H reaction, so that constant is rooted here.
%%%
%%% A thermodynamic ratio is unavailable when a reaction requires a
%%% dissolved product that the model does not represent. Those phases use
%%% far-from-equilibrium dissolution kinetics in the callers; NaN prevents
%%% the placeholder from being reported as a saturation index.
%%%
%%% MINERALOGY CAVEAT: the ratio returned here is why albite and quartz
%%% never dissolve at Kantola (no silica sink, so the pore water becomes
%%% supersaturated with both and their affinity max(0, 1 - Omega) is zero),
%%% and the NaN returned for epidote, chlorite and the amphibole is why
%%% those dissolve at the far-from-equilibrium rate. The plain-language
%%% account, and what it means for Ca+Mg release and CDR, is the block
%%% headed MINERALOGY CAVEAT in ERW_function.m, above the call to this
%%% function.

if iscell(mineral)
    if numel(mineral) ~= 1
        error('ERW:invalidParam','Mineral name must be scalar.');
    end
    mineral = mineral{1};
end
mineral = lower(strtrim(char(string(mineral))));

unsupported = {'labradorite','augite','apatite','chlorite','epidote', ...
    'ferroactinolite','actinolite','inert'};
if any(strcmp(mineral,unsupported))
    Omega = NaN;
    thermodynamic_defined = false;
    return
end

if ~(isfinite(conv_mol) && conv_mol > 0) || ~(isfinite(H) && H > 0)
    error('ERW:invalidState','Mineral saturation requires positive H and conv_mol.');
end
dissolved = [K,Ca,Mg,Na,Si,Al,PO4];
if any(~isfinite(dissolved)) || any(dissolved < 0)
    error('ERW:invalidState', ...
        'Mineral saturation requires finite nonnegative dissolved concentrations.');
end
if ~(isfinite(K_sp) && K_sp > 0)
    error('ERW:invalidParam','Mineral equilibrium constant must be positive and finite.');
end
if nargin < 13
    error('ERW:invalidParam', ...
        ['Al_CEC_scale and I_str are required; pass the solver''s conv_Al and ' ...
        'the step''s pore-water ionic strength.']);
end
if ~(isscalar(Al_CEC_scale) && isfinite(Al_CEC_scale) && Al_CEC_scale > 0)
    error('ERW:invalidParam','Al_CEC_scale must be a positive finite scalar.');
end
if ~(isscalar(I_str) && isfinite(I_str) && I_str >= 0)
    error('ERW:invalidState','Ionic strength must be a finite nonnegative scalar.');
end

gamma1 = ERW_Activity_Coefficient(1,I_str);
gamma2 = ERW_Activity_Coefficient(2,I_str);
gamma3 = ERW_Activity_Coefficient(3,I_str);
ca = gamma2*Ca/conv_mol;
mg = gamma2*Mg/conv_mol;
na = gamma1*Na/conv_mol;
si = Si/conv_mol;
al = gamma3*Al/(conv_mol*Al_CEC_scale);
h = gamma1*H/conv_mol;

switch mineral
    case 'forsterite'
        Omega = mg^(2/4)*si^(1/4)/(h*K_sp);
    case 'wollastonite'
        Omega = ca^(1/2)*si^(1/2)/(h*K_sp);
    case 'diopside'
        Omega = ca^(1/4)*mg^(1/4)*si^(2/4)/(h*K_sp);
    case 'albite'
        Omega = na^(1/4)*al^(1/4)*si^(3/4)/(h*K_sp);
    case 'anorthite'
        Omega = ca^(1/8)*al^(2/8)*si^(2/8)/(h*K_sp);
    case 'zoisite'
        Omega = ca^(2/13)*al^(3/13)*si^(3/13)/(h*K_sp^(1/13));
    case 'quartz'
        Omega = si/K_sp;
    otherwise
        error('ERW:invalidParam','No mineral reaction is defined for %s.',mineral);
end

if ~(isfinite(Omega) && Omega >= 0)
    error('ERW:invalidState','Mineral saturation ratio is not finite and nonnegative.');
end
thermodynamic_defined = true;
end
