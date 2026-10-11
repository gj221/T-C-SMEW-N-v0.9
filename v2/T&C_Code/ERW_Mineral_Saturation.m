function [Omega, thermodynamic_defined] = ERW_Mineral_Saturation( ...
    mineral, K, Ca, Mg, Na, Si, Al, PO4, H, K_sp, conv_mol)
%%% Primary-mineral saturation ratio from the complete modeled reaction.
%%%
%%% Dissolved concentrations are in umol/L and gamma_i=1. K_sp is the
%%% proton-normalized equilibrium constant for forsterite, wollastonite,
%%% diopside, albite and anorthite. Thermoddem reports the zoisite constant
%%% for the complete 13-H reaction, so that constant is rooted here.
%%%
%%% A thermodynamic ratio is unavailable when a reaction requires a
%%% dissolved product that the model does not represent. Those phases use
%%% far-from-equilibrium dissolution kinetics in the callers; NaN prevents
%%% the placeholder from being reported as a saturation index.

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

ca = Ca/conv_mol;
mg = Mg/conv_mol;
na = Na/conv_mol;
si = Si/conv_mol;
al = Al/conv_mol;
h = H/conv_mol;

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
