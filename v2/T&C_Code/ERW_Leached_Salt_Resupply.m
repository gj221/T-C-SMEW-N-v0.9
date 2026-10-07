function [I_Ca,I_Mg] = ERW_Leached_Salt_Resupply(f_Lw,An_tot,Na_tot,Ca,Mg)
%%% Return the leached residual anion as the salt it left as.
%%%
%%% The legacy steady closure (keyword_add = 1) tops the residual anion and
%%% Na inventories back up each day by what leached or was taken up, as
%%% bare ions. That is alkalinity-neutral on its own: the ERW books its own
%%% An and Na leaching and uptake, and the re-supply cancels them. The
%%% cations that balanced the anion in the leachate are T&C's Ca and Mg,
%%% and while their leaching coefficient was 0.4 the leachate carried less
%%% cation than anion charge, so nothing left with the anion. With
%%% dissolved Ca and Mg leaving at their concentration the leachate is
%%% neutral: the anion charge Na did not carry out left on Ca and Mg, and
%%% returning the anion alone would be an acid input of about 0.65 eq m-2
%%% yr-1 in the control (budget of package m, 1 October 2026). The leached
%%% part of the closure therefore comes back with those cations, split by
%%% their pore-water charge share, through the I_Ca and I_Mg boundary slots
%%% ERW_function books in dAlk_boundary. The uptake part keeps the legacy
%%% booking, which has no cation counterpart in the ERW.
%%%
%%% f_Lw        fraction of the active-layer water drained today   [d-1]
%%% An_tot      residual anion inventory                           [umol charge m-2]
%%% Na_tot      sodium inventory                                   [umol m-2]
%%% Ca, Mg      pore-water concentrations                          [umol L-1]
%%% I_Ca, I_Mg  boundary supply, as ERW_function expects them      [umol m-2 d-1]
x = [f_Lw An_tot Na_tot Ca Mg];
if any(~isfinite(x)) || any(x < 0)
    error('ERW:invalidParam','Leached-salt inputs must be finite and nonnegative.');
end
charge = f_Lw*max(0,An_tot-Na_tot);
if Ca+Mg > 0
    share_Ca = Ca/(Ca+Mg);
else
    share_Ca = 0.5;
end
I_Ca = share_Ca*charge/2;
I_Mg = (1-share_Ca)*charge/2;
end
