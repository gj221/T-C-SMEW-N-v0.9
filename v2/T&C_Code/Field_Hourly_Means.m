function m = Field_Hourly_Means(h,w)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Field means of one hour's fluxes, for the tower comparisons (package y).
%%
%% h holds the hour's values of every square, one column per square:
%%   1 x N   Rn H QE G Ts EG ELitter ESN ESN_In EWAT EICE EIn_urb EIn_rock
%%           Lk Rh Rd
%%   cc x N  T_H T_L EIn_H EIn_L An_H An_L Rdark_H Rdark_L Ccrown
%%   nl x N  O Tdp (the active-layer soil layers)
%% w is 1 x N, the weights Domain uses (each square's share of the field).
%%
%%   ET   as the 1-D frame defines it: transpiration and interception from
%%        both vegetation layers plus soil, litter, snow, water, ice, urban
%%        and rock evaporation [mm h-1]
%%   T    transpiration T_H + T_L [mm h-1]
%%   GPP  sum over PFT slots of Ccrown x (An + Rdark), both layers
%%        [umol CO2 m-2 s-1]; An and Rdark are per unit crown area
%%   Runoff  Rh + Rd [mm per step]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fm = @(x) sum(x.*w,2).';
ET = sum(h.T_H+h.EIn_H,1) + sum(h.T_L+h.EIn_L,1) + h.EG + h.ELitter + ...
    h.ESN + h.ESN_In + h.EWAT + h.EICE + h.EIn_urb + h.EIn_rock;
m.Rn = fm(h.Rn);
m.H = fm(h.H);
m.QE = fm(h.QE);
m.G = fm(h.G);
m.Ts = fm(h.Ts);
m.ET = fm(ET);
m.T = fm(sum(h.T_H,1) + sum(h.T_L,1));
m.EG = fm(h.EG);
m.GPP = fm(sum(h.Ccrown.*(h.An_H + h.Rdark_H + h.An_L + h.Rdark_L),1));
m.O = fm(h.O);
m.Tdp = fm(h.Tdp);
m.Lk = fm(h.Lk);
m.Runoff = fm(h.Rh + h.Rd);
end
