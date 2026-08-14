function W = ERW_Window_Interp(dtw, day_idx, res_per_h, ...
    CO2_h_day, s_pair, temp_pair, rhet_pair, raut_pair, Tw_day, Lw_day, Iw_day)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ERW_Window_Interp
%%
%% Build the ERW sub-hourly forcing window for ONE day, without building a
%% timetable over the whole run.
%%
%% WHY THIS EXISTS
%%   The single-column MAIN_FRAME re-runs synchronize() over the entire
%%   simulation on every simulated day, then uses only that day's rows.
%%   For a 5-year run at 10-minute resolution that is a ~263,000-row
%%   timetable rebuilt 1825 times.  On one column it is merely wasteful;
%%   with one column per sub-plot cell it becomes the dominant cost and
%%   makes the multi-column run intractable.  This routine reproduces the
%%   same interpolation for the current day only.
%%
%% EQUIVALENCE TO THE TIMETABLE PATH
%%   MAIN_FRAME synchronizes with:
%%     'linear'   for CO2_atm (hourly) and for s, temp_soil, r_het, r_aut
%%                (daily)
%%     'previous' for Tw, Lw, Iw (hourly values retimed to daily sums)
%%   The same two rules are applied here: linear ramps across the day for
%%   the daily-resolution drivers, linear within the hour for CO2, and a
%%   zero-order hold for the daily water fluxes.
%%
%% Inputs
%%   dtw        ERW timestep [day] (= 1/(24*res_per_h))
%%   day_idx    index of the day being solved (used only for reporting)
%%   res_per_h  ERW sub-steps per hour (MOD_PARAM t_res; 6 -> 10 minutes)
%%   CO2_h_day  [25 x 1] hourly atmospheric CO2 spanning this day plus the
%%              first hour of the next [umol-conv/L]
%%   s_pair     [2 x 1] soil saturation at the start and end of the day [-]
%%   temp_pair  [2 x 1] soil temperature, start and end of day [degC]
%%   rhet_pair  [2 x 1] heterotrophic respiration, start and end [umol/m2 d]
%%   raut_pair  [2 x 1] autotrophic root respiration, start and end
%%   Tw_day     scalar daily transpiration within Zbio [m/d]
%%   Lw_day     scalar daily leakage [m/d]
%%   Iw_day     scalar daily infiltration [m/d]
%%
%% Output
%%   W  struct of [n_sub x 1] columns: CO2_atm, s, temp_soil, r_het,
%%      r_aut, Tw, Lw, Iw  -- the same fields the ERW solver reads off
%%      TT_m, in the same units.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_sub = round(1/dtw);              %% sub-steps in one day (24*res_per_h)
tau   = ((1:n_sub)' - 1) / n_sub;  %% fractional position through the day

%% Daily-resolution drivers: linear ramp from the day's start to its end.
W.s         = s_pair(1)    + tau .* (s_pair(2)    - s_pair(1));
W.temp_soil = temp_pair(1) + tau .* (temp_pair(2) - temp_pair(1));
W.r_het     = rhet_pair(1) + tau .* (rhet_pair(2) - rhet_pair(1));
W.r_aut     = raut_pair(1) + tau .* (raut_pair(2) - raut_pair(1));

%% CO2 is hourly: linear within each hour of the day.
h_pos   = tau * 24;                       %% hours into the day, [0,24)
h_lo    = floor(h_pos);
frac    = h_pos - h_lo;
idx_lo  = min(h_lo + 1, numel(CO2_h_day) - 1);
W.CO2_atm = CO2_h_day(idx_lo) + frac .* (CO2_h_day(idx_lo + 1) - CO2_h_day(idx_lo));

%% Daily water fluxes: zero-order hold, matching 'previous'.
W.Tw = Tw_day * ones(n_sub, 1);
W.Lw = Lw_day * ones(n_sub, 1);
W.Iw = Iw_day * ones(n_sub, 1);

W.n_sub   = n_sub;
W.day_idx = day_idx;

end
