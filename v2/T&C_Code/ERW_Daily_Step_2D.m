function E = ERW_Daily_Step_2D(E, cfg, p, i, j, NN, NNd, Datam, ~, Ca, ...
    O_c, Lk_c, f_c, T_H_c, T_L_c, Tdp_bio, VSUM, VTSUM, ...
    R_microbe_c, R_litter_c, R_ew_c, R_litter_sur_c, Rmr_H_c, Rmr_L_c, Ccrown, ...
    P_c, MOD_PARAM_input, ERW_PARAM_input, conv_mol, conv_Al, N_col)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ERW_Daily_Step_2D
%%
%% One day of the enhanced-rock-weathering solver for ONE sub-plot column.
%% Holds that column's ERW state in E between calls, so columns weather
%% independently and can receive different amendments.
%%
%% TREATMENT PLACEMENT
%%   cfg.apply_rock decides whether this column receives the amendment.
%%   A column with apply_rock = false still runs the full carbonate and
%%   cation-exchange chemistry -- it simply gets no rock added, which is
%%   what makes an untreated field interior comparable with a treated
%%   margin inside one simulation.
%%
%% FIDELITY TO THE SINGLE-COLUMN FRAME
%%   The 1-D frame rebuilds synchronised timetables spanning the whole
%%   run on every simulated day and then reads one day out of them:
%%   ~263,000 rows rebuilt 1825 times for a 5-year run, and once per
%%   column here.  This version reproduces those timetables exactly
%%   without rebuilding them:
%%
%%     - the regular 10-minute time grid and the CO2 column depend only
%%       on forcing that is known upfront, so they are built once at
%%       initialisation and cached;
%%     - the daily drivers (s, temp_soil, r_het, r_aut) are linearly
%%       interpolated onto the grid with interp1, which is what
%%       synchronize(...,'linear') does, evaluated only on the rows the
%%       current day actually reads;
%%     - the daily water sums (Tw, Lw, Iw) are held zero-order from the
%%       calendar day containing each row, which is what
%%       synchronize(...,'previous') does over retime(...,'daily','sum').
%%
%%   An earlier attempt at this shortcut silently changed the answers:
%%   infiltration was taken as the current hour scaled by 24 instead of
%%   summed over the day, and the day-boundary sub-step alignment was
%%   off by one.  Both are fixed here, and
%%   test_2d_driver_equivalence checks the result bit for bit against
%%   the 1-D frame rather than taking the equivalence on trust.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

E.has_update = false;

%% ---------------------------------------------------------------------
%% One-off initialisation for this column
%% ---------------------------------------------------------------------
P_prev_for_erw_delta = P_c(j-1,:);
if ~E.init
    E.opt_TC_biochem = ERW_PARAM_input.opt_TC_biochem;
    E.dtd = 1;
    E.t_end = NN/24;
    E.dtw = 1/(24*MOD_PARAM_input.t_res);
    E.t = E.dtw:E.dtw:E.t_end;

    C.n        = sum(p.Bio_Zs.*p.Osat);
    C.rsd_bio  = sum(p.Bio_Zs.*p.rsd);
    C.PH_in    = p.PH;
    C.Zbio     = p.Zbio;
    C.ms       = p.ms;
    C.conv_mol = conv_mol;
    C.conv_Al  = conv_Al;
    if isfield(p,'Porg');   C.Porg  = p.Porg;   end
    if isfield(p,'Pcla');   C.Pcla  = p.Pcla;   end
    if isfield(p,'rain_pH_vwm'); C.rain_pH_vwm = p.rain_pH_vwm; end
    if isfield(p,'bulk_pH_params'); C.bulk_pH_params = p.bulk_pH_params; end
    E.TC_const = C;

    %% Per-column series, filled in as the run proceeds
    E.TC.s         = zeros(NNd,1);
    E.TC.temp_soil = zeros(NNd,1);
    E.TC.r_het     = zeros(NNd,1);
    E.TC.r_aut     = zeros(NNd,1);
    E.TC.f_Lw_bio  = zeros(NNd,1);
    E.TC.f_T_bio   = zeros(NNd,1);
    E.init_done_j2 = false;

    %% Time grids, in hours since Date(1).  These mirror the timetable
    %% axes: tt_h is hourly with an extra trailing hour, tt_d is every
    %% 24th hourly stamp, and TT_m is a regular 10-minute grid.
    res = MOD_PARAM_input.t_res;
    E.n_t = round((NN/24)/(1/(24*res)));
    E.th  = (0:NN)';                       %% hourly stamps  [h]
    E.td  = E.th(1:24:end);                %% daily stamps   [h]
    E.tm  = ((0:E.n_t-1)')/res;            %% 10-minute grid [h]
    E.res = res;
    %% Keep the hourly CO2 series. The daily window below uses the same
    %% timetable interpolation as the single-column frame because direct
    %% interp1 rounding can select a different nonlinear chemistry branch.
    E.CO2_h = [Ca(1); Ca(:)]/22.41;
end
C = E.TC_const;

%% ---------------------------------------------------------------------
%% Water fluxes for this column  (full hourly series, as in the 1-D frame)
%% ---------------------------------------------------------------------
cfTL = (sum((ones(length(Ccrown),1)*(p.Bio_Zs>0)).*p.RfL_Zs,2))';
cfTH = (sum((ones(length(Ccrown),1)*(p.Bio_Zs>0)).*p.RfH_Zs,2))';
E.TC.Tw = sum(cfTH.*T_H_c + cfTL.*T_L_c, 2)/1000;   %% [m/h]
E.TC.Iw = f_c/1000;                                  %% [m/h]
E.TC.Lw = Lk_c/1000;                                 %% [m/h]
%% The single-column frame carries the same two statements in the opposite
%% order (MAIN_FRAME_BG_ERW_10mins ~700):
%%     ERW_TC.Lw(1) = Lk(2)/1000;   <- overwritten on the next line
%%     ERW_TC.Lw    = Lk/1000;
%% so there the override is dead and Lw(1) keeps Lk(1). Here it used to
%% run, which shifted the first daily drainage bin: the first ERW call saw
%% Lk(2)+Lk(3)+Lk(4) against the frame's Lk(2). ERW state is cumulative, so
%% that day-0 difference never washed out; it biased Al by 3%, H by 1.4%
%% and Ca and Mg by 0.6% across every subsequent call.
%% Matching the frame, because the frame produced the reference results and
%% the 2-D driver is required to reproduce it at N_col = 1.
%% Both are papering over Lk(1), which is structurally 0 because the
%% hydrology loop starts at i = 2. That is worth fixing on its own terms,
%% in both drivers together, rather than silently here.


pdind = max(1,i-24):i-1;

%% Leaching fraction of the active-layer store
Lk_tot = sum(Lk_c(pdind));
Lk_tot(Lk_tot > VSUM) = VSUM;
if VSUM == 0
    E.TC.f_Lw_bio(j) = 0;
else
    E.TC.f_Lw_bio(j) = Lk_tot/VSUM;
end
%% Passive uptake fraction from transpiration in the active layers
T_tot = sum(cfTH.*sum(T_H_c(pdind,:),1) + cfTL.*sum(T_L_c(pdind,:),1), 2);
T_tot(T_tot > VTSUM) = VTSUM;
if VTSUM == 0
    E.TC.f_T_bio(j) = 0;
else
    E.TC.f_T_bio(j) = T_tot/VTSUM;
end

r_het_now = (R_microbe_c(j) + R_litter_c(j) + R_ew_c(j) - R_litter_sur_c(j))/12*conv_mol;
r_aut_now = (Rmr_H_c(j,:) + Rmr_L_c(j,:))*Ccrown'/12*conv_mol;

if ~E.init_done_j2
    E.TC.s(1:2,:)         = erw_pore_saturation(p,O_c,pdind,C.n,cfg.saturation_exponent);
    E.TC.temp_soil(1:2,:) = Tdp_bio;
    E.TC.r_het(1:2,:)     = (R_microbe_c(2) + R_litter_c(2) + R_ew_c(2) - R_litter_sur_c(2))/12*conv_mol;
    E.TC.r_aut(1:2,:)     = (Rmr_H_c(2,:) + Rmr_L_c(2,:))*Ccrown'/12*conv_mol;
else
    E.TC.s(j,:)         = erw_pore_saturation(p,O_c,pdind,C.n,cfg.saturation_exponent);
    E.TC.temp_soil(j,:) = Tdp_bio;
    E.TC.r_het(j,:)     = r_het_now;
    E.TC.r_aut(j,:)     = r_aut_now;
end

E.TC.CO2_atm = Ca(:)/22.41;      %% [ppm] -> [umol-conv/L]

cov_moll = C.n * C.Zbio/1000 * E.TC.s(j,:) * 1000;
E.TC_const.cov_moll = cov_moll;

%% Solution fractions for the BG -> ERW handoff, identical to the 1-D
%% frame (MAIN_FRAME_BG_ERW_10mins).  P(:,43), P(:,32) and P(:,31) are
%% labile MINERAL pools; only the dissolved part carries pore-water
%% charge, and these are the model's own solubility coefficients, the
%% same ones Biogeo_Leakage and Biogeo_uptake apply to those pools.
BiogeoPar_erw_handoff = Biogeochemistry_Parameter(1);
erw_PO4_fraction = BiogeoPar_erw_handoff.aP;
erw_NO3_fraction = BiogeoPar_erw_handoff.aNO3;
erw_NH4_fraction = BiogeoPar_erw_handoff.aNH4;

%% ---------------------------------------------------------------------
%% Seed dissolved and exchangeable pools on the first ERW day
%% ---------------------------------------------------------------------
if ~E.init_done_j2
    bg.Ca  = P_c(1,60)/40  /cov_moll*1e6;
    bg.Mg  = P_c(1,68)/24.3/cov_moll*1e6;
    bg.K   = P_c(1,52)/39.1/cov_moll*1e6;
    bg.Si  = P_c(2,76)/28.1/cov_moll*1e6;
    bg.PO4 = erw_PO4_fraction * P_c(2,43)/31  /cov_moll*1e6;
    bg.NO3 = erw_NO3_fraction * P_c(2,32)/14  /cov_moll*1e6;
    bg.NH4 = erw_NH4_fraction * P_c(2,31)/14  /cov_moll*1e6;
    bg.excCa = P_c(1,61)/40  /cov_moll*1e6;
    bg.excMg = P_c(1,69)/24.3/cov_moll*1e6;
    bg.excK  = P_c(1,53)/39.1/cov_moll*1e6;
    E.bg_seed = bg;
end

%% ---------------------------------------------------------------------
%% Rock application schedule for THIS column
%% ---------------------------------------------------------------------
t   = E.t;
dtw = E.dtw;
[kinetic_mass_fraction, direct_carbonate_mass_fraction, rock_f_in_erw] = ...
    ERW_Feedstock_Pathway(ERW_PARAM_input.mineral, ERW_PARAM_input.rock_f_in);
if ERW_PARAM_input.diss_f > 0 && cfg.apply_rock
    %% Every 10 November from the first field application (2016) to the
    %% column's last application year: 2019 is the field record, a
    %% projection extends it. Dates outside the run are dropped by
    %% ERW_Application_Indices.
    application_years = 2016:cfg.last_application_year;
    time_application = days(datetime(application_years,11,10) - ...
        datetime([Datam(1,:) 0 0]));
    M_rock_add = zeros(length(t),1);
    M_rock_add_CaCO3 = zeros(length(t),1);
    M_rock_add(1) = ERW_PARAM_input.M_rock_in*kinetic_mass_fraction;
    M_rock_add_CaCO3(1) = ERW_PARAM_input.M_rock_in*direct_carbonate_mass_fraction/100.087*1e6;
    idx = ERW_Application_Indices(time_application,dtw,length(t));
    M_rock_add(idx) = M_rock_add(idx) + ...
        ERW_PARAM_input.M_rock_event*kinetic_mass_fraction;
    M_rock_add_CaCO3(idx) = M_rock_add_CaCO3(idx) + ...
        ERW_PARAM_input.M_rock_event*direct_carbonate_mass_fraction/100.087*1e6;
else
    M_rock_add = zeros(length(t),1);
    M_rock_add_CaCO3 = zeros(length(t),1);
    if cfg.apply_rock
        M_rock_add(1) = ERW_PARAM_input.M_rock_in*kinetic_mass_fraction;
        M_rock_add_CaCO3(1) = ERW_PARAM_input.M_rock_in*direct_carbonate_mass_fraction/100.087*1e6;
    end
end
ERW_PARAM_input.num_rock_app = sum(M_rock_add > 0 | M_rock_add_CaCO3 > 0);
E.TC_const.M_rock_add = M_rock_add;
E.TC_const.M_rock_add_CaCO3 = M_rock_add_CaCO3;

%% ---------------------------------------------------------------------
%% Window on the ERW time axis, evaluated without rebuilding timetables
%% ---------------------------------------------------------------------
%% Index arithmetic identical to the timetable lookup it replaces:
%%   dm_start is the grid row at hour (i-25), dm_end the last row
%%   strictly before hour (i-1), so the window is exactly 24*t_res rows.
if i == NN
    h_first = find(Datam(:,4)==1,1,'last') - 1;
    h_last  = NN;                       %% length(tt_h) - 1 in hours
else
    h_first = max(1,i-24) - 1;
    h_last  = i - 1;
end
dm_start = find(E.tm == h_first, 1, 'first');
dm_end   = find(E.tm <  h_last,  1, 'last');
if isempty(dm_start); dm_start = 1; end
if isempty(dm_end);   dm_end = min(E.n_t, dm_start); end
dm_end = min(dm_end, E.n_t);

%% Rows this day reads, plus one earlier row for the s_pre lookup.
rows = (max(dm_start-1,1):dm_end)';
tw_r = E.tm(rows);

%% Daily drivers: use the exact timetable interpolation operator from the
%% single-column frame. Restricting it to the bracketing daily interval is
%% bitwise identical to the full synchronization without rebuilding the
%% full sub-daily series for every column and day.
h0 = 24*floor(tw_r(1)/24);
h1 = 24*ceil(tw_r(end)/24);
hidx = find(E.th >= h0 & E.th <= h1);
didx = find(E.td >= h0 & E.td <= h1);
TT_h_window = timetable(hours(E.th(hidx)), E.CO2_h(hidx), ...
    'VariableNames', {'CO2_atm'});
TT_d_window = timetable(hours(E.td(didx)), E.TC.s(didx,:), ...
    E.TC.temp_soil(didx,:), E.TC.r_het(didx,:), E.TC.r_aut(didx,:), ...
    'VariableNames', {'s','temp_soil','r_het','r_aut'});
TT_linear_window = synchronize(TT_h_window, TT_d_window, 'regular', ...
    'linear', 'TimeStep', minutes(60/E.res));
window_first_row = round(h0*E.res) + 1;
window_rows = rows - window_first_row + 1;
W.s         = TT_linear_window.s(window_rows,:);
W.temp_soil = TT_linear_window.temp_soil(window_rows,:);
W.r_het     = TT_linear_window.r_het(window_rows,:);
W.r_aut     = TT_linear_window.r_aut(window_rows,:);
W.CO2_atm   = TT_linear_window.CO2_atm(window_rows,:);

%% Daily water sums held zero-order from the calendar day of each row,
%% as synchronize(...,'previous') over retime(...,'daily','sum').
Twv = [E.TC.Tw(1); E.TC.Tw];
Lwv = [E.TC.Lw(1); E.TC.Lw];
Iwv = [E.TC.Iw(1); E.TC.Iw];
day_of_row = floor(tw_r/24);
W.Tw = zeros(numel(rows),1); W.Lw = W.Tw; W.Iw = W.Tw;
for d_u = unique(day_of_row)'
    k0 = 24*d_u + 1;
    k1 = min(24*d_u + 24, numel(Twv));
    sel = (day_of_row == d_u);
    W.Tw(sel) = sum(Twv(k0:k1));
    W.Lw(sel) = sum(Lwv(k0:k1));
    W.Iw(sel) = sum(Iwv(k0:k1));
end
%% Map a grid row index onto its position within the window arrays.
row_of = @(dm) dm - rows(1) + 1;

%% ---------------------------------------------------------------------
%% First ERW day: build the parameter and state structures
%% ---------------------------------------------------------------------
if ~E.init_done_j2
    %% ERW_MOD_PARAM_excel preallocates on length(t) and reads the first
    %% element of each forcing series plus max(r_het+r_aut); supply the
    %% full-length columns it expects.
    TT_h_init = timetable(hours(E.th), E.CO2_h, ...
        'VariableNames', {'CO2_atm'});
    TT_d_init = timetable(hours(E.td), E.TC.s, E.TC.temp_soil, ...
        E.TC.r_het, E.TC.r_aut, ...
        'VariableNames', {'s','temp_soil','r_het','r_aut'});
    TT_init = synchronize(TT_h_init, TT_d_init, 'regular', 'linear', ...
        'TimeStep', minutes(60/E.res));
    TT_init = TT_init(1:E.n_t,:);
    %% Daily water sums, built exactly as the window block below builds
    %% W.Tw and W.Lw. These used to be passed as zeros on the assumption
    %% that the seed only reads the first element of each series plus
    %% max(r_het+r_aut). It also reads Lw, and that feeds the initial
    %% chemistry, so the 2-D driver seeded itself with no drainage at all:
    %% sum(Lw) was 0 against the frame's 0.00647, which set the ERW state
    %% off from the first solve and never recovered.
    Twv_a = [E.TC.Tw(1); E.TC.Tw];
    Lwv_a = [E.TC.Lw(1); E.TC.Lw];
    day_a = floor(E.tm/24);
    TT_init.Tw = zeros(E.n_t,1);
    TT_init.Lw = zeros(E.n_t,1);
    for d_a = unique(day_a)'
        ka0 = 24*d_a + 1;
        ka1 = min(24*d_a + 24, numel(Lwv_a));
        sel_a = (day_a == d_a);
        TT_init.Tw(sel_a) = sum(Twv_a(ka0:ka1));
        TT_init.Lw(sel_a) = sum(Lwv_a(ka0:ka1));
    end
    [E.bg, E.plant, E.weath, E.const, E.var, ...
     E.bg_names, E.weath_names, ~, ~, ~] = ...
        ERW_MOD_PARAM_excel(t, TT_init, E.TC_const, E.bg_seed, ERW_PARAM_input, 1);
    E.P_initial = ERW_Initialize_TC_Cation_Pools( ...
        P_c(1,:), E.bg, E.const, TT_init.s(1));
    P_c(1,:) = E.P_initial;
    E.init_done_j2 = true;
    E.init = true;
end

%% Change in the biogeochemical totals since yesterday.  ERW Ca_tot,
%% Mg_tot, and K_tot include the CEC inventory, so internal transfer between
%% the legacy solution and exchange pools must cancel here.
apply_tc_daily_delta = ERW_Should_Apply_TC_Daily_Delta(j,i,NN,Datam);
if apply_tc_daily_delta
    dCa  = ((P_c(j,60) + P_c(j,61)) - (P_prev_for_erw_delta(60) + P_prev_for_erw_delta(61)))/40  *1e6;
    dMg  = ((P_c(j,68) + P_c(j,69)) - (P_prev_for_erw_delta(68) + P_prev_for_erw_delta(69)))/24.3*1e6;
    dK   = ((P_c(j,52) + P_c(j,53)) - (P_prev_for_erw_delta(52) + P_prev_for_erw_delta(53)))/39.1*1e6;
    dSi  = (P_c(j,76) - P_prev_for_erw_delta(76))/28.1*1e6;
    dPO4 = erw_PO4_fraction * (P_c(j,43) - P_prev_for_erw_delta(43))/31  *1e6;
    dNO3 = erw_NO3_fraction * (P_c(j,32) - P_prev_for_erw_delta(32))/14  *1e6;
    dNH4 = erw_NH4_fraction * (P_c(j,31) - P_prev_for_erw_delta(31))/14  *1e6;
else
    dCa=0; dMg=0; dK=0; dSi=0; dPO4=0; dNO3=0; dNH4=0;
end

%% Background input
if ERW_PARAM_input.keyword_add == 1
    if j == 2
        bg_row = 1;
    else
        bg_row = E.dm_end;
    end
    E.const.I_Na = mean(E.TC.f_Lw_bio(j)+E.TC.f_T_bio(j)) * E.bg.Na_tot(bg_row);
    E.const.I_An = mean(E.TC.f_Lw_bio(j)+E.TC.f_T_bio(j)) * E.bg.An_tot(bg_row);
    %% The leached part of that closure returns as the salt it left as, so
    %% the bare anion does not acidify the box (ERW_Leached_Salt_Resupply).
    [E.const.I_Ca,E.const.I_Mg] = ERW_Leached_Salt_Resupply(E.TC.f_Lw_bio(j), ...
        E.bg.An_tot(bg_row),E.bg.Na_tot(bg_row),E.bg.Ca(bg_row),E.bg.Mg(bg_row));
end

%% ---------------------------------------------------------------------
%% Sub-daily ERW integration over this day's window
%% ---------------------------------------------------------------------
bg = E.bg; weath = E.weath; cst = E.const; var_ = E.var;
%% The particle-class histories (rows x classes x applications) of an
%% amended column are each about 3.6 GB over 35 years. Written in place
%% here they would be copied whole on every daily call, since the frame
%% still holds the column state, so their rows collect in daily blocks and
%% the frame assembles them after the run (ERW_Particle_History).
%% ERW_function reads exactly the rows it read from the full arrays.
if ~isfield(E,'hist')
    [E.hist,weath] = ERW_Particle_History('start',weath,E.weath_names);
end
hist = E.hist;
day = ERW_Particle_History('open',hist,max(dm_start,2),max(0,dm_end-max(dm_start,2)+1));
bulk_pH_params_2d = cst.bulk_pH_params;
bulk_pH_params_2d.suppress_no_converge_warning = true;
erw_window_days = (dm_end-max(dm_start,2)+1)*dtw;
d_erw = ERW_Normalize_TC_Daily_Deltas( ...
    [dCa; dMg; dK; dSi; dPO4; dNO3; dNH4],erw_window_days);
dCa=d_erw(1,:); dMg=d_erw(2,:); dK=d_erw(3,:); dSi=d_erw(4,:);
dPO4=d_erw(5,:); dNO3=d_erw(6,:); dNH4=d_erw(7,:);

%% NH4, NO3 and dissolved PO4 are biogeochemical boundary conditions. The
%% PO4 partition follows the same daily equilibrium fraction used by T&C
%% uptake and leakage. The final flush reuses P(j) without its delta again.
if j >= 2
    anchor_row = max(dm_start,2) - 1;
    anchor_state_row = j - double(apply_tc_daily_delta);
    bg.PO4_tot(anchor_row,:) = ...
        erw_PO4_fraction * P_c(anchor_state_row,43)/31 * 1e6;
    bg.NH4_tot(anchor_row,:) = ...
        erw_NH4_fraction * P_c(anchor_state_row,31)/14 * 1e6;
    bg.NO3_tot(anchor_row,:) = ...
        erw_NO3_fraction * P_c(anchor_state_row,32)/14 * 1e6;
    %% Ca/Mg/K/Si_tot must be re-anchored the same way, but ONLY when
    %% lateral transport can actually move mass into or out of this column
    %% (N_col > 1). In the 2-D driver, lateral transport (Biogeo_Lateral_
    %% Pools, called every hour) keeps writing into
    %% P(j,60/61/68/69/52/53/76,:) for all 24 hours of a day, but this
    %% function only runs once, at hour 1. Relying on P_c(j-1,:) as the
    %% delta baseline (below) captures only hour 1's lateral change:
    %% P_c(j-1,:) is exactly the value BIOGEO_UNIT copied forward into
    %% P_c(j,:) before hour 1's lateral step ran, so it cancels out of
    %% (P_c(j,:) - P_c(j-1,:)) and hours 2-24 of the PREVIOUS day's lateral
    %% transport -- already folded into that same P_c(j-1,:) snapshot --
    %% disappear from the signal entirely. Anchoring the absolute total to
    %% P_c(anchor_state_row,:) every day, exactly as PO4/NH4/NO3 do above,
    %% discards no history: each day's ERW state starts from T&C's true
    %% accumulated total instead of ERW's own dCa-integrated total.
    %%
    %% Unlike the PO4/NH4/NO3 anchors, this one is not a bit-exact
    %% round-trip: it reads back TWO independently-rounded stored columns
    %% (e.g. P_c(60)+P_c(61)) and rescales them, where PO4/NH4/NO3 read a
    %% single column with no addition. At N_col = 1 nothing external ever
    %% touches these pools between ERW calls -- the lateral block above is
    %% itself gated on N_col > 1 -- so re-deriving Ca_tot from P_c there is
    %% a mathematical no-op that still costs a few ULPs relative to ERW's
    %% own already-exact running total, which then amplifies through the
    %% same nonlinear-solver sensitivity already documented above for pH
    %% (confirmed against the 1-D frame: unconditional, this broke
    %% test_2d_driver_equivalence's bit-exact fields by day 6). Gating on
    %% N_col > 1 keeps the fix exactly where the 23/24-hour loss actually
    %% happens and leaves the N_col = 1 case, which has no loss to fix,
    %% bit-identical to before.
    if N_col > 1
        bg.Ca_tot(anchor_row,:) = ...
            (P_c(anchor_state_row,60) + P_c(anchor_state_row,61))/40 * 1e6;
        bg.Mg_tot(anchor_row,:) = ...
            (P_c(anchor_state_row,68) + P_c(anchor_state_row,69))/24.3 * 1e6;
        bg.K_tot(anchor_row,:) = ...
            (P_c(anchor_state_row,52) + P_c(anchor_state_row,53))/39.1 * 1e6;
        bg.Si_tot(anchor_row,:) = ...
            P_c(anchor_state_row,76)/28.1 * 1e6;
    end
end
for dt_m = max(dm_start,2):dm_end
    cst.s_pre = W.s(row_of(dt_m-1));
    rock_app  = sum(M_rock_add(1:dt_m) > 0);
    prev = struct();
    for nm = {'psd','d','lamb','SSA'}
        if any(strcmp(hist.names,nm{1}))
            prev.(nm{1}) = ERW_Particle_History('row',hist,day,weath,nm{1},dt_m-1);
        else
            prev.(nm{1}) = weath.(nm{1})(dt_m-1,:,:);
        end
    end
    if any(strcmp(hist.names,'psd'))
        psd_0 = ERW_Particle_History('first',hist,day,'psd');
    else
        psd_0 = weath.psd;
    end
    if any(strcmp(hist.names,'lamb'))
        lamb_0 = ERW_Particle_History('first',hist,day,'lamb');
    else
        lamb_0 = weath.lamb;
    end

    [~, ~, bg_it, weath_it, ~, ~, cst_it] = ERW_function( ...
        cst.n, cst.s_pre, W.s(row_of(dt_m)), W.Lw(row_of(dt_m-1)), W.Tw(row_of(dt_m-1)), W.Iw(row_of(dt_m)), ...
        cst.Zrw, W.r_het(row_of(dt_m-1)), W.r_aut(row_of(dt_m-1)), W.temp_soil(row_of(dt_m)), ...
        cst.K_CEC, bg.CEC_tot(dt_m-1,:), ...
        cst.M_rock_in, cst.k_diss_H, cst.k_diss_w, cst.k_diss_OH, ...
        cst.E_H, cst.E_w, cst.E_OH, cst.R, cst.T_ref(dt_m-1,:), ...
        cst.K1, cst.K2, cst.K3, cst.K4, ...
        cst.rho_rock, cst.b, cst.K_CaCO3, cst.K_MgCO3, cst.MM_min, ...
        cst.min_st, dtw, cst.Al_CEC_scale, conv_mol, ...
        bg.Ca(dt_m-1,:), bg.Mg(dt_m-1,:), bg.Na(dt_m-1,:), bg.K(dt_m-1,:), bg.An(dt_m-1,:), ...
        var_.fsolve_options, cst.Dw_0, cst.D_0, cst.pH_rain, W.CO2_atm(row_of(dt_m)), ...
        bg.CO2_air(dt_m-1,:), cst.Z_CO2, ...
        bg.Ca_tot(dt_m-1,:), bg.Mg_tot(dt_m-1,:), bg.K_tot(dt_m-1,:), ...
        bg.Na_tot(dt_m-1,:), bg.An_tot(dt_m-1,:), bg.Al_tot(dt_m-1,:), bg.Si_tot(dt_m-1,:), ...
        bg.Al(dt_m-1,:), bg.AlOH4(dt_m-1,:), bg.Si(dt_m-1,:), ...
        bg.R_alk(dt_m-1,:), bg.H(dt_m-1,:), ...
        weath.EW(dt_m-1,:), bg.IC_tot(dt_m-1,:), weath.Wr(dt_m-1,:), ...
        prev.psd, prev.d, prev.lamb, prev.SSA, ...
        bg.Fs(dt_m-1,:), bg.DIC(dt_m-1,:), ...
        cst.I_Ca, cst.I_Mg, cst.I_Na, cst.I_K, cst.I_An, cst.I_Si, ...
        weath.W_CaCO3(dt_m-1), weath.W_MgCO3(dt_m-1), 1, ...
        bg.CaCO3(dt_m-1,:), bg.MgCO3(dt_m-1,:), ...
        bg.f_CEC(dt_m-1,:), bg.f_Ca(dt_m-1,:), bg.f_Mg(dt_m-1,:), bg.f_K(dt_m-1,:), ...
        bg.f_Na(dt_m-1,:), bg.f_Al(dt_m-1,:), bg.f_H(dt_m-1,:), ...
        cst.r_CaCO3, cst.r_MgCO3, cst.tau_CaCO3, cst.tau_MgCO3, ...
        bg.M_min(dt_m-1,:), cst.M_iner, cst.a, psd_0, cst.d_in, lamb_0, ...
        dt_m, var_.mineral, cst.K_sp, cst.diss_f, cst.n_H, cst.n_OH, ...
        E.opt_TC_biochem, ...
        dCa, dMg, dK, dSi, ...
        cov_moll, M_rock_add(dt_m,:), cst.psd_perc_in, rock_f_in_erw, ...
        cst.SSA_in, rock_app, ...
        E.TC.f_Lw_bio(j), E.TC.f_T_bio(j), ...
        bg.PO4_tot(dt_m-1,:), bg.PO4(dt_m-1,:), dPO4, ...
        bg.NO3_tot(dt_m-1,:), bg.NO3(dt_m-1,:), dNO3, ...
        bg.NH4_tot(dt_m-1,:), bg.NH4(dt_m-1,:), dNH4, ...
        ERW_PARAM_input.keyword_add, ...
        M_rock_add_CaCO3(dt_m,:), ...
        cst.K_AlOH3, bulk_pH_params_2d);

    for v = 1:size(E.bg_names,1)
        nm = cell2mat(E.bg_names(v));
        if isfield(bg_it, nm)
            bg.(nm)(dt_m,:) = bg_it.(nm);
        end
    end
    for v = 1:size(E.weath_names,1)
        nm = cell2mat(E.weath_names(v));
        if isfield(weath_it, nm) && ~any(strcmp(hist.names,nm))
            if size(weath.(nm),3) > 1
                weath.(nm)(dt_m,:,:) = weath_it.(nm);
            else
                weath.(nm)(dt_m,:) = weath_it.(nm);
            end
        end
    end
    if ~isempty(hist.names)
        day = ERW_Particle_History('write',hist,day,weath_it,dt_m);
        hist = ERW_Particle_History('track',hist,day,weath_it,dt_m,weath);
    end
    cst.M_iner = cst_it.M_iner;
end
hist = ERW_Particle_History('close',hist,day);
E.hist = hist;

%% Forward extraction of the solved pore-water state to the air-dried
%% bulk extract pH (Porewater_to_Bulk_pH), as in the 1-D frame.
pH_bulk_now = nan(size(bg.pH(dm_end,:)));
for col = 1:numel(bg.pH(dm_end,:))
    bulk_pH_params_now = cst.bulk_pH_params;
    s_col = W.s(row_of(dm_end));
    bulk_pH_params_now.field_water_L_m2 = cst.n * s_col * cst.Zrw * 1000;
    bulk_pH_params_now.DIC_field_umol_L = bg.DIC(dm_end,col);
    %% The paste is made from the plot's soil, so it carries the plot's
    %% solid calcite: applied grains not yet dissolved plus the pedogenic
    %% pool. M_min is g m-2 and MM_min g per umol, so the sum is umol m-2.
    solid_calcite_now = max(bg.CaCO3(dm_end,col),0);
    calcite_names = strtrim(string(var_.mineral));
    is_calcite = strcmpi(calcite_names(:).',"calcite");
    if any(is_calcite)
        M_min_now = reshape(bg.M_min(dm_end,:),1,[]);
        MM_min_now = reshape(cst.MM_min,1,[]);
        solid_calcite_now = solid_calcite_now + ...
            sum(M_min_now(is_calcite)./MM_min_now(is_calcite));
    end
    bulk_pH_params_now.solid_CaCO3_umol_m2 = solid_calcite_now;
    %% It carries the plot's solid Mg carbonate, the nesquehonite pool
    %% (umol m-2), on the same terms.
    bulk_pH_params_now.solid_MgCO3_umol_m2 = max(bg.MgCO3(dm_end,col),0);
    [pH_bulk_now(col),bulk_pH_diag_now] = Porewater_to_Bulk_pH( ...
        bg.H(dm_end, col), ...
        bg.Ca(dm_end, col), bg.Mg(dm_end, col), ...
        bg.K(dm_end, col), bg.Na(dm_end, col), ...
        bg.Al_w(dm_end, col), bg.Al(dm_end, col), ...
        bg.An(dm_end, col), ...
        bg.NO3(dm_end, col), bg.PO4(dm_end, col), ...
        bg.NH4(dm_end, col), ...
        W.CO2_atm(row_of(dm_end)), ...
        bg.k1(dm_end, col), bg.k2(dm_end, col), ...
        bg.k_w(dm_end, col), bg.k_nh3(dm_end, col), bg.k_H(dm_end, col), ...
        cst.K1, cst.K2, cst.K3, cst.K4, ...
        bg.f_Ca(dm_end, col), bg.f_Mg(dm_end, col), ...
        bg.f_K(dm_end, col), bg.f_Na(dm_end, col), ...
        bg.f_Al(dm_end, col), bg.f_H(dm_end, col), ...
        bg.CEC_tot(dm_end, col), cst.rho_bulk, ...
        cst.Zrw, conv_mol, bulk_pH_params_now, ...
        cst.Al_CEC_scale);
    cst.bulk_pH_params.H_prev = bulk_pH_diag_now.H;
    cst.bulk_pH_params.H_prev_valid = true;
end
E.PH = bg.pH(dm_end,:);
E.PH_pore = bg.pH(dm_end,:);
E.PH_bulk = pH_bulk_now;
bg.pH_bulk(dm_end,:) = pH_bulk_now;
E.bg = bg; E.weath = weath; E.const = cst;
E.dm_end = dm_end;
E.current_s = W.s(row_of(dm_end));
cov_moll_end = cst.n * cst.Zrw * E.current_s * 1000;

%% ---------------------------------------------------------------------
%% Hand the updated pore-water chemistry back to the biogeochemistry
%% ---------------------------------------------------------------------
E.P60 = bg.Ca(dm_end,:) * 40   * cov_moll_end / 1e6;
E.P61 = bg.f_Ca(dm_end,:) .* bg.CEC_tot(dm_end,:) / 2 * 40 / 1e6;
E.P68 = bg.Mg(dm_end,:) * 24.3 * cov_moll_end / 1e6;
E.P69 = bg.f_Mg(dm_end,:) .* bg.CEC_tot(dm_end,:) / 2 * 24.3 / 1e6;
E.P52 = bg.K(dm_end,:)  * 39.1 * cov_moll_end / 1e6;
E.P53 = bg.f_K(dm_end,:) .* bg.CEC_tot(dm_end,:) * 39.1 / 1e6;
E.P76 = bg.Si(dm_end,:) * 28.1 * cov_moll_end / 1e6;
[E.P43,P_weathering_release_umol] = ERW_Return_Weathered_Phosphorus( ...
    P_c(j,43),bg.PO4_tot(dm_end,:),erw_PO4_fraction, ...
    weath.P_weathering_release(max(dm_start,2):dm_end,:));
E.P_weathering_release_umol = P_weathering_release_umol;
E.has_update = true;

end

function s = erw_pore_saturation(p,O_c,pdind,n,exponent)
%% Active-layer water content over its porosity: the saturation the ERW
%% module weathers, breathes and dilutes with. The uncertainty ensemble
%% raises it to an exponent (sqrt(s) to s^2, the legacy T&C-SMEW range,
%% which perturbed the module's whole s input the same way). Exponent 1
%% leaves the value untouched, bit for bit.
s = sum(p.Bio_Zs.*mean(O_c(pdind,:)))/n;
if exponent ~= 1
    s = s.^exponent;
end
end
