function E = ERW_Daily_Step_2D(E, cfg, p, i, j, NN, NNd, Datam, ~, Ca, ...
    O_c, Lk_c, f_c, T_H_c, T_L_c, Tdp_bio, VSUM, VTSUM, ...
    R_microbe_c, R_litter_c, R_ew_c, R_litter_sur_c, Rmr_H_c, Rmr_L_c, Ccrown, ...
    P_c, MOD_PARAM_input, ERW_PARAM_input, conv_mol, conv_Al)
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
    %% CO2 never changes during the run, so interpolate it once.
    CO2v  = [Ca(1); Ca(:)]/22.41;
    E.CO2_m = interp1(E.th, CO2v, E.tm, 'linear');
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
E.TC.Lw(1) = Lk_c(2)/1000;

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
    E.TC.s(1:2,:)         = sum(p.Bio_Zs.*mean(O_c(pdind,:)))/C.n;
    E.TC.temp_soil(1:2,:) = Tdp_bio;
    E.TC.r_het(1:2,:)     = (R_microbe_c(2) + R_litter_c(2) + R_ew_c(2) - R_litter_sur_c(2))/12*conv_mol;
    E.TC.r_aut(1:2,:)     = (Rmr_H_c(2,:) + Rmr_L_c(2,:))*Ccrown'/12*conv_mol;
else
    E.TC.s(j,:)         = sum(p.Bio_Zs.*mean(O_c(pdind,:)))/C.n;
    E.TC.temp_soil(j,:) = Tdp_bio;
    E.TC.r_het(j,:)     = r_het_now;
    E.TC.r_aut(j,:)     = r_aut_now;
end

E.TC.CO2_atm = Ca(:)/22.41;      %% [ppm] -> [umol-conv/L]

cov_moll = C.n * C.Zbio/1000 * E.TC.s(j,:) * 1000;
E.TC_const.cov_moll = cov_moll;

%% ---------------------------------------------------------------------
%% Seed dissolved and exchangeable pools on the first ERW day
%% ---------------------------------------------------------------------
if ~E.init_done_j2
    bg.Ca  = P_c(1,60)/40  /cov_moll*1e6;
    bg.Mg  = P_c(1,68)/24.3/cov_moll*1e6;
    bg.K   = P_c(1,52)/39.1/cov_moll*1e6;
    bg.Si  = P_c(2,76)/28.1/cov_moll*1e6;
    bg.PO4 = P_c(2,43)/31  /cov_moll*1e6;
    bg.NO3 = P_c(2,32)/14  /cov_moll*1e6;
    bg.NH4 = P_c(2,31)/14  /cov_moll*1e6;
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
if ERW_PARAM_input.diss_f > 0 && cfg.apply_rock
    ta1 = days(datetime([2016 11 10 0 0 0])-datetime([Datam(1,:) 0 0]));
    ta2 = days(datetime([2017 11 10 0 0 0])-datetime([Datam(1,:) 0 0]));
    ta3 = days(datetime([2018 11 10 0 0 0])-datetime([Datam(1,:) 0 0]));
    ta4 = days(datetime([2019 11 10 0 0 0])-datetime([Datam(1,:) 0 0]));
    time_application = [ta1 ta2 ta3 ta4];
    M_rock_add = zeros(length(t),1);
    M_rock_add_CaCO3 = zeros(length(t),1);
    M_rock_add(1) = ERW_PARAM_input.M_rock_in;
    idx = time_application/dtw;
    idx = idx(idx >= 1 & idx <= length(t));
    if contains(ERW_PARAM_input.ExperimentalDescription,"Beerling")
        M_rock_add(idx) = 5000*(1-0.026);
        M_rock_add_CaCO3(idx) = 5000*0.026/100.087*1e6;   %% umol
    else
        M_rock_add(idx) = 5000;
    end
else
    M_rock_add = zeros(length(t),1);
    M_rock_add_CaCO3 = zeros(length(t),1);
    if cfg.apply_rock
        M_rock_add(1) = ERW_PARAM_input.M_rock_in;
    end
end
ERW_PARAM_input.num_rock_app = sum(M_rock_add > 0);
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

%% Daily drivers: linear interpolation, as synchronize(...,'linear').
W.s         = interp1(E.td, E.TC.s,         tw_r, 'linear');
W.temp_soil = interp1(E.td, E.TC.temp_soil, tw_r, 'linear');
W.r_het     = interp1(E.td, E.TC.r_het,     tw_r, 'linear');
W.r_aut     = interp1(E.td, E.TC.r_aut,     tw_r, 'linear');
W.CO2_atm   = E.CO2_m(rows);

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
    TT_init.temp_soil = interp1(E.td, E.TC.temp_soil, E.tm, 'linear');
    TT_init.CO2_atm   = E.CO2_m;
    TT_init.s         = interp1(E.td, E.TC.s,     E.tm, 'linear');
    TT_init.r_het     = interp1(E.td, E.TC.r_het, E.tm, 'linear');
    TT_init.r_aut     = interp1(E.td, E.TC.r_aut, E.tm, 'linear');
    TT_init.Tw        = zeros(E.n_t,1);
    TT_init.Lw        = zeros(E.n_t,1);
    [E.bg, E.plant, E.weath, E.const, E.var, ...
     E.bg_names, E.weath_names, ~, ~, ~] = ...
        ERW_MOD_PARAM_excel(t, TT_init, E.TC_const, E.bg_seed, ERW_PARAM_input, 1);
    E.init_done_j2 = true;
    E.init = true;
end

%% Change in the biogeochemical totals since yesterday
if j > 2
    dCa  = (P_c(j,60) - P_c(j-1,60))/40  *1e6;
    dMg  = (P_c(j,68) - P_c(j-1,68))/24.3*1e6;
    dK   = (P_c(j,52) - P_c(j-1,52))/39.1*1e6;
    dSi  = (P_c(j,76) - P_c(j-1,76))/28.1*1e6;
    dPO4 = (P_c(j,43) - P_c(j-1,43))/31  *1e6;
    dNO3 = (P_c(j,32) - P_c(j-1,32))/14  *1e6;
    dNH4 = (P_c(j,31) - P_c(j-1,31))/14  *1e6;
else
    dCa=0; dMg=0; dK=0; dSi=0; dPO4=0; dNO3=0; dNH4=0;
end

%% Background input
if ERW_PARAM_input.keyword_add == 1
    if j == 2
        E.const.I_Na = mean(E.TC.f_Lw_bio(j)+E.TC.f_T_bio(j)) * E.bg.Na_tot(1);
        E.const.I_An = mean(E.TC.f_Lw_bio(j)+E.TC.f_T_bio(j)) * E.bg.An_tot(1);
    else
        E.const.I_Na = mean(E.TC.f_Lw_bio(j)+E.TC.f_T_bio(j)) * E.bg.Na_tot(E.dm_end);
        E.const.I_An = mean(E.TC.f_Lw_bio(j)+E.TC.f_T_bio(j)) * E.bg.An_tot(E.dm_end);
    end
end

%% ---------------------------------------------------------------------
%% Sub-daily ERW integration over this day's window
%% ---------------------------------------------------------------------
bg = E.bg; weath = E.weath; cst = E.const; var_ = E.var;
for dt_m = max(dm_start,2):dm_end
    cst.s_pre = W.s(row_of(dt_m-1));
    rock_app  = sum(M_rock_add(1:dt_m) > 0);

    [~, ~, bg_it, weath_it, ~, ~, cst_it] = ERW_function( ...
        cst.n, cst.s_pre, W.s(row_of(dt_m)), W.Lw(row_of(dt_m-1)), W.Tw(row_of(dt_m-1)), W.Iw(row_of(dt_m)), ...
        cst.Zrw, W.r_het(row_of(dt_m-1)), W.r_aut(row_of(dt_m-1)), W.temp_soil(row_of(dt_m)), ...
        cst.K_CEC, bg.CEC_tot(dt_m-1,:), ...
        cst.M_rock_in, cst.k_diss_H, cst.k_diss_w, cst.k_diss_OH, ...
        cst.E_H, cst.E_w, cst.E_OH, cst.R, cst.T_ref(dt_m-1,:), ...
        cst.K1, cst.K2, cst.K3, cst.K4, ...
        cst.rho_rock, cst.b, cst.K_CaCO3, cst.K_MgCO3, cst.MM_min, ...
        cst.min_st, dtw, conv_Al, conv_mol, ...
        bg.Ca(dt_m-1,:), bg.Mg(dt_m-1,:), bg.Na(dt_m-1,:), bg.K(dt_m-1,:), bg.An(dt_m-1,:), ...
        var_.fsolve_options, cst.Dw_0, cst.D_0, cst.Alk_rain, W.CO2_atm(row_of(dt_m)), ...
        bg.CO2_air(dt_m-1,:), cst.Z_CO2, ...
        bg.Ca_tot(dt_m-1,:), bg.Mg_tot(dt_m-1,:), bg.K_tot(dt_m-1,:), ...
        bg.Na_tot(dt_m-1,:), bg.An_tot(dt_m-1,:), bg.Al_tot(dt_m-1,:), bg.Si_tot(dt_m-1,:), ...
        bg.Al(dt_m-1,:), bg.AlOH4(dt_m-1,:), bg.Si(dt_m-1,:), ...
        bg.R_alk(dt_m-1,:), bg.H(dt_m-1,:), ...
        weath.EW(dt_m-1,:), bg.IC_tot(dt_m-1,:), weath.Wr(dt_m-1,:), ...
        weath.psd(dt_m-1,:,:), weath.d(dt_m-1,:,:), weath.lamb(dt_m-1,:,:), weath.SSA(dt_m-1,:,:), ...
        bg.Fs(dt_m-1,:), bg.DIC(dt_m-1,:), ...
        cst.I_Ca, cst.I_Mg, cst.I_Na, cst.I_K, cst.I_An, cst.I_Si, ...
        weath.W_CaCO3(dt_m-1), weath.W_MgCO3(dt_m-1), 1, ...
        bg.CaCO3(dt_m-1,:), bg.MgCO3(dt_m-1,:), ...
        bg.f_CEC(dt_m-1,:), bg.f_Ca(dt_m-1,:), bg.f_Mg(dt_m-1,:), bg.f_K(dt_m-1,:), ...
        bg.f_Na(dt_m-1,:), bg.f_Al(dt_m-1,:), bg.f_H(dt_m-1,:), ...
        cst.r_CaCO3, cst.r_MgCO3, cst.tau_CaCO3, cst.tau_MgCO3, ...
        bg.M_min(dt_m-1,:), cst.M_iner, cst.a, weath.psd, cst.d_in, weath.lamb, ...
        dt_m, var_.mineral, cst.K_sp, cst.diss_f, cst.n_H, cst.n_OH, ...
        E.opt_TC_biochem, ...
        dCa, dMg, dK, dSi, ...
        cov_moll, M_rock_add(dt_m,:), cst.psd_perc_in, cst.rock_f_in, ...
        cst.SSA_in, rock_app, ...
        E.TC.f_Lw_bio(j), E.TC.f_T_bio(j), ...
        bg.PO4_tot(dt_m-1,:), bg.PO4(dt_m-1,:), dPO4, ...
        bg.NO3_tot(dt_m-1,:), bg.NO3(dt_m-1,:), dNO3, ...
        bg.NH4_tot(dt_m-1,:), bg.NH4(dt_m-1,:), dNH4, ...
        ERW_PARAM_input.keyword_add, ...
        M_rock_add_CaCO3(dt_m,:));

    for v = 1:size(E.bg_names,1)
        nm = cell2mat(E.bg_names(v));
        bg.(nm)(dt_m,:) = bg_it.(nm);
    end
    for v = 1:size(E.weath_names,1)
        nm = cell2mat(E.weath_names(v));
        if size(weath.(nm),3) > 1
            weath.(nm)(dt_m,:,:) = weath_it.(nm);
        else
            weath.(nm)(dt_m,:) = weath_it.(nm);
        end
    end
    cst.M_iner = cst_it.M_iner;
end
E.bg = bg; E.weath = weath; E.const = cst;
E.dm_end = dm_end;
E.current_s = W.s(row_of(dm_end));

%% ---------------------------------------------------------------------
%% Hand the updated chemistry back to the biogeochemistry
%% ---------------------------------------------------------------------
E.PH  = bg.pH(dm_end,:);
E.P60 = bg.Ca(dm_end,:) * 40   * cov_moll / 1e6;
E.P68 = bg.Mg(dm_end,:) * 24.3 * cov_moll / 1e6;
E.P52 = bg.K(dm_end,:)  * 39.1 * cov_moll / 1e6;
E.P76 = bg.Si(dm_end,:) * 28.1 * cov_moll / 1e6;
E.P43 = bg.PO4(dm_end,:)* 31   * cov_moll / 1e6;
E.has_update = true;

end
