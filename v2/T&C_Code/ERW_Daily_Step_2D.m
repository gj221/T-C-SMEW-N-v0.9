function E = ERW_Daily_Step_2D(E, col, cfg, p, i, j, NN, Datam, Date, Ca, ...
    O_pd, Lk_pd, f_now, T_H_pd, T_L_pd, Tdp_bio, ...
    R_microbe_c, R_litter_c, R_ew_c, R_litter_sur_c, Rmr_H_c, Rmr_L_c, Ccrown, ...
    P_c, MOD_PARAM_input, ERW_PARAM_input, conv_mol, conv_Al, Parameter_Directory)
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
%%   margin in the same simulation.
%%
%% DIFFERENCE FROM THE SINGLE-COLUMN ERW BLOCK
%%   The 1-D frame rebuilds a timetable spanning the whole simulation on
%%   every simulated day and then reads one day out of it.  That is
%%   ~263,000 rows rebuilt 1825 times for a 5-year run; multiplied by one
%%   column per cell it dominates everything else.  Here the sub-hourly
%%   forcing window is built for the current day only, by
%%   ERW_Window_Interp, using the same interpolation rules.
%%
%% Returns E with .has_update set when PH and the pool updates are valid.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

E.has_update = false;
res_per_h = MOD_PARAM_input.t_res;
dtw   = 1/(24*res_per_h);          %% [day] ERW timestep
t_end = NN/24;                     %% [day] simulated days
n_t   = round(t_end/dtw);          %% ERW steps over the whole run

%% ---------------------------------------------------------------------
%% One-off initialisation for this column
%% ---------------------------------------------------------------------
if ~E.init
    t = dtw:dtw:t_end;
    E.t = t;
    E.dtw = dtw;
    E.opt_TC_biochem = ERW_PARAM_input.opt_TC_biochem;

    C.n        = sum(p.Bio_Zs.*p.Osat);      %% porosity [-]
    C.rsd_bio  = sum(p.Bio_Zs.*p.rsd);
    C.PH_in    = p.PH;
    C.Zbio     = p.Zbio;
    C.ms       = p.ms;
    C.conv_mol = conv_mol;
    C.conv_Al  = conv_Al;

    %% Rock application schedule for THIS column.  Mass is per unit area,
    %% so a column receives the same areal rate as the field trial; the
    %% area weighting happens when results are aggregated.
    C.M_rock_add       = zeros(n_t,1);
    C.M_rock_add_CaCO3 = zeros(n_t,1);
    if cfg.apply_rock
        C.M_rock_add(1) = ERW_PARAM_input.M_rock_in;
        if strcmp(char(MOD_PARAM_input.Location),'Kantola') && ERW_PARAM_input.diss_f > 0
            ta = zeros(1,4);
            for y = 1:4
                ta(y) = days(datetime([2015+y 11 10 0 0 0]) - datetime([Datam(1,:) 0 0]));
            end
            idx = round(ta/dtw);
            idx = idx(idx >= 1 & idx <= n_t);
            if contains(ERW_PARAM_input.ExperimentalDescription,"Beerling")
                C.M_rock_add(idx)       = 5000*(1-0.026);
                C.M_rock_add_CaCO3(idx) = 5000*0.026/100.087*1e6;   %% umol
            else
                C.M_rock_add(idx) = 5000;
            end
        end
    end
    ERW_PARAM_input.num_rock_app = sum(C.M_rock_add > 0);
    E.TC_const = C;

    %% Seed dissolved and exchangeable pools from the biogeochemistry.
    s0 = sum(p.Bio_Zs.*mean(O_pd,1))/C.n;
    cov_moll = C.n * p.Zbio/1000 * s0 * 1000;      %% [mol/l] -> [mol]
    E.TC_const.cov_moll = cov_moll;

    B = P_c(1,:);
    bg.Ca  = B(60)/40  /cov_moll*1e6;
    bg.Mg  = B(68)/24.3/cov_moll*1e6;
    bg.K   = B(52)/39.1/cov_moll*1e6;
    bg.Si  = P_c(2,76)/28.1/cov_moll*1e6;
    bg.PO4 = P_c(2,43)/31  /cov_moll*1e6;
    bg.NO3 = P_c(2,32)/14  /cov_moll*1e6;
    bg.NH4 = P_c(2,31)/14  /cov_moll*1e6;
    bg.excCa = B(61)/40  /cov_moll*1e6;
    bg.excMg = B(69)/24.3/cov_moll*1e6;
    bg.excK  = B(53)/39.1/cov_moll*1e6;

    %% ERW_MOD_PARAM_excel preallocates on length(t) and reads the FIRST
    %% element of each forcing series, plus max(r_het+r_aut) for the
    %% initial soil-CO2 estimate.  Supply full-length series seeded from
    %% the opening state, which reproduces those reads exactly.
    r_het0 = (R_microbe_c(2) + R_litter_c(2) + R_ew_c(2) - R_litter_sur_c(2))/12*conv_mol;
    r_aut0 = (Rmr_H_c(2,:) + Rmr_L_c(2,:))*Ccrown'/12*conv_mol;
    one = ones(n_t,1);
    TT0.temp_soil = one * Tdp_bio;
    TT0.CO2_atm   = one * (Ca(1)/22.41);
    TT0.Tw        = zeros(n_t,1);
    TT0.Lw        = zeros(n_t,1);
    TT0.s         = one * s0;
    TT0.r_het     = one * r_het0;
    TT0.r_aut     = one * r_aut0;

    [E.bg, E.plant, E.weath, E.const, E.var, ...
     E.bg_names, E.weath_names, ~, ~, ~, ~] = ...
        ERW_MOD_PARAM_excel(t, TT0, E.TC_const, bg, ERW_PARAM_input, 1);

    E.prev = struct('s', s0, 'temp', TT0.temp_soil(1), ...
                    'r_het', r_het0, 'r_aut', r_aut0);
    E.dm_end = 1;
    E.init = true;
end

%% ---------------------------------------------------------------------
%% Today's forcing for this column
%% ---------------------------------------------------------------------
C = E.TC_const;

%% Transpiration drawn from within the biogeochemically active depth
cfTL = sum((ones(numel(Ccrown),1)*(p.Bio_Zs>0)).*p.RfL_Zs, 2)';
cfTH = sum((ones(numel(Ccrown),1)*(p.Bio_Zs>0)).*p.RfH_Zs, 2)';
Tw_day = sum(cfTH.*sum(T_H_pd,1) + cfTL.*sum(T_L_pd,1)) / 1000;   %% [m/d]
Lw_day = sum(Lk_pd)/1000;                                          %% [m/d]
Iw_day = f_now/1000*24;                                            %% [m/d]

s_now    = sum(p.Bio_Zs.*mean(O_pd,1))/C.n;
temp_now = Tdp_bio;
r_het_now = (R_microbe_c(j) + R_litter_c(j) + R_ew_c(j) - R_litter_sur_c(j))/12*conv_mol;
r_aut_now = (Rmr_H_c(j,:) + Rmr_L_c(j,:))*Ccrown'/12*conv_mol;

%% Leaching and passive-uptake fractions of the active-layer water store
VSUM  = sum(p.Bio_Zs.*mean(O_pd,1)) * p.Zbio;
VTSUM = VSUM;
Lk_tot = min(sum(Lk_pd), VSUM);
T_tot  = min(sum(cfTH.*sum(T_H_pd,1) + cfTL.*sum(T_L_pd,1)), VTSUM);
if VSUM  == 0; f_Lw = 0; else; f_Lw = Lk_tot/VSUM;  end
if VTSUM == 0; f_T  = 0; else; f_T  = T_tot/VTSUM;  end

%% Hourly CO2 spanning today plus the first hour of tomorrow
h0 = max(1, i-24);
h1 = min(numel(Ca), h0+24);
CO2_h_day = Ca(h0:h1)/22.41;
if numel(CO2_h_day) < 25
    CO2_h_day(end+1:25) = CO2_h_day(end);
end

W = ERW_Window_Interp(dtw, j, res_per_h, CO2_h_day, ...
    [E.prev.s; s_now], [E.prev.temp; temp_now], ...
    [E.prev.r_het; r_het_now], [E.prev.r_aut; r_aut_now], ...
    Tw_day, Lw_day, Iw_day);

cov_moll = C.n * p.Zbio/1000 * s_now * 1000;
E.TC_const.cov_moll = cov_moll;

%% Change in the biogeochemical totals since yesterday, handed to ERW
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

if ERW_PARAM_input.keyword_add == 1
    E.const.I_Na = (f_Lw + f_T) * E.bg.Na_tot(max(E.dm_end,1));
    E.const.I_An = (f_Lw + f_T) * E.bg.An_tot(max(E.dm_end,1));
end

%% ---------------------------------------------------------------------
%% Sub-daily ERW integration
%% ---------------------------------------------------------------------
dm_start = max(2, E.dm_end + 1);
dm_end   = min(dm_start + W.n_sub - 1, numel(E.t));
if dm_end < dm_start
    return
end

bg = E.bg; weath = E.weath; cst = E.const; var_ = E.var;
for dt_m = dm_start:dm_end
    wi = dt_m - dm_start + 1;
    cst.s_pre = W.s(max(wi-1,1));
    rock_app = sum(C.M_rock_add(1:dt_m) > 0);

    [~, ~, bg_it, weath_it, ~, ~, cst_it] = ERW_function( ...
        cst.n, cst.s_pre, W.s(wi), W.Lw(max(wi-1,1)), W.Tw(max(wi-1,1)), W.Iw(wi), ...
        cst.Zrw, W.r_het(max(wi-1,1)), W.r_aut(max(wi-1,1)), W.temp_soil(wi), ...
        cst.K_CEC, bg.CEC_tot(dt_m-1,:), ...
        cst.M_rock_in, cst.k_diss_H, cst.k_diss_w, cst.k_diss_OH, ...
        cst.E_H, cst.E_w, cst.E_OH, cst.R, cst.T_ref(dt_m-1,:), ...
        cst.K1, cst.K2, cst.K3, cst.K4, ...
        cst.rho_rock, cst.b, cst.K_CaCO3, cst.K_MgCO3, cst.MM_min, ...
        cst.min_st, dtw, conv_Al, conv_mol, ...
        bg.Ca(dt_m-1,:), bg.Mg(dt_m-1,:), bg.Na(dt_m-1,:), bg.K(dt_m-1,:), bg.An(dt_m-1,:), ...
        var_.fsolve_options, cst.Dw_0, cst.D_0, cst.Alk_rain, W.CO2_atm(wi), ...
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
        cov_moll, C.M_rock_add(dt_m,:), cst.psd_perc_in, cst.rock_f_in, ...
        cst.SSA_in, rock_app, f_Lw, f_T, ...
        bg.PO4_tot(dt_m-1,:), bg.PO4(dt_m-1,:), dPO4, ...
        bg.NO3_tot(dt_m-1,:), bg.NO3(dt_m-1,:), dNO3, ...
        bg.NH4_tot(dt_m-1,:), bg.NH4(dt_m-1,:), dNH4, ...
        ERW_PARAM_input.keyword_add, ...
        C.M_rock_add_CaCO3(dt_m,:));

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

E.prev = struct('s', s_now, 'temp', temp_now, 'r_het', r_het_now, 'r_aut', r_aut_now);

end

