function budget(res_dir)
% Alkalinity budget of runs 10 (control) and 11 (basalt): what the rock
% released, where it went (exchanger, pore water, drainage, crop) and the
% exchanger state each August, when the soil cores were taken (Beerling et
% al. 2024, Methods). Observed base saturation is the digitised Fig. 1F
% (Data_Beerling.mat, CEC_Exp). Companion of phobs.m; both read the saved
% ERW state and the T&C leaching and uptake series of the result file.
%
% Units: ERW sub-step fluxes (dAlk_*) are already per step; EW is a rate
% in umol d^-1 and is integrated with the step length; the T&C series are
% daily g m^-2 d^-1.
BS_obs.control = [79.5 80.79 85.15 79.11 72.87; 80.99 77.13 76.53 76.14 82.87];
BS_obs.basalt = [80.69 83.17 89.6 87.82 83.37; 83.27 79.6 81.09 80.59 86.53];
years = 2016:2020;
for run = [10 11]
    d = dir(fullfile(res_dir,sprintf('TC_2D_run%02d_*.mat',run)));
    if numel(d) ~= 1
        error('budget:ResultFile','Expected one result file for run %02d in %s, found %d.', ...
            run,res_dir,numel(d));
    end
    S = load(fullfile(d(1).folder,d(1).name),'ERW','Datam','LEAK_Ca','LEAK_Mg','LEAK_NO3', ...
        'Cauptake_H','Cauptake_L','Mguptake_H','Mguptake_L','PH_bulk_d');
    bg = S.ERW{1}.bg; weath = S.ERW{1}.weath; cst = S.ERW{1}.const;
    t0 = datenum(S.Datam(1,1),S.Datam(1,2),S.Datam(1,3));
    n_t = numel(bg.pH); n_d = numel(S.PH_bulk_d);
    t_sub = t0 + (1:n_t)'/(n_t/n_d);
    v_sub = datevec(t_sub);
    t_day = t0 + (0:n_d-1)'; v_day = datevec(t_day);
    rel = cumsum(bg.dAlk_sil + bg.dAlk_carb)/1e6;
    dAlk_all = cumsum(bg.dAlk_total)/1e6;
    leach = cumsum(2*S.LEAK_Ca(:,1)/40.08 + 2*S.LEAK_Mg(:,1)/24.31);
    upt = cumsum(2*sum(S.Cauptake_H(:,:,1),2)/40.08 + 2*sum(S.Cauptake_L(:,:,1),2)/40.08 + ...
        2*sum(S.Mguptake_H(:,:,1),2)/24.31 + 2*sum(S.Mguptake_L(:,:,1),2)/24.31);
    no3 = cumsum(S.LEAK_NO3(:,1)/14.01);
    BS = 100*(bg.f_Ca + bg.f_Mg + bg.f_K + bg.f_Na);
    fprintf('\n===== run %02d: %s =====\n',run,d(1).name);
    fprintf('CEC %.1f eq m^-2; exchangeable acidity at start %.1f eq m^-2\n', ...
        bg.CEC_tot(1)/1e6,bg.CEC_tot(1)*(bg.f_H(1)+bg.f_Al(1))/1e6);
    fprintf('%-34s','August means'); fprintf('%10d',years); fprintf('\n');
    rows = {'pore-water pH',bg.pH,1; 'bulk (extract) pH',S.PH_bulk_d(:),0; ...
        'base saturation, %%',BS,1; 'f_H, %%',100*bg.f_H,1; 'f_Al, %%',100*bg.f_Al,1; ...
        'pore Alk, meq L^-1',bg.Alk/1e3,1; 'pore Ca, mmol L^-1',bg.Ca/1e3,1; ...
        'pore Mg, mmol L^-1',bg.Mg/1e3,1; 'pore DIC, mmol L^-1',bg.DIC/1e3,1; ...
        'soil pCO2, %%',100*bg.pCO2_atm,1; 'ionic strength, mmol L^-1',1e3*bg.ionic_str,1};
    for r = 1:size(rows,1)
        x = rows{r,2}; if rows{r,3}, v = v_sub; else, v = v_day; end
        fprintf(['  ' rows{r,1} repmat(' ',1,max(0,32-numel(rows{r,1})))]);
        for y = years, w = v(:,1)==y & v(:,2)==8; fprintf('%10.3g',mean(x(w))); end
        fprintf('\n');
    end
    nm = 'control'; if run == 11, nm = 'basalt'; end
    fprintf('  obs base saturation 0-10 cm, %%     '); fprintf('%10.1f',BS_obs.(nm)(1,:)); fprintf('\n');
    fprintf('  obs base saturation 10-30 cm, %%    '); fprintf('%10.1f',BS_obs.(nm)(2,:)); fprintf('\n');
    fprintf('%-34s','Cumulative to 31 Aug, eq m^-2'); fprintf('%10d',years); fprintf('\n');
    cum = {'alkalinity released by rock',rel,1; 'change in total alkalinity',dAlk_all,1; ...
        'Ca+Mg leached (T&C)',leach,0; 'Ca+Mg crop uptake (T&C)',upt,0; ...
        'NO3 leached, mol N m^-2',no3,0};
    for r = 1:size(cum,1)
        x = cum{r,2}; if cum{r,3}, t = t_sub; else, t = t_day; end
        fprintf(['  ' cum{r,1} repmat(' ',1,max(0,32-numel(cum{r,1})))]);
        for y = years, k = find(t <= datenum(y,8,31),1,'last'); fprintf('%10.3g',x(k)); end
        fprintf('\n');
    end
    if run == 11
        EWc = cumsum(weath.EW,1)/1e6*mean(diff(t_sub));
        fprintf('  per mineral, mol m^-2 dissolved by 31 Aug 2020, and %% of mass applied by then:\n');
        MM = cst.MM_min*1e6; rock_f = cst.rock_f_in;
        mass_applied = 4*5000*rock_f;   % four 50 t/ha applications, Nov 2016-2019
        k = find(t_sub <= datenum(2020,8,31),1,'last');
        for m = 1:numel(rock_f)
            fprintf('    mineral %d: %.3g mol m^-2, %.1f %% of %.0f g m^-2\n',m,EWc(k,m), ...
                100*EWc(k,m)*MM(m)/mass_applied(m),mass_applied(m));
        end
        fprintf('  field (Beerling et al. 2024): 32 %% of basalt Ca and 12 %% of Mg lost from 0-10 cm by 2020\n');
    end
end
end
