function [T,fig] = Analyse_Sensitivity(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Summarise the nitrogen sensitivity batch produced by Sensitivity_Run.
%%
%% Reports, per parameter set and for the equilibrium years only, the
%% control values and the change the feedstock produces. Year one is a
%% spin-up transient: NO3 starts near 0.737 gN/m2 and roughly doubles
%% before settling, so including it understates every flux.
%%
%% The rN2O partition sweep is done here rather than by running the model.
%% The nitrate balance carries -DENITflx whatever the N2O/N2 split, so
%% pH0_rN2O and k_rN2O cannot change NO3, and their effect on N2O is
%% recovered exactly by re-splitting the saved DENITflx.
%%
%% Options
%%   Directory   sensitivity root         [SIM/Sensitivity]
%%   Runs        treatment group          [10:12]
%%   NDays       shortened-run tag        []
%%   Visible     figure visibility        ['on']
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sens_dir = fileparts(mfilename('fullpath'));
v2 = fileparts(sens_dir);
addpath(v2,sens_dir,fullfile(v2,'T&C_Code'), ...
    fullfile(v2,'T&C_Code','Parameters'));

opt = struct('Directory',fullfile(v2,'SIM','Sensitivity'), ...
    'Runs',10:12,'NDays',[],'Visible','on');
for k = 1:2:numel(varargin); opt.(char(varargin{k})) = varargin{k+1}; end

[~,sets] = Sensitivity_Parameters();
RM = Run_Matrix_12();
BP = Biogeochemistry_Parameter(1);
labels = {'control','basalt','calcite'};

suffix = '';
if ~isempty(opt.NDays); suffix = sprintf('__%dd',opt.NDays); end

row = 0; T = table();
keep = false(numel(sets),1);
store = cell(numel(sets),1);
for s = 1:numel(sets)
    sdir = fullfile(opt.Directory,sets{s}.tag);
    ok = true; D = struct();
    for i = 1:3
        id = opt.Runs(i);
        f = fullfile(sdir,sprintf('TC_2D_run%02d_%s%s.mat',id,RM(id).name,suffix));
        if ~isfile(f); ok = false; break; end
        S = load(f,'P','N2flx','N2Oflx','LEAK_NO3','PH_d', ...
            'NH4uptake_H','NH4uptake_L','NO3uptake_H','NO3uptake_L');
        P = S.P; if ndims(P) == 3; P = P(:,:,1); end
        nd = size(P,1); sl = min(nd,round(365.25)+1):nd;
        g = 10/(numel(sl)/365.25);
        n2o = col(S,'N2Oflx',nd); n2 = col(S,'N2flx',nd);
        D(i).n2o   = sum(n2o(sl),'omitnan')*g;
        D(i).den   = sum(n2o(sl)+n2(sl),'omitnan')*g;
        D(i).leach = sum(subcol(S,'LEAK_NO3',nd,sl),'omitnan')*g;
        D(i).up    = uptake(S,nd,sl)*g;
        ph = col(S,'PH_d',nd); D(i).pH = mean(ph(sl),'omitnan');
        D(i).ph_series = ph(sl); D(i).n2_series = n2(sl); D(i).n2o_series = n2o(sl);
        D(i).no3pool = mean(P(sl,32));
    end
    if ~ok
        fprintf('  skipping %-20s (results missing)\n',sets{s}.tag); continue
    end
    keep(s) = true; store{s} = D;
    row = row+1;
    T(row,:) = table(string(sets{s}.tag), ...
        D(1).leach, D(1).n2o, D(1).den, D(1).up, D(1).pH, ...
        D(2).leach-D(1).leach, D(2).n2o-D(1).n2o, ...
        100*(D(1).leach-D(2).leach)/max(D(1).leach,eps), ...
        100*(D(1).n2o-D(2).n2o)/max(D(1).n2o,eps), ...
        'VariableNames',{'set','ctl_leach','ctl_N2O','ctl_denit','ctl_uptake', ...
        'ctl_pH','d_leach','d_N2O','leach_cut_pct','N2O_cut_pct'});
end
if isempty(T); error('Analyse_Sensitivity:NoResults','No completed sets found.'); end

fprintf('\n=== control values and feedstock effect, equilibrium years ===\n');
disp(T);

%% Exact rN2O sweep on the baseline set, which needs no extra simulation.
idx = find(keep);
kept_tags = cellfun(@(x) x.tag,sets(idx),'UniformOutput',false);
bi = find(strcmp(kept_tags,'baseline'),1);
if true
    if ~isempty(bi)
        Db = store{idx(bi)};
        fprintf('\n=== rN2O partition sweep, baseline set, no extra runs ===\n');
        fprintf('%6s %5s %12s %12s %11s\n','pH0','k','control N2O','basalt N2O','cut %%');
        for pH0 = [5.0 5.5 6.06 6.2]
            for kk = [1.5 2.5]
                v = zeros(1,2);
                for i = 1:2
                    r0 = 1./(1+exp(BP.k_rN2O*(Db(i).ph_series-BP.pH0_rN2O)));
                    DEN = Db(i).n2_series./max(1-r0,eps);
                    nit = Db(i).n2o_series - r0.*DEN;
                    rr = 1./(1+exp(kk*(Db(i).ph_series-pH0)));
                    n = numel(Db(i).ph_series);
                    v(i) = sum(nit+rr.*DEN,'omitnan')*10/(n/365.25);
                end
                fprintf('%6.2f %5.1f %12.1f %12.1f %10.1f%%\n', ...
                    pH0,kk,v(1),v(2),100*(v(1)-v(2))/max(v(1),eps));
            end
        end
    end
end

%% Spread across sets: the point of the exercise is the range, so plot the
%% control value and the feedstock effect side by side rather than ranking.
fig = figure('Color','w','Visible',char(opt.Visible), ...
    'Name','Nitrogen sensitivity');
fig.Position(3:4) = [1150 620];
tl = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
bars = {'ctl_leach','NO_3 leaching, control','kg N ha^{-1} yr^{-1}'; ...
        'ctl_N2O','N_2O, control','kg N ha^{-1} yr^{-1}'; ...
        'leach_cut_pct','NO_3 leaching cut by feedstock','%'; ...
        'N2O_cut_pct','N_2O cut by feedstock','%'};
for p = 1:4
    ax = nexttile(tl);
    y = T.(bars{p,1});
    barh(ax,y,'FaceColor',[0 0.25 0.60],'EdgeColor','none');
    set(ax,'YTick',1:height(T),'YTickLabel',T.set,'YDir','reverse', ...
        'Color','w','XColor','k','YColor','k','TickLabelInterpreter','none');
    grid(ax,'off'); box(ax,'on');
    xlabel(ax,bars{p,3},'Color','k');
    text(ax,0.98,0.04,bars{p,2},'Units','normalized', ...
        'HorizontalAlignment','right','FontWeight','bold','Color','k');
end
title(tl,'Nitrogen parameter sensitivity, equilibrium years','Color','k');
out = fullfile(opt.Directory,'sensitivity_summary.png');
exportgraphics(fig,out,'Resolution',200);
writetable(T,fullfile(opt.Directory,'sensitivity_summary.csv'));
fprintf('\nSaved: %s\n',out);
fprintf('Saved: %s\n',fullfile(opt.Directory,'sensitivity_summary.csv'));

end

%% =======================================================================
function v = col(S,name,nd)
if isfield(S,name); v = S.(name); v = v(:); v = v(1:nd); else; v = zeros(nd,1); end
end
function v = subcol(S,name,nd,sl)
v = col(S,name,nd); v = v(sl);
end
function u = uptake(S,nd,sl)
u = 0;
for f = {'NH4uptake_H','NH4uptake_L','NO3uptake_H','NO3uptake_L'}
    if isfield(S,f{1})
        x = S.(f{1}); if size(x,1) >= nd; x = sum(x(1:nd,:),2); u = u + sum(x(sl)); end
    end
end
end
