function [fig,plot_file,stats] = Plot_1D_Group(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Weathering signal and nitrogen partitioning for one treatment group.
%%
%% Complements Plot_Subfield_12 rather than repeating it. That figure
%% already carries cumulative N2O, cumulative NO3 leaching, the soil NO3
%% pool, pore-water pH, NEE and crop yield, so none of those appear here.
%% What is missing there and shown here:
%%
%%   base cations   Ca, Mg and Si solution pools, which carry the
%%                  enhanced weathering signal
%%   denitrification total N2 + N2O, where Plot_Subfield_12 shows only the
%%                  N2O component
%%   WFPS           water-filled pore space, the control on how much of the
%%                  mineral N is denitrified at all
%%
%% The denitrification N2O fraction rN2O is still computed, and reported
%% flux-weighted in the printed statistics, but is not plotted: as a daily
%% curve it is just the partition function re-evaluated on that day's pH,
%% which reads as a larger treatment effect than the run actually produces.
%%
%% Colours, line weights, axis treatment, application markers and legend
%% placement follow Plot_Subfield_12 so the two figures read as one set.
%%
%% Options
%%   ResultDirectory  result folder              [SIM/Subfield_12]
%%   Runs             one group of three runs    [10:12]
%%   NDays            shortened-run filename tag []
%%   Visible          figure visibility          ['on']
%%   OutputFile       output PNG and companion FIG [in the result folder]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

here = fileparts(mfilename('fullpath'));
opt = struct('ResultDirectory',fullfile(here,'SIM','Subfield_12'), ...
    'Runs',10:12,'NDays',[],'Visible','on','OutputFile','');
if mod(numel(varargin),2) ~= 0
    error('Plot_1D_Group:BadOptions','Options must be name/value pairs.');
end
for k = 1:2:numel(varargin)
    name = char(varargin{k});
    if ~isfield(opt,name)
        error('Plot_1D_Group:BadOptions','Unknown option %s.',name);
    end
    opt.(name) = varargin{k+1};
end
run_ids = opt.Runs(:).';
if numel(run_ids) ~= 3
    error('Plot_1D_Group:BadRuns', ...
        'Expected one complete treatment group of three runs, got %d.', ...
        numel(run_ids));
end
addpath(here,fullfile(here,'T&C_Code'), ...
    fullfile(here,'T&C_Code','Parameters'));
RM = Run_Matrix_12();

%% Shared with Plot_Subfield_12 so the figures match.
treatment_labels = {'Control','Basalt','Calcite'};
treatment_colors = [0 0 0; 0 0.25 0.60; 0.35 0.70 0.90];
application_color = [0.45 0.45 0.45];

%% Pool column indices in P.
iCa = 60; iMg = 68; iSi = 76;

suffix = '';
if ~isempty(opt.NDays); suffix = sprintf('__%dd',opt.NDays); end

Date_d = datetime.empty(0,1);
for i = 1:3
    id = run_ids(i);
    result_file = fullfile(opt.ResultDirectory, ...
        sprintf('TC_2D_run%02d_%s%s.mat',id,RM(id).name,suffix));
    if ~isfile(result_file)
        error('Plot_1D_Group:MissingResult','Missing result: %s',result_file);
    end
    S = load(result_file);
    P = S.P; if ndims(P) == 3; P = P(:,:,1); end
    nd = size(P,1);
    if i == 1
        Date_d = datetime(S.Datam(1,1),S.Datam(1,2),S.Datam(1,3)) + ...
            caldays((0:nd-1).');
        ca = zeros(nd,3); mg = zeros(nd,3); si = zeros(nd,3);
        den = zeros(nd,3); n2o = zeros(nd,3); wfps = zeros(nd,3);
        ph = zeros(nd,3); n2 = zeros(nd,3);
    end
    ca(:,i) = P(:,iCa); mg(:,i) = P(:,iMg); si(:,i) = P(:,iSi);
    n2o(:,i) = daily(S,'N2Oflx',nd);
    n2(:,i)  = daily(S,'N2flx',nd);
    den(:,i) = n2(:,i) + n2o(:,i);
    wfps(:,i) = water_filled_pore_space(S,nd);
    ph(:,i) = daily(S,'PH_d',nd);
end

%% Denitrification N2O fraction, evaluated from the model's own partition.
%%
%% Do not compute this as N2Oflx/(N2Oflx+N2flx) from the saved fluxes. That
%% ratio is not the denitrification partition, because N2Oflx carries a
%% second source:
%%
%%   N2Oflx = f_N2O_nit*NO3flx  +  rN2O*DENITflx
%%
%% The first term is nitrification N2O and has nothing to do with the N2O
%% against N2 split. On days when denitrification is near zero but
%% nitrification continues, that ratio runs to 1.0, which reads as complete
%% suppression of N2O reductase when it is only a mixed numerator. The saved
%% output does not separate the two terms, so reconstruct rN2O directly from
%% daily pH using the same expression BIOGEOCHEMISTRY_DYNAMIC3 uses. It is
%% then defined on every day, needs no window and no threshold.
BiogeoPar = Biogeochemistry_Parameter(1);
rN2O = 100./(1+exp(BiogeoPar.k_rN2O*(ph-BiogeoPar.pH0_rN2O)));
rN2O(ph <= 0) = NaN;

%% Denitrification-weighted rN2O over the equilibrium years, used to label
%% panel 5. The daily curve is the partition function evaluated on that
%% day's pH; what the nitrogen actually experiences is this curve weighted
%% by the denitrification flux, and the two differ because denitrification
%% concentrates on wet days that do not carry the mean pH. Without the
%% label the curve reads as a larger treatment effect than the run produces.
nd_all = size(den,1);
eq = min(nd_all,round(365.25)+1):nd_all;
rN2O_fw = zeros(1,3);
for i = 1:3
    rN2O_fw(i) = rN2O_fluxweighted(rN2O(eq,i),n2(eq,i));
end

time_value = {ca,mg,si,cumsum(den,1),wfps};
time_label = {'Ca solution','Mg solution','Si solution', ...
    'Cumulative denitrification (N_2 + N_2O)','Water-filled pore space'};
time_units = {'\mu mol','\mu mol','\mu mol','g N m^{-2}','%'};
time_tag = {'CaSolution','MgSolution','SiSolution', ...
    'CumulativeDenitrification','WFPS'};

application_dates = datetime([2016 11 10; 2017 11 10; ...
    2018 11 10; 2019 11 10]);
application_dates = application_dates( ...
    application_dates >= Date_d(1) & application_dates <= Date_d(end));

group_label = strrep(RM(run_ids(1)).layout,'_','\_');
fig = figure('Color','w','Visible',char(opt.Visible), ...
    'Name','T&C-SMEW-N one-group weathering and nitrogen partitioning');
fig.Position(3:4) = [1150 720];
%% Two rows: the three base cations that carry the weathering signal, then
%% the two nitrogen panels. A six-column grid with spans lets a row of three
%% sit above a row of two without leaving an empty tile.
tl = tiledlayout(fig,2,6,'TileSpacing','compact','Padding','compact');
tile_span = {[1 2],[1 2],[1 2],[1 3],[1 3]};

time_axes = gobjects(5,1);
legend_lines = gobjects(3,1);
application_legend = gobjects(1);
for measurement = 1:5
    ax = nexttile(tl,tile_span{measurement});
    time_axes(measurement) = ax;
    ax.Tag = sprintf('TimeMetric%d',measurement);
    ln = plot(ax,Date_d,time_value{measurement},'LineWidth',1.15);
    for treatment = 1:3
        ln(treatment).Color = treatment_colors(treatment,:);
    end
    set(ln,'Tag',time_tag{measurement});
    if ~isempty(application_dates)
        marker = xline(ax,application_dates,':','Color',application_color, ...
            'LineWidth',1.1,'Tag','FeedstockApplication');
        if measurement == 1; application_legend = marker(1); end
    end
    format_time_axis(ax);
    ylabel(ax,time_units{measurement},'Color','k');
    if measurement >= 4; xlabel(ax,'Date','Color','k'); end
    text(ax,0.02,0.96,time_label{measurement},'Units','normalized', ...
        'VerticalAlignment','top','Color','k','FontWeight','bold');
    if measurement == 1; legend_lines = ln; end
end
linkaxes(time_axes,'x');
set(time_axes(1:3),'XTickLabel',[]);

legend_handles = legend_lines;
legend_labels = treatment_labels;
if isgraphics(application_legend)
    legend_handles(end+1) = application_legend;
    legend_labels{end+1} = 'Feedstock application';
end
lgd = legend(legend_handles,legend_labels,'Orientation','horizontal', ...
    'Box','off','TextColor','k','Color','w','EdgeColor','none');
lgd.Layout.Tile = 'south';
title(tl,sprintf('%s: weathering signal and nitrogen partitioning', ...
    group_label),'Color','k');

if isempty(opt.OutputFile)
    plot_file = fullfile(opt.ResultDirectory, ...
        sprintf('%s_weathering_nitrogen%s.png',RM(run_ids(1)).layout,suffix));
else
    plot_file = char(opt.OutputFile);
end
exportgraphics(fig,plot_file,'Resolution',200);
[plot_folder,plot_name] = fileparts(plot_file);
savefig(fig,fullfile(plot_folder,[plot_name '.fig']));
fprintf('Saved: %s\n',plot_file);
fprintf('Saved: %s\n',fullfile(plot_folder,[plot_name '.fig']));

%% Equilibrium statistics, excluding year one. The spin-up starts NO3 near
%% 0.737 gN/m2 and roughly doubles before settling, so a first-year mean
%% reports the transient rather than the model's own steady state.
nd = size(den,1);
sl = min(nd,round(365.25)+1):nd;
f = 10/(numel(sl)/365.25);          %% gN/m2 over the window -> kg N/ha/yr
stats = struct();
fprintf('\nyears 2 onward, kg N/ha/yr\n');
fprintf('%-8s %7s %8s %10s %8s %9s\n', ...
    'run','N2O','denit','rN2O','WFPS','Ca end');
for i = 1:3
    st = struct('treatment',treatment_labels{i}, ...
        'N2O',sum(n2o(sl,i),'omitnan')*f, ...
        'denit',sum(den(sl,i),'omitnan')*f, ...
        'rN2O',rN2O_fw(i), ...
        'wfps',mean(wfps(sl,i),'omitnan'), ...
        'Ca_end',ca(end,i),'Mg_end',mg(end,i),'Si_end',si(end,i));
    stats.(lower(treatment_labels{i})) = st;
    fprintf('%-8s %7.1f %8.1f %9.1f%% %7.1f%% %9.2f\n', ...
        st.treatment,st.N2O,st.denit,st.rN2O,st.wfps,st.Ca_end);
end

end

%% =======================================================================
function r = rN2O_fluxweighted(share_pct,n2_flux)
%% Denitrification-weighted mean of the N2O fraction [%].
%%
%% The plain mean over days is not what the nitrogen actually experiences.
%% Denitrification is concentrated on wet days, and those days do not carry
%% the mean pH, so weighting by the flux changes the answer materially: on
%% the calcite run the unweighted mean reads 7.9% against 9.6% weighted.
%% Recover total denitrification from N2flx = (1-rN2O)*DENITflx.
f = share_pct/100;
den = n2_flux./max(1-f,eps);
r = 100*sum(f.*den,'omitnan')/max(sum(den,'omitnan'),eps);
end

%% =======================================================================
function v = daily(S,name,nd)
%% Daily flux as a column, or zeros when the run did not save it.
if isfield(S,name); v = S.(name); v = v(:); v = v(1:nd); else; v = zeros(nd,1); end
end

%% =======================================================================
function w = water_filled_pore_space(S,nd)
%% Depth-weighted WFPS of the biogeochemical layer [%], daily mean.
p = S.p;
O = S.O; if ndims(O) == 3; O = O(:,:,1); end
sph = size(O,1)/nd;
Osat = sum(p.Bio_Zs.*p.Soil_Param.Osat);
w = zeros(nd,1);
for j = 1:nd
    rows = round((j-1)*sph)+1:round(j*sph);
    w(j) = sum(p.Bio_Zs.*mean(O(rows,:),1))/Osat*100;
end
end

%% =======================================================================
function format_time_axis(ax)
%% Identical to the local formatter in Plot_Subfield_12.
ax.Color = 'w';
ax.XColor = 'k';
ax.YColor = 'k';
ax.XTickLabelRotation = 20;
grid(ax,'off');
end
