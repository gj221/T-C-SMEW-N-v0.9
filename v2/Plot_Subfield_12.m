function [summary,fig,plot_file] = Plot_Subfield_12(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot the four scenario groups by treatment class.
%%
%% Options
%%   ResultDirectory  result folder              [SIM/Subfield_12]
%%   NDays            shortened-run filename tag []
%%   Visible          figure visibility          ['on']
%%   OutputFile       output PNG                  [summary in result folder]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

here = fileparts(mfilename('fullpath'));
opt = struct('ResultDirectory',fullfile(here,'SIM','Subfield_12'), ...
    'NDays',[],'Visible','on','OutputFile','');
if mod(numel(varargin),2) ~= 0
    error('Plot_Subfield_12:BadOptions','Options must be name/value pairs.');
end
for k = 1:2:numel(varargin)
    name = char(varargin{k});
    if ~isfield(opt,name)
        error('Plot_Subfield_12:BadOptions','Unknown option %s.',name);
    end
    opt.(name) = varargin{k+1};
end
if ~isempty(opt.NDays) && ~(isnumeric(opt.NDays) && isscalar(opt.NDays) && ...
        isfinite(opt.NDays) && opt.NDays == fix(opt.NDays) && opt.NDays >= 1)
    error('Plot_Subfield_12:BadOptions','NDays must be a positive integer.');
end
if ~ismember(char(opt.Visible),{'on','off'})
    error('Plot_Subfield_12:BadOptions','Visible must be on or off.');
end
if ~isfolder(opt.ResultDirectory)
    error('Plot_Subfield_12:MissingDirectory', ...
        'Result directory not found: %s',opt.ResultDirectory);
end

RM = Run_Matrix_12();
n = numel(RM);
n2o = nan(n,1);
no3_leach = nan(n,1);
no3_final = nan(n,1);
ph_final = nan(n,1);
for k = 1:n
    tag = RM(k).name;
    if ~isempty(opt.NDays); tag = sprintf('%s__%dd',tag,opt.NDays); end
    result_file = fullfile(opt.ResultDirectory, ...
        sprintf('TC_2D_run%02d_%s.mat',k,tag));
    if ~isfile(result_file)
        error('Plot_Subfield_12:MissingResult','Result not found: %s',result_file);
    end
    S = load(result_file,'R','Domain');
    if ~isfield(S,'R') || ~isfield(S,'Domain') || ...
            S.R.id ~= k || ~strcmp(S.R.name,RM(k).name)
        error('Plot_Subfield_12:WrongResult', ...
            'Result does not match run %d: %s',k,result_file);
    end
    required = {'N2Oflx','LEAK_NO3','NO3_pool','PH_d'};
    if ~all(isfield(S.Domain,required)) || ...
            isempty(S.Domain.NO3_pool) || isempty(S.Domain.PH_d)
        error('Plot_Subfield_12:IncompleteResult', ...
            'Domain summary is incomplete in %s.',result_file);
    end
    if ~isnumeric(S.Domain.N2Oflx) || ~isnumeric(S.Domain.LEAK_NO3) || ...
            ~isnumeric(S.Domain.NO3_pool) || ~isnumeric(S.Domain.PH_d) || ...
            ~any(isfinite(S.Domain.N2Oflx(:))) || ...
            ~any(isfinite(S.Domain.LEAK_NO3(:))) || ...
            ~isfinite(S.Domain.NO3_pool(end)) || ~isfinite(S.Domain.PH_d(end))
        error('Plot_Subfield_12:NonfiniteMeasurement', ...
            'Run %d does not contain finite values for every plotted measurement.',k);
    end
    n2o(k) = sum(S.Domain.N2Oflx(:),'omitnan');
    no3_leach(k) = sum(S.Domain.LEAK_NO3(:),'omitnan');
    no3_final(k) = S.Domain.NO3_pool(end);
    ph_final(k) = S.Domain.PH_d(end);
    if any(~isfinite([n2o(k),no3_leach(k),no3_final(k),ph_final(k)]))
        error('Plot_Subfield_12:NonfiniteMeasurement', ...
            'Run %d does not contain finite values for every plotted measurement.',k);
    end
end

summary = table([RM.id].',string({RM.layout}).',string({RM.treatment}).', ...
    n2o,no3_leach,no3_final,ph_final, ...
    'VariableNames',{'RunID','ScenarioGroup','Treatment', ...
    'CumulativeN2O_gN_m2','CumulativeNO3Leaching_gN_m2', ...
    'FinalNO3Pool_gN_m2','FinalSoilPH'});

group_labels = {'Whole field','Margin with crop','Margin with prairie','Vertical 1-D'};
treatment_labels = {'Control','Basalt','Calcite'};
fig = figure('Color','w','Visible',char(opt.Visible), ...
    'Name','T&C-SMEW-N twelve-scenario summary');
fig.Position(3:4) = [1400 800];
tl = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
b = draw_panel(nexttile(tl),n2o,'Cumulative N_2O','g N m^{-2}', ...
    group_labels);
draw_panel(nexttile(tl),no3_leach,'Cumulative NO_3 leaching','g N m^{-2}', ...
    group_labels);
draw_panel(nexttile(tl),no3_final,'Final soil NO_3 pool','g N m^{-2}', ...
    group_labels);
draw_panel(nexttile(tl),ph_final,'Final soil pH','pH', ...
    group_labels);
lgd = legend(b,treatment_labels,'Orientation','horizontal','Box','off', ...
    'TextColor','k','Color','w');
lgd.Layout.Tile = 'north';
title(tl,'Treatment response across four scenario groups','Color','k');

if isempty(opt.OutputFile)
    suffix = '';
    if ~isempty(opt.NDays); suffix = sprintf('__%dd',opt.NDays); end
    plot_file = fullfile(opt.ResultDirectory, ...
        ['Subfield_12_summary' suffix '.png']);
else
    plot_file = char(opt.OutputFile);
end
exportgraphics(fig,plot_file,'Resolution',200);
fprintf('Saved: %s\n',plot_file);

end

function b = draw_panel(ax,value,panel_title,y_label,group_labels)
b = bar(ax,reshape(value,3,4).','grouped');
ax.Color = 'w';
ax.XColor = 'k';
ax.YColor = 'k';
ax.GridColor = [0.75 0.75 0.75];
ax.XTick = 1:4;
ax.XTickLabel = group_labels;
ax.XTickLabelRotation = 18;
grid(ax,'on');
title(ax,panel_title,'Color','k');
ylabel(ax,y_label,'Color','k');
end
