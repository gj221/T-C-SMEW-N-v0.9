function [validation,fig,plot_file] = Plot_Kantola_Validation(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compare annual modeled grain carbon with Kantola field observations.
%% The initial daily state is excluded before calendar-year aggregation.
%%
%% Options
%%   ResultDirectory  result folder            [SIM/Subfield_12]
%%   Runs             control and basalt pair  [10 11]
%%   ObservationFile  block-level harvest CSV  [Illinois Data Bank]
%%   Visible          figure visibility        ['on']
%%   OutputFile       output PNG                [result folder]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

here = fileparts(mfilename('fullpath'));
opt = struct('ResultDirectory',fullfile(here,'SIM','Subfield_12'), ...
    'Runs',[10 11],'ObservationFile','', ...
    'Visible','on','OutputFile','');
if mod(numel(varargin),2) ~= 0
    error('Plot_Kantola_Validation:BadOptions', ...
        'Options must be name/value pairs.');
end
for k = 1:2:numel(varargin)
    name = char(varargin{k});
    if ~isfield(opt,name)
        error('Plot_Kantola_Validation:BadOptions','Unknown option %s.',name);
    end
    opt.(name) = varargin{k+1};
end
if ~ismember(char(opt.Visible),{'on','off'})
    error('Plot_Kantola_Validation:BadOptions','Visible must be on or off.');
end
if ~(ischar(opt.ResultDirectory) || ...
        (isstring(opt.ResultDirectory) && isscalar(opt.ResultDirectory))) || ...
        ~(ischar(opt.ObservationFile) || ...
        (isstring(opt.ObservationFile) && isscalar(opt.ObservationFile))) || ...
        ~(ischar(opt.OutputFile) || ...
        (isstring(opt.OutputFile) && isscalar(opt.OutputFile)))
    error('Plot_Kantola_Validation:BadOptions', ...
        'File and directory options must be character vectors or string scalars.');
end
if ~isfolder(opt.ResultDirectory)
    error('Plot_Kantola_Validation:MissingDirectory', ...
        'Result directory not found: %s',char(opt.ResultDirectory));
end

run_ids = reshape(opt.Runs,1,[]);
if ~(isnumeric(run_ids) && isreal(run_ids) && numel(run_ids) == 2 && ...
        all(isfinite(run_ids)) && all(run_ids == fix(run_ids)) && ...
        all(run_ids >= 1 & run_ids <= 12))
    error('Plot_Kantola_Validation:BadRuns', ...
        'Runs must be one control and basalt pair from the same scenario group.');
end
RM_all = Run_Matrix_12();
RM = RM_all(run_ids);
if ~strcmp(RM(1).treatment,'control') || ...
        ~strcmp(RM(2).treatment,'basalt') || ...
        ~strcmp(RM(1).layout,RM(2).layout)
    error('Plot_Kantola_Validation:BadRuns', ...
        'Runs must be one control and basalt pair from the same scenario group.');
end

if isempty(opt.ObservationFile)
    observations = Load_Kantola_Harvest_Observations();
else
    observations = Load_Kantola_Harvest_Observations(opt.ObservationFile);
end
years = unique(observations.Year,'stable');
model = nan(numel(years),numel(RM));
required = {'R','Datam','RB_L','Ccrown_t','topo','ColCfg'};
for i = 1:numel(RM)
    result_file = fullfile(char(opt.ResultDirectory), ...
        sprintf('TC_2D_run%02d_%s.mat',RM(i).id,RM(i).name));
    if ~isfile(result_file)
        error('Plot_Kantola_Validation:MissingResult', ...
            'Result not found: %s',result_file);
    end
    saved = whos('-file',result_file);
    if ~all(ismember(required,{saved.name}))
        error('Plot_Kantola_Validation:IncompleteCarbonResult', ...
            'Required harvest-carbon state is incomplete in %s.',result_file);
    end
    S = load(result_file,required{:});
    if ~isstruct(S.R) || ~isfield(S.R,'id') || ~isfield(S.R,'name') || ...
            S.R.id ~= RM(i).id || ~strcmp(S.R.name,RM(i).name)
        error('Plot_Kantola_Validation:WrongResult', ...
            'Result does not match run %d: %s',RM(i).id,result_file);
    end
    model(:,i) = annual_grain_c(S,years,result_file);
end

n_obs = height(observations);
RunID = nan(n_obs,1);
ScenarioGroup = repmat(string(RM(1).layout),n_obs,1);
ModelGrainC_gC_m2 = nan(n_obs,1);
for row = 1:n_obs
    treatment = find(string({RM.treatment}) == observations.Treatment(row));
    year_id = find(years == observations.Year(row));
    if ~isscalar(treatment) || ~isscalar(year_id)
        error('Plot_Kantola_Validation:ObservationMismatch', ...
            'Each observation must match one modeled year and treatment.');
    end
    RunID(row) = RM(treatment).id;
    ModelGrainC_gC_m2(row) = model(year_id,treatment);
end
GrainCResidual_gC_m2 = ModelGrainC_gC_m2-observations.Mean_gC_m2;
validation = table(RunID,ScenarioGroup,observations.Year, ...
    observations.Crop,observations.Treatment,ModelGrainC_gC_m2, ...
    observations.Mean_gC_m2,observations.SEM_gC_m2, ...
    GrainCResidual_gC_m2,observations.SampleSize, ...
    observations.SourceKind,observations.SourceFile, ...
    'VariableNames',{'RunID','ScenarioGroup','Year','Crop','Treatment', ...
    'ModelGrainC_gC_m2','ObservedMean_gC_m2','ObservedSEM_gC_m2', ...
    'GrainCResidual_gC_m2','ObservedSampleSize','SourceKind','SourceFile'});

treatment_names = ["control" "basalt"];
treatment_labels = {'Control','Basalt'};
treatment_colors = [0 0 0; 0 0.25 0.60];
fig = figure('Color','w','Visible',char(opt.Visible), ...
    'Name','Kantola crop validation');
try
    fig.Position(3:4) = [1200 520];
    tl = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
    ax1 = nexttile(tl);
    ax1.Tag = 'KantolaGrainC';
    [legend_handles,legend_labels] = draw_comparison(ax1,validation, ...
        years,treatment_names,treatment_labels,treatment_colors);

    ax2 = nexttile(tl);
    ax2.Tag = 'KantolaGrainCResidual';
    draw_residual(ax2,validation,years,treatment_names,treatment_colors);

    lgd = legend(legend_handles,legend_labels,'Orientation','horizontal', ...
        'Box','off','TextColor','k','Color','w','EdgeColor','none');
    lgd.Layout.Tile = 'south';
    title(tl,sprintf('Kantola crop validation: %s', ...
        strrep(RM(1).layout,'_',' ')),'Color','k');
    if isempty(opt.OutputFile)
        plot_file = fullfile(char(opt.ResultDirectory),sprintf( ...
            'Kantola_validation_%s.png',RM(1).layout));
    else
        plot_file = char(opt.OutputFile);
    end
    exportgraphics(fig,plot_file,'Resolution',200);
catch ME
    if isgraphics(fig); close(fig); end
    rethrow(ME)
end
fprintf('Saved: %s\n',plot_file);

end

function annual = annual_grain_c(S,years,result_file)
if ~isnumeric(S.RB_L) || ~isnumeric(S.Ccrown_t) || ...
        ~isnumeric(S.Datam) || ~isstruct(S.topo) || ...
        ~isfield(S.topo,'cell_area_m2') || ~isstruct(S.ColCfg)
    error('Plot_Kantola_Validation:IncompleteCarbonResult', ...
        'Harvest-carbon state has invalid types in %s.',result_file);
end

n_time = size(S.RB_L,1);
cc = size(S.RB_L,2);
area = S.topo.cell_area_m2(:);
n_col = numel(area);
if n_time < 2 || cc < 1 || n_col < 1 || numel(S.ColCfg) ~= n_col || ...
        any(~isfinite(area)) || any(area <= 0)
    error('Plot_Kantola_Validation:IncompleteCarbonResult', ...
        'Harvest-carbon dimensions or cell metadata are invalid in %s.', ...
        result_file);
end
if size(S.RB_L,3) ~= 7 || size(S.RB_L,4) ~= n_col || ...
        numel(S.RB_L) ~= n_time*cc*7*n_col || ...
        size(S.Ccrown_t,1) ~= n_time || size(S.Ccrown_t,2) ~= cc || ...
        size(S.Ccrown_t,3) ~= n_col || ...
        numel(S.Ccrown_t) ~= n_time*cc*n_col
    error('Plot_Kantola_Validation:CarbonTimeLengthMismatch', ...
        'Harvest-carbon arrays do not share the daily dimensions in %s.', ...
        result_file);
end
RB_L = reshape(S.RB_L,n_time,cc,7,n_col);
Ccrown_t = reshape(S.Ccrown_t,n_time,cc,n_col);
if any(~isfinite(RB_L(:))) || any(~isfinite(Ccrown_t(:))) || ...
        any(Ccrown_t(:) < 0 | Ccrown_t(:) > 1) || ...
        any(sum(Ccrown_t(2:end,:,:),2) > 1+1e-9,'all') || ...
        any(RB_L(:,:,5,:) < -1e-9,'all')
    error('Plot_Kantola_Validation:UnphysicalCarbonState', ...
        'Grain removal or crown fraction is unphysical in %s.',result_file);
end

if size(S.Datam,2) < 4 || any(~isfinite(S.Datam(:)))
    error('Plot_Kantola_Validation:IncompleteCarbonResult', ...
        'Hourly forcing dates are incomplete in %s.',result_file);
end
daily_rows = find(S.Datam(:,4) == 1);
if numel(daily_rows) ~= n_time-1
    error('Plot_Kantola_Validation:CarbonTimeLengthMismatch', ...
        'Daily process rows do not match hourly forcing dates in %s.',result_file);
end
try
    process_dates = datetime(S.Datam(daily_rows,1), ...
        S.Datam(daily_rows,2),S.Datam(daily_rows,3));
catch
    error('Plot_Kantola_Validation:IncompleteCarbonResult', ...
        'Hourly forcing dates are invalid in %s.',result_file);
end
if any(diff(process_dates) ~= days(1))
    error('Plot_Kantola_Validation:CarbonTimeLengthMismatch', ...
        'Daily process dates are missing, duplicated or out of order in %s.', ...
        result_file);
end
model_years = year(process_dates);
if ~all(ismember(years,model_years))
    error('Plot_Kantola_Validation:ObservationOutsideModelPeriod', ...
        'Observation years fall outside the saved model period in %s.', ...
        result_file);
end
for k = 1:numel(years)
    dates_in_year = process_dates(model_years == years(k));
    if dates_in_year(1) ~= datetime(years(k),1,1) || ...
            dates_in_year(end) < datetime(years(k),12,30)
        error('Plot_Kantola_Validation:IncompleteModelYear', ...
            'The saved result does not cover the full harvest year %d in %s.', ...
            years(k),result_file);
    end
end

crop_cell = false(n_col,1);
for c = 1:n_col
    if ~isfield(S.ColCfg(c),'PFT_OVERRIDE') || ...
            ~(ischar(S.ColCfg(c).PFT_OVERRIDE) || ...
            (isstring(S.ColCfg(c).PFT_OVERRIDE) && ...
            isscalar(S.ColCfg(c).PFT_OVERRIDE)))
        error('Plot_Kantola_Validation:IncompleteCarbonResult', ...
            'Crop configuration is incomplete in %s.',result_file);
    end
    crop_cell(c) = isempty(char(S.ColCfg(c).PFT_OVERRIDE));
end
if ~any(crop_cell)
    error('Plot_Kantola_Validation:UnphysicalCarbonState', ...
        'The result does not contain a positive cropped area in %s.',result_file);
end

daily_cell = zeros(n_time-1,n_col);
for c = find(crop_cell).'
    for k = 1:cc
        fruit = max(0,reshape(RB_L(2:end,k,5,c),n_time-1,1));
        crown = reshape(Ccrown_t(2:end,k,c),n_time-1,1);
        daily_cell(:,c) = daily_cell(:,c) + crown.*fruit;
    end
end
crop_area = sum(area(crop_cell));
daily = daily_cell(:,crop_cell)*area(crop_cell)/crop_area;
annual = zeros(numel(years),1);
for k = 1:numel(years)
    annual(k) = sum(daily(year(process_dates) == years(k)));
end

end

function [handles,labels] = draw_comparison(ax,T,years,treatment_names, ...
        treatment_labels,treatment_colors)
hold(ax,'on');
set(ax,'Color','w','XColor','k','YColor','k');
handles = gobjects(4,1);
labels = cell(4,1);
offsets = [-0.08 0.08];
for treatment = 1:2
    rows = T.Treatment == treatment_names(treatment);
    [plot_years,order] = sort(T.Year(rows));
    x = plot_years+offsets(treatment);
    model_value = T.ModelGrainC_gC_m2(rows);
    observed_mean = T.ObservedMean_gC_m2(rows);
    observed_sem = T.ObservedSEM_gC_m2(rows);
    handles(2*treatment-1) = plot(ax,x,model_value(order),'s-', ...
        'Color',treatment_colors(treatment,:),'LineWidth',1.4, ...
        'MarkerFaceColor',treatment_colors(treatment,:),'Tag','ModelGrainC');
    handles(2*treatment) = errorbar(ax,x,observed_mean(order), ...
        observed_sem(order),'o','Color',treatment_colors(treatment,:), ...
        'LineWidth',1.2,'MarkerFaceColor','w','CapSize',7, ...
        'Tag','ObservedGrainC');
    labels{2*treatment-1} = ['Model ' treatment_labels{treatment}];
    labels{2*treatment} = ['Observed ' treatment_labels{treatment}];
end
xlabel(ax,'Harvest year','Color','k');
ylabel(ax,'g C m^{-2}','Color','k');
sample_size = unique(T.ObservedSampleSize);
if isscalar(sample_size)
    panel_title = sprintf('Annual harvested grain C (raw blocks, n = %d)', ...
        sample_size);
else
    panel_title = 'Annual harvested grain C (raw block data)';
end
title(ax,panel_title,'Color','k');
ax.XTick = years;
ax.XLim = [min(years)-0.5 max(years)+0.5];
box(ax,'on');
end

function draw_residual(ax,T,years,treatment_names,treatment_colors)
hold(ax,'on');
set(ax,'Color','w','XColor','k','YColor','k');
residual = nan(numel(years),2);
for treatment = 1:2
    for year_id = 1:numel(years)
        row = T.Year == years(year_id) & ...
            T.Treatment == treatment_names(treatment);
        if nnz(row) ~= 1
            error('Plot_Kantola_Validation:ObservationMismatch', ...
                'Each year must contain one control and one basalt observation.');
        end
        residual(year_id,treatment) = T.GrainCResidual_gC_m2(row);
    end
end
bars = bar(ax,years,residual,'grouped');
for treatment = 1:2
    bars(treatment).FaceColor = treatment_colors(treatment,:);
    bars(treatment).Tag = 'GrainCResidual';
end
yline(ax,0,'k-');
xlabel(ax,'Harvest year','Color','k');
ylabel(ax,'Model - observed (g C m^{-2})','Color','k');
title(ax,'Harvest-C residual','Color','k');
ax.XTick = years;
ax.XLim = [min(years)-0.5 max(years)+0.5];
box(ax,'on');
end
