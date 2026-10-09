function summary = Summarize_Uncertainty(sources,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Collect the three parameter sets of the Kantola uncertainty ensemble
%% (Uncertainty_Members: the original, the lower and the upper set) into
%% one summary file, per run, as field (domain) means on the daily axis,
%% with the yearly statistics, the member envelopes and the observations.
%%
%%   Summarize_Uncertainty({original_dir, ensemble_dir})
%%   Summarize_Uncertainty({orig_1, ens_1; orig_2, ens_2}, 'Runs',1:9, ...
%%       'Output',file)
%%
%% Each row of sources is one campaign: original_dir holds the original
%% set's results (TC_2D_runNN_*.mat of a Base, SSP or field SSP campaign)
%% and ensemble_dir the lower/ and upper/ directories of its ensemble. Rows
%% may hold different runs on different grids (the controls on the control
%% field, the treatments on the basalt field); a run must come from one row
%% only, and every result must share one daily axis. Options:
%%   'Runs'    runs to look for (default 1:12); runs no row holds are dropped
%%   'Output'  the summary file (default uncertainty_summary.mat in the
%%             first ensemble_dir); a yearly and a ranking CSV go beside it
%%   'Label'   what the summary covers, for the figures' titles (default:
%%             Kantola 1-D, 2-D plot or field grids, with the scenario
%%             and case the campaign tags name)
%% A missing member result is NaN throughout and listed.
%%
%% Every result is checked against its set (Uncertainty_Result_Matches: the
%% knobs per square, the horizon, the scenario and calcite case in the file
%% name), and the three sets of a run against each other (squares, grid
%% file and its SHA-256, start date).
%%
%% WHAT IS READ. A multi-square result of package y keeps no hourly or
%% ten-minute records: daily water content and drainage come from its
%% daily stores (O_d, Lk_d), the feedstock release from its ERW day
%% summaries (ERW_daily, the minerals' EW totals with the stoichiometry of
%% ERW{c}.const), and nothing is needed from Final_state, which holds only
%% the last days and hours of what the stores already cover. A result
%% with full records (one column, or a package before y) is read from the
%% hourly O and Lk and the ten-minute weath.EW instead, binned into the
%% same rows. The pools come from P, the fluxes and pH from Domain (the
%% frame's field means), NEE from the crown-weighted NPP and respiration
%% stores; per-square quantities are averaged with Domain's weights
%% (Field_Aggregation_Areas: the polygon overlap on a field grid, the cell
%% area on the plot).
%%
%% Daily rows are dated as Plot_1D_Then_All_12 dates them: row j holds the
%% day from 01:00 of start+(j-3), row 2 is the lone first hour dated the day
%% before the start, and row 1 the initial state. Calendar-year statistics
%% use rows 3 onward (complete days).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'runtime','Uncertainty_1D'),fullfile(here,'runtime','Base_1D'), ...
    fullfile(here,'v2','run'),fullfile(here,'v2','T&C_Code'),fullfile(here,'v2','plotting'));
if ischar(sources) || isstring(sources) || ~iscell(sources) || size(sources,2) ~= 2 || ...
        isempty(sources) || ~all(cellfun(@(x) ischar(x) || (isstring(x) && isscalar(x)),sources(:)))
    error('Summarize_Uncertainty:BadSources', ...
        'sources must be a cell array of {original_dir, ensemble_dir} rows.');
end
sources = cellfun(@(x) regexprep(char(x),'[/]+$',''),sources,'UniformOutput',false);
opt = struct('Runs',1:12,'Output','','Label','');
if mod(numel(varargin),2) ~= 0
    error('Summarize_Uncertainty:BadOptions','Options must be name/value pairs.');
end
for k = 1:2:numel(varargin)
    name = char(varargin{k});
    if ~isfield(opt,name)
        error('Summarize_Uncertainty:BadOptions','Unknown option %s.',name);
    end
    opt.(name) = varargin{k+1};
end
if isempty(opt.Output)
    opt.Output = fullfile(sources{1,2},'uncertainty_summary.mat');
end

P = Campaign_Parameters;
members = Uncertainty_Members;
n_members = numel(members);
RM = Run_Matrix_12;

%% Which row holds each run, and the campaign each row's tags name.
run_ids = opt.Runs(:).';
row_of_run = zeros(size(run_ids));
for r = 1:numel(run_ids)
    for s = 1:size(sources,1)
        if ~isempty(find_results(sources{s,1},run_ids(r))) || ...
                ~isempty(find_results(fullfile(sources{s,2},'lower'),run_ids(r))) || ...
                ~isempty(find_results(fullfile(sources{s,2},'upper'),run_ids(r)))
            if row_of_run(r) ~= 0
                error('Summarize_Uncertainty:RunInTwoSources', ...
                    'Run %d is held by more than one row of sources.',run_ids(r));
            end
            row_of_run(r) = s;
        end
    end
end
if ~any(row_of_run)
    error('Summarize_Uncertainty:NoResults','No result for runs %s.',mat2str(run_ids));
end
if any(row_of_run == 0)
    fprintf('Runs without any result, dropped: %s\n',mat2str(run_ids(row_of_run == 0)));
end
run_ids = run_ids(row_of_run > 0);
row_of_run = row_of_run(row_of_run > 0);
n_runs = numel(run_ids);
treatments = arrayfun(@(id) RM([RM.id] == id).treatment,run_ids,'UniformOutput',false);
layouts = arrayfun(@(id) RM([RM.id] == id).layout,run_ids,'UniformOutput',false);

daily_fields = {'CaSolution','MgSolution','SiSolution', ...
    'CumulativeGaseousNLoss','WFPS','Drainage_mm','CumulativeN2O','CumulativeNO3Leaching', ...
    'NO3Pool','BulkPH','PoreWaterPH','NEE','DailyN2O_gN_ha_d', ...
    'CaRelease_mol_m2','MgRelease_mol_m2'};
summary_fields = {'CumulativeN2O_gN_m2','CumulativeNO3Leaching_gN_m2', ...
    'FinalNO3Pool_gN_m2','FinalBulkPH','CumulativeNEE_gC_m2'};

daily = struct(); Date = []; n_days = 0; first_date = [];
summary_value = struct();
result_file = cell(n_members,n_runs);
grid = repmat(struct('n_columns',NaN,'grid',''),1,n_runs);
missing = {};
for r = 1:n_runs
    s = row_of_run(r);
    for m = 1:n_members
        if members(m).reused
            folder = sources{s,1};
        else
            folder = fullfile(sources{s,2},members(m).tag);
        end
        [~,campaign] = fileparts(sources{s,1});
        if ~members(m).reused
            [~,campaign] = fileparts(sources{s,2});
        end
        files = find_results(folder,run_ids(r));
        if isempty(files)
            missing{end+1} = sprintf('%s run %d (%s)',members(m).tag,run_ids(r),folder); %#ok<AGROW>
            continue
        end
        if numel(files) > 1
            error('Summarize_Uncertainty:AmbiguousResult', ...
                '%d results for run %d in %s.',numel(files),run_ids(r),folder);
        end
        file = fullfile(files(1).folder,files(1).name);
        S = load(file,'run_provenance','Datam','Domain','topo');
        check_member(S.run_provenance,files(1).name,members(m),P,campaign,sources(s,:),file);
        nd = numel(S.Domain.N2Oflx);
        start = datenum(S.Datam(1,1),S.Datam(1,2),S.Datam(1,3));
        if isempty(Date)
            n_days = nd; first_date = start;
            Date = start + (0:nd-1).' - 2;
            for f = daily_fields
                daily.(f{1}) = nan(n_days,n_members,n_runs);
            end
            for f = summary_fields
                summary_value.(f{1}) = nan(n_members,n_runs);
            end
        elseif nd ~= n_days || start ~= first_date
            error('Summarize_Uncertainty:TimeAxisMismatch', ...
                '%s is not on the summary''s daily axis (%d days from %s).', ...
                file,n_days,datestr(first_date,'yyyy-mm-dd'));
        end
        this_grid = grid_identity(S.run_provenance,S.topo);
        if isnan(grid(r).n_columns)
            grid(r) = this_grid;
        elseif ~isequal(grid(r),this_grid)
            error('Summarize_Uncertainty:GridMismatch', ...
                'The sets of run %d did not run on the same squares: %s.',run_ids(r),file);
        end
        result_file{m,r} = file;
        d = result_daily(file,S,nd);
        for f = daily_fields
            daily.(f{1})(:,m,r) = d.(f{1});
        end
        summary_value.CumulativeN2O_gN_m2(m,r) = sum(S.Domain.N2Oflx(:));
        summary_value.CumulativeNO3Leaching_gN_m2(m,r) = sum(S.Domain.LEAK_NO3(:));
        summary_value.FinalNO3Pool_gN_m2(m,r) = d.NO3Pool(end);
        summary_value.FinalBulkPH(m,r) = d.BulkPH(end);
        summary_value.CumulativeNEE_gC_m2(m,r) = sum(d.NEE,'omitnan');
        fprintf('  %-8s run %2d  %s (%s)\n',members(m).tag,run_ids(r),files(1).name,d.layout);
    end
end
if isempty(Date)
    error('Summarize_Uncertainty:NoResults','No result was readable.');
end
if ~isempty(missing)
    fprintf('Missing results (NaN in the summary): %s\n',strjoin(missing,', '));
end

%% The three sets span the ensemble: every member is in the envelope, the
%% original (member 1) is its reference.
in_envelope = true(1,n_members);
envelope = struct();
for f = daily_fields
    x = daily.(f{1});
    envelope.(f{1}) = struct('min',reshape(min(x,[],2,'omitnan'),n_days,n_runs), ...
        'max',reshape(max(x,[],2,'omitnan'),n_days,n_runs), ...
        'base',reshape(x(:,1,:),n_days,n_runs));
end
ctl = find(strcmp(treatments,'control'),1); bas = find(strcmp(treatments,'basalt'),1);
if ~isempty(ctl) && ~isempty(bas)
    sig = daily.BulkPH(:,:,bas)-daily.BulkPH(:,:,ctl);
    envelope.signal_bulk_pH = struct('min',min(sig,[],2,'omitnan'), ...
        'max',max(sig,[],2,'omitnan'),'base',sig(:,1));
end
summary_range = struct();
for f = summary_fields
    x = summary_value.(f{1});
    summary_range.(f{1}) = struct('min',min(x,[],1,'omitnan'), ...
        'max',max(x,[],1,'omitnan'),'base',x(1,:),'members',x);
end

%% Calendar years of complete days (rows 3 onward).
dv = datevec(Date);
complete = (1:n_days).' >= 3;
years = unique(dv(complete,1)).';
n_years = numel(years);
yearly = struct('years',years);
yearly.pH_season_mean = nan(n_members,n_years,n_runs);
yearly.N2O_cumulative_gN_m2 = nan(n_members,n_years,n_runs);
yearly.NO3_leaching_cumulative_gN_m2 = nan(n_members,n_years,n_runs);
yearly.Ca_release_mol_m2 = nan(n_members,n_years,n_runs);
yearly.Mg_release_mol_m2 = nan(n_members,n_years,n_runs);
yearly.drainage_mm = nan(n_members,n_years,n_runs);
daily_n2o = daily.DailyN2O_gN_ha_d/1e4;
daily_no3 = [daily.CumulativeNO3Leaching(1,:,:); diff(daily.CumulativeNO3Leaching,1,1)];
for y = 1:n_years
    in_year = complete & dv(:,1) == years(y);
    season = in_year & dv(:,2) >= 5 & dv(:,2) <= 9;
    last = find(in_year,1,'last');
    yearly.pH_season_mean(:,y,:) = mean(daily.BulkPH(season,:,:),1);
    yearly.N2O_cumulative_gN_m2(:,y,:) = sum(daily_n2o(in_year,:,:),1);
    yearly.NO3_leaching_cumulative_gN_m2(:,y,:) = sum(daily_no3(in_year,:,:),1);
    yearly.drainage_mm(:,y,:) = sum(daily.Drainage_mm(in_year,:,:),1);
    yearly.Ca_release_mol_m2(:,y,:) = daily.CaRelease_mol_m2(last,:,:);
    yearly.Mg_release_mol_m2(:,y,:) = daily.MgRelease_mol_m2(last,:,:);
end
%% Field-paper CDR potential: 2 x (Ca + Mg) released, as CO2, t/ha.
yearly.CDR_potential_tCO2_ha = 2*(yearly.Ca_release_mol_m2+yearly.Mg_release_mol_m2)*44e-6*1e4;
if ~isempty(ctl) && ~isempty(bas)
    yearly.pH_season_signal = yearly.pH_season_mean(:,:,bas)-yearly.pH_season_mean(:,:,ctl);
end

%% Which set gives each edge of the spread in the last year.
ranked = {'pH_season_mean','N2O_cumulative_gN_m2','NO3_leaching_cumulative_gN_m2', ...
    'drainage_mm','Ca_release_mol_m2','Mg_release_mol_m2','CDR_potential_tCO2_ha'};
ranking = repmat(struct('metric','','run_id',0,'treatment','','year',0, ...
    'min_member','','max_member','','original',NaN,'min',NaN,'max',NaN),0,1);
for f = ranked
    x = yearly.(f{1});
    for r = 1:n_runs
        v = x(:,end,r);
        if all(isnan(v)); continue; end
        [vmin,imin] = min(v); [vmax,imax] = max(v);
        ranking(end+1,1) = struct('metric',f{1},'run_id',run_ids(r), ...
            'treatment',treatments{r},'year',years(end),'min_member',members(imin).tag, ...
            'max_member',members(imax).tag,'original',v(1),'min',vmin,'max',vmax); %#ok<AGROW>
    end
end

member_rows = struct('id',{members.id},'tag',{members.tag}, ...
    'parameter','all sixteen','level',{members.tag},'envelope',true, ...
    'KNmaxd',[],'AppeloAlpha',[],'result_file',[]);
for m = 1:n_members
    K = Uncertainty_Knobs(members(m),P);
    member_rows(m).KNmaxd = K.KNmaxd;
    member_rows(m).AppeloAlpha = K.AppeloAlpha;
    member_rows(m).result_file = result_file(m,:);
end
obs = observations(here);
label = char(opt.Label);
if isempty(label)
    label = summary_label(sources,grid);
end
summary = struct('label',label,'sources',{sources},'members',member_rows,'in_envelope',in_envelope, ...
    'run_ids',run_ids,'treatments',{treatments},'layouts',{layouts},'grid',grid, ...
    'daily_fields',{daily_fields},'Date',Date,'daily',daily, ...
    'envelope',envelope,'summary_fields',{summary_fields}, ...
    'summary_value',summary_value,'summary_range',summary_range, ...
    'yearly',yearly,'ranking',ranking,'obs',obs,'missing',{missing});
output_dir = fileparts(opt.Output);
if ~isempty(output_dir) && ~isfolder(output_dir); mkdir(output_dir); end
save(opt.Output,'-struct','summary','-v7.3');
write_tables(opt.Output,summary);
fprintf('Summary: %d sets x %d runs over %d days -> %s\n',n_members,n_runs,n_days,opt.Output);
for k = 1:numel(ranking)
    fprintf('  %-30s run %d %-8s %d: original %.4g, min %.4g (%s), max %.4g (%s)\n', ...
        ranking(k).metric,ranking(k).run_id,ranking(k).treatment,ranking(k).year, ...
        ranking(k).original,ranking(k).min,ranking(k).min_member,ranking(k).max, ...
        ranking(k).max_member);
end
end

%% =======================================================================
function files = find_results(folder,run_id)
files = dir(fullfile(folder,sprintf('TC_2D_run%02d_*.mat',run_id)));
files = files(~[files.isdir]);
end

%% =======================================================================
function check_member(p,name,member,P,campaign,source,file)
%% The result was run as this set, under its campaign: Uncertainty_Result_Matches
%% on the result's own horizon, squares and rock schedule, with the scenario
%% and calcite case the campaign's tag names. The validation stamp is not
%% compared: each campaign has its own.
K = Uncertainty_Knobs(member,P);
scenario = char(regexp(campaign,'(?<=_)ssp[0-9]{3}(?=_|$)','match','once'));
case_token = char(regexp(campaign,'_case[0-9]+(?=__|$)','match','once'));
calcite_case = [];
if ~isempty(case_token); calcite_case = str2double(case_token(6:end)); end
fields = {'n_days','n_columns','last_application_year'};
if ~all(isfield(p,fields))
    error('Summarize_Uncertainty:Provenance','%s carries no run provenance.',file);
end
C = struct('scenario',scenario,'calcite_case',calcite_case,'n_days',p.n_days, ...
    'n_columns',p.n_columns,'last_application_year',p.last_application_year);
if ~Uncertainty_Result_Matches(p,name,C,K,campaign,'')
    error('Summarize_Uncertainty:Provenance', ...
        '%s was not run as set %s of %s (sources %s, %s).',file,member.tag,campaign, ...
        source{1},source{2});
end
end

%% =======================================================================
function g = grid_identity(p,topo)
g = struct('n_columns',p.n_columns,'grid','');
if isfield(p,'field_topology_file_sha256')
    g.grid = p.field_topology_file_sha256;
elseif isstruct(topo) && isfield(topo,'cell_area_m2')
    g.grid = sprintf('cells:%s',mat2str(topo.cell_area_m2(:).',12));
end
end

%% =======================================================================
function label = summary_label(sources,grid)
%% Kantola 1-D, 2-D plot or field grids, and the scenario and calcite case
%% the first campaign's tag names.
if all([grid.n_columns] == 1)
    label = 'Kantola 1-D';
elseif any(cellfun(@(g) ~isempty(g) && ~startsWith(g,'cells:'),{grid.grid}))
    label = 'Kantola field grids';
else
    label = 'Kantola 2-D plot';
end
[~,campaign] = fileparts(sources{1,1});
scenario = char(regexp(campaign,'(?<=_)ssp[0-9]{3}(?=_|$)','match','once'));
case_token = char(regexp(campaign,'(?<=_)case[0-9]+(?=__|$)','match','once'));
if ~isempty(scenario); label = [label ', ' scenario]; end
if ~isempty(case_token); label = [label ', ' case_token]; end
end

%% =======================================================================
function d = result_daily(file,S,nd)
%% The daily series of one result as field (domain) means.
dom = S.Domain;
[~,w] = Field_Aggregation_Areas(S.topo);
w = w(:);
n_col = numel(w);
contents = whos('-file',file);
stores = all(ismember({'O_d','Lk_d','ERW_daily'},{contents.name}));
d.layout = 'daily stores';
if ~stores; d.layout = 'full records'; end
L = load(file,'PAR','NPP_H','NPP_L','Ccrown_t','R_litter','R_microbe','R_ew','ERW');
m = matfile(file);
if numel(size(m,'P')) == 2
    pools = m.P(:,60:8:76);
else
    pools = m.P(:,60:8:76,:);
end
pools = reshape(double(pools),nd,3,n_col);
d.CaSolution = reshape(pools(:,1,:),nd,n_col)*w;
d.MgSolution = reshape(pools(:,2,:),nd,n_col)*w;
d.SiSolution = reshape(pools(:,3,:),nd,n_col)*w;
d.CumulativeGaseousNLoss = cumsum(dom.N2flx(:)+dom.N2Oflx(:));
d.CumulativeN2O = cumsum(dom.N2Oflx(:));
d.CumulativeNO3Leaching = cumsum(dom.LEAK_NO3(:));
d.NO3Pool = dom.NO3_pool(:);
d.DailyN2O_gN_ha_d = 1e4*dom.N2Oflx(:);
if isfield(dom,'PH_bulk_d') && isfield(dom,'PH_pore_d')
    d.BulkPH = dom.PH_bulk_d(:); d.PoreWaterPH = dom.PH_pore_d(:);
else
    d.BulkPH = dom.PH_d(:); d.PoreWaterPH = dom.PH_d(:);
end
d.NEE = net_ecosystem_exchange(L,nd,n_col)*w;

%% Water of the biogeochemical layer and drainage, per square.
wfps = nan(nd,n_col);
drainage = nan(nd,n_col);
if stores
    H = load(file,'O_d','Lk_d');
    O_d = reshape(double(H.O_d),nd,[],n_col);
    Lk_d = reshape(double(H.Lk_d),nd,n_col);
    n_layer = size(O_d,2);
    for c = 1:n_col
        bio = double(L.PAR{c}.Bio_Zs(:).');
        osat = double(L.PAR{c}.Soil_Param.Osat(:).');
        wfps(:,c) = 100*O_d(:,:,c)*bio(1:n_layer).'/sum(bio.*osat);
    end
    drainage = Lk_d;
else
    H = load(file,'O','Lk');
    O = double(H.O); Lk = double(H.Lk);
    for c = 1:n_col
        wfps(:,c) = Daily_WFPS_From_Hourly(O(:,:,c),L.PAR{c},nd);
        %% Lk is a rate [mm h-1] on the hourly step (dth = 1 h).
        drainage(:,c) = hourly_to_daily(Lk(:,1,c),nd);
    end
end
d.WFPS = wfps*w;
d.Drainage_mm = drainage*w;

%% Feedstock Ca and Mg release [mol m-2], cumulative on the daily axis.
ca = zeros(nd,n_col); mg = zeros(nd,n_col);
if stores
    D = load(file,'ERW_daily');
    for c = 1:n_col
        [ca(:,c),mg(:,c)] = release_from_days(D.ERW_daily{c},L.ERW{c},nd);
    end
else
    for c = 1:n_col
        [ca(:,c),mg(:,c)] = release_from_rows(L.ERW{c},nd);
    end
end
d.CaRelease_mol_m2 = ca*w;
d.MgRelease_mol_m2 = mg*w;
end

%% =======================================================================
function nee = net_ecosystem_exchange(L,nd,n_col)
%% Plot_1D_Then_All_12's carbon_dynamics per square: heterotrophic
%% respiration (litter + microbial + weathering zone) minus crown-weighted
%% NPP, positive to the atmosphere.
cc = numel(L.NPP_L)/(nd*n_col);
NPP_H = reshape(double(L.NPP_H),nd,cc,n_col);
NPP_L = reshape(double(L.NPP_L),nd,cc,n_col);
Ccrown_t = reshape(double(L.Ccrown_t),nd,cc,n_col);
respiration = reshape(double(L.R_litter)+double(L.R_microbe)+double(L.R_ew),nd,n_col);
nee = respiration-reshape(sum(Ccrown_t.*(NPP_H+NPP_L),2),nd,n_col);
end

%% =======================================================================
function x = hourly_to_daily(hourly,nd)
%% Hourly amounts summed into the daily rows BIOGEO fills: row 1 none,
%% row 2 the first hour, row j the 24 hours ending at its 01:00.
x = nan(nd,1);
x(2) = hourly(1);
x(3:nd) = sum(reshape(hourly(2:24*(nd-2)+1),24,nd-2),1).';
end

%% =======================================================================
function [ca,mg] = release_from_days(D,E,nd)
%% Each ERW day's mineral dissolution [umol mineral m-2] times the release
%% stoichiometry, on the daily row the day was stored with.
ca = zeros(nd,1); mg = zeros(nd,1);
if isempty(D) || ~isfield(D,'total') || ~isfield(D.total,'EW') || D.n == 0
    return
end
st = stoichiometry(E,size(D.total.EW,2));
if isempty(st); return; end
row = min(max(round(D.j(:)),1),nd);
ca = cumsum(accumarray(row,D.total.EW*st(:,1),[nd 1]))*1e-6;
mg = cumsum(accumarray(row,D.total.EW*st(:,2),[nd 1]))*1e-6;
end

%% =======================================================================
function [ca,mg] = release_from_rows(E,nd)
%% The ten-minute EW [umol mineral m-2 d-1] of a full record, times dtw,
%% binned as the ERW days run: rows 2-6 the first call (daily row 2), then
%% 144 rows from 01:00 to 00:50 per daily row, the flush on the last row.
ca = zeros(nd,1); mg = zeros(nd,1);
if ~isfield(E,'weath') || ~isfield(E.weath,'EW') || isempty(E.weath.EW)
    return
end
EW = E.weath.EW;
st = stoichiometry(E,size(EW,2));
if isempty(st); return; end
steps_per_hour = round(1/(24*E.dtw));
steps_per_day = 24*steps_per_hour;
n_rows = size(EW,1);
rows = (1:n_rows).';
row = 3 + floor((rows-(steps_per_hour+1))/steps_per_day);
row(rows <= steps_per_hour) = 2;
row(1) = 1;
row = min(max(row,1),nd);
ca = cumsum(accumarray(row,EW*st(:,1)*E.dtw,[nd 1]))*1e-6;
mg = cumsum(accumarray(row,EW*st(:,2)*E.dtw,[nd 1]))*1e-6;
end

%% =======================================================================
function st = stoichiometry(E,n_minerals)
st = [];
if ~isfield(E,'const') || ~isfield(E.const,'min_st') || isempty(E.const.min_st)
    return
end
st = E.const.min_st;
if size(st,1) ~= n_minerals
    error('Summarize_Uncertainty:Stoichiometry', ...
        'min_st (%d rows) does not match EW (%d minerals).',size(st,1),n_minerals);
end
end

%% =======================================================================
function obs = observations(root)
%% Observations with their standard errors, on the field record's years.
obs = struct('years',2016:2020);
csv = fullfile(root,'v2','Experimental observations','Web plot digitizer', ...
    'Beerling2024_Fig1E_vector_extraction.csv');
if ~isfile(csv)
    warning('Summarize_Uncertainty:NoObservations','No pH observations at %s.',csv);
else
    opts = detectImportOptions(csv,'CommentStyle','#','Delimiter',',');
    opts = setvartype(opts,{'treatment','depth_cm','period'},'string');
    opts = setvartype(opts,{'pH','SEM'},'double');
    T = readtable(csv,opts);
    periods = ["pre" "2017" "2018" "2019" "2020"];
    for treatment = {'control','basalt'}
        for depth = {'0-10','10-30'}
            rows = T.treatment == treatment{1} & T.depth_cm == depth{1};
            value = nan(1,5); sem = nan(1,5);
            for k = 1:5
                pick = find(rows & T.period == periods(k));
                if isscalar(pick)
                    value(k) = T.pH(pick); sem(k) = T.SEM(pick);
                end
            end
            obs.pH.(sprintf('%s_%s',treatment{1},strrep(depth{1},'-','_'))) = ...
                struct('value',value,'sem',sem);
        end
    end
end
%% Beerling et al. 2024 Dataset S06 as the legacy benchmark carried it:
%% cumulative potential CDR (t CO2/ha) with its SEM, and the weathered Ca
%% and Mg (ug per g of soil-basalt mix, mean and SD/sqrt(n)) converted to
%% mmol m-2 with the mix mass of each year's cumulative application.
obs.CDR_potential_tCO2_ha = struct('value',[0 3.82 6.75 11.49 10.51], ...
    'sem',[0 1.8261 2.5109 3.2609 3.7174]);
mix_g = [0, (1:4)*5000./[0.023595168 0.046102539 0.067595642 0.088141529]];
ca_ug_g = [0 429.93 1159.85 1791.11 1218.10];
mg_ug_g = [0 237.25 156.23 343.32 540.11];
ca_sem_ug_g = [0 1464.62/20 1351.67/20 1334.53/21 1954.46/16];
mg_sem_ug_g = [0 545.48/20 1241.83/20 1624.39/21 1391.67/16];
obs.Ca_release_mmol_m2 = struct('value',ca_ug_g.*mix_g*1e-3/40, ...
    'sem',ca_sem_ug_g.*mix_g*1e-3/40);
obs.Mg_release_mmol_m2 = struct('value',mg_ug_g.*mix_g*1e-3/24, ...
    'sem',mg_sem_ug_g.*mix_g*1e-3/24);
end

%% =======================================================================
function write_tables(output,s)
[folder,stem] = fileparts(output);
n = numel(s.members)*numel(s.run_ids)*numel(s.yearly.years);
rows = cell(n,13);
k = 0;
for m = 1:numel(s.members)
    for r = 1:numel(s.run_ids)
        for y = 1:numel(s.yearly.years)
            k = k+1;
            rows(k,:) = {s.members(m).tag,s.run_ids(r),s.treatments{r},s.layouts{r}, ...
                s.yearly.years(y),s.yearly.pH_season_mean(m,y,r), ...
                s.yearly.N2O_cumulative_gN_m2(m,y,r),s.yearly.NO3_leaching_cumulative_gN_m2(m,y,r), ...
                s.yearly.drainage_mm(m,y,r),s.yearly.Ca_release_mol_m2(m,y,r), ...
                s.yearly.Mg_release_mol_m2(m,y,r),s.yearly.CDR_potential_tCO2_ha(m,y,r), ...
                s.grid(r).n_columns};
        end
    end
end
T = cell2table(rows,'VariableNames',{'member','run_id','treatment','layout','year', ...
    'pH_season_mean','N2O_gN_m2','NO3_leaching_gN_m2','drainage_mm', ...
    'Ca_release_mol_m2','Mg_release_mol_m2','CDR_potential_tCO2_ha','squares'});
writetable(T,fullfile(folder,[strrep(stem,'summary','yearly') '.csv']));
if ~isempty(s.ranking)
    writetable(struct2table(s.ranking),fullfile(folder,[strrep(stem,'summary','ranking') '.csv']));
end
end
