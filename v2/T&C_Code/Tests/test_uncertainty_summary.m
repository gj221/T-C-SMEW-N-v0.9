function test_uncertainty_summary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Summarize_Uncertainty (package root) on synthetic results: the three
%% sets of a control run on one grid and a basalt run on another, 400 days
%% across a year boundary, written once with package y's daily stores
%% (O_d, Lk_d, ERW_daily) and once with full hourly and ten-minute records.
%% Both layouts must give the same summary; the field means must be the
%% Domain-weighted means; the yearly statistics, envelopes and checks must
%% hold; and Plot_Uncertainty must draw the summary.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_dir = fileparts(mfilename('fullpath'));
code_dir = fileparts(test_dir);
repo_root = fileparts(fileparts(code_dir));
addpath(repo_root,fullfile(repo_root,'runtime','Uncertainty_1D'), ...
    fullfile(repo_root,'runtime','Base_1D'),fullfile(repo_root,'v2','run'), ...
    code_dir,fullfile(repo_root,'v2','plotting'));
P = Campaign_Parameters;
members = Uncertainty_Members;
root = tempname;
mkdir(root);
cleanup = onCleanup(@() tidy(root));

n_days = 400;
grids = struct('name',{'grid_a','grid_b'},'n_col',{3,2},'run',{1,2}, ...
    'run_name',{'control__whole_field','basalt__whole_field'});
for layout = {'stores','full'}
    for g = 1:numel(grids)
        original = fullfile(root,layout{1},['Field_SSP_2D_35yr_ssp245__' grids(g).name]);
        ensemble = fullfile(root,layout{1},['Field_Uncertainty_2D_35yr_ssp245__' grids(g).name]);
        for m = 1:numel(members)
            if members(m).reused
                folder = original;
            else
                folder = fullfile(ensemble,members(m).tag);
            end
            [~,campaign] = fileparts(original);
            if ~members(m).reused; [~,campaign] = fileparts(ensemble); end
            write_result(folder,grids(g),members(m),Uncertainty_Knobs(members(m),P), ...
                campaign,n_days,strcmp(layout{1},'stores'),10*g+m);
        end
    end
end
sources = @(layout) { ...
    fullfile(root,layout,'Field_SSP_2D_35yr_ssp245__grid_a'), ...
    fullfile(root,layout,'Field_Uncertainty_2D_35yr_ssp245__grid_a'); ...
    fullfile(root,layout,'Field_SSP_2D_35yr_ssp245__grid_b'), ...
    fullfile(root,layout,'Field_Uncertainty_2D_35yr_ssp245__grid_b')};
out_stores = fullfile(root,'stores_summary.mat');
out_full = fullfile(root,'full_summary.mat');
S = Summarize_Uncertainty(sources('stores'),'Runs',1:3,'Output',out_stores);
F = Summarize_Uncertainty(sources('full'),'Runs',1:3,'Output',out_full);

%% Shape and the three sets.
assert(isequal(S.run_ids,[1 2]) && isequal(S.treatments,{'control','basalt'}) && ...
    isequal(S.in_envelope,true(1,3)) && isequal({S.members.tag},{'original','lower','upper'}) && ...
    isequal([S.grid.n_columns],[3 2]),'Runs 1 and 2 come from their own grids, with three sets each.');
assert(strcmp(S.label,'Kantola field grids, ssp245'),'The summary names what it covers.');
assert(numel(S.Date) == n_days+1 && S.Date(3) == datenum(2016,1,1) && ...
    S.Date(2) == datenum(2015,12,31),'Rows are dated as Plot_1D_Then_All_12 dates them.');
assert(isfile(out_stores) && isfile(fullfile(root,'stores_yearly.csv')) && ...
    isfile(fullfile(root,'stores_ranking.csv')),'The summary and its tables are written.');

%% The two layouts give the same summary.
for f = S.daily_fields
    a = S.daily.(f{1}); b = F.daily.(f{1});
    assert(isequal(isnan(a),isnan(b)) && max(abs(a(:)-b(:)),[],'omitnan') <= ...
        1e-12*max(1,max(abs(a(:)),[],'omitnan')), ...
        'The daily %s differs between daily stores and full records.',f{1});
end
for f = fieldnames(S.yearly).'
    if strcmp(f{1},'years'); continue; end
    a = S.yearly.(f{1}); b = F.yearly.(f{1});
    assert(max(abs(a(:)-b(:)),[],'omitnan') <= 1e-12*max(1,max(abs(a(:)),[],'omitnan')), ...
        'The yearly %s differs between the layouts.',f{1});
end

%% The field means are Domain's weights applied to each square.
X = load(result_path(fullfile(root,'stores','Field_SSP_2D_35yr_ssp245__grid_a'),1));
[~,w] = Field_Aggregation_Areas(X.topo);
ca = reshape(X.P(:,60,:),n_days+1,[])*w(:);
assert(max(abs(S.daily.CaSolution(:,1,1)-ca)) < 1e-12,'CaSolution is the weighted P(:,60).');
nee = reshape(X.R_litter+X.R_microbe+X.R_ew,n_days+1,[]) - ...
    reshape(sum(X.Ccrown_t.*(X.NPP_H+X.NPP_L),2),n_days+1,[]);
assert(max(abs(S.daily.NEE(:,1,1)-nee*w(:))) < 1e-12,'NEE is the weighted crown balance.');
bio = X.PAR{1}.Bio_Zs; osat = X.PAR{1}.Soil_Param.Osat;
wfps = zeros(n_days+1,X.topo.N_cells);
for c = 1:X.topo.N_cells
    wfps(:,c) = 100*X.O_d(:,:,c)*bio(1:3).'/sum(bio.*osat);
end
assert(isequaln(S.daily.WFPS(:,1,1),wfps*w(:)),'WFPS is the weighted store of water content.');
assert(isequaln(S.daily.Drainage_mm(:,1,1),reshape(X.Lk_d,n_days+1,[])*w(:)), ...
    'Drainage is the weighted daily store.');
assert(isequal(S.daily.CumulativeN2O(:,1,1),cumsum(X.Domain.N2Oflx)) && ...
    isequal(S.daily.BulkPH(:,1,1),X.Domain.PH_bulk_d),'Fluxes and pH are Domain''s.');
%% Feedstock release: each ERW day's EW totals times the stoichiometry.
D = X.ERW_daily{1};
ca_square = cumsum(accumarray(D.j,D.total.EW*X.ERW{1}.const.min_st(:,1),[n_days+1 1]))*1e-6;
B = load(result_path(fullfile(root,'stores','Field_SSP_2D_35yr_ssp245__grid_b'),2));
[~,wb] = Field_Aggregation_Areas(B.topo);
release = zeros(n_days+1,B.topo.N_cells);
for c = 1:B.topo.N_cells
    Dc = B.ERW_daily{c};
    release(:,c) = cumsum(accumarray(Dc.j,Dc.total.EW*B.ERW{c}.const.min_st(:,1),[n_days+1 1]))*1e-6;
end
assert(max(abs(S.daily.CaRelease_mol_m2(:,1,2)-release*wb(:))) <= 1e-12*max(release(:)) && ...
    any(ca_square > 0), ...
    'The release is the weighted cumulative ERW-day dissolution.');

%% Yearly statistics on complete days.
assert(isequal(S.yearly.years,[2016 2017]),'Calendar years start with the first complete day.');
dv = datevec(S.Date);
season = (1:n_days+1).' >= 3 & dv(:,1) == 2016 & dv(:,2) >= 5 & dv(:,2) <= 9;
in_2016 = (1:n_days+1).' >= 3 & dv(:,1) == 2016;
assert(abs(S.yearly.pH_season_mean(2,1,2)-mean(S.daily.BulkPH(season,2,2))) < 1e-12 && ...
    abs(S.yearly.N2O_cumulative_gN_m2(3,1,1)-sum(S.daily.DailyN2O_gN_ha_d(in_2016,3,1))/1e4) < 1e-12 && ...
    abs(S.yearly.Ca_release_mol_m2(1,1,2)-S.daily.CaRelease_mol_m2(find(in_2016,1,'last'),1,2)) < 1e-15, ...
    'Season pH, calendar N2O and year-end release follow their definitions.');
assert(max(abs(S.yearly.CDR_potential_tCO2_ha(:)-2*(S.yearly.Ca_release_mol_m2(:)+ ...
    S.yearly.Mg_release_mol_m2(:))*44e-6*1e4)) < 1e-12,'CDR potential is 2 x (Ca + Mg) as CO2.');
assert(isequal(size(S.yearly.pH_season_signal),[3 2]) && ...
    max(abs(S.yearly.pH_season_signal(:)-reshape(S.yearly.pH_season_mean(:,:,2)- ...
    S.yearly.pH_season_mean(:,:,1),[],1))) < 1e-12,'The pH signal is basalt minus control.');
x = S.daily.BulkPH;
assert(isequaln(S.envelope.BulkPH.min,reshape(min(x,[],2),n_days+1,2)) && ...
    isequaln(S.envelope.BulkPH.base,reshape(x(:,1,:),n_days+1,2)), ...
    'The envelope spans the three sets with the original as its base.');

%% The checks: a set run with the wrong knobs, a missing set, a run held
%% twice, a result on another axis.
lower_a = fullfile(root,'stores','Field_Uncertainty_2D_35yr_ssp245__grid_a','lower');
bad = load(result_path(lower_a,1));
bad.run_provenance.KNmaxd = P.KNmaxd;
bad_root = fullfile(root,'bad');
copyfile(fullfile(root,'stores'),bad_root);
save(result_path(fullfile(bad_root,'Field_Uncertainty_2D_35yr_ssp245__grid_a','lower'),1), ...
    '-struct','bad','-v7.3');
bad_sources = strrep(sources('stores'),fullfile(root,'stores'),bad_root);
assert_error(@() Summarize_Uncertainty(bad_sources,'Runs',1:2,'Output',fullfile(root,'bad.mat')), ...
    'Summarize_Uncertainty:Provenance');
delete(result_path(fullfile(bad_root,'Field_Uncertainty_2D_35yr_ssp245__grid_a','lower'),1));
M = Summarize_Uncertainty(bad_sources,'Runs',1:2,'Output',fullfile(root,'missing.mat'));
assert(all(isnan(M.daily.BulkPH(:,2,1))) && ~any(isnan(M.daily.BulkPH(3:end,1,1))) && ...
    isscalar(M.missing) && contains(M.missing{1},'lower run 1'), ...
    'A missing set is NaN throughout and listed.');
twice = [sources('stores'); sources('stores')];
assert_error(@() Summarize_Uncertainty(twice,'Runs',1,'Output',fullfile(root,'twice.mat')), ...
    'Summarize_Uncertainty:RunInTwoSources');
short = fullfile(root,'short','Field_SSP_2D_35yr_ssp245__grid_c');
write_result(short,struct('name','grid_c','n_col',2,'run',3,'run_name','calcite__whole_field'), ...
    members(1),Uncertainty_Knobs(members(1),P),'Field_SSP_2D_35yr_ssp245__grid_c',n_days-10,true,99);
mixed = [sources('stores'); {short,fullfile(root,'short','Field_Uncertainty_2D_35yr_ssp245__grid_c')}];
assert_error(@() Summarize_Uncertainty(mixed,'Runs',1:3,'Output',fullfile(root,'mixed.mat')), ...
    'Summarize_Uncertainty:TimeAxisMismatch');
assert_error(@() Summarize_Uncertainty('not a cell'),'Summarize_Uncertainty:BadSources');

%% The figures of the summary.
figure_dir = fullfile(root,'figures');
[fig,plot_file] = Plot_Uncertainty(out_stores,figure_dir,'off');
assert(numel(fig) == 3 && all(cellfun(@isfile,plot_file)), ...
    'Plot_Uncertainty must draw the three-set summary.');
layout = findall(fig(1),'Type','tiledlayout');
assert(isscalar(layout) && startsWith(layout.Title.String,'Kantola field grids, ssp245:'), ...
    'The figures take their title from the summary''s label.');
close(fig);
fprintf('PASS: the uncertainty summary reads daily stores and full records alike.\n');
end

%% =======================================================================
function write_result(folder,grid,member,K,campaign,n_days,stores,seed)
%% One synthetic result as Run_One_Subfield saves it, with the fields the
%% summary reads, in the daily-store or the full-record layout.
rng(seed);
if ~isfolder(folder); mkdir(folder); end
n_col = grid.n_col;
NN = 24*n_days; NNd = n_days+1;
dates = datenum(2016,1,1) + (0:NN-1).'/24;
[YE,MO,DA,HO] = datevec(dates);
Datam = [YE MO DA HO];
topo = struct('N_cells',n_col,'cell_area_m2',100*ones(n_col,1), ...
    'cell_overlap_area_m2',100*(1:n_col).'/n_col);
w = topo.cell_overlap_area_m2/sum(topo.cell_overlap_area_m2);
PAR = repmat({struct('Bio_Zs',[100 100 100 0], ...
    'Soil_Param',struct('Osat',[0.5 0.45 0.4 0.4]))},n_col,1);
P = zeros(NNd,79,n_col);
P(:,[32 60 68 76],:) = rand(NNd,4,n_col);
cc = 2;
NPP_H = rand(NNd,cc,n_col); NPP_L = rand(NNd,cc,n_col); Ccrown_t = rand(NNd,cc,n_col);
R_litter = rand(NNd,1,n_col); R_microbe = rand(NNd,1,n_col); R_ew = rand(NNd,1,n_col);
cell_daily = @(x) reshape(x,NNd,n_col)*w;
Domain = struct();
Domain.N2Oflx = cell_daily(rand(NNd,1,n_col)*1e-3);
Domain.N2flx = cell_daily(rand(NNd,1,n_col)*1e-3);
Domain.LEAK_NO3 = cell_daily(rand(NNd,1,n_col)*1e-2);
Domain.NO3_pool = reshape(P(:,32,:),NNd,n_col)*w;
Domain.PH_bulk_d = 6+rand(NNd,1);
Domain.PH_pore_d = 6.5+rand(NNd,1);
O = 0.2+0.2*rand(NN,4,n_col);
Lk = 0.01*rand(NN,1,n_col);
dtw = 1/144;
n_t = 6*NN;
ERW = cell(n_col,1);
ERW_daily = cell(n_col,1);
for c = 1:n_col
    EW = 1e3*rand(n_t,2);
    EW(1,:) = 0;
    ERW{c} = struct('dtw',dtw,'const',struct('min_st',[1 0.5 0 0 0 1 0; 0 1 0 0 0 1 0]), ...
        'weath',struct('EW',EW));
    %% The ERW days as ERW_Daily_Step_2D tiles the grid: rows 2-6 on daily
    %% row 2, then 144 rows per daily row, the flush on the last.
    first = [2; 7+144*(0:n_days-1).'];
    last = [6; 6+144*(1:n_days-1).'; n_t];
    j = [2; (3:NNd).'; NNd];
    total = zeros(numel(first),2);
    for k = 1:numel(first)
        total(k,:) = sum(EW(first(k):last(k),:),1)*dtw;
    end
    ERW_daily{c} = struct('n',numel(first),'rows',[first last],'j',j, ...
        'total',struct('EW',total));
end
run_provenance = struct('campaign',campaign,'validation_fingerprint','synthetic', ...
    'n_days',n_days,'n_columns',n_col,'last_application_year',2050, ...
    'KNmaxd',K.KNmaxd,'Kbot',K.Kbot,'biogeo_overrides',K.BiogeoOverrides, ...
    'dissolution_factor',K.DissolutionFactor,'saturation_exponent',K.SaturationExponent, ...
    'AppeloAlpha',repmat(K.AppeloAlpha,n_col,1),'AppeloLogKH',repmat(K.AppeloLogKH,n_col,1), ...
    'field_topology_file_sha256',grid.name);
vars = {'run_provenance','Datam','Domain','topo','PAR','P','NPP_H','NPP_L','Ccrown_t', ...
    'R_litter','R_microbe','R_ew','ERW'};
if stores
    %% The frame's daily stores: the hours each daily step averaged.
    O_d = nan(NNd,3,n_col); Lk_d = nan(NNd,1,n_col);
    O_d(2,:,:) = O(1,1:3,:); Lk_d(2,:,:) = Lk(1,:,:);
    for jd = 3:NNd
        hours = 24*(jd-3)+2:24*(jd-2)+1;
        O_d(jd,:,:) = mean(O(hours,1:3,:),1);
        Lk_d(jd,:,:) = sum(Lk(hours,:,:),1);
    end
    for c = 1:n_col
        ERW{c} = rmfield(ERW{c},'weath');
    end
    vars = [vars {'O_d','Lk_d','ERW_daily'}];
else
    vars = [vars {'O','Lk'}];
end
name = sprintf('TC_2D_run%02d_%s__ssp245__field_%s__%s__%dd__from20160101.mat', ...
    grid.run,grid.run_name,grid.name,member.tag,n_days);
save(fullfile(folder,name),vars{:},'-v7.3');
end

%% =======================================================================
function file = result_path(folder,run_id)
files = dir(fullfile(folder,sprintf('TC_2D_run%02d_*.mat',run_id)));
assert(isscalar(files),'Expected one result for run %d in %s.',run_id,folder);
file = fullfile(files.folder,files.name);
end

%% =======================================================================
function tidy(root)
close(findall(groot,'Type','figure'));
if isfolder(root); rmdir(root,'s'); end
end

%% =======================================================================
function assert_error(fun,identifier)
try
    fun();
catch ME
    assert(strcmp(ME.identifier,identifier), ...
        'Expected %s, got %s: %s',identifier,ME.identifier,ME.message);
    return
end
error('Expected error %s was not raised.',identifier);
end
