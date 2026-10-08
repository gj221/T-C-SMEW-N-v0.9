function test_field_cell_maps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The field figures (v2/plotting: Plot_Field_Layout, Plot_Field_Maps and
%% their helpers) read each cell's series as the model saves it, reduce
%% them over time windows correctly, weigh field means as Domain does, and
%% write the layout, the maps, the CSV and the animation, for the 12-cell
%% plot grid and for a polygon grid with a DEM and a projection.
%%
%% Synthetic results with known values stand in for model runs: N2O grows
%% with the cell number, NEE is -0.2 g C m-2 every day, pH is constant per
%% cell, the NO3 pool rises by row, and so on, so every window value has a
%% closed form.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

v2 = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(v2,'plotting'),fullfile(v2,'T&C_Code'));

out_dir = tempname;
mkdir(out_dir);
before = findall(groot,'Type','figure');
cleanup = onCleanup(@() tidy(out_dir,before));
fig_dir = fullfile(out_dir,'maps');

%% --- the 12-cell plot: prairie margin with basalt, crop interior -------
topo = Build_Cell_Topology(60,5,2,2);
N = topo.N_cells;
is_main = strcmp(topo.cell_region,'main');
ColCfg = repmat(struct('PARAM_IC',fullfile('Parameters','MOD_PARAM_Kantola_Maize_Soy.m'), ...
    'PFT_OVERRIDE','','apply_rock',false),N,1);
for c = find(~is_main(:)).'
    ColCfg(c).PFT_OVERRIDE = 'PFT_Prairie_Margin';
    ColCfg(c).apply_rock = true;
end
R = struct('name','basalt__margin_prairie','treatment','basalt');
n_day = 400;
file = fullfile(out_dir,'plot_result.mat');
write_result(file,topo,ColCfg,R,n_day,true);

res = Field_Cell_Result(file);
assert(res.N_cells == 12 && res.NNd == n_day+1 && res.NN == 24*n_day, ...
    'The result dimensions are read from the file.');
assert(res.Date_d(3) == datetime(2016,11,8) && res.day(2) == datetime(2016,11,7), ...
    'Rows are dated as Plot_1D_Then_All_12 dates them.');
assert(isequal(res.pft(is_main),repmat({'Crop (Kantola Maize Soy)'},4,1)) && ...
    isequal(res.pft(~is_main),repmat({'Prairie Margin'},8,1)), ...
    'Vegetation labels come from PARAM_IC and PFT_OVERRIDE.');
assert(isequal(res.rock(~is_main),repmat({'Basalt'},8,1)) && ...
    all(strcmp(res.rock(is_main),'None')),'Feedstock labels follow apply_rock.');
assert(abs(sum(res.weight)-1) < 1e-12,'Weights are the area shares.');

%% --- windows ------------------------------------------------------------
[win,tag] = Field_Cell_Windows(res.day,'year');
assert(strcmp(tag,'year') && numel(win) == 2 && isequal({win.label},{'2016','2017'}), ...
    'Calendar years from 8 Nov 2016 to 11 Dec 2017.');
assert(win(1).first_row == 2 && win(1).last_row == 56 && ...
    win(2).first_row == 57 && win(2).last_row == n_day+1 && win(1).prev_row == 1, ...
    'Window rows exclude the initial state and split at the new year.');
assert(numel(Field_Cell_Windows(res.day,'season')) == 6, 'Six seasons.');
assert(numel(Field_Cell_Windows(res.day,'month')) == 14, 'Fourteen months.');
assert(numel(Field_Cell_Windows(res.day,30)) == 14, 'Fourteen 30-day blocks.');
gs = Field_Cell_Windows(res.day,'growing_season');
assert(isscalar(gs) && strcmp(gs.label,'May-Oct 2017'),'One growing season.');
snap = Field_Cell_Windows(res.day,'day','Snapshots',datetime(2017,[3 8],1));
assert(numel(snap) == 2 && all([snap.n_days] == 1),'Snapshots are single days.');
clipped = Field_Cell_Windows(res.day,'month','Period',{'2017-01-01','2017-06-30'});
assert(numel(clipped) == 6 && strcmp(clipped(1).label,'Jan 2017'),'Period clips the months.');
n_rows = arrayfun(@(w) sum(w.rows),win);

%% --- maps, with every value checked -----------------------------------
keys = {'N2O','NO3_leaching','NEE','pH_bulk','NO3_pool','WFPS','crop_harvest', ...
    'LAI','drainage','lateral_water'};
out = Plot_Field_Maps(file,'Variables',keys,'Interval','year', ...
    'FigureDirectory',fig_dir,'Visible','off','TraceCells',[1 5]);
value = @(key) out.variables(strcmp({out.variables.key},key)).values;
cell_id = 1:N;
assert(max(abs(value('N2O') - 10*1e-3*n_rows(:)*cell_id),[],'all') < 1e-12, ...
    'N2O is the window total in kg N ha-1.');
expected = zeros(2,N); expected(:,1) = 10*2e-3*n_rows(:);
assert(max(abs(value('NO3_leaching')-expected),[],'all') < 1e-12,'NO3 leaching totals.');
assert(max(abs(value('NEE') - (-0.2)*n_rows(:)*ones(1,N)),[],'all') < 1e-10, ...
    'NEE is respiration minus crown-weighted NPP.');
ph = 6+0.02*cell_id;
assert(max(abs(value('pH_bulk') - repmat(ph,2,1)),[],'all') < 1e-12,'Constant pH is kept.');
fm = out.variables(strcmp({out.variables.key},'pH_bulk')).field_mean;
assert(max(abs(fm - pH_Area_Weighted(ph,res.weight))) < 1e-12, ...
    'The pH field mean weighs hydrogen-ion activity by area, as Domain does.');
assert(max(abs(value('NO3_pool') - 10*0.1*[win.last_row].'*cell_id),[],'all') < 1e-10, ...
    'Pools are read at the window end.');
wfps = (100*(0.2+0.01*cell_id)+200*0.3)/150*100;
assert(max(abs(value('WFPS') - repmat(wfps,2,1)),[],'all') < 1e-10, ...
    'WFPS is the depth-weighted saturation of the biogeochemical layer.');
series = Field_Cell_Series(res,'WFPS');
L = load(file,'O','PAR');
for c = [1 7 12]
    reference = Daily_WFPS_From_Hourly(L.O(:,:,c),L.PAR{c},res.NNd);
    assert(isequal(isnan(series(:,c)),isnan(reference)) && ...
        max(abs(series(2:end,c)-reference(2:end))) < 1e-10, ...
        'Daily WFPS matches Daily_WFPS_From_Hourly row for row.');
end
harvest = zeros(2,N); harvest(2,is_main) = 1.5;
assert(isequal(value('crop_harvest'),harvest), ...
    'Harvest falls in 2017 on the crop cells only.');
assert(max(abs(value('LAI') - 1),[],'all') < 1e-12,'LAI is crown-weighted.');
drain = [0.01+0.24*(n_rows(1)-1); 0.24*n_rows(2)];
assert(max(abs(value('drainage') - drain*ones(1,N)),[],'all') < 1e-9, ...
    'Hourly drainage rates sum into the daily rows BIOGEO fills.');
lateral = zeros(2,N); lateral(:,1) = n_rows(:);
lateral(:,2) = -n_rows(:)*topo.cell_area_m2(1)/topo.cell_area_m2(2);
assert(max(abs(value('lateral_water')-lateral),[],'all') < 1e-9, ...
    'The lateral ledger is summed per window.');
nee = out.variables(strcmp({out.variables.key},'NEE'));
assert(nee.diverging && nee.clim(1) == -nee.clim(2),'NEE has a scale centred on zero.');
for k = 1:numel(keys)
    assert(isfile(fullfile(fig_dir,[keys{k} '_year.png'])) && ...
        isfile(fullfile(fig_dir,[keys{k} '_year.fig'])),'Each variable writes PNG and FIG.');
end
assert(numel(findobj(out.figures(1),'Type','axes')) == 3, ...
    'A variable figure holds the series and one map per window.');
T = readtable(out.csv);
assert(height(T) == numel(keys)*2*N && all(ismember({'zone','vegetation', ...
    'feedstock','value','field_mean','window_start'},T.Properties.VariableNames)), ...
    'The CSV holds every value with the cell attributes.');
assert(isfile(out.mat),'The values are saved as MAT.');

%% --- statistic override, reference, anomaly, combined, animation ------
o_end = Plot_Field_Maps(file,'Variables','N2O','Interval','year','Statistic','end', ...
    'Export',false,'Visible','off','TimeSeries',false);
assert(max(abs(o_end.variables.values - 10*1e-3*ones(2,1)*cell_id),[],'all') < 1e-12, ...
    'Statistic end reads the last day of the window.');
o_ref = Plot_Field_Maps(file,'Variables',{'N2O','pH_bulk'},'Interval','year', ...
    'Reference',file,'FigureDirectory',fig_dir,'Visible','off');
assert(all([o_ref.variables.values] == 0,'all') && ...
    isequal(o_ref.variables(1).clim,[-1 1]), ...
    'A run minus itself is zero on a centred scale.');
assert(all(ismember({'run_value','reference_value'},o_ref.table.Properties.VariableNames)), ...
    'The reference values are kept.');
o_anom = Plot_Field_Maps(file,'Variables','N2O','Interval','year','Anomaly',true, ...
    'Export',false,'Visible','off');
assert(max(abs(o_anom.variables.values*res.weight(:))) < 1e-12, ...
    'Anomalies average to zero over the field.');
Plot_Field_Maps(file,'Variables',{'N2O','pH_bulk'},'Interval','year', ...
    'Layout','combined','FigureDirectory',fig_dir,'Visible','off');
assert(isfile(fullfile(fig_dir,'maps_year.png')),'The combined figure is written.');
o_gif = Plot_Field_Maps(file,'Variables','pH_bulk','Interval','season', ...
    'Static',false,'Animation',true,'FigureDirectory',fig_dir,'Visible','off');
assert(isfile(o_gif.variables.gif) && numel(imfinfo(o_gif.variables.gif)) == 6, ...
    'The animation has a frame per season.');
assert_error(@() Plot_Field_Maps(file,'Variables','N2O','Interval','day', ...
    'Export',false,'Visible','off'),'Plot_Field_Maps:TooManyWindows');
assert_error(@() Field_Cell_Variables('nonsense'),'Field_Cell_Variables:UnknownVariable');
plain = fullfile(out_dir,'plain_result.mat');
write_result(plain,topo,ColCfg,R,10,false);
assert_error(@() Field_Cell_Series(Field_Cell_Result(plain),'lateral_water'), ...
    'Field_Cell_Series:NoLateralLedger');

%% --- the layout of the plot result -------------------------------------
[fig_l,files] = Plot_Field_Layout(file,'FigureDirectory',fig_dir,'Visible','off');
assert(isfile(files.png) && numel(findobj(fig_l,'Type','axes')) == 3, ...
    'A flat run shows zones, vegetation and feedstock.');

%% --- a polygon grid with a DEM and a projection ------------------------
field = polyshape([0 31 34 25 4 -2],[0 -1 14 24 22 9]);
[grid_polygon,~,offset] = Place_Field_Polygon_Grid(field,5);
ptopo = Build_Polygon_Grid_Topology(grid_polygon,5,'majority',offset);
ptopo.cell_region(ptopo.outer_face_width_m > 0) = {'buffer'};
xs = floor(min(ptopo.cell_x_lo))-10:2.5:ceil(max(ptopo.cell_x_hi))+10;
ys = floor(min(ptopo.cell_y_lo))-10:2.5:ceil(max(ptopo.cell_y_hi))+10;
[gx,gy] = meshgrid(xs,ys);
Z = 100 + 0.05*(max(ys)-gy) + 0.3*sin(gx/8);
ptopo = Attach_Topography_To_Topology(ptopo,struct('Z',Z,'x',xs,'y',ys));
G = Field_Cell_Geometry(ptopo,'Coordinates','local');
drawn = 0;
for k = 1:size(G.faces,1)
    f = G.faces(k,isfinite(G.faces(k,:)));
    drawn = drawn + polyarea(G.vertices(f,1),G.vertices(f,2));
end
assert(abs(drawn-sum(ptopo.cell_overlap_area_m2)) < 1e-6*sum(ptopo.cell_overlap_area_m2), ...
    'The clipped faces cover exactly the in-field area the model weighs.');
ptopo.grid_rotation_degrees = 30;
ptopo.projected_grid_origin_m = [397000 4435000];
Gp = Field_Cell_Geometry(ptopo);
[px,py] = Grid_Local_To_Projected(ptopo.cell_xc(:),ptopo.cell_yc(:),ptopo);
assert(strcmp(Gp.coordinates,'projected') && isequal(Gp.north,[0 1]) && ...
    max(abs(Gp.centre_x(:)+Gp.offset(1)-px)) < 1e-6 && ...
    max(abs(Gp.centre_y(:)+Gp.offset(2)-py)) < 1e-6, ...
    'Projected maps place the cells where the grid lies, north up.');
Gl = Field_Cell_Geometry(ptopo,'Coordinates','local');
assert(max(abs(Gl.north-[sind(30) cosd(30)])) < 1e-12,'Local maps carry the north arrow.');
[fig_p,files_p] = Plot_Field_Layout(ptopo,'FigureDirectory',fig_dir, ...
    'Name','polygon_layout','Visible','off');
assert(isfile(files_p.png) && numel(findobj(fig_p,'Type','axes')) >= 3 && ...
    ~isempty(findobj(fig_p,'Tag','DEMContours')) && ...
    ~isempty(findobj(fig_p,'Tag','OverlandFlow')), ...
    'A DEM grid shows elevation with contours, flow paths and slope.');

pColCfg = repmat(struct('PARAM_IC','MOD_PARAM_Kantola_Maize_Soy.m', ...
    'PFT_OVERRIDE','','apply_rock',false),ptopo.N_cells,1);
pfile = fullfile(out_dir,'polygon_result.mat');
write_result(pfile,ptopo,pColCfg,R,40,true);
o_poly = Plot_Field_Maps(pfile,'Variables',{'N2O','lateral_water'},'Interval','month', ...
    'FigureDirectory',fig_dir,'Visible','off');
pres = Field_Cell_Result(pfile);
n2o = o_poly.variables(1);
assert(isequal(size(n2o.values),[2 ptopo.N_cells]) && ...
    max(abs(n2o.field_mean - n2o.values*pres.weight(:))) < 1e-12 && ...
    abs(sum(pres.area_m2)-sum(ptopo.cell_overlap_area_m2)) < 1e-9, ...
    'Field means on a polygon grid weigh cells by their overlap with the field.');
fprintf('test_field_cell_maps: PASS\n');
end

%% =======================================================================
function write_result(file,topo,ColCfg,R,n_day,with_ledger) %#ok<INUSL>
N = topo.N_cells;
NN = 24*n_day; NNd = n_day+1; ms = 4; cc = 1; dth = 1; %#ok<NASGU>
Datam = datevec(datetime(2016,11,8,1,0,0) + hours(0:NN-1).'); %#ok<NASGU>
cell_id = reshape(1:N,1,1,N);
row = (1:NNd).';
N2Oflx = 1e-3*(row > 1).*cell_id;
N2Oflx_nit = 0.25*N2Oflx; %#ok<NASGU>
LEAK_NO3 = 2e-3*(row > 1).*(cell_id == 1); %#ok<NASGU>
PH_bulk_d = repmat(6+0.02*cell_id,NNd,1);
PH_pore_d = PH_bulk_d-0.3; %#ok<NASGU>
P = zeros(NNd,79,N);
P(:,32,:) = row.*(0.1*cell_id); %#ok<NASGU>
NPP_H = zeros(NNd,cc,N); NPP_L = 2*ones(NNd,cc,N); Ccrown_t = 0.5*ones(NNd,cc,N); %#ok<NASGU>
LAI_H = ones(NNd,cc,N); LAI_L = ones(NNd,cc,N); %#ok<NASGU>
R_litter = 0.5*ones(NNd,1,N); R_microbe = 0.25*ones(NNd,1,N); R_ew = 0.05*ones(NNd,1,N); %#ok<NASGU>
RB_L = zeros(NNd,cc,7,N);
if NNd >= 200; RB_L(200,1,5,:) = 3; end %#ok<NASGU>
O = 0.3*ones(NN,ms,N);
O(:,1,:) = repmat(0.2+0.01*cell_id,NN,1); %#ok<NASGU>
PAR = repmat({struct('Bio_Zs',[100 200 0 0],'Soil_Param',struct('Osat',[0.5 0.5 0.45 0.45]))},N,1); %#ok<NASGU>
Lk = 0.01*ones(NN,1,N); q_runon = zeros(NN,1,N); Rh = zeros(NN,1,N); %#ok<NASGU>
Rd = zeros(NN,1,N); ZWT = 1500*ones(NN,1,N); f = zeros(NN,1,N); %#ok<NASGU>
names = {'topo','ColCfg','R','Datam','NN','NNd','dth','N2Oflx','N2Oflx_nit', ...
    'LEAK_NO3','PH_bulk_d','PH_pore_d','P','NPP_H','NPP_L','Ccrown_t','LAI_H', ...
    'LAI_L','R_litter','R_microbe','R_ew','RB_L','O','PAR','Lk','q_runon','Rh', ...
    'Rd','ZWT','f'};
if with_ledger
    water = zeros(NNd,N);
    water(2:end,1) = 1;
    water(2:end,2) = -topo.cell_area_m2(1)/topo.cell_area_m2(2);
    lateral_diag = struct('water_cell_net_mm',water, ...
        'NO3_cell_net_g_m2',1e-3*water); %#ok<NASGU>
    names{end+1} = 'lateral_diag';
end
save(file,names{:},'-v7.3');
end

function assert_error(fcn,identifier)
try
    fcn();
catch ME
    assert(strcmp(ME.identifier,identifier), ...
        'Expected %s, got %s: %s',identifier,ME.identifier,ME.message);
    return
end
error('test_field_cell_maps:NoError','Expected the error %s.',identifier);
end

function tidy(out_dir,before)
opened = setdiff(findall(groot,'Type','figure'),before);
if ~isempty(opened); close(opened); end
if isfolder(out_dir); rmdir(out_dir,'s'); end
end
