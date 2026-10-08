function spec = Field_Cell_Variables(key)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Field_Cell_Variables
%%
%% The per-cell quantities Plot_Field_Maps can draw, one row per key.
%%
%%   spec = Field_Cell_Variables            every row (struct array)
%%   spec = Field_Cell_Variables('N2O')     one row (keys ignore case)
%%   Field_Cell_Variables('list')           print the keys
%%
%% Fields of a row
%%   key          name passed to Plot_Field_Maps('Variables',...)
%%   label        title text (TeX subscripts, no backslashes)
%%   units        units of a WINDOW value, after scale
%%   source       'daily'    saved [NNd x 1 x N_col] daily series
%%                'pool'     columns of the daily pool matrix P, summed
%%                'hourly'   saved [NN x 1 x N_col] step series, or the
%%                           daily store <name>_d a multi-square run of
%%                           package y keeps in its place
%%                'lateral'  per-cell daily ledger in lateral_diag
%%                'derived'  computed in Field_Cell_Series
%%   name         saved variable, or the derived quantity
%%   index        P columns ('pool')
%%   hourly       'sum' or 'mean': how the steps of a day combine
%%   rate         true when the step series is a rate per hour (x dth)
%%   stat         the default window statistic: 'sum' (fluxes), 'mean',
%%                'hmean' (mean hydrogen-ion activity, for pH), 'end',
%%                'start', 'change' (end minus start), 'max', 'min'
%%   scale        multiplies every value drawn (10: g m-2 -> kg ha-1)
%%   diverging    zero is the natural centre of the colour scale
%%   description  one line for the list
%%
%% A struct with these fields can be passed to Plot_Field_Maps in place of
%% a key, to map a saved series this table does not list; missing fields
%% take the defaults of a daily flux (sum, scale 1).
%%
%% Units as saved: N, C and element fluxes g m-2 d-1, pools g m-2 (P),
%% water mm per step (Rh, Rd) or mm h-1 (Lk, q_runon, f), WTD mm.
%% NEE follows the T&C sign: positive is carbon released to the
%% atmosphere. R_ew is earthworm respiration (BIOGEOCHEMISTRY_DYNAMIC3),
%% part of the heterotrophic respiration, not a weathering flux.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent rows
if isempty(rows)
    rows = build_rows();
end
if nargin < 1 || isempty(key)
    spec = rows;
    return
end
if isstruct(key)
    spec = complete_row(key);
    return
end
if ~(ischar(key) || (isstring(key) && isscalar(key)))
    error('Field_Cell_Variables:BadKey', ...
        'A variable is a key (text) or a struct row.');
end
key = char(key);
if strcmpi(key,'list')
    print_rows(rows);
    if nargout > 0; spec = rows; end
    return
end
k = find(strcmpi({rows.key},key),1);
if isempty(k)
    error('Field_Cell_Variables:UnknownVariable', ...
        'Unknown variable "%s". Known keys: %s.',key,strjoin({rows.key},', '));
end
spec = rows(k);
end

%% =======================================================================
function rows = build_rows()
kgN = 'kg N ha^{-1}';
r = {};
%% Nitrogen gases, transformations and inputs [g N m-2 d-1]
r{end+1} = row('N2O','N_2O emission',kgN,'daily','N2Oflx', ...
    'Soil N2O emission, nitrification and denitrification');
r{end+1} = row('N2O_nitrification','N_2O from nitrification',kgN,'daily','N2Oflx_nit', ...
    'The nitrification share of the N2O emission');
r{end+1} = row('N2','N_2 emission',kgN,'daily','N2flx','Denitrification N2');
r{end+1} = row('NH3','NH_3 volatilisation',kgN,'daily','VOL','Ammonia volatilisation');
r{end+1} = row('nitrification','Nitrification',kgN,'daily','NITflx','Gross nitrification');
r{end+1} = row('denitrification','Denitrification',kgN,'daily','DENITflx','Gross denitrification');
r{end+1} = row('N_mineralisation','Net N mineralisation',kgN,'daily','Min_N', ...
    'Net mineralisation (negative is immobilisation)');
r{end+1} = row('N_fixation','Biological N fixation',kgN,'daily','BfixN','Symbiotic and free N fixation');
r{end+1} = row('N_fertiliser','Fertiliser N',kgN,'daily','FertN_d','Fertiliser N applied');
r{end+1} = row('N_deposition','N deposition',kgN,'daily','DepN_d','Atmospheric N deposition');
%% Leaching at the base of the biogeochemical layer [g m-2 d-1]
r{end+1} = row('NO3_leaching','NO_3^- leaching',kgN,'daily','LEAK_NO3','Nitrate leached');
r{end+1} = row('NH4_leaching','NH_4^+ leaching',kgN,'daily','LEAK_NH4','Ammonium leached');
r{end+1} = row('DON_leaching','DON leaching',kgN,'daily','LEAK_DON','Dissolved organic N leached');
r{end+1} = row('DOC_leaching','DOC leaching','kg C ha^{-1}','daily','LEAK_DOC','Dissolved organic C leached');
r{end+1} = row('Ca_leaching','Ca leaching','kg Ca ha^{-1}','daily','LEAK_Ca','Calcium leached');
r{end+1} = row('Mg_leaching','Mg leaching','kg Mg ha^{-1}','daily','LEAK_Mg','Magnesium leached');
r{end+1} = row('Si_leaching','Si leaching','kg Si ha^{-1}','daily','LEAK_Si','Silicon leached');
r{end+1} = row('K_leaching','K leaching','kg K ha^{-1}','daily','LEAK_K','Potassium leached');
r{end+1} = row('P_leaching','P leaching','kg P ha^{-1}','daily','LEAK_P','Phosphorus leached');
%% Carbon [g C m-2 d-1], from the crown-weighted vegetation and the soil
r{end+1} = row('NEE','NEE','g C m^{-2}','derived','NEE', ...
    'Net ecosystem exchange, positive to the atmosphere','scale',1,'diverging',true);
r{end+1} = row('NPP','Net primary production','g C m^{-2}','derived','NPP', ...
    'Crown-weighted NPP of both vegetation layers','scale',1);
r{end+1} = row('Rh','Heterotrophic respiration','g C m^{-2}','derived','Rh', ...
    'Litter, microbial and earthworm respiration','scale',1);
r{end+1} = row('crop_harvest','Crop harvest','g C m^{-2}','derived','crop_harvest', ...
    'Fruit carbon removed at harvest (crop cells)','scale',1);
r{end+1} = row('LAI','Leaf area index','m^2 m^{-2}','derived','LAI', ...
    'Crown-weighted LAI of both vegetation layers','scale',1,'stat','mean');
%% Soil chemistry
r{end+1} = row('pH_bulk','Bulk soil pH','pH','daily','PH_bulk_d', ...
    'Bulk (extract) pH','scale',1,'stat','hmean');
r{end+1} = row('pH_pore','Pore-water pH','pH','daily','PH_pore_d', ...
    'Pore-water pH','scale',1,'stat','hmean');
r{end+1} = row('NO3_pool','Soil NO_3^-',kgN,'pool','P', ...
    'Nitrate in the biogeochemical layer','index',32,'stat','end');
r{end+1} = row('NH4_pool','Soil NH_4^+',kgN,'pool','P', ...
    'Ammonium in the biogeochemical layer','index',31,'stat','end');
r{end+1} = row('SOC','Soil organic C','t C ha^{-1}','pool','P', ...
    'Particulate (lignin, cellulose) and mineral-associated SOM carbon', ...
    'index',[9 10 11],'stat','end','scale',0.01);
r{end+1} = row('POC','Particulate organic C','t C ha^{-1}','pool','P', ...
    'SOM POC, lignin and cellulose','index',[9 10],'stat','end','scale',0.01);
r{end+1} = row('MAOC','Mineral-associated organic C','t C ha^{-1}','pool','P', ...
    'SOM MOC','index',11,'stat','end','scale',0.01);
r{end+1} = row('SON','Soil organic N','t N ha^{-1}','pool','P', ...
    'SOM nitrogen','index',26,'stat','end','scale',0.01);
r{end+1} = row('microbial_C','Microbial biomass C','g C m^{-2}','pool','P', ...
    'Bacteria and fungi','index',[18 19],'stat','end','scale',1);
r{end+1} = row('Ca_exchangeable','Exchangeable Ca','g Ca m^{-2}','pool','P', ...
    'Calcium on the exchanger','index',61,'stat','end','scale',1);
r{end+1} = row('Ca_solution','Soil-solution Ca','g Ca m^{-2}','pool','P', ...
    'Calcium in the mineral solution','index',60,'stat','end','scale',1);
%% Water
r{end+1} = row('WFPS','Water-filled pore space','%','derived','WFPS', ...
    'Depth-weighted WFPS of the biogeochemical layer','scale',1,'stat','mean');
r{end+1} = row('soil_moisture','Soil moisture','m^3 m^{-3}','derived','theta_bio', ...
    'Depth-weighted water content of the biogeochemical layer','scale',1,'stat','mean');
r{end+1} = row('water_table','Water-table depth','m','hourly','ZWT', ...
    'Depth to the water table','hourly','mean','scale',0.001,'stat','mean');
r{end+1} = row('runoff','Surface runoff','mm','derived','runoff', ...
    'Infiltration- and saturation-excess runoff generated in the cell','scale',1);
r{end+1} = row('runon','Overland run-on','mm','hourly','q_runon', ...
    'Routed runoff received from upslope cells (DEM runs)','hourly','sum','rate',true,'scale',1);
r{end+1} = row('infiltration','Infiltration','mm','hourly','f', ...
    'Water entering the soil surface','hourly','sum','rate',true,'scale',1);
r{end+1} = row('drainage','Drainage at the column base','mm','hourly','Lk', ...
    'Leakage through the bottom of the soil column','hourly','sum','rate',true,'scale',1);
%% Subsurface lateral exchange between cells (field grids and DEM runs)
r{end+1} = row('lateral_water','Net lateral inflow','mm','lateral','water_cell_net_mm', ...
    'Subsurface (Darcy) water gained from neighbouring cells','scale',1,'diverging',true);
r{end+1} = row('lateral_NO3','Net lateral NO_3^- inflow',kgN,'lateral','NO3_cell_net_g_m2', ...
    'Nitrate gained with the subsurface lateral flow','diverging',true);
rows = [r{:}];
end

%% =======================================================================
function s = row(key,label,units,source,name,description,varargin)
s = struct('key',key,'label',label,'units',units,'source',source, ...
    'name',name,'index',[],'hourly','sum','rate',false,'stat','sum', ...
    'scale',10,'diverging',false,'description',description);
for k = 1:2:numel(varargin)
    s.(varargin{k}) = varargin{k+1};
end
end

%% =======================================================================
function s = complete_row(s)
%% A caller's own row: fill what it leaves out with a daily flux.
if ~isscalar(s) || ~isfield(s,'name') || isempty(s.name)
    error('Field_Cell_Variables:BadRow', ...
        'A custom variable needs at least a saved variable name (.name).');
end
default = row(char(s.name),char(s.name),'','daily',char(s.name),'');
default.scale = 1;
f = fieldnames(default);
for k = 1:numel(f)
    if ~isfield(s,f{k}) || isempty(s.(f{k}))
        if ~(strcmp(f{k},'index') && isfield(s,'index'))
            s.(f{k}) = default.(f{k});
        end
    end
end
s = orderfields(rmfield(s,setdiff(fieldnames(s),f)),default);
sources = {'daily','pool','hourly','lateral','derived'};
stats = {'sum','mean','hmean','end','start','change','max','min'};
if ~ismember(s.source,sources)
    error('Field_Cell_Variables:BadRow','Source must be one of: %s.',strjoin(sources,', '));
end
if ~ismember(s.stat,stats)
    error('Field_Cell_Variables:BadRow','Statistic must be one of: %s.',strjoin(stats,', '));
end
if strcmp(s.source,'pool') && isempty(s.index)
    error('Field_Cell_Variables:BadRow','A pool row needs the P columns (.index).');
end
end

%% =======================================================================
function print_rows(rows)
fprintf('%-18s %-14s %-7s %s\n','key','units','stat','description');
for k = 1:numel(rows)
    units = strrep(strrep(rows(k).units,'^{',''),'}','');
    fprintf('%-18s %-14s %-7s %s\n',rows(k).key,units,rows(k).stat,rows(k).description);
end
end
