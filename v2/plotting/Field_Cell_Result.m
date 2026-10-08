function res = Field_Cell_Result(source)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Field_Cell_Result
%%
%% Open a saved 2-D result, a saved field grid or a topology struct for
%% the field figures. Nothing large is read here: Field_Cell_Series reads
%% the series later, slice by slice.
%%
%%   res = Field_Cell_Result(result_file)
%%   res = Field_Cell_Result(topology_file)     a FieldTopologyFile (topo)
%%   res = Field_Cell_Result(topo)
%%
%% res carries
%%   file, name        the file and its stem ('' and 'topology' for a struct)
%%   topo, N_cells     the cell topology the run used
%%   has_series        true for a model result
%%   NN, NNd, dth      steps, daily rows, step length [h]
%%   Date_d, day       the date of each daily row, as Plot_1D_Then_All_12
%%                     dates it: first day + (row - 2) days. Row 1 is the
%%                     initial state and in no window; row 2 is the
%%                     BIOGEO step that closes at the first hour, dated
%%                     the day before the start, and counts as a day of
%%                     flux as it does in the Domain totals
%%   region, is_main   each cell's zone ('main' or the margin)
%%   pft, rock         each cell's vegetation and feedstock labels
%%   area_m2, weight   in-field area and aggregation weight per cell,
%%                     from Field_Aggregation_Areas (polygon overlaps on a
%%                     field grid, whole cells on the plot grid)
%%   has_topography    a DEM is attached to the topology
%%   run_label         R.name when saved, else the file stem
%%   short_name        the stem cut to 60 characters with a checksum, for
%%                     figure folders (result names run to 228 characters)
%%   variables         whos of the file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

here = fileparts(mfilename('fullpath'));
if exist('Field_Aggregation_Areas','file') ~= 2
    addpath(fullfile(fileparts(here),'T&C_Code'));
end

res = struct('file','','name','topology','topo',[],'N_cells',0, ...
    'has_series',false,'NN',0,'NNd',0,'dth',1,'steps_per_day',24, ...
    'Date_d',datetime.empty(0,1),'day',datetime.empty(0,1), ...
    'region',{{}},'is_main',[],'pft',{{}},'rock',{{}},'crop_cell',[], ...
    'area_m2',[],'weight',[],'has_topography',false,'run_label','', ...
    'short_name','topology', ...
    'variables',struct('name',{},'size',{}),'ColCfg',[],'R',[]);

if isstruct(source)
    res.topo = source;
else
    file = char(source);
    if ~isfile(file)
        error('Field_Cell_Result:MissingFile','File not found: %s',file);
    end
    res.file = file;
    [~,res.name] = fileparts(file);
    info = whos('-file',file);
    res.variables = struct('name',{info.name},'size',{info.size});
    names = {info.name};
    if ~ismember('topo',names)
        error('Field_Cell_Result:NoTopology', ...
            'No cell topology (topo) in %s.',file);
    end
    wanted = intersect({'topo','ColCfg','R','dth','NN','NNd'},names);
    L = load(file,wanted{:});
    res.topo = L.topo;
    if isfield(L,'ColCfg'); res.ColCfg = L.ColCfg; end
    if isfield(L,'R'); res.R = L.R; end
    if isfield(L,'dth') && isscalar(L.dth) && L.dth > 0; res.dth = double(L.dth); end
    res.has_series = all(ismember({'Datam','P'},names));
    if res.has_series
        res.NNd = saved_count(L,'NNd',info,'P',1);
        res.NN = saved_count(L,'NN',info,'Datam',1);
        res.steps_per_day = round(24/res.dth);
        m = matfile(file);
        first = double(m.Datam(1,1:3));
        start = datetime(first(1),first(2),first(3));
        res.Date_d = start + caldays((0:res.NNd-1).'-2);
        res.day = res.Date_d;
    end
end

topo = res.topo;
if ~isstruct(topo) || ~isfield(topo,'N_cells') || ~isfield(topo,'cell_x_lo')
    error('Field_Cell_Result:BadTopology', ...
        'The topology must be a Build_*_Topology struct with cell rectangles.');
end
N = topo.N_cells;
res.N_cells = N;
res.region = cellstr(topo.cell_region(:));
res.is_main = strcmp(res.region,'main');
[res.area_m2,res.weight] = Field_Aggregation_Areas(topo);
res.has_topography = isfield(topo,'topography_attached') && ...
    topo.topography_attached && isfield(topo,'cell_z_m');

%% Vegetation and feedstock, as the run configured each column.
res.pft = repmat({''},N,1);
res.rock = repmat({''},N,1);
res.crop_cell = false(N,1);
if ~isempty(res.ColCfg)
    if numel(res.ColCfg) ~= N
        error('Field_Cell_Result:BadConfiguration', ...
            'ColCfg has %d columns for %d cells.',numel(res.ColCfg),N);
    end
    treatment = '';
    if isstruct(res.R) && isfield(res.R,'treatment')
        treatment = char(res.R.treatment);
    end
    for c = 1:N
        cfg = res.ColCfg(c);
        [res.pft{c},res.crop_cell(c)] = pft_label(cfg);
        res.rock{c} = rock_label(cfg,treatment);
    end
end
if isstruct(res.R) && isfield(res.R,'name') && ~isempty(res.R.name)
    res.run_label = char(res.R.name);
else
    res.run_label = res.name;
end
res.short_name = short_name(res.name);
end

%% =======================================================================
function name = short_name(name)
if numel(name) > 60
    h = mod(sum(double(name).*(1:numel(name)).^2),2^32);
    name = sprintf('%s_%08x',name(1:51),h);
end
end

%% =======================================================================
function n = saved_count(L,scalar_name,info,array_name,dim)
if isfield(L,scalar_name) && isscalar(L.(scalar_name))
    n = double(L.(scalar_name));
    return
end
k = find(strcmp({info.name},array_name),1);
n = info(k).size(dim);
end

%% =======================================================================
function [label,is_crop] = pft_label(cfg)
override = '';
if isfield(cfg,'PFT_OVERRIDE') && ~isempty(cfg.PFT_OVERRIDE)
    override = char(cfg.PFT_OVERRIDE);
end
is_crop = isempty(override);
if is_crop
    label = 'Crop';
    if isfield(cfg,'PARAM_IC') && ~isempty(cfg.PARAM_IC)
        [~,param] = fileparts(char(cfg.PARAM_IC));
        param = regexprep(param,'^MOD_PARAM_','');
        label = ['Crop (' strrep(param,'_',' ') ')'];
    end
else
    [~,override] = fileparts(override);
    label = strrep(regexprep(override,'^PFT_',''),'_',' ');
end
end

%% =======================================================================
function label = rock_label(cfg,treatment)
label = 'None';
if isfield(cfg,'apply_rock') && cfg.apply_rock
    if isempty(treatment) || strcmpi(treatment,'control')
        label = 'Feedstock';
    else
        label = [upper(treatment(1)) treatment(2:end)];
    end
end
end
