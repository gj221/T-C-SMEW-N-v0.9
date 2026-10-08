function value = Field_Cell_Series(res,spec)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Field_Cell_Series
%%
%% Daily per-cell series [NNd x N_cells] of one Field_Cell_Variables row,
%% read from a result opened with Field_Cell_Result, in the units the model
%% saves; the row's scale is applied later, to the values drawn.
%%
%%   value = Field_Cell_Series(res,'N2O')
%%   value = Field_Cell_Series(res,Field_Cell_Variables('pH_bulk'))
%%
%% Row 1 is the initial state. Step series are combined into the daily
%% rows the way BIOGEO fills them (Daily_WFPS_From_Hourly): row 2 holds
%% the first step, row j > 2 the 24 steps ending at its 01:00, and the
%% steps after the last full row are left out. Large arrays are read
%% through matfile, a slice of pools or a block of cells at a time, so a
%% field grid of thousands of cells does not have to fit in memory whole.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~res.has_series
    error('Field_Cell_Series:NoSeries', ...
        'This is a grid, not a model result: there are no series to read.');
end
spec = Field_Cell_Variables(spec);
switch spec.source
    case 'daily'
        value = daily_variable(res,spec.name);
    case 'pool'
        value = pool_sum(res,spec.index);
    case 'hourly'
        value = steps_to_days(res,spec.name,1,[],spec.hourly,spec.rate);
    case 'lateral'
        value = lateral_ledger(res,spec.name);
    case 'derived'
        value = derived(res,spec.name);
    otherwise
        error('Field_Cell_Series:BadSource','Unknown source "%s".',spec.source);
end
if ~isequal(size(value),[res.NNd res.N_cells])
    error('Field_Cell_Series:ShapeMismatch', ...
        '%s gave a %s series; expected %d days x %d cells.', ...
        spec.key,mat2str(size(value)),res.NNd,res.N_cells);
end
end

%% =======================================================================
function require(res,names)
names = cellstr(names);
missing = setdiff(names,{res.variables.name});
if ~isempty(missing)
    error('Field_Cell_Series:MissingVariable', ...
        'The result %s does not hold %s.',res.name,strjoin(missing,', '));
end
end

%% =======================================================================
function sz = saved_size(res,name)
k = find(strcmp({res.variables.name},name),1);
sz = res.variables(k).size;
end

%% =======================================================================
function X = daily_variable(res,name)
%% A [NNd x 1 x N] (or [NNd x N]) daily series as [NNd x N].
require(res,name);
L = load(res.file,name);
X = double(L.(name));
if size(X,1) ~= res.NNd || numel(X) ~= res.NNd*res.N_cells
    error('Field_Cell_Series:ShapeMismatch', ...
        '%s is %s; expected %d days x 1 x %d cells.', ...
        name,mat2str(size(X)),res.NNd,res.N_cells);
end
X = reshape(X,res.NNd,res.N_cells);
end

%% =======================================================================
function X = pool_sum(res,index)
%% Sum of P columns, one column read at a time.
require(res,'P');
m = matfile(res.file);
one_cell = numel(saved_size(res,'P')) == 2;
X = zeros(res.NNd,res.N_cells);
for k = index(:).'
    if one_cell
        X = X + double(m.P(:,k));
    else
        X = X + reshape(double(m.P(:,k,:)),res.NNd,res.N_cells);
    end
end
end

%% =======================================================================
function X = lateral_ledger(res,field)
names = {res.variables.name};
if ismember('lateral_diag',names)
    L = load(res.file,'lateral_diag');
    lat = L.lateral_diag;
elseif ismember('Domain',names)
    L = load(res.file,'Domain');
    lat = L.Domain.lateral;
else
    lat = struct();
end
if ~isfield(lat,field)
    error('Field_Cell_Series:NoLateralLedger', ...
        ['%s holds no per-cell lateral ledger (lateral_diag.%s): it is ' ...
        'kept on field grids and with topography, from package x on.'],res.name,field);
end
X = double(lat.(field));
end

%% =======================================================================
function X = derived(res,name)
N = res.N_cells;
switch name
    case {'NEE','NPP','Rh'}
        if ~strcmp(name,'Rh')
            npp = crown_sum(res,'NPP_H','NPP_L');
        end
        if ~strcmp(name,'NPP')
            rh = daily_variable(res,'R_litter') + ...
                daily_variable(res,'R_microbe') + daily_variable(res,'R_ew');
        end
        switch name
            case 'NEE'; X = rh - npp;
            case 'NPP'; X = npp;
            case 'Rh';  X = rh;
        end
    case 'LAI'
        X = crown_sum(res,'LAI_H','LAI_L');
    case 'crop_harvest'
        %% Both Kantola crops set fract_left_fr = 0, so all fruit removal
        %% (RB_L class 5) is yield (Plot_1D_Then_All_12's carbon_dynamics).
        require(res,{'RB_L','Ccrown_t'});
        m = matfile(res.file);
        sz = saved_size(res,'RB_L');
        crown = load(res.file,'Ccrown_t');
        crown = reshape(double(crown.Ccrown_t),res.NNd,[],N);
        if numel(sz) == 3
            fruit = double(m.RB_L(:,:,5));
        else
            fruit = double(m.RB_L(:,:,5,:));
        end
        fruit = reshape(fruit,res.NNd,sz(2),N);
        X = reshape(sum(crown.*max(0,fruit),2),res.NNd,N);
        if isempty(res.crop_cell) || ~any(res.crop_cell)
            if ~isempty(res.ColCfg); X(:) = 0; end
        else
            X(:,~res.crop_cell) = 0;
        end
    case {'WFPS','theta_bio'}
        X = biogeochemical_layer_water(res,strcmp(name,'WFPS'));
    case 'runoff'
        X = steps_to_days(res,'Rh',1,[],'sum',false) + ...
            steps_to_days(res,'Rd',1,[],'sum',false);
    otherwise
        error('Field_Cell_Series:UnknownDerived','No derived quantity "%s".',name);
end
end

%% =======================================================================
function X = crown_sum(res,high,low)
%% Sum over crown classes of Ccrown x (high + low vegetation) [NNd x N].
require(res,{high,low,'Ccrown_t'});
L = load(res.file,high,low,'Ccrown_t');
crown = reshape(double(L.Ccrown_t),res.NNd,[],res.N_cells);
both = reshape(double(L.(high))+double(L.(low)),res.NNd,[],res.N_cells);
X = reshape(sum(crown.*both,2),res.NNd,res.N_cells);
end

%% =======================================================================
function X = biogeochemical_layer_water(res,as_wfps)
%% Depth-weighted water content [-] or WFPS [%] of the biogeochemical
%% layer from the hourly moisture O [NN x ms x N], as
%% Daily_WFPS_From_Hourly computes it for one column.
require(res,{'O','PAR'});
L = load(res.file,'PAR');
PAR = L.PAR;
if ~iscell(PAR) || numel(PAR) ~= res.N_cells
    error('Field_Cell_Series:BadParameters', ...
        'PAR must hold one parameter set per cell.');
end
weights = cell(res.N_cells,1);
n_layer = 0;
for c = 1:res.N_cells
    bio = double(PAR{c}.Bio_Zs(:).');
    if as_wfps
        osat = double(PAR{c}.Soil_Param.Osat(:).');
        weights{c} = 100*bio/sum(bio.*osat);
    else
        weights{c} = bio/sum(bio);
    end
    n_layer = max(n_layer,find(bio > 0,1,'last'));
end
weights = cellfun(@(w) w(1:n_layer),weights,'UniformOutput',false);
X = steps_to_days(res,'O',1:n_layer,weights,'mean',false);
end

%% =======================================================================
function X = steps_to_days(res,name,layers,weights,how,is_rate)
%% Read a step series a block of cells at a time, reduce its layers
%% (weights{c}, or the single layer) and combine the steps of each daily
%% row by sum or mean. Saved shapes: [NN x layers x N], or for one cell
%% [NN x layers], or a single-layer [NN x N].
require(res,name);
sz = saved_size(res,name);
NN = sz(1);
N = res.N_cells;
if numel(sz) >= 3
    layout = 'steps_layers_cells';
    ok = prod(sz(3:end)) == N && sz(2) >= layers(end);
elseif N == 1
    layout = 'steps_layers';
    ok = sz(2) >= layers(end);
else
    layout = 'steps_cells';
    ok = sz(2) == N && isscalar(layers) && layers == 1;
end
if ~ok
    error('Field_Cell_Series:ShapeMismatch', ...
        '%s is %s; expected steps x layers x %d cells.',name,mat2str(sz),N);
end
[A,count] = day_map(res,NN);
m = matfile(res.file);
budget = 2^25;                                   %% values per read
block = max(1,floor(budget/(NN*numel(layers))));
X = zeros(res.NNd,N);
for c1 = 1:block:N
    c2 = min(N,c1+block-1);
    switch layout
        case 'steps_layers_cells'
            raw = double(m.(name)(:,layers(1):layers(end),c1:c2));
        case 'steps_layers'
            raw = double(m.(name)(:,layers(1):layers(end)));
        case 'steps_cells'
            raw = reshape(double(m.(name)(:,c1:c2)),NN,1,c2-c1+1);
    end
    raw = reshape(raw,NN,numel(layers),c2-c1+1);
    steps = zeros(NN,c2-c1+1);
    for c = c1:c2
        if isempty(weights)
            steps(:,c-c1+1) = raw(:,1,c-c1+1);
        else
            steps(:,c-c1+1) = raw(:,:,c-c1+1)*weights{c}(:);
        end
    end
    if is_rate
        steps = steps*res.dth;
    end
    X(:,c1:c2) = A*steps;
end
if strcmp(how,'mean')
    X = X./count;
end
X(count == 0,:) = NaN;
end

%% =======================================================================
function [A,count] = day_map(res,NN)
%% Sparse [NNd x NN] map of steps to the daily rows they fill.
spd = res.steps_per_day;
row = 2 + ceil(((1:NN).'-1)/spd);
keep = row <= res.NNd;
A = sparse(row(keep),find(keep),1,res.NNd,NN);
count = full(sum(A,2));
end
