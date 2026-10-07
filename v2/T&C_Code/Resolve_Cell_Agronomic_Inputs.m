function resolved = Resolve_Cell_Agronomic_Inputs(topo,input_table)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Resolve measured soil properties or their interpolated cell baselines.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isstruct(topo) || ~isscalar(topo) || ~isfield(topo,'N_cells') || ...
        ~isfield(topo,'cell_id')
    error('Resolve_Cell_Agronomic_Inputs:InvalidTopology', ...
        'The topology must contain N_cells and cell_id.');
end
if ~istable(input_table)
    error('Resolve_Cell_Agronomic_Inputs:InvalidInput', ...
        'Cell agronomic inputs must be supplied as a table.');
end
required = {'cell_id','baseline_bulk_pH','measured_bulk_pH', ...
    'use_measured_bulk_pH','baseline_sand_fraction', ...
    'baseline_clay_fraction','measured_sand_fraction', ...
    'measured_clay_fraction','use_measured_texture', ...
    'baseline_OM_fraction','measured_OM_fraction','use_measured_OM', ...
    'baseline_CEC_cmolc_kg','measured_CEC_cmolc_kg','use_measured_CEC'};
missing = required(~ismember(required,input_table.Properties.VariableNames));
if ~isempty(missing)
    error('Resolve_Cell_Agronomic_Inputs:MissingColumns', ...
        'Cell agronomic input is missing: %s.',strjoin(missing,', '));
end

cell_id = input_table.cell_id;
valid_id = isnumeric(cell_id) && isreal(cell_id) && iscolumn(cell_id) && ...
    all(isfinite(cell_id)) && all(cell_id == fix(cell_id)) && all(cell_id >= 1);
if ~valid_id || numel(unique(cell_id)) ~= height(input_table)
    error('Resolve_Cell_Agronomic_Inputs:InvalidCellId', ...
        'cell_id must contain unique positive integers.');
end
topology_id = topo.cell_id(:);
if topo.N_cells ~= numel(topology_id) || height(input_table) ~= topo.N_cells || ...
        ~isequal(sort(cell_id),sort(topology_id))
    error('Resolve_Cell_Agronomic_Inputs:CellSetMismatch', ...
        'The agronomic table must contain one row for every topology cell.');
end
[~,order] = ismember(topology_id,cell_id);
T = input_table(order,:);

[cell_x,cell_y] = topology_centres(topo);
has_x = ismember('cell_x_projected_m',T.Properties.VariableNames);
has_y = ismember('cell_y_projected_m',T.Properties.VariableNames);
if xor(has_x,has_y)
    error('Resolve_Cell_Agronomic_Inputs:MissingColumns', ...
        'Projected cell coordinates must supply both x and y columns.');
end
if has_x
    x = T.cell_x_projected_m;
    y = T.cell_y_projected_m;
    if any(~isfinite(x)) || any(~isfinite(y)) || ...
            any(abs(x-cell_x) > 1e-6) || any(abs(y-cell_y) > 1e-6)
        error('Resolve_Cell_Agronomic_Inputs:CellCoordinateMismatch', ...
            'Projected cell coordinates do not match the selected grid topology.');
    end
end

use_pH = logical_flag(T.use_measured_bulk_pH,'use_measured_bulk_pH');
use_texture = logical_flag(T.use_measured_texture,'use_measured_texture');
use_OM = logical_flag(T.use_measured_OM,'use_measured_OM');
use_CEC = logical_flag(T.use_measured_CEC,'use_measured_CEC');

baseline_names = {'baseline_bulk_pH','baseline_sand_fraction', ...
    'baseline_clay_fraction','baseline_OM_fraction','baseline_CEC_cmolc_kg'};
for k = 1:numel(baseline_names)
    values = T.(baseline_names{k});
    if ~isnumeric(values) || ~isreal(values) || ~iscolumn(values) || ...
            any(~isfinite(values))
        error('Resolve_Cell_Agronomic_Inputs:MissingBaseline', ...
            '%s must contain one finite value for every simulated cell.',baseline_names{k});
    end
end
require_measurement(T.measured_bulk_pH,use_pH,'measured_bulk_pH');
require_measurement(T.measured_sand_fraction,use_texture,'measured_sand_fraction');
require_measurement(T.measured_clay_fraction,use_texture,'measured_clay_fraction');
require_measurement(T.measured_OM_fraction,use_OM,'measured_OM_fraction');
require_measurement(T.measured_CEC_cmolc_kg,use_CEC,'measured_CEC_cmolc_kg');

bulk_pH = choose(T.baseline_bulk_pH,T.measured_bulk_pH,use_pH);
sand = choose(T.baseline_sand_fraction,T.measured_sand_fraction,use_texture);
clay = choose(T.baseline_clay_fraction,T.measured_clay_fraction,use_texture);
OM = choose(T.baseline_OM_fraction,T.measured_OM_fraction,use_OM);
CEC = choose(T.baseline_CEC_cmolc_kg,T.measured_CEC_cmolc_kg,use_CEC);
silt = 1-sand-clay-OM;

if any(bulk_pH <= 0 | bulk_pH > 14)
    error('Resolve_Cell_Agronomic_Inputs:InvalidPH', ...
        'Resolved bulk pH must be greater than zero and no greater than 14.');
end
if any(sand < 0 | sand >= 1 | clay < 0 | clay >= 1 | OM < 0 | OM >= 1) || ...
        any(silt < -1e-12)
    error('Resolve_Cell_Agronomic_Inputs:TextureClosure', ...
        'Sand, clay, organic matter and derived silt fractions must be nonnegative and sum to one.');
end
if any(CEC <= 0)
    error('Resolve_Cell_Agronomic_Inputs:InvalidCEC', ...
        'Resolved CEC must be positive in cmol charge kg-1.');
end
silt = max(silt,0);

%% Optional: exchangeable base saturation (Ca, Mg, K as fractions of CEC),
%% measured or a baseline.  There is no pedotransfer baseline for the
%% exchanger -- the spin-up's own state is the baseline -- so the columns
%% are optional and a cell with neither simply leaves the exchangeable
%% pools to the spin-up.  A site may instead give a baseline (all three
%% baseline_bs_* columns, one finite value per cell, as the other
%% baselines): its representative exchanger until measurements replace
%% it.  Where used, the exchangeable stock handed to the spin-up
%% constraint is fraction x resolved CEC [cmol_c kg-1].
bs_cols = {'measured_bs_Ca_fraction','measured_bs_Mg_fraction', ...
    'measured_bs_K_fraction'};
base_bs_cols = {'baseline_bs_Ca_fraction','baseline_bs_Mg_fraction', ...
    'baseline_bs_K_fraction'};
exch = nan(topo.N_cells,3);
use_bs = false(topo.N_cells,1);
if all(ismember([bs_cols,{'use_measured_base_saturation'}],T.Properties.VariableNames))
    use_bs = logical_flag(T.use_measured_base_saturation,'use_measured_base_saturation');
    F = [T.(bs_cols{1}) T.(bs_cols{2}) T.(bs_cols{3})];
    for k = 1:3
        require_measurement(F(:,k),use_bs,bs_cols{k});
    end
    exch(use_bs,:) = base_saturation(F(use_bs,:),'Measured').*CEC(use_bs);
end
has_base_bs = ismember(base_bs_cols,T.Properties.VariableNames);
if any(has_base_bs) && ~all(has_base_bs)
    error('Resolve_Cell_Agronomic_Inputs:MissingColumns', ...
        'A base-saturation baseline needs all of %s.',strjoin(base_bs_cols,', '));
end
use_base_bs = false(topo.N_cells,1);
if all(has_base_bs)
    B = T{:,base_bs_cols};
    if ~isnumeric(B) || any(~isfinite(B(:)))
        error('Resolve_Cell_Agronomic_Inputs:MissingBaseline', ...
            '%s must contain one finite value for every simulated cell.', ...
            strjoin(base_bs_cols,', '));
    end
    use_base_bs = ~use_bs;
    exch(use_base_bs,:) = base_saturation(B(use_base_bs,:),'Baseline').*CEC(use_base_bs);
end
bs_source = repmat("none",topo.N_cells,1);
bs_source(use_bs) = "measured";
bs_source(use_base_bs) = "baseline";

%% Optional: land cover (Cell_Land_Cover, from a land-cover map or set by
%% hand) -- the cell's ground divided among tree crown (with the low
%% vegetation beneath it), open low vegetation, urban, water, bare soil and
%% rock.  All six columns or none; every row sums to one.  Absent, each
%% cell keeps the parameter row's own cover.
cover_cols = {'cover_tree','cover_grass','cover_urban','cover_water', ...
    'cover_bare','cover_rock'};
has_cover = ismember(cover_cols,T.Properties.VariableNames);
if any(has_cover) && ~all(has_cover)
    error('Resolve_Cell_Agronomic_Inputs:MissingColumns', ...
        'Land cover needs all of %s.',strjoin(cover_cols,', '));
end
cover = nan(topo.N_cells,numel(cover_cols));
if all(has_cover)
    cover = T{:,cover_cols};
    if ~isnumeric(cover) || any(~isfinite(cover(:))) || ...
            any(cover(:) < -1e-9 | cover(:) > 1+1e-9) || any(abs(sum(cover,2)-1) > 1e-6)
        error('Resolve_Cell_Agronomic_Inputs:InvalidLandCover', ...
            'Land-cover fractions must lie in [0,1] and sum to one in every cell.');
    end
end

pH_source = source_label(use_pH);
texture_source = source_label(use_texture);
OM_source = source_label(use_OM);
CEC_source = source_label(use_CEC);
input_identity = strings(topo.N_cells,1);
for c = 1:topo.N_cells
    input_identity(c) = sprintf([ ...
        'cell=%d;x=%.17g;y=%.17g;pH=%.17g:%s;sand=%.17g;clay=%.17g:' ...
        '%s;OM=%.17g:%s;CEC=%.17g:%s'],topology_id(c),cell_x(c),cell_y(c), ...
        bulk_pH(c),pH_source(c),sand(c),clay(c),texture_source(c), ...
        OM(c),OM_source(c),CEC(c),CEC_source(c));
    %% Appended only when used, so identities of cells without base
    %% saturation are unchanged from before this input existed.
    if bs_source(c) ~= "none"
        input_identity(c) = sprintf('%s;exch=%.17g,%.17g,%.17g:%s', ...
            input_identity(c),exch(c,1),exch(c,2),exch(c,3),bs_source(c));
    end
    %% Likewise appended only when land cover is given.
    if all(has_cover)
        input_identity(c) = sprintf('%s;cover=%s',input_identity(c), ...
            strjoin(compose('%.17g',cover(c,:)),','));
    end
end

resolved = table(topology_id,cell_x,cell_y,bulk_pH,sand,clay,silt,OM,CEC, ...
    exch(:,1),exch(:,2),exch(:,3), ...
    pH_source,texture_source,OM_source,CEC_source,bs_source,input_identity, ...
    'VariableNames',{'cell_id','cell_x_projected_m','cell_y_projected_m', ...
    'bulk_pH','sand_fraction','clay_fraction','silt_fraction','OM_fraction', ...
    'CEC_cmolc_kg','exch_Ca_cmolc_kg','exch_Mg_cmolc_kg','exch_K_cmolc_kg', ...
    'bulk_pH_source','texture_source','OM_source','CEC_source', ...
    'base_saturation_source','input_identity'});
resolved = [resolved array2table(cover,'VariableNames',cover_cols)];
end

function F = base_saturation(F,kind)
%% Base-saturation fractions: each in (0,1), their sum below one (the
%% remainder is the acid saturation).
if any(F(:) <= 0 | F(:) >= 1) || any(sum(F,2) >= 1)
    error('Resolve_Cell_Agronomic_Inputs:InvalidBaseSaturation', ...
        ['%s base saturation fractions must each lie in (0,1) and ' ...
         'sum to less than one (the remainder is the acid saturation).'],kind);
end
end

function [x,y] = topology_centres(topo)
if isfield(topo,'cell_xc_projected_m') && isfield(topo,'cell_yc_projected_m')
    x = topo.cell_xc_projected_m(:);
    y = topo.cell_yc_projected_m(:);
else
    x = nan(topo.N_cells,1);
    y = nan(topo.N_cells,1);
end
if numel(x) ~= topo.N_cells || numel(y) ~= topo.N_cells
    error('Resolve_Cell_Agronomic_Inputs:InvalidTopology', ...
        'Projected topology centres must contain one coordinate per cell.');
end
end

function flag = logical_flag(values,name)
valid = (islogical(values) || isnumeric(values)) && isreal(values) && ...
    iscolumn(values) && all(isfinite(double(values))) && ...
    all(double(values) == 0 | double(values) == 1);
if ~valid
    error('Resolve_Cell_Agronomic_Inputs:InvalidToggle', ...
        '%s must contain only logical or numeric zero/one values.',name);
end
flag = logical(values);
end

function require_measurement(values,use_measurement,name)
valid_type = isnumeric(values) && isreal(values) && iscolumn(values);
if ~valid_type || any(~isfinite(values(use_measurement)))
    error('Resolve_Cell_Agronomic_Inputs:MissingMeasurement', ...
        '%s must be finite wherever its measurement toggle is enabled.',name);
end
end

function value = choose(baseline,measured,use_measurement)
value = baseline;
value(use_measurement) = measured(use_measurement);
end

function label = source_label(use_measurement)
label = repmat("baseline",numel(use_measurement),1);
label(use_measurement) = "measured";
end
