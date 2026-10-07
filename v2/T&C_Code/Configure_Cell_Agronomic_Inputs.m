function CellInputs = Configure_Cell_Agronomic_Inputs(MOD_base,ERW_base,resolved,CellBGStates)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Apply resolved agronomic properties to each cell's model input rows.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~istable(MOD_base) || ~istable(ERW_base) || ~istable(resolved)
    error('Configure_Cell_Agronomic_Inputs:InvalidInput', ...
        'MOD, ERW and resolved agronomic inputs must be tables.');
end
n = height(resolved);
if ~any(height(MOD_base) == [1 n]) || ~any(height(ERW_base) == [1 n])
    error('Configure_Cell_Agronomic_Inputs:RowCount', ...
        'MOD and ERW inputs must contain one row or one row per cell.');
end
required_MOD = {'PH','Psan','Pcla','Porg'};
required_ERW = {'CEC_tot_in'};
if ~all(ismember(required_MOD,MOD_base.Properties.VariableNames)) || ...
        ~all(ismember(required_ERW,ERW_base.Properties.VariableNames))
    error('Configure_Cell_Agronomic_Inputs:MissingColumns', ...
        'MOD requires PH, Psan, Pcla and Porg; ERW requires CEC_tot_in.');
end
if nargin < 4
    CellBGStates = [];
end
if ~isempty(CellBGStates) && (~isstruct(CellBGStates) || numel(CellBGStates) ~= n)
    error('Configure_Cell_Agronomic_Inputs:StateCount', ...
        'CellBGStates must contain one state per resolved cell.');
end

CellInputs = repmat(struct('cell_id',0,'MOD_PARAM',table(), ...
    'ERW_PARAM',table(),'BG_STATE',struct(),'agronomic_input',table(), ...
    'agronomic_identity',''),n,1);
for c = 1:n
    im = min(c,height(MOD_base));
    ie = min(c,height(ERW_base));
    M = MOD_base(im,:);
    E = ERW_base(ie,:);
    M.PH = resolved.bulk_pH(c);
    M.Psan = resolved.sand_fraction(c);
    M.Pcla = resolved.clay_fraction(c);
    M.Porg = resolved.OM_fraction(c);
    E.CEC_tot_in = resolved.CEC_cmolc_kg(c);
    %% Base saturation (measured, or the site's baseline) becomes this
    %% cell's ERW exchanger seed, f_CEC_in = [Ca Mg K Na Al H].
    %% Invert_Bulk_pH treats a supplied seed as a conserved exchangeable
    %% inventory that it combines with the dissolved cations and
    %% re-partitions by the daily solver's own Gaines-Thomas law at the
    %% pore-water pH that reproduces the measured bulk pH, so the given
    %% stock reaches the handoff exchanger as a fixed point of the model's
    %% exchange law.  The acidity share
    %% (1 - Ca - Mg - K) starts as exchangeable H+; Na and Al exchange
    %% start at zero and follow that equilibrium.  Absent -> the row's own
    %% f_CEC_in stands and the handoff derives the exchanger from pore water
    %% and bulk pH alone, exactly as before.
    if ismember('base_saturation_source',resolved.Properties.VariableNames) && ...
            any(resolved.base_saturation_source(c) == ["measured" "baseline"])
        cec = resolved.CEC_cmolc_kg(c);
        f = [resolved.exch_Ca_cmolc_kg(c) resolved.exch_Mg_cmolc_kg(c) ...
            resolved.exch_K_cmolc_kg(c)]/cec;
        if ~(all(isfinite(f)) && all(f > 0) && sum(f) < 1)
            error('Configure_Cell_Agronomic_Inputs:InvalidBaseSaturation', ...
                'Cell %d: %s base saturation fractions %s cannot seed the exchanger.', ...
                resolved.cell_id(c),resolved.base_saturation_source(c),mat2str(f,6));
        end
        E.f_CEC_in = string(sprintf('[%.12g %.12g %.12g 0 0 %.12g]',f(1),f(2),f(3),1-sum(f)));
        %% The workbook's Kantola convention prescribes its own 80 % base
        %% saturation in ERW_pre_processing and would overwrite this seed;
        %% v0.9's f_CEC_inventory convention keeps supplied fractions as the
        %% conserved initial inventory (same full-charge anion closure).
        if ismember('fX_conc_convention',E.Properties.VariableNames)
            E.fX_conc_convention = "f_CEC_inventory";
        end
    end
    %% Land cover, where given (Cell_Land_Cover): crown 1 is the tree
    %% canopy with the low vegetation beneath it, crown 2 the open low
    %% vegetation; the rest of the ground is bare, urban, water or rock.
    %% The base parameter file reads Ccrown and Cbare from the row, and
    %% Curb, Cwat and Crock when the row carries them.
    if ismember('cover_tree',resolved.Properties.VariableNames) && ...
            isfinite(resolved.cover_tree(c))
        M.Ccrown = string(sprintf('[%.12g %.12g]', ...
            resolved.cover_tree(c),resolved.cover_grass(c)));
        M.Cbare = resolved.cover_bare(c);
        M.Curb = resolved.cover_urban(c);
        M.Cwat = resolved.cover_water(c);
        M.Crock = resolved.cover_rock(c);
    end
    CellInputs(c).cell_id = resolved.cell_id(c);
    CellInputs(c).MOD_PARAM = M;
    CellInputs(c).ERW_PARAM = E;
    CellInputs(c).agronomic_input = resolved(c,:);
    CellInputs(c).agronomic_identity = char(resolved.input_identity(c));
    if ~isempty(CellBGStates)
        if CellBGStates(c).cell_id ~= resolved.cell_id(c) || ...
                ~strcmp(char(string(CellBGStates(c).cell_agronomic_identity)), ...
                CellInputs(c).agronomic_identity)
            error('Configure_Cell_Agronomic_Inputs:StateIdentity', ...
                'Spin-up state %d does not match resolved cell %d.',c,resolved.cell_id(c));
        end
        CellInputs(c).BG_STATE = CellBGStates(c);
    end
end
end
