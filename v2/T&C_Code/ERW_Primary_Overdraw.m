function [max_overdraw,overdraw,reacted_mass] = ERW_Primary_Overdraw( ...
    EW,dtw,MM_min,mineral,available_mass)

if ~isscalar(dtw) || ~isfinite(dtw) || dtw <= 0
    error('ERW:PrimaryOverdrawTimeStep', ...
        'The weathering time step must be finite and positive.');
end
if ~isequal(size(EW),size(available_mass))
    error('ERW:PrimaryOverdrawSize', ...
        'Weathering flux and primary inventory must have the same size.');
end
if any(~isfinite(EW),'all') || any(EW < 0,'all')
    error('ERW:PrimaryOverdrawFlux', ...
        'Primary weathering flux must be finite and nonnegative.');
end
if any(~isfinite(available_mass),'all') || any(available_mass < 0,'all')
    error('ERW:PrimaryOverdrawInventory', ...
        'Primary mineral inventory must be finite and nonnegative.');
end

if isempty(mineral)
    if any(EW ~= 0,'all') || any(available_mass ~= 0,'all')
        error('ERW:PrimaryOverdrawInactiveInventory', ...
            'A run without primary minerals has nonzero flux or inventory.');
    end
    reacted_mass = zeros(size(EW));
else
    MM = ERW_Mineral_Molar_Mass_Row(MM_min,mineral);
    if size(EW,2) ~= numel(MM)
        error('ERW:PrimaryOverdrawSize', ...
            'Weathering columns must follow the primary-mineral ordering.');
    end
    reacted_mass = EW.*MM*dtw;
end

overdraw = reacted_mass-available_mass;
max_overdraw = max([overdraw(:);0]);
end
