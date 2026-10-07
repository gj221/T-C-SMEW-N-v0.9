function MM = ERW_Mineral_Molar_Mass_Row(MM_min,mineral)

MM = reshape(MM_min,1,[]);
if ischar(mineral)
    n_min = size(mineral,1);
else
    n_min = numel(mineral);
end
if numel(MM) ~= n_min
    error('ERW:MineralMolarMassSize', ...
        'Expected one molar mass for each of %d minerals; received %d.', ...
        n_min,numel(MM));
end
end
