function nu_C = ERW_Mineral_DIC_Stoichiometry(mineral)
% Inorganic carbon released per mole of weathered feedstock mineral.

if ischar(mineral)
    names = cellstr(mineral);
elseif isstring(mineral)
    names = cellstr(mineral(:));
elseif iscell(mineral)
    names = mineral(:);
else
    error('ERW:InvalidMineralNames','Mineral names must be text.');
end

names = strtrim(names);
nu_C = zeros(numel(names),1);
nu_C(strcmpi(names,'calcite')) = 1;
end
