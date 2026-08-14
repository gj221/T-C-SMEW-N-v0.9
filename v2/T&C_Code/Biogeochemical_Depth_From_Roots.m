function [Zbio,Bio_Zs] = Biogeochemical_Depth_From_Roots( ...
    Zs,ZR95_H,ZR95_L)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set the active biogeochemical depth from the deepest active PFT root.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~(isnumeric(Zs) && isvector(Zs) && numel(Zs) >= 2 && ...
        all(isfinite(Zs)) && Zs(1) == 0 && all(diff(Zs) > 0)) || ...
        any(~isfinite([ZR95_H(:);ZR95_L(:)])) || ...
        any([ZR95_H(:);ZR95_L(:)] < 0)
    error('Biogeochemical_Depth_From_Roots:BadInput', ...
        'The layer mesh and root depths must be finite and nonnegative.');
end
root_depth = [ZR95_H(:);ZR95_L(:)];
root_depth = root_depth(root_depth > 0);
if isempty(root_depth)
    error('Biogeochemical_Depth_From_Roots:NoRoots', ...
        'At least one active PFT root depth is required.');
end
Zbio = max(root_depth);
if Zbio > Zs(end)
    error('Biogeochemical_Depth_From_Roots:RootBelowMesh', ...
        'Root depth %.0f mm exceeds the %.0f mm soil mesh.',Zbio,Zs(end));
end
Bio_Zs = Evaporation_layers(Zs,Zbio);
if abs(sum(Bio_Zs)-1) > 1e-12
    error('Biogeochemical_Depth_From_Roots:BadWeights', ...
        'Biogeochemical layer weights must sum to one.');
end

end
