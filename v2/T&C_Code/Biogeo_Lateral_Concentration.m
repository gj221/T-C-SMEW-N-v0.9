function C_aq = Biogeo_Lateral_Concentration( ...
    mobile_pool,cell_area_m2,theta_cell,Zlayer_m,active_fraction)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Concentration of T&C mobile pools within the biogeochemical depth.
%% mobile_pool is gX/m2 and C_aq is gX/m3 water.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[N_cells,N_species] = size(mobile_pool);
ms = numel(Zlayer_m);
if ~isequal(size(theta_cell),[N_cells ms]) || ...
        ~isequal(size(active_fraction),[N_cells ms]) || ...
        numel(cell_area_m2) ~= N_cells
    error('Biogeo_Lateral_Concentration:BadSize', ...
        'Pool, water, active-depth and area arrays must have matching cell dimensions.');
end
if any(~isfinite(mobile_pool(:))) || any(mobile_pool(:) < 0) || ...
        any(~isfinite(cell_area_m2(:))) || any(cell_area_m2(:) <= 0) || ...
        any(~isfinite(theta_cell(:))) || any(theta_cell(:) < 0) || ...
        any(~isfinite(Zlayer_m(:))) || any(Zlayer_m(:) <= 0) || ...
        any(~isfinite(active_fraction(:))) || ...
        any(active_fraction(:) < 0) || any(active_fraction(:) > 1)
    error('Biogeo_Lateral_Concentration:BadState', ...
        'Pools and physical state must be finite, nonnegative and within bounds.');
end

cell_area_m2 = cell_area_m2(:);
Zlayer_m = Zlayer_m(:).';
water_volume = sum(theta_cell.*active_fraction.*Zlayer_m,2).*cell_area_m2;
total_mass = mobile_pool.*cell_area_m2;
if any(water_volume <= 0 & any(total_mass > 0,2))
    error('Biogeo_Lateral_Concentration:DryActiveZone', ...
        'A cell with mobile mass must contain water in its active depth.');
end

concentration = zeros(N_cells,N_species);
wet = water_volume > 0;
concentration(wet,:) = total_mass(wet,:)./water_volume(wet);
C_aq = zeros(N_cells,ms,N_species);
for c = 1:N_cells
    C_aq(c,active_fraction(c,:) > 0,:) = repmat( ...
        reshape(concentration(c,:),1,1,N_species), ...
        [1 nnz(active_fraction(c,:) > 0) 1]);
end

end
