function [Q_active,edge_fraction] = Restrict_Lateral_Flux_To_Active_Depth( ...
    topo,Q_edge,active_fraction)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Restrict a layer flux to the active-depth overlap at each cell face.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ms = size(Q_edge,2);
if ~isequal(size(Q_edge),[topo.N_edges ms]) || ...
        ~isequal(size(active_fraction),[topo.N_cells ms])
    error('Restrict_Lateral_Flux_To_Active_Depth:BadSize', ...
        'Flux and active-depth arrays must match the topology.');
end
if any(~isfinite(Q_edge(:))) || any(~isfinite(active_fraction(:))) || ...
        any(active_fraction(:) < 0) || any(active_fraction(:) > 1)
    error('Restrict_Lateral_Flux_To_Active_Depth:BadState', ...
        'Fluxes must be finite and active fractions must lie within [0,1].');
end

edge_fraction = min(active_fraction(topo.edge_a,:), ...
    active_fraction(topo.edge_b,:));
if isfield(topo,'edge_active')
    edge_fraction = edge_fraction.*topo.edge_active(:);
end
Q_active = Q_edge.*edge_fraction;

end
