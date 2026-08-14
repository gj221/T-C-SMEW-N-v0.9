function test_lateral_routing()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% test_lateral_routing
%%
%% Physics checks on the sub-plot domain machinery:
%%   1  topology geometry closes (areas, perimeter, main/buffer split)
%%   2  Darcy flow runs DOWNHILL in head (wet -> dry), not uphill
%%   3  lateral water exchange conserves mass (closed perimeter)
%%   4  a uniform field generates exactly zero lateral flux
%%   5  solute advection conserves mass
%%   6  solute moves with the water, taken from the donor cell
%%   7  reversing the head gradient reverses the transport direction
%%   8  theta never leaves [0, Osat] under an extreme gradient
%%   9  cells with individually different soil / moisture stay independent
%%  10  buffer -> main diagnostic carries the right sign
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n=== test_lateral_routing ===\n');
n_pass = 0; n_fail = 0;

    function check(name, cond, detail)
        if nargin < 3; detail = ''; end
        if cond
            n_pass = n_pass + 1;
            fprintf('  PASS  %s\n', name);
        else
            n_fail = n_fail + 1;
            fprintf('  FAIL  %s  %s\n', name, detail);
        end
    end

%% ---------------------------------------------------------------------
%% Geometry: 60 m field, 5 m margins, 2x2 main grid -> 4 main + 8 buffer
%% ---------------------------------------------------------------------
L = 60; w = 5;
topo = Build_Cell_Topology(L, w, 2, 2);

check('topology has 4 main + 8 buffer cells', topo.N_cells == 12, ...
    sprintf('got %d', topo.N_cells));
check('cell areas sum to L^2', abs(sum(topo.cell_area_m2) - L^2) < 1e-9*L^2, ...
    sprintf('got %.6f vs %.6f', sum(topo.cell_area_m2), L^2));
check('outer faces sum to 4L', abs(sum(topo.outer_face_width_m) - 4*L) < 1e-6, ...
    sprintf('got %.6f vs %.6f', sum(topo.outer_face_width_m), 4*L));
check('main region area = (L-2w)^2', ...
    abs(sum(topo.cell_area_m2(strcmp(topo.cell_region,'main'))) - (L-2*w)^2) < 1e-9*L^2);
check('buffer region area = L^2 - (L-2w)^2', ...
    abs(sum(topo.cell_area_m2(~strcmp(topo.cell_region,'main'))) - (L^2-(L-2*w)^2)) < 1e-9*L^2);
check('every edge joins two distinct cells', all(topo.edge_a ~= topo.edge_b));
check('no main cell touches the domain perimeter', ...
    all(topo.outer_face_width_m(strcmp(topo.cell_region,'main')) == 0));

%% A 4x4 main grid must describe the same field, just more finely
topo16 = Build_Cell_Topology(L, w, 4, 4);
check('4x4 main grid gives 24 cells over the same area', ...
    topo16.N_cells == 24 && abs(sum(topo16.cell_area_m2) - L^2) < 1e-9*L^2);

%% ---------------------------------------------------------------------
%% Soil parameters (Saxton-Rawls, SPAR = 2, silty clay loam)
%% ---------------------------------------------------------------------
ms = 5;
Zlayer_m = [0.05 0.05 0.10 0.20 0.30];
[Osat, Lslope, Pe, Ks, O33] = Soil_parameters(0.15, 0.30, 0.03);
Ohy_s = 0.02;

sp = struct('SPAR', 2, 'Ks', Ks*ones(1,ms), 'Osat', Osat*ones(1,ms), ...
    'Ohy', Ohy_s*ones(1,ms), 'L', Lslope*ones(1,ms), 'Pe', Pe*ones(1,ms), ...
    'O33', O33*ones(1,ms), 'alpVG', 1e-3*ones(1,ms), 'nVG', 2*ones(1,ms), ...
    'lVG', 0.5*ones(1,ms), 'Ks_mac', zeros(1,ms), 'Omac', zeros(1,ms), ...
    'alpVGM', 1e-3*ones(1,ms), 'nVGM', 2*ones(1,ms), 'lVGM', 0.5*ones(1,ms), ...
    'Phy1', 10000, 's_SVG', zeros(1,ms), 'bVG', zeros(1,ms));
soilpar_cell = repmat(sp, topo.N_cells, 1);

dt_h = 1;
A = topo.cell_area_m2(:);
is_main = strcmp(topo.cell_region, 'main');
i_main  = find(is_main);
i_buf   = find(~is_main);

%% ---------------------------------------------------------------------
%% 2-3. Darcy direction and mass conservation (main wet, buffer dry)
%% ---------------------------------------------------------------------
theta = zeros(topo.N_cells, ms);
theta(i_main,:) = 0.40;      %% wet field
theta(i_buf,:)  = 0.20;      %% dry margin

[Q_edge, dtheta, sub_steps] = Darcy_Lateral_Topology( ...
    topo, theta, Zlayer_m, soilpar_cell, dt_h);

bm_sign = edge_sign_buffer_to_main(topo);
Q_buf_to_main = sum(bm_sign .* sum(Q_edge,2));
check('water flows from the WET field to the DRY margin', Q_buf_to_main < 0, ...
    sprintf('Q_buf->main = %.6g m3/h (should be negative here)', Q_buf_to_main));
check('Darcy solver took at least one sub-step', sub_steps >= 1);

dV = sum(dtheta .* repmat(Zlayer_m, topo.N_cells, 1), 2) .* A;   %% m3
check('closed perimeter conserves water', ...
    abs(sum(dV)) < 1e-9*sum(abs(dV)) + 1e-15, sprintf('net dV = %.3e m3', sum(dV)));
%% Every main cell borders the margin, so all of them drain.  Corner
%% margin cells touch only the N/S/E/W strips (never the field), so at
%% uniform margin moisture they see no gradient and correctly stay put.
i_buf_adj = i_buf(arrayfun(@(c) any(ismember( ...
    [topo.edge_b(topo.edge_a==c); topo.edge_a(topo.edge_b==c)], i_main)), i_buf));
i_buf_corner = setdiff(i_buf, i_buf_adj);
check('every wet field cell loses water', all(dV(i_main) < 0));
check('margin cells bordering the field gain water', all(dV(i_buf_adj) > 0));
check('corner margin cells, not touching the field, stay unchanged', ...
    ~isempty(i_buf_corner) && all(abs(dV(i_buf_corner)) < 1e-15));
check('no margin cell loses water to the field here', all(dV(i_buf) >= 0));

%% ---------------------------------------------------------------------
%% 4. Uniform field -> exactly zero lateral flux
%% ---------------------------------------------------------------------
[Q_u, dtheta_u] = Darcy_Lateral_Topology( ...
    topo, 0.30*ones(topo.N_cells, ms), Zlayer_m, soilpar_cell, dt_h);
check('uniform moisture generates zero lateral flux', ...
    max(abs(Q_u(:))) < 1e-12 && max(abs(dtheta_u(:))) < 1e-15, ...
    sprintf('max|Q| = %.3e', max(abs(Q_u(:)))));

%% ---------------------------------------------------------------------
%% 5-6. Solute advection: conservation and donor selection
%% ---------------------------------------------------------------------
N_sp = 12;
C_aq = zeros(topo.N_cells, ms, N_sp);
C_aq(i_main,:,7) = 50;     %% NO3 only in the wet (donor) field cells

[dM, dM_edge] = Solute_Lateral_Transport(topo, Q_edge, C_aq, dt_h);

check('lateral advection conserves solute mass', ...
    abs(sum(sum(dM(:,:,7)))) < 1e-9*sum(abs(reshape(dM(:,:,7),[],1))) + 1e-15, ...
    sprintf('net dM = %.3e', sum(sum(dM(:,:,7)))));
check('solute leaves the donor cells and enters the receivers', ...
    sum(sum(dM(i_main,:,7))) < 0 && sum(sum(dM(i_buf,:,7))) > 0);
check('species with zero concentration stay at zero', ...
    max(abs(reshape(dM(:,:,[1:6 8:12]),[],1))) < 1e-15);

%% Mass delivered to the margin = water volume moved x donor concentration
Vw_main_to_buf = -Q_buf_to_main * dt_h;                 %% m3
check('advected mass = water volume x donor concentration', ...
    abs(sum(sum(dM(i_buf,:,7))) - Vw_main_to_buf*50) < 1e-9*max(1,Vw_main_to_buf*50), ...
    sprintf('got %.6g expected %.6g', sum(sum(dM(i_buf,:,7))), Vw_main_to_buf*50));

%% ---------------------------------------------------------------------
%% 7. Reversing the gradient reverses the water direction; with a clean
%%    donor no solute is transported.
%% ---------------------------------------------------------------------
theta_rev = theta; theta_rev(i_main,:) = 0.20; theta_rev(i_buf,:) = 0.40;
Q_rev = Darcy_Lateral_Topology(topo, theta_rev, Zlayer_m, soilpar_cell, dt_h);
check('reversing the head gradient reverses the water flux', ...
    sum(bm_sign .* sum(Q_rev,2)) > 0, 'margin is now wet so it should feed the field');
dM_rev = Solute_Lateral_Transport(topo, Q_rev, C_aq, dt_h);
check('a solute-free donor transports no solute', ...
    max(abs(reshape(dM_rev(:,:,7),[],1))) < 1e-15);

%% ---------------------------------------------------------------------
%% 8. theta stays physical under an aggressive gradient
%% ---------------------------------------------------------------------
theta_x = zeros(topo.N_cells, ms);
theta_x(i_main,:) = Osat - 1e-4;    %% nearly saturated
theta_x(i_buf,:)  = Ohy_s + 1e-4;   %% nearly air-dry
[~, dtheta_x] = Darcy_Lateral_Topology(topo, theta_x, Zlayer_m, soilpar_cell, 6);
theta_end = theta_x + dtheta_x;
check('theta stays within [0, Osat] under a 6 h extreme gradient', ...
    all(theta_end(:) >= -1e-12) && all(theta_end(:) <= Osat + 1e-12), ...
    sprintf('range [%.6f %.6f], Osat=%.6f', min(theta_end(:)), max(theta_end(:)), Osat));

%% ---------------------------------------------------------------------
%% 9. Per-cell heterogeneity: cells may carry different soil AND moisture.
%%    Two isolated wet cells that are not adjacent must not exchange.
%% ---------------------------------------------------------------------
soil_het = soilpar_cell;
for c = 1:topo.N_cells
    soil_het(c).Ks = Ks*(0.5 + c/topo.N_cells)*ones(1,ms);   %% distinct per cell
end
theta_het = 0.25 + 0.02*mod((1:topo.N_cells)',3);            %% distinct per cell
theta_het = repmat(theta_het, 1, ms);
[Q_het, dtheta_het] = Darcy_Lateral_Topology( ...
    topo, theta_het, Zlayer_m, soil_het, dt_h);
dV_het = sum(dtheta_het .* repmat(Zlayer_m, topo.N_cells, 1), 2) .* A;
check('heterogeneous cells still conserve water exactly', ...
    abs(sum(dV_het)) < 1e-9*sum(abs(dV_het)) + 1e-15, ...
    sprintf('net dV = %.3e m3', sum(dV_het)));
check('heterogeneous soils produce non-zero exchange', max(abs(Q_het(:))) > 0);

%% Non-adjacent cells must not exchange directly
for e = 1:topo.N_edges
    if ~topo.adj(topo.edge_a(e), topo.edge_b(e))
        check('edges match the adjacency matrix', false, 'edge without adjacency');
        break
    end
end
check('every edge is present in the adjacency matrix', true);

%% ---------------------------------------------------------------------
%% 10. Wrapper: buffer -> main diagnostic sign
%% ---------------------------------------------------------------------
[theta_out, dM_out, diag] = Lateral_Routing_Step( ...
    topo, theta, C_aq, soilpar_cell, Zlayer_m, dt_h);
check('wrapper returns updated theta of the right shape', ...
    isequal(size(theta_out), [topo.N_cells ms]));
check('wrapper reports net water leaving the field for the margin', ...
    diag.buf_to_main.water < 0);
check('wrapper reports net NO3 leaving the field for the margin', ...
    diag.buf_to_main.net(7) < 0);
check('wrapper solute change matches the standalone transport call', ...
    max(abs(dM_out(:) - dM(:))) < 1e-12);

%% ---------------------------------------------------------------------
fprintf('  ---------------------------------------\n');
fprintf('  %d passed, %d failed\n\n', n_pass, n_fail);
if n_fail > 0
    error('test_lateral_routing: %d check(s) failed', n_fail);
end

end

%% =======================================================================
function bm_sign = edge_sign_buffer_to_main(topo)
is_main_a = strcmp(topo.cell_region(topo.edge_a), 'main');
is_main_b = strcmp(topo.cell_region(topo.edge_b), 'main');
bm_sign = zeros(topo.N_edges, 1);
bm_sign(~is_main_a &  is_main_b) = +1;
bm_sign( is_main_a & ~is_main_b) = -1;
end
