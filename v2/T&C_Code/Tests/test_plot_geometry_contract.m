function test_plot_geometry_contract
% Plot and margin areas must be explicit, finite and positive.

code_dir = fileparts(fileparts(mfilename('fullpath')));
addpath(code_dir);

topo = Build_Cell_Topology(60,5,2,2);
is_main = strcmp(topo.cell_region,'main');
assert(abs(sum(topo.cell_area_m2)-3600) < 1e-12, ...
    'The complete modeled square must remain 3600 m2.');
assert(abs(sum(topo.cell_area_m2(is_main))-2500) < 1e-12, ...
    'An internal 5 m ring leaves a 50 m by 50 m main plot.');
assert(abs(sum(topo.cell_area_m2(~is_main))-1100) < 1e-12, ...
    'The internal margin area must close the plot-area identity.');

bad_values = [0,-1,Inf,NaN];
for k = 1:numel(bad_values)
    did_error = false;
    try
        Build_Single_Cell_Topology(bad_values(k));
    catch ME
        did_error = strcmp(ME.identifier,'Build_Single_Cell_Topology:BadArea');
    end
    assert(did_error,'Single-cell area must be finite and positive.');
end

fprintf('PASS: plot geometry and margin areas have explicit boundaries.\n');
end
