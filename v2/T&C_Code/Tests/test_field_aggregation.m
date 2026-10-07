function test_field_aggregation()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% In-field totals use polygon overlap, while full cells remain simulated.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_directory = fileparts(mfilename('fullpath'));
code_directory = fileparts(test_directory);
addpath(code_directory);

topo.cell_area_m2 = [100;100;100];
topo.cell_overlap_area_m2 = [75;25;0];
cell_values = [1 2 100; 3 4 200];

[field_mean,field_total] = Field_Aggregate_Areal(cell_values,topo);
[accounting_area,weights,fraction] = Field_Aggregation_Areas(topo);

assert(max(abs(accounting_area-[75;25;0])) < 1e-12);
assert(max(abs(weights-[0.75 0.25 0])) < 1e-12);
assert(max(abs(fraction-[0.75;0.25;0])) < 1e-12);
assert(max(abs(field_mean-[1.25;3.25])) < 1e-12);
assert(max(abs(field_total-[125;325])) < 1e-12);

field_pH = pH_Area_Weighted([6 7 9],weights);
expected_pH = -log10(0.75*10^-6+0.25*10^-7);
assert(abs(field_pH-expected_pH) < 1e-12);

%% A majority-rule grid leaves field land in its dropped squares; the
%% totals carry it at the field mean, the mean and weights are unchanged.
majority = topo;
majority.field_missed_area_m2 = 20;
[majority_mean,majority_total] = Field_Aggregate_Areal(cell_values,majority);
[~,majority_weights,~,majority_area] = Field_Aggregation_Areas(majority);
assert(abs(majority_area-120) < 1e-12);
assert(max(abs(majority_weights-weights)) < 1e-12);
assert(max(abs(majority_mean-field_mean)) < 1e-12);
assert(max(abs(majority_total-field_mean*120)) < 1e-12);
[~,~,~,overlap_area] = Field_Aggregation_Areas(topo);
assert(overlap_area == 100);

legacy.cell_area_m2 = [3;1];
[legacy_mean,legacy_total] = Field_Aggregate_Areal([2 6],legacy);
assert(abs(legacy_mean-3) < 1e-12);
assert(abs(legacy_total-12) < 1e-12);

bad = topo;
bad.cell_overlap_area_m2 = [0;0;0];
expect_error(@() Field_Aggregation_Areas(bad), ...
    'Field_Aggregation_Areas:EmptyField');
bad.cell_overlap_area_m2 = [101;0;0];
expect_error(@() Field_Aggregation_Areas(bad), ...
    'Field_Aggregation_Areas:BadOverlapArea');
expect_error(@() Field_Aggregate_Areal(ones(2,2),topo), ...
    'Field_Aggregate_Areal:SizeMismatch');

%% This tree's driver takes the polygon-grid weights from the same function
%% (its domain means are Field_Aggregate_Areal's), and keeps the plot's own
%% area weights, bit for bit, on a grid without overlap areas.
driver_source = fileread(fullfile(code_directory,'MAIN_FRAME_BG_ERW_2D.m'));
assert(contains(driver_source, ...
    '[A_field,w_col,~,A_field_total] = Field_Aggregation_Areas(topo);'));
assert(contains(driver_source,'w_col = reshape(A_cell/A_tot, 1, N_col);'));
assert(contains(driver_source,'Domain.cell_overlap_area_m2'));
assert(contains(driver_source,'Domain.exit_water_field_total_m3'));
[plot_mean,plot_total] = Field_Aggregate_Areal(cell_values,struct('cell_area_m2',[100;100;100]));
assert(isequal(plot_mean,cell_values*([100;100;100]/300)) && ...
    max(abs(plot_total-plot_mean*300)) < 1e-9, ...
    'Without overlap areas the weights are the plot''s area weights.');

fprintf('test_field_aggregation: all checks passed\n');

end

function expect_error(fcn,identifier)
try
    fcn();
catch ME
    assert(strcmp(ME.identifier,identifier), ...
        'Expected %s, got %s.',identifier,ME.identifier);
    return
end
error('Expected error %s was not raised.',identifier);
end
