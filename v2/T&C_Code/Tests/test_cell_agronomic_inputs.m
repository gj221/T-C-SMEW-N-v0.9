function test_cell_agronomic_inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cell measurements override interpolated baselines only when enabled.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n=== test_cell_agronomic_inputs ===\n');
topo = struct('N_cells',2,'cell_id',[1;2], ...
    'cell_xc_projected_m',[10;20],'cell_yc_projected_m',[30;40]);

T = table([2;1],[20;10],[40;30], ...
    [6.0;5.5],[7.0;6.5],logical([0;1]), ...
    [0.30;0.25],[0.35;0.40],[0.40;0.45],[0.20;0.15],logical([1;0]), ...
    [0.030;0.025],[0.040;0.035],logical([0;1]), ...
    [8;9],[10;11],logical([1;0]), ...
    'VariableNames',{'cell_id','cell_x_projected_m','cell_y_projected_m', ...
    'baseline_bulk_pH','measured_bulk_pH','use_measured_bulk_pH', ...
    'baseline_sand_fraction','baseline_clay_fraction', ...
    'measured_sand_fraction','measured_clay_fraction','use_measured_texture', ...
    'baseline_OM_fraction','measured_OM_fraction','use_measured_OM', ...
    'baseline_CEC_cmolc_kg','measured_CEC_cmolc_kg','use_measured_CEC'});

R = Resolve_Cell_Agronomic_Inputs(topo,T);
assert(isequal(R.cell_id,[1;2]));
assert(isequal(R.bulk_pH,[6.5;6.0]));
assert(isequal(R.sand_fraction,[0.25;0.40]));
assert(isequal(R.clay_fraction,[0.40;0.20]));
assert(isequal(R.OM_fraction,[0.035;0.030]));
assert(isequal(R.CEC_cmolc_kg,[9;10]));
assert(max(abs(R.silt_fraction-[0.315;0.37])) < 1e-12);
assert(isequal(R.bulk_pH_source,["measured";"baseline"]));
assert(isequal(R.texture_source,["baseline";"measured"]));

MOD = table(6.1,0.1,0.35,0.025,300,1295, ...
    'VariableNames',{'PH','Psan','Pcla','Porg','Zbio','rsd'});
ERW = table(8.6,'VariableNames',{'CEC_tot_in'});
S = repmat(struct('cell_id',0,'cell_agronomic_identity','', ...
    'B_out',zeros(1,79)),2,1);
for c = 1:2
    S(c).cell_id = c;
    S(c).cell_agronomic_identity = char(R.input_identity(c));
    S(c).B_out(c) = c;
end
C = Configure_Cell_Agronomic_Inputs(MOD,ERW,R,S);
assert(C(1).MOD_PARAM.PH == 6.5 && C(2).MOD_PARAM.PH == 6.0);
assert(C(1).MOD_PARAM.Psan == 0.25 && C(2).MOD_PARAM.Pcla == 0.20);
assert(C(1).MOD_PARAM.Porg == 0.035 && C(2).ERW_PARAM.CEC_tot_in == 10);
assert(C(2).BG_STATE.B_out(2) == 2);

T_bad = T;
T_bad.baseline_bulk_pH(T_bad.cell_id == 2) = NaN;
assert_error(@() Resolve_Cell_Agronomic_Inputs(topo,T_bad), ...
    'Resolve_Cell_Agronomic_Inputs:MissingBaseline');
T_bad = T;
T_bad.measured_bulk_pH(T_bad.cell_id == 1) = NaN;
assert_error(@() Resolve_Cell_Agronomic_Inputs(topo,T_bad), ...
    'Resolve_Cell_Agronomic_Inputs:MissingMeasurement');
T_bad = T;
T_bad.baseline_sand_fraction(T_bad.cell_id == 1) = 0.7;
T_bad.baseline_clay_fraction(T_bad.cell_id == 1) = 0.4;
assert_error(@() Resolve_Cell_Agronomic_Inputs(topo,T_bad), ...
    'Resolve_Cell_Agronomic_Inputs:TextureClosure');
T_bad = T;
T_bad.cell_x_projected_m(T_bad.cell_id == 1) = 11;
assert_error(@() Resolve_Cell_Agronomic_Inputs(topo,T_bad), ...
    'Resolve_Cell_Agronomic_Inputs:CellCoordinateMismatch');

%% Base saturation: a site baseline in every cell, a measurement where
%% used; both seed the cell's exchanger and enter its identity.
B = T;
B.baseline_bs_Ca_fraction = [0.20;0.30]; B.baseline_bs_Mg_fraction = [0.05;0.06];
B.baseline_bs_K_fraction = [0.02;0.03];
B.measured_bs_Ca_fraction = [0.50;NaN]; B.measured_bs_Mg_fraction = [0.10;NaN];
B.measured_bs_K_fraction = [0.04;NaN]; B.use_measured_base_saturation = logical([1;0]);
RB = Resolve_Cell_Agronomic_Inputs(topo,B);   % rows in grid order: cell 1 is B's row 2
assert(isequal(RB.base_saturation_source,["baseline";"measured"]));
assert(max(abs([RB.exch_Ca_cmolc_kg RB.exch_Mg_cmolc_kg RB.exch_K_cmolc_kg] - ...
    [0.30 0.06 0.03; 0.50 0.10 0.04].*[9;10]),[],'all') < 1e-12);
assert(endsWith(RB.input_identity(1),':baseline') && endsWith(RB.input_identity(2),':measured'));
assert(~contains(R.input_identity(1),'exch='));   % no base saturation, identity as before
CB = Configure_Cell_Agronomic_Inputs(MOD,table(8.6,"Kantola",NaN, ...
    'VariableNames',{'CEC_tot_in','fX_conc_convention','f_CEC_in'}),RB);
assert(CB(1).ERW_PARAM.fX_conc_convention == "f_CEC_inventory" && ...
    max(abs(str2num(char(CB(1).ERW_PARAM.f_CEC_in)) - [0.30 0.06 0.03 0 0 0.61])) < 1e-12); %#ok<ST2NM>
B_bad = B; B_bad.baseline_bs_K_fraction = [];
assert_error(@() Resolve_Cell_Agronomic_Inputs(topo,B_bad), ...
    'Resolve_Cell_Agronomic_Inputs:MissingColumns');
B_bad = B; B_bad.baseline_bs_Ca_fraction(1) = NaN;
assert_error(@() Resolve_Cell_Agronomic_Inputs(topo,B_bad), ...
    'Resolve_Cell_Agronomic_Inputs:MissingBaseline');
B_bad = B; B_bad.baseline_bs_Ca_fraction(2) = 0.95;   % cell 1, which uses its baseline
assert_error(@() Resolve_Cell_Agronomic_Inputs(topo,B_bad), ...
    'Resolve_Cell_Agronomic_Inputs:InvalidBaseSaturation');

fprintf('  cell override, fallback and validation checks passed\n');
end

function assert_error(f,expected_id)
try
    f();
catch ME
    assert(strcmp(ME.identifier,expected_id), ...
        'Expected %s, got %s: %s',expected_id,ME.identifier,ME.message);
    return
end
error('Expected error %s was not raised.',expected_id);
end
