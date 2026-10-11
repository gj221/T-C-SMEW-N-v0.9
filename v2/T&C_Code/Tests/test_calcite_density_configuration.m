function tests = test_calcite_density_configuration
tests = functiontests(localfunctions);
end

function testPureCalciteUsesCalciteDensity(testCase)
code_dir = fileparts(fileparts(mfilename('fullpath')));
input_file = fullfile(code_dir,'Inputs','ERW_MOD_PARAM_parameters.xlsx');
T = readtable(input_file,'Sheet','ERW_MOD','VariableNamingRule','preserve');
is_calcite = strcmp(string(T.ExperimentalDescription),'crop_calcite');

verifyEqual(testCase,sum(is_calcite),1);
verifyEqual(testCase,T.rho_rock(is_calcite),2.71e6,'AbsTol',1e3);
end

function testBeerlingUsesFiniteCompositeDensity(testCase)
code_dir = fileparts(fileparts(mfilename('fullpath')));
input_file = fullfile(code_dir,'Inputs','ERW_MOD_PARAM_parameters.xlsx');
T = readtable(input_file,'Sheet','ERW_MOD','VariableNamingRule','preserve');
is_beerling = strcmp(string(T.ExperimentalDescription),'crop_Beerling');

verifyEqual(testCase,sum(is_beerling),1);
verifyEqual(testCase,T.rho_rock(is_beerling),3e6,'AbsTol',1e3);
end
