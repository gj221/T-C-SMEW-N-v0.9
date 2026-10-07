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

function testCalciteCasesKeepTheirApplicationRates(testCase)
code_dir = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(fileparts(code_dir),'run'));
RM = Run_Matrix_12();
is_calcite = strcmp({RM.treatment},'calcite');
verifyEqual(testCase,[RM(is_calcite).case_number],[7 7 7 7]);

T = readtable(fullfile(code_dir,'Inputs','ERW_MOD_PARAM_parameters.xlsx'), ...
    'Sheet','ERW_MOD','VariableNamingRule','preserve');
verifyEqual(testCase,T.M_rock_event(T.Number == 7),5000);
verifyEqual(testCase,T.M_rock_event(T.Number == 8),500);
end

function testBeerlingUsesFiniteCompositeDensity(testCase)
code_dir = fileparts(fileparts(mfilename('fullpath')));
input_file = fullfile(code_dir,'Inputs','ERW_MOD_PARAM_parameters.xlsx');
T = readtable(input_file,'Sheet','ERW_MOD','VariableNamingRule','preserve');
is_beerling = strcmp(string(T.ExperimentalDescription),'crop_Beerling');

verifyEqual(testCase,sum(is_beerling),1);
verifyEqual(testCase,T.rho_rock(is_beerling),3e6,'AbsTol',1e3);
end

function testLewisUsesFiniteCompositeDensity(testCase)
code_dir = fileparts(fileparts(mfilename('fullpath')));
input_file = fullfile(code_dir,'Inputs','ERW_MOD_PARAM_parameters.xlsx');
T = readtable(input_file,'Sheet','ERW_MOD','VariableNamingRule','preserve');
is_lewis = strcmp(string(T.ExperimentalDescription),'crop_Lewis');
verifyEqual(testCase,sum(is_lewis),1);
verifyEqual(testCase,T.rho_rock(is_lewis),3e6,'AbsTol',1e3);
end
