function tests = test_plant_exports_standing_leaf_nutrients
tests = functiontests(localfunctions);
end

function test_near_complete_turnover_is_nonnegative_and_conservative(testCase)
St = Veg_Stoichiometric_Parameter(20);
NuLITtm1 = 1./leaf_ratios(St);
[NuLIT, exported] = run_case(St,0.01,0.99,0,NuLITtm1);

verifyGreaterThanOrEqual(testCase,NuLIT,zeros(1,6));
verifyEqual(testCase,NuLIT,0.01*NuLITtm1,'AbsTol',1e-12);
verifyEqual(testCase,NuLIT + exported,NuLITtm1,'AbsTol',1e-12);
end

function test_zero_inventory_does_not_create_nutrients(testCase)
St = Veg_Stoichiometric_Parameter(20);
[NuLIT, exported] = run_case(St,0.01,0.99,0,zeros(1,6));

verifyEqual(testCase,NuLIT,zeros(1,6),'AbsTol',1e-12);
verifyEqual(testCase,exported,zeros(1,6),'AbsTol',1e-12);
end

function test_complete_turnover_exports_available_nutrients(testCase)
St = Veg_Stoichiometric_Parameter(20);
NuLITtm1 = 1./leaf_ratios(St);
Sll = 0.2;
[NuLIT, exported] = run_case(St,0,1 + Sll,Sll,NuLITtm1);
NuLITadd = 0.4*Sll./leaf_ratios(St);

verifyEqual(testCase,NuLIT,zeros(1,6),'AbsTol',1e-12);
verifyEqual(testCase,exported,NuLITtm1 + NuLITadd,'AbsTol',1e-12);
end

function [NuLIT, exported] = run_case(St,Bstanding,Slf,Sll,NuLITtm1)
Mpar = Vegetation_Management_Parameter();
Btm1 = zeros(1,8);
Btm1([4 7]) = [1 Bstanding + Slf - Sll];
B = Btm1;
B(7) = Bstanding;
RB = zeros(1,8);
out = cell(1,14);
[out{:}] = Plant_Exports(B,Btm1,NuLITtm1,Slf,0,0,Sll,0,zeros(1,3), ...
    St,Mpar,0.5,0.5,RB,0,0,0,0,0,0,1,1,1,1,1,1,0,1,'Kantola');
NuLIT = out{14};
exported = [out{2} out{3} out{4} out{5} out{6} out{7}];
end

function ratios = leaf_ratios(St)
ratios = [St.Nl St.Phol St.Kpotl St.Call St.Magl St.Sill];
end
