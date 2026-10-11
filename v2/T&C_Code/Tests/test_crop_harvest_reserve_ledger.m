function test_crop_harvest_reserve_ledger
% Harvested reserves must accumulate once on a plot-area basis.

code_dir = fileparts(fileparts(mfilename('fullpath')));
addpath(code_dir);

one_d = fileread(fullfile(code_dir,'MAIN_FRAME_BG_ERW_10mins.m'));
two_d = fileread(fullfile(code_dir,'MAIN_FRAME_BG_ERW_2D.m'));
for source = {one_d,two_d}
    assert(contains(source{1},'AgrHarNut_cc'), ...
        'Crop calls must not overwrite the daily harvested-reserve ledger.');
    assert(contains(source{1},'Plant_Harvest_Reserve_Export'), ...
        'Crop reserve export must be converted from crown to plot area.');
end

initial = [0.1 0.2 0.3 0.4 0.5 0.6];
reserve = [1 2 3 4 5 6];
actual = Plant_Harvest_Reserve_Export(initial,reserve,0.25);
expected = initial + 0.25*reserve;
assert(max(abs(actual-expected)) < 10*eps, ...
    'Harvested reserve mass must accumulate without loss or duplication.');

St = Veg_Stoichiometric_Parameter(20);
Mpar = Vegetation_Management_Parameter();
reserve = [0.068574062635 0.025 0.03 0.014478941890 0.001 0.03];
fab = 0.4;
out_zero = run_harvest_case(St,Mpar,fab,reserve,zeros(1,8));
expected_reserve = (1-fab)*reserve;
assert(max(abs([out_zero{15:20}]-expected_reserve)) < 1e-12, ...
    'Harvest must retain only the non-litter reserve share.');
assert(max(abs([out_zero{2:7}]-fab*reserve)) < 1e-12, ...
    'Harvest must route the litter reserve share through Tex.');
assert(max(abs([out_zero{2:7}] + [out_zero{15:20}]-reserve)) < 1e-12, ...
    'Harvest reserve routing must close for every represented element.');

St = Veg_Stoichiometric_Parameter(12,'Soybean');
reserve = [2 1 2 2 1 2];
r = [1.4 0.7 1.2 0.6 1.3 0.8];
Btm1 = [4 3 2 1 5 6 0 0];
ratios = [St.Nl St.Ns St.Nr St.Ns St.Nf St.Nh; ...
    St.Phol St.Phos St.Phor St.Phos St.Phof St.Phoh; ...
    St.Kpotl St.Kpots St.Kpotr St.Kpots St.Kpotf St.Kpoth; ...
    St.Call St.Cals St.Calr St.Cals St.Calf St.Calh; ...
    St.Magl St.Mags St.Magr St.Mags St.Magf St.Magh; ...
    St.Sill St.Sils St.Silr St.Sils St.Silf St.Silh];
initial = sum(Btm1(1:6)./ratios,2)' + reserve;

RB = [0.4 0.3 0.2 0.1 0.5 0.6 0 0];
for ManI = [0 -1 -5]
    [transfer,state] = run_management_case(St,Mpar,Btm1-RB,Btm1,RB, ...
        reserve,r,ManI);
    assert(max(abs(transfer+state-initial)) < 1e-12, ...
        'Partial management must conserve every represented element.');
end

[transfer,state] = run_management_case(St,Mpar,zeros(1,8),Btm1,Btm1, ...
    reserve,r,-2);
assert(all(transfer >= 0) && all(state >= 0), ...
    'Complete harvest must retain finite nonnegative element inventories.');
assert(max(abs(transfer+state-initial)) < 1e-12, ...
    'Complete harvest must conserve base tissue inventory and reserve.');

Btm1 = [0 0 10 0 0 0 0 0];
reserve = 2*ones(1,6);
r = 1.5*ones(1,6);
initial = Btm1(3)./ratios(:,3)' + reserve;
[transfer,state,litter] = run_management_case( ...
    St,Mpar,zeros(1,8),Btm1,Btm1,reserve,r,-2);
organ_export = transfer-litter;
assert(all(organ_export >= -1e-12), ...
    'Complete harvest litter must not exceed the physical element transfer.');
assert(max(abs(litter+organ_export+state-initial)) < 1e-12, ...
    'Complete harvest litter, export and retained state must close.');

overdraw_failed = false;
try
    B = [0 0 5 0 0 0 0 0];
    RB = [0 0 5 0 0 0 0 0];
    run_management_case(St,Mpar,B,Btm1,RB,zeros(1,6),r,0);
catch ME
    overdraw_failed = strcmp(ME.identifier,'Plant_Exports:ReserveOverdraw');
end
assert(overdraw_failed, ...
    'Ordinary plant removal must reject a nutrient-reserve overdraw.');

Mpar.fire_eff = 0.5;
fab = 0.74;
RB = [1 0 0 0 0 0 0 0];
out = cell(1,27);
[out{:}] = Plant_Exports(Btm1-RB,Btm1,zeros(1,6),0,0,0,0,0, ...
    zeros(1,3),St,Mpar,fab,1-fab,RB,reserve(1),reserve(2), ...
    reserve(3),reserve(4),reserve(5),reserve(6),1,1,1,1,1,1, ...
    -5,1,'Kantola');
expected_reserve = (1-Mpar.fire_eff*fab)*reserve;
assert(max(abs([out{15:20}]-expected_reserve)) < 1e-12, ...
    'The fire indicator must reach the reserve-export boundary.');

fprintf('PASS: crop reserve harvest is a plot-area boundary export.\n');
end

function out = run_harvest_case(St,Mpar,fab,reserve,RB)
B = zeros(1,8);
B(4) = 1;
out = cell(1,27);
[out{:}] = Plant_Exports(B,B,zeros(1,6),0,0,0,0,0,zeros(1,3), ...
    St,Mpar,fab,1-fab,RB,reserve(1),reserve(2),reserve(3), ...
    reserve(4),reserve(5),reserve(6),1,1,1,1,1,1,-2,1,'Kantola');
end

function [transfer,state,litter] = run_management_case(St,Mpar,B,Btm1,RB, ...
    reserve,r,ManI)
out = cell(1,27);
[out{:}] = Plant_Exports(B,Btm1,zeros(1,6),0,0,0,0,0,zeros(1,3), ...
    St,Mpar,1,0,RB,reserve(1),reserve(2),reserve(3),reserve(4), ...
    reserve(5),reserve(6),r(1),r(2),r(3),r(4),r(5),r(6), ...
    ManI,1,'Kantola');
transfer = [out{2:7}];
state = [out{8:13}] + [out{15:20}];
soil_input = out{27};
litter = [sum(soil_input(10:12)) sum(soil_input(13:15)) ...
    sum(soil_input(16:18)) sum(soil_input(19:21)) ...
    sum(soil_input(22:24)) sum(soil_input(25:27))];
end
