function test_crop_harvest_scale_reserve_export
St = Veg_Stoichiometric_Parameter(20);
Mpar = Vegetation_Management_Parameter();
reserve = [1.1 0.22 0.7 0.1 0.08 0.15];
fab = 0.5;

Btm1 = zeros(1,8);
Btm1(1:6) = [60 25 46 39 37 10];
pre = run_case(St,Mpar,Btm1,Btm1,zeros(1,8),reserve,fab,0);
harvest = run_case(St,Mpar,zeros(1,8),Btm1,Btm1,reserve,fab,-2);

pre_inventory = cell2mat(pre(8:13))+reserve;
transfer = cell2mat(harvest(2:7));
reserve_after = cell2mat(harvest(15:20));
if any(abs(transfer+reserve_after-pre_inventory)>1e-10)
    error('test_crop_harvest_scale_reserve_export:Closure', ...
        'The complete-harvest organ and reserve boundary does not close.')
end
if any(abs(reserve_after-reserve)>1e-12)
    error('test_crop_harvest_scale_reserve_export:Reserve', ...
        'Plant_Exports must retain the whole reserve for the off-field AgrHarNut export.')
end
IS = harvest{27};
indices = {10:12,13:15,16:18,19:21,22:24,25:27};
soil_litter = cellfun(@(idx) sum(IS(idx)),indices);
organ_export = transfer-soil_litter;
if any(~isfinite([transfer reserve_after soil_litter organ_export])) || ...
        any([transfer reserve_after soil_litter organ_export] < -1e-12)
    error('test_crop_harvest_scale_reserve_export:Physicality', ...
        'Harvest transfer, litter, export and reserve must be nonnegative.')
end

fprintf('PASS: complete harvest closes without reserve overdraw.\n');
end

function out = run_case(St,Mpar,B,Btm1,RB,reserve,fab,ManI)
out = cell(1,27);
[out{:}] = Plant_Exports(B,Btm1,zeros(1,6),0,0,0,0,0,zeros(1,3), ...
    St,Mpar,fab,1-fab,RB,reserve(1),reserve(2),reserve(3), ...
    reserve(4),reserve(5),reserve(6),0.3,0.3,0.3,0.3,0.3,0.3, ...
    ManI,1,'Kantola');
end
