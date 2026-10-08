function test_kantola_bulk_pH_calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Beerling bulk-pH calibration target and launch contract.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n=== test_kantola_bulk_pH_calibration ===\n');
n_pass = 0;
v2 = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(v2,'run'),fullfile(v2,'calibration'));

result_dir = fullfile(v2,'SIM','Subfield_12','calibration_isolation_test');
log_dir = fullfile(v2,'logs','Subfield_12','calibration_isolation_test');
report = Calibrate_Kantola_Bulk_pH('DryRun',true, ...
    'AlphaGrid',[0.4 0.8 1.3],'NDays',679, ...
    'ResultDirectory',result_dir,'LogDirectory',log_dir);
check(numel(report.candidates) == 3, ...
    'screen retains every requested alpha candidate');
check(isequal(size(report.targets.depth_cm),[2 2]) && ...
    isequal(report.targets.depth_cm,[0 10;10 30]), ...
    'observations retain the measured depth intervals');
check(numel(report.targets.pre_bulk_pH) == 2 && ...
    numel(report.targets.control_bulk_pH) == 2 && ...
    numel(report.targets.basalt_bulk_pH) == 2, ...
    'bulk-pH targets remain separate by depth');
check(max(abs(report.targets.delta_bulk_pH-[0.32525951 0.11805556])) < 1e-6, ...
    'calibration retains both observed pH treatment responses');
check(numel(report.targets.delta_base_sat_pct) == 2 && ...
    all(isfinite(report.targets.delta_base_sat_pct)), ...
    'base-saturation responses remain separate without unmeasured CEC weights');
check(all(contains({report.candidates.control_file},'run10_control__vertical_1d')), ...
    'every candidate includes the untreated trajectory');
check(all(contains({report.candidates.basalt_file},'run11_basalt__vertical_1d')), ...
    'every candidate includes the Beerling basalt trajectory');
check(all(contains({report.candidates.control_file},'__alpha')), ...
    'candidate output names encode alpha');
check(all(startsWith({report.candidates.control_file},result_dir)), ...
    'calibration outputs are isolated from production outputs');
check(strcmp(report.options.ResultDirectory,result_dir) && ...
    strcmp(report.options.LogDirectory,log_dir), ...
    'calibration retains its explicit result and log roots');
check(isfield(report.options,'ObservationBoundary') && ...
    strcmp(report.options.ObservationBoundary,'pre_application'), ...
    'calibration declares the observation side of the second-application boundary');
boundary_fields = {'control_pre_application_bulk_pH', ...
    'basalt_pre_application_bulk_pH','control_post_application_bulk_pH', ...
    'basalt_post_application_bulk_pH','delta_pre_application_bulk_pH', ...
    'delta_post_application_bulk_pH'};
check(all(isfield(report.candidates,boundary_fields)), ...
    'calibration reports pre- and post-second-application bulk pH separately');

assert_error(@() Calibrate_Kantola_Bulk_pH('DryRun',true, ...
    'AlphaGrid',[-0.1 0.8]),'Calibrate_Kantola_Bulk_pH:BadOptions');
n_pass = n_pass + 1;
assert_error(@() Calibrate_Kantola_Bulk_pH('DryRun',true, ...
    'NDays',400),'Calibrate_Kantola_Bulk_pH:BadOptions');
n_pass = n_pass + 1;
assert_error(@() Run_All_12_Parallel('DryRun',true,'Runs',[10 11], ...
    'NDays',679),'Run_All_12_Parallel:MissingAppeloAlpha');
n_pass = n_pass + 1;
assert_error(@() Calibrate_Kantola_Bulk_pH('DryRun',true, ...
    'LogDirectory',3),'Calibrate_Kantola_Bulk_pH:BadOptions');
n_pass = n_pass + 1;

fixture = [tempname '.mat'];
cleanup = onCleanup(@() delete_if_present(fixture));
Date = datenum(2016,11,8)+(0:47)/24;
E = struct;
E.dm_end = 288;
E.tm = (0:E.dm_end-1)'/6;
E.bg.pH_bulk = nan(E.dm_end,1);
E.bg.pH_bulk([1 144 288]) = [6.1 6.2 6.3];
E.bg.f_Na = 0.001*ones(E.dm_end,1);
E.bg.f_K = 0.040*ones(E.dm_end,1);
E.bg.f_Ca = 0.500*ones(E.dm_end,1);
E.bg.f_Mg = 0.100*ones(E.dm_end,1);
E.bg.f_H = 0.300*ones(E.dm_end,1);
E.bg.f_Al = 0.059*ones(E.dm_end,1);
ERW = {E};
save(fixture,'Date','ERW');
summary = Kantola_Bulk_pH_Model_Summary( ...
    fixture,datenum(2016,11,10),'pre_application');
check(abs(summary.pre_application_bulk_pH-6.2) < 1e-12 && ...
    abs(summary.post_application_bulk_pH-6.3) < 1e-12, ...
    'calibration samples sparse daily bulk pH on the ERW time grid');

%% The same column as a multi-square result of package y keeps it: no
%% ten-minute records, the t = 0 row and each ERW call's summary.
fixture_daily = [tempname '.mat'];
cleanup_daily = onCleanup(@() delete_if_present(fixture_daily));
names = {'pH_bulk','f_Na','f_K','f_Ca','f_Mg','f_H','f_Al'};
calls = [2 144; 145 288];
Ey = struct('res',6,'dm_end',144);
D = struct('n',2,'rows',calls,'j',[2; 3]);
f_sum = E.bg.f_Ca+E.bg.f_Mg+E.bg.f_K+E.bg.f_Na+E.bg.f_Al+E.bg.f_H;
D.max.f_sum = [max(abs(f_sum(2:144)-1)); max(abs(f_sum(145:288)-1))];
for k = 1:numel(names)
    Ey.daily_initial.bg.(names{k}) = E.bg.(names{k})(1);
    D.end.(names{k}) = E.bg.(names{k})(calls(:,2));
end
ERW = {Ey}; ERW_daily = {D};
save(fixture_daily,'Date','ERW','ERW_daily');
summary_daily = Kantola_Bulk_pH_Model_Summary( ...
    fixture_daily,datenum(2016,11,10),'pre_application');
same = true;
for f = setdiff(fieldnames(summary),{'exchange_closure_max'})'
    same = same && isequal(summary_daily.(f{1}),summary.(f{1}));
end
check(same && abs(summary_daily.exchange_closure_max - ...
    summary.exchange_closure_max) < 1e-15, ...
    'a daily-store result gives the same calibration samples');

fprintf('  ---------------------------------------\n');
fprintf('  %d passed, 0 failed\n',n_pass);

    function check(ok,label)
        assert(ok,'test_kantola_bulk_pH_calibration:Failed','%s',label);
        n_pass = n_pass + 1;
        fprintf('  PASS  %s\n',label);
    end

end

function delete_if_present(file)
if isfile(file); delete(file); end
end

function assert_error(fn,identifier)
try
    fn();
catch ME
    assert(strcmp(ME.identifier,identifier), ...
        'Expected %s, received %s.',identifier,ME.identifier);
    return
end
error('test_kantola_bulk_pH_calibration:ExpectedError', ...
    'Expected %s.',identifier);
end
