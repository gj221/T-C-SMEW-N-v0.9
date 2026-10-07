function [decision,results,allFailures] = Audit_Post_Emergence_Checkpoints(config)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% One-pass audit of immutable post-emergence checkpoints.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

config = validate_config(config);
rootDir = config.RootDirectory;
tag = config.Tag;
labels = config.Labels;
runs = config.Runs;
expectedIterations = config.ExpectedIterations;
expectedHashes = config.ExpectedHashes;
minimumIteration = config.MinimumIteration;
auditDir = config.AuditDirectory;
files = config.CheckpointFiles;
reportFile = fullfile(auditDir,'checkpoint_audit_report.txt');
decisionFile = fullfile(auditDir,'checkpoint_audit_decision.txt');
resultFile = fullfile(auditDir,'checkpoint_audit_result.mat');
inventoryFile = fullfile(auditDir,'matlab_checkpoint_inventory.csv');
preflightFile = fullfile(auditDir,'checkpoint_preflight.csv');
markerFile = fullfile(auditDir,'CONSUMED_LOADED_ONCE.txt');
unitsFile = fullfile(auditDir,'dimensional_analysis.csv');

assert(isfolder(auditDir),'Audit:AuditDirectory', ...
    'The audit directory does not exist: %s',auditDir);
assert(~isfile(markerFile),'Audit:ConsumedCheckpoint', ...
    'This checkpoint generation was already consumed by an audit.');
addpath(fullfile(rootDir,'v2','T&C_Code'));
addpath(fullfile(rootDir,'v2','T&C_Code','Tests'));
preflight = Audit_Checkpoint_Preflight(files,expectedHashes, ...
    expectedIterations,minimumIteration);
writetable(preflight,preflightFile);

fid = fopen(reportFile,'wt');
assert(fid >= 0,'Audit:ReportOpen','Could not open the checkpoint audit report.');
reportCleanup = onCleanup(@() fclose(fid));
fprintf(fid,'TC-SMEW-N post-emergence checkpoint audit\n');
fprintf(fid,'tag=%s\n',tag);
fprintf(fid,'minimum_iteration_exclusive=%d\n',minimumIteration);
fprintf(fid,'soil_state=P(:,1:79) [g element m^-2]\n');
fprintf(fid,'ground_state=sum(per_crown_state.*Ccrown_t)\n');
fprintf(fid,'ground_flux=sum(per_crown_flux.*Ccrown_t)\n');
fprintf(fid,'harvest_ledger=reserve_before-export-litter-reserve_after\n');
fprintf(fid,'N_ledger=inputs-outputs-storage_change [g N m^-2]\n');
fprintf(fid,'microbial_pH=PH_d=PH_bulk_d; chemistry_pH=PH_pore_d\n');
fprintf(fid,'bulk_pH_method=closed 1:2 soil:DI-water operator at 25 C\n');
fprintf(fid,'assumptions=one spatial column; daily row d+1 follows 24 hourly steps\n');
fprintf(fid,['failure_points=daily/hourly row alignment; per-crown/ground-area ' ...
    'conversion; management discontinuity mistaken for numerical instability\n\n']);

inventory = table();
schemas = cell(1,3);
for k = 1:3
    assert(isfile(files{k}),'Missing checkpoint %s.',files{k});
    info = whos('-file',files{k});
    schemas{k} = sort({info.name});
    sizeText = strings(numel(info),1);
    for q = 1:numel(info)
        sizeText(q) = mat2str(info(q).size);
    end
    block = table(repmat(runs(k),numel(info),1), ...
        repmat(string(files{k}),numel(info),1),string({info.name})', ...
        string({info.class})',sizeText,double([info.bytes])', ...
        'VariableNames',{'Run','Checkpoint','Variable','Class','Size','Bytes'});
    inventory = [inventory;block]; %#ok<AGROW>
    fprintf(fid,'%s inventory=%d expected_sha256=%s\n', ...
        labels{k},numel(info),expectedHashes{k});
end
assert(isequal(schemas{1},schemas{2}), ...
    'Checkpoint schemas differ before load.');
writetable(inventory,inventoryFile);
write_units(unitsFile);

marker = fopen(markerFile,'wt');
assert(marker >= 0,'Could not create the consumption marker.');
fprintf(marker,'tag=%s\nstate=LOAD_IN_PROGRESS\n',tag);
fclose(marker);

results = cell(3,1);
globalFailures = {};
allFailures = {};
decision = 'FAIL_KEEP_STOPPED';
loadedCount = 0;
try
    for k = 1:3
        wrapper = load(files{k});
        loadedCount = loadedCount+1;
        marker = fopen(markerFile,'at');
        fprintf(marker,'run%d_load_count=1 utc=%s\n',runs(k),datestr(now,31));
        fclose(marker);
        try
            S = find_state_struct(wrapper,0);
            Audit_Checkpoint_State_Gate(double(S.i),expectedIterations(k), ...
                minimumIteration,labels{k});
            results{k} = audit_run(S,labels{k},runs(k),auditDir,fid);
        catch runError
            globalFailures{end+1} = sprintf('%s audit error %s: %s', ...
                labels{k},runError.identifier,runError.message); %#ok<AGROW>
            fprintf(fid,'RUN %d %s AUDIT_ERROR=%s: %s\n\n', ...
                runs(k),labels{k},runError.identifier,runError.message);
        end
        clear wrapper S
    end

    assert(all(~cellfun(@isempty,results)),'Audit:IncompleteRuns', ...
        'One or more continuation runs could not be audited.');
    commonN = min(cellfun(@(x) x.n,results));
    bulkCommon = zeros(commonN,3);
    poreCommon = zeros(commonN,3);
    nh4Common = zeros(commonN,3);
    no3Common = zeros(commonN,3);
    laiCommon = zeros(commonN,3);
    plantCCommon = zeros(commonN,3);
    for k = 1:3
        bulkCommon(:,k) = results{k}.bulk(1:commonN);
        poreCommon(:,k) = results{k}.pore(1:commonN);
        nh4Common(:,k) = results{k}.nh4(1:commonN);
        no3Common(:,k) = results{k}.no3(1:commonN);
        laiCommon(:,k) = results{k}.lai(1:commonN);
        plantCCommon(:,k) = results{k}.plantC(1:commonN);
    end
    commonDates = results{1}.dates(1:commonN);
    dateAligned = all(cellfun(@(x) isequal(x.dates(1:commonN),commonDates), ...
        results));
    if ~dateAligned
        globalFailures{end+1} = ...
            'treatment state dates are not aligned'; %#ok<AGROW>
    else
        firstAmendment = datetime(2016,11,10);
        commonSpread = [Audit_Preapplication_Treatment_Spread( ...
            commonDates,bulkCommon,firstAmendment); ...
            Audit_Preapplication_Treatment_Spread( ...
            commonDates,poreCommon,firstAmendment); ...
            Audit_Preapplication_Treatment_Spread( ...
            commonDates,nh4Common,firstAmendment); ...
            Audit_Preapplication_Treatment_Spread( ...
            commonDates,no3Common,firstAmendment); ...
            Audit_Preapplication_Treatment_Spread( ...
            commonDates,laiCommon,firstAmendment); ...
            Audit_Preapplication_Treatment_Spread( ...
            commonDates,plantCCommon,firstAmendment)];
        spreadNames = {'bulk pH','pore pH','NH4','NO3','LAI','plant C'};
        for q = 1:numel(commonSpread)
            if ~isfinite(commonSpread(q)) || commonSpread(q) > 1e-10
                globalFailures{end+1} = sprintf( ...
                    'pre-application %s treatment spread is %.12g', ...
                    spreadNames{q},commonSpread(q)); %#ok<AGROW>
            end
        end
    end

    allFailures = globalFailures;
    for k = 1:3
        for q = 1:numel(results{k}.failures)
            allFailures{end+1} = sprintf('%s: %s',labels{k}, ...
                results{k}.failures{q}); %#ok<AGROW>
        end
    end
    marker = fopen(markerFile,'at');
    fprintf(marker,'state=CONSUMED_LOADED_ONCE\n');
    fprintf(marker,'note=Never reload these checkpoint versions.\n');
    fclose(marker);

    decision = 'PASS_RESUME_ALLOWED';
    if ~isempty(allFailures)
        decision = 'FAIL_KEEP_STOPPED';
    end
    write_decision(decisionFile,decision,allFailures,results);
    fprintf(fid,'\nGLOBAL_CHECKS\n');
    for q = 1:numel(commonSpread)
        fprintf(fid,'common_%s_max_spread=%.12g\n', ...
            strrep(spreadNames{q},' ','_'),commonSpread(q));
    end
    fprintf(fid,'pH_operator_runtime_identity=checked\n');
    fprintf(fid,'accepted_limitation=shared legacy first-root P/Si uptake pulse\n');
    fprintf(fid,'FAILURE_COUNT=%d\n',numel(allFailures));
    for q = 1:numel(allFailures)
        fprintf(fid,'FAILURE %d: %s\n',q,allFailures{q});
    end
    fprintf(fid,'DECISION=%s\n',decision);
    save(resultFile,'results','allFailures','decision','expectedHashes', ...
        'expectedIterations','minimumIteration','preflight','-v7');
    fprintf('CHECKPOINT_AUDIT decision=%s failures=%d report=%s\n', ...
        decision,numel(allFailures),reportFile);
catch ME
    marker = fopen(markerFile,'at');
    if marker >= 0
        fprintf(marker,'state=CONSUMED_LOADED_ONCE_AUDIT_FAILED\n');
        fprintf(marker,'loaded_count=%d\n',loadedCount);
        fprintf(marker,'error=%s: %s\n',ME.identifier,ME.message);
        fprintf(marker,'note=Never rerun or reload any checkpoint already counted above.\n');
        fclose(marker);
    end
    failures = {sprintf('%s: %s',ME.identifier,ME.message)};
    allFailures = failures;
    decision = 'FAIL_KEEP_STOPPED';
    write_decision(decisionFile,decision,failures,results);
    fprintf(fid,'FATAL_AUDIT_ERROR=%s: %s\n',ME.identifier,ME.message);
    fprintf(fid,'DECISION=FAIL_KEEP_STOPPED\n');
    fprintf('CHECKPOINT_AUDIT_FAILED %s: %s\n',ME.identifier,ME.message);
end
clear reportCleanup
end

function R = audit_run(S,label,runNumber,auditDir,fid)
required = {'i','P','Ccrown_t','PH_d','PH_bulk_d','PH_pore_d', ...
    'B_H','B_L','LAI_H','LAI_L','PHE_S_L','dflo_L','NPP_H','NPP_L', ...
    'RA_H','RA_L','Nuptake_H','Nuptake_L','FertN_d','Min_N', ...
    'N2Oflx','N2Oflx_nit','N2flx','DENITflx','Datam','AgrHarNut'};
missing = required(~isfield_many(S,required));
assert(isempty(missing),'Audit:MissingFields','Missing fields: %s', ...
    strjoin(missing,', '));

nd = floor(double(S.i)/24);
n = nd+1;
startDate = datetime(double(S.Datam(1,1)),double(S.Datam(1,2)), ...
    double(S.Datam(1,3)));
dates = startDate+caldays((0:nd)');
R = struct();
R.label = label;
R.run = runNumber;
R.iteration = double(S.i);
R.completed_days = nd;
R.current_date = dates(end);
R.n = n;
R.failures = {};

P = pool_matrix(S.P,n);
[poolNames,poolElements] = pool_metadata();
R = require_check(R,size(P,2) == 79,'retained state does not contain 79 pools');
R = require_check(R,all(isfinite(P),'all'),'soil pools are nonfinite');
R = require_check(R,min(P,[],'all') >= -1e-10,'soil pool is negative');
write_pool_stats(P,poolNames,poolElements, ...
    fullfile(auditDir,sprintf('run_%d_%s_soil_pool_stats.csv',runNumber,label)));
[groupNames,groupTotals] = pool_group_totals(P);
write_group_stats(groupTotals,groupNames, ...
    fullfile(auditDir,sprintf('run_%d_%s_element_group_stats.csv',runNumber,label)));

bulk = column_series(S.PH_bulk_d,n,'PH_bulk_d');
pore = column_series(S.PH_pore_d,n,'PH_pore_d');
response = column_series(S.PH_d,n,'PH_d');
R.bulk = bulk;
R.pore = pore;
R.bulk_pH_current = bulk(end);
R.pore_pH_current = pore(end);
R.bulk_pH_min = min(bulk);
R.bulk_pH_max = max(bulk);
R.pore_pH_min = min(pore);
R.pore_pH_max = max(pore);
R.bulk_pH_max_daily_change = max(abs(diff(bulk)));
R.pore_pH_max_daily_change = max(abs(diff(pore)));
R.pH_response_identity_error = max(abs(response-bulk));
R = require_check(R,all(isfinite([bulk;pore;response])),'daily pH is nonfinite');
R = require_check(R,all(bulk > 3 & bulk < 9), ...
    sprintf('bulk pH outside [3,9]: [%.6g %.6g]',min(bulk),max(bulk)));
R = require_check(R,all(pore > 2.5 & pore < 10), ...
    sprintf('pore pH outside [2.5,10]: [%.6g %.6g]',min(pore),max(pore)));
R = require_check(R,R.pH_response_identity_error < 1e-12, ...
    sprintf('microbial response pH differs from bulk pH by %.6g', ...
    R.pH_response_identity_error));
R = require_check(R,R.bulk_pH_max_daily_change < 1.5, ...
    sprintf('bulk pH daily jump is %.6g',R.bulk_pH_max_daily_change));
R = require_check(R,R.pore_pH_max_daily_change < 1.5, ...
    sprintf('pore pH daily jump is %.6g',R.pore_pH_max_daily_change));

nh4 = P(:,31);
no3 = P(:,32);
R.nh4 = nh4;
R.no3 = no3;
R.dates = dates;
R.NH4 = series_stats(nh4);
R.NO3 = series_stats(no3);
R.mineralN = series_stats(nh4+no3);
R.NH4_current_kgNha = 10*nh4(end);
R.NO3_current_kgNha = 10*no3(end);
R = require_check(R,all(isfinite([nh4;no3])) && min([nh4;no3]) >= -1e-12, ...
    'NH4 or NO3 is nonfinite or negative');
R = require_check(R,max(abs(diff(nh4))) <= 5, ...
    sprintf('NH4 daily change exceeds 5 g N m^-2: %.6g',max(abs(diff(nh4)))));
R = require_check(R,max(abs(diff(no3))) <= 5, ...
    sprintf('NO3 daily change exceeds 5 g N m^-2: %.6g',max(abs(diff(no3)))));

fert = daily_total(S.FertN_d,n);
mineralization = daily_total(S.Min_N,n);
is_2016 = year(dates) == 2016;
R.fertilizer_gNm2 = sum(fert(is_2016));
R.mineralization_gNm2 = sum(mineralization(2:end));
R = require_check(R,abs(R.fertilizer_gNm2) < 1e-12, ...
    sprintf('2016 soybean received %.6g g N m^-2 fertilizer',R.fertilizer_gNm2));

w = double(S.Ccrown_t(1:n,:,:));
R = require_check(R,all(isfinite(w),'all') && min(w,[],'all') >= -1e-12 && ...
    max(w,[],'all') <= 1+1e-12,'Ccrown_t is outside [0,1]');
lai = weighted_scalar(S.LAI_H,w,n)+weighted_scalar(S.LAI_L,w,n);
BH = weighted_components(S.B_H,w,n);
BL = weighted_components(S.B_L,w,n);
B = BH+BL;
plantC = sum(B(:,1:6),2);
leafC = B(:,1);
stemC = B(:,2);
rootC = B(:,3);
reserveC = B(:,4);
reproductiveC = B(:,5);
R.lai = lai;
R.plantC = plantC;
R.LAI_current = lai(end);
[growthEnd,harvestRow,harvested,harvestFailure] = ...
    Audit_Managed_Crop_Season(dates,lai,plantC,datetime(2016,9,12));
R = require_check(R,isempty(harvestFailure),harvestFailure);
R.scheduled_harvest_observed = harvested;
R.growth_endpoint_date = dates(growthEnd);
R.harvest_state_date = NaT;
if harvested
    R.harvest_state_date = dates(harvestRow);
end
R.LAI_growth_endpoint = lai(growthEnd);
R.plant_C_growth_endpoint_gCm2 = plantC(growthEnd);
[R.LAI_max,peakLai] = max(lai(1:growthEnd));
[R.plant_C_max_gCm2,peakPlant] = max(plantC(1:growthEnd));
R.LAI_peak_date = dates(peakLai);
R.plant_C_peak_date = dates(peakPlant);
R.plant_C_current_gCm2 = plantC(end);
R.root_C_max_gCm2 = max(rootC);
R.reproductive_C_max_gCm2 = max(reproductiveC);
R.reserve_C_max_gCm2 = max(reserveC);

emerge = find(plantC > 1e-12,1,'first');
leafOn = find(lai > 1e-12,1,'first');
R = require_check(R,~isempty(emerge),'managed soybean did not emerge');
R = require_check(R,~isempty(leafOn),'managed soybean LAI did not become positive');
if isempty(emerge); emerge = n; end
if isempty(leafOn); leafOn = n; end
R.emergence_state_date = dates(emerge);
R.leaf_on_state_date = dates(leafOn);
R = require_check(R,R.emergence_state_date >= datetime(2016,6,9) && ...
    R.emergence_state_date <= datetime(2016,6,10), ...
    sprintf('retained emergence row is %s',datestr(R.emergence_state_date,29)));
R = require_check(R,R.leaf_on_state_date >= datetime(2016,6,9) && ...
    R.leaf_on_state_date <= datetime(2016,6,11), ...
    sprintf('retained LAI leaf-on row is %s',datestr(R.leaf_on_state_date,29)));

activePft = active_low_pft_at(S.LAI_L,w,n,peakLai);
R.active_low_pft = activePft;
R.PHE_S_at_emergence = pft_value_at(S.PHE_S_L,activePft,n,emerge);
R.dflo_at_emergence = pft_value_at(S.dflo_L,activePft,n,emerge);
R.dflo_max = max_pft_value(S.dflo_L,activePft,n);
R = require_check(R,activePft == 2, ...
    sprintf('peak managed low-canopy PFT is %d rather than soybean PFT 2',activePft));
R = require_check(R,R.PHE_S_at_emergence >= 2, ...
    sprintf('PHE_S at emergence is %.6g',R.PHE_S_at_emergence));
R = require_check(R,R.dflo_at_emergence >= 1 && R.dflo_max >= 1, ...
    sprintf('dflo failed to initiate/progress: %.6g/%.6g', ...
    R.dflo_at_emergence,R.dflo_max));
R = require_check(R,R.LAI_max >= 0.05 && R.LAI_max <= 8, ...
    sprintf('post-emergence LAI outside [0.05,8]: %.6g',R.LAI_max));
R = require_check(R,R.LAI_growth_endpoint > 0 && ...
    R.LAI_growth_endpoint >= 0.5*R.LAI_max, ...
    sprintf('LAI collapsed before harvest: endpoint %.6g peak %.6g', ...
    R.LAI_growth_endpoint,R.LAI_max));
R = require_check(R,R.LAI_peak_date >= datetime(2016,6,10) && ...
    R.LAI_peak_date <= R.growth_endpoint_date, ...
    sprintf('LAI peak date is %s',datestr(R.LAI_peak_date,29)));
R = require_check(R,nnz(diff(lai(leafOn:growthEnd)) > 1e-8) >= 5, ...
    'LAI did not develop across at least five retained daily states');
R = require_check(R,leafC(min(emerge+1,n)) > 0, ...
    'reserve-to-leaf transfer did not create leaf carbon');
R = require_check(R,stemC(growthEnd) > 0,'stem carbon is absent before harvest');
R = require_check(R,rootC(emerge) > 0 && rootC(growthEnd) > rootC(emerge), ...
    'root carbon was absent or failed to develop');
R = require_check(R,reserveC(emerge) > 0,'reserve carbon was absent at emergence');
R = require_check(R,R.plant_C_growth_endpoint_gCm2 >= 5 && ...
    R.plant_C_max_gCm2 <= 2000, ...
    sprintf('pre-harvest plant C outside [5,2000] g C m^-2: %.6g', ...
    R.plant_C_growth_endpoint_gCm2));
R = require_check(R,R.plant_C_growth_endpoint_gCm2 >= ...
    0.5*R.plant_C_max_gCm2, ...
    sprintf('plant carbon collapsed before harvest: endpoint %.6g peak %.6g', ...
    R.plant_C_growth_endpoint_gCm2,R.plant_C_max_gCm2));
R = require_check(R,max([0;abs(diff(plantC(1:growthEnd)))]) <= 50, ...
    sprintf('plant C daily change exceeds 50 g C m^-2: %.6g', ...
    max([0;abs(diff(plantC(1:growthEnd)))])));
R = require_check(R,all(isfinite(B),'all') && min(B,[],'all') >= -1e-10, ...
    'live or reserve biomass is nonfinite or negative');

deadLai = weighted_scalar(S.LAIdead_H,w,n)+weighted_scalar(S.LAIdead_L,w,n);
rootLength = weighted_scalar(S.Rrootl_H,w,n)+weighted_scalar(S.Rrootl_L,w,n);
anpp = weighted_scalar(S.ANPP_H,w,n)+weighted_scalar(S.ANPP_L,w,n);
nppi = weighted_scalar(S.NPPI_H,w,n)+weighted_scalar(S.NPPI_L,w,n);
fixation = daily_total(S.BfixN,n);
psi = [double(S.Psi_s_H(:));double(S.Psi_s_L(:)); ...
    double(S.Psi_x_H(:));double(S.Psi_x_L(:)); ...
    double(S.Psi_l_H(:));double(S.Psi_l_L(:))];
R.dead_LAI_current = deadLai(end);
R.root_length_current = rootLength(end);
R.root_length_growth_endpoint = rootLength(growthEnd);
R.root_length_max = max(rootLength(1:growthEnd));
R.ANPP_current = anpp(end);
R.NPPI_current = nppi(end);
R.cumulative_fixation_gNm2 = sum(fixation(2:end));
R.water_potential_min_MPa = min(psi);
R.water_potential_max_MPa = max(psi);
R = require_check(R,all(isfinite([deadLai;rootLength;anpp;nppi;fixation])) && ...
    min([deadLai;rootLength;fixation]) >= -1e-10, ...
    'senescence, root, productivity or fixation diagnostic is nonphysical');
R = require_check(R,R.root_length_growth_endpoint > 0 && ...
    R.root_length_growth_endpoint >= rootLength(emerge), ...
    'fine-root length is absent or declining below its emergence state');
R = require_check(R,all(isfinite(psi)) && min(psi) >= -100 && max(psi) <= 1, ...
    sprintf('plant water potential outside [-100,1] MPa: [%.6g %.6g]', ...
    min(psi),max(psi)));
R = require_check(R,isfield(S,'RfL_Zs') && isfield(S,'ZR95_L'), ...
    'root depth or layer-distribution state is missing');
if isfield(S,'RfL_Zs') && isfield(S,'ZR95_L')
    rootFraction = double(S.RfL_Zs);
    rootDepth = double(S.ZR95_L);
    R = require_check(R,all(isfinite(rootFraction),'all') && ...
        min(rootFraction,[],'all') >= -1e-12 && ...
        sum(rootFraction(activePft,:)) > 0.99 && ...
        sum(rootFraction(activePft,:)) < 1.01, ...
        'soybean root-layer fractions are nonphysical or do not sum to one');
    R = require_check(R,all(isfinite(rootDepth)) && ...
        rootDepth(activePft) > 0, ...
        'soybean root depth is nonfinite or nonpositive');
    R.root_depth_code_units = rootDepth(activePft);
end

npp = weighted_scalar(S.NPP_H,w,n)+weighted_scalar(S.NPP_L,w,n);
ra = weighted_scalar(S.RA_H,w,n)+weighted_scalar(S.RA_L,w,n);
gpp = npp+ra;
transpiration = hourly_transpiration(S,w,nd,n);
R.GPP_max_gCm2d = max(gpp);
R.NPP_max_gCm2d = max(npp);
R.cumulative_transpiration_mm = sum(transpiration(2:end));
R = require_check(R,all(isfinite([npp;ra;gpp;transpiration])), ...
    'vegetation C or transpiration flux is nonfinite');
R = require_check(R,R.GPP_max_gCm2d > 0 && R.GPP_max_gCm2d <= 50, ...
    sprintf('peak GPP outside (0,50] g C m^-2 d^-1: %.6g',R.GPP_max_gCm2d));
R = require_check(R,R.NPP_max_gCm2d > 0 && R.NPP_max_gCm2d <= 40, ...
    sprintf('peak NPP outside (0,40] g C m^-2 d^-1: %.6g',R.NPP_max_gCm2d));
R = require_check(R,min(ra) >= -1e-10,'autotrophic respiration is negative');
R = require_check(R,R.cumulative_transpiration_mm > 0 && ...
    R.cumulative_transpiration_mm <= 1000, ...
    sprintf('seasonal transpiration outside (0,1000] mm: %.6g', ...
    R.cumulative_transpiration_mm));

[stores,reserves,uptake,totalTransfer,soilLitterTransfer,organExport, ...
    elementNames,stateFailure] = ...
    plant_element_dynamics(S,w,n);
R.plant_store_current = stores(end,:);
R.plant_store_max = max(stores,[],1);
R.plant_reserve_current = reserves(end,:);
R.plant_reserve_max = max(reserves,[],1);
R.cumulative_uptake = sum(uptake(2:end,:),1);
R.cumulative_total_plant_transfer = sum(totalTransfer(2:end,:),1);
R.cumulative_soil_litter_transfer = sum(soilLitterTransfer(2:end,:),1);
R.cumulative_organ_export = sum(organExport(2:end,:),1);
harvestExport = double(S.AgrHarNut(1:n,1:6));
R.cumulative_reserve_export = sum(harvestExport(2:end,:),1);
R.cumulative_harvest_export = R.cumulative_organ_export+ ...
    R.cumulative_reserve_export;
R = require_check(R,isempty(stateFailure),stateFailure);
R = require_check(R,all(isfinite(uptake),'all') && min(uptake,[],'all') >= -1e-10, ...
    'plant nutrient uptake is nonfinite or negative');
R = require_check(R,all(isfinite(harvestExport),'all') && ...
    min(harvestExport,[],'all') >= -1e-10, ...
    'harvest nutrient export is nonfinite or negative');
R = require_check(R,all(isfinite(organExport),'all') && ...
    min(organExport,[],'all') >= -1e-10, ...
    'organ nutrient export is nonfinite or negative');
storeScreens = [50 20 50 50 20 50];
for q = 1:6
    R = require_check(R,R.plant_store_max(q) <= storeScreens(q), ...
        sprintf('plant %s store exceeds broad screen: %.6g g m^-2', ...
        elementNames{q},R.plant_store_max(q)));
end
R = require_check(R,R.cumulative_uptake(1) > 0 && R.cumulative_uptake(1) <= 50, ...
    sprintf('seasonal N uptake outside (0,50] g N m^-2: %.6g', ...
    R.cumulative_uptake(1)));

[N,nLedger] = audit_n_ledger(S,w,P,dates,n);
R.N_closure_residual_uncorrected_gNm2 = N.uncorrected_residual;
R.N_closure_residual_gNm2 = N.corrected_residual;
R.N_daily_closure_max_abs_gNm2 = N.max_abs_daily_residual;
R.N_implicit_organ_export_gNm2 = N.organ_export;
R.N_reserve_export_gNm2 = N.reserve_export;
R.N_uptake_ledger_gNm2 = sum(uptake(2:end,1));
R.NO3_leak_gNm2 = N.no3_leak;
R.NH4_leak_gNm2 = N.nh4_leak;
denitrification = daily_total(S.DENITflx,n);
R.denitrification_gNm2 = sum(denitrification(2:end));
R = require_check(R,abs(R.N_closure_residual_gNm2) <= 1e-8, ...
    sprintf('corrected full N closure residual is %.12g g N m^-2', ...
    R.N_closure_residual_gNm2));
R = require_check(R,R.N_daily_closure_max_abs_gNm2 <= 1e-8, ...
    sprintf('maximum daily N closure residual is %.12g g N m^-2', ...
    R.N_daily_closure_max_abs_gNm2));
R = require_check(R,abs(R.N_implicit_organ_export_gNm2- ...
    R.cumulative_organ_export(1)) <= 1e-10, ...
    sprintf('N organ-export equation mismatch is %.12g g N m^-2', ...
    R.N_implicit_organ_export_gNm2-R.cumulative_organ_export(1)));
R = require_check(R,abs(R.N_uptake_ledger_gNm2-R.cumulative_uptake(1)) <= 1e-10, ...
    sprintf('N uptake crown conversion mismatch is %.12g g N m^-2', ...
    R.N_uptake_ledger_gNm2-R.cumulative_uptake(1)));

rows = 2:n;
denitN2O = double(S.N2Oflx(rows,:))-double(S.N2Oflx_nit(rows,:));
partition = double(S.N2flx(rows,:))+denitN2O-double(S.DENITflx(rows,:));
gas = [double(S.N2Oflx(rows,:));double(S.N2Oflx_nit(rows,:)); ...
    double(S.N2flx(rows,:));double(S.DENITflx(rows,:));denitN2O];
R.cumulative_N2O_gNm2 = sum(double(S.N2Oflx(rows,:)),'all');
R.max_N_partition_residual_gNm2d = max(abs(partition),[],'all');
R = require_check(R,all(isfinite(gas),'all') && min(gas,[],'all') >= -1e-12, ...
    'N2/N2O flux component is nonfinite or negative');
R = require_check(R,R.max_N_partition_residual_gNm2d < 1e-10, ...
    sprintf('N2 plus denitrification N2O partition residual is %.12g', ...
    R.max_N_partition_residual_gNm2d));

dailyFile = fullfile(auditDir,sprintf('run_%d_%s_daily.csv',runNumber,label));
daily = table(dates,bulk,pore,response,nh4,no3,10*nh4,10*no3,lai, ...
    plantC,leafC,stemC,rootC,reserveC,reproductiveC,gpp,npp,ra, ...
    transpiration,deadLai,rootLength,fixation,fert,mineralization,stores(:,1),stores(:,2),stores(:,3), ...
    stores(:,4),stores(:,5),stores(:,6),uptake(:,1),uptake(:,2),uptake(:,3), ...
    uptake(:,4),uptake(:,5),uptake(:,6), ...
    'VariableNames',{'Date','Bulk_pH','Pore_pH','Microbial_Response_pH', ...
    'NH4_gN_m2','NO3_gN_m2','NH4_kgN_ha','NO3_kgN_ha','LAI_m2_m2', ...
    'PlantC_gC_m2','LeafC_gC_m2','StemC_gC_m2','RootC_gC_m2', ...
    'ReserveC_gC_m2','ReproductiveC_gC_m2','GPP_gC_m2_d','NPP_gC_m2_d', ...
    'AutotrophicResp_gC_m2_d','Transpiration_mm_d','DeadLAI_m2_m2', ...
    'FineRootLength_code_units','Fixation_gN_m2_d','Fertilizer_gN_m2_d', ...
    'Mineralization_gN_m2_d','PlantN_gN_m2','PlantP_gP_m2','PlantK_gK_m2', ...
    'PlantCa_gCa_m2','PlantMg_gMg_m2','PlantSi_gSi_m2','NUptake_gN_m2_d', ...
    'PUptake_gP_m2_d','KUptake_gK_m2_d','CaUptake_gCa_m2_d', ...
    'MgUptake_gMg_m2_d','SiUptake_gSi_m2_d'});
writetable(daily,dailyFile);
ledgerFile = fullfile(auditDir,sprintf('run_%d_%s_N_ledger_daily.csv', ...
    runNumber,label));
writetable(nLedger,ledgerFile);
elementFile = fullfile(auditDir,sprintf('run_%d_%s_element_flux_daily.csv', ...
    runNumber,label));
elementDaily = table(dates,totalTransfer(:,1),totalTransfer(:,2), ...
    totalTransfer(:,3),totalTransfer(:,4),totalTransfer(:,5), ...
    totalTransfer(:,6),soilLitterTransfer(:,1),soilLitterTransfer(:,2), ...
    soilLitterTransfer(:,3),soilLitterTransfer(:,4), ...
    soilLitterTransfer(:,5),soilLitterTransfer(:,6),organExport(:,1), ...
    organExport(:,2),organExport(:,3),organExport(:,4),organExport(:,5), ...
    organExport(:,6),harvestExport(:,1),harvestExport(:,2), ...
    harvestExport(:,3),harvestExport(:,4),harvestExport(:,5), ...
    harvestExport(:,6),'VariableNames',{'Date','N_TotalPlantTransfer', ...
    'P_TotalPlantTransfer','K_TotalPlantTransfer','Ca_TotalPlantTransfer', ...
    'Mg_TotalPlantTransfer','Si_TotalPlantTransfer','N_SoilLitterTransfer', ...
    'P_SoilLitterTransfer','K_SoilLitterTransfer','Ca_SoilLitterTransfer', ...
    'Mg_SoilLitterTransfer','Si_SoilLitterTransfer','N_OrganExport', ...
    'P_OrganExport','K_OrganExport','Ca_OrganExport','Mg_OrganExport', ...
    'Si_OrganExport','N_ReserveExport','P_ReserveExport','K_ReserveExport', ...
    'Ca_ReserveExport','Mg_ReserveExport','Si_ReserveExport'});
writetable(elementDaily,elementFile);

fprintf(fid,'RUN %d %s\n',runNumber,label);
fprintf(fid,'  iteration=%d completed_days=%d date=%s\n', ...
    R.iteration,nd,datestr(R.current_date,29));
fprintf(fid,['  bulk_pH initial/current/min/max/max_delta=' ...
    '%.9g/%.9g/%.9g/%.9g/%.9g\n'],bulk(1),bulk(end),min(bulk),max(bulk), ...
    R.bulk_pH_max_daily_change);
fprintf(fid,['  pore_pH initial/current/min/max/max_delta=' ...
    '%.9g/%.9g/%.9g/%.9g/%.9g response_identity_error=%.3g\n'], ...
    pore(1),pore(end),min(pore),max(pore),R.pore_pH_max_daily_change, ...
    R.pH_response_identity_error);
fprintf(fid,['  NH4 initial/current/min/p05/p50/p95/max/max_delta=' ...
    '%.9g/%.9g/%.9g/%.9g/%.9g/%.9g/%.9g/%.9g gN_m2\n'], ...
    R.NH4.initial,R.NH4.current,R.NH4.minimum,R.NH4.p05,R.NH4.p50, ...
    R.NH4.p95,R.NH4.maximum,R.NH4.max_daily_change);
fprintf(fid,['  NO3 initial/current/min/p05/p50/p95/max/max_delta=' ...
    '%.9g/%.9g/%.9g/%.9g/%.9g/%.9g/%.9g/%.9g gN_m2\n'], ...
    R.NO3.initial,R.NO3.current,R.NO3.minimum,R.NO3.p05,R.NO3.p50, ...
    R.NO3.p95,R.NO3.maximum,R.NO3.max_daily_change);
fprintf(fid,['  soybean_PFT=%d emergence_state=%s leaf_on_state=%s ' ...
    'PHE_S=%.9g dflo_emergence/max=%.9g/%.9g\n'],activePft, ...
    datestr(R.emergence_state_date,29),datestr(R.leaf_on_state_date,29), ...
    R.PHE_S_at_emergence,R.dflo_at_emergence,R.dflo_max);
fprintf(fid,['  LAI current/peak=%.9g/%.9g peak_date=%s ' ...
    'plantC current/peak=%.9g/%.9g peak_date=%s\n'], ...
    R.LAI_current,R.LAI_max,datestr(R.LAI_peak_date,29), ...
    R.plant_C_current_gCm2,R.plant_C_max_gCm2,datestr(R.plant_C_peak_date,29));
fprintf(fid,['  root_length_current/growth_endpoint/max=%.9g/%.9g/%.9g ' ...
    'dead_LAI=%.9g ' ...
    'fixation=%.9g water_potential_min/max=%.9g/%.9g MPa\n'], ...
    R.root_length_current,R.root_length_growth_endpoint, ...
    R.root_length_max,R.dead_LAI_current, ...
    R.cumulative_fixation_gNm2,R.water_potential_min_MPa, ...
    R.water_potential_max_MPa);
fprintf(fid,['  rootC_peak=%.9g reproductiveC_peak=%.9g reserveC_peak=%.9g ' ...
    'GPP_peak=%.9g NPP_peak=%.9g transpiration=%.9g\n'], ...
    R.root_C_max_gCm2,R.reproductive_C_max_gCm2,R.reserve_C_max_gCm2, ...
    R.GPP_max_gCm2d,R.NPP_max_gCm2d,R.cumulative_transpiration_mm);
for q = 1:6
    fprintf(fid,['  %s store_current/max=%.9g/%.9g reserve_current/max=' ...
        '%.9g/%.9g uptake=%.9g soil_litter=%.9g organ_export=%.9g ' ...
        'reserve_export=%.9g total_harvest_export=%.9g g_m2\n'], ...
        elementNames{q},R.plant_store_current(q),R.plant_store_max(q), ...
        R.plant_reserve_current(q),R.plant_reserve_max(q), ...
        R.cumulative_uptake(q),R.cumulative_soil_litter_transfer(q), ...
        R.cumulative_organ_export(q),R.cumulative_reserve_export(q), ...
        R.cumulative_harvest_export(q));
end
fprintf(fid,['  fertilizer=%.9g mineralization=%.9g NO3leak=%.9g ' ...
    'NH4leak=%.9g denit=%.9g N2O=%.9g Nclosure_uncorrected/corrected=' ...
    '%.12g/%.12g daily_max=%.12g organ_export=%.12g partition=%.12g\n'], ...
    R.fertilizer_gNm2,R.mineralization_gNm2,R.NO3_leak_gNm2, ...
    R.NH4_leak_gNm2,R.denitrification_gNm2,R.cumulative_N2O_gNm2, ...
    R.N_closure_residual_uncorrected_gNm2,R.N_closure_residual_gNm2, ...
    R.N_daily_closure_max_abs_gNm2,R.N_implicit_organ_export_gNm2, ...
    R.max_N_partition_residual_gNm2d);
fprintf(fid,'  peak_uncorrected_N_ledger_date=%s\n', ...
    datestr(N.peak_uncorrected_date,29));
fprintf(fid,'  STATUS=%s failures=%d\n\n', ...
    ternary(isempty(R.failures),'PASS','FAIL'),numel(R.failures));
end

function P = pool_matrix(value,n)
value = double(value(1:n,:,:));
sz = size(value);
assert(sz(2) >= 79,'Audit:PoolShape','Pool state has fewer than 79 columns.');
P = sum(reshape(value,n,sz(2),[]),3);
end

function [names,elements] = pool_metadata
names = { ...
    'C_litter_above_metabolic','C_litter_above_structural_nonlignin', ...
    'C_litter_above_structural_lignin','C_woody_above_nonlignin', ...
    'C_woody_above_lignin','C_litter_below_metabolic', ...
    'C_litter_below_structural_nonlignin','C_litter_below_structural_lignin', ...
    'C_SOM_POC_lignin','C_SOM_POC_cellulose','C_SOM_MOC','C_DOC_bacterial', ...
    'C_DOC_fungal','C_enzyme_POC_bacterial','C_enzyme_POC_fungal', ...
    'C_enzyme_MOC_bacterial','C_enzyme_MOC_fungal','C_bacteria','C_fungi', ...
    'C_AM_mycorrhiza','C_EM_mycorrhiza','C_earthworm', ...
    'N_litter_above','N_woody_above','N_litter_below','N_SOM','N_bacteria', ...
    'N_fungi','N_AM_mycorrhiza','N_EM_mycorrhiza','N_NH4','N_NO3','N_DON', ...
    'N_earthworm','P_litter_above','P_woody_above','P_litter_below','P_SOM', ...
    'P_bacteria','P_fungi','P_AM_mycorrhiza','P_EM_mycorrhiza','P_mineral', ...
    'P_primary','P_secondary','P_occluded','P_DOP','K_litter_above', ...
    'K_woody_above','K_litter_below','K_SOM','K_solution','K_exchangeable', ...
    'K_fixed','K_primary','Ca_litter_above','Ca_woody_above','Ca_litter_below', ...
    'Ca_SOM','Ca_solution','Ca_exchangeable','Ca_fixed','Ca_primary', ...
    'Mg_litter_above','Mg_woody_above','Mg_litter_below','Mg_SOM','Mg_solution', ...
    'Mg_exchangeable','Mg_fixed','Mg_primary','Si_litter_above','Si_woody_above', ...
    'Si_litter_below','Si_SOM','Si_solution','Si_exchangeable','Si_fixed','Si_primary'};
elements = [repmat({'C'},1,22),repmat({'N'},1,12),repmat({'P'},1,13), ...
    repmat({'K'},1,8),repmat({'Ca'},1,8),repmat({'Mg'},1,8),repmat({'Si'},1,8)];
assert(numel(names) == 79 && numel(elements) == 79,'Pool metadata is incomplete.');
end

function write_pool_stats(P,names,elements,fileName)
nPool = size(P,2);
initial = P(1,:)';
current = P(end,:)';
minimum = min(P,[],1)';
maximum = max(P,[],1)';
p05 = zeros(nPool,1); p50 = p05; p95 = p05; maxChange = p05;
for q = 1:nPool
    s = series_stats(P(:,q));
    p05(q) = s.p05; p50(q) = s.p50; p95(q) = s.p95;
    maxChange(q) = s.max_daily_change;
end
T = table((1:nPool)',string(names(:)),string(elements(:)), ...
    repmat("g element m^-2",nPool,1),initial,current,minimum,p05,p50,p95, ...
    maximum,maxChange,'VariableNames',{'Index','Pool','Element','Actual_Units', ...
    'Initial','Current','Minimum','P05','P50','P95','Maximum','Max_Daily_Change'});
writetable(T,fileName);
end

function [names,totals] = pool_group_totals(P)
names = {'C','N','P','K','Ca','Mg','Si'};
ranges = {1:22,23:34,35:47,48:55,56:63,64:71,72:79};
totals = zeros(size(P,1),numel(names));
for q = 1:numel(names)
    totals(:,q) = sum(P(:,ranges{q}),2);
end
end

function write_group_stats(values,names,fileName)
nGroup = size(values,2);
initial = values(1,:)'; current = values(end,:)';
minimum = min(values,[],1)'; maximum = max(values,[],1)';
maxChange = max(abs(diff(values)),[],1)';
T = table(string(names(:)),repmat("g element m^-2",nGroup,1), ...
    initial,current,minimum,maximum,maxChange, ...
    'VariableNames',{'Element','Actual_Units','Initial','Current', ...
    'Minimum','Maximum','Max_Daily_Change'});
writetable(T,fileName);
end

function x = column_series(value,n,name)
x = reshape(double(value(1:n,:,:)),n,[]);
assert(size(x,2) == 1,'Audit:ColumnShape','%s is not one-column data.',name);
end

function x = daily_total(value,n)
x = sum(reshape(double(value(1:n,:,:)),n,[]),2);
end

function x = weighted_scalar(value,w,n)
value = double(value(1:n,:,:));
szW = size(w); nc = szW(2); ncell = numel(w)/(n*nc);
w = reshape(w,n,nc,ncell);
value = reshape(value,n,nc,ncell);
x = zeros(n,1);
for d = 1:n
    x(d) = sum(reshape(value(d,:,:),nc,ncell).* ...
        reshape(w(d,:,:),nc,ncell),'all');
end
end

function x = weighted_components(value,w,n)
value = double(value(1:n,:,:,:));
sz = size(value); nc = sz(2); ncomp = sz(3);
ncell = numel(value)/(n*nc*ncomp);
value = reshape(value,n,nc,ncomp,ncell);
w = reshape(w,n,nc,ncell);
x = zeros(n,ncomp);
for d = 1:n
    for c = 1:ncell
        block = reshape(value(d,:,:,c),nc,ncomp);
        weight = reshape(w(d,:,c),nc,1);
        x(d,:) = x(d,:)+sum(block.*weight,1);
    end
end
end

function x = hourly_transpiration(S,w,nd,n)
x = zeros(n,1);
assert(isfield(S,'T_H') && isfield(S,'T_L') && size(S.T_H,1) >= 24*nd, ...
    'Audit:Transpiration','Hourly transpiration is unavailable.');
nc = size(w,2); ncell = numel(w)/(n*nc);
w = reshape(w,n,nc,ncell);
TH = reshape(double(S.T_H(1:24*nd,:,:)),24*nd,nc,ncell);
TL = reshape(double(S.T_L(1:24*nd,:,:)),24*nd,nc,ncell);
for d = 1:nd
    rows = (d-1)*24+(1:24);
    weight = reshape(w(d+1,:,:),1,nc,ncell);
    x(d+1) = sum((TH(rows,:,:)+TL(rows,:,:)).*weight,'all');
end
end

function pft = active_low_pft_at(lai,w,n,index)
lai = double(lai(1:n,:,:)); nc = size(w,2);
ncell = numel(w)/(n*nc);
lai = reshape(lai,n,nc,ncell); w = reshape(w,n,nc,ncell);
current = zeros(1,nc);
for c = 1:ncell
    current = current+reshape(lai(index,:,c),1,nc).* ...
        reshape(w(index,:,c),1,nc);
end
[~,pft] = max(current);
end

function value = pft_value_at(x,pft,n,index)
x = reshape(double(x(1:n,:,:)),n,size(x,2),[]);
value = mean(reshape(x(index,pft,:),1,[]));
end

function value = max_pft_value(x,pft,n)
x = reshape(double(x(1:n,:,:)),n,size(x,2),[]);
value = max(reshape(x(:,pft,:),n,[]),[],'all');
end

function [stores,reserves,uptake,totalTransfer,soilLitterTransfer, ...
    organExport,names,failure] = ...
    plant_element_dynamics(S,w,n)
names = {'N','P','K','Ca','Mg','Si'};
totalNames = {'TNIT','TPHO','TPOT','TCAL','TMAG','TSIL'};
reserveNames = {'Nreserve','Preserve','Kreserve','Careserve','Mgreserve','Sireserve'};
uptakeNames = {'Nuptake','Puptake','Kuptake','Cauptake','Mguptake','Siuptake'};
transferNames = {'TexN','TexP','TexK','TexCa','TexMg','TexSi'};
isoilIndex = {10:12,13:15,16:18,19:21,22:24,25:27};
stores = zeros(n,6); reserves = stores; uptake = stores;
totalTransfer = stores; soilLitterTransfer = stores; organExport = stores;
failure = '';
for q = 1:6
    fields = {[totalNames{q} '_H'],[totalNames{q} '_L'], ...
        [reserveNames{q} '_H'],[reserveNames{q} '_L'], ...
        [uptakeNames{q} '_H'],[uptakeNames{q} '_L'], ...
        [transferNames{q} '_H'],[transferNames{q} '_L']};
    if ~all(isfield_many(S,fields))
        failure = sprintf('NOT ASSESSABLE: missing plant %s field',names{q});
        return
    end
    totalH = double(S.(fields{1})(1:n,:,:));
    totalL = double(S.(fields{2})(1:n,:,:));
    reserveH = double(S.(fields{3})(1:n,:,:));
    reserveL = double(S.(fields{4})(1:n,:,:));
    rawState = [totalH(:);totalL(:);reserveH(:);reserveL(:)];
    if any(~isfinite(rawState)) || min(rawState) < -1e-10
        failure = sprintf('plant %s organ or reserve pool is nonfinite or negative',names{q});
        return
    end
    stores(:,q) = weighted_scalar(totalH,w,n)+weighted_scalar(totalL,w,n)+ ...
        weighted_scalar(reserveH,w,n)+weighted_scalar(reserveL,w,n);
    reserves(:,q) = weighted_scalar(reserveH,w,n)+weighted_scalar(reserveL,w,n);
    uptake(:,q) = weighted_scalar(S.(fields{5}),w,n)+ ...
        weighted_scalar(S.(fields{6}),w,n);
    totalTransfer(:,q) = weighted_scalar(S.(fields{7}),w,n)+ ...
        weighted_scalar(S.(fields{8}),w,n);
    soilLitterTransfer(:,q) = weighted_pending_element( ...
        S.ISOIL_H,w,n,isoilIndex{q})+weighted_pending_element( ...
        S.ISOIL_L,w,n,isoilIndex{q});
    organExport(:,q) = totalTransfer(:,q)-soilLitterTransfer(:,q);
end
end

function [A,T] = audit_n_ledger(S,w,P,dates,n)
soil = sum(P(:,23:34),2);
plant = weighted_scalar(S.TNIT_H,w,n)+weighted_scalar(S.TNIT_L,w,n);
reserve = weighted_scalar(S.Nreserve_H,w,n)+ ...
    weighted_scalar(S.Nreserve_L,w,n);
pending = weighted_pending_element(S.ISOIL_H,w,n,10:12)+ ...
    weighted_pending_element(S.ISOIL_L,w,n,10:12);
storage = soil+plant+reserve+pending;

fertilizer = daily_total(S.FertN_d,n);
deposition = daily_total(S.DepN_d,n);
fixation = daily_total(S.BfixN,n);
seed = daily_total(S.SeedInput(1:n,2),n);
inputs = fertilizer+deposition+fixation+seed;

n2o = daily_total(S.N2Oflx,n);
n2 = daily_total(S.N2flx,n);
volatilization = daily_total(S.VOL,n);
no3Leak = daily_total(S.LEAK_NO3,n);
nh4Leak = daily_total(S.LEAK_NH4,n);
donLeak = daily_total(S.LEAK_DON,n);
reserveExport = daily_total(S.AgrHarNut(1:n,1),n);
baseOutputs = n2o+n2+volatilization+no3Leak+nh4Leak+donLeak+reserveExport;

totalTransfer = weighted_scalar(S.TexN_H,w,n)+ ...
    weighted_scalar(S.TexN_L,w,n);
organExport = totalTransfer-pending;
assert(all(isfinite(organExport)) && min(organExport) >= -1e-10, ...
    'Audit:OrganExport','Derived organ N export is nonfinite or negative.');
outputs = baseOutputs+organExport;

uncorrectedDaily = zeros(n,1);
correctedDaily = zeros(n,1);
uncorrectedDaily(2:end) = inputs(2:end)-baseOutputs(2:end)-diff(storage);
correctedDaily(2:end) = inputs(2:end)-outputs(2:end)-diff(storage);
uncorrectedCumulative = cumsum(uncorrectedDaily);
correctedCumulative = cumsum(correctedDaily);

A.uncorrected_residual = uncorrectedCumulative(end);
A.corrected_residual = correctedCumulative(end);
A.max_abs_daily_residual = max(abs(correctedDaily));
A.organ_export = sum(organExport(2:end));
A.reserve_export = sum(reserveExport(2:end));
A.no3_leak = sum(no3Leak(2:end));
A.nh4_leak = sum(nh4Leak(2:end));
[~,peakRow] = max(abs(uncorrectedDaily));
A.peak_uncorrected_date = dates(peakRow);

allFlux = [inputs;baseOutputs;organExport;outputs];
assert(all(isfinite(allFlux)) && min(allFlux) >= -1e-10, ...
    'Audit:NFlux','N boundary flux is nonfinite or negative.');
T = table(dates,soil,plant,reserve,pending,storage,fertilizer,deposition, ...
    fixation,seed,inputs,n2o,n2,volatilization,no3Leak,nh4Leak,donLeak, ...
    reserveExport,totalTransfer,organExport,baseOutputs,outputs, ...
    uncorrectedDaily,correctedDaily,uncorrectedCumulative, ...
    correctedCumulative,'VariableNames',{'Date','SoilN','PlantOrganN', ...
    'PlantReserveN','PendingLitterN','TotalStorageN','FertilizerN', ...
    'DepositionN','FixationN','SeedN','TotalInputN','N2O','N2', ...
    'VolatilizationN','NO3LeakN','NH4LeakN','DONLeakN','ReserveExportN', ...
    'TotalPlantTransferN','OrganHarvestExportN','BaseOutputN', ...
    'CorrectedOutputN','UncorrectedDailyResidualN', ...
    'CorrectedDailyResidualN','UncorrectedCumulativeResidualN', ...
    'CorrectedCumulativeResidualN'});
end

function value = weighted_pending_element(x,w,n,indices)
x = double(x(1:n,:,indices,:));
sz = size(x); nc = sz(2); ncomp = numel(indices);
ncell = numel(x)/(n*nc*ncomp);
x = reshape(x,n,nc,ncomp,ncell);
w = reshape(double(w),n,nc,ncell);
value = zeros(n,1);
for d = 1:n
    for c = 1:ncell
        block = reshape(x(d,:,:,c),nc,ncomp);
        weight = reshape(w(d,:,c),nc,1);
        value(d) = value(d)+sum(sum(block,2).*weight);
    end
end
end

function S = series_stats(x)
x = double(x(:)); sx = sort(x);
S.initial = x(1); S.current = x(end); S.minimum = min(x); S.maximum = max(x);
S.p05 = sample_quantile(sx,0.05); S.p50 = sample_quantile(sx,0.50);
S.p95 = sample_quantile(sx,0.95); S.max_daily_change = max(abs(diff(x)));
end

function q = sample_quantile(x,p)
position = 1+(numel(x)-1)*p; lo = floor(position); hi = ceil(position);
if lo == hi
    q = x(lo);
else
    q = x(lo)+(position-lo)*(x(hi)-x(lo));
end
end

function y = row_range(x)
y = max(x,[],2)-min(x,[],2);
end

function mask = isfield_many(S,names)
mask = false(size(names));
for k = 1:numel(names)
    mask(k) = isfield(S,names{k});
end
end

function R = require_check(R,condition,message)
if ~condition
    R.failures{end+1} = message;
end
end

function X = find_state_struct(value,depth)
assert(depth <= 8,'Audit:CheckpointNesting','Checkpoint nesting exceeds eight levels.');
if isstruct(value) && isscalar(value) && isfield(value,'P') && isfield(value,'i')
    X = value;
    return
end
assert(isstruct(value) && isscalar(value),'Audit:CheckpointBranch', ...
    'Checkpoint branch is not a scalar struct.');
names = fieldnames(value);
for k = 1:numel(names)
    child = value.(names{k});
    if isstruct(child) && isscalar(child)
        try
            X = find_state_struct(child,depth+1);
            return
        catch ME
            if ~startsWith(ME.identifier,'Audit:Checkpoint')
                rethrow(ME)
            end
        end
    end
end
error('Audit:CheckpointState','Checkpoint state containing P and i was not found.');
end

function write_units(fileName)
variable = ["pH";"LAI";"Ccrown_t";"soil nutrient pool"; ...
    "daily nutrient flux";"plant biomass";"plant nutrient store"; ...
    "per-crown nutrient uptake";"ground-area nutrient uptake"; ...
    "harvest nutrient export";"transpiration";"fine-root length"; ...
    "root layer fraction";"root depth";"plant water potential"; ...
    "biological N fixation"];
expectedSI = ["dimensionless";"m2 leaf m-2 ground";"dimensionless"; ...
    "kg element m-2";"kg element m-2 s-1";"kg C m-2 ground"; ...
    "kg element m-2 ground";"kg element m-2 crown s-1"; ...
    "kg element m-2 ground s-1";"kg element m-2 event-1"; ...
    "kg water m-2 s-1";"m root m-2 crown";"dimensionless"; ...
    "m";"Pa";"kg N m-2 s-1"];
actual = ["dimensionless";"m2 leaf m-2 ground";"fraction [0,1]"; ...
    "g element m-2; multiply by 1e-3 for SI"; ...
    "g element m-2 d-1; multiply by 1e-3/86400 for SI"; ...
    "g C m-2 crown; Ccrown-weighted and multiply by 1e-3"; ...
    "g element m-2 crown; Ccrown-weighted and multiply by 1e-3"; ...
    "g element m-2 crown d-1; multiply by 1e-3/86400"; ...
    "g element m-2 d-1 after Ccrown weighting; multiply by 1e-3/86400"; ...
    "g element m-2 event-1; multiply by 1e-3 for SI"; ...
    "mm h-1; 1 mm water m-2 equals 1 kg water m-2"; ...
    "m root m-2 crown";"fraction [0,1]"; ...
    "mm; multiply by 1e-3 for SI";"MPa; multiply by 1e6 for SI"; ...
    "g N m-2 d-1; multiply by 1e-3/86400 for SI"];
T = table(variable,expectedSI,actual, ...
    'VariableNames',{'Variable','Expected_SI_Units','Actual_Units_Provided_In_Code'});
writetable(T,fileName);
end

function write_decision(fileName,decision,failures,results)
fid = fopen(fileName,'wt');
assert(fid >= 0,'Could not write the audit decision.');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid,'decision=%s\n',decision);
fprintf(fid,'failure_count=%d\n',numel(failures));
for k = 1:numel(failures)
    fprintf(fid,'failure_%d=%s\n',k,failures{k});
end
for k = 1:numel(results)
    if ~isempty(results{k}) && isfield(results{k},'label')
        fprintf(fid,'%s_status=%s\n',results{k}.label, ...
            ternary(isempty(results{k}.failures),'PASS','FAIL'));
    end
end
clear cleanup
end

function value = ternary(condition,yesValue,noValue)
if condition
    value = yesValue;
else
    value = noValue;
end
end

function config = validate_config(config)
required = {'RootDirectory','Tag','AuditDirectory','Labels','Runs', ...
    'ExpectedIterations','ExpectedHashes','CheckpointFiles','MinimumIteration'};
assert(isstruct(config) && isscalar(config),'Audit:Config', ...
    'The audit configuration must be a scalar struct.');
missing = required(~isfield_many(config,required));
assert(isempty(missing),'Audit:Config', ...
    'Missing configuration fields: %s',strjoin(missing,', '));

for name = {'RootDirectory','Tag','AuditDirectory'}
    value = config.(name{1});
    assert(ischar(value) || (isstring(value) && isscalar(value)), ...
        'Audit:Config','%s must be text.',name{1});
    config.(name{1}) = char(value);
end

assert(iscellstr(config.Labels) || isstring(config.Labels), ...
    'Audit:Config','Labels must be text values.');
assert(iscellstr(config.ExpectedHashes) || isstring(config.ExpectedHashes), ...
    'Audit:Config','ExpectedHashes must be text values.');
assert(iscellstr(config.CheckpointFiles) || isstring(config.CheckpointFiles), ...
    'Audit:Config','CheckpointFiles must be text values.');
config.Labels = cellstr(string(config.Labels(:)'));
config.ExpectedHashes = cellstr(string(config.ExpectedHashes(:)'));
config.CheckpointFiles = cellstr(string(config.CheckpointFiles(:)'));
config.Runs = double(config.Runs(:)');
config.ExpectedIterations = double(config.ExpectedIterations(:)');

assert(numel(config.Runs) == 3 && numel(config.Labels) == 3 && ...
    numel(config.ExpectedIterations) == 3 && ...
    numel(config.ExpectedHashes) == 3 && ...
    numel(config.CheckpointFiles) == 3,'Audit:Config', ...
    'The audit requires control, basalt and calcite configuration entries.');
assert(all(isfinite(config.Runs)) && all(config.Runs == fix(config.Runs)), ...
    'Audit:Config','Run identifiers must be finite integers.');
end
