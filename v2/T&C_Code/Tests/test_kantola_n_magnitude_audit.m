function test_kantola_n_magnitude_audit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% test_kantola_n_magnitude_audit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

S.P = zeros(4,79);
S.P(:,31) = [1; 1.9; 2.8; 3.8];
S.TNIT_H = zeros(4,1);
S.TNIT_L = zeros(4,1);
S.Nreserve_H = zeros(4,1);
S.Nreserve_L = zeros(4,1);
S.FertN_d = [0; ones(3,1)];
S.DepN_d = zeros(4,1);
S.BfixN = zeros(4,1);
S.N2Oflx = zeros(4,1);
S.N2flx = zeros(4,1);
S.VOL = zeros(4,1);
S.LEAK_NO3 = [0; 0.1*ones(3,1)];
S.LEAK_NH4 = zeros(4,1);
S.LEAK_DON = zeros(4,1);
S.NN = 72;
S.NNd = 4;
S.Lk = (10/24)*ones(72,1);
S.Lk(1) = 10; %% the first block's window is the single hour before it
S.Date = datenum(2020,1,1) + (0:71)'/24;
S.AgrHarNut = zeros(4,6);
S.SeedInput = zeros(4,7);
S.SeedInput(2,2) = 0.1;
S.Nuptake_H = [0; 0.2; 0.3; 0.4];
S.Nuptake_L = zeros(4,1);
S.DENITflx = [0; 0.01; 0.02; 0.03];
S.Min_N = [0; 0.4; 0.5; 0.6];
S.Ccrown_t = ones(4,1);
S.ISOIL_H = zeros(4,1,27);
S.ISOIL_L = zeros(4,1,27);
S.TexN_H = zeros(4,1);
S.TexN_L = zeros(4,1);

A = Kantola_N_Magnitude_Audit(S,'control',3);

assert(abs(A.standing_NH4_gNm2-3.8) < 1e-12);
assert(abs(A.standing_NH4_kgNha-38) < 1e-12);
assert(abs(A.cumulative_NO3_leak_gNm2-0.3) < 1e-12);
assert(abs(A.cumulative_NO3_leak_kgNha-3) < 1e-12);
assert(abs(A.maximum_NO3_mobile_mgNL-10) < 1e-12);
assert(abs(A.plant_uptake_gNm2-0.9) < 1e-12);
assert(abs(A.denitrification_gNm2-0.06) < 1e-12);
assert(abs(A.fertilizer_gNm2-3) < 1e-12);
assert(abs(A.mineralization_gNm2-1.5) < 1e-12);
assert(abs(A.seed_input_gNm2-0.1) < 1e-12);
assert(abs(A.N_closure_residual_gNm2) < 1e-12);
assert(A.observed_2017_NO3_leak_kgNha == 67.06);
assert(isnan(A.modeled_2017_NO3_leak_kgNha));

S_calendar = calendar_case();
A_calendar = Kantola_N_Magnitude_Audit(S_calendar,'control',731);
assert(abs(A_calendar.modeled_2017_NO3_leak_kgNha-36.5)<1e-12, ...
    'Calendar 2017 leakage was not selected from the retained forcing dates.');

S_area = S;
area_day = size(S_area.P,1)-1;
S_area.Ccrown_t = 0.5*ones(size(S_area.TNIT_H));
S_area.TNIT_H(:) = 0;
S_area.TNIT_L(:) = 0;
S_area.Nreserve_H(:) = 0;
S_area.Nreserve_L(:) = 0;
plant_crown = 0.4;
plant_ground = 0.5*plant_crown;
S_area.TNIT_H(area_day+1,1) = plant_crown;
S_area.P(area_day+1,31) = S_area.P(area_day+1,31)-plant_ground;
A_area = Kantola_N_Magnitude_Audit(S_area,'control',area_day);
assert(abs(A_area.N_closure_residual_gNm2)<1e-12, ...
    'Per-crown plant N was not converted once to ground-area storage.');
assert(abs(A_area.plant_uptake_gNm2-0.45)<1e-12, ...
    'Per-crown uptake was not converted once to ground-area flux.');

S_pending = S;
S_pending.Ccrown_t = 0.5*ones(size(S_pending.TNIT_H));
S_pending.ISOIL_H(end,1,10) = 0.2;
S_pending.TexN_H(end,1) = 0.2;
S_pending.P(end,31) = S_pending.P(end,31)-0.1;
A_pending = Kantola_N_Magnitude_Audit(S_pending,'control',3);
assert(abs(A_pending.N_closure_residual_gNm2)<1e-12, ...
    'Pending litter N was not retained until the next soil update.');

S_export = S;
S_export.TexN_H(end) = 0.2;
S_export.P(end,31) = S_export.P(end,31)-0.2;
A_export = Kantola_N_Magnitude_Audit(S_export,'control',3);
assert(abs(A_export.organ_harvest_export_gNm2-0.2)<1e-12 && ...
    abs(A_export.N_closure_residual_gNm2)<1e-12, ...
    'Organ harvest N was omitted from the external N boundary ledger.');

S_imm = S;
S_imm.Min_N(3) = -0.2;
A_imm = Kantola_N_Magnitude_Audit(S_imm,'control',3);
assert(abs(A_imm.mineralization_gNm2-0.8)<1e-12 && ...
    abs(A_imm.net_immobilization_gNm2-0.2)<1e-12 && ...
    A_imm.net_immobilization_days == 1 && ...
    abs(A_imm.N_closure_residual_gNm2)<1e-12, ...
    'Net immobilization days must be accepted as signed net mineralization.');
assert(A.net_immobilization_gNm2 == 0 && A.net_immobilization_days == 0, ...
    'A trajectory without immobilization must report none.');
S_dilute = S;
S_dilute.TNIT_H(end) = 0.5;
S_dilute.Nreserve_H(end) = -0.2; %% tissue diluted to 0.6 of nominal
S_dilute.P(end,31) = S_dilute.P(end,31)-0.3;
A_dilute = Kantola_N_Magnitude_Audit(S_dilute,'control',3);
assert(abs(A_dilute.N_closure_residual_gNm2)<1e-12, ...
    'A reserve deficit is diluted tissue and must enter plant N signed.');
S_dilute.Nreserve_H(end) = -0.6;
S_dilute.P(end,31) = S_dilute.P(end,31)+0.4;
try
    Kantola_N_Magnitude_Audit(S_dilute,'control',3);
    error('Kantola_N_Magnitude_Audit:TestFailure', ...
        'A deficit larger than the tissue N was accepted.');
catch ME
    assert(strcmp(ME.identifier,'Kantola_N_Magnitude_Audit:InvalidState'), ...
        'Negative plant N must be refused as InvalidState.');
end
S_shift = S;
S_shift.Lk(26:49) = 0; %% the third block's window (hours 26 to 49) is dry
try
    Kantola_N_Magnitude_Audit(S_shift,'control',3);
    error('Kantola_N_Magnitude_Audit:TestFailure', ...
        'Leaching after a dry previous day was accepted.');
catch ME
    assert(strcmp(ME.identifier,'Kantola_N_Magnitude_Audit:DryLeakage'), ...
        'Leaching on the day after a dry day must be refused as DryLeakage.');
end
S_shift.LEAK_NO3(4) = 0;
S_shift.P(4,31) = S_shift.P(4,31)+0.1;
A_shift = Kantola_N_Magnitude_Audit(S_shift,'control',3);
assert(abs(A_shift.N_closure_residual_gNm2)<1e-12 && ...
    abs(A_shift.maximum_NO3_mobile_mgNL-10)<1e-12 && ...
    abs(A_shift.latest_NO3_mobile_mgNL-10)<1e-12, ...
    'Leaching must be matched with the drainage window of its daily block.');
S_window = S;
S_window.Lk(:) = 0;
S_window.Lk(25) = 10; %% 00:00 of day 2: inside the window of block 2 only
S_window.LEAK_NO3 = [0; 0; 0.1; 0];
S_window.P(:,31) = [1; 2; 3; 4];
A_window = Kantola_N_Magnitude_Audit(S_window,'control',3);
assert(abs(A_window.N_closure_residual_gNm2)<1e-12 && ...
    abs(A_window.maximum_NO3_mobile_mgNL-10)<1e-12, ...
    'A midnight-hour drainage belongs to the following daily block.');
S_calendar_shift = S_window;
S_calendar_shift.Lk(:) = 0;
S_calendar_shift.Lk(60) = 10; %% 11:00 of day 3: after the last block, outside every window
try
    Kantola_N_Magnitude_Audit(S_calendar_shift,'control',3);
    error('Kantola_N_Magnitude_Audit:TestFailure', ...
        'Leaching driven by no block-window drainage was accepted.');
catch ME
    assert(strcmp(ME.identifier,'Kantola_N_Magnitude_Audit:DryLeakage'), ...
        'Leaching without block-window drainage must be refused.');
end
S_init = S;
S_init.Nreserve_L(1) = 1000;
S_init.Ccrown_t(1) = 1;
A_init = Kantola_N_Magnitude_Audit(S_init,'control',3);
assert(abs(A_init.N_closure_residual_gNm2)<1e-12 && ...
    abs(A_init.row1_placeholder_plant_N_gNm2-1000)<1e-12, ...
    'The row-1 legacy seed reserve must not enter the boundary ledger.');
assert(abs(A.row1_placeholder_plant_N_gNm2)<1e-12, ...
    'A fixture without a row-1 placeholder must report none.');
S_live = S_init;
S_live.Nreserve_L(2) = 1000;
S_live.Ccrown_t(2) = 1;
try
    Kantola_N_Magnitude_Audit(S_live,'control',3);
    error('Kantola_N_Magnitude_Audit:TestFailure', ...
        'A live legacy reserve under a present crop was accepted.');
catch ME
    assert(strcmp(ME.identifier,'Kantola_N_Magnitude_Audit:LegacyReserveLive'), ...
        'A live legacy reserve under a present crop must be refused.');
end
S_nan = S;
S_nan.Min_N(3) = NaN;
try
    Kantola_N_Magnitude_Audit(S_nan,'control',3);
    error('Kantola_N_Magnitude_Audit:TestFailure', ...
        'A nonfinite net mineralization day was accepted.');
catch ME
    assert(strcmp(ME.identifier,'Kantola_N_Magnitude_Audit:InvalidFlux'), ...
        'A nonfinite net mineralization day must be refused as InvalidFlux.');
end

code_dir = fileparts(fileparts(mfilename('fullpath')));
parameter_source = fileread(fullfile(code_dir,'Parameters', ...
    'MOD_PARAM_Kantola_Maize_Soy.m'));
calendar = regexp(parameter_source,'rotation_years = (\d+):(\d+);','tokens','once');
assert(~isempty(calendar) && str2double(calendar{1}) == 2009 && ...
    str2double(calendar{2}) >= 2050, ...
    'The Kantola crop calendar must run from 2009 to at least 2050.');
assert(contains(parameter_source,'soybean_year = mod(rotation_years-2010,3) == 0;') && ...
    contains(parameter_source,'Mpar_L(1).Crop_crown = double(~soybean_year);') && ...
    contains(parameter_source,'Mpar_L(2).Crop_crown = double(soybean_year);'), ...
    'The Kantola rotation must be maize-maize-soybean, soybean every third year from 2010.');
rotation_years = 2009:2050;
soybean_year = mod(rotation_years-2010,3) == 0;
assert(isequal(soybean_year(1:12),logical([0 1 0 0 1 0 0 1 0 0 1 0])), ...
    'The calendar must keep the 2009-2020 record: soybean in 2010, 2013, 2016 and 2019.');
assert(contains(parameter_source,'B_IO.FertN_years = rotation_years(~soybean_year);'), ...
    'The maize fertilizer years must follow the crop calendar.');

drivers = {fileread(fullfile(code_dir,'MAIN_FRAME_BG_ERW_10mins.m')), ...
    fileread(fullfile(code_dir,'MAIN_FRAME_BG_ERW_2D.m'))};
for k = 1:numel(drivers)
    assert(contains(drivers{k},'FertN_years)') && ...
        isempty(regexp(drivers{k},'ismember\(Datam\(i,1\),\[\d','once')), ...
        'The drivers must take the maize fertilizer years from the crop calendar.');
    assert(contains(drivers{k},'FertN(81:100) = 1.01') || ...
        contains(drivers{k},'FertN(81:100)=1.01') || ...
        contains(drivers{k},'FertN(81:100) = 1.01;'), ...
        'The legacy 202 kg N ha^-1 fertilizer amount or timing changed.');
end

fprintf('PASS: Kantola N audit preserves units, flux sums and full N closure.\n');
end

function S = calendar_case()
nd = 731;
leak = [0.02*ones(366,1); 0.01*ones(365,1)];
S.P = zeros(nd+1,79);
S.P(:,32) = [100; 100-cumsum(leak)];
S.TNIT_H = zeros(nd+1,1);
S.TNIT_L = zeros(nd+1,1);
S.Nreserve_H = zeros(nd+1,1);
S.Nreserve_L = zeros(nd+1,1);
S.FertN_d = zeros(nd+1,1);
S.DepN_d = zeros(nd+1,1);
S.BfixN = zeros(nd+1,1);
S.N2Oflx = zeros(nd+1,1);
S.N2flx = zeros(nd+1,1);
S.VOL = zeros(nd+1,1);
S.LEAK_NO3 = [0; leak];
S.LEAK_NH4 = zeros(nd+1,1);
S.LEAK_DON = zeros(nd+1,1);
S.Lk = ones(24*nd,1);
S.Date = datenum(2016,1,1) + (0:24*nd-1)'/24;
S.AgrHarNut = zeros(nd+1,6);
S.SeedInput = zeros(nd+1,7);
S.Nuptake_H = zeros(nd+1,1);
S.Nuptake_L = zeros(nd+1,1);
S.DENITflx = zeros(nd+1,1);
S.Min_N = zeros(nd+1,1);
S.Ccrown_t = ones(nd+1,1);
S.ISOIL_H = zeros(nd+1,1,27);
S.ISOIL_L = zeros(nd+1,1,27);
S.TexN_H = zeros(nd+1,1);
S.TexN_L = zeros(nd+1,1);
end
