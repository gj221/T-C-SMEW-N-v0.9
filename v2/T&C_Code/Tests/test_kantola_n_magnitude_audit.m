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

code_dir = fileparts(fileparts(mfilename('fullpath')));
parameter_source = fileread(fullfile(code_dir,'Parameters', ...
    'MOD_PARAM_Kantola_Maize_Soy.m'));
assert(~isempty(regexp(parameter_source, ...
    'Mpar_L\(1\)\.Crop_crown\s*=\s*\[1 0 1 1 0 1 1 0 1 1 0 1\]', ...
    'once')), ...
    'The Kantola crop calendar must retain maize in 2017 and 2018, with soybean in 2016.');
assert(~isempty(regexp(parameter_source, ...
    'Mpar_L\(2\)\.Crop_crown\s*=\s*\[0 1 0 0 1 0 0 1 0 0 1 0\]', ...
    'once')), ...
    'The Kantola crop calendar must retain soybean in 2016 and 2019.');

fertilizer_years = '\[2009 2011 2012 2014 2015 2017 2018 2020\]';
drivers = {fileread(fullfile(code_dir,'MAIN_FRAME_BG_ERW_10mins.m')), ...
    fileread(fullfile(code_dir,'MAIN_FRAME_BG_ERW_2D.m'))};
for k = 1:numel(drivers)
    assert(~isempty(regexp(drivers{k},fertilizer_years,'once')), ...
        'The maize fertilizer years are not aligned with the crop calendar.');
    assert(contains(drivers{k},'FertN(81:100) = 0.84') || ...
        contains(drivers{k},'FertN(81:100)=0.84') || ...
        contains(drivers{k},'FertN(81:100) = 0.84;'), ...
        'The legacy 168 kg N ha^-1 fertilizer amount or timing changed.');
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
