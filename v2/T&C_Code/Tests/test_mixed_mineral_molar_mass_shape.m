function test_mixed_mineral_molar_mass_shape
% Mixed feedstocks retain one molar mass per mineral.

test_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(test_dir,'..'));

mineral = char('albite','ferroactinolite','epidote', ...
    'chlorite','quartz','calcite');
MM_col = [263;865.8;474;620;60;100.087] / 1e6;
MM = ERW_Mineral_Molar_Mass_Row(MM_col,mineral);

assert(isequal(size(MM),[1 6]),'Molar masses must follow mineral columns.');
assert(isequal(MM,MM_col.'),'Mineral molar masses changed during reshaping.');

M0 = [233 119 178 340 90 26];
EW = [2 3 5 7 11 13] * 1e3;
dtw = 1/144;
M1 = M0 - EW.*MM*dtw;
expected = M0 - EW.*MM_col.'*dtw;
assert(max(abs(M1-expected)) < 1e-12, ...
    'Mixed-mineral mass loss used the wrong molar mass.');

rho_rock = 3e6;
Wr = EW / 10;
rock_f = M0 / 1000;
d_shrink = sum(Wr.*rock_f.*MM/rho_rock);
expected_shrink = sum(Wr.*rock_f.*MM_col.'/rho_rock);
assert(abs(d_shrink-expected_shrink) < 1e-18, ...
    'Particle shrinkage used the wrong mineral molar mass or phase area.');

unweighted_shrink = sum(Wr.*MM/rho_rock);
assert(d_shrink < unweighted_shrink, ...
    'A minor fast phase must not retreat the entire mixed particle surface.');

source = fileread(fullfile(test_dir,'..','ERW_function.m'));
assert(contains(source, ...
    'Wr(j_only,:).*rock_f(j_only,:).*MM_min_row/rho_rock'), ...
    'Production particle retreat must use the same phase-area weighting.');
assert(contains(source, ...
    'psd_0_in(ind_check).*(d(j_only,ind_check,t_appl)./d_0(ind_check)).^3.*'), ...
    'Each PSD class must use elementwise cubic particle shrinkage.');

assert_error(@() ERW_Mineral_Molar_Mass_Row(MM_col(1:5),mineral), ...
    'ERW:MineralMolarMassSize');

fprintf('PASS: mixed-mineral molar masses preserve mineral columns.\n');
end

function assert_error(fn,id)
try
    fn();
catch ME
    assert(strcmp(ME.identifier,id), ...
        'Expected %s, received %s.',id,ME.identifier);
    return
end
error('Expected error %s was not raised.',id);
end
