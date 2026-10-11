function test_erw_cation_partition_ownership()
% TEST_ERW_CATION_PARTITION_OWNERSHIP Verify ERW owns cation partitioning.

test_dir = fileparts(mfilename('fullpath'));
code_dir = fileparts(test_dir);
old_path = path;
cleanup_path = onCleanup(@() path(old_path));
addpath(code_dir);

rates = [1/2, 1/2, 1/600, 1/13500];

legacy_rates = Biogeo_Cation_Partition_Rates(rates, 0, 'appelo');
assert(isequal(legacy_rates, rates), ...
    'Standalone T&C-SMEW must retain its legacy cation partition rates.');

gt_rates = Biogeo_Cation_Partition_Rates(rates, 1, 'TC_GT');
assert(isequal(gt_rates, [0 0 rates(3:4)]), ...
    'TC_GT must retain the legacy fixed-pool kinetics in ERW mode.');

appelo_rates = Biogeo_Cation_Partition_Rates(rates, 1, 'appelo');
assert(isequal(appelo_rates, zeros(size(rates))), ...
    'Appelo must own all cation partitioning in ERW mode.');

kantola_gt = struct();
kantola_gt.excel_input = struct('fX_conc_convention', "Kantola");
kantola_gt.conv_Al = 1;
kantola_gt.ERW_TC_const_input = struct('bulk_pH_params', ...
    struct('scheme', 'TC_GT'));
assert_error(@() ERW_pre_processing(fieldnames(kantola_gt), kantola_gt), ...
    'ERW:KantolaRequiresAppelo');

kantola_missing = kantola_gt;
kantola_missing.ERW_TC_const_input.bulk_pH_params = struct();
assert_error(@() ERW_pre_processing(fieldnames(kantola_missing), kantola_missing), ...
    'ERW:KantolaRequiresAppelo');

kantola_array = kantola_gt;
kantola_array.ERW_TC_const_input.bulk_pH_params.scheme = ["appelo", "TC_GT"];
assert_error(@() ERW_pre_processing(fieldnames(kantola_array), kantola_array), ...
    'ERW:KantolaRequiresAppelo');

kantola_struct_array = kantola_gt;
kantola_struct_array.ERW_TC_const_input.bulk_pH_params = ...
    repmat(struct('scheme', 'appelo'), 1, 2);
assert_error(@() ERW_pre_processing(fieldnames(kantola_struct_array), ...
    kantola_struct_array), 'ERW:KantolaRequiresAppelo');

kantola_whitespace = kantola_gt;
kantola_whitespace.ERW_TC_const_input.bulk_pH_params.scheme = ' appelo ';
assert_error(@() ERW_pre_processing(fieldnames(kantola_whitespace), ...
    kantola_whitespace), 'ERW:KantolaRequiresAppelo');

source = fileread(fullfile(code_dir,'BIOGEOCHEMISTRY_DYNAMIC3.m'));
assert(count(source,'Biogeo_Cation_Partition_Rates(') == 3, ...
    'K, Ca and Mg must all delegate partition ownership in ERW mode.');

unit_source = fileread(fullfile(code_dir,'BIOGEO_UNIT.m'));
frame_1d = fileread(fullfile(code_dir,'MAIN_FRAME_BG_ERW_10mins.m'));
frame_2d = fileread(fullfile(code_dir,'MAIN_FRAME_BG_ERW_2D.m'));
frame_tc = fileread(fullfile(code_dir,'MAIN_FRAME_BG.m'));
assert(contains(unit_source,'WFPS,exchange_scheme,KNmaxd_override') && ...
    contains(frame_1d,'biogeo_exchange_scheme,[]') && ...
    contains(frame_2d,'biogeo_exchange_scheme{c},biogeo_KNmaxd{c}') && ...
    contains(frame_tc,"WFPS_bio,'TC_GT',[]"), ...
    'The selected exchange scheme is not wired through every driver.');

fprintf('test_erw_cation_partition_ownership: PASS\n');

    function assert_error(f, expected_id)
        try
            f();
        catch ME
            assert(strcmp(ME.identifier, expected_id), ...
                'Expected %s, got %s: %s', expected_id, ME.identifier, ME.message);
            return
        end
        error('Expected error %s was not raised.', expected_id);
    end

end
