function test_kantola_reserve_initialization
% Bare post-spin-up crops must not contain unaccounted nutrient reserves.

test_dir = fileparts(mfilename('fullpath'));
code_dir = fileparts(test_dir);
parameter_file = fullfile(code_dir,'Parameters','MOD_PARAM_Kantola_Maize_Soy.m');
code = fileread(parameter_file);

reserve_names = {'Nreserve_H','Preserve_H','Kreserve_H', ...
    'Careserve_H','Mgreserve_H','Sireserve_H', ...
    'Nreserve_L','Preserve_L','Kreserve_L', ...
    'Careserve_L','Mgreserve_L','Sireserve_L'};

for k = 1:numel(reserve_names)
    pattern = sprintf('%s\\(1,:\\)\\s*=\\s*\\[([^\\]]+)\\]',reserve_names{k});
    token = regexp(code,pattern,'tokens','once');
    assert(~isempty(token),'Missing initialization for %s.',reserve_names{k});
    values = sscanf(token{1},'%f').';
    assert(~isempty(values) && all(values == 0), ...
        '%s must start at zero g element m^-2 PFT.',reserve_names{k});
end

reserve = zeros(1,2);
crown_old = [1 1];
crown_new = [0.05 0.02];
scaled = reserve.*crown_old./crown_new;
assert(all(scaled == 0), ...
    'Crown rescaling must not create nutrient mass from a bare crop state.');

fprintf('PASS: Kantola crop reserves start at zero and remain zero under crown rescaling.\n');
end
