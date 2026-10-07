function test_kantola_reserve_initialization
% Legacy T&C initialises crops with large non-limiting nutrient reserves, so
% growth is carbon- and water-limited rather than nutrient-limited. This test
% pins those legacy values rather than requiring zero.

test_dir = fileparts(mfilename('fullpath'));
code_dir = fileparts(test_dir);
parameter_file = fullfile(code_dir,'Parameters','MOD_PARAM_Kantola_Maize_Soy.m');
code = fileread(parameter_file);

reserve_names = {'Nreserve_H','Preserve_H','Kreserve_H', ...
    'Careserve_H','Mgreserve_H','Sireserve_H', ...
    'Nreserve_L','Preserve_L','Kreserve_L', ...
    'Careserve_L','Mgreserve_L','Sireserve_L'};
expected = [0 0 0 1000 1000 1000 1000 100 1000 1000 1000 1000];

for k = 1:numel(reserve_names)
    % Built by concatenation, not sprintf: sprintf consumes backslash escapes
    % and would silently mangle the regular expression.
    pattern = [reserve_names{k} '\(1,:\)\s*=\s*\[([^\]]+)\]'];
    token = regexp(code,pattern,'tokens','once');
    assert(~isempty(token),'Missing initialization for %s.',reserve_names{k});
    values = sscanf(token{1},'%f').';
    assert(~isempty(values) && all(values == expected(k)), ...
        '%s must start at the legacy %g g element m^-2 PFT.', ...
        reserve_names{k},expected(k));
end

reserve = zeros(1,2);
crown_old = [1 1];
crown_new = [0.05 0.02];
scaled = reserve.*crown_old./crown_new;
assert(all(scaled == 0), ...
    'Crown rescaling must not create nutrient mass from a bare crop state.');

fprintf('PASS: Kantola crop reserves match the legacy initialization.\n');
end
