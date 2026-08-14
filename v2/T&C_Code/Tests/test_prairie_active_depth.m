function test_prairie_active_depth
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Prairie biogeochemical depth follows the deepest active root profile.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n=== test_prairie_active_depth ===\n');
n_pass = 0;
Zs = [0 10 50 100 150 200 300 400 600 800 1000 1250 1500 ...
      1750 2000 2500 3000 3500 4000 4500 5000];

[Zbio,Bio_Zs] = Biogeochemical_Depth_From_Roots( ...
    Zs,[2500 0],[1200 0]);
check(Zbio == 2500, 'active depth equals the 2500 mm big bluestem root depth');
check(abs(sum(Bio_Zs)-1) < 1e-12, ...
    'active-layer weights sum to one');
layer_thickness = Bio_Zs*Zbio;
check(abs(sum(layer_thickness)-Zbio) < 1e-12 && ...
      all(layer_thickness(Zs(1:end-1) >= Zbio) == 0), ...
    'active thickness ends exactly at the root-depth boundary');

assert_error(@() Biogeochemical_Depth_From_Roots( ...
    Zs,[5500 0],[1200 0]), ...
    'Biogeochemical_Depth_From_Roots:RootBelowMesh');
n_pass = n_pass + 1;

fprintf('  ---------------------------------------\n');
fprintf('  %d passed, 0 failed\n',n_pass);

    function check(ok,label)
        assert(ok,'test_prairie_active_depth:Failed','%s',label);
        n_pass = n_pass + 1;
        fprintf('  PASS  %s\n',label);
    end

end

function assert_error(f,expected_id)
try
    f();
catch ME
    assert(strcmp(ME.identifier,expected_id), ...
        'Expected %s, got %s: %s',expected_id,ME.identifier,ME.message);
    fprintf('  PASS  rejects a root depth below the soil mesh\n');
    return
end
error('Expected error %s was not raised.',expected_id);
end
