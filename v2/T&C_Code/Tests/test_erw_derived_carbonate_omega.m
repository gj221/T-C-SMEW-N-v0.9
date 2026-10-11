function test_erw_derived_carbonate_omega
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Carbonate saturation reconstruction when Omega is not retained.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n=== test_erw_derived_carbonate_omega ===\n');
n_pass = 0;

here = fileparts(mfilename('fullpath'));
addpath(fileparts(here));

Ca = [2;4];
Mg = [1;2];
CO3 = [5;10];
[Omega_CaCO3,Omega_MgCO3] = ERW_Derived_Carbonate_Omega( ...
    Ca,Mg,CO3,10,5);
check(isequal(Omega_CaCO3,[1;4]) && isequal(Omega_MgCO3,[1;4]), ...
    'Omega = cation*CO3/K on the retained concentration basis');
check(isequal(size(Omega_CaCO3),size(Ca)) && ...
    isequal(size(Omega_MgCO3),size(Mg)), ...
    'derived saturation preserves the retained chemistry dimensions');

assert_error(@() ERW_Derived_Carbonate_Omega([1;-1],[1;1],[1;1],1,1), ...
    'ERW_Derived_Carbonate_Omega:BadConcentration');
n_pass = n_pass + 1;
assert_error(@() ERW_Derived_Carbonate_Omega([1;1],1,[1;1],1,1), ...
    'ERW_Derived_Carbonate_Omega:BadInput');
n_pass = n_pass + 1;
assert_error(@() ERW_Derived_Carbonate_Omega(1,1,1,0,1), ...
    'ERW_Derived_Carbonate_Omega:BadSolubilityProduct');
n_pass = n_pass + 1;

fprintf('  ---------------------------------------\n');
fprintf('  %d passed, 0 failed\n',n_pass);

    function check(ok,label)
        assert(ok,'test_erw_derived_carbonate_omega:Failed','%s',label);
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
    fprintf('  PASS  rejects invalid input with %s\n',expected_id);
    return
end
error('Expected error %s was not raised.',expected_id);
end
