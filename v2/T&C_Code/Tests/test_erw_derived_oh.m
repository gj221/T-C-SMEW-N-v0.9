function test_erw_derived_oh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Hydroxide reconstruction when OH is not retained in an ERW snapshot.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n=== test_erw_derived_oh ===\n');
n_pass = 0;

here = fileparts(mfilename('fullpath'));
addpath(fileparts(here));

H = [1e-2;1;10];
k_w = [100;100;100];
OH = ERW_Derived_OH(H,k_w);
check(isequal(OH,[1e4;100;10]), ...
    'OH = k_w/H without requiring a retained OH field');
check(isequal(size(OH),size(H)), ...
    'derived OH preserves the retained chemistry dimensions');

assert_error(@() ERW_Derived_OH([1;0],[100;100]), ...
    'ERW_Derived_OH:BadHydrogen');
n_pass = n_pass + 1;
assert_error(@() ERW_Derived_OH([1;2],100), ...
    'ERW_Derived_OH:BadInput');
n_pass = n_pass + 1;
assert_error(@() ERW_Derived_OH(1,-100), ...
    'ERW_Derived_OH:BadIonProduct');
n_pass = n_pass + 1;

fprintf('  ---------------------------------------\n');
fprintf('  %d passed, 0 failed\n',n_pass);

    function check(ok,label)
        assert(ok,'test_erw_derived_oh:Failed','%s',label);
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
