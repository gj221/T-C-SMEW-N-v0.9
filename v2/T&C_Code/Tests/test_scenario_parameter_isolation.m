function test_scenario_parameter_isolation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Treatment rows may change feedstock inputs, not the base T&C setup.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n=== test_scenario_parameter_isolation ===\n');
n_pass = 0;
MOD = table([3;5;7],[2020;2016;2020],[300;300;300], ...
    'VariableNames',{'Number','Year','Zbio'});
ERW = table([3;5;7],repmat("Kantola",3,1), ...
    ["basalt";"";"calcite"],["Beerling";"control";"calcite"], ...
    [1;NaN;1],ones(3,1), ...
    'VariableNames',{'Number','Location','mineral', ...
    'ExperimentalDescription','diss_f','opt_TC_biochem'});

[M,E,C] = Select_Subfield_Parameter_Rows(MOD,ERW,3);
check(M.Number == 5 && M.Year == 2016, ...
    'base T&C parameters always come from the control row');
check(E.Number == 3 && E.mineral == "basalt", ...
    'selected treatment row supplies the feedstock');
check(C.Number == 5 && C.mineral == "", ...
    'untreated cells use the control chemistry row');

ERW.Location(1) = "OtherSite";
assert_error(@() Select_Subfield_Parameter_Rows(MOD,ERW,3), ...
    'Select_Subfield_Parameter_Rows:NonFeedstockDifference');
n_pass = n_pass + 1;

fprintf('  ---------------------------------------\n');
fprintf('  %d passed, 0 failed\n',n_pass);

    function check(ok,label)
        assert(ok,'test_scenario_parameter_isolation:Failed','%s',label);
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
    fprintf('  PASS  rejects a non-feedstock treatment difference\n');
    return
end
error('Expected error %s was not raised.',expected_id);
end
