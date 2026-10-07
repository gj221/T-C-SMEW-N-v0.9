function Audit_Checkpoint_State_Gate(actual_iteration,expected_iteration, ...
    minimum_iteration,label)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Confirm the loaded state is the recorded post-emergence checkpoint.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

values = [actual_iteration expected_iteration minimum_iteration];
if any(~isfinite(values)) || any(values ~= fix(values)) || ...
        actual_iteration < 0 || expected_iteration < 0 || minimum_iteration < 0
    error('Audit:IterationState','Checkpoint iterations must be nonnegative integers.');
end
if ~(ischar(label) || (isstring(label) && isscalar(label)))
    error('Audit:IterationLabel','The checkpoint label must be text.');
end
label = char(label);
if actual_iteration ~= expected_iteration
    error('Audit:IterationMismatch', ...
        '%s loaded iteration %d differs from recorded iteration %d.', ...
        label,actual_iteration,expected_iteration);
end
if actual_iteration <= minimum_iteration
    error('Audit:BeforeEmergenceGate', ...
        '%s checkpoint iteration %d has not crossed %d.', ...
        label,actual_iteration,minimum_iteration);
end
end
