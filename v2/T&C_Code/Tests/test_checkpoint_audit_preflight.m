function test_checkpoint_audit_preflight
% Checkpoint audits enforce immutable identity and the emergence boundary.

test_dir = fileparts(mfilename('fullpath'));
run_dir = fullfile(test_dir,'..','..','audit');
addpath(run_dir);

assert_error(@() Audit_Post_Emergence_Checkpoints(struct),'Audit:Config');

tmp = tempname;
mkdir(tmp);
cleanup = onCleanup(@() rmdir(tmp,'s')); %#ok<NASGU>
checkpoint = fullfile(tmp,'candidate.mat');
write_text(checkpoint,'abc');
known_hash = 'BA7816BF8F01CFEA414140DE5DAE2223B00361A396177A9CB410FF61F20015AD';

inventory = Audit_Checkpoint_Preflight({checkpoint},{known_hash},3841,3840);
assert(height(inventory) == 1 && inventory.ExpectedIteration(1) == 3841 && ...
    inventory.SHA256(1) == string(known_hash), ...
    'Valid post-emergence checkpoint metadata was not retained.');
Audit_Checkpoint_State_Gate(3841,3841,3840,'control');

assert_error(@() Audit_Checkpoint_Preflight( ...
    {checkpoint},{known_hash},3840,3840),'Audit:BeforeEmergenceConfig');
assert_error(@() Audit_Checkpoint_Preflight( ...
    {checkpoint},{repmat('0',1,64)},3841,3840),'Audit:HashMismatch');
assert_error(@() Audit_Checkpoint_State_Gate(3840,3840,3840,'control'), ...
    'Audit:BeforeEmergenceGate');
assert_error(@() Audit_Checkpoint_State_Gate(3841,3842,3840,'control'), ...
    'Audit:IterationMismatch');

source = fileread(fullfile(run_dir,'Audit_Post_Emergence_Checkpoints.m'));
assert(~contains(source,'7000') && ...
    contains(source,'config.ExpectedIterations') && ...
    contains(source,'Audit_Checkpoint_Preflight') && ...
    contains(source,'Audit_Checkpoint_State_Gate') && ...
    ~contains(source,'run10_checkpoint=consumed once') && ...
    ~contains(source,'MPa\\n'), ...
    'The reusable audit retained stale one-use assumptions.');

fprintf('PASS: checkpoint audit preflight and emergence boundary.\n');
end

function assert_error(action,identifier)
try
    action();
catch ME
    assert(strcmp(ME.identifier,identifier), ...
        'Expected %s but received %s.',identifier,ME.identifier);
    return
end
error('Expected error %s was not raised.',identifier);
end

function write_text(file_name,text)
fid = fopen(file_name,'wb');
assert(fid > 0,'Could not create checkpoint audit fixture.');
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
fwrite(fid,uint8(text),'uint8');
end
