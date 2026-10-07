function inventory = Audit_Checkpoint_Preflight(files,expected_hashes, ...
    expected_iterations,minimum_iteration)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Validate immutable checkpoint candidates before any state is loaded.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isstring(files); files = cellstr(files); end
if isstring(expected_hashes); expected_hashes = cellstr(expected_hashes); end
files = files(:);
expected_hashes = expected_hashes(:);
expected_iterations = double(expected_iterations(:));
n = numel(files);

if n < 1 || numel(expected_hashes) ~= n || ...
        numel(expected_iterations) ~= n
    error('Audit:PreflightShape', ...
        'Files, hashes and expected iterations must have equal nonzero lengths.');
end
if ~(isscalar(minimum_iteration) && isfinite(minimum_iteration) && ...
        minimum_iteration == fix(minimum_iteration) && minimum_iteration >= 0)
    error('Audit:MinimumIteration', ...
        'The minimum iteration must be a nonnegative integer.');
end
if any(~isfinite(expected_iterations)) || ...
        any(expected_iterations ~= fix(expected_iterations)) || ...
        any(expected_iterations <= minimum_iteration)
    error('Audit:BeforeEmergenceConfig', ...
        'Every recorded checkpoint iteration must be strictly above %d.', ...
        minimum_iteration);
end

bytes = zeros(n,1);
modified = NaT(n,1);
actual_hashes = strings(n,1);
for k = 1:n
    file_path = files{k};
    expected_hash = expected_hashes{k};
    if ~(ischar(file_path) && isrow(file_path) && ~isempty(file_path))
        error('Audit:CheckpointPath','Checkpoint paths must be nonempty text.');
    end
    if ~(ischar(expected_hash) && ...
            ~isempty(regexp(expected_hash,'^[0-9A-Fa-f]{64}$','once')))
        error('Audit:ExpectedHash','Expected SHA-256 values must contain 64 hex digits.');
    end
    before = dir(file_path);
    if numel(before) ~= 1 || before.isdir
        error('Audit:CheckpointMissing','Missing checkpoint %s.',file_path);
    end
    actual_hash = file_sha256(file_path);
    after = dir(file_path);
    if numel(after) ~= 1 || after.bytes ~= before.bytes || ...
            after.datenum ~= before.datenum
        error('Audit:CheckpointChanged', ...
            'Checkpoint metadata changed while hashing %s.',file_path);
    end
    if ~strcmpi(actual_hash,expected_hash)
        error('Audit:HashMismatch', ...
            'Checkpoint SHA-256 does not match the recorded value for %s.',file_path);
    end
    bytes(k) = before.bytes;
    modified(k) = datetime(before.datenum,'ConvertFrom','datenum');
    actual_hashes(k) = upper(string(actual_hash));
end

inventory = table(string(files),expected_iterations,bytes,modified, ...
    actual_hashes,'VariableNames',{'Checkpoint','ExpectedIteration', ...
    'Bytes','Modified','SHA256'});
end

function hex = file_sha256(file_path)
fid = fopen(file_path,'rb');
if fid < 0
    error('Audit:CheckpointRead','Cannot read checkpoint %s.',file_path);
end
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
digest = javaMethod('getInstance','java.security.MessageDigest','SHA-256');
while true
    bytes = fread(fid,1024*1024,'*uint8');
    if isempty(bytes); break; end
    digest.update(typecast(bytes,'int8'));
end
raw = typecast(digest.digest(),'uint8');
hex = upper(reshape(dec2hex(raw,2).',1,[]));
end
