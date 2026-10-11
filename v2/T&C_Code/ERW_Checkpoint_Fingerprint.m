%#########################################################################
%   Function ERW_CHECKPOINT_FINGERPRINT                                  #
%#########################################################################
%%% Identity string for a checkpoint, used to decide whether a saved
%%% trajectory may be resumed.
%%%
%%% A checkpoint holds the state the solver produced up to that iteration.
%%% Replaying it under changed chemistry would splice two model versions
%%% into one trajectory and the join would be invisible in the output, so
%%% the fingerprint covers both the run configuration (horizon, column
%%% count) and the source files that decide the chemistry. Files are
%%% identified by content so a same-size edit cannot reuse an incompatible
%%% checkpoint.
%%%
%%% INPUTS
%%%   NN        number of timesteps in the run
%%%   N_col     number of columns
%%%   runtime_identity exact effective runtime chemistry identity
%%%   varargin  full paths of the source files to fingerprint
%%% OUTPUT
%%%   fp        char row vector; compare with strcmp
function fp = ERW_Checkpoint_Fingerprint(NN, N_col, runtime_identity, varargin)

if ~(ischar(runtime_identity) || (isstring(runtime_identity) && isscalar(runtime_identity)))
    error('ERW:checkpointRuntimeIdentity', ...
        'runtime_identity must be a char vector or string scalar.');
end
parts = {sprintf('NN=%d;N_col=%d', NN, N_col), ...
    ['runtime=' char(runtime_identity)]};
for k = 1:numel(varargin)
    f = varargin{k};
    d = dir(f);
    d = d(~[d.isdir]);
    if isempty(d)
        parts{end+1} = sprintf('%s:missing', f); %#ok<AGROW>
    else
        paths = fullfile({d.folder},{d.name});
        [paths,order] = sort(paths);
        d = d(order);
        for j = 1:numel(d)
            parts{end+1} = sprintf('%s:%s',paths{j},file_sha256(paths{j})); %#ok<AGROW>
        end
    end
end
fp = strjoin(parts, '|');

end

function hex = file_sha256(file_path)
fid = fopen(file_path,'rb');
if fid < 0
    error('ERW:checkpointFingerprintRead', ...
        'Cannot read checkpoint dependency %s.',file_path);
end
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
bytes = fread(fid,Inf,'*uint8');
digest = javaMethod('getInstance','java.security.MessageDigest','SHA-256');
digest.update(typecast(bytes,'int8'));
raw = typecast(digest.digest(),'uint8');
hex = lower(reshape(dec2hex(raw,2).',1,[]));
end
