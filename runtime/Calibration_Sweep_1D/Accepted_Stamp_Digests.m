function [accepted,lineage,rejected] = Accepted_Stamp_Digests( ...
    stored_stamp,logs_root,repo_root)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Validation-stamp digests whose results may be scored under the current
%% stamp.
%%
%% The current digest is always accepted. Each digest named in
%% TC_ACCEPTED_RESULT_DIGESTS (separated by ':' or ',') is accepted only when its
%% archived stamp <logs_root>/validation_stamps/<digest>.txt digests to that
%% name, shares the current stamp's horizon and numerical runtime identity,
%% and differs from the current stamp only in the scoring-side files named
%% by allowed_relative_files. None of those files can change a trajectory,
%% so a result produced under such a prior stamp comes from the same
%% validated model state as the code that scores it now.
%%
%% A requested digest that fails those checks is not an error here: it is
%% returned in rejected with the reason, and the summary refuses a result
%% only when the result actually carries a digest that is not accepted.
%%
%% OUTPUTS
%%   accepted  cell row of digests: the current one first
%%   lineage   struct array, one per accepted prior digest, with fields
%%             digest and changed (cell of the scoring-side files that differ)
%%   rejected  struct array, one per refused digest, with fields digest,
%%             identifier and reason
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stored_stamp = strtrim(char(stored_stamp));
current_digest = Fingerprint_Digest(stored_stamp);
accepted = {current_digest};
lineage = struct('digest',{},'changed',{});
rejected = struct('digest',{},'identifier',{},'reason',{});

requested = strtrim(getenv('TC_ACCEPTED_RESULT_DIGESTS'));
if isempty(requested)
    return
end
requested = strtrim(strsplit(requested,{':',','}));
requested = requested(~cellfun(@isempty,requested));

[current_header,current_map] = parse_stamp(stored_stamp,'current');
allowed = allowed_absolute_files(repo_root);

for k = 1:numel(requested)
    digest = lower(requested{k});
    if isempty(regexp(digest,'^[0-9a-f]{64}$','once'))
        error('Bouchet_Joint_Sweep:LineageDigest', ...
            'TC_ACCEPTED_RESULT_DIGESTS entry %s is not a SHA-256 digest.', ...
            requested{k});
    end
    if strcmp(digest,current_digest) || any(strcmp(digest,accepted)) || ...
            any(strcmp(digest,{rejected.digest}))
        continue
    end
    archive_file = fullfile(logs_root,'validation_stamps',[digest '.txt']);
    if ~isfile(archive_file)
        rejected(end+1) = refusal(digest, ...
            'Bouchet_Joint_Sweep:LineageArchiveMissing', ...
            sprintf('the archived validation stamp %s does not exist', ...
            archive_file)); %#ok<AGROW>
        continue
    end
    prior_stamp = strtrim(fileread(archive_file));
    if ~strcmp(Fingerprint_Digest(prior_stamp),digest)
        rejected(end+1) = refusal(digest, ...
            'Bouchet_Joint_Sweep:LineageArchiveMismatch', ...
            sprintf('the archived validation stamp %s does not digest to its name', ...
            archive_file)); %#ok<AGROW>
        continue
    end
    try
        [prior_header,prior_map] = parse_stamp(prior_stamp,archive_file);
    catch ME
        rejected(end+1) = refusal(digest,ME.identifier,ME.message); %#ok<AGROW>
        continue
    end
    if ~isequal(prior_header,current_header)
        rejected(end+1) = refusal(digest, ...
            'Bouchet_Joint_Sweep:LineageRuntime', ...
            'its stamp has a different horizon or numerical runtime identity than the current stamp'); %#ok<AGROW>
        continue
    end
    paths = union(keys(current_map),keys(prior_map));
    changed = {};
    blocked = '';
    for j = 1:numel(paths)
        path = paths{j};
        if isKey(current_map,path) && isKey(prior_map,path) && ...
                strcmp(current_map(path),prior_map(path))
            continue
        end
        if ~any(strcmp(path,allowed))
            blocked = path;
            break
        end
        changed{end+1} = path; %#ok<AGROW>
    end
    if ~isempty(blocked)
        rejected(end+1) = refusal(digest, ...
            'Bouchet_Joint_Sweep:LineageRejected', ...
            sprintf(['its stamp differs from the current stamp in %s, ' ...
            'which is not a scoring-side file'],blocked)); %#ok<AGROW>
        continue
    end
    accepted{end+1} = digest; %#ok<AGROW>
    lineage(end+1) = struct('digest',digest,'changed',{changed}); %#ok<AGROW>
end

end

function r = refusal(digest,identifier,reason)
r = struct('digest',digest,'identifier',identifier,'reason',reason);
end

function rel = allowed_relative_files()
% Files that score, validate, submit or document results but never run
% inside a trajectory; the campaign fingerprint on Bouchet also hashes the
% controller script. Keep this list identical to rescore_allowed_files in
% submit_sweep.sh; the configuration test compares them.
rel = { ...
    'v2/T&C_Code/Tests/Kantola_N_Magnitude_Audit.m', ...
    'v2/T&C_Code/Tests/test_kantola_n_magnitude_audit.m', ...
    'v2/T&C_Code/Tests/test_bouchet_1d_joint_sweep_configuration.m', ...
    'runtime/Calibration_Sweep_1D/Summarize_Sweep.m', ...
    'runtime/Calibration_Sweep_1D/Validate_Sweep.m', ...
    'runtime/Calibration_Sweep_1D/Assert_Sweep_Provenance.m', ...
    'runtime/Calibration_Sweep_1D/Accepted_Stamp_Digests.m', ...
    'runtime/Calibration_Sweep_1D/README.md', ...
    'runtime/Calibration_Sweep_1D/submit_sweep.sh'};
end

function absolute = allowed_absolute_files(repo_root)
rel = allowed_relative_files();
absolute = cell(size(rel));
for k = 1:numel(rel)
    segments = strsplit(rel{k},'/');
    absolute{k} = fullfile(repo_root,segments{:});
end
end

function [header,map] = parse_stamp(text,label)
% A stamp is 'NN=...;N_col=...|runtime=...|<path>:<sha256>|...' as written
% by ERW_Checkpoint_Fingerprint. The two leading parts describe the horizon
% and the numerical runtime; every other part names one fingerprinted file.
parts = strsplit(text,'|');
header = {};
map = containers.Map('KeyType','char','ValueType','char');
for k = 1:numel(parts)
    part = parts{k};
    if startsWith(part,'NN=') || startsWith(part,'runtime=')
        header{end+1} = part; %#ok<AGROW>
        continue
    end
    tokens = regexp(part,'^(.+):([0-9a-f]{64})$','tokens','once');
    if isempty(tokens)
        error('Bouchet_Joint_Sweep:LineageUnreadable', ...
            'Validation stamp %s has an entry without a file digest: %s', ...
            label,part);
    end
    if isKey(map,tokens{1})
        error('Bouchet_Joint_Sweep:LineageUnreadable', ...
            'Validation stamp %s lists %s twice.',label,tokens{1});
    end
    map(tokens{1}) = tokens{2};
end
if numel(header) ~= 2
    error('Bouchet_Joint_Sweep:LineageUnreadable', ...
        'Validation stamp %s does not start with the horizon and runtime parts.', ...
        label);
end
end
