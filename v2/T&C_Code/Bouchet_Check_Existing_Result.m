function state = Bouchet_Check_Existing_Result( ...
    results_dir,run_id,validation_fingerprint,n_columns,n_days)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Decide whether a final MAT already in the campaign directory may stand.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Matching a filename wildcard proves nothing about which code produced
%% the file. A campaign that is validated, partly run, revalidated after a
%% source change and then resumed would silently mix two model versions:
%% the completed runs stay old while the missing ones use the new code.
%% So the saved provenance is loaded and its validation fingerprint,
%% column count and horizon are compared with the current campaign gate.
%%
%% INPUTS
%%   results_dir            campaign result directory
%%   run_id                 treatment number
%%   validation_fingerprint stamp the current job validated against
%%   n_columns              columns this campaign run must have
%%   n_days                 horizon in days this campaign run must have
%% OUTPUT
%%   state                  'none' if the run still has to be executed,
%%                          'complete' if a matching result already exists

if ~(ischar(validation_fingerprint) || ...
        (isstring(validation_fingerprint) && isscalar(validation_fingerprint))) || ...
        strlength(string(validation_fingerprint)) == 0
    error('Bouchet:ExistingResultFingerprint', ...
        'The validation fingerprint must be nonempty text.');
end
results_dir = char(results_dir);
validation_fingerprint = char(validation_fingerprint);

existing = dir(fullfile(results_dir,sprintf('TC_2D_run%02d_*.mat',run_id)));
existing = existing(~[existing.isdir]);
if isempty(existing)
    state = 'none';
    return
end
if numel(existing) > 1
    error('Bouchet:AmbiguousResult', ...
        'Run %d has %d final MAT files in %s. Keep exactly one.', ...
        run_id,numel(existing),results_dir);
end
result_file = fullfile(existing.folder,existing.name);

try
    contents = whos('-file',result_file);
catch ME
    error('Bouchet:UnreadableResult', ...
        'Run %d has an unreadable final MAT %s: %s', ...
        run_id,result_file,ME.message);
end
if ~any(strcmp({contents.name},'run_provenance'))
    error('Bouchet:UnprovenancedResult', ...
        ['Run %d has a final MAT without run_provenance: %s. It predates ' ...
         'campaign provenance and cannot be accepted.'],run_id,result_file);
end
try
    saved = load(result_file,'run_provenance');
catch ME
    error('Bouchet:UnreadableResult', ...
        'Run %d has an unreadable final MAT %s: %s', ...
        run_id,result_file,ME.message);
end
prov = saved.run_provenance;

if ~isfield(prov,'validation_fingerprint') || ...
        ~strcmp(char(prov.validation_fingerprint),validation_fingerprint)
    error('Bouchet:StaleResult', ...
        ['Run %d was produced under a different validation than the one ' ...
         'this job passed: %s. Move it out of the campaign directory or ' ...
         'rerun the whole campaign under one validation.'], ...
        run_id,result_file);
end
if ~isfield(prov,'n_columns') || prov.n_columns ~= n_columns
    error('Bouchet:WrongResultColumns', ...
        'Run %d has a final MAT with the wrong column count: %s.', ...
        run_id,result_file);
end
if ~isfield(prov,'n_days') || prov.n_days ~= n_days
    error('Bouchet:WrongResultHorizon', ...
        'Run %d has a final MAT that does not cover %d days: %s.', ...
        run_id,n_days,result_file);
end

state = 'complete';
fprintf('Run %d already has a matching final MAT: %s\n',run_id,result_file);

end
