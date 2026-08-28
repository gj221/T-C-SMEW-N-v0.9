function hex = Bouchet_Fingerprint_Digest(fingerprint)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Short, comparison-safe identity for a campaign validation fingerprint.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The full fingerprint names and hashes every dependency and runs to tens
%% of kilobytes: 49595 characters for this campaign. MATLAB's setenv
%% rejects anything over 32766, so the stamp cannot travel to the solver as
%% environment text.
%%
%% The campaign gate only ever tests the stamp for equality, and a SHA-256
%% digest is exactly as strict for that. The full text stays on disk in
%% validation_stamp.txt, so a mismatch can still be diagnosed by diffing
%% the stamp against a freshly computed fingerprint.
%%
%% INPUT
%%   fingerprint  full validation fingerprint text
%% OUTPUT
%%   hex          lowercase SHA-256 hex digest, 64 characters

if ~(ischar(fingerprint) || ...
        (isstring(fingerprint) && isscalar(fingerprint))) || ...
        strlength(string(fingerprint)) == 0
    error('Bouchet:FingerprintDigestInput', ...
        'The fingerprint must be nonempty text.');
end
bytes = unicode2native(char(fingerprint),'UTF-8');
digest = javaMethod('getInstance','java.security.MessageDigest','SHA-256');
digest.update(typecast(bytes(:),'int8'));
raw = typecast(digest.digest(),'uint8');
hex = lower(reshape(dec2hex(raw,2).',1,[]));

end
