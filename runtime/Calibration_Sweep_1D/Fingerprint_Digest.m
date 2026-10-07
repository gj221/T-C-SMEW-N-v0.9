function hex = Fingerprint_Digest(fingerprint)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Short, comparison-safe identity for a campaign validation fingerprint.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
