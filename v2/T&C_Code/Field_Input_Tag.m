function tag = Field_Input_Tag(agronomic_file,state_file)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Identify the exact per-cell input and spin-up state file contents.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

digest = javaMethod('getInstance','java.security.MessageDigest','SHA-256');
files = {agronomic_file,state_file};
for k = 1:numel(files)
    fid = fopen(files{k},'rb');
    if fid < 0
        error('Field_Input_Tag:FileRead', ...
            'Cannot read cell input dependency: %s',files{k});
    end
    cleanup = onCleanup(@() fclose(fid));
    bytes = fread(fid,Inf,'*uint8');
    digest.update(typecast(bytes,'int8'));
    clear cleanup
end
raw = typecast(digest.digest(),'uint8');
hex = lower(reshape(dec2hex(raw,2).',1,[]));
tag = ['cellsoil_' hex(1:12)];
end
