function Bouchet_Write_Validation_Stamp(stamp_file,fingerprint)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Atomically replace a Bouchet validation stamp.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~(ischar(stamp_file) || (isstring(stamp_file) && isscalar(stamp_file))) || ...
        ~(ischar(fingerprint) || (isstring(fingerprint) && isscalar(fingerprint))) || ...
        strlength(string(fingerprint)) == 0
    error('Bouchet:ValidationStampInput', ...
        'The stamp path and fingerprint must be nonempty text.');
end
stamp_file = char(stamp_file);
stamp_dir = fileparts(stamp_file);
if ~exist(stamp_dir,'dir'); mkdir(stamp_dir); end
tmp_file = [tempname(stamp_dir) '.txt'];
fid = fopen(tmp_file,'wt');
if fid < 0
    error('Bouchet:ValidationStampWrite', ...
        'Cannot open temporary validation stamp %s.',tmp_file);
end
try
    fprintf(fid,'%s\n',char(fingerprint));
    fclose(fid);
    fid = -1;
    movefile(tmp_file,stamp_file,'f');
catch ME
    if fid >= 0; fclose(fid); end
    if isfile(tmp_file); delete(tmp_file); end
    rethrow(ME)
end
end
