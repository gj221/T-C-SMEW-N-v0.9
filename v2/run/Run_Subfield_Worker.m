function Run_Subfield_Worker(run_id, status_file, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Child-process entry point used by the parallel scenario launchers.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

write_status(status_file,'running','','');
try
    result_file = Run_One_Subfield(run_id,varargin{:});
    write_status(status_file,'complete',result_file,'');
catch ME
    write_status(status_file,'failed',ME.identifier,ME.message);
    rethrow(ME)
end

end

function write_status(status_file, state, detail, message)
[status_dir,~,~] = fileparts(status_file);
if ~isempty(status_dir) && ~exist(status_dir,'dir'); mkdir(status_dir); end
tmp_file = [status_file '.tmp'];
fid = fopen(tmp_file,'w');
if fid < 0
    error('Run_Subfield_Worker:StatusWrite','Cannot write %s.',tmp_file);
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid,'state=%s\n',state);
fprintf(fid,'detail=%s\n',strrep(char(detail),newline,' '));
fprintf(fid,'message=%s\n',strrep(char(message),newline,' '));
fprintf(fid,'finished=%s\n',datestr(now,30));
clear cleanup
movefile(tmp_file,status_file,'f');
end
