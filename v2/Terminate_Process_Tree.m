function terminated = Terminate_Process_Tree(process)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Stop a launcher process and its descendants.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

terminated = true;
if isempty(process); return; end

try
    if ispc && isa(process,'System.Diagnostics.Process')
        if process.HasExited; return; end
        [~,~] = system(sprintf( ...
            'taskkill /PID %d /T /F >NUL 2>&1',process.Id));
        process.WaitForExit(5000);
        terminated = process.HasExited;
        if ~terminated
            process.Kill();
            process.WaitForExit(5000);
            terminated = process.HasExited;
        end
    elseif process.isAlive()
        process.destroyForcibly();
        terminated = ~process.isAlive();
    end
catch ME
    terminated = false;
    warning('Terminate_Process_Tree:Failed','%s',ME.message);
end

if ~terminated
    warning('Terminate_Process_Tree:StillRunning', ...
        'A launched worker process could not be stopped.');
end

end
