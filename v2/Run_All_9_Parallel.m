function plan = Run_All_9_Parallel(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Launch the original nine scenarios through the twelve-case launcher.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if mod(numel(varargin),2) ~= 0
    error('Run_All_9_Parallel:BadOptions','Options must be name/value pairs.');
end
names = varargin(1:2:end);
idx = find(cellfun(@(x) strcmp(char(x),'Runs'),names),1);
if isempty(idx)
    varargin = [varargin {'Runs',1:9}];
else
    runs = varargin{2*idx};
    if ~(isnumeric(runs) && isvector(runs) && all(isfinite(runs)) && ...
            all(runs == fix(runs)) && all(ismember(runs,1:9)) && ...
            numel(unique(runs)) == numel(runs))
        error('Run_All_9_Parallel:BadRuns', ...
            'Runs must contain unique integer ids from 1 to 9.');
    end
end

try
    plan = Run_All_12_Parallel(varargin{:});
catch ME
    id = strrep(ME.identifier,'Run_All_12_Parallel','Run_All_9_Parallel');
    if strcmp(id,ME.identifier); rethrow(ME); end
    throwAsCaller(MException(id,'%s',ME.message));
end

end
