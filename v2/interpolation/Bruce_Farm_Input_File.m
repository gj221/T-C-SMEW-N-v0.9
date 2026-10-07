function input_file = Bruce_Farm_Input_File(relative_path,repository_path)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Resolve a Bruce Farm input, preferring the copy tracked in this repository.
%%
%% Shapefiles are mirrored under v2/Spatial inputs so a run is reproducible
%% from a clone alone. The sample workbooks are too large to track, so they
%% are read from the survey folder in OneDrive. Set BRUCE_FARM_ROOT to point
%% at that folder somewhere else.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2; repository_path = ''; end
script_directory = fileparts(mfilename('fullpath'));
v2_directory = fileparts(script_directory);

if ~isempty(repository_path)
    tracked_file = fullfile(v2_directory,repository_path);
    if isfile(tracked_file)
        input_file = tracked_file;
        return
    end
end

survey_root = getenv('BRUCE_FARM_ROOT');
if isempty(survey_root)
    survey_root = fullfile(getenv('USERPROFILE'), ...
        'OneDrive - Yale University','Documents','Bruce Farm_Virginia');
end
input_file = fullfile(survey_root,relative_path);
if ~isfile(input_file)
    error('Bruce_Farm_Input_File:MissingInput', ...
        ['%s was not found. Set BRUCE_FARM_ROOT to the folder holding ' ...
         'the Bruce Farm survey data.'],input_file);
end

end
