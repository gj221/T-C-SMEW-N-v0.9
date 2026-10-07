function test_calcite_case_campaigns
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A calcite-case campaign (TC_CALCITE_CASE, such as 8, the agronomic
%% 5 t/ha rate) runs only the calcite treatments, 1-D run 12 and 2-D runs
%% 3, 6 and 9, under its own tag: on the validation forcing through the
%% Base controller and on a scenario forcing through the SSP one.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(fileparts(test_dir)));
runtime_dir = fullfile(repo_root,'runtime');
one_d_dir = fullfile(runtime_dir,'Base_1D');
old_path = addpath(one_d_dir);
old_case = getenv('TC_CALCITE_CASE');
cleanup = onCleanup(@() restore(old_path,old_case)); %#ok<NASGU>

%% Unset is the standard campaign; a positive integer names a workbook row.
setenv('TC_CALCITE_CASE','');
[calcite_case,suffix] = Calcite_Case;
assert(isempty(calcite_case) && isempty(suffix), ...
    'An unset TC_CALCITE_CASE must select the standard campaign.');
setenv('TC_CALCITE_CASE','8');
[calcite_case,suffix] = Calcite_Case;
assert(isequal(calcite_case,8) && strcmp(suffix,'_case8'), ...
    'TC_CALCITE_CASE=8 must select workbook case 8 under the _case8 tag.');
for bad = {'0','-8','8.5','x'}
    setenv('TC_CALCITE_CASE',bad{1});
    refused = false;
    try
        Calcite_Case;
    catch err
        refused = strcmp(err.identifier,'Bouchet_Kantola_CEC13:BadCalciteCase');
    end
    assert(refused,'TC_CALCITE_CASE=%s must be refused.',bad{1});
end

%% The Base controller gives the case campaign its own tags, job names and
%% trees, and submits only the calcite runs.
assert_tokens(fullfile(one_d_dir,'submit_all_12.sh'),{ ...
    'calcite_case="${TC_CALCITE_CASE:-}"', ...
    'case_suffix="_case${calcite_case}"', ...
    'tag_1d="Base_1D${case_suffix}"', ...
    'tag_2d="Base_2D${case_suffix}"', ...
    'prod_1d_ntasks=1; prod_1d_runs=''12''', ...
    'prod_2d_array=''3,6,9''', ...
    'results_1d="$repo_root/v2/SIM/Subfield_12/$tag_1d"', ...
    'logs_2d="$repo_root/v2/logs/Subfield_12/$tag_2d"', ...
    'TC_CALCITE_CASE=$calcite_case', ...
    '--job-name="Base1D${case_suffix}"', ...
    '--job-name="Base2D${case_suffix}"', ...
    '--ntasks="$prod_1d_ntasks"', ...
    '--array="$prod_2d_array"'});
assert_tokens(fullfile(one_d_dir,'submit_base_1d.sbatch'),{ ...
    'runs=(12)', ...
    'export TC_CAMPAIGN_TAG="Base_1D${case_suffix}"', ...
    'for run_id in "${runs[@]}"; do', ...
    '"${SLURM_NTASKS:-0}" != "${#runs[@]}"', ...
    'results_dir="$repo_root/v2/SIM/Subfield_12/$TC_CAMPAIGN_TAG"'});
assert_tokens(fullfile(runtime_dir,'Base_2D','submit_base_2d.sbatch'),{ ...
    'export TC_CAMPAIGN_TAG="Base_2D${case_suffix}"'});

%% Every task passes the case to the model and refuses a non-calcite run,
%% after the validation branch where there is one: a 2-D validation is
%% array element 1, run 1.
tasks = { ...
    fullfile(one_d_dir,'Run_Base_1D_Task.m'),'run_id ~= 12'; ...
    fullfile(runtime_dir,'Base_2D','Run_Base_2D_Task.m'),'~ismember(run_id,[3 6 9])'; ...
    fullfile(runtime_dir,'SSP_1D','Run_SSP_1D_Task.m'),'run_id ~= 12'; ...
    fullfile(runtime_dir,'SSP_2D','Run_SSP_2D_Task.m'),'~ismember(run_id,[3 6 9])'};
for k = 1:size(tasks,1)
    text = fileread(tasks{k,1});
    [~,name] = fileparts(tasks{k,1});
    guard = strfind(text,['~isempty(calcite_case) && ' tasks{k,2}]);
    assert(isscalar(guard),'%s must refuse a non-calcite run once.',name);
    assert(contains(text,'{''CalciteCase'',calcite_case}'), ...
        '%s must pass the calcite case to the model.',name);
    validation = strfind(text,'if strcmp(validate_only,''1'')');
    assert(isempty(validation) || validation < guard, ...
        '%s must refuse a non-calcite run after its validation branch.',name);
end
for name = {'Run_Base_1D_Task.m','Validate_Base_1D.m'}
    assert(contains(fileread(fullfile(one_d_dir,name{1})), ...
        '[''Base_1D'' case_suffix]'),'%s must tag the campaign with the case.',name{1});
end
two_d_task = fileread(fullfile(runtime_dir,'Base_2D','Run_Base_2D_Task.m'));
assert(contains(two_d_task,'[''Base_2D'' case_suffix]') && ...
    contains(two_d_task,'fullfile(one_d_dir,''Calcite_Case.m'')'), ...
    'The 2-D task must tag the campaign with the case and gate on the helper.');
fprintf('PASS: calcite-case campaigns run the calcite treatments under their own tag.\n');
end

function assert_tokens(file,tokens)
text = fileread(file);
missing = tokens(~cellfun(@(t) contains(text,t),tokens));
[~,name,ext] = fileparts(file);
assert(isempty(missing),'%s omits: %s',[name ext],strjoin(missing,' | '));
end

function restore(old_path,old_case)
path(old_path);
setenv('TC_CALCITE_CASE',old_case);
end
