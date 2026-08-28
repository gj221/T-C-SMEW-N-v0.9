function test_bouchet_linux_portability(runtime_dirs)
% Deployed MATLAB path strings must be portable to Linux.
%
% Called with no argument the test scans every deployed campaign package.
% A campaign gate passes its own directory so that validating one campaign
% cannot be satisfied by the portability of a different one.

test_dir = fileparts(mfilename('fullpath'));
code_dir = fileparts(test_dir);
v2_dir = fileparts(code_dir);
repo_root = fileparts(v2_dir);
run_dir = fullfile(v2_dir,'run');
if nargin < 1 || isempty(runtime_dirs)
    runtime_dirs = { ...
        fullfile(repo_root,'runtime', ...
            'Kantola_5yr_1D_Runs10to12_ScheduleFix_R2025a_Bouchet'), ...
        fullfile(repo_root,'runtime', ...
            'Kantola_5yr_2D_Runs01to09_DenitPH7p8_Alpha1p3_Bouchet')};
elseif ~iscell(runtime_dirs)
    runtime_dirs = {char(runtime_dirs)};
end

files = [dir(fullfile(code_dir,'*.m')); ...
    dir(fullfile(code_dir,'Parameters','*.m')); ...
    dir(fullfile(run_dir,'*.m'))];
for k = 1:numel(runtime_dirs)
    package = char(runtime_dirs{k});
    assert(exist(package,'dir') == 7, ...
        'Runtime package %s does not exist.',package);
    files = [files; dir(fullfile(package,'*.m'))]; %#ok<AGROW>
end
pattern = "(?:""[^""\r\n]*""|'(?:[^']|'')*')";
failures = strings(0,1);
for k = 1:numel(files)
    file = fullfile(files(k).folder,files(k).name);
    literals = regexp(fileread(file),pattern,'match');
    for j = 1:numel(literals)
        %% Escapes that are not path separators: the printf set and
        %% the regexp whitespace class, plus the one LaTeX label.
        literals{j} = strrep(literals{j},[char(92) 'Psi'],'');
        literal = literals{j};
        for escape = {'\n','\r','\t','\s'}
            literal = strrep(literal,escape{1},'');
        end
        if contains(literal,'\')
            failures(end+1,1) = string(file) + ": " + string(literals{j}); %#ok<AGROW>
        end
    end
end
assert(isempty(failures), ...
    'Windows-only path string(s) remain:\n%s',strjoin(failures,newline));

%% A deployed wrapper that does not parse fails only after the job has
%% waited in the queue and started MATLAB, so check it here instead.
syntax_ids = {'STRIN','SNTXERR','BADCH','MDOTM','NOPTS','ENDCT'};
syntax_failures = strings(0,1);
for k = 1:numel(files)
    file = fullfile(files(k).folder,files(k).name);
    messages = checkcode(file,'-id','-struct');
    for j = 1:numel(messages)
        if any(strcmp(messages(j).id,syntax_ids))
            syntax_failures(end+1,1) = ...
                string(file) + ":" + messages(j).line + ": " + ...
                string(messages(j).message); %#ok<AGROW>
        end
    end
end
assert(isempty(syntax_failures), ...
    'Deployed MATLAB source does not parse:\n%s', ...
    strjoin(syntax_failures,newline));

equivalence_source = fileread(fullfile(code_dir,'ERW_Equivalence_Fingerprint.m'));
assert(contains(equivalence_source, ...
       "fullfile(code_dir,'Inputs','ERW_MOD_PARAM_parameters.xlsx')") && ...
       contains(equivalence_source, ...
       "fullfile(code_dir,'Inputs','blueridge_CDF_data.json')") && ...
       ~contains(equivalence_source, ...
       "fullfile(code_dir,'ERW_MOD_PARAM_parameters.xlsx')"), ...
    'The equivalence fingerprint must use the deployed workbook and rock data.');
assert(contains(equivalence_source,"version('-release')"), ...
    'The equivalence fingerprint must include the MATLAB release.');

baseline_launcher = fileread(fullfile(run_dir,'ERW_run.m'));
assert(contains(baseline_launcher,'ERW_Equivalence_Fingerprint'), ...
    'ERW_run must fingerprint a generated equivalence baseline.');

addpath(code_dir);
fingerprint = ERW_Equivalence_Fingerprint(v2_dir,10);
assert(contains(fingerprint,'ERW_MOD_PARAM_parameters.xlsx:') && ...
       contains(fingerprint,'blueridge_CDF_data.json:') && ...
       ~contains(fingerprint,':missing'), ...
    'The equivalence fingerprint must hash every deployed dependency.');
if ispc
    alternate_v2_dir = strrep(v2_dir,char(92),'/');
    alternate_fingerprint = ERW_Equivalence_Fingerprint(alternate_v2_dir,10);
    assert(strcmp(fingerprint,alternate_fingerprint), ...
        'Equivalent Windows path separators must give the same fingerprint.');
end

fprintf('PASS: deployed MATLAB paths and baseline dependencies are Linux-safe.\n');
end
