function test_ssp_controller_files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Every file the SSP controller requires before it submits must exist in
%% the checkout. Its pin list once pointed at runtime/SSP_1D, where the
%% list has never been, and the controller refused every campaign.
%% Entries named after the scenario (the forcing, deployed apart from the
%% package) are the controller's own runtime check, not the package's.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(fileparts(test_dir)));
script = fileread(fullfile(repo_root,'runtime','SSP_1D','submit_ssp.sh'));

%% Top-level names: the literal ones, then those built from them.
assigned = [regexp(script,'(?m)^(\w+)=''([^'']*)''\s*$','tokens') ...
    regexp(script,'(?m)^(\w+)="([^"]*)"\s*$','tokens')];
block = regexp(script,'required_files=\((.*?)\n\)','tokens','once');
assert(~isempty(assigned) && ~isempty(block), ...
    'The SSP controller no longer declares its required files.');
entries = regexp(block{1},'["'']([^"'']+)["'']','tokens');
entries = cellfun(@(e) e{1},entries,'UniformOutput',false);
assert(any(strcmp(entries,'$pin_list')),'The pin list must be a required file.');
n_checked = 0;
for k = 1:numel(entries)
    path = entries{k};
    for pass = 1:2
        for a = 1:numel(assigned)
            path = strrep(path,['$' assigned{a}{1}],assigned{a}{2});
        end
    end
    if contains(path,'${TC_SSP}'); continue; end
    assert(~contains(path,'$'),'The required file %s did not resolve.',entries{k});
    assert(isfile(fullfile(repo_root,path)), ...
        'The SSP controller requires %s, which the checkout does not hold.',path);
    n_checked = n_checked+1;
end
assert(n_checked == numel(entries)-1, ...
    'Every required file but the scenario forcing must be checked.');
%% The controllers check their pin lists strictly, so a file a list names
%% and the package no longer holds refuses every submission.
for list = {'shared_pins.sha256','ssp_pins.sha256'}
    pins = regexp(fileread(fullfile(repo_root,'runtime','Base_1D',list{1})), ...
        '(?m)^[0-9a-f]{64} \*([^\r\n]+)','tokens');
    paths = cellfun(@(t) t{1},pins,'UniformOutput',false);
    missing = paths(~cellfun(@(q) isfile(fullfile(repo_root,q)),paths));
    assert(numel(paths) > 40 && isempty(missing), ...
        '%s names files the checkout does not hold: %s',list{1},strjoin(missing,', '));
end
fprintf('PASS: every file the SSP controller requires or pins is in the checkout.\n');
end
