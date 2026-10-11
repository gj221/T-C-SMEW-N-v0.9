function test_biogeochemistry_io_path_portability
% Deposition inputs must resolve independently of the MATLAB working folder.

test_dir = fileparts(mfilename('fullpath'));
code_dir = fileparts(test_dir);
input_dir = fullfile(code_dir,'Inputs');
addpath(code_dir);

location = load(fullfile(input_dir,'Data_Kantola.mat'),'Lat','Lon');
temporary_dir = tempname;
mkdir(temporary_dir);
old_dir = pwd;
cleanup = onCleanup(@() restore_folder(old_dir,temporary_dir));
cd(temporary_dir);

B = Biogeochemistry_IO([0;300],location.Lat,location.Lon,0.01,0);
fields = {'DepN','DepP','DepK','DepCa','DepMg','DepSi'};
for k = 1:numel(fields)
    value = B.(fields{k});
    assert(~isempty(value) && all(isfinite(value(:))), ...
        '%s must be finite after loading the deposition input.',fields{k});
end

fprintf('PASS: deposition inputs resolve outside the code directory.\n');
end

function restore_folder(old_dir,temporary_dir)
cd(old_dir);
rmdir(temporary_dir,'s');
end
