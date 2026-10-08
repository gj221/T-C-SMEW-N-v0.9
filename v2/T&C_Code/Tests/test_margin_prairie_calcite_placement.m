function test_margin_prairie_calcite_placement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Native prairie and calcite stay in the margin cells for run 9.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_dir = fileparts(mfilename('fullpath'));
v2 = fileparts(fileparts(test_dir));
addpath(fullfile(v2,'run'),fullfile(v2,'T&C_Code'));

result_dir = tempname;
mkdir(result_dir);
old_checkpoint_minutes = getenv('TC_CKPT_MINUTES');
cleanup = onCleanup(@() restore_test_environment( ...
    result_dir,old_checkpoint_minutes)); %#ok<NASGU>
setenv('TC_CKPT_MINUTES','0');
%% These checks read the run's full hourly and ten-minute records, which
%% a multi-square run keeps only when asked (package y).
old_records = getenv('TC_FULL_RECORDS');
records_cleanup = onCleanup(@() setenv('TC_FULL_RECORDS',old_records)); %#ok<NASGU>
setenv('TC_FULL_RECORDS','1');

result_file = Run_One_Subfield(9,'StartDate','2016-11-09','NDays',3, ...
    'AppeloAlpha',1.3,'ParallelWorkers',6,'ResultDirectory',result_dir);
S = load(result_file,'R','topo','ColCfg','PAR', ...
    'ERW_PARAM_input','ERW_PARAM_control','ERW');

is_main = strcmp(S.topo.cell_region,'main');
is_margin = ~is_main;
assert(sum(is_main) == 4 && sum(is_margin) == 8, ...
    'The 2 x 2 field grid must contain four center and eight margin cells.');
assert(abs(sum(S.topo.cell_area_m2(is_main))-2500) < 1e-10 && ...
       abs(sum(S.topo.cell_area_m2(is_margin))-1100) < 1e-10, ...
    'The field-center and margin areas are incorrect.');
assert(strcmp(S.R.treatment,'calcite') && ...
       strcmp(S.R.layout,'margin_prairie') && ...
       ~S.R.apply_main && S.R.apply_margin, ...
    'Run 9 is not the calcite-treated native-prairie margin scenario.');

apply_rock = [S.ColCfg.apply_rock].';
pft = string({S.ColCfg.PFT_OVERRIDE}).';
assert(~any(apply_rock(is_main)) && all(apply_rock(is_margin)), ...
    'Calcite must be applied to the margin and not the field center.');
assert(all(pft(is_main) == "") && ...
       all(pft(is_margin) == "PFT_Prairie_Margin"), ...
    'Native prairie must replace the crop only in the margin.');

for c = find(is_main).'
    assert(isequaln(S.ColCfg(c).ERW_PARAM,S.ERW_PARAM_control), ...
        'Field-center cell %d did not retain control chemistry.',c);
    assert(sum(S.ERW{c}.TC_const.M_rock_add) == 0, ...
        'Field-center cell %d received calcite.',c);
end
for c = find(is_margin).'
    assert(isequaln(S.ColCfg(c).ERW_PARAM,S.ERW_PARAM_input), ...
        'Margin cell %d did not receive calcite chemistry.',c);
    assert(isequal(S.PAR{c}.Ccrown,[1 0]) && ...
           S.PAR{c}.ZR95_H(1) == 2500 && S.PAR{c}.ZR95_L(1) == 1200, ...
        'Margin cell %d did not load the native prairie PFT.',c);
    application_row = ERW_Application_Indices(1,1/144,3*144);
    assert(application_row == 145 && ...
        all(S.ERW{c}.TC_const.M_rock_add(1:application_row-1) == 0) && ...
        S.ERW{c}.TC_const.M_rock_add(application_row) == 5000, ...
        'Margin cell %d did not apply calcite at the midnight boundary.',c);
    assert(all(S.ERW{c}.bg.M_min(1:application_row-1,:) == 0,'all') && ...
        any(S.ERW{c}.bg.M_min(application_row,:) > 0), ...
        'Margin cell %d has primary mineral on the wrong side of the event.',c);
end

applied_g_m2 = cellfun(@(x) sum(x.TC_const.M_rock_add),S.ERW);
assert(all(applied_g_m2(is_margin) == 5000), ...
    'Each margin cell must receive the prescribed 5000 g m^-2 event.');
applied_mass_g = sum(applied_g_m2(:).*S.topo.cell_area_m2(:));
assert(abs(applied_mass_g-5000*1100) < 1e-8*applied_mass_g, ...
    'The margin calcite mass does not close to rate times treated area.');

fprintf(['PASS: run 9 keeps calcite and native prairie in the 1100 m2 ' ...
    'margin; the 2500 m2 field center remains unamended crop.\n']);
end

function restore_test_environment(result_dir,checkpoint_minutes)
setenv('TC_CKPT_MINUTES',checkpoint_minutes);
pool = gcp('nocreate');
if ~isempty(pool); delete(pool); end
if isfolder(result_dir); rmdir(result_dir,'s'); end
end
