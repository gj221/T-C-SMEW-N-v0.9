function Validate_Bouchet_1D_Preapplication_Identity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Same-node discriminator for pre-application treatment state.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

campaign_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(campaign_dir));
v2_dir = fullfile(repo_root,'v2');
addpath(fullfile(v2_dir,'run'),fullfile(v2_dir,'T&C_Code'));

root = tempname(tempdir);
mkdir(root);
old_ckpt = getenv('TC_CKPT_MINUTES');
old_force = getenv('TC_CKPT_FORCE_RESUME');
old_workers = getenv('TC_PARALLEL_WORKERS');
cleanup = onCleanup(@() restore_environment( ...
    root,old_ckpt,old_force,old_workers));
setenv('TC_CKPT_MINUTES','0');
setenv('TC_CKPT_FORCE_RESUME','');
setenv('TC_PARALLEL_WORKERS','0');

run_ids = [10 11 12 12];
labels = {'control','basalt','calcite_a','calcite_b'};
variables = {'Date','Pr','Ta','Ws','ea','Pre','PARB','PARD', ...
    'P','O','V','Ts','Tdp','Lk','T_H','T_L','EG','LAI_L','B_L', ...
    'PH_d','VOL','NITflx','DENITflx','N2Oflx_nit','N2flx','N2Oflx', ...
    'LEAK_NO3','LEAK_NH4','NavlI','ERW'};
validation_workers = 4;
if isempty(ver('parallel')) || ~license('test','Distrib_Computing_Toolbox')
    error('Preapplication:ParallelToolboxUnavailable', ...
        'Parallel Computing Toolbox is required for concurrent validation.');
end
pool = parpool('Processes',validation_workers);
pctRunOnAll maxNumCompThreads(1)
result_files = cell(size(run_ids));
parfor (k = 1:numel(run_ids),validation_workers)
    result_dir = fullfile(root,labels{k});
    mkdir(result_dir);
    result_files{k} = Run_One_Subfield(run_ids(k), ...
        'StartDate',datenum(2016,1,1),'NDays',315, ...
        'AppeloAlpha',1.3,'ParallelWorkers',0, ...
        'SaveResults',true,'ResultDirectory',result_dir);
end
delete(pool);
clear pool

states = cell(size(run_ids));
for k = 1:numel(run_ids)
    states{k} = load(result_files{k},variables{:});
end

forcing = {'Date','Pr','Ta','Ws','ea','Pre','PARB','PARD'};
for k = 2:3
    for q = 1:numel(forcing)
        name = forcing{q};
        assert(isequaln(states{1}.(name),states{k}.(name)), ...
            'Preapplication:ForcingMismatch', ...
            '%s differs between control and run %d.',name,run_ids(k));
    end
end

first_amendment = datenum(2016,11,10);
hour_mask = states{1}.Date(:) < first_amendment;
daily_dates = states{1}.Date(1) + (0:size(states{1}.P,1)-1)';
daily_mask = daily_dates < first_amendment;
state_names = {'P','O','V','Ts','Tdp','Lk','T_H','T_L','EG', ...
    'LAI_L','B_L','PH_d','VOL','NITflx','DENITflx','N2Oflx_nit', ...
    'N2flx','N2Oflx','LEAK_NO3','LEAK_NH4','NavlI'};
for k = 2:3
    for q = 1:numel(state_names)
        name = state_names{q};
        a = before_amendment(states{1}.(name),hour_mask,daily_mask,name);
        b = before_amendment(states{k}.(name),hour_mask,daily_mask,name);
        assert(isequaln(a,b),'Preapplication:StateMismatch', ...
            '%s differs before amendment between control and run %d.', ...
            name,run_ids(k));
    end
end

application_index = ERW_Application_Indices(314,1/144,315*144);
assert(application_index == 45217, ...
    'Preapplication:ApplicationIndex', ...
    'The 10-Nov-2016 amendment must occupy ERW row 45217.');
erw_fields = {'pH','Ca','Mg','K','Na','An','Al','Si','PO4', ...
    'NO3','NH4','H','DIC','CaCO3','MgCO3','pH_bulk'};
for k = 2:3
    for q = 1:numel(erw_fields)
        name = erw_fields{q};
        a = states{1}.ERW{1}.bg.(name);
        b = states{k}.ERW{1}.bg.(name);
        a = a(1:application_index-1,:);
        b = b(1:application_index-1,:);
        assert(isequaln(a,b),'Preapplication:ERWStateMismatch', ...
            'ERW %s differs before amendment for run %d.',name,run_ids(k));
    end
end

control_schedule = states{1}.ERW{1}.TC_const.M_rock_add;
assert(all(control_schedule == 0), ...
    'Preapplication:ControlAmendment', ...
    'The control schedule contains rock.');
for k = 2:3
    schedule = states{k}.ERW{1}.TC_const.M_rock_add;
    primary = states{k}.ERW{1}.bg.M_min;
    assert(all(schedule(1:application_index-1) == 0) && ...
        schedule(application_index) == 5000, ...
        'Preapplication:AmendmentBoundary', ...
        'Run %d does not apply 5000 g m^-2 at row 45217.',run_ids(k));
    assert(all(primary(1:application_index-1,:) == 0,'all'), ...
        'Preapplication:PrematureMineral', ...
        'Run %d contains primary mineral before amendment.',run_ids(k));
end

for q = 1:numel(state_names)
    name = state_names{q};
    assert(isequaln(states{3}.(name),states{4}.(name)), ...
        'Preapplication:RepeatMismatch', ...
        'Repeated calcite runs differ in %s.',name);
end
for q = 1:numel(erw_fields)
    name = erw_fields{q};
    assert(isequaln(states{3}.ERW{1}.bg.(name), ...
        states{4}.ERW{1}.bg.(name)), ...
        'Preapplication:RepeatERWMismatch', ...
        'Repeated calcite runs differ in ERW %s.',name);
end

fprintf(['PASS: control, basalt and calcite are exact before amendment; ' ...
    'the repeated calcite trajectory is exact on one node.\n']);
end

function selected = before_amendment(value,hour_mask,daily_mask,name)
if size(value,1) == numel(hour_mask)
    mask = hour_mask;
elseif size(value,1) == numel(daily_mask)
    mask = daily_mask;
else
    error('Preapplication:StateLength', ...
        '%s has an unsupported leading dimension.',name);
end
subs = repmat({':'},1,ndims(value));
subs{1} = mask;
selected = value(subs{:});
end

function restore_environment(root,ckpt,force_resume,workers)
setenv('TC_CKPT_MINUTES',ckpt);
setenv('TC_CKPT_FORCE_RESUME',force_resume);
setenv('TC_PARALLEL_WORKERS',workers);
pool = gcp('nocreate');
if ~isempty(pool); delete(pool); end
if isfolder(root) && startsWith(root,tempdir)
    rmdir(root,'s');
end
end
