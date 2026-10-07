function [state,result_file] = Uncertainty_Result_State(results_dir,run_id,C,K,campaign,digest)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% What a result directory holds for one run: 'none', 'complete' (one
%% result, the expected one; see Uncertainty_Result_Matches) or 'stale'
%% (one result of anything else, or without provenance). More than one
%% final MAT is an error, since nothing says which of them to trust.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

state = 'none';
result_file = '';
existing = dir(fullfile(results_dir,sprintf('TC_2D_run%02d_*.mat',run_id)));
existing = existing(~[existing.isdir]);
if isempty(existing); return; end
if numel(existing) > 1
    error('Kantola_Uncertainty:AmbiguousResult', ...
        'Run %d has %d final MAT files in %s; keep exactly one.', ...
        run_id,numel(existing),results_dir);
end
result_file = fullfile(existing.folder,existing.name);
state = 'stale';
contents = whos('-file',result_file);
if ~any(strcmp({contents.name},'run_provenance')); return; end
S = load(result_file,'run_provenance');
if Uncertainty_Result_Matches(S.run_provenance,existing.name,C,K,campaign,digest)
    state = 'complete';
end
end
