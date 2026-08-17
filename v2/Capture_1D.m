%% Capture ERW_function inputs from the single-column frame.
%% Arm with: TC_ERW_CAPTURE=1 TC_CASE=3 TC_SCENARIO=<scen> TC_NDAYS=<n>
%% Writes the capture to TC_CAP_OUT.
ERW_run;                       %% opens with clear all; the hook self-arms
%% Resolve the output path AFTER the run: ERW_run's clear all discards any
%% variable set before it, so this must be read from the environment here.
out = getenv('TC_CAP_OUT');
if isempty(out); out = fullfile(cd,'SIM','cap_1d.mat'); end
global ERW_CAP %#ok<GVMIS,TLEV>
cap = ERW_CAP;
cap.vals = cap.vals(1:min(cap.k,cap.max),:);
cap.source = '1D MAIN_FRAME_BG_ERW_10mins';
cap.init = ERW_CAPI_get();
fprintf('\n[Capture_1D] %d ERW_function calls, %d scalar terms\n', ...
    cap.k,numel(cap.names));
save(out,'cap','-v7.3');
fprintf('[Capture_1D] saved %s\n',out);
