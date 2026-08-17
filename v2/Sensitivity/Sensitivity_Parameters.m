function [names,sets] = Sensitivity_Parameters()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sweepable nitrogen parameters, and the parameter sets to run.
%%
%% names  cell row of BiogeoPar field names the sweep is allowed to move
%% sets   cell column, one struct per run set: .tag .note and one field per
%%        parameter it changes. Absent fields keep the model default. A cell
%%        array rather than a struct array because the sets deliberately
%%        carry different fields, which a struct array cannot hold.
%%
%% Only two of these need a simulation. rN2O re-partitions DENITflx after
%% the nitrate balance has already taken it, so pH0_rN2O and k_rN2O cannot
%% move NO3 at all and their effect on N2O is recovered exactly in
%% post-processing by Analyse_Sensitivity. They are listed as sweepable for
%% completeness but should be left out of the run sets: sweeping them here
%% multiplies the compute for an answer arithmetic already gives.
%%
%% Ranges below are bracketing values, not recommendations. The point of
%% the exercise is the spread they produce, not picking the member that
%% flatters the treatment.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

names = {'KNmaxn','KNmaxd','Kc','f_N2O_nit','k_rN2O','pH0_rN2O', ...
    'pH_nit_0','pH_nit_s','pH_den_min','pH_den_rng'};

if nargout < 2; return; end

%% Baseline first so every comparison has an anchor in the same batch.
sets = {};

sets{end+1,1} = mk('baseline','model defaults, anchor for the batch');

%% Full factorial over the two pH-response shape levers, at two bracketing
%% levels each, with the rate constants at stock. Denitrification responds
%% to the feedstock pH rise through the span of its linear pH response
%% (Heinen 2006 reviews linear forms with differing thresholds and spans),
%% and the nitrification arctan response (Parton et al 1996) is governed by
%% where its inflection sits: at 5.0 it is nearly saturated across the
%% working pH of a fertilised midwest maize soil, at 6.0 it is mid-ramp
%% where the feedstock raises pH. The two levers act on opposite ends of
%% the same N cycle, so the full factorial belongs in the batch, the
%% interaction is the point. A third lever would double the set count, so
%% the factorial is kept at the minimum that answers the question.
sets{end+1,1} = mk('den_sharp_nit_ph5','pH_den_rng 2.0 with pH_nit_0 5.0', ...
    'pH_den_rng',2.0,'pH_nit_0',5.0);
sets{end+1,1} = mk('den_sharp_nit_ph6','pH_den_rng 2.0 with pH_nit_0 6.0', ...
    'pH_den_rng',2.0,'pH_nit_0',6.0);
sets{end+1,1} = mk('den_flat_nit_ph5','pH_den_rng 4.0 with pH_nit_0 5.0', ...
    'pH_den_rng',4.0,'pH_nit_0',5.0);
sets{end+1,1} = mk('den_flat_nit_ph6','pH_den_rng 4.0 with pH_nit_0 6.0', ...
    'pH_den_rng',4.0,'pH_nit_0',6.0);

end

%% =======================================================================
function s = mk(tag,note,varargin)
s = struct('tag',tag,'note',note);
for k = 1:2:numel(varargin)
    s.(varargin{k}) = varargin{k+1};
end
end
