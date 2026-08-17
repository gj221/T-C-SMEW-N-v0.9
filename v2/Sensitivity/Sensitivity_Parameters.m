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

%% Denitrification pH sensitivity. Heinen 2006 reviews linear forms with
%% differing thresholds and spans; a narrower span makes denitrification
%% respond more sharply to the feedstock pH rise, which is the mechanism
%% that converts nitrate to gas rather than leachate.
sets{end+1,1} = mk('den_pH_sharp','denitrification pH span 3.0 -> 2.0', ...
    'pH_den_rng',2.0);
sets{end+1,1} = mk('den_pH_flat','denitrification pH span 3.0 -> 4.0', ...
    'pH_den_rng',4.0);

%% Denitrification rate constant. Dickinson et al 2002 give 0.217 /d.
sets{end+1,1} = mk('den_fast','KNmaxd 0.217 -> 0.325 /d','KNmaxd',0.325);
sets{end+1,1} = mk('den_slow','KNmaxd 0.217 -> 0.145 /d','KNmaxd',0.145);

%% Nitrification. Slower nitrification leaves N as NH4, which neither
%% leaches nor denitrifies appreciably, so it cuts both loss terms at once.
sets{end+1,1} = mk('nit_slow','KNmaxn 0.086 -> 0.057 /d','KNmaxn',0.057);
sets{end+1,1} = mk('nit_pH_sharp','nitrification inflection 5.0 -> 6.0', ...
    'pH_nit_0',6.0);

%% One combination, to test whether the two mechanisms compose or cancel.
sets{end+1,1} = mk('den_sharp_nit_slow','pH_den_rng 2.0 with KNmaxn 0.057', ...
    'pH_den_rng',2.0,'KNmaxn',0.057);

end

%% =======================================================================
function s = mk(tag,note,varargin)
s = struct('tag',tag,'note',note);
for k = 1:2:numel(varargin)
    s.(varargin{k}) = varargin{k+1};
end
end
