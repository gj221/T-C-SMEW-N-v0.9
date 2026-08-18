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
%% Known infeasible corner: this set crossed with basalt (run 11) drives
%% lsqnonlin in ERW_function past convergence from roughly iteration 37000
%% onward, thousands of consecutive failures with no successes. The N
%% levers here change nitrification H+ production, which moves soil pH,
%% which is an input to the geochemical equilibrium solve, so they can push
%% the chemistry into a region with no attainable solution. It fails
%% quietly: ERW_function only disp()s a warning and then clamps the
%% unconverged answer, so the run continues on unreliable pH, alkalinity
%% and Ca and slows about 16x. Treat the cell as not simulable rather than
%% re-running it; res_threshold or the x0 restart needs revisiting first.
%% The same set with calcite, and basalt at stock parameters, are clean.
sets{end+1,1} = mk('den_sharp_nit_ph5','pH_den_rng 2.0 with pH_nit_0 5.0', ...
    'pH_den_rng',2.0,'pH_nit_0',5.0);
sets{end+1,1} = mk('den_sharp_nit_ph6','pH_den_rng 2.0 with pH_nit_0 6.0', ...
    'pH_den_rng',2.0,'pH_nit_0',6.0);
%% Infeasible corner, second occurrence: this set crossed with basalt
%% (run 11) fails the same way as den_sharp_nit_ph5 x basalt, in the
%% warning-and-clamp mode documented there. It ran clean through ~64%
%%
%% (step ~150000) before the ERW_function solver left convergence from
%% roughly step 212419 onward with thousands of consecutive failures and
%% no successes; audited 2026-08-17, killed and recorded failed rather
%% than left to crawl at ~16x on clamped pH/alkalinity/Ca. The common
%% trigger with the sharp corner is pH_nit_0 = 5.0 crossed with basalt:
%% the pH 6.0 corners of the same span are clean. Treat as not simulable
%% unless res_threshold or the x0 restart is revisited.
sets{end+1,1} = mk('den_flat_nit_ph5','pH_den_rng 4.0 with pH_nit_0 5.0', ...
    'pH_den_rng',4.0,'pH_nit_0',5.0);
sets{end+1,1} = mk('den_flat_nit_ph6','pH_den_rng 4.0 with pH_nit_0 6.0', ...
    'pH_den_rng',4.0,'pH_nit_0',6.0);

%% Third level on both levers, completing the 3x3 grid. These were added to
%% support a quadratic fit locating an optimum for the 2D runs, and that
%% purpose does not hold: the optimum is not identifiable here. Choosing
%% these parameters needs N2O or nitrate-leaching observations to fit
%% against and the repository has none, the measured data being GPP, NEE
%% and Re for Kantola and CEC, Ca, Mg and pH from Beerling. The chemistry
%% observations cannot discriminate either, because these levers move the
%% response to pH rather than pH itself: control mean pH came out at 6.049,
%% 6.060 and 6.057 across sets. Absolute N2O is additionally degenerate
%% with the rN2O partition, which is not swept. The 2D runs therefore use
%% the published defaults, pH_nit_0 5.0 with pH_den_rng 3.0, which is the
%% baseline set above.
%%
%% They are kept because they ran, and because they widen the spread the
%% sweep exists to measure. They also settle where the solver fails: every
%% one of these completed clean, basalt included, so the infeasibility
%% documented above belongs to pH_nit_0 5.0 crossed with a perturbed span,
%% not to the 5.5 level and not to basalt generally. Seven of the nine
%% basalt cells are simulable; control and calcite complete the full grid.
sets{end+1,1} = mk('den_sharp_nit_ph55','pH_den_rng 2.0 with pH_nit_0 5.5', ...
    'pH_den_rng',2.0,'pH_nit_0',5.5);
sets{end+1,1} = mk('den_mid_nit_ph55','pH_den_rng 3.0 with pH_nit_0 5.5', ...
    'pH_den_rng',3.0,'pH_nit_0',5.5);
sets{end+1,1} = mk('den_mid_nit_ph6','pH_den_rng 3.0 with pH_nit_0 6.0', ...
    'pH_den_rng',3.0,'pH_nit_0',6.0);
sets{end+1,1} = mk('den_flat_nit_ph55','pH_den_rng 4.0 with pH_nit_0 5.5', ...
    'pH_den_rng',4.0,'pH_nit_0',5.5);
end

%% =======================================================================
function s = mk(tag,note,varargin)
s = struct('tag',tag,'note',note);
for k = 1:2:numel(varargin)
    s.(varargin{k}) = varargin{k+1};
end
end
