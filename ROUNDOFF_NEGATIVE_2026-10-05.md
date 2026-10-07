# A rounding-level negative in the ionic strength, 5 October 2026

Package w is package v with one change to `ERW_Ionic_Strength.m` and its
test (`test_erw_activity_coefficient.m`). The comparison utility at the
package root (`compare_u_v.m`, outside the model and its audit) also
handles two results that both keep daily particle rows, so it can compare
w with v.

## What happened

The 35-year 2-D uncertainty run ssp245, upper set, run 7 (prairie margin
control; element 28314649_16 of package v) stopped after 20.5 hours, at
59 % of its run (model date October 2036):

    state=failed
    detail=ERW:negativeConcentration
    message=Ionic strength inputs must be nonnegative concentrations;
            most negative is -4.54747e-16.

One of the sixteen pore-water concentrations handed to the ionic strength
was -4.5e-16 umol-conv/L. That is a few units in the last place of numbers
of order one, a cancellation residue rather than a negative amount: the
Al(OH)4- concentration, for one, is the dissolved Al left after its other
species (`AlOH4 = Al_w - (Al + AlOH + AlOH2 + AlOH3)` in `ERW_function`),
and a pool emptied by a flux of its own size can end a few ulps below zero.
The check refused any negative value, so the run stopped. The model is
deterministic: resuming from its checkpoint stops at the same step.

## The change

Inputs above -1e-10 umol-conv/L (1e-16 mol/L, far below any concentration
the model resolves) count as zero; anything more negative still stops the
run with the same error. Nonnegative inputs take exactly the same path as
before, so every run that got past this check gives the same numbers in v
and in w. The test pins both sides: -4.5e-16 gives the ionic strength of
zero exactly, and -2e-10 is refused.

## Using it

w's code differs from v's, so a run of v cannot continue from its
checkpoint in w; the failed run starts again from 2016 (about 36 hours on
the day partition). RUNBOOK section 0w deploys w, runs the full audit and
the one-year comparison against v (which must pass, as v's did against u),
and reruns the failed element from a scratch copy of w. Runs of v that
have not stopped need nothing.
