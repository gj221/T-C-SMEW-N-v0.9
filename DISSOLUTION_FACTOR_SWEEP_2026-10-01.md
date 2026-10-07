# Package p: the dissolution factor joins the sweep, 1 October 2026

Package p is package o with one change of scope and no change of physics:
the 1-D calibration sweep gains a fourth dimension, the dissolution factor
(`DissolutionFactor`, Run_One_Subfield's option, the workbook `diss_f` of
every applied feedstock mineral), and the other three grids are trimmed so
the sweep keeps its earlier cost. Every model file is o's.

(Packages p and q called this quantity the release factor and this note
`RELEASE_FACTOR_SWEEP_2026-10-01.md`. From package r on it is the
dissolution factor, the name of the option and of the uncertainty
ensemble's parameter. Nothing but the wording changed.)

(Package t, 2 October 2026, takes the dissolution factor out of the sweep
again and holds it at 1. The excess release this note describes was traced
to one rate constant, the amphibole's, which t corrects:
`AMPHIBOLE_KINETICS_2026-10-02.md`. The grid, counts and smoke cases below
are those of p to s; t runs o's. The factor remains the uncertainty
ensemble's parameter, 0.8 and 1.2 around 1.)

## Why

Package o's basalt run releases 12.6 eq m-2 of alkalinity by August 2020
against an exchanger acidity of 10.3 eq m-2 (f_H 0.20 of a 50.5 eq m-2
CEC). The exchanger saturates in 2019 and the bulk pH climbs to 7.78 while
the field's 0-30 cm August pH is 6.65 (treatment effect +1.94 against
+0.49). This is a mass balance: no alpha, logKH, KNmaxd or Kbot changes the
size of the acidity reservoir or the amount released, so the o sweep could
move the pH at a given base saturation but could not keep the exchanger
from saturating. The field's own soil chemistry (Beerling et al. 2024 pH
and base saturation 2017-2020; Pihlap et al. 2025, pH 7.2 in 1:1 water and
exchangeable H 2 % in October 2022 after two further applications) is
consistent with roughly 0.3 to 0.6 of the model's release, and bicarbonate
export cannot absorb the difference below pore pH 7 (alkalinity under
0.5 meq L-1 at 0.4-3 % soil CO2). The REE/TiCAT release estimates are not
used as a target: they are contested (Reershemius and Suhrhoff 2024) and
cannot be reconciled with the soil pH and base saturation.

The grid: alpha 1.3, 2.35, 3.4, 4.5, 6; KNmaxd 0.0032, 0.0048, 0.0064 d-1;
Kbot free drainage and 0.03 mm h-1; dissolution factor 1, 0.7, 0.5, 0.35,
0.25; runs 10 to 12. 150 combinations, 450 trajectories, about 1,125
CPU-hours, about 24 hours at 48 concurrent slots. Case 1 to 3 remain the
reference tuple with factor 1, so the one-year control gate repeats o's
(residual -0.0248) and the all-12 campaign label is
`a1p3_kn0p0048_kbfree_fd1`. The last case (alpha 6, KNmaxd 0.0064, Kbot
0.03, factor 0.25, calcite) is the pilot.

## What is deliberately not changed

(This note was written for package p. Package q changes one of the items
below, the paste reading of plots that hold solid calcite; see
`PASTE_CALCITE_2026-10-01.md`.)

- The export operator. The biogeochemistry and the ERW leach the 0-300 mm
  pools with the leakage out of the bottom of the 5 m column over the whole
  column's water (`Biogeo_Leakage(Ptm1,Lk,V,...)` with `V = sum(V)`;
  `f_Lw_bio = Lk_tot/VSUM`). That is inherited unchanged from T&C-BG and the
  published T&C-SMEW, whose Kantola rows already carried the 5 m mesh with
  Zbio 300 mm. Measured on the 2016 control: bottom leakage 333 mm, flux
  across 300 mm about 425 mm, column water 9.9 times the active layer's,
  applied fraction 0.3 of the pool per year against 4.4. By decision of
  1 October 2026 it stays as inherited for this package; the dissolution factor
  the sweep selects therefore absorbs the missing export as well as the
  lab-to-field rate gap, and the 35-year export of alkalinity (the carbon
  removal) keeps the inherited rate. A corrected operator (the net Darcy
  flux across the active-layer base, prototyped and smoke-tested in the
  session scratchpad) is held back for a later package.
- The calcite chemistry. Applied calcite is primary, dissolves under the
  Palandri-Kharaka law with (1 - Omega) on the activity basis and a per-step
  closed-water equilibrium limiter, and the secondary CaCO3 pool
  precipitates above saturation. 50 t/ha therefore saturates the pore
  water within days and the pH sits at the calcite-CO2 ceiling for the
  prevailing soil pCO2 and Ca (pore 6.9-7.35, bulk 7.4-7.9 under the closed
  Kanzaki-2024 extract); the basalt run reaches the same ceiling in 2020.
  More calcite cannot raise it; the factor changes calcite dissolution only
  while the water is undersaturated. `calcite_check.m` at the package root
  reads the carbonate state of any result file; it sits outside every
  fingerprinted directory, so running it or copying it elsewhere does not
  touch a validation stamp. Never add a file under `runtime/<campaign>` or
  `v2/` of a validated root.
- The extract operator and the seed's background Ca, which set where that
  ceiling reads (a calcareous soil reads 7.8-8.3 in an open paste). Changing
  either moves the seed inversion and the whole bulk-pH scale.

## Files changed from o

`runtime/Calibration_Sweep_1D/Sweep_Cases.m` (grid and tags),
`Run_Sweep_Task.m` (option and log line), `Validate_Sweep.m` (smoke cases,
option, provenance and filename checks), `Summarize_Sweep.m` (case count,
reference check, DissolutionFactor column in both row builders and the
combination table), `Assert_Sweep_Provenance.m` (`dissolution_factor`
required and compared), `submit_sweep.sh` (array sizes 450 and 150,
rescore counts), `README.md`; `v2/T&C_Code/Tests/test_bouchet_1d_joint_sweep_configuration.m`
(grids, counts, tag pattern, tuple key, array literals);
`runtime/Base_1D/Campaign_Parameters.m` (combination label). The pins are
regenerated. No file under `v2/T&C_Code` other than that test changes, so
the model and its 93 other audit tests are o's.

## Validation

Validation runs on Bouchet through the sweep's own `SweepCheck` job (the
full post-legacy audit and the three-case smoke under the campaign
identity), then the one-year control gate. Locally, as a sanity check of
the edits only: `checkcode` on every edited file; the tests whose text or
inputs the edits touch (`test_bouchet_1d_joint_sweep_configuration`,
`test_all_12_configuration`, `test_scepter_cec_alignment`,
`test_pH_conversion`, `test_bouchet_linux_portability`); and a replica of
the three-case validation smoke, three 3-day runs with the new option
(result file names 168 to 205 characters). The audit and every five-year
trajectory run on Bouchet.

## Reading the result

`joint_sweep_combination_summary.csv` gains a `DissolutionFactor` column.
Read `Front_pH_difference_base_sat_difference` first, as before, and
require the chosen combination to hold the August 2019 and 2020 effects
(+0.50, +0.49) with base saturation 83-87 %. Then run `budget.m` on its
control and basalt results and `calcite_check.m` on its calcite result. A
combination adopted for production goes into
`runtime/Base_1D/Campaign_Parameters.m` (alpha, KNmaxd, Kbot and the
factor) in a rebuilt package, and the 5-year, 2-D and 35-year campaigns run
in that root.

## Relation to the uncertainty ensemble

(Added in package r.) The uncertainty ensemble (`runtime/Uncertainty_1D`,
after Table S2 of the T&C-SMEW field analysis) moves the dissolution factor
too: 0.8 in the lower set and 1.2 in the upper set, against the workbook's
1 in the original. Those are absolute values, handed to the run as
`DissolutionFactor`. Two things follow.

The ensemble spreads the result around the central value; it does not
choose it. Taken alone, a factor of 0.8 still delivers 10.1 eq m-2 by
August 2020 against 10.3 eq m-2 of exchanger acidity. (The lower set also
slows weathering through the saturation exponent and moves fourteen other
parameters, so its own outcome is not that number.) The sweep is what
looks for the central value.

If a dissolution factor other than 1 is adopted, the ensemble has to move
with it. The Base and SSP tasks must pass the adopted value, and the lower
and upper sets must become 0.8 and 1.2 times it (`Uncertainty_Members.m`,
`Uncertainty_Knobs.m`, the table of `runtime/Uncertainty_1D/README.md`,
`test_uncertainty_campaign`). Left as they are, both sets would run at 0.8
and 1.2 of the workbook's rate whatever the original uses. None of this is
changed in this package, because no factor has been adopted yet.
