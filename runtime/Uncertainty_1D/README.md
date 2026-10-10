# Kantola uncertainty ensemble

The parameter-uncertainty analysis of the T&C-SMEW field paper (Zhang et
al., GCB; supplement section 4, Table S2), carried into the coupled model
with the two parameters this model adds. Three parameter sets are
compared: the original (calibrated) values and two moved sets, each
moving all sixteen parameters to the end that pushes the basalt pH
response its way (the "lower" and the "upper" response; the bounds of
9 October 2026, `UNCERTAINTY_BOUNDS_REVISION_2026-10-09.md`). Each set
covers every treatment in 1-D and 2-D, at the standard rates and in the
5 t/ha calcite arm, over the five validation years and over the 35
projection years of each SSP.

The legacy analysis ran the ERW module standalone on the resulting
respiration range. Here the whole model runs, so the pH, nitrogen and
weathering responses carry the parameter spread together.

## Parameter sets

The sets since package ac (the first run's sets, Table S2's ends, are in
the last table of this section):

| parameter | meaning | lower (low response) | original | upper (high response) | Table S2 range |
|---|---|---|---|---|---|
| Vm | max. specific decomposition rate of MOC | 200 | 100 | 50 | 1.2-528 |
| Vpc | max. specific decomposition rate of POC cellulose | 400 | 200 | 100 | 4.8-792 |
| Km | half-saturation constant, MOC decomposition | 0.5 | 0.25 | 0.125 | 0.057-0.95 |
| rem | turnover rate of the MOC enzyme pools | 0.036 | 0.018 | 0.009 | 0.0075-0.075 |
| rep | turnover rate of the POC enzyme pools | 0.036 | 0.018 | 0.009 | 0.0075-0.075 |
| pepf | maintenance-respiration fraction to POC enzymes (fungi) | 0.012 | 0.006 | 0.0031 (S2 floor) | 0.0031-0.031 |
| mrb | maintenance rate of bacteria | 0.020 | 0.010 | 0.005 | 0.0012-0.03 |
| mrf | maintenance rate of saprotrophic fungi | 0.010 | 0.005 | 0.0025 | 0.0012-0.03 |
| fd | fraction of decomposed POC to DOC | 0.5 | 0.4 | 0.3 | 0.2-0.8 |
| Vpl | max. specific decomposition rate of POC lignin | 46 | 23 | 11.5 | 4.8-792 |
| Vdb | max. DOC uptake rate, bacteria | 0.06 (S2 ceiling) | 0.04 | 0.02 | 0.0024-0.06 |
| Vdf | max. DOC uptake rate, saprotrophic fungi | 0.04 | 0.02 | 0.01 | 0.0024-0.06 |
| FD | dissolution factor on the weathering rate law | 0.8 | 1 | 1.2 | 0.8-1.2 |
| exponent on the ERW pore saturation | s^2 slow, sqrt(s) fast | 2 | 1 | 0.5 | |
| KNmaxd | denitrification rate scale, d^-1 | 0.002 | 0.0048 | 0.0064 | |
| alpha | Appelo exchange slope, constants slid with it | 1.3 | 3.4 | 6 | |

The original values are what the model ships: the SI "original" column,
except mrb and mrf, which keep the legacy code's 0.010 and 0.005 (the SI
lists 0.005 and 0.002) because `Data_Kantola_BG.mat` was spun up with
them. The workbook's dissolution factor is 1 on every feedstock row, the
SI's original.

**The twelve literals** move a factor of two around the shipped values
(fd, a fraction, by 0.1), clipped to Table S2 where the factor would
leave it: pepf's half, 0.003, becomes the table's floor 0.0031, and
Vdb's double, 0.08, its ceiling 0.06. Table S2's ends are the cross-site
literature priors of the MEND-type microbial model (Wang et al. 2013,
Waring and Averill 2013), two to three orders of magnitude wide; they
describe ignorance across soils, not the uncertainty left after the
model has been calibrated to this field's tower NEE, respiration and
chambers. The first run of the ensemble (package u) put all twelve at
those ends on pools spun up at the shipped values, and the set at the
tops (Vm x5.3, Vpl x34 ...) mineralised about a decade of litter in its
first year: +865 g C/m2 of net respiration, 145 kg N/ha of nitrate
leached and 36 kg N2O-N/ha emitted in 2016, before any fertiliser, and
the control pH went from 6.1 to 4.9 within months. A factor of two is
the usual post-calibration bracket for microbial kinetic constants and
keeps the potential decomposition within 2x of the state the pools were
spun up with, which that run shows is the condition for no burst (x4-34
burst; x0.01-0.4 did not). The range comments in
`Biogeochemistry_Parameter.m` still give Table S2.

**KNmaxd** is this model's addition, 0.0048 d^-1 calibrated to the
multi-year chamber total (N2O/obs = 1). That total hides two
compensating per-year errors recorded in `Biogeochemistry_Parameter.m`:
2017 modelled 2.4x the chambers, 2018 0.75x, so the per-year fits
bracket the constant at 0.0048/2.4 = 0.002 and 0.0048/0.75 = 0.0064.
Both ends have literature support: Hou et al. (2024,
doi:10.3389/fpls.2024.1482292) measure first-order soil denitrification
constants of 0.002-0.003 d^-1 in the upper 40 cm; CropSyst's first-order
nitrate loss (Confalonieri et al. 2006, doi:10.1051/agro:2006022) has a
0.01 d^-1 reference at 15 C (Heinen 2003, Alterra report 690; Heinen
2006 reviews the family), so 0.0064 is inside the modelled range, and
the first run's 0.01 is dropped as the top of a cross-model bracket
rather than a value fitted to this site. 0.0064 ran in the package-t
sweep grid (`kn0p0064`), so its exchange behaviour is known.

**Alpha** 3.4 is the Appelo (1994) proton-loading slope, the value
SCEPTER carries as its default and the campaign's since package u. The
range is the ends of the sweep grid, 1.3 (the value fitted on 26
September) and 6, both formed as the sweep forms them (the first point of
`linspace(1.3,3.4,5)` and the literal 6). The field's own pH-versus-base-
saturation pairs (Beerling et al. 2024 / Pihlap, 0.063-0.080 pH per
point of base saturation against the model's 0.056 at alpha 3.4)
correspond to alpha 4.1-5.8 at the seed's f_H of 0.20, so 6 brackets the
field slope with a margin. Moving alpha slides the exchange constants by the campaign's rule,
`logKH = reference_logKH + charge x (alpha - 1.3) x 0.20`, which keeps the
seed's measured 80 % base saturation. A set's constants are given to
twelve significant figures, the precision a result name carries: at alpha
6 the Mg constant of the rule has thirteen (10.29854858734, run as
10.2985485873) and the K constant falls a rounding error away from its
twelve-figure decimal, and either would otherwise add its hexadecimal form
to the name and take the longest 35-year 2-D name past the 255 characters
a file name may have. The sweep's own alpha 6 cases carry the unrounded
constants, so the upper set's result names are not theirs:

| alpha | Na | K | Ca | Mg | Al |
|---|---|---|---|---|---|
| 1.3 | 4.716 | 3.616 | 8.103 | 8.419 | 12.919 |
| 3.4 | 5.136 | 4.036 | 8.943 | 9.259 | 14.179 |
| 6 | 5.656 | 4.556 | 9.983 | 10.299 | 15.739 |

Alpha sets how far the bulk pH moves for a given change in base
saturation, (0.434/f_H + alpha)/100 pH units per point at an H fraction
f_H: 0.035, 0.056 and 0.082 at the seed's f_H of 0.20 for the three
values.

**Which end sits in which set.** "Lower" and "upper" name the ends of
the basalt pH response, the quantity the paper reports, not the ends of
the parameter ranges: each parameter's end sits in the set it pushes
that response towards. Faster mineralisation (literals x2) means more
nitrification acidity and nitrate export, so a lower pH; less
denitrification (KNmaxd 0.002) consumes less acidity and nitrate; a low
alpha moves the pH less per point of base saturation; a low dissolution
factor and s^2 release less. So the lower set is the fast-nitrogen,
low-alpha, slow-weathering set and the upper set the slow-nitrogen,
high-alpha, fast-weathering one:

| parameter | lower (low response) | upper (high response) |
|---|---|---|
| twelve literals | x2 | x0.5 |
| KNmaxd | 0.002 | 0.0064 |
| alpha | 1.3 | 6 |
| dissolution factor | 0.8 | 1.2 |
| saturation exponent | 2 (s^2) | 0.5 (sqrt s) |

Relative to the first run this flips only the literals' side. In that
run both moved sets were the ends of the parameter ranges, so
parameters that push the pH opposite ways shared a member and cancelled:
the only member above the central control line was the set with alpha
1.3, and alpha 6, which gives the largest basalt effect, sat in the
acidified member where it could not lift the band. Alpha alone cannot do
it: it scales the control's acid drift and the basalt gain alike (0.035,
0.056 and 0.082 pH per point of base saturation at 1.3, 3.4 and 6), so
the band's upper edge for the control is set by the nitrogen cycle and
the basalt plot's upper edge needs high alpha in that same slow-nitrogen
member. Other outputs need not follow the sets: N2O pairs fast
mineralisation with slow denitrification and vice versa, so its band is
narrow (roughly 10-45 kg N/ha over five years, inside the Midwest maize
range of Hoben et al. 2011), and the original can fall outside a pair.

The first run's sets, for the record (package u to ab; every parameter
at a Table S2 end, KNmaxd 0.002-0.01, alpha 1.3 in the lower set and 6
in the upper):

| parameter | lower | upper |
|---|---|---|
| Vm, Vpc, Km | 1.2, 4.8, 0.057 | 528, 792, 0.95 |
| rem, rep, pepf | 0.0075, 0.0075, 0.0031 | 0.075, 0.075, 0.031 |
| mrb, mrf, fd | 0.0012, 0.0012, 0.2 | 0.03, 0.03, 0.8 |
| Vpl, Vdb, Vdf | 4.8, 0.0024, 0.0024 | 792, 0.06, 0.06 |
| FD, exponent, KNmaxd, alpha | 0.8, 2, 0.002, 1.3 | 1.2, 0.5, 0.01, 6 |

## Arms, horizons and where the results are

In package u the original set is the Base and SSP campaigns of the same
root (`TC_UC_ORIGINAL_ROOT="$PWD"`, no `TC_UC_WITH_ORIGINAL`), so it is run
once. Results of another package cannot stand in for it: the model check
below refuses a checkout whose source or inputs differ. `TC_UC_WITH_ORIGINAL=1`
remains for an ensemble kept in a root of its own, as m ran it. The lower
and upper sets run here, each over every run of an arm, with everything
else as the original's campaign ran it: the campaign's set
(`Campaign_Parameters.m`), the horizon, forcing and application years,
and for 2-D the worker count (6 at five years as `Base_2D`, 12 at 35 as
`SSP_2D`).

| horizon | arm | runs | original | lower and upper |
|---|---|---|---|---|
| 5 yr, legacy forcing | 1-D | 10, 11, 12 | `Base_1D` | `Uncertainty_1D` |
| | 1-D, 5 t/ha | 12 | `Base_1D_case8` | `Uncertainty_1D_case8` |
| | 2-D | 1-9 | `Base_2D` | `Uncertainty_2D` |
| | 2-D, 5 t/ha | 3, 6, 9 | `Base_2D_case8` | `Uncertainty_2D_case8` |
| 35 yr, SSP forcing | 1-D | 10, 11, 12 | `SSP_1D_35yr_<ssp>` | `Uncertainty_1D_35yr_<ssp>` |
| | 1-D, 5 t/ha | 12 | `SSP_1D_35yr_<ssp>_case8` | `Uncertainty_1D_35yr_<ssp>_case8` |
| | 2-D | 1-9 | `SSP_2D_35yr_<ssp>` | `Uncertainty_2D_35yr_<ssp>` |
| | 2-D, 5 t/ha | 3, 6, 9 | `SSP_2D_35yr_<ssp>_case8` | `Uncertainty_2D_35yr_<ssp>_case8` |

All are directories under `v2/SIM/Subfield_12/`. The table's original
column gives the matching Base/SSP campaign directories, which hold the
original set in package u. `lower/` and `upper/` sit under each ensemble
tag (with `original/` beside them only when a call adds
`TC_UC_WITH_ORIGINAL=1`), with the 2-D logs, status files and run locks
under the same path in `v2/logs/Subfield_12/`.

The three sets compare only if they ran the same model, so the controller
first compares this checkout with the original's (`same_model_code.sh`):
`v2/T&C_Code` without its tests, the forcing and every other input in it,
`v2/run`, and the calibrated set must be byte for byte the same. Nothing
is submitted otherwise.

## Files

| file | role |
|---|---|
| `Uncertainty_Members.m` | the three sets above |
| `Uncertainty_Campaign.m` | reads `TC_SSP`, `TC_HORIZON` and `TC_CALCITE_CASE` into the arm, horizon, tags, run list, workers and array elements (member-major) |
| `Uncertainty_Knobs.m` | a set's values and its `Run_One_Subfield` options; for the original, exactly what the Base and SSP tasks pass |
| `Uncertainty_Result_Matches.m`, `Uncertainty_Result_State.m` | whether a saved result is the expected run: campaign, stamp, horizon, columns, every knob per column, and the scenario and calcite case in its name |
| `Run_Uncertainty_Task.m`, `submit_uncertainty.sbatch` | one 1-D element: one CPU, one thread |
| `../Uncertainty_2D/Run_Uncertainty_2D_Task.m`, `../Uncertainty_2D/submit_uncertainty_2d.sbatch` | one 2-D element: the SSP 2-D campaign's path, with a run lock and the lateral conservation check |
| `submit_uncertainty.sh` | the controller: the model check, then the 1-D and 2-D arrays of one horizon and case |
| `same_model_code.sh` | the byte comparison with the original's checkout |
| `requeue_near_limit.sh` | sourced by both batch files: a run requeues itself before the 24-hour limit and continues from its checkpoint |

## Submitting

`RUNBOOK_2026-09-24.md`, section 0u, gives the order for package u: the
equivalence baseline, the campaign's own 1-D runs with their validation,
then the ensemble. Submit the ensemble only after the R2025a validation of
the root has passed (the `Base1Dcheck` job). Section 0ac gives the order
for the revised sets: the two moved sets run in 1-D first
(`TC_UC_DIMS=1d`, runs 10-12, five years) and are accepted against the
checks listed there (`precheck_1d.m` at the package root prints them)
before any 2-D or 35-year array is submitted.

One controller call per horizon and case, eight in all. A call submits a
1-D and a 2-D array; `TC_UC_DIMS=1d` or `TC_UC_DIMS=2d` submits that array
alone, so the five-year 1-D sets can run and be looked at first:

```bash
for c in '' 8; do TC_BOUCHET_ACCOUNT=pi_np363 TC_UC_ORIGINAL_ROOT="$PWD" TC_UC_DIMS=1d TC_SSP= TC_HORIZON= TC_CALCITE_CASE=$c bash runtime/Uncertainty_1D/submit_uncertainty.sh || break; done
```

```bash
for c in '' 8; do TC_BOUCHET_ACCOUNT=pi_np363 TC_UC_ORIGINAL_ROOT="$PWD" TC_UC_DIMS=2d TC_SSP= TC_HORIZON= TC_CALCITE_CASE=$c bash runtime/Uncertainty_1D/submit_uncertainty.sh || break; done
```

```bash
for c in '' 8; do for s in ssp126 ssp245 ssp585; do TC_BOUCHET_ACCOUNT=pi_np363 TC_UC_ORIGINAL_ROOT="$PWD" TC_SSP=$s TC_HORIZON=35yr TC_CALCITE_CASE=$c bash runtime/Uncertainty_1D/submit_uncertainty.sh || break 2; done; done
```

Resources per element, set by the controller:

| horizon | 1-D | 2-D |
|---|---|---|
| 5 yr | 1 CPU, 12 GiB, `day`, 12 h | 6 CPUs, 48 GiB, `day`, 24 h |
| 35 yr | 1 CPU, 56 GiB, `day`, 24 h per piece | 12 CPUs, 256 GiB, `day`, 24 h per piece |

Everything runs on the `day` partition. The per-user limits there, 1500
CPUs and 20 TB, are far wider than on `week` (96 CPUs and 1.5 TB, where
a 256 GiB element runs six at a time), so the ensemble does not queue
behind the Base and SSP campaigns. A day job stops at 24 hours. A 35-year
run longer than a day requeues itself when Slurm signals
the limit half an hour ahead, and continues from its last checkpoint
under the same job id, which keeps its run lock (`requeue_near_limit.sh`,
at most five restarts). Each restart loses at most the half hour since
that checkpoint. Interrupted 1-D and twelve-column 2-D runs were checked
bit for bit in k; the daily state replay and reload checks of l are recorded in
`PARTICLE_HISTORY_REVIEW_2026-09-28.md`. One
diagnostic does not carry over a restart: the plant carbon-overdraw count
a run prints and saves at its end (`carbon_deficit`) covers its last
piece only.

Without originals, each horizon has 8 one-dimensional and 24
two-dimensional elements (the standard and 5 t/ha arms together), 128
elements across both horizons and all SSPs. The six 35-year calls include
72 two-dimensional and 24 one-dimensional elements. Allow space for
their results and checkpoints, including temporary checkpoint copies:
a 35-year 2-D result is about 11 GB for a control and tens of GB for
whole-field basalt, so the 35-year sets belong in a scratch copy of the
checkout (RUNBOOK section 0v). The controller refuses a 5-year call below
100 GiB free and a 35-year call below 600 GiB, a floor only. The group
quota can be tighter than `df` shows. The 35-year 2-D elements dominate
the cost. Package u's 35-year runs measured what sets their memory: the
particle-class histories, 17 GiB per treated column held whole, put a
whole-field basalt run near 400 GiB; package v keeps those histories
daily (`PARTICLE_RECORDS_DAILY_2026-10-04.md`), which puts it near
150-200 GiB, hence 256 GiB (`TC_UC_2D_MEM` overrides it). Read the
MaxRSS of a finished whole-field basalt run before going lower. An element that runs out of
memory stops (OUT_OF_MEMORY) with its checkpoint; resubmitting the call
with more memory, once its arrays have left the queue, continues it.

`TC_UC_WITH_ORIGINAL=1` on a controller call runs the original set here
as well, first in each array, with the campaign's knobs and under the
same checks, into `original/` beside `lower/` and `upper/`. Package u does
not need it. The array grows by half; a call already submitted
without it only needs the
same call with it added once its jobs have left the queue, since
finished lower and upper elements exit at once.

After a TIMEOUT or node loss, resubmit the same command. An element whose
result already carries its set's provenance exits at once (a 2-D element
first rereads its lateral conservation), the rest continue from their
checkpoints, and a result of anything else is refused rather than
overwritten. The controller refuses while any of that call's jobs are
queued, and lists stale 2-D run locks to remove.

## Summary and figures

From package z, `Summarize_Uncertainty.m` at the package root gathers the
three sets of an arm: the original from the Base or SSP campaign
directory of the table above, `lower/` and `upper/` under the ensemble
tag, as pairs `{original_dir, ensemble_dir}`. It writes per run the
domain means on the daily axis, the yearly statistics, the envelope of
the three sets and the observations, which `v2/plotting/Plot_Uncertainty.m`
draws, and reads one-column results from their full records and 2-D
results of package y from their daily stores. It lives outside every
campaign's stamp, so it can change while elements run; nothing under
`runtime/`, `v2/plotting`, `v2/audit` or `v2/T&C_Code/Tests` of a root
may, since every element's stamp covers those files and a change discards
the runs in flight. For example, the five-year 2-D arm:

```bash
sbatch --account=pi_np363 --export=ALL,UC_SOURCES="v2/SIM/Subfield_12/Base_2D|v2/SIM/Subfield_12/Uncertainty_2D",UC_RUNS="1:9",UC_OUTPUT=v2/SIM/Subfield_12/uncertainty_2d_5yr.mat summarize_uncertainty.sbatch
```

The field-grid ensemble is `runtime/Field_Uncertainty_2D`.
