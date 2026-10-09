# The 35-year runs on the 10 m field grids, 8 October 2026

Package aa is package z (package y, `FIELD_RECORDS_DAILY_2026-10-08.md`,
with the ERW chemistry kernels, `ERW_KERNELS_2026-10-08.md`) with the
campaign side of the 108 35-year runs on the 10 m Kantola grids: 36
central runs and the 72 lower and upper runs of the uncertainty ensemble.
No equation, parameter, input or solver setting changes. One model-tree
file changes (`v2/run/Run_One_Subfield.m`), in how it names a result or
checkpoint that would not fit in a file name; every name that fits is
unchanged. The runs use z's fast kernels (`TC_ERW_KERNEL=fast`), which
compute every value as y's code does, bit for bit.

## 1. The runs

| | runs | grid | calls |
|---|---|---|---|
| central, standard | 1, 4, 7 | `kantola_control_dem_10m` (369 squares) | 3 scenarios |
| central, standard | 2, 3, 5, 6, 8, 9 | `kantola_basalt_dem_10m` (340 squares) | 3 scenarios |
| central, 5 t/ha calcite | 3, 6, 9 | `kantola_basalt_dem_10m` | 3 scenarios |
| lower and upper sets | the same 36 | the same | 9 calls |

36 central runs (three scenarios, nine standard runs and three calcite
runs each) and 72 ensemble runs. The table puts the controls on the
control field and the basalt and calcite runs on the basalt field, as the
two fields of Kantola et al. (2023) were run; that assignment is a choice
in the commands only (`TC_FIELD_RUNS`), not in the code, and any run can
go on either grid.

* **Central** (`runtime/Field_SSP_2D`): runs of `Run_Matrix_12` on a field
  grid under the scenario forcing, 2016-01-01 to 2050-12-31, rock every
  10 November to 2050, with the calibrated set (`Campaign_Parameters`),
  as `runtime/SSP_2D` runs them on the 12-column plot with the field grid
  in place of the 2x2 grid. These are the original set of the ensemble.
* **Ensemble** (`runtime/Field_Uncertainty_2D`): the lower and upper sets
  (`Uncertainty_Members`: the sixteen parameters at the ends of their
  ranges, the dissolution factor 0.8 and 1.2, alpha 1.3 and 6) on the same
  runs and grids, as the original's campaign runs them.

Both campaigns pass Field_2D's gates: the grid file beside the deployed
root (its content in the stamp), `TC_FIELD_WORKERS` workers for the
squares' daily, hourly and ERW work, a validation element (the complete
audit, the Linux checks, one day of the grid across the first
application under the campaign's scenario, case and options, saved and
read back), the stamp checked before and recomputed after each run, the
run lock, two-hourly checkpoints, the calibrated (or the set's)
provenance, the grid file's SHA-256 in the result, and the lateral
conservation read back.

## 2. File names past 255 characters

`Run_One_Subfield` names a result or checkpoint by its run, scenario,
case and grid, then one token per knob passed, then the span. On the field
grids the grid's tag (`field_kantola_control_10m__cell10m__margin20__dem10`,
51 characters) takes the place of `grid2x2` (7), and every lower and upper
name, result and checkpoint, of every run, grid, scenario and case comes
out at 258 to 272 characters (the longest: run 9 of the 5 t/ha arm under
ssp585 on the control grid). A Linux file name may have 255. A run would have checkpointed into a
temporary file that could not be moved into place, warned
(`ERW:CheckpointFailed`) and carried on with no checkpoint, then failed to
move its finished result into place after 35 years and deleted it.

`v2/run/Bounded_File_Tag.m` keeps a name within 255 characters: a name
that would not fit keeps its run, scenario, case, grid and span tokens and
has its knob tokens replaced by `__knobs` and twelve hex digits of their
SHA-256 (`run_provenance` carries every knob, and the ensemble checks a
result by its provenance). A name that fits is returned as it is, so the
central field runs (at most 238 characters), the plot ensemble (227) and
every earlier campaign keep their names. Run 9 of the upper set reads
`TC_2D_run09_calcite__margin_prairie__ssp585__case8__field_kantola_control_10m__cell10m__margin20__dem10__knobs<12 hex>__12784d__from20160101.mat`
(148 characters).

## 3. Memory, time and disk

* **Memory.** Package y holds every daily array of every square whole from
  the start (672 values a square and day), with the 12 hourly-derived
  daily stores and the 79 values of each ERW day's summary: 6.1 kB a square
  and day, 26.8 GiB on 369 squares over 12,785 daily rows (24.7 GiB on
  340), resident by the end of the run. With about 4 GiB for MATLAB, the
  pool and the per-square state, a task refuses a job below 30.8 GiB on
  369 squares or 28.7 GiB on 340 (`Field_Campaign_Grid('memory')`). The
  first ERW day also builds each square's whole-run ten-minute timetable
  on its worker, a passing peak of about one timetable per worker that no
  measurement yet bounds; the pilot's MaxRSS covers it.
* **Time.** No 35-year field run has been timed on Bouchet. Package x's
  one-column five-year run took 2 h 14 min, 4.4 s a square and day, about
  nine tenths of it the ERW chemistry; y does that work for every square in
  its three parallel loops, so with y's own kernels a run is about 5,800
  CPU-hours on 369 squares (7.5 days on 32 workers if the loops scale) and
  the 108 runs about 0.6 million. z's fast kernels run the ERW step 11
  times faster on Bouchet (15 times on y locally), a one-column year 4.8
  times faster and a 12-square three-day run on four workers 2.6 times, so
  a field run should take roughly 1.5 to 3 days on 32 workers and the 108
  runs 0.1 to 0.2 million CPU-hours. The pilot measures it.
* **Partition.** The `day` partition by default: 1500 CPUs and 20 TB per
  user, so about 46 runs of 32 workers at once, against 96 CPUs per user on
  `week`. Pieces are 24 hours; a run requeues itself half an hour before
  the limit (`requeue_near_limit.sh`), at most `TC_FIELD_RESTARTS` times
  (default 12, about 12.5 days of computing), and continues from its
  two-hourly checkpoint under the same job id and lock, repeating up to
  two hours each time. `TC_FIELD_PARTITION=week` gives pieces of up to 7
  days.
* **Disk.** A checkpoint holds the daily stores, ERW day summaries and
  state uncompressed, about 9 GB on 369 squares, and a temporary copy
  while it is written; a result is about 8.3 GB before compression. The
  central controller asks 300 GiB free and the ensemble controller 600 GiB.
  The 108 results are under 1 TB, so both campaigns can run in the
  project root.
* **Pilot.** Run 2 on the basalt grid (every square treated) alone, with
  192 G and 32 workers. After a few hours `sstat` gives its MaxRSS and the
  worker log its pace (consecutive checkpoint lines are two hours apart).
  `TC_FIELD_MEM` is then that MaxRSS plus the daily arrays still to be
  touched, with 15 % to spare, and `TC_FIELD_RESTARTS` the pieces 12,784
  days need at that pace, plus one. The other runs of the pilot's
  campaign go in beside it (`TC_FIELD_STAGE=resume`: the controller
  refuses only queued elements of the runs it is asked for), and the
  other grids, scenarios and cases as fresh calls.

## 4. The uncertainty summary

`Summarize_Uncertainty.m` (package root, with `summarize_uncertainty.sbatch`)
gathers the three sets of an arm into one summary file: per run, the
field means on the daily axis (Ca, Mg and Si in solution, cumulative
gaseous N loss, WFPS, drainage, cumulative N2O and NO3 leaching, the NO3
pool, bulk and pore-water pH, NEE, daily N2O, cumulative feedstock Ca and
Mg release), the yearly statistics (May-September bulk pH, calendar-year
N2O, NO3 leaching and drainage, year-end release, potential CDR), the
envelope of the three sets around the original, the basalt-minus-control
pH signal, the observations, and a yearly and a ranking CSV. The file
follows the contract `Plot_Uncertainty` and `Plot_1D_Then_All_12` read.

It reads a multi-square result of y through what y keeps: water content
and drainage from the daily stores `O_d` and `Lk_d`, feedstock release
from `ERW_daily` (each ERW day's mineral dissolution with the
stoichiometry of `ERW{c}.const`), fluxes and pH from `Domain`, pools from
`P` and NEE from the crown-weighted NPP and respiration stores, every
per-square quantity averaged with Domain's weights (the polygon overlap on
a field grid). `Final_state` holds only the last days and hours of what
the stores cover, and the summary needs none of it. A result with full
records (one column, or a package before y) is read from its hourly `O`
and `Lk` and ten-minute `weath.EW`, binned into the same daily rows; the
two give the same summary (`test_uncertainty_summary`). Rows are dated as
`Plot_1D_Then_All_12` dates them, and calendar years use complete days.
Every result is checked against its set (`Uncertainty_Result_Matches`),
and the three sets of a run against each other (squares, grid file,
start date). A run may come from either grid, so one summary holds a
layout's control, basalt and calcite runs.

It lives at the package root, outside every campaign's stamp, so it can
be corrected while campaigns run.

## 5. Checks

Run locally (R2026a, this PC; no model runs):

* `checkcode` on every new and edited MATLAB file: none of the messages
  `test_bouchet_linux_portability` rejects.
* New tests, both added to the audit list:
  `test_field_projection_campaigns` (the name bound on the Kantola grids'
  longest names, the grid helper, the ensemble's campaign definition, the
  tasks', batch files' and controllers' gates, the Linux checks of both
  directories) and `test_uncertainty_summary` (synthetic three-set results
  on two grids in both layouts: hand-computed field means, layouts equal,
  yearly statistics, envelopes, the provenance, missing-set, duplicate-run
  and time-axis checks, and `Plot_Uncertainty` on the output). Both pass.
* Tests that read the edited files pass: `test_uncertainty_run_options`,
  `test_field_switches`, `test_cec_pH_dependence`,
  `test_checkpoint_fingerprint_dependencies`, `test_uncertainty_campaign`,
  `test_ssp_controller_files`, `test_calcite_case_campaigns`,
  `test_bouchet_linux_portability`, `test_plot_uncertainty`.
* Both controllers were driven in a copy of the package with Slurm stubbed
  (submissions, overlap and lock refusals, every bad input, the job
  environment each element receives). This found one defect, fixed
  before the build: the ensemble's run list, which holds commas, was
  passed inside `--export`, whose commas separate variables, so every
  element would have read run 2 alone. It now reaches the jobs through the
  environment `--export=ALL` passes on.

The campaign side was written on y and moved onto z (`git rebase` of the
six commits onto z's 3436327: the driver and the audit list merged
without conflict; the README, RUNBOOK and pin lists keep both packages'
entries), then given the kernel gate; the checks above were run again on
the merged tree (`ERW_KERNELS_2026-10-08.md` for the kernels themselves).

On Bouchet, before production: z's checks (RUNBOOK section 0z: the MEX
build, the kernel comparison and the comparison with x under the fast
kernels, which replaces y's 0y), the equivalence baseline for this root,
then the pilot; every campaign's validation runs the complete audit with
the two new tests and z's kernel test, and one day of its grid under the
kernel it will run.

## 6. Files

New: `runtime/Field_SSP_2D/` (`submit_field_ssp_2d.sh`,
`submit_field_ssp_2d.sbatch`, `Run_Field_SSP_2D_Task.m`,
`Field_Campaign_Grid.m`, `Check_Lateral_Conservation.m` as a byte copy,
`README.md`), `runtime/Field_Uncertainty_2D/` (`submit_field_uncertainty.sh`,
`submit_field_uncertainty.sbatch`, `Run_Field_Uncertainty_Task.m`,
`Field_Uncertainty_Campaign.m`, `README.md`), `v2/run/Bounded_File_Tag.m`,
`Summarize_Uncertainty.m` and `summarize_uncertainty.sbatch` at the
package root, and the two tests.

Changed: `v2/run/Run_One_Subfield.m` (the bound, before both names are
formed), `v2/T&C_Code/Tests/run_post_legacy_equation_audit_tests.m` (the
two tests), `v2/plotting/Plot_Uncertainty.m` (titles from the summary's
label), both pin lists (the new campaign files, the tests and
`Bounded_File_Tag.m`), `runtime/Uncertainty_1D/README.md` (the summary),
this note, `README.md` and `RUNBOOK_2026-09-24.md` (section 0aa).

The kernel gate: both field controllers take `TC_ERW_KERNEL` (matlab, the
default, or fast; mex is refused before anything is submitted, since
MATLAB runs no MEX function on a thread pool) and pass it to every
element; the batch files check it again, and the tasks log it and check
that each result, and each validation smoke, records it in
`run_provenance.erw_kernel`. As z decided, the kernel is in neither the
validation stamp nor the checkpoint identity: the kernels compute the same
numbers, a run may resume under the other one, and
`run_provenance.execution_segments` keeps each segment's. `Run_Field_2D_Task`
(the five-year field runs) has no kernel gate yet.
