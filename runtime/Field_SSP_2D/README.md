# Kantola 35-year SSP runs on a field grid (Bouchet, R2025a)

Campaign tag: `Field_SSP_2D_35yr_<ssp>[_case8]__<grid name>`

Runs 1 to 9 of `Run_Matrix_12.m` (control, basalt and calcite on the `whole_field`, `margin_crop` and `margin_prairie` layouts) on a field grid (`runtime/Field_2D`, the 10 m Kantola grids), over the 35 projection years: the scenario forcing `Data_Kantola_<ssp>.mat`, 2016-01-01 to 2050-12-31 (12784 days), rock every 10 November to 2050, as `runtime/SSP_2D` runs them on the 12-column plot. `TC_CALCITE_CASE=8` is the 5 t/ha calcite arm (runs 3, 6 and 9). Everything else is the calibrated set of `runtime/Base_1D` (`Campaign_Parameters.m`, `Assert_Base_Provenance.m`) with the horizon moved to 35 years, and the run options are `SSP_2D`'s with the field grid in place of the 2x2 plot grid. A result here is the original member of the field-grid uncertainty ensemble (`runtime/Field_Uncertainty_2D`).

The model is package y's (`FIELD_RECORDS_DAILY_2026-10-08.md`): the squares' daily vegetation and biogeochemistry, hourly physics and daily ERW chemistry share `TC_FIELD_WORKERS` workers, and a result keeps per-square daily stores, hourly field means and the final state instead of every hourly and ten-minute row.

## The grid

As for `runtime/Field_2D`: a MAT file holding `topo`, saved by `Build_Field_Grid_Topology`, **beside** the deployed root (a validated root never gains files), named by at most 32 letters, digits and underscores; the name enters the campaign tag and the content the validation stamp, so a grid replaced after validation stops production with `StaleValidation`. The Kantola grids are `kantola_basalt_dem_10m.mat` (340 squares) and `kantola_control_dem_10m.mat` (369 squares) from `Downloads/kantola_field_grids_10m` (`SHA256SUMS` beside them).

One call submits the runs of one grid, scenario and case. Which runs go on which grid is the call's `TC_FIELD_RUNS`: the RUNBOOK's commands put the controls (1, 4, 7) on the control field and the basalt and calcite runs (2, 3, 5, 6, 8, 9) on the basalt field, 36 runs over the three scenarios and the 5 t/ha arm.

## Resources

* **Memory.** Package y holds every daily array of every square whole from the start: 672 values a square and day, the 12 hourly-derived daily stores and the 79 values of an ERW day's summary, 6.1 kB a square and day. Over 12,785 daily rows that is 26.8 GiB on 369 squares and 24.7 GiB on 340, all resident by the end of the run. Add MATLAB, the thread pool and the per-square state (about 4 GiB, an estimate): the task refuses a job whose `TC_FIELD_MEM` is below that floor (`Field_Campaign_Grid('memory',...)`: 30.8 GiB on 369 squares, 28.7 GiB on 340). The first ERW day also builds each square's whole-run ten-minute timetable on its worker, a passing peak of about `TC_FIELD_WORKERS` timetables that the pilot measures.
* **Time.** No 35-year field-grid run of y has been timed. Package x's cost model (2.1 h of ERW chemistry and 0.16 h of hourly physics per square per five years) gives about 4 days for 369 squares on 64 workers and 8 days on 32 if y's three parallel loops scale; the field rows of `compare_x_y.sbatch` (RUNBOOK section 0y) and the pilot replace it with measurements.
* **Partition.** The `day` partition (default): per-user limits of 1500 CPUs and 20 TB, pieces of 24 hours. A run requeues itself half an hour before the limit (`runtime/Uncertainty_1D/requeue_near_limit.sh`), at most `TC_FIELD_RESTARTS` times (default 12, about 12.5 days of computing), and continues from its two-hourly checkpoint under the same job id and run lock; each restart repeats up to two hours. The `week` partition (`TC_FIELD_PARTITION=week`, pieces of up to 7 days) holds 96 CPUs per user, one 64-worker run at a time.
* **Checkpoints** every two hours (`TC_CKPT_MINUTES=120`): the daily stores, ERW day summaries and state, about 9 GB on 369 squares, written uncompressed, plus a temporary copy while each is written. A result is about 8.3 GB before compression.
* **Pilot first.** Submit run 2 alone (basalt on the whole field: every square treated, the heaviest) with `TC_FIELD_RUNS=2 TC_FIELD_WORKERS=32 TC_FIELD_MEM=192G`. After its first hours, `sstat -a -j <job>_2 --format=JobID,MaxRSS,AveCPU` gives the peak so far (the first ERW day's timetables included) and the worker log's checkpoint lines the pace: consecutive `[checkpoint] iter` lines are two hours apart, so their difference in iterations is simulated hours per two wall hours. Then size the other runs: `TC_FIELD_MEM` = the pilot's MaxRSS plus the floor's daily arrays not yet touched (up to 27 GiB on 369 squares), rounded up with 15 % to spare; `TC_FIELD_RESTARTS` = 12,784 days divided by the days one 23.5-hour piece simulates, rounded up, plus one. Submit them while the pilot runs (`TC_FIELD_STAGE=resume TC_FIELD_RUNS=1,3-9`, or the other grids, scenarios and cases as fresh calls).

## Submission

From the checkout root, after the equivalence baseline has been built for this deployment and the scenario forcing copied in:

```bash
TC_FIELD_TOPOLOGY=/full/path/beside/the/root/kantola_basalt_dem_10m.mat \
TC_SSP=ssp245 TC_FIELD_RUNS=2 TC_FIELD_WORKERS=32 TC_FIELD_MEM=192G \
MATLAB_MODULE=MATLAB/2025a TC_BOUCHET_ACCOUNT=<pi account> \
bash runtime/Field_SSP_2D/submit_field_ssp_2d.sh
```

* `TC_FIELD_RUNS` (default `1-9`; `3,6,9` with a case): the Slurm array of runs, such as `2`, `1,4,7`, or `2-9%4` to hold four at a time.
* `TC_FIELD_WORKERS` (default 32, 2 to 64): the pool and the CPUs of every element. It enters the stamp, so every call of one campaign must keep it.
* `TC_FIELD_MEM` (no default), `TC_FIELD_PARTITION` (`day`), `TC_FIELD_TIME` (the partition's limit), `TC_FIELD_RESTARTS` (12).
* `TC_FIELD_STAGE`: `fresh` (default) needs empty result and log trees; `resume` adds runs or continues them. Either refuses while a validation of the campaign is queued, or an element of a requested run, or while a requested run holds a lock (it prints the `rm` to run when no job holds it); elements of other runs keep running.

The controller submits one validation element and the production array.
* **Validation:** `day`, `TC_FIELD_WORKERS` CPUs, 96 GB, 6 h, never requeued. It runs the complete audit, the Linux portability checks and one day of the grid across the first application (run 5, or 6 with a case) under the campaign's scenario, case and calibrated options, saved and read back (grid provenance, lateral conservation), prints the memory floor and writes the stamp.
* **Production:** the runs of `TC_FIELD_RUNS`, eligible after the validation passes. Finished runs are kept; a run whose result already carries the current stamp exits at once.

## Files

| file | role |
|---|---|
| `submit_field_ssp_2d.sh` | controller: preflight, overlap checks, validation element, production array |
| `submit_field_ssp_2d.sbatch` | array batch file: validation mode, and production with the requeue |
| `Run_Field_SSP_2D_Task.m` | one element: gates, memory floor, stamp, lock, existing-result check, `Run_Subfield_Worker` with `FieldTopologyFile` and `Scenario`, the gate recomputed, provenance (calibrated set, grid file and its SHA-256, scenario, case, last application year), lateral conservation |
| `Field_Campaign_Grid.m` | the grid file's checks and provenance, run lists, the memory floor; shared with `runtime/Field_Uncertainty_2D` |
| `Check_Lateral_Conservation.m` | byte copy of the 2-D campaigns' whole-run check |

## Outputs

* Results: `v2/SIM/Subfield_12/<tag>/TC_2D_runNN_<treatment>__<layout>__<ssp>[__case8]__field_<shapefile>__cell10m__margin20__dem10__alpha3p4__logkh...__knmaxd0p0048__kbotfree__12784d__from20160101.mat` (at most 238 characters on the Kantola grids).
* Logs, status, locks and stamp: `v2/logs/Subfield_12/<tag>/`; the worker log of a run gathers every piece.
* Slurm output: `runtime/Field_SSP_2D/slurm_<array job>_<element>.out|err`, appended by each piece.

The field maps (`runtime/Field_2D/README.md`, Figures) read these results from their daily stores. The ensemble's summary is `Summarize_Uncertainty.m` at the package root (`runtime/Field_Uncertainty_2D/README.md`).
