# Kantola uncertainty ensemble on a field grid (Bouchet, R2025a)

Campaign tag: `Field_Uncertainty_2D_35yr_<ssp>[_case8]__<grid name>`, with `lower/` and `upper/` under it

The lower and upper parameter sets of the Kantola ensemble (`runtime/Uncertainty_1D/README.md`: all sixteen parameters at one end of their ranges, alpha 1.3 and 6 around the campaign's 3.4) on the runs of the 35-year SSP campaign on a field grid (`runtime/Field_SSP_2D`). The original set is that campaign's own result and is not run here. Each set runs as the original's campaign runs: the scenario forcing, 12784 days, rock to 2050, the field grid, `TC_FIELD_WORKERS` workers, the 5 t/ha calcite arm with `TC_CALCITE_CASE=8`. The members, their knobs and the provenance a result must carry are the plot ensemble's own (`Uncertainty_Members.m`, `Uncertainty_Knobs.m`, `Uncertainty_Result_Matches.m`, `Uncertainty_Result_State.m`), so a field set differs from the plot's only by the grid.

## Names

Every field-grid set's result and checkpoint name is 258 to 272 characters as `Run_One_Subfield` strings its tokens together (run, scenario, case, the grid's tag, then alpha, five exchange constants, KNmaxd, the overrides, the dissolution factor, the saturation exponent, Kbot, the span), past the 255 a Linux file name may have. `Bounded_File_Tag` (`v2/run`) therefore replaces the knob tokens of such a name by `__knobs` and twelve hex digits of their SHA-256; the run, scenario, case, grid and span stay readable, and `run_provenance` carries every knob. Any name that fits is left as it was, so no earlier result or campaign name changes. The validation element writes one day of each set under such a name and reads it back.

## The model check

The three sets compare only if they ran the same model, so the controller compares this checkout with the original's (`runtime/Uncertainty_1D/same_model_code.sh`: `v2/T&C_Code` without its tests, its inputs, `v2/run` and the calibrated set, byte for byte) before it submits anything. With both campaigns in one root, `TC_UC_ORIGINAL_ROOT="$PWD"`.

## Submission

From the checkout root, after the field SSP campaign of the same grid, scenario and case has been sized by its pilot:

```bash
TC_UC_ORIGINAL_ROOT="$PWD" TC_ERW_KERNEL=fast \
TC_FIELD_TOPOLOGY=/full/path/beside/the/root/kantola_basalt_dem_10m.mat \
TC_SSP=ssp245 TC_FIELD_RUNS=2,3,5,6,8,9 TC_FIELD_WORKERS=32 TC_FIELD_MEM=<from the pilot> \
MATLAB_MODULE=MATLAB/2025a TC_BOUCHET_ACCOUNT=<pi account> \
bash runtime/Field_Uncertainty_2D/submit_field_uncertainty.sh
```

* `TC_FIELD_RUNS` (default every run of the arm): the runs, as the field SSP controller takes them; a `%N` suffix holds N elements at once. The array is member-major: elements 1 to n are the lower set over the runs, n+1 to 2n the upper.
* `TC_FIELD_WORKERS`, `TC_FIELD_MEM`, `TC_FIELD_PARTITION`, `TC_FIELD_TIME`, `TC_FIELD_RESTARTS` and `TC_ERW_KERNEL` as for `runtime/Field_SSP_2D` (README, Resources); the memory and pace of the original's pilot hold for both sets, which move parameters, not arrays.
* The controller refuses while any job of the call's campaign is queued and lists stale run locks to remove. After a time-out or node loss, resubmit the same command once its jobs have left the queue: finished elements exit at once (after rereading their grid and conservation), the rest continue from their checkpoints, and a result of anything else is refused rather than overwritten.

The controller submits one validation element (`day`, `TC_FIELD_WORKERS` CPUs, 96 GB, 6 h, never requeued: the complete audit, the Linux checks of this directory and the code it borrows, and one day of the grid under each set, saved and read back) and the production array, eligible after it passes. Each production element requeues itself before the `day` partition's 24-hour limit, at most `TC_FIELD_RESTARTS` times, and continues from its two-hourly checkpoint under the same job id and run lock.

## Files

| file | role |
|---|---|
| `submit_field_uncertainty.sh` | controller: preflight, the model check, validation element, production array |
| `submit_field_uncertainty.sbatch` | array batch file: validation mode, and production with the requeue |
| `Run_Field_Uncertainty_Task.m` | one element: gates, memory floor, stamp, lock, the set's existing result, `Run_Subfield_Worker` with the set's knobs, the gate recomputed, the set's and the grid's provenance, lateral conservation |
| `Field_Uncertainty_Campaign.m` | the campaign from the environment: tags, runs, options, workers, dependencies, elements |

## Outputs and the summary

* Results: `v2/SIM/Subfield_12/<tag>/<lower|upper>/TC_2D_runNN_<treatment>__<layout>__<ssp>[__case8]__field_<shapefile>__cell10m__margin20__dem10__knobs<hex>__12784d__from20160101.mat`.
* Logs, status and locks: `v2/logs/Subfield_12/<tag>/<lower|upper>/`; the stamp: `v2/logs/Subfield_12/<tag>/validation_stamp.txt`.
* Slurm output: `runtime/Field_Uncertainty_2D/slurm_<array job>_<element>.out|err`.

`Summarize_Uncertainty.m` at the package root gathers the three sets of an arm into one summary (field means on the daily axis, yearly statistics, the envelope of the three sets, the observations; `v2/plotting/Plot_Uncertainty.m` draws it). It reads the daily stores, the ERW day summaries and Domain of these results, and takes the runs of one arm from both grids at once:

```bash
o=v2/SIM/Subfield_12; sbatch --account=pi_np363 --export=ALL,UC_SOURCES="$o/Field_SSP_2D_35yr_ssp245__kantola_control_dem_10m|$o/Field_Uncertainty_2D_35yr_ssp245__kantola_control_dem_10m;$o/Field_SSP_2D_35yr_ssp245__kantola_basalt_dem_10m|$o/Field_Uncertainty_2D_35yr_ssp245__kantola_basalt_dem_10m",UC_RUNS="[1 2 3]",UC_OUTPUT=$o/uncertainty_ssp245_whole_field.mat summarize_uncertainty.sbatch
```
