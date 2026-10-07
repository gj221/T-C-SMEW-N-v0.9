# Kantola five-year 1-D runs 10 to 12, CEC 13 set at alpha 3.4 (Bouchet, R2025a)

Campaign tag: `Base_1D`

This directory holds the controller for the whole calibrated all-12 campaign and the 1-D half of it (runs 10, 11 and 12: control, basalt and calcite on the single `vertical_1d` column of `Run_Matrix_12.m`). The companion 2-D half (runs 1 to 9) is `runtime/Base_2D`.

## Parameter set

Alpha 3.4 was adopted on 2026-10-02 (package u) in place of the 1.3 fitted on 26 September. It is the Appelo (1994) proton-loading slope, the value SCEPTER carries as its default, and the uncertainty ensemble runs 1.3 and 6 around it (`runtime/Uncertainty_1D`). KNmaxd 0.0048 and free drainage are kept from the joint sweep of 26 September (`Calibration_Sweep_1D`, arrays 27559013 to 27559017, root `tcsmewn_v09_kantola_postlegacy_20260925`; `POSTLEGACY_FIXES_2026-09-25.md` records that choice): KNmaxd 0.0048 reproduced the total chamber N2O at alpha 1.3 (0.97 of observed) and is carried to alpha 3.4 unchanged, so the first runs check it there, and free drainage is the boundary under which `Data_Kantola_BG.mat` was spun up (no new spin-up). In the sweep's case table the set is combination `a3p4_kn0p0048_kbfree_fd1` (combination id 65, sweep cases 193 to 195):

| parameter | value | how it reaches the model |
|---|---|---|
| Appelo alpha | 3.4, formed as `linspace(1.3,3.4,5)(5)` exactly as the sweep forms it (hex `400b333333333333`) | `AppeloAlpha` option, with the matched `AppeloLogKH` vector (5.136, 4.036, 8.943, 9.259, 14.179 for Na, K, Ca, Mg, Al) |
| KNmaxd | 0.0048 d^-1 | `KNmaxd` option |
| Kbot | free drainage (NaN, `Leakage_Bottom`: Lk = Ko), as in the spin-up | `Kbot` option |
| CEC | 13 cmolc kg^-1 | pinned workbook `v2/T&C_Code/Inputs/ERW_MOD_PARAM_parameters.xlsx` (not an option) |
| exchange scheme | Appelo | pinned `v2/T&C_Code/Parameters/MOD_PARAM_Kantola_Maize_Soy.m` (not an option) |

`Campaign_Parameters.m` is the only file that spells these numbers; every wrapper, validator and assertion reads it. `Assert_Base_Provenance.m` checks the saved `run_provenance` of every result against it (alpha, logKH, KNmaxd, Kbot, CEC 13 in every column, Appelo, constant CEC, 1826 days, column count, worker count, numerical identity, validation digest, campaign tag). Because alpha is bit-identical to the sweep's, runs 10 to 12 produce the same result filenames as sweep cases 193 to 195 and the 1-D task compares their daily states with the sweep results under this root (`sweep_identity_runNN.txt` in the logs; informational). No sweep has run in a package u root, so a fresh root reports the sweep result as not present. The sweep of package t ran the same options with the paste closed to the air, so its bulk pH, and everything that responds to it, differs from this package's.

## Files

| file | role |
|---|---|
| `submit_all_12.sh` | controller: preflights, submits both validations and both productions, held and released together |
| `submit_base_1d.sbatch` | 1-D batch file: validation mode or production of runs 10 to 12 as three single-CPU steps on one node |
| `Run_Base_1D_Task.m` | one production run: gates, lock, existing-result check, `Run_Subfield_Worker` with the calibrated options, post-run checks, provenance assertion |
| `Validate_Base_1D.m` | validation: full audit suite, Linux portability of both campaign dirs, configuration test, pre-application identity with the calibrated hydrology, one- and twelve-column smoke; writes `validation_stamp.txt` |
| `Validate_Preapplication_Identity.m` | four 315-day runs (control, basalt, calcite twice) must be exact before the first application and repeatable on one node |
| `Validate_All12_Smoke.m` | three-day run 10 (one column) and run 1 (twelve columns) with provenance and lateral-conservation checks |
| `test_all_12_configuration.m` | pins the parameter identity with the sweep, workbook CEC 13, Appelo scheme, wrapper options, Slurm resources, controller graph and README sentinels |
| `Fingerprint_Digest.m` | byte copy of the sweep's digest helper |
| `Calcite_Case.m` | reads `TC_CALCITE_CASE` for the tasks and the validator: the workbook row and the campaign-tag suffix |
| `phobs.m` | August bulk pH of runs 10 and 11 against both observed depth intervals, with the daily range as the sampling-day envelope |
| `budget.m` | alkalinity budget of runs 10 and 11: release, exchanger state, pore water, Ca+Mg leached and taken up, per-mineral dissolution, against the digitised Fig. 1F base saturation |
| `shared_pins.sha256` | SHA-256 of the shared code, workbook, parameter script and old campaign files the audit suite reads; the controller refuses any other checkout state |

## Submission

From the checkout root, after `deploy.sh` has installed both directories:

```bash
MATLAB_MODULE=MATLAB/2025a TC_BOUCHET_ACCOUNT=pi_np363 \
bash runtime/Base_1D/submit_all_12.sh
```

`TC_CALCITE_CASE=8` in front of the same command submits the calcite treatments alone (run 12 here, runs 3, 6 and 9 in `Base_2D`) with ERW workbook case 8, the agronomic 5 t/ha rate, as a campaign of their own: tags `Base_1D_case8` and `Base_2D_case8`, job names ending `_case8`, their own validation stamps and result trees, and result names carrying `__case8`. It runs beside the standard campaign.

`TC_ALL12_DIMS=1d` in front of the same command submits the 1-D half alone: the 1-D validation, then runs 10 to 12 (run 12 for a calcite case). The 2-D half follows once those results have been looked at, from the same command with `TC_ALL12_STAGE=resume` and no `TC_ALL12_DIMS`: both validations run again, the finished 1-D runs are kept, and runs 1 to 9 start.

Do not call this batch file directly with `sbatch`; it refuses without the controller token. The controller never rebuilds the equivalence baseline `v2/SIM/EQUIV_1D_TMP/TC_ERW_results_Kantola_crop_control.mat` (it is hashed into the finished CEC 13 sweep's stamp) and refuses if it is missing.

Job graph and resources:

| job | name | resources | dependency |
|---|---|---|---|
| 1-D validation | `Base1Dcheck` | day, 1 task, 6 CPUs, 40 GB, 4 h (audit ~10 min, identity ~30 min on a 4-process pool, smoke ~5 min) | none |
| 2-D validation | `Base2Dcheck` | day, array element 1, 6 CPUs, 32 GB, 4 h (configuration and portability tests) | none |
| 1-D production | `Base1D` | day, one node, 3 tasks x 1 CPU, 24 GB, 6 h (each run took 2 h 14 min and 3.1 GB in the sweep) | both validations |
| 2-D production | `Base2D` | week, array 1-9, 6 CPUs, 48 GB, 4 days per element | both validations |

Stages: `TC_ALL12_STAGE=fresh` (default) requires empty result and log trees for both campaign tags; `TC_ALL12_STAGE=resume` resubmits the same graph after a TIMEOUT or a cancellation: validation rewrites an identical stamp when nothing changed, finished runs are kept (`Bouchet_Check_Existing_Result`), interrupted runs continue from their checkpoints. A node failure needs no resume: both batch files carry `--requeue`, the requeued element keeps its job id, reclaims its own lock and resumes its checkpoint. Both stages refuse while any campaign job is queued or a `locks/run_NN.lock` remains (the controller prints the `rm -r` line for each stale lock).

Sequencing with the CEC 13 sweep: the sweep's rescore packages replace `v2/T&C_Code/Tests/test_bouchet_1d_joint_sweep_configuration.m`, which this campaign's stamp also covers. Deploy such a package before this campaign is submitted or after all twelve runs have finished, including any resume: a resume revalidates, and with a changed stamp every finished run is refused as `Bouchet:StaleResult` and the campaign has to start again in fresh trees. `deploy.sh` refuses while any job runs from the checkout.

## Outputs

* Results: `v2/SIM/Subfield_12/<tag>/TC_2D_runNN_<treatment>__vertical_1d__alpha3p4__logkh..._14p178822881__knmaxd0p0048__kbotfree__1826d__from20160101.mat` with `<treatment>` control, basalt or calcite for runs 10, 11, 12 (about 250 MB each); `checkpoints/` is emptied on success. The logKH values of the seed-anchored rule print at 12 digits, so the longest campaign filename (2-D margin_prairie) is 178 characters against the filesystem's 255 limit.
* Logs, status, locks and stamp: `v2/logs/Subfield_12/<tag>/` (`validation_stamp.txt`, `run_NN.status`, `worker_runNN.log`, `srun_runNN_<job>.out|err`, `sweep_identity_runNN.txt`).
* Slurm output: `runtime/<tag>/slurm_<job>.out|err` (not fingerprinted).

Post-processing: `Plot_1D_Then_All_12` and `pull_bouchet_campaign.sh` expect result names of the form `<treatment>__<layout>__1826d`; copy the MATs and rename them by keeping the first two `__` fields, for example `awk -F'__' '{print $1"__"$2"__1826d.mat"}'`. Do not edit the plotting code (it is fingerprinted); the pull script is not fingerprinted, but keep it unchanged as well.
