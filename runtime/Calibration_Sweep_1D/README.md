# Kantola 1-D alpha, KNmaxd and Kbot sweep on Bouchet

This sweep covers 112 parameter combinations across runs 10 to 12. Each trajectory runs 1,826 days from 1 January 2016. Its daily rows are the lone first hour, dated 31 December 2015, then 1,825 complete days, so the reported 2020 period is 1 January through 29 December. The forcing does not contain the last six hours needed for a complete 1,827-day window.

| Parameter | Values | Code units |
|---|---|---|
| fixed CEC | 13 | cmolc kg^-1 |
| Appelo alpha | 1.3, 1.825, 2.35, 2.875, 3.4, 4.5, 6 | dimensionless |
| KNmaxd | 0.0016, 0.0032, 0.0048, 0.0064 | d^-1 |
| Kbot | free drainage, 0.07, 0.03, 0.01 | mm h^-1 |
| dissolution factor (DissolutionFactor, workbook diss_f) | 1 (held) | dimensionless |
| run | 10, 11, 12 | control, basalt, calcite |

The fixed value is the best site-specific untreated-control proxy. Pihlap et al. (2025, doi:10.1016/j.isci.2025.114232) reported 13 +/- 2 cmolc kg^-1 for Energy Farm control plots sampled in October 2022 at 1 to 6 and 15 to 20 cm using the A&L/NCERA-13 soil analysis. Basalt and lime plots were 14 +/- 1 cmolc kg^-1 and were not significantly different from control, so the same capacity is used for every treatment. The earlier value of 8.6 was prescribed rather than measured. The 19 to 20 cmolc kg^-1 CEC7 evidence is retained as an upper structural sensitivity.

SCEPTER assigns a fixed mass-specific CEC to each exchanger solid, although aggregate cell CEC can change if solid masses change. This lumped TC Appelo treatment has no evolving exchanger-solid inventory, so soil CEC stays fixed at 13 cmolc kg^-1 and the pH response comes from H-site-dependent apparent selectivity. A separate positive linear `dCEC_dpH` is therefore rejected in Appelo mode. The seven alpha values (the five earlier ones, then 4.5 and 6, added when the budget of package m showed the exchanger saturating where the field's did not) use matched code-input `logKH` values:

```text
logKH_i(alpha) = logKH_i(alpha=1.3) + z_i*(alpha - 1.3)*0.20
```

Here `z_i = [1, 1, 2, 2, 3]` for Na, K, Ca, Mg and Al. SCEPTER accepts alpha and `logKH` independently. Selectivity is `10^(-logKH + alpha*z*f_H)`, so this coupling holds it fixed at the seed's `f_H = 0.20` and every alpha starts from the measured 80 % base saturation; alpha then only sets how the exchanger responds once pH moves. The earlier coupling, `z_i*(alpha - 1.3)`, held selectivity fixed at `f_H = 1` and let the seed fall to 69 % base saturation at alpha 2.35, so results from it are not comparable.

The pore-to-bulk and bulk-to-pore operators use the same charge-conserving equilibria of a paste open to the laboratory air (`PASTE_OPEN_TO_AIR_2026-10-02.md`), so their inverse/forward round trip is internally consistent. SCEPTER 1.0.2 instead runs a separate 25 C laboratory-extraction simulation with transferred solid and exchanger inputs and its configured aqueous activity model. TC does not reproduce that full extraction run or Davies activity coefficients. This is a reduced SCEPTER-like Appelo exchange treatment.

The legacy moisture response remains `fmoist = Se^(1/0.5) = Se^2`. `KNmaxd` multiplies the remaining dimensionless denitrification response terms. A blank `Kbot`, represented by `NaN`, retains free drainage.

| Variable | Expected SI Units | Actual Units provided in the code |
|---|---|---|
| CEC | molc kg^-1 | cmolc kg^-1 |
| Appelo alpha | 1 | 1 |
| Appelo fractions | 1 | 1 |
| KNmaxd | s^-1 | d^-1 |
| Kbot | m s^-1 | mm h^-1 |
| GPP, NEE and Re | kg C m^-2 s^-1 | g C m^-2 d^-1 |
| cumulative N2O and NO3 loss | kg N m^-2 | g N m^-2 |

The conversions are `CEC_SI = CEC*0.01`, `KNmaxd_SI = KNmaxd/86400`, `Kbot_SI = Kbot*1e-3/3600`, and `Cflux_SI = Cflux*1e-3/86400`.

## Runtime plan

One completed trajectory from the earlier sweep took 2:29:23, used about 95% of one CPU, and reached 3.12 GiB MaxRSS. More CPUs inside one trajectory would sit idle because this is a one-column, single-thread calculation. More RAM would not make it faster either.

The default submission is a one-year control gate, not the 336-case sweep. It runs case 1 for all 366 days of 2016 using alpha 1.3, KNmaxd 0.0048 d^-1, free drainage and fixed CEC 13 cmolc kg^-1. Scaling the measured five-year runtime gives about 30 minutes for this trajectory after validation and queueing. The result is stored under `pre_control_gate`, separate from five-year task 1.

Only 1 January through 9 November 2016 is compared with the manuscript `pre` control observation. This gives 314 daily states. The model and the two observation depths are each averaged in ideal hydrogen-concentration space. The output reports the signed model-minus-observation discrepancy and its absolute value. With one comparison, bias equals the signed discrepancy and RMSE equals the absolute error. This is a drift and stability check because the initial pH was conditioned to the pretreatment soil. It is not a multi-year fit or a test of the basalt response. N2O, NO3, GPP, NEE and Re observations begin in 2017, so they are not scored by this gate.

The controller uses one CPU and 8 GiB per trajectory. It allows up to 48 concurrent trajectories across available nodes, reducing the ideal production time from about 25 hours at 24 slots to about 12.5 hours at 48 slots. Queueing still applies. Set `TC_SWEEP_MAX_CONCURRENT` to an integer from 1 through 48 to use a lower cap.

The default `independent` layout has 336 array elements. `TC_SWEEP_LAYOUT=packed` has 112 array elements. Each packed allocation reserves three CPUs and 24 GiB, then starts the control, basalt and calcite cases as three separate single-thread MATLAB processes. Every child step has an explicit 8 GiB limit, its own CPU, scratch directory, lock, checkpoint, status and result. The packed pilot first proves that three exclusive Slurm steps overlap on distinct CPU masks. A failed process makes the packed element fail, while completed siblings remain valid and are reused after that element is resubmitted. Packed mode reduces scheduler entries, but it does not reduce core-hours and may be harder to backfill. Use it when Bouchet is limiting active array elements rather than available cores.

Before the production array, a full extreme-case pilot runs at alpha 6, KNmaxd 0.0064 d^-1 and Kbot 0.01 mm h^-1. In independent mode this is task 336, the calcite case, with one CPU and 8 GiB. In packed mode this is array element 112, which runs tasks 334 to 336 with three CPUs and 24 GiB. Production has an `afterok` dependency on the pilot and reuses its validated results.

No numerical solver was altered for speed. Shorter root scans look promising, but they require MATLAB profiler evidence and trajectory-equivalence tests first.

## Submission graph

With no stage option, the controller submits three or four jobs in order:

1. Only when `v2/SIM/EQUIV_1D_TMP/TC_ERW_results_Kantola_crop_control.mat` is missing, or `TC_EQUIVALENCE_BASELINE=rebuild` is given: a fresh R2025a one-dimensional equivalence baseline.
2. The full post-legacy audit and treatment smoke runs, which write the validation stamp.
3. One 366-day control trajectory.
4. The checked manuscript pre-pH comparison, followed by a stop.

The controller prints `Equivalence baseline: reused ...` or the baseline job ID. The baseline MAT file is a dependency of the validation fingerprint and, like every MAT file, embeds its creation timestamp, so rebuilding it changes the fingerprint even when its numbers are identical. Every result produced under the previous fingerprint then reads as stale: `Bouchet_Check_Existing_Result` refuses it and the pre-control recheck fails. Rebuild only when the baseline is missing or the source changed, and only through the exploratory stage.

The fingerprint (`Bouchet_Campaign_Validation_Fingerprint` as installed on Bouchet, read back from the archived stamp) hashes exactly these files: `v2/T&C_Code/*.m`, `v2/T&C_Code/Parameters/*.m`, `v2/T&C_Code/Tests/*.m`, `v2/T&C_Code/Tests/Fixtures/*.mat`, the six inputs `Data_Kantola.mat`, `Data_Kantola_BG.mat`, `Ca_Data.mat`, `All_deposition_data.mat`, `ERW_MOD_PARAM_parameters.xlsx` and `blueridge_CDF_data.json` under `v2/T&C_Code/Inputs`, `v2/run/*.m`, `v2/audit/*.m`, `v2/plotting/*.m`, the equivalence baseline `v2/SIM/EQUIV_1D_TMP/TC_ERW_results_Kantola_crop_control.mat`, the campaign `*.m`, `*.sbatch`, `*.sh` and this README, and the four observation files under `v2/Experimental observations/Web plot digitizer`. Plotting and audit scripts and the controller shell script are all included. Installing a package that touches any of them re-stamps the campaign on its next validation, and none of them may be edited, added or removed between the exploratory run and the end of production (restoring the exact previous bytes makes the stamp match again). The controller reads the stamp on disk, re-hashes every file it lists, scans the globs above for files the stamp does not list, and refuses to submit when a fingerprinted file changed, disappeared or appeared, when a one-year result exists that the next validation would not accept, or when a `run_NN.lock` from a killed job remains under the campaign logs (it prints the exact `rm -r` command; it has already confirmed with `squeue` that no sweep job is active).

No pilot or production array is submitted in the default stage. After reviewing the one-row CSV, the full stage keeps the equivalence baseline (it refuses `TC_EQUIVALENCE_BASELINE=rebuild`), checks that the reviewed one-year result carries the SHA-256 digest of the validation stamp on disk, reruns validation, rechecks the exploratory result against that stamp in MATLAB, then submits the extreme pilot, production array and full summary. The full stage requires the exact acknowledgement `TC_EXPLORATORY_REVIEWED=YES`.

### Rescoring a complete sweep after a scoring fix

Scoring-side files (`Kantola_N_Magnitude_Audit.m`, `Summarize_Sweep.m`, `Validate_Sweep.m`, `Assert_Sweep_Provenance.m`, `Accepted_Stamp_Digests.m`, the controller `submit_sweep.sh`, this README and the two test files that cover them) are fingerprinted, so a fix changes the validation stamp although no trajectory could be affected. The `rescore` stage scores the 336 finished results under the corrected code without recomputing anything:

```bash
MATLAB_MODULE=MATLAB/2025a \
TC_BOUCHET_ACCOUNT=pi_np363 \
TC_SWEEP_STAGE=rescore \
TC_ACCEPT_PRIOR_STAMP=YES \
bash runtime/Calibration_Sweep_1D/submit_sweep.sh
```

The controller refuses unless all 336 results and complete status files exist and every file that differs from the stamp on disk is in the scoring-side list above (nine fingerprinted files; both the controller and `Accepted_Stamp_Digests` carry that list, and the configuration test keeps them identical). Every stage archives the stamp it reads as `validation_stamps/<digest>.txt` before any job can overwrite it, and validation archives every stamp it writes. The rescore submits validation, which writes the new stamp, then the summary with `TC_ACCEPTED_RESULT_DIGESTS` naming every archived stamp other than the one it validates against (`:`-separated, because `sbatch --export` splits on commas). The summary accepts a result carrying a prior digest only after re-deriving, from the archived stamp, that the prior and current stamps share the horizon and numerical runtime identity and differ solely in scoring-side files; an archived stamp that fails that test is merely not accepted, and a result carrying it is refused with `Bouchet_Joint_Sweep:ResultStampNotAccepted` and the reason. Every trajectory row records `scored_under_validation_fingerprint` and `ERW_steps_per_day` (the weathering module retains its closure state per step, 144 per day in this campaign; the closure tolerances are applied to every retained step), and `joint_sweep_summary.status` records `validation_fingerprint` and `accepted_prior_digests`. When nothing changed since the stamp, the controller submits the summary alone, still offering the archived stamps, so rerunning the rescore after a failed or interrupted summary is safe.

Lineage acceptance applies to the full summary only. After a rescore the stamp on disk is new, so the exploratory and full stages refuse the earlier one-year gate and the 336 results as stale and, when finished production results exist, say so and point back to the rescore; a sweep that needs any trajectory recomputed is a new campaign.

### Recovering a stale one-year gate

If either stage refuses with `different validation stamp`, `fingerprinted file changed` or `would fail with Bouchet:StaleResult`, or a control task fails with `Bouchet:StaleResult`, the one-year result no longer matches the stamp the next validation will write. Move it aside and rerun the exploratory stage from the repository root:

```bash
stale=$(date -u +%Y%m%dT%H%M%SZ)
mkdir -p v2/SIM/Subfield_12/stale_gates v2/logs/Subfield_12/stale_gates
mv v2/SIM/Subfield_12/Calibration_Sweep_1D/pre_control_gate "v2/SIM/Subfield_12/stale_gates/pre_control_gate_$stale"
mv v2/logs/Subfield_12/Calibration_Sweep_1D/pre_control_gate "v2/logs/Subfield_12/stale_gates/pre_control_gate_$stale"
```

The moved directories are never scanned by the campaign. The model is deterministic on one thread and one CPU model, so the summary log line `Model minus observation` is expected to repeat the reviewed value, and the full-precision `pH_residual_model_minus_observed` in the new CSV and status file is expected to be identical to the reviewed run's; a difference in the leading six decimals indicates a changed input or runtime and should be understood before the full stage is released.

Production retains MATLAB R2025a, `cpumodel:9575f`, one computational thread, per-case locks, rolling checkpoints, and the validation fingerprint. Result provenance must match CEC 13, zero `dCEC_dpH`, the selected alpha and its matched `logKH`, KNmaxd, Kbot, run ID, dates, runtime identity and array task.

## Storage

The exploratory stage requires 10 GiB free according to `df -Pk`. Final MAT files from the full stage are expected to use about 47,600 MiB, or 46.5 GiB. Forty-eight live tasks can also write checkpoints and temporary saves, so the full stage requires at least 85 GiB free.

```bash
R=/nfs/roberts/project/pi_np363/gj326/tcsmewn_v09_20260926j
df -Pk "$R"
quota -s
quota -gs pi_np363
```

Do not submit unless the filesystem and applicable quota meet the requirement for the selected stage.

## Submit

The controller refuses to start while any of the sweep's job names is active.

Run the exploratory gate first:

```bash
cd /nfs/roberts/project/pi_np363/gj326/tcsmewn_v09_20260926j

test -r 'runtime/Calibration_Sweep_1D/submit_sweep.sh'

MATLAB_MODULE=MATLAB/2025a \
TC_BOUCHET_ACCOUNT=pi_np363 \
TC_SWEEP_STAGE=exploratory \
bash runtime/Calibration_Sweep_1D/submit_sweep.sh
```

The command prints `Baseline: reused` and three job IDs (four when a baseline job was submitted). Acceptance requires `COMPLETED 0:0` for every printed job, the R2025a validation PASS sentinel, and this final summary sentinel:

```text
PASS: one-year control trajectory and manuscript pre-pH comparison
STOP: no pilot or production sweep was submitted
```

Review:

```text
v2/logs/Subfield_12/Calibration_Sweep_1D/pre_control_gate/pre_control_pH_discrepancy.csv
```

Only after that review, release the full independent sweep with:

```bash
MATLAB_MODULE=MATLAB/2025a \
TC_BOUCHET_ACCOUNT=pi_np363 \
TC_SWEEP_STAGE=full \
TC_EXPLORATORY_REVIEWED=YES \
TC_SWEEP_MAX_CONCURRENT=48 \
bash runtime/Calibration_Sweep_1D/submit_sweep.sh
```

Add `TC_SWEEP_LAYOUT=packed` to that full-stage command to use 112 three-process array elements.

Run this from the repository root. Do not invoke the sbatch file from inside the campaign directory.

The controller uses an atomic `.submission_controller.lock` directory at `v2/logs/Subfield_12/Calibration_Sweep_1D/.submission_controller.lock`. If a stopped controller leaves the empty lock behind, first confirm that no old or new sweep job is active, i.e. that this prints nothing:

```bash
squeue -h -u "$USER" --name=SweepBase,SweepCheck,SweepExplore,SweepExploreSum,SweepPilot,Sweep1D,SweepSummary -o '%A|%j|%T'
```

Only then remove that empty lock with `rmdir`.

## Check progress

For the exploratory stage, replace the names below with the printed IDs, omitting `BASELINE_ID` when the controller printed `Baseline: reused`:

```bash
sacct -j BASELINE_ID,VALIDATION_ID,EXPLORATORY_ID,EXPLORATORY_SUMMARY_ID -X \
  --format=JobID%24,JobName%22,State%12,ExitCode,Elapsed,MaxRSS,NodeList%20

squeue -j BASELINE_ID,VALIDATION_ID,EXPLORATORY_ID,EXPLORATORY_SUMMARY_ID \
  -o '%.18i %.22j %.10T %.10M %.10l %.35R'
```

The validation log must finish with `PASS: Bouchet R2025a joint 1-D parameter sweep validation`.

The exploratory inventory must contain one result and one `state=complete` status. After a separately released full production stage, the exact inventory must contain 336 results and 336 complete status files. The full summary compares pH, NO3, N2O, GPP, NEE and Re. It does not submit the final all-12 campaign.

### Scoring objectives

`joint_sweep_combination_summary.csv` holds one row per alpha, KNmaxd and Kbot tuple (runs 10, 11 and 12 are its control, basalt and calcite trajectories). `Pareto_nondominated` minimises five objectives at once:

1. `N2O_cumulative_log_ratio_RMS`: for control 2017, control 2018, basalt 2017 and basalt 2018, the observed cumulative emission is the trapezoid over the chamber dates of that year (the convention of the field studies) and the model's daily flux is integrated by the same rule over the same date window; the objective is the RMS of ln(model/obs) over the four, so 0 is exact and 0.69 is a factor of two. The model does not resolve the wetting pulses that dominate the chamber record, so a day-matched score measures pulse timing rather than emission. In the pinned observation file the campaigns differ in length: control 24 April to 20 September 2017 (14 dates) and 27 April to 31 August 2018 (17 dates), basalt 12 May to 21 June 2017 (9 dates) and 10 May to 27 June 2018 (10 dates). Each treatment-year is scored over its own window, so the basalt totals are May to June totals and are not comparable with the control totals in the same row; the windows and observed totals are carried beside the modelled totals and ratios (`N2O_cumulative_<treatment>_<year>_window`, `_observed_kgN_ha`, `_kgN_ha`, `_ratio`).
2. `pH_season_discrete_RMSE`: the arithmetic mean of daily bulk pH over August of 2016 (pre-application) and of 2017 to 2020, the month the soil cores were taken (Beerling et al. 2024, Methods; SI Fig. S7), scored against the observed 0 to 10 cm value and the observed 10 to 30 cm value of that year separately and pooled: both treatments, twenty comparisons. The model is one lumped 0 to 30 cm box (`Zbio` = 300 mm) and the field holds a gradient the box cannot, so the two intervals bracket where a lumped value can fall; `pH_season_0_10_RMSE` and `pH_season_10_30_RMSE` (each with control and basalt variants) show which interval the model tracks. The cores were taken inside the peak-biomass quadrats in August, with no day recorded.
3. `pH_difference_discrete_RMSE`: the modelled basalt-minus-control August pH against the observed difference at 0 to 10 cm and at 10 to 30 cm, pooled (ten comparisons); `pH_difference_0_10_RMSE`, `pH_difference_10_30_RMSE` and `pH_difference_2020_model_over_observed_0_10` sit beside it. Sampling and laboratory variability common to both plots cancels; this is the treatment signal.
4. `base_sat_difference_discrete_RMSE`: the modelled basalt-minus-control August base saturation of the lumped box (Ca, Mg, K and Na sites, per cent of CEC) against the observed difference at 0 to 10 cm and at 10 to 30 cm from the digitised Fig. 1F (`Data_Beerling.mat`, `CEC_Exp`), pooled (ten comparisons); `base_sat_season_discrete_RMSE` scores the absolute values the same way as the pH, and the 0 to 30 cm value is the soil-mass average of the two intervals. The field reached pH 6.9 at 83 to 90 % base saturation; an exchanger that saturates to reach that pH fails this objective however well it matches the pH.
5. `NO3_2017_relative_RMSE` and 6. `carbon_GPP_NEE_NRMSE`, unchanged.

`Front_N2O_pH_season`, `Front_N2O_pH_difference` and `Front_pH_difference_base_sat_difference` flag the two-objective fronts. Diagnostics that remain in the combination CSV: the day-matched `N2O_NRMSE`, the annual `pH_RMSE` (pre period 1 January to 9 November 2016, then calendar years, ideal hydrogen-concentration means), the 0 to 30 cm hydrogen-activity-average scores `pH_season_RMSE` and `pH_difference_RMSE` (the single target of objective sets v9 and earlier) with their per-treatment `pH_season_control_RMSE` and `pH_season_basalt_RMSE`, the per-period `pH_difference_pre` to `pH_difference_2020`, and the pooled `pH_envelope_miss_RMS` with its interval forms `pH_envelope_miss_0_10_RMS`, `pH_envelope_miss_10_30_RMS` and `pH_envelope_miss_discrete_RMS` (the timing-agnostic scores: the soil samples carry no date within their year, so an observation inside the model's August daily range counts as matched) (RMS over the ten comparisons of the distance of an observation outside the model's growing-season daily range, zero when inside). Per trajectory, `joint_sweep_trajectory_summary.csv` adds the growing-season means and daily range (`bulk_pH_season_*`, `bulk_pH_season_min_*`, `bulk_pH_season_max_*`), `pH_season_RMSE`, `pH_envelope_miss_SSE`, and the per-year `N2O_cumulative_<year>_window`, `_observation_count`, `_observed_kgN_ha`, `_modeled_kgN_ha` and `_log_ratio`. `joint_sweep_summary.status` records `objective_set=cumulative_N2O_discrete_depth_pH_base_sat_v11` and the two front counts. No winner is selected automatically.

## Restart behavior

The rolling checkpoint interval is 30 minutes. A requeue keeps the job ID, so the same task can resume its checkpoint. Completed results are accepted only when their validation fingerprint and complete parameter provenance still match.

If one case fails, the summary still runs and names it. Rerunning the same full-stage command recovers every completed case and recomputes only the missing ones, provided no fingerprinted file changed. A case killed by the scheduler (time limit, `scancel`, node loss without requeue) leaves its `locks/run_NN.lock` behind; the controller then refuses and prints the exact `rm -r` command for each stale lock. Never clear the lock tree while production is running.

## Dissolution factor (held at 1 from package t, 2 October 2026)

Packages p to s swept Run_One_Subfield's `DissolutionFactor` (the workbook
`diss_f` of every applied feedstock mineral) as a fourth dimension, because
at factor 1 the basalt run released 12.6 eq m-2 of alkalinity by August
2020 against 10.3 eq m-2 of exchanger acidity. Package t traces that excess
to one rate constant, the amphibole's, corrects it
(`AMPHIBOLE_KINETICS_2026-10-02.md`) and returns the sweep to package o's
three grids. Every case passes `DissolutionFactor` 1, so result tags keep
`_fd1`, the provenance check still compares the factor, and
`joint_sweep_combination_summary.csv` keeps its `DissolutionFactor` column.
Case 1 to 3 are the reference tuple, so the one-year control gate and the
all-12 campaign label `a1p3_kn0p0048_kbfree_fd1` are unchanged. The factor
is still the uncertainty ensemble's parameter (0.8 and 1.2). The export
operator (the bottom leakage over the whole column's water that T&C-BG and
the published T&C-SMEW use) is kept as inherited, by decision of 1 October
2026. The calcite run's bulk pH is read from a paste in which the plot's
carbonates dissolve (`PASTE_CALCITE_2026-10-01.md`,
`PASTE_MG_CARBONATE_2026-10-01.md`).

Score the basalt effect on `Front_pH_difference_base_sat_difference`, and
read the control's `pH_season_control_RMSE` on the same rows: a steeper
alpha raises the basalt effect and deepens the control's acid offset.
