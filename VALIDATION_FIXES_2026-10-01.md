# Validation fixes for the Bouchet deployment, 1 October 2026

Package o. The model is package n's, unchanged: no parameter, rate, pool,
driver or input differs from n. Five files change, all in the validation
path, plus the runbook and this note.

## What happened to n

Both 1-D validations of n (`Base1Dcheck` 28003062 and `Base1Dcheck_case8`
28003066) failed after 12:21 and 11:25. In `Validate_Base_1D` that is the
end of the audit suite; the pre-application identity stage, which takes
about 30 minutes, and the smoke, which follows it, were not reached. The
2-D validations run only the four configuration tests.

The audit's 89th test, `test_2d_driver_equivalence`, refuses the saved
equivalence baseline unless `ERW_Equivalence_Fingerprint` of the checkout
matches the fingerprint stored in it. That fingerprint hashed every
`Inputs/*.mat`, and the three scenario forcings `Data_Kantola_ssp*.mat`
are copied into a root apart from the package. A baseline built before that
copy, or a root whose forcings were replaced after the build, fails there
although the baseline reads none of those files. The regen job was run
twice in the n root (27999185 and 28002844), which fits that reading; the
Slurm log of 28003062 confirms it (runbook, section 0).

## What o changes

1. `v2/T&C_Code/ERW_Equivalence_Fingerprint.m`: the baseline fingerprint
   names the four inputs the 1-D control run reads (`Data_Kantola.mat`,
   `Data_Kantola_BG.mat`, `Ca_Data.mat`, `All_deposition_data.mat`, the
   same four the checkpoint identity and the campaign stamp list) instead
   of `Inputs/*.mat`. The workbook, the rock file, the code, the parameter
   scripts, the run entry points and the MATLAB release stay in it.
2. `v2/T&C_Code/Tests/test_2d_driver_equivalence.m`: a refused baseline
   prints the fingerprint entries that differ, as `baseline had` and
   `checkout has` lines with the file path and the first 12 hex digits of
   its hash, before raising `test_2d_driver_equivalence:StaleBaseline`.
   The file keeps its CRLF line endings.
3. `runtime/Base_1D/Validate_Base_1D.m` and
   `runtime/SSP_1D/Validate_SSP_1D.m`: a timestamped banner
   (`[yyyy-mm-dd HH:MM:SS] Validation stage: ...`) before each stage, so
   `slurm_<job>.out` shows where a validation stopped. The stages and
   their order are unchanged.
4. `v2/T&C_Code/Tests/test_erw_lateral_cation_daily_delta.m`: the one
   audit test n's salt closure broke. It compares the ERW `Ca_tot`,
   `Mg_tot`, `K_tot` and `Si_tot` at every daily boundary with the T&C
   pools. The last boundary precedes the final flush, which revisits the
   same P row without its delta (`ERW_Should_Apply_TC_Daily_Delta`), so
   that anchor is the pre-flush total while the saved P row holds the
   post-flush state. The flush adds the boundary supply for every row
   after the boundary; with n's leached-salt return that supply is no
   longer zero for Ca and Mg, and the saved total exceeded the anchor by
   exactly `I_Ca` times the flush length (440.87 umol m-2 against
   460.04 umol m-2 d-1 over 138 ten-minute rows in the test's 6-day run
   7; relative 2.8e-5 for Ca, 5.6e-5 for Mg; days 1 to 5 and K and Si
   exact). The test now removes that supply from the last expectation,
   using the saved `E.const.I_*`, `E.dm_end` and `E.dtw`, and keeps its
   1e-10 tolerance. The model, the write-back and the final-pool check
   were already consistent; only the expectation predated the supply.
   On Bouchet this test is number 94 of 98, after the equivalence test.
5. `RUNBOOK_2026-09-24.md`: the o root, a baseline freshness check to run
   before any submission, and the n log commands.

The pin lists (`runtime/Base_1D/shared_pins.sha256`, `ssp_pins.sha256`)
are regenerated for the two changed pinned files. Every campaign stamp and
the baseline change with the source, so o needs its own R2025a baseline and
validation, as the runbook orders.

## What was run before packaging (PC, MATLAB R2026a, Windows)

R2025a is not installed on the PC, so these runs establish that the tests
and the model execute, not the R2025a numerical identity, which only the
Bouchet validation can.

* Package n, unchanged: the equivalence baseline built (10 days, 133 s);
  the complete list of `run_post_legacy_equation_audit_tests` (95 tests
  and the plot regression) run with every test wrapped so that all
  failures would be collected, then `test_bouchet_linux_portability` over
  every deployed campaign directory,
  `test_biogeochemistry_io_path_portability` and
  `test_all_12_configuration`: 98 entries, 96 passed. The two failures:
  `test_erw_lateral_cation_daily_delta` (fix 4 above) and
  `test_plot_configuration_figures`, which died in `exportgraphics`
  inside R2026a's web figure renderer after 650 s (a JavaScript stack
  trace and "Could not load TeX fonts") with the plotting code unchanged
  since m; rerun alone on the o tree it passed (885 s), so it was the
  PC's renderer under memory pressure, not the package, and Bouchet
  R2025a renders through the classic path that passed m.
* The three-day smoke's one-column run passed `Assert_Base_Provenance`
  with the Slurm variables emulated; its twelve-column run was stopped by
  the PC's memory limit (commit charge 65 of 69 GB with three MATLAB
  sessions open), not by the model.
* Package o: the baseline rebuilt with the new fingerprint (268 s); the
  configuration, portability, pin-list, calcite-case, campaign-stamp,
  checkpoint-dependency, launcher, resume-provenance and equivalence
  tests rerun (15, all passed); the stale-baseline report exercised
  against a baseline whose stored fingerprint was altered (it named the
  file and raised `StaleBaseline`; the restored baseline matched); the
  runbook's freshness one-liner printed `baseline fresh`; the corrected
  `test_erw_lateral_cation_daily_delta` rerun on the same 6-day run 7
  (13 of 13 checks passed; worst daily shadow error 6.4e-15 for Ca and
  5.8e-15 for Mg, from 2.8e-5 and 5.6e-5 before the correction).
