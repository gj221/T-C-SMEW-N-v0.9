# Field-grid records kept daily, squares in parallel, 8 October 2026

Package y is package x with changes to how a multi-square run holds,
computes and saves its state, for the 108 35-year runs on the 10 m Kantola
grids (340 and 369 squares). No equation, parameter, input or solver
setting changes: every number the model computes is x's, bit for bit. What
changes is what a multi-square result keeps, and that the squares' work
runs on the thread pool. One-column runs keep x's full records, so every
1-D tool, campaign and check is unchanged.

Package x needed more than 1.5 TB of memory and 0.3-1 TB per result for a
35-year field-grid run (`field_profile` job, 8 October): it held every
hourly and ten-minute row of every square, copied them into the daily
loops, and saved the whole workspace.

## 1. What a run keeps

**One column** (every 1-D campaign, the 1-D sweeps, the equivalence
test): x's full records, unchanged.

**Several squares** (the 12-column plot and the field grids): the list
Greg approved on 8 October (`PACKAGE_Y_OUTPUT_LISTS.md` in the session
notes), in the result file:

- **Daily stores**, every square and day, under x's names and shapes
  (`DAILY_STORE_NAMES`): `P` (all 79 pools), `N2Oflx`, `N2Oflx_nit`,
  `N2flx`, `VOL`, `NITflx`, `DENITflx`, `Min_N`, `BfixN`, `FertN_d`,
  `DepN_d`, the ten `LEAK_*`, `PH_bulk_d`, `PH_pore_d`, `R_litter`,
  `R_microbe`, `R_ew`, `NPP_H`, `NPP_L`, `LAI_H`, `LAI_L`, `Ccrown_t`,
  `RB_L`, with `Domain` and `lateral_diag` as in x. `PH_d` repeats
  `PH_bulk_d` row for row and is not stored.
- **From the hourly records**, per square and model day (the hours each
  daily step averages: row 1 none, row 2 the first hour alone, then 01:00
  to 00:00): the mean water content of the active layers `O_d` (the layers
  with `Bio_Zs > 0`) and water-table depth `ZWT_d`, and the totals [mm] of
  drainage `Lk_d`, infiltration `f_d`, run-on `q_runon_d` and runoff
  `Rh_d`, `Rd_d`. The hours after the last daily step are in `Hourly_tail`,
  so run totals close (`Field_Daily_From_Hourly`).
- **From the ten-minute chemistry**, per square and ERW day (the rows one
  daily call computes: 01:00 to 00:50, five rows on the first call, the
  final flush its own entry, NDays+1 entries): `ERW_daily{c}`
  (`ERW_Daily_Record`), with the end-of-day values, the means, the totals
  and the largest closure residuals of the approved list, each call's
  first and last ten-minute row and its daily row. Each call is
  summarised from the values its rows finally hold (after the next hours'
  lateral moves and the next call's anchors), as x's arrays hold them at
  the end. `ERW{c}.daily_initial` is the t = 0 row.
- **Hourly field means** of the tower fluxes (`Field_hourly`,
  `Field_Hourly_Means`), weighted as `Domain` weighs: `Rn`, `H`, `QE`,
  `G`, `Ts`, `ET` (as the 1-D frame defines it), `T`, `EG`, `GPP`, `O` and
  `Tdp` of the active layers, `Lk`, `Runoff` (`Rh + Rd`).
- **The final state** (`Final_state`): the last 50 hours of every hourly
  array and the last two rows of every daily array, in order, and every
  square's ERW state (its last day's block of ten-minute rows); the
  carried per-square state and the run's inputs as in x.
- **Particle-size histories**: final state only (the rows set at
  initialisation and the last row of each of the last three ERW days),
  for every amended square, one application or several
  (`weath.particle_rows` gives their rows on the ten-minute grid). A
  full-record run keeps them as x does.

`TC_FULL_RECORDS=1` (or `0`) forces either layout, whatever the number of
squares; the tests and the comparison job use it. It enters the
checkpoint identity (`records=full|daily`).

## 2. Memory

- **Hourly arrays** keep a window of 50 hours (`HOURLY_WINDOW`): hour h in
  row `mod(h-1,50)+1`. The model reads hour i-1, the 24 hours a daily step
  averages, and the ERW water bins, which reach back at most to hour i-49.
  Every hourly array is written every hour, so a row is always the hour it
  stands for; every read is of an hour inside the window.
- **Ten-minute chemistry** (`ERW_Daily_Step_2D`): `E.bg` and `E.weath`
  hold one call's block. Block row 1 is the previous call's last row, and
  the call's own rows follow; the day's code indexes the block, so it reads
  exactly the rows the full arrays gave it. The anchored row is held twice
  (the old block's last row and the new block's first) and both take the
  anchors before the old call is summarised. `ERW_MOD_PARAM_excel` takes an
  optional allocation length (two rows here; the 1-D frame and every other
  caller allocate the whole grid as before). The rock schedule is built
  once, on the whole grid, exactly as x builds it, and kept as the rows
  that carry an addition.
- **Daily arrays** keep their full length in memory (the model reads rows
  j-1 and j; the stores are saved whole, the others in the final state).
- **No full-length series per day**: x rebuilt the water series `Tw`,
  `Lw`, `Iw`, the CO2 series and the rock schedule over the whole run in
  every square every day, and passed each square's whole hourly and daily
  history into the ERW loop. y passes the window: the hours from i-49, two
  daily rows of `P`, the CO2 of the window (the whole series on the first
  call, which seeds the state from the full-length timetable as x does).
- **Measured** on the 3-day local runs (section 9): per square and day the
  daily arrays hold 672 values (5.4 kB), all kept in memory, of which a
  result saves the 140 of the stores, plus 79 per ERW day. x also held
  2.9 kB per square-hour of hourly records and 1.5 kB per square per
  ten-minute row. For 35 years on 369 squares that is, for y, about 25 GB
  of daily arrays and 3 GB of ERW summaries in memory and about 8.3 GB per
  result before compression (5.3 GB stores, 3.0 GB ERW summaries); for x,
  about 330 GB of hourly and 1,050 GB of ten-minute records.

## 3. Speed

- **The squares' work runs on the thread pool**, square by square, in
  three `parfor` loops bounded by `ParallelWorkers` (`TC_FIELD_WORKERS`
  in the field campaign): the daily vegetation and biogeochemistry, the
  hourly vertical physics (`HYDROLOGIC_UNIT`), and the daily ERW chemistry
  as in x. Nothing a square computes reads another square inside these
  loops; the lateral transfer of water, heat, solutes and runoff couples
  the squares and stays serial between them. The per-hour radiation
  projection on the slopes (`Topographic_Radiation_Projection`) is made
  once per hour before the loop.
- **Carbon-reserve events** (`Plant_Carbon_Deficit`) are taken on the
  worker after each square's day and replayed on the client in square
  order, so the run's report is the serial loop's, bit for bit. The
  nutrient-warning throttle of `Nutrients_Available` is per worker; it
  decides only which warnings print.
- **Chemistry housekeeping**: the `lsqnonlin` options of the ERW solver
  are built once per worker instead of at every step (`ERW_function`);
  `ERW_Appelo_State_At_H` forms its pH-only terms once per trial pH
  instead of at every bisection step, and evaluates the aluminium
  equilibrium (dissolved Al3+, the solid and the inventory residual) only
  at the converged exchanger state, since the bisection reads only the
  site sum (`ERW_Al_Equilibrium_Target` is pure arithmetic, so skipping
  its unused calls changes nothing). This is the extra Greg did not veto.
- **Measured** on the 3-day local runs of 12 squares (a shared laptop, so
  only the ratios mean much): the time loop took 0.220 h in x (ERW pool
  off), 0.205 h in y serial and 0.068 h in y on 4 workers, 3.0 times
  faster. The field grids have 28-31 times as many squares to share out;
  the comparison job's grid pair gives the first Bouchet figure.
- **Not in y**: the restructured ERW kernels of the session "Compile the
  ERW chemistry step" (`TC_ERW_KERNEL=fast`, bit-identical, about 11 times
  faster per ERW step on Bouchet) are a patch on x that applies to y with
  small manual merges in `ERW_function.m` and `ERW_Daily_Step_2D.m`.

## 4. Checkpoints and results

A daily-store run checkpoints every variable but the daily arrays outside
the stores, whose last two rows go into `ckpt_daily_rows`, uncompressed
(`-nocompression`), and resumes from it exactly. Its result holds the
stores, the hourly field means, the final state and the run's inputs and
metadata, not the hourly windows under the hourly names (a reader would
take a window for the whole run). Checkpoints of x, and of the other
records layout, are refused by the checkpoint identity, as any code change
refuses them.

## 5. Readers

- `Field_Cell_Series` (`Plot_Field_Maps`, `Plot_Field_Layout`) reads the
  hourly-derived maps (WFPS, soil moisture, water table, runoff, run-on,
  infiltration, drainage) from the daily stores when the result has no
  hourly record. It reads the hourly record when there is one (x, one
  column), as before. The two agree to rounding: WFPS from daily mean
  water contents rather than the daily mean of hourly WFPS, sums in a
  different order.
- `Kantola_Bulk_pH_Model_Summary` reads `ERW_daily` when the result has
  no ten-minute records, for any square (fourth argument, default 1). The
  closure check then uses each day's largest residual and the sign check
  the day-end fractions.
- `test_subfield_12_smoke_results` reads the final state, the daily stores
  and the tail of a daily-store result.
- `Plot_1D_Then_All_12`, `budget`, `phobs`, `readout_1d`, `calcite_check`,
  `Summarize_Sweep`, the N audit and `Check_Lateral_Conservation` read
  one-column results or what both layouts keep; unchanged.

## 6. Checks

- **x against y on Bouchet** (`compare_x_y.sbatch`, RUNBOOK section 0y):
  run 2 from 1 November 2016 in both roots: 40 days on the 12-column plot
  across the rock application of 10 November (daily stores, and again
  with `TC_FULL_RECORDS=1`), 400 days on the plot across the applications
  of 2016 and 2017 (the particle histories of two applications in daily
  blocks), and 40 days on the 10 m basalt DEM grid. Every variable both
  results hold (paths into the two trees read with x's root as y's), every
  hourly-derived store, the hourly field means, every square's ERW day
  summaries and t = 0 row, and the final state must equal the same
  reduction of x's records with `isequaln` (`Daily_Store_Check`, which
  `compare_x_y.m` prints).
- **Tests**: `test_2d_parfor_equivalence` (run on demand, as before; not
  in the audit list) now runs run 2 for three days from 9 November 2016
  serial, on four workers and with full records: the parallel result must
  equal the serial one in every saved variable but timing and worker
  count, and the serial daily stores must reproduce the full records
  (`Daily_Store_Check`). The comparison job holds y's pool runs to x's
  serial hourly loop as well. `test_2d_driver_equivalence` holds
  the one-column run to the 1-D frame bit for bit, as before; the
  equivalence baseline must be rebuilt for the new code
  (`regen_equiv_1d.sbatch`).
- Local runs on this PC: section 9.

## 7. Also in y

`Load_Field_DEM` read the EPSG code of the DEM with a regular expression
written with backslashes, which `test_bouchet_linux_portability` rejects,
so package x fails the gate `Run_Field_2D_Task` runs before the field
campaign. The expression now uses one-character sets for the brackets and
`[0-9]` for the digits, which match the same text.

## 8. Files

New: `v2/T&C_Code/ERW_Daily_Record.m`, `Field_Daily_From_Hourly.m`,
`Field_Hour_Values.m`, `Field_Hourly_Means.m`, `Daily_Store_Check.m`;
`compare_x_y.m`, `compare_x_y.sbatch` at the package root.

Changed: `MAIN_FRAME_BG_ERW_2D.m` (records, stores, the three loops,
checkpoints, final state), `ERW_Daily_Step_2D.m` (windows, blocks, the
call summaries), `ERW_Particle_History.m` (`keep_rows`: no rows of dropped days, single
applications kept too; block length),
`Parameters/ERW_MOD_PARAM_excel.m` (allocation length), `ERW_function.m`
(the `lsqnonlin` options), `ERW_Appelo_State_At_H.m`,
`Plant_Carbon_Deficit.m` (`take`, `replay`), `Load_Field_DEM.m`,
`v2/run/Run_One_Subfield.m` (the result's variables),
`v2/plotting/Field_Cell_Series.m`, `Field_Cell_Variables.m` (a comment),
`v2/calibration/Kantola_Bulk_pH_Model_Summary.m`.

Tests: `test_2d_parfor_equivalence` (rewritten), `test_parallel_erw_configuration`
and `test_bouchet_2d_parallel_configuration` (three loops),
`test_erw_exchange_pool_coupling` (the two daily rows of `P_c`),
`test_checkpoint_omits_cleanup` (the daily-store save),
`test_kantola_bulk_pH_calibration` (a daily-store fixture),
`test_subfield_12_smoke_results`; `test_2d_amendment_conservation`,
`test_erw_handoff_consistency`, `test_erw_lateral_cation_daily_delta`,
`test_margin_prairie_calcite_placement`, `test_parallel_erw_equivalence`
and `test_2d_topographic_smoke` read 12-column records and run with full
records.
