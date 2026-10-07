# Particle records kept daily, 4 October 2026

Package v is package u with two changes to how a run holds and saves its
state. Neither touches the model's equations, inputs or parameters, and
every value the model computes is u's.

## 1. Particle-class histories: the last three days whole, then daily

`psd`, `d`, `lamb`, `SSA` and `delta_d` are rows x size classes x
applications for every amended column: one row per ten-minute ERW step,
one slice per scheduled rock application. Over 35 years (1,840,896 rows,
35 applications, 7 classes) each is 3.6 GB, so the five took 16.8 GiB in
every amended column. In package u they were held in daily blocks during
the run and assembled into full arrays at its end, which for a moment held
them twice.

The model reads two things from them (`ERW_Daily_Step_2D`):

- the previous row, `dt_m-1`, of `psd`, `d`, `lamb` and `SSA`;
- for `psd` and `lamb`, each application's first row whose last size class
  is positive, its initial distribution, which the helper caches as it is
  written (`'first'`).

`delta_d` is never read back. Nothing else reads them: not the frame, the
lateral transport, the bulk-pH operator, the campaign checks, `budget`,
`phobs`, `readout_1d` or the plotting code. Only the parameter script
allocates them and `test_erw_particle_history` tests the helper.

`ERW_Particle_History` now keeps the blocks of the last three ERW days
whole and, of every older day, only its last written row, at its latest
value. Every previous-row read comes from today's block or the day before
it, as the daily windows follow one another without overlapping, so it
returns exactly what u returned (a window repeating the day before would
be served exactly too). A read of any other row of a dropped day stops
the run (`RowNotKept`) instead of returning a wrong value, and so does the
one path that would search the whole history again (`RescanNotKept`: a
rewrite of an application's first positive row with a non-positive value,
which cannot happen while each row is written once).

What the results now hold, per amended column: `weath.psd`, `weath.d`,
`weath.lamb`, `weath.SSA` and `weath.delta_d` at the rows set at
initialisation and at the last row of every ERW day, and
`weath.particle_rows`, each kept row's index on the ten-minute grid. For a
35-year run that is about 12,800 rows instead of 1,840,896. Columns
without these histories (controls, one application) are unchanged and get
no `particle_rows`.

Memory, from package u's 35-year runs (ssp126, 4 October): the whole-field
basalt run had reached 205 GB by half way and was expected to end at
365-405 GB against its 384 G request; the finished 35-year 1-D runs peaked
at 7 GB (control) and 48 GB (basalt), 34 GB of the difference being these
histories held twice. Without them a whole-field basalt run is expected
near 150-200 GB and a control near 105 GB.

## 2. SSP 2-D checkpoints every two hours

Each checkpoint writes the whole frame state, compressed, at about 25 MB a
second. Measured in package u's 35-year runs at a 30-minute interval (the
count of checkpoint lines in each worker log against its elapsed time):

| Run (ssp126) | Checkpoints | Average write | Share of run time |
|---|---|---|---|
| 1 control, whole field | 62 | 3.2 min | 10 % |
| 2 basalt, whole field | 35 | 7.8 min | 21 % |
| 3 calcite, whole field | 38 | 4.1 min | 12 % |
| 4 control, margin | 31 | 2.0 min | 6 % |

The writes grow with the run; the basalt run's took about 15 minutes by
half way. The SSP 2-D campaign now checkpoints every 120 minutes
(`runtime/SSP_2D/submit_ssp_2d.sbatch`, enforced in `Run_SSP_2D_Task`),
and the checkpoints are smaller by the particle histories. A job stopped
part-way loses up to two hours. The uncertainty ensemble keeps 30 minutes:
its runs requeue themselves every 24 hours, half an hour after Slurm's
warning.

## 3. Memory requests

The SSP 2-D request for 35 years and the ensemble's 2-D default are 256 G
(were 384 G). The week partition lets one user hold 1.5 TB at once, so six
runs go together at 256 G against four at 384 G; controls can be lowered
to 160 G with `scontrol update` before they start. Read the MaxRSS of the
first finished whole-field basalt run before lowering anything else.

## 4. What it means for the campaign

Expected per 35-year 2-D run: whole-field basalt about 50 h in u, about
36 h in v; whole-field calcite about 56 h and 45 h; controls about 38 h
and 34 h. With six to eight runs at a time instead of four, the 36
central runs take about ten days instead of about nineteen.

## 5. Checks

Edited files: `v2/T&C_Code/ERW_Particle_History.m`,
`v2/T&C_Code/Tests/test_erw_particle_history.m` (rewritten for the kept
rows: reads, a repeated day, a rewritten kept row, both refusals, the
v7.3 round trip and the assembled rows), `v2/T&C_Code/MAIN_FRAME_BG_ERW_2D.m`
(one comment), `runtime/SSP_2D/Run_SSP_2D_Task.m` and
`submit_ssp_2d.sbatch`, `runtime/SSP_1D/submit_ssp.sh`,
`runtime/Uncertainty_1D/submit_uncertainty.sh`,
`runtime/Uncertainty_2D/submit_uncertainty_2d.sbatch`,
`v2/T&C_Code/Tests/test_uncertainty_campaign.m` (the pinned 256G line),
the README and the RUNBOOK (section 0v). Nothing was run locally; the
audit runs in the campaign validations on Bouchet. Before the campaign,
RUNBOOK section 0v runs one year of 1-D basalt across two applications in
u and in v and compares every saved variable: all must be identical except
the particle histories, which must equal u's at v's kept rows.
