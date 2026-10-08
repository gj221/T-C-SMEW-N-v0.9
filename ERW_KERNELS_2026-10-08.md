# Faster ERW chemistry, bit for bit, 8 October 2026

The ten-minute ERW step (`ERW_function`) can now run through restructured
kernels that compute the same floating-point operations in the same order,
so every output is bit-identical to the original code. A switch selects
them; it defaults to the original code, so nothing changes unless a run
asks for it.

| `TC_ERW_KERNEL` | What runs | Runs on thread-pool workers |
|---|---|---|
| unset or `matlab` | the original code, unchanged | yes |
| `fast` | restructured MATLAB kernels | yes |
| `mex` | the same kernels compiled with MATLAB Coder | **no** |

`bulk_pH_params.erw_kernel`, when set, overrides the environment variable.
MATLAB (R2025a and R2026a) refuses MEX functions on `parpool('Threads')`
workers, and the 2-D frame runs its daily ERW columns on such a pool. So
the field-grid campaign would use `fast`. `mex` serves 1-D runs and serial
2-D runs (`TC_PARALLEL_WORKERS` 0 or 1).

## Package z: the kernels on package y

The change was built and checked on package x (the sections below). Package
z puts it on package y, whose squares share the thread pool in every loop.
The kernel commit applies to y as it is, except two files merged by hand.
Every line added or removed is byte for byte the one x's commit adds or
removes:

- **`ERW_function.m`:** y changed only `solve_erw_branch`, where the
  `lsqnonlin` options are now built once per process. That is clear of the
  three dispatch points and the helper. The file mixes CRLF and LF lines,
  and the merge keeps each line's ending.
- **`ERW_Daily_Step_2D.m`:** y rewrote it for daily blocks. The three lines
  still follow `suppress_no_converge_warning`, before the step's one
  `ERW_function` call.
- **`ERW_Appelo_State_At_H.m`:** y has its own bit-identical edit (the
  pH-only terms once, the Al(OH)3 equilibrium only at the converged f_H).
  On z, `test_erw_chemistry_kernels` compares `fast` with that edited
  original, so a pass checks y's edit as well.

z adds:

- **`run_provenance.erw_kernel`** (`Run_One_Subfield`): `matlab`, `fast` or
  `mex`, read from `TC_ERW_KERNEL` when the run starts. `mex` with a thread
  pool stops there (`Run_One_Subfield:MexOnThreadPool`). A `mex` run checks
  its MEX files there too, not at the first step.
- **`compare_erw_kernels.m`:** two runs of one root under two kernels,
  every saved value compared, inside structs and cells too. A daily-store
  result keeps most of its numbers in `ERW_daily`, `Final_state` and
  `Field_hourly`. `compare_u_v` lists a struct that differs without failing.
- **`compare_erw_kernels.sbatch`, for y:** the runs take the Base_1D set,
  and run 2 keeps daily stores as the field runs will. With `X_ROOT` set,
  x runs both as well and each fast run is compared with x's
  (`compare_x_y.m`).
- **`compare_x_y.sbatch`:** `CMP_Y_KERNEL=fast` runs y's side with the
  fast kernels.

Checked locally on z (R2026a, 8 October; the laptop ran three other MATLAB
jobs, so times are rough):

- **`checkcode`:** no new message in the edited files.
- **`test_erw_chemistry_kernels`:** passed against y's code, for `fast`
  and, after `Build_ERW_Chemistry_Mex` in a test copy (MinGW 13.2, 43 s),
  for `mex`. That is 600 reduced solves with roots at pore pH 3.50-9.45
  (12 out-of-scan errors raised identically), 200 charge scans, all seven
  outputs of the 23 fixture variants, and three 12-step chains.
- **`mex` in a run:** `Run_One_Subfield` stops a 4-worker run with
  `Run_One_Subfield:MexOnThreadPool` before it starts, and runs a serial
  one-column day under `mex`.
- **13 more tests** that read or call the edited files pass, five of them
  again under `fast`. `test_campaign_resume_provenance` took 131 s, and
  68 s under `fast`.
- **Run 2** on the 12-column plot, 3 days from 9 November 2016, daily
  stores, 4 thread workers, `fast`:
  - against x (serial), `compare_x_y`: PASS; every daily output and the
    final state are identical;
  - against y's original code on 4 workers, `compare_erw_kernels`: 332 of
    342 variables identical and the provenance alike. Of the rest, 8 hold
    run times and result paths, and one is the checkpoint fingerprint,
    which hashes the two trees' different code;
  - time loop 0.026 h, against 0.068 h for y's code on 4 workers and
    0.220 h for x's, serial.
- **Run 11** (one column, full records), 3 days from 9 November 2016,
  `fast`, against x: `compare_x_y` PASS, 558 of 571 variables identical
  and the other 13 different by design (times, paths, provenance, code
  identity). Time loop 0.007 h, against 0.020 h for x.
- **`compare_erw_kernels` catches** a one-ulp change planted in `ERW`,
  `ERW_daily` and `Final_state`.

On Bouchet, package z follows `RUNBOOK_2026-09-24.md`, section 0z.

## Where the time went

This is a profile of one step in isolation: the step-22053 fixture of
`test_erw_appelo_exchange_floor`, a Kantola column under Appelo exchange
with the Al(OH)3 pin. No model was run.

- **62 %: the reduced Appelo solve** (`ERW_Appelo_Mass_Conserving_State`).
  It scans 57 pH points, then bisects. Every trial H runs a fixed 60-step
  bisection on f_H through three nested functions. That is about 5,800
  inner evaluations and 17,000 function calls per step.
- **31 %: building `ERW_data`.** This is `who` plus about 350 `eval`s, once
  per activity pass. Neither frame keeps this output: the 2-D step receives
  it as `~`, and the 1-D frame never reads `ERW_data_ittn`.
- **About 1 %: the 281-point pore-water charge scan.**
- **lsqnonlin did not run.** The reduced solve already satisfies the system.

## What `fast` does

- **`ERW_Appelo_Seed_Kernel`.** It evaluates the original's expressions on
  the same operands. It removes only work whose result the original
  discards or repeats:
  - the Al(OH)3 inventory split runs once per trial H instead of at all 61
    values of f_H;
  - terms that depend on H alone are formed once;
  - the 57 scan points run their f_H bisections side by side.

  MATLAB's elementwise `+ - * / .^ 10.^` and its `sum` of four or six terms
  round exactly as the scalar left-to-right forms. This was checked on
  4,000,000 random operands.
- **`ERW_Charge_Root_Kernel`** does the same for the 281-point scan and its
  bisection.
- **`ERW_data`**:
  - It is built with one `eval` into the identical struct, field for field.
    That includes the partial copy of itself that the original loop stores
    in `ERW_data.ERW_data`.
  - It is not built at all when the caller sets
    `bulk_pH_params.skip_erw_data`. Both frames now set it.
- **A step with no root inside the scan.** There the original stops with
  `MATLAB:badsubscript`. The original then runs, so the error is its own.

## Bit-identity, local (R2026a, Windows, MinGW 13.2)

`test_erw_chemistry_kernels` compares everything bit for bit: class, size,
the IEEE bits of every float, field order, and function handles by text and
captured values. It passed for `fast` and `mex`:

- **600 reduced solves.** Roots range over pore pH 3.50-9.45. Cases also
  vary temperature 2-30 C, ionic strength 1e-3 to 0.1, alpha 1.3/3.4/6,
  CEC x0.5-2, water, untreated and treated pools, and an Al(OH)3 inventory
  that is none, inventory-limited or saturated. 12 steps outside the scan
  raised identical errors.
- **200 charge scans.**
- **All seven `ERW_function` outputs** on 23 variants of the fixture:
  - pore pH about 5 to 11;
  - the Al pin on and off;
  - 2, 12 and 30 C;
  - saturation 0.35-0.95;
  - alpha 1.3/3.4/6;
  - treated pools.

  Variants the original refuses raise identical errors.
- **Three 12-step chains**, each step fed the state the previous one
  returned.

The 31 isolated tests that read or call the edited files pass under all
three modes.

## Bouchet: bit-identical, and the speed on the production CPU

`build_erw_mex.sbatch` ran on 8 October 2026 as job 28759434 on a1114u18n01
(cpumodel 9575f), with R2025a Update 1 and gcc 11.5.0. It ran from a scratch
test copy, package x plus this change. The MEX built in 80 s and passed its
reference check. `test_erw_chemistry_kernels(3000,48)` passed for both `fast`
and `mex`, every value bit for bit:

- 3000 reduced solves, roots at pore pH 3.50-13.91; 57 steps outside the
  scan raised identical errors;
- 1000 charge scans;
- all seven `ERW_function` outputs on the 23 fixture variants;
- three chains of 48 steps.

ERW_function step time on one computational thread, median of 25
interleaved rounds:

| pore pH | activity passes | original | fast | mex |
|---|---|---|---|---|
| 10.9 | 1 | 26.9 ms | 3.43 ms (7.9x) | 1.92 ms (14.0x) |
| 7.8 | 2 | 53.8 ms | 6.39 ms (8.4x) | 3.06 ms (17.6x) |
| 7.1 | 2 | 52.6 ms | 5.67 ms (9.3x) | 2.83 ms (18.6x) |
| 6.5 | 2 | 90.4 ms | 5.45 ms (16.6x) | 2.62 ms (34.5x) |
| 5.8 | 2 | 46.9 ms | 3.31 ms (14.2x) | 2.16 ms (21.7x) |
| sum | | 270.7 ms | 24.3 ms (11.2x) | 12.6 ms (21.5x) |

The original's 27-50 ms per step agrees with the README cost model:
2.1 h per square per five years is about 29 ms per ten-minute step.
lsqnonlin ran in none of these steps.

## Bouchet: the full model, bit for bit

`compare_erw_kernels.sbatch` ran on 8 October 2026 as job 28761479 on
a1132u37n02 (no CPU pin), from the same scratch copy. It compared every
saved variable with `compare_u_v`:

| comparison | result |
|---|---|
| run 11 (basalt, one column, 375 days, two applications), `fast` against `matlab` | no numeric model variable differs; ERW column identical; PASS |
| run 11, `mex` against `matlab` | the same |
| run 2 (basalt on the 12-column plot, 40 days, 6 thread-pool workers), `fast` against `matlab` | no numeric model variable differs; all 12 ERW columns identical; PASS |

The five runs shared the node:

| run | `matlab` | `fast` | `mex` |
|---|---|---|---|
| run 11 | 2105 s | 435 s (4.8x) | 396 s (5.3x) |
| run 2 | 732 s | 362 s (2.0x) | – |

- **Run 11's 4.8x** is what an 11x faster ERW step gives a whole run when
  the step is about 88 % of it.
- **Run 2 gains less in package x.** Its hourly T&C loop is serial over the
  12 columns, while the ERW columns share 6 workers. With package y's
  parallel square loops, a field run's gain approaches run 11's.

## Speed, local

These are medians of 15 interleaved rounds (`erw_kernel_benchmark`). The
laptop was shared with other jobs, so its ratios ran higher than Bouchet's
(above).

| pore pH | original | fast | mex |
|---|---|---|---|
| 10.9 (1 activity pass) | 477 ms | 25.8 ms (18x) | 17.4 ms (27x) |
| 7.8 | 271 ms | 15.9 ms (17x) | 10.1 ms (27x) |
| 7.1 | 284 ms | 15.5 ms (18x) | 10.1 ms (28x) |
| 6.5 | 245 ms | 14.0 ms (18x) | 8.9 ms (28x) |
| 5.8 | 226 ms | 12.8 ms (18x) | 8.6 ms (26x) |

These are control-column steps. A treated column also runs the mineral
kinetics and particle-size updates, which are unchanged, so its gain is
smaller. `compare_erw_kernels.sbatch` measures it.

## Files

New, in `v2/T&C_Code`:

- `ERW_Chemistry_Kernel_Mode.m` (the switch)
- `ERW_Appelo_Seed.m` (dispatch and error fallback)
- `ERW_Appelo_Seed_Kernel.m`
- `ERW_Charge_Root_Kernel.m`
- `ERW_Chemistry_Mex_Check.m` (in `mex` mode, refuses a missing or
  non-identical MEX on first use)
- `Build_ERW_Chemistry_Mex.m`
- `Tests/test_erw_chemistry_kernels.m`, added to the post-legacy audit

Edited:

- `ERW_function.m`: three dispatch points and a local helper. With the
  mode at `matlab` the original statements run unchanged.
- `ERW_Daily_Step_2D.m` and `MAIN_FRAME_BG_ERW_10mins.m`: the
  `skip_erw_data` flag, three lines each.
- `v2/run/Run_One_Subfield.m` (package z): `run_provenance.erw_kernel`,
  and the refusal of `mex` with a thread pool.

At the package root:

- `build_erw_mex.sbatch`
- `compare_erw_kernels.sbatch` (in z, for y's records and with `X_ROOT`)
- `compare_erw_kernels.m` (package z)
- `erw_kernel_benchmark.m`
- `compare_x_y.sbatch` (package y; z adds `CMP_Y_KERNEL`)

The campaign fingerprints hash `*.m`, so the new sources enter them. The
MEX files (`*_mex.mexa64`) and their manifest do not. Checkpoints hash all
of `v2/T&C_Code/*.m` and not the mode, so a run may resume across modes.

## On Bouchet

1. **Feasibility: done on 8 October.** MATLAB R2025a Update 1; the Coder
   license test and checkout succeed; MEX selects gcc 11 (system gcc
   11.5.0). The one-line check, on `devel`, which starts within a minute:

   ```bash
   module load MATLAB/2025a && srun --account=pi_np363 -p devel -t 00:15:00 -c 1 --mem=8G matlab -batch "disp(version); fprintf('Coder license %d\n',license('test','MATLAB_Coder')); [ok,msg]=license('checkout','MATLAB_Coder'); fprintf('checkout %d %s\n',ok,msg); cc=mex.getCompilerConfigurations('C','Selected'); if isempty(cc), disp('no C compiler selected'), else, fprintf('%s %s\n',cc.Name,cc.Version), end" && gcc --version | head -n 1
   ```

2. **Build and test in the deployed root, once per root, before
   validation.** This passed in the scratch test copy (above). A MEX
   belongs to its root, so the package that carries this change runs it
   again in its own root. Look for `PASS: ERW chemistry MEX files built
   ...`, the four `PASS` lines of `test_erw_chemistry_kernels`, and the
   timing table. The `[ERW cap]` and `[ERW fail]` lines are the original
   code's prints for the fixture variants it refuses. To read past them:
   `grep -v "^\[ERW" build_erw_mex_<job>.out | tail -n 40`.

   ```bash
   sbatch --account=pi_np363 build_erw_mex.sbatch
   ```

3. **Full model, before switching production: passed on 8 October**
   (above), on x with `compare_u_v`. It is a model run, so a later package
   repeats it only with your approval. In package z each comparison ends
   with `PASS: every saved value is identical under both kernels.`, each
   comparison with x (`X_ROOT`) with `PASS: y reproduces x`, and the job
   with `ALL PASS`; the `finished in` lines give the wall times
   (`RUNBOOK_2026-09-24.md`, section 0z, has the commands).

   ```bash
   sbatch --account=pi_np363 --export=ALL,X_ROOT=<x root> compare_erw_kernels.sbatch
   ```

4. **Production.** Run `export TC_ERW_KERNEL=fast` before the campaign
   controller; it submits with `--export=ALL`. Package z records the mode
   in `run_provenance.erw_kernel`; asserting it in the campaign task is left
   to the campaign integration.

## If a later release does not reproduce bit for bit

Bouchet's R2025a reproduced the original bit for bit at the function level
(above), and `compare_erw_kernels.sbatch` checks the full frame. A future
MATLAB release or another platform could change how MATLAB's elementwise
math rounds, or how Coder builds. `test_erw_chemistry_kernels` would then
fail in the audit, and `fast` or `mex` would need this validation before
use:

- the 35-year 1-D central configurations, old against new;
- daily pore and bulk pH within 0.001;
- exchange fractions within 1e-6;
- cumulative weathering and CDR within 0.1 %;
- N2O and NO3 leaching within 0.1 %.
