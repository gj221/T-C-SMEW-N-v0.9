# Particle history storage, 28 September 2026

Claude session: `ddd36c34-e868-4960-84ce-35090b8b4ab0`, titled
"Resume recent codex session". The last implementation was the particle
history helper in this tree. Package k is the comparison copy.

## Before changing the helper

The next steps are to remove duplicate history storage, test reload and
replayed windows, compare treated runs with k, then check the full audit
and package pins. The solver and management equations stay as written.

| Failure point | Check required |
|---|---|
| A replay replaces the first positive size-class row | Compare the cached row with a direct search after every write, including overwrites |
| Reload loses blocks or reads the wrong previous row | Save and reload the live state, then continue the same sequence |
| Full arrays coexist with blocks and exhaust memory | Remove the original arrays from the live weathering state; assemble them only for output |

## Equations and assumptions

1. A global row r belongs to a daily window starting at f at local index
   i = r - f + 1. A row is valid in that window when 1 <= i <= n.
   These are integer indices. They carry no physical units.

2. Read a history row from the latest written window covering r, otherwise
   from its initialized value. The initial arrays are zero except for
   explicitly populated rows. Keep those populated rows exactly, including
   any future row, rather than assuming that only row 1 can be populated.

3. For each application a, the initial distribution read by ERW_function is
   F_a = X(r_a,:,a), where r_a = min{r: X(r,end,a) > 0}. If that set is empty,
   the solver uses its existing current-step fallback. Overwriting r_a
   requires another search if its last class is no longer positive. Read
   cache updates from the row after assignment: before the first rock
   application the solver returns scalar zeros, which MATLAB expands
   across the size classes and applications.

4. Final assembly starts from the initialized values and applies the
   windows in write order. Later writes win. This is the same assignment
   rule used by the full arrays, including the final repeated window.

The physical equations consuming those rows are unchanged:

    delta_d = [d_1, diff(d)]
    lambda = a * d^b
    SSA = 6 * lambda / (d * rho_rock)
    psd_initial = class_fraction * M_added / delta_d
    d_shrink = sum(Wr * rock_fraction * MM_min / rho_rock)
    d_loss = min(d, 2 * d_shrink * lambda * dtw)
    psd = psd_initial * (d/d_initial)^3 * (lambda_initial/lambda)
    SA = sum(SSA * psd * delta_d)

The diameter update is bounded by the diameter already present. SSA and
the distribution update only divide by a positive remaining diameter.
The initialized diameter classes and their widths must be positive; the
existing particle-distribution regression checks that input contract.
The factor 6 is the existing sphere surface-to-volume factor and lambda
is the existing surface-area correction. No coefficient or boundary
condition is being fitted here.

| Variable | Expected SI Units | Actual Units provided in the code |
|---|---|---|
| r, f, i, n, application index | dimensionless integer | MATLAB indices |
| d, delta_d | m | m |
| b, lambda, class_fraction, rock_fraction | dimensionless | dimensionless |
| a | m^-b | m^-b |
| rho_rock | kg m^-3 | g m^-3; needed with the gram-based mass inventory |
| M_added, mineral inventory | kg m^-2 ground | g m^-2 ground |
| psd | kg m^-3, per diameter and ground area | g m^-3; inherited [g/m] comments omit the ground-area normalization |
| SSA | m^2 kg^-1 | m^2 g^-1 |
| MM_min | kg mol^-1 | g umol^-1, from conv_mol = 1e6 |
| Wr | mol m^-2 rock s^-1 | umol m^-2 rock d^-1; inherited comments saying mol or seconds do not describe the internal basis |
| dtw | s | d; 1/(24*t_res), with t_res = 6 in this campaign |
| d_shrink | m s^-1 | m d^-1 |
| d_loss | m | m |
| SA | m^2 rock / m^2 ground | same area ratio; inherited [m2] comments omit the ground-area normalization |

The gram, micromole and day bases cancel consistently in d_shrink and
d_loss. The comment discrepancies above are documentation issues, not a
change to those equations. Storage and retrieval preserve all five
histories in double precision.

A 35-year run has 306816 hourly T&C rows but 1840896 ten-minute ERW rows.
With seven classes and 35 applications, each double history is 3608156160
bytes. Keeping five original arrays alongside their populated blocks adds
about 18.04 GB per treated column. The previous hourly-size benchmark
underestimated this footprint by a factor of six. Initialization still
allocates the full arrays before the block store starts. Final assembly
also needs both the blocks and the output arrays for the column being
assembled, so this change does not establish a lower Slurm memory request.

## Edge cases checked

| Case | Expected result |
|---|---|
| The solver returns scalar zeros before the first application | Expand them across the stored row before indexing the cache by application |
| A final window has no writable rows | No extra block; a write reports `RowOutsideDay` |
| A replay makes the first positive last class zero or negative | Rescan agrees with the original full-history search |
| A previous-row request is zero or past the history end | Report `RowOutsideHistory`; no invalid indexing |

No physical division is added by the helper. The existing solver only
updates SSA and the PSD ratio for positive diameters. Single-application
arrays keep their existing two-dimensional form. Block allocation and
assembly keep the array type; the campaign inputs and saved histories are
double precision.

## Validation

The tests first failed when the helper retained the full arrays, then
again when scalar zeros were indexed by application. Both regressions
pass after their fixes. The latter cache update reads the assigned block
row; array-valued solver outputs keep the same values as before.

| Check | Result |
|---|---|
| Storage, overlaps, first-row replacement, bounds and v7.3 reload | Pass |
| Scalar zeros before the first application | Failed with an index error before the fix; passes after it |
| Checkpoint source dependencies and Linux paths | Pass after the final helper change |
| Post-legacy numerical audit | All 92 entries pass; final PASS sentinel recorded |
| Plot configuration regression | Pass in a fresh MATLAB process, exit 0, without WebGL or font warnings |
| Eight-day, twelve-column amendment test | All 25 checks pass, including water, carbon, nutrient and lateral transport closure |
| Fresh ten-day local baseline and driver comparison after the scalar fix | All 26 checks pass; every reported difference is zero |
| Basalt, calcite at 5 and 50 t/ha, and a run before the first application | Selected physical outputs and the full ERW state match k bit for bit |
| Real daily solves with an added inactive future application | Replay twice and v7.3 reload match the original daily routine in all four cases |
| Hash manifests, shell syntax and original-root gate | 55 shared pins and 376 SSP pins pass; l accepts itself and rejects k |
| Final-helper history profile, 262944 rows, seven classes, four applications | Full-array writes and searches 8.019521 s; blocks 0.379734 s; ratio 21.12 |

The profile measures history handling with the MATLAB profiler enabled.
It excludes solver time and the one-off conversion to blocks. Whole-run
speed still needs measuring.

The scalar fix was found during the audit; the affected helper,
fingerprint, portability and driver checks were rerun after it. The
combined audit process was stopped after the numerical PASS sentinel
because the plot export stopped advancing following WebGL and TeX-font
warnings. The plotting regression then passed in a fresh process without
those warnings. Its configuration and uncertainty-band checks use
synthetic fixtures, not campaign results.

The validation scripts, logs and profiler output are retained in
`_validation`, outside the source archive. The local runs use MATLAB
R2026a. They do not establish Bouchet R2025a compatibility, full 35-year
runtime, peak memory or an interrupted 35-year campaign's resume behavior.
The runbook requires a fresh R2025a baseline and validation with the same
SSP inputs before release. Existing j/k checkpoints and stamps cannot be
reused; the byte-for-byte original-root gate remains in place.
