# Energy Farm (Kantola) five-year 2-D runs 1-9 on Bouchet

This Energy Farm campaign runs the nine retained 2-D configurations from
`Run_Matrix_12`:

| Runs | Layout | Treatments |
|---|---|---|
| 1-3 | whole field | control, basalt, calcite |
| 4-6 | cropped margin | control, basalt, calcite |
| 7-9 | prairie margin | control, basalt, calcite |

Each Slurm array element runs one configuration. Every run uses a 2 x 2 main grid, 12 columns, 20 lateral edges, 1826 days from 2016-01-01, a 300 mm active biogeochemical depth and `AppeloAlpha = 1.3`.

The `DenitPH7p8` tag identifies the configured gross-denitrification pH
sensitivity. It does not mean that total denitrification is a nitrate
Michaelis-Menten rate. The packaged solver uses the Dickinson and legacy
T&C-SMEW first-order nitrate sink, with temperature, moisture, reconstructed
bulk-pH and volumetric DOC factors. `Kc = 17.5 g C m^-3`; the volumetric
`Kn = 82.5 g N m^-3` is retained for provenance but is inactive. The 7.8 pH
plateau adapts the Blanc-Betes response shape using the soil-specific evidence
from Simek et al. and is not a universal optimum. N2 and
denitrification-derived N2O sum to the gross sink. Volatilization retains the
legacy first-order NH4 term with no direct pH multiplier. Its pH effect is
indirect through NH4 availability. The equations, units and citations are in
[`v2/POST_LEGACY_EQUATION_AUDIT.md`](../../v2/POST_LEGACY_EQUATION_AUDIT.md),
which is included in the deployment manifest.

Run 9 applies calcite to the 1100 m2 prairie margin only. The 2500 m2 cropped center receives no direct calcite. Each prescribed 5000 g m^-2 event therefore adds 5.5 t to the margin, with four events adding 22 t over the simulation. This is a high-dose model treatment, not a field recommendation. Compare run 9 with run 7 to isolate calcite under the same prairie vegetation. Compare run 9 with run 6 to isolate prairie vegetation under the same margin-only calcite placement.

The job reserves 6 CPUs per run. MATLAB uses a six-worker thread pool only for the independent daily ERW chemistry calculations. Daily biology, hourly column physics, lateral routing and state write-back remain serial. BLAS and OpenMP are restricted to one thread to prevent oversubscription.

Do not call this batch file directly for the twelve-run campaign. From a new,
empty deployment root, run the coordinated controller after `slurm_checkup.sh`
confirms the account:

```bash
MATLAB_MODULE=MATLAB/2025a TC_BOUCHET_ACCOUNT=pi_np363 bash runtime/submit_bouchet_all_12.sh
```

The controller rebuilds the ten-day 1-D equivalence baseline under R2025a.
After it succeeds, the 1-D and 2-D validation allocations may run concurrently.
The 1-D allocation uses four MATLAB process workers for the 315-day control,
basalt, calcite and repeated-calcite trajectories. The separate 2-D validation
uses array element 1 with 6 CPUs and 32 GB. Production becomes eligible only
if both jobs succeed. The 2-D validation must end with `PASS: Bouchet R2025a
parallel 2-D validation`. It writes a content and numerical-platform stamp.
Production stops if the stamp is missing or stale.

`TC_VALIDATE_ONLY` is set to `1` only for validation. It stays unset for
production.

After each production run, the wrapper checks the sum of absolute ERW mass residuals against gross lateral transport. It checks soil-atmosphere inorganic carbon exchange the same way. Opposing errors on different hours cannot cancel and pass the run.

Production uses Bouchet's `day` partition with a 23:30:00 walltime, 32 GB per
array element and the `cpumodel:9575f` constraint. The job records the node,
CPU details and source hashes before MATLAB starts. A
matching checkpoint is loaded when an unfinished array task is submitted
again.

The controller submits this array and the three-task 1-D allocation held, builds
the complete dependency graph, then releases it. Both production allocations
depend on both R2025a validations and become eligible together. Slurm controls
their actual start times. Runs 10 to 12 start as steps only after the 1-D
allocation starts. This array permits nine active elements, so all nine 2-D
runs and all three 1-D steps can run concurrently when Bouchet has the
resources. Eligibility does not guarantee that all twelve start together.

A final MAT only stops a resubmission after the run has finished. While a run
is live, exclusion comes from the lock directory each array element takes
under `v2/logs/Subfield_12/<campaign>/locks/`, so a second submission of the
same element fails immediately instead of sharing its status file, checkpoint
and result path. A requeued job keeps its job id and reclaims its own orphaned
lock, so `--requeue` and the rolling checkpoint still work unattended. A lock
left by any other job stops the run: remove it by hand only after confirming
with `squeue` that no job is still running that element.

An existing final MAT is accepted only when its saved provenance carries the
same validation fingerprint, column count and horizon as the job now
starting. A result produced under an earlier validation stops the campaign
rather than being silently mixed with newer runs. Move it out of the campaign
directory, or rerun every treatment under one validation.

This campaign requires MATLAB R2025a. The coordinated controller checks that
the module exists before creating any job and stops if it is unavailable.

Results are written to:

```text
v2/SIM/Subfield_12/Kantola_5yr_2D_Runs01to09_DenitPH7p8_Alpha1p3_Bouchet
```

Worker logs and atomic status files are written to:

```text
v2/logs/Subfield_12/Kantola_5yr_2D_Runs01to09_DenitPH7p8_Alpha1p3_Bouchet
```

Slurm startup logs remain in this runtime directory. Periodic checkpoints are written every ten minutes beneath the result directory. Normal resubmission can load a matching checkpoint. Forced resume remains disabled.

Useful checks:

```bash
squeue -u "$USER" -n Kantola2D
jobstats <job-id>
```

Cluster guidance:

- [YCRC MATLAB jobs](https://docs.ycrc.yale.edu/clusters-at-yale/guides/matlab/)
- [YCRC CPU and memory requests](https://docs.ycrc.yale.edu/clusters-at-yale/job-scheduling/resource-requests/)
- [YCRC job arrays](https://docs.ycrc.yale.edu/clusters-at-yale/job-scheduling/dsq/)
