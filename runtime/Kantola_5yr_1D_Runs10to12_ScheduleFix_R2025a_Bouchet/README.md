# Corrected Energy Farm (Kantola) five-year 1-D runs 10-12 on Bouchet

This Energy Farm campaign reruns control, basalt and calcite after correcting
the ERW application-state index. The old files should be kept as archive
evidence but not used for the treatment comparison.

The packaged solver also contains the separate denitrification correction
documented in
[`v2/POST_LEGACY_EQUATION_AUDIT.md`](../../v2/POST_LEGACY_EQUATION_AUDIT.md).
Gross denitrification is first-order in nitrate, with temperature, moisture,
bulk-pH and volumetric DOC factors. `Kc = 17.5 g C m^-3`; the volumetric
`Kn = 82.5 g N m^-3` remains inactive because this model has no DNDC
denitrifier-capacity state. The gross pH response reaches its configured
plateau at 7.8. This is a treatment sensitivity, not the pH 7.0 breakpoint in
the original Blanc-Betes equation or a universal optimum. Volatilization
retains `VOL = 8.6e-5*NH4` per day, with pH acting indirectly through the NH4
pool rather than through an explicit NH3 equilibrium. The deployment manifest
includes the audit so the archived release retains this rationale.

Production starts three separate MATLAB processes on one physical Bouchet
node. Each exclusive step gets one CPU and 8 GB from the 24 GB allocation.
The explicit step memory prevents Slurm from treating each `srun` as if it
needed the complete allocation and serializing the three runs. MATLAB, BLAS
and OpenMP use one numerical thread. In completed job 23970418, the largest
step used 3174016K. The 8 GB step limit retains more than 2.5 times that
measured peak.

Each treatment takes a lock directory under `v2/logs/Subfield_12/<campaign>/
locks/` before it touches the status file, checkpoint or result path, so a
duplicate submission fails immediately instead of racing. A requeued job keeps
its job id and reclaims its own orphaned lock, so `--requeue` and the rolling
checkpoint still work unattended. A lock left by any other job stops the run:
remove it by hand only after confirming with `squeue` that no job is still
running that treatment.

Do not call this batch file directly for the twelve-run campaign. From a new,
empty deployment root, run the coordinated controller after `slurm_checkup.sh`
confirms the account:

```bash
MATLAB_MODULE=MATLAB/2025a TC_BOUCHET_ACCOUNT=pi_np363 bash runtime/submit_bouchet_all_12.sh
```

The controller rebuilds the ten-day baseline under R2025a. After that succeeds,
the 1-D and 2-D validation allocations may run concurrently. The 1-D
allocation reserves 6 CPUs and 40 GB. Four MATLAB process workers run the
315-day control, basalt, calcite and repeated-calcite identity trajectories.
Production becomes eligible only if both validation jobs succeed. The 1-D
validation log must end with:

```text
PASS: Bouchet R2025a corrected 1-D validation
```

Production is not a Slurm array. All three treatments share one node. The
controller submits this allocation and the nine-element 2-D array held, builds
the complete dependency graph, then releases it. Both production jobs become
eligible together. The 1-D allocation then starts three exclusive MATLAB steps,
while the 2-D array permits up to nine active elements. Bouchet can therefore
run all twelve simulations at once when resources are available. Slurm does not
guarantee an identical start time.

Results are written to:

```text
v2/SIM/Subfield_12/Kantola_5yr_1D_Runs10to12_ScheduleFix_R2025a_Bouchet
```

Worker logs and status files are written to the matching directory under
`v2/logs/Subfield_12`. Validation writes a content and numerical-platform
stamp there. Production stops if that stamp is missing or stale. Checkpoints
are written every 30 minutes and removed
after a successful save. Resubmission leaves a completed MAT unchanged only
when its saved provenance carries the same validation fingerprint, column count
and horizon as the job now starting; a result from an earlier validation stops
the campaign instead of being mixed with newer runs. Unfinished runs resume
only when their code, input and forcing-window fingerprint still matches.
Forced resume remains disabled. Each run recomputes the campaign fingerprint
after it finishes, so a source edit made while the job was running is caught
rather than silently accepted.

Useful checks:

```bash
squeue -u "$USER" -n Kantola1Dfix
sacct -j JOB_ID --format=JobID,State,ExitCode,Start,End,Elapsed,MaxRSS
```
