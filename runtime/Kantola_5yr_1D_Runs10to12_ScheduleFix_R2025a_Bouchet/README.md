# Corrected Kantola five-year 1-D runs 10-12 on Bouchet

This campaign reruns control, basalt and calcite after correcting the ERW
application-state index. The old files should be kept as archive evidence but
not used for the treatment comparison.

Production starts three separate MATLAB processes on one physical Bouchet
node. Each process gets one CPU. The `--mem=24G` request is a single limit
shared by all three steps, not a per-process reservation, so one greedy
process can starve the other two; raise the request rather than assuming 8 GB
each. MATLAB, BLAS and OpenMP use one numerical thread. This avoids the mixed
CPU and thread configuration used by the earlier array while keeping the
treatments parallel.

Each treatment takes a lock directory under `v2/logs/Subfield_12/<campaign>/
locks/` before it touches the status file, checkpoint or result path, so a
duplicate submission fails immediately instead of racing. A requeued job keeps
its job id and reclaims its own orphaned lock, so `--requeue` and the rolling
checkpoint still work unattended. A lock left by any other job stops the run:
remove it by hand only after confirming with `squeue` that no job is still
running that treatment.

Run these commands from the repository root. First rebuild the ten-day
equivalence baseline with Bouchet R2025a:

```bash
sbatch --ntasks=1 --cpus-per-task=1 --mem=16G --time=01:00:00 --export=ALL,TC_REBUILD_EQUIVALENCE=1 runtime/Kantola_5yr_1D_Runs10to12_ScheduleFix_R2025a_Bouchet/submit_kantola_1d_runs_10_to_12.sbatch
```

After it completes, run the full audit and the 315-day same-node identity
check:

```bash
sbatch --ntasks=1 --cpus-per-task=6 --mem=32G --time=04:00:00 --export=ALL,TC_VALIDATE_ONLY=1 runtime/Kantola_5yr_1D_Runs10to12_ScheduleFix_R2025a_Bouchet/submit_kantola_1d_runs_10_to_12.sbatch
```

Do not submit the production job unless the validation log ends with:

```text
PASS: Bouchet R2025a corrected 1-D validation
```

The production job is not a Slurm array. All treatments share one node:

```bash
sbatch runtime/Kantola_5yr_1D_Runs10to12_ScheduleFix_R2025a_Bouchet/submit_kantola_1d_runs_10_to_12.sbatch
```

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
sacct -j JOB_ID --format=JobID,State,ExitCode,Elapsed,MaxRSS
```
