# Kantola five-year 2-D runs 1-9 on Bouchet

This campaign runs the nine retained 2-D configurations from `Run_Matrix_12`:

| Runs | Layout | Treatments |
|---|---|---|
| 1-3 | whole field | control, basalt, calcite |
| 4-6 | cropped margin | control, basalt, calcite |
| 7-9 | prairie margin | control, basalt, calcite |

Each Slurm array element runs one configuration. Every run uses a 2 x 2 main grid, 12 columns, 20 lateral edges, 1826 days from 2016-01-01, a 300 mm active biogeochemical depth and `AppeloAlpha = 1.3`.

Run 9 applies calcite to the 1100 m2 prairie margin only. The 2500 m2 cropped center receives no direct calcite. Each prescribed 5000 g m^-2 event therefore adds 5.5 t to the margin, with four events adding 22 t over the simulation. This is a high-dose model treatment, not a field recommendation. Compare run 9 with run 7 to isolate calcite under the same prairie vegetation. Compare run 9 with run 6 to isolate prairie vegetation under the same margin-only calcite placement.

The job reserves 6 CPUs per run. MATLAB uses a six-worker thread pool only for the independent daily ERW chemistry calculations. Daily biology, hourly column physics, lateral routing and state write-back remain serial. BLAS and OpenMP are restricted to one thread to prevent oversubscription.

Rebuild the ten-day 1-D equivalence baseline under Bouchet R2025a after
installing a code patch:

```bash
sbatch --partition=day --array=1 --time=01:00:00 --cpus-per-task=1 --mem=16G --job-name=Kantola2Dbase --export=ALL,TC_REBUILD_EQUIVALENCE=1 runtime/Kantola_5yr_2D_Runs01to09_DenitPH7p8_Alpha1p3_Bouchet/submit_kantola_2d_runs_01_to_09.sbatch
```

Run the R2025a validation job only after that job completes:

```bash
sbatch --partition=day --array=1 --time=04:00:00 --job-name=Kantola2Dcheck --export=ALL,TC_VALIDATE_ONLY=1 runtime/Kantola_5yr_2D_Runs01to09_DenitPH7p8_Alpha1p3_Bouchet/submit_kantola_2d_runs_01_to_09.sbatch
```

Do not submit the production array unless this job ends with `PASS: Bouchet R2025a parallel 2-D validation`. Validation writes a content and numerical-platform stamp. Production stops if the stamp is missing or stale.

Submit from the repository root:

```bash
sbatch runtime/Kantola_5yr_2D_Runs01to09_DenitPH7p8_Alpha1p3_Bouchet/submit_kantola_2d_runs_01_to_09.sbatch
```

Production uses Bouchet's `week` partition with a four-day walltime and the
`cpumodel:9575f` constraint. The job records the node, CPU details and source
hashes before MATLAB starts. A
matching checkpoint is loaded when an unfinished array task is submitted
again.

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

Bouchet currently provides MATLAB R2025a. Override the module only if the cluster module name changes:

```bash
MATLAB_MODULE=MATLAB/2025a sbatch runtime/Kantola_5yr_2D_Runs01to09_DenitPH7p8_Alpha1p3_Bouchet/submit_kantola_2d_runs_01_to_09.sbatch
```

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
