# Check, archive, compare, rerun: the Bouchet twelve-run campaign

The order matters. The submission controller refuses to start a campaign in a
root that already holds results, so the archive step is not housekeeping that
can be deferred — it is the step that unblocks the rerun.

Everything below assumes the Yale VPN is up. Bouchet is not reachable without
it.

## 1. Check the run that just finished

Log in and look at the accounting record rather than at the queue. A finished
job leaves the queue, so `squeue` alone cannot tell a completed campaign from
one that was cancelled.

```bash
ssh NETID@bouchet.ycrc.yale.edu
```

```bash
sacct -u "$USER" --starttime now-3days \
  --name=KantolaBase,Kantola1Dcheck,Kantola2Dcheck,Kantola1Dfix,Kantola2D \
  --format=JobID%20,JobName%18,State,ExitCode,Submit,Start,End,Elapsed,MaxRSS,NodeList%20
```

Read it as the dependency graph in `submit_bouchet_all_12.sh`:

| JobName | Stage | Wall clock requested |
|---|---|---|
| `KantolaBase` | R2025a ten-day equivalence baseline | 1 h |
| `Kantola1Dcheck` | 1-D validation, four 315-day identity trajectories | 4 h |
| `Kantola2Dcheck` | 2-D validation, array element 1 | 4 h |
| `Kantola1Dfix` | production runs 10-12, three steps on one node | 6 h |
| `Kantola2D` | production runs 1-9, nine-element array | 23:30 h |

A roughly five-hour elapsed time is the 1-D production allocation
(`Kantola1Dfix`) or a validation job that used most of its four hours. Confirm
which from the `JobName` and `Elapsed` columns before assuming the whole
campaign is done — the nine 2-D elements sit in a separate, much longer job.

Every stage must show `COMPLETED` with exit code `0:0`. A production job that
shows `TIMEOUT`, `NODE_FAIL` or `CANCELLED` may still have written a usable
checkpoint; resubmit it rather than archiving a partial result.

Then confirm what the runs themselves recorded, from the deployment root:

```bash
cd /path/to/deployment/root
c1d=Kantola_5yr_1D_Runs10to12_ScheduleFix_R2025a_Bouchet
c2d=Kantola_5yr_2D_Runs01to09_DenitPH7p8_Alpha1p3_Bouchet

grep -H state v2/logs/Subfield_12/$c1d/run_*.status
grep -H state v2/logs/Subfield_12/$c2d/run_*.status
ls -1 v2/SIM/Subfield_12/$c1d/TC_2D_run1[0-2]_*__1826d.mat | wc -l   # expect 3
ls -1 v2/SIM/Subfield_12/$c2d/TC_2D_run0[1-9]_*__1826d.mat | wc -l   # expect 9
```

The validation logs carry the two gate lines the production jobs depend on.
Both must be present:

```bash
grep -r 'PASS: Bouchet R2025a corrected 1-D validation'  runtime/$c1d/slurm_*.out
grep -r 'PASS: Bouchet R2025a parallel 2-D validation'   runtime/$c2d/slurm_*.out
```

Check the conservation report each production run writes after it finishes —
ERW mass residual against gross lateral transport, and soil-atmosphere
inorganic carbon exchange:

```bash
grep -iE 'residual|conservation|PASS|FAIL' v2/logs/Subfield_12/$c2d/worker_run*.log | tail -40
```

Leftover lock directories mean a run died without releasing its lock. Remove
one only after `squeue` shows nothing running that element:

```bash
ls v2/logs/Subfield_12/$c1d/locks v2/logs/Subfield_12/$c2d/locks 2>/dev/null
```

## 2. Archive the outputs to OneDrive

Run this on the workstation in Git Bash, from the repository root. It refuses
to archive while a Kantola job is still queued or running, and refuses if any
of the twelve results is missing or any status file does not say
`state=complete`.

```bash
cd "/c/Users/Gregory/OneDrive/Yale/TC-SMEW-N/T&C-SMEW-N v0.9"
bash runtime/pull_bouchet_campaign.sh \
  --netid NETID \
  --remote-root '/path/to/deployment/root' \
  --label denitph7p8_alpha1p3
```

What it does, in order: opens one ssh control connection so Duo is answered
once; runs the completeness checks above; reports the transfer size and asks
before moving anything; takes SHA-256 sums on Bouchet; streams the results and
logs down through `tar` over ssh, leaving checkpoints and locks behind;
verifies every file against its Bouchet checksum; and only then moves the
finished folder into OneDrive.

The staging folder lives outside OneDrive on purpose. A partially transferred
multi-gigabyte MAT must never start syncing, and a failed transfer must not
leave half an archive in the tree that the next campaign will read.

The archive lands beside the existing ones as
`v2/SIM/Subfield_12/archive_<UTC stamp>_<label>/`:

```text
RUN_INFO.txt        provenance in the existing house format
MANIFEST.sha256     Bouchet-side checksums, verified after transfer
results/            all twelve TC_2D_run*__1826d.mat, flattened
logs/1D_runs10to12/ status files and worker logs
logs/2D_runs01to09/
slurm/              sacct.txt and the Slurm startup logs
figures/            filled in by step 3
```

The twelve MATs are flattened into one `results/` folder because
`Plot_1D_Then_All_12` reads a single `ResultDirectory`. Campaign provenance is
not lost: it stays in `RUN_INFO.txt` and in the two separate log trees.

Size is the thing to watch. The 1-D 1826-day results ran 160-275 MB each in
the August campaign, and the 2-D results are larger. The script prints the
total before it transfers anything. If OneDrive should not carry the full set,
use `--logs-only` to bring back logs and Slurm output alone and leave the MATs
on Bouchet until they are needed.

Do not delete the deployment root on Bouchet until `sha256sum -c` has passed
locally. The script does that verification for you and stops on the first
mismatch.

## 3. Compare the twelve runs

`Plot_1D_Then_All_12` is the comparison. It runs the vertical 1-D weathering
and nitrogen figure first, then produces summary and time-dynamics figures for
each metric, comparing control, basalt and calcite within every layout it
finds.

```bash
cd "/c/Users/Gregory/OneDrive/Yale/TC-SMEW-N/T&C-SMEW-N v0.9"
ARCHIVE='v2/SIM/Subfield_12/archive_<stamp>_<label>'
matlab -batch "cd('v2/plotting'); NDays=1826; ResultDirectory='$PWD/$ARCHIVE/results'; FigureDirectory='$PWD/$ARCHIVE/figures'; Plot_1D_Then_All_12"
```

`NDays=1826` is required, and the files must have been renamed on the way in.
`Run_One_Subfield` tags a result with every option that was set, so Bouchet
writes:

```text
TC_2D_run01_control__whole_field__grid2x2__alpha1p3__1826d__from20160101.mat
```

`Plot_1D_Then_All_12` builds only `<treatment>__<layout>__1826d`, and raises
`MissingResult` on anything else. `pull_bouchet_campaign.sh` trims the extra
tags after verifying checksums against the delivered names, and records the
mapping in `RENAMED.txt`. Downloading by hand means doing that trim yourself:

```bash
cd "$ARCHIVE/results" && for f in TC_2D_run*.mat; do n="$(sed -E 's/__grid[0-9]+x[0-9]+//; s/__alpha[0-9p]+//; s/__from[0-9]{8}//' <<<"$f")"; [ "$n" = "$f" ] || mv -n "$f" "$n"; done
```

Verify checksums *before* renaming — the manifest carries the delivered names.

The design being compared is three treatments crossed with four layouts:

|  | whole_field | margin_crop | margin_prairie | vertical_1d |
|---|---|---|---|---|
| control | 1 | 4 | 7 | 10 |
| basalt | 2 | 5 | 8 | 11 |
| calcite | 3 | 6 | 9 | 12 |

Runs 1-9 are 2-D with a 2x2 main grid, 12 columns, 20 lateral edges and a 5 m
margin. Runs 10-12 are a single 3600 m2 column with no subdomains and no
lateral transport, which is why they are the cleanest treatment contrast and
why the script plots them first.

Two contrasts worth reading deliberately, from the 2-D campaign README: run 9
against run 7 isolates calcite under the same prairie vegetation; run 9 against
run 6 isolates prairie vegetation under the same margin-only calcite placement.
Run 9's 5000 g m^-2 events land on the 1100 m2 margin only, 22 t over the
simulation — a high-dose model treatment, not a field rate.

`DryRun` prints the plan without plotting, which is the fast way to confirm the
script sees all twelve files:

```bash
matlab -batch "cd('v2/plotting'); NDays=1826; ResultDirectory='$PWD/$ARCHIVE/results'; DryRun=true; Plot_1D_Then_All_12"
```

Supporting figures, same `ResultDirectory` convention: `Plot_Kantola_Benchmark`
against the chamber observations, `Plot_Nflux_Diagnostics` for the nitrogen
fluxes, `Plot_Biogeo_Pools` for the pools.

Stage 2 is skipped with a warning if 2-D runs are missing, so the script is
still usable while the nine-element array is in flight. Stage 1 raises if a 1-D
result is missing.

## 4. Rerun with the updated nitrogen response functions

The package in `Downloads` already carries them. It was built on 2026-09-01 at
14:50 from this working tree, its checksum verifies, and nothing in the
manifest has been edited since:

```text
tc-smewn-v09-energy-farm-p-return-bouchet.tar.gz
tc-smewn-v09-energy-farm-p-return-bouchet.tar.gz.sha256.txt
```

If any manifest file changes before you deploy, rebuild it with the GNU tar
block in the repository `README.md` and regenerate the sidecar. A source edit
made while a job is running is caught by the post-run fingerprint check rather
than silently accepted, so a stale package fails loudly rather than quietly.

Copy both files to Bouchet, into the **parent** of a new, empty root. Do not
extract over the root you just archived. That tree holds results under an older
campaign fingerprint, and the controller stops rather than mixing them:

```bash
scp tc-smewn-v09-energy-farm-p-return-bouchet.tar.gz* NETID@bouchet.ycrc.yale.edu:~/project/
```

On Bouchet, confirm the account first, then submit from inside the new root:

```bash
slurm_checkup.sh
mkdir -p ~/project/tc-smewn-v09-nresponse && cd ~/project/tc-smewn-v09-nresponse
set -euo pipefail
(cd .. && sha256sum -c tc-smewn-v09-energy-farm-p-return-bouchet.tar.gz.sha256.txt)
tar --no-same-owner -xzf ../tc-smewn-v09-energy-farm-p-return-bouchet.tar.gz
test -r 'v2/T&C_Code/Inputs/Data_Kantola.mat'
test -r 'v2/T&C_Code/Inputs/Data_Kantola_BG.mat'
test -r 'v2/T&C_Code/Tests/run_post_legacy_equation_audit_tests.m'
test -r 'v2/POST_LEGACY_EQUATION_AUDIT.md'
MATLAB_MODULE=MATLAB/2025a TC_BOUCHET_ACCOUNT=pi_np363 bash runtime/submit_bouchet_all_12.sh
```

The controller submits all five jobs held, builds the dependency graph, then
releases them. It cancels what it already submitted if any step fails. Confirm
`pi_np363` against what `slurm_checkup.sh` reports rather than trusting the
value written here.

Watch it with:

```bash
squeue -u "$USER"
squeue -u "$USER" -n Kantola2D
jobstats <job-id>
```

Then return to step 1 for the new campaign. Keep the old archive: the corrected
nitrogen response is a treatment-comparison change, and the previous archive is
the before-side of that comparison.
