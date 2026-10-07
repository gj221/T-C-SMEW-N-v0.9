#!/usr/bin/env bash
# Submit the lower and upper sets of the Kantola uncertainty ensemble for
# one horizon and calcite case: a 1-D array over runs 10-12 and a 2-D array
# over runs 1-9 (run 12 and runs 3, 6, 9 for a calcite case). The original
# set is not run: it is the Base or SSP campaign result in the checkout
# TC_UC_ORIGINAL_ROOT, whose model must be byte for byte this one's.
#
# Usage, from the checkout root:
#   TC_BOUCHET_ACCOUNT=pi_np363 TC_UC_ORIGINAL_ROOT=/path/to/original \
#   bash runtime/Uncertainty_1D/submit_uncertainty.sh
# Without TC_SSP the ensemble is the 5-year validation on the legacy
# forcing (originals Base_1D and Base_2D); TC_SSP=ssp126|ssp245|ssp585 with
# TC_HORIZON=35yr is that projection (originals SSP_1D_35yr_<ssp> and
# SSP_2D_35yr_<ssp>). Add TC_CALCITE_CASE=8 for the 5 t/ha arm.
# TC_UC_2D_MEM overrides the 256G asked for each 35-year 2-D element, the
# SSP 2-D campaign's figure (see its batch file). TC_UC_WITH_ORIGINAL=1 runs
# the original set here too, first in each array, under the same checks.
# TC_UC_DIMS=1d or 2d submits that array alone (default all: both).

set -euo pipefail

repo_root="$(pwd -P)"
campaign_dir='runtime/Uncertainty_1D'
sbatch_1d="$campaign_dir/submit_uncertainty.sbatch"
sbatch_2d='runtime/Uncertainty_2D/submit_uncertainty_2d.sbatch'
token='UNCERTAINTY_CAMPAIGN_V2'
run_members=2   # lower and upper (Uncertainty_Members.m); the original is reused
if [[ ! -f "$repo_root/v2/run/Run_One_Subfield.m" || ! -f "$repo_root/$sbatch_1d" || ! -f "$repo_root/$sbatch_2d" ]]; then
    echo "Run this from the T&C-SMEW-N v0.9 checkout root." >&2
    exit 2
fi
if [[ ! -f "$repo_root/v2/SIM/EQUIV_1D_TMP/TC_ERW_results_Kantola_crop_control.mat" ]]; then
    echo "The equivalence baseline is missing; copy it from the original's checkout." >&2
    exit 6
fi
: "${TC_BOUCHET_ACCOUNT:?Set TC_BOUCHET_ACCOUNT to the pi account reported by slurm_checkup.sh}"
: "${TC_UC_ORIGINAL_ROOT:?Set TC_UC_ORIGINAL_ROOT to the checkout whose Base and SSP campaigns hold the original set}"

ssp="${TC_SSP:-}"
horizon="${TC_HORIZON:-}"
if [[ -z "$ssp" ]]; then
    if [[ -n "$horizon" && "$horizon" != 5yr ]]; then
        echo "Without TC_SSP the ensemble is the 5-year validation; TC_HORIZON=$horizon is not it." >&2
        exit 2
    fi
    horizon='5yr'
else
    case "$ssp" in ssp126|ssp245|ssp585) ;; *) echo "Bad TC_SSP: $ssp" >&2; exit 2 ;; esac
    if [[ "$horizon" != 35yr ]]; then
        echo "A scenario ensemble is the 35-year projection; TC_HORIZON=$horizon is not it." >&2
        exit 2
    fi
    if [[ ! -f "$repo_root/v2/T&C_Code/Inputs/Data_Kantola_${ssp}.mat" ]]; then
        echo "The scenario forcing v2/T&C_Code/Inputs/Data_Kantola_${ssp}.mat is missing." >&2
        exit 6
    fi
fi
with_original="${TC_UC_WITH_ORIGINAL:-}"
case "$with_original" in
    ''|0) ;;
    1) run_members=3 ;;   # the original as well, as Uncertainty_Campaign.m orders them
    *) echo "TC_UC_WITH_ORIGINAL must be unset, 0 or 1 (got '$with_original')." >&2; exit 2 ;;
esac
dims="${TC_UC_DIMS:-all}"
case "$dims" in
    all|1d|2d) ;;
    *) echo "TC_UC_DIMS must be all, 1d or 2d (got '$dims')." >&2; exit 2 ;;
esac
calcite_case="${TC_CALCITE_CASE:-}"
case_suffix=''
runs_1d=3; runs_2d=9
if [[ -n "$calcite_case" ]]; then
    if [[ ! "$calcite_case" =~ ^[1-9][0-9]*$ ]]; then
        echo "Bad TC_CALCITE_CASE: $calcite_case" >&2; exit 2
    fi
    case_suffix="_case${calcite_case}"
    runs_1d=1; runs_2d=3
fi
# Per element: the 1-D memory is the SSP 1-D campaign's per-run figure; a
# 2-D element has the memory and workers of its original's campaign
# (Base_2D at 5 years, SSP_2D at 35), whose worker count the task
# requires. Everything runs on the day partition, whose per-user limits
# (1500 CPUs, 20 TB) are far wider than the week partition's (96 CPUs,
# 1.5 TB); a run past the time limit requeues itself and continues from
# its checkpoint.
case "$horizon" in
    5yr)  time_1d=12:00:00;   mem_1d=12G;  workers_2d=6;  mem_2d=48G ;;
    35yr) time_1d=1-00:00:00; mem_1d=56G;  workers_2d=12; mem_2d="${TC_UC_2D_MEM:-256G}" ;;
esac
# A 35-year 2-D element holds a result of about 11 GB for a control and
# more for a treatment, tens of GB for whole-field basalt, a checkpoint of
# about the same and a temporary copy while each checkpoint is written;
# eighteen of them at once can need a terabyte. The check below is only a
# floor: a scratch copy of the checkout is the place for the 35-year sets
# (RUNBOOK section 0v).
minimum_free_gib=100
[[ "$horizon" == 35yr ]] && minimum_free_gib=600
if ! available_kib="$(LC_ALL=C df -Pk "$repo_root" | awk 'NR == 2 {print $4}')" || \
        [[ ! "$available_kib" =~ ^[0-9]+$ ]]; then
    echo "Cannot read the free space of $repo_root." >&2; exit 4
fi
if (( available_kib < minimum_free_gib * 1024 * 1024 )); then
    echo "Only $(( available_kib / 1024 / 1024 )) GiB free; this call needs $minimum_free_gib GiB. No jobs were submitted." >&2
    exit 4
fi

bash "$campaign_dir/same_model_code.sh" "$TC_UC_ORIGINAL_ROOT" || {
    echo 'No jobs were submitted.' >&2; exit 6; }

suffix="${ssp:-legacy}_${horizon}${case_suffix}"
tag_2d="Uncertainty_2D${ssp:+_${horizon}_${ssp}}${case_suffix}"
if ! active_jobs="$(squeue -h -u "$USER" --name="UC1D_$suffix,UC2D_$suffix" -o '%A|%j|%T')"; then
    echo 'Cannot check for active ensemble jobs.' >&2; exit 5
fi
if [[ -n "$active_jobs" ]]; then
    echo "Ensemble jobs for $suffix are still queued or running. No jobs were submitted." >&2
    printf '%s\n' "$active_jobs" >&2
    exit 5
fi
stale_locks="$(find "$repo_root/v2/logs/Subfield_12/$tag_2d" -type d -name 'run_*.lock' 2>/dev/null || true)"
if [[ -n "$stale_locks" ]]; then
    echo 'Stale run locks remain from killed jobs. Remove them and rerun:' >&2
    while IFS= read -r lock_dir; do
        printf '  rm -r -- %q\n' "$lock_dir" >&2
    done <<< "$stale_locks"
    exit 5
fi

unset TC_VALIDATE_ONLY TC_CKPT_FORCE_RESUME
for submit_env in "${!SBATCH_@}"; do unset "$submit_env"; done
export_common="ALL,TC_UC_CONTROLLER_TOKEN=$token,TC_SSP=$ssp,TC_HORIZON=$horizon,TC_CALCITE_CASE=$calcite_case,TC_UC_WITH_ORIGINAL=$with_original"
job_1d=''
job_2d=''
if [[ "$dims" != 2d ]]; then
    job_1d="$(sbatch --parsable --account="$TC_BOUCHET_ACCOUNT" \
        --partition=day --time="$time_1d" --mem="$mem_1d" --cpus-per-task=1 \
        --array="1-$(( run_members * runs_1d ))" --job-name="UC1D_$suffix" \
        --export="$export_common" "$sbatch_1d")"
    job_1d="${job_1d%%;*}"
fi
if [[ "$dims" != 1d ]]; then
    if ! job_2d="$(sbatch --parsable --account="$TC_BOUCHET_ACCOUNT" \
        --partition=day --time=1-00:00:00 --mem="$mem_2d" --cpus-per-task="$workers_2d" \
        --array="1-$(( run_members * runs_2d ))" --job-name="UC2D_$suffix" \
        --export="$export_common" "$sbatch_2d")"; then
        if [[ -n "$job_1d" ]]; then
            echo "The 2-D submission failed; cancelling the 1-D array $job_1d." >&2
            scancel "$job_1d"
        else
            echo 'The 2-D submission failed.' >&2
        fi
        exit 1
    fi
    job_2d="${job_2d%%;*}"
fi

printf 'Uncertainty ensemble %s (original from %s)\n' "$suffix" "$TC_UC_ORIGINAL_ROOT"
monitor=''
if [[ -n "$job_1d" ]]; then
    printf '  1-D %s: %d elements, day, %s, %s each\n' "$job_1d" "$(( run_members * runs_1d ))" "$time_1d" "$mem_1d"
    monitor="$job_1d"
fi
if [[ -n "$job_2d" ]]; then
    printf '  2-D %s: %d elements, day, 1-00:00:00, %s and %d CPUs each\n' "$job_2d" "$(( run_members * runs_2d ))" "$mem_2d" "$workers_2d"
    monitor="${monitor:+$monitor,}$job_2d"
fi
printf 'Monitor: sacct -j %s -X --format=JobID%%16,JobName%%24,State%%12,Elapsed,MaxRSS\n' "$monitor"
