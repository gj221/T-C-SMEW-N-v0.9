#!/usr/bin/env bash
#SBATCH --job-name=SSP1D
#SBATCH --partition=week
#SBATCH --constraint=cpumodel:9575f
#SBATCH --nodes=1
#SBATCH --ntasks=3
#SBATCH --cpus-per-task=1
#SBATCH --mem=192G
#SBATCH --time=1-12:00:00
#SBATCH --requeue
#SBATCH --output=runtime/SSP_1D/slurm_%j.out
#SBATCH --error=runtime/SSP_1D/slurm_%j.err
#
# Runs 10-12 of one SSP campaign, three one-column treatments in parallel on
# one node. TC_SSP and TC_HORIZON select scenario and run length.
#
# The directives above size the 35-year case; the controller overrides
# --partition, --time, --ntasks, --cpus-per-task and --mem for validation and
# for 5yr. A one-column run holds roughly a twelfth of the 2-D state, so the
# five-year 8 G per task scales to well under 56 G at 35 years; the per-task
# figure below is set from TC_HORIZON with generous headroom. MaxRSS was
# never recorded for the five-year 1-D runs -- only that they completed
# inside 8 G -- so 56 G assumes the worst case that they used all of it and
# then scales by seven. Cheap insurance on a 2251 GB node.

set -euo pipefail
umask 027

repo_root="${SLURM_SUBMIT_DIR:?Submit this job from the repository root}"
if [[ ! -f "$repo_root/v2/run/Run_One_Subfield.m" ]]; then
    echo "Submit from the T&C-SMEW-N v0.9 checkout root." >&2
    exit 2
fi
cd "$repo_root"

if [[ "${TC_SSP_CONTROLLER_TOKEN:-}" != "SSP_CAMPAIGN_V1" ]]; then
    echo "Submit through runtime/SSP_1D/submit_ssp.sh, not sbatch." >&2
    exit 7
fi
if [[ "${TC_REBUILD_EQUIVALENCE:-0}" != "0" ]]; then
    echo "This campaign never rebuilds the equivalence baseline." >&2
    exit 7
fi

: "${TC_SSP:?TC_SSP must be ssp126, ssp245 or ssp585}"
: "${TC_HORIZON:?TC_HORIZON must be 5yr or 35yr}"
case "$TC_SSP" in ssp126|ssp245|ssp585) ;; *) echo "Bad TC_SSP: $TC_SSP" >&2; exit 2 ;; esac
case "$TC_HORIZON" in
    5yr)  run_mem='12G' ;;
    35yr) run_mem='56G' ;;
    *) echo "Bad TC_HORIZON: $TC_HORIZON" >&2; exit 2 ;;
esac

module purge
module load "${MATLAB_MODULE:-MATLAB/2025a}"

node_record="$(scontrol show node "$SLURMD_NODENAME")"
printf '%s\n' "$node_record"
if [[ "$node_record" != *"cpumodel:9575f"* ]]; then
    echo "The job is not running on the required cpumodel:9575f node." >&2
    exit 3
fi
hostname
module list

export TC_PARALLEL_WORKERS=0
# A calcite-case campaign (TC_CALCITE_CASE) runs the calcite treatment only,
# under its own tag.
case_suffix=''
runs=(10 11 12)
if [[ -n "${TC_CALCITE_CASE:-}" ]]; then
    case_suffix="_case${TC_CALCITE_CASE}"
    runs=(12)
fi
export TC_CAMPAIGN_TAG="SSP_1D_${TC_HORIZON}_${TC_SSP}${case_suffix}"
export TC_NUMERICAL_IDENTITY='cpumodel=9575f;matlab_release=2025a;single_comp_thread=1;parallel_workers=0'
export TC_CKPT_MINUTES=30
unset TC_CKPT_FORCE_RESUME || true
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

scratch_parent="${TMPDIR:-/tmp/$USER}"
mkdir -p "$scratch_parent"
scratch_parent="$(cd "$scratch_parent" && pwd -P)"
scratch_root="$scratch_parent/tc_smewn_ssp1d_${SLURM_JOB_ID}_${SLURM_RESTART_COUNT:-0}"
mkdir -p "$scratch_root"
cleanup_scratch() {
    case "$scratch_root" in
        "$scratch_parent"/tc_smewn_ssp1d_*) rm -rf -- "$scratch_root" ;;
        *) echo "Refusing to remove unexpected scratch path: $scratch_root" >&2 ;;
    esac
}
trap cleanup_scratch EXIT

campaign='SSP_1D'

# The validator needs a single task with at least six CPUs; production needs
# three tasks of one CPU. The controller asks for the right shape, and these
# checks refuse anything else rather than failing deep inside MATLAB.
if [[ "${TC_VALIDATE_ONLY:-0}" == "1" ]]; then
    if [[ "${SLURM_NTASKS:-0}" != "1" || "${SLURM_CPUS_PER_TASK:-0}" -lt 6 ]]; then
        echo "Validation requires --ntasks=1 and --cpus-per-task=6." >&2
        exit 4
    fi
    export MATLAB_PREFDIR="$scratch_root/prefs_validation"
    export TMPDIR="$scratch_root/tmp_validation"
    mkdir -p "$MATLAB_PREFDIR" "$TMPDIR"
    export TMP="$TMPDIR"
    export TEMP="$TMPDIR"
    srun --cpu-bind=cores matlab -singleCompThread -batch \
        "addpath(fullfile(pwd,'runtime','$campaign')); Validate_SSP_1D"
    exit 0
fi

if [[ "${SLURM_JOB_NUM_NODES:-0}" != "1" || \
      "${SLURM_NTASKS:-0}" != "${#runs[@]}" || \
      "${SLURM_CPUS_PER_TASK:-0}" != "1" ]]; then
    echo "Production requires one node, ${#runs[@]} task(s) and one CPU per task." >&2
    exit 5
fi

logs_dir="$repo_root/v2/logs/Subfield_12/$TC_CAMPAIGN_TAG"
results_dir="$repo_root/v2/SIM/Subfield_12/$TC_CAMPAIGN_TAG"
mkdir -p "$logs_dir" "$results_dir"

pids=()
launched_runs=()
# Whether an existing MAT may stand is a validation question, not a filename
# question, so the decision is left to the MATLAB wrapper.
for run_id in "${runs[@]}"; do
    run_scratch="$scratch_root/run_${run_id}"
    mkdir -p "$run_scratch/prefs" "$run_scratch/tmp"
    (
        export TC_RUN_ID="$run_id"
        export MATLAB_PREFDIR="$run_scratch/prefs"
        export TMPDIR="$run_scratch/tmp"
        export TMP="$TMPDIR"
        export TEMP="$TMPDIR"
        srun --exclusive --nodes=1 --ntasks=1 --cpus-per-task=1 \
            --mem="$run_mem" \
            --cpu-bind=cores \
            --output="$logs_dir/srun_run${run_id}_%j.out" \
            --error="$logs_dir/srun_run${run_id}_%j.err" \
            matlab -singleCompThread -batch \
            "addpath(fullfile(pwd,'runtime','$campaign')); Run_SSP_1D_Task"
    ) &
    pids+=("$!")
    launched_runs+=("$run_id")
done

failed=0
for k in "${!pids[@]}"; do
    if ! wait "${pids[$k]}"; then
        echo "Run ${launched_runs[$k]} failed." >&2
        failed=1
    fi
done
if [[ "$failed" -ne 0 ]]; then
    exit 6
fi
echo "SSP 1-D ${TC_SSP} ${TC_HORIZON}${case_suffix} treatments ${runs[*]} completed."
