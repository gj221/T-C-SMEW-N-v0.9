#!/usr/bin/env bash
#SBATCH --job-name=Sweep1D
#SBATCH --partition=day
#SBATCH --constraint=cpumodel:9575f
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=05:00:00
#SBATCH --array=1
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --output=runtime/Calibration_Sweep_1D/slurm_%A_%a.out
#SBATCH --error=runtime/Calibration_Sweep_1D/slurm_%A_%a.err

set -Eeuo pipefail
umask 027

if [[ "${TC_SWEEP_CONTROLLER_TOKEN:-}" != "CEC13_STAGED_V4" ]]; then
    echo 'Submit this batch file through submit_sweep.sh.' >&2
    exit 6
fi

repo_root="${SLURM_SUBMIT_DIR:?Submit this job from the repository root}"
if [[ ! -f "$repo_root/v2/run/Run_One_Subfield.m" ]]; then
    echo "Submit from the T&C-SMEW-N v0.9 repository root." >&2
    exit 2
fi
cd "$repo_root"

campaign_name='Calibration_Sweep_1D'
campaign_dir="runtime/$campaign_name"

module purge
module load "${MATLAB_MODULE:-MATLAB/2025a}"

node_record="$(scontrol show node "$SLURMD_NODENAME")"
printf '%s\n' "$node_record"
if [[ "$node_record" != *"cpumodel:9575f"* ]]; then
    echo "The job is not running on the required cpumodel:9575f node." >&2
    exit 3
fi
hostname
lscpu
module list
sha256sum \
    'v2/T&C_Code/BIOGEOCHEMISTRY_DYNAMIC3.m' \
    'v2/T&C_Code/MAIN_FRAME_BG_ERW_2D.m' \
    'v2/run/Run_One_Subfield.m' \
    'v2/run/Run_All_12_Parallel.m' \
    "$campaign_dir/Sweep_Cases.m" \
    "$campaign_dir/Run_Sweep_Task.m" \
    "$campaign_dir/Validate_Sweep.m" \
    "$campaign_dir/Summarize_Sweep.m"

export TC_PARALLEL_WORKERS=0
export TC_CAMPAIGN_TAG="$campaign_name"
export TC_NUMERICAL_IDENTITY='cpumodel=9575f;matlab_release=2025a;single_comp_thread=1;parallel_workers=0'
export TC_CKPT_MINUTES=30
unset TC_CKPT_FORCE_RESUME || true
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

sweep_scope="${TC_SWEEP_SCOPE:-pre_control}"
summary_scope="${TC_SWEEP_SUMMARY_SCOPE:-pre_control}"
case "$sweep_scope" in
    full|pre_control) ;;
    *) echo 'TC_SWEEP_SCOPE must be full or pre_control.' >&2; exit 6 ;;
esac
case "$summary_scope" in
    full|pre_control) ;;
    *) echo 'TC_SWEEP_SUMMARY_SCOPE must be full or pre_control.' >&2; exit 6 ;;
esac

scratch_parent="${TMPDIR:-/tmp/$USER}"
mkdir -p "$scratch_parent"
scratch_parent="$(cd "$scratch_parent" && pwd -P)"
scratch_dir="$scratch_parent/tc_smewn_joint1d_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID:-0}_${SLURM_RESTART_COUNT:-0}"
mkdir -p "$scratch_dir/prefs" "$scratch_dir/tmp"
cleanup_scratch() {
    case "$scratch_dir" in
        "$scratch_parent"/tc_smewn_joint1d_*) rm -rf -- "$scratch_dir" ;;
        *) echo "Refusing to remove unexpected scratch path: $scratch_dir" >&2 ;;
    esac
}
trap cleanup_scratch EXIT
export MATLAB_PREFDIR="$scratch_dir/prefs"
export TMPDIR="$scratch_dir/tmp"
export TMP="$TMPDIR"
export TEMP="$TMPDIR"

if [[ "${TC_REBUILD_EQUIVALENCE:-0}" == "1" ]]; then
    export TC_CKPT_MINUTES=0
    export TC_CASE=5
    export TC_SCENARIO=EQUIV_1D_TMP
    export TC_NDAYS=10
    srun --kill-on-bad-exit=1 --cpu-bind=cores \
        matlab -singleCompThread -batch \
        "assert(strcmp(version('-release'),'2025a'),'The equivalence baseline requires MATLAB R2025a.'); cd(fullfile(pwd,'v2','run')); ERW_run"
    exit 0
fi

if [[ "${TC_VALIDATE_ONLY:-0}" == "1" ]]; then
    if [[ "${SLURM_NTASKS:-0}" != "1" || "${SLURM_CPUS_PER_TASK:-0}" -lt 6 ]]; then
        echo "Validation requires --ntasks=1 and --cpus-per-task=6." >&2
        exit 4
    fi
    export TC_CKPT_MINUTES=0
    srun --kill-on-bad-exit=1 --cpu-bind=cores \
        matlab -singleCompThread -batch \
        "addpath(fullfile(pwd,'$campaign_dir')); Validate_Sweep"
    exit 0
fi

if [[ "${TC_SUMMARIZE_ONLY:-0}" == "1" ]]; then
    if [[ "${SLURM_NTASKS:-0}" != "1" || "${SLURM_CPUS_PER_TASK:-0}" != "1" ]]; then
        echo "Summarization requires one task and one CPU." >&2
        exit 5
    fi
    export TC_CKPT_MINUTES=0
    srun --kill-on-bad-exit=1 --cpu-bind=cores \
        matlab -singleCompThread -batch \
        "addpath(fullfile(pwd,'$campaign_dir')); Summarize_Sweep('$summary_scope')"
    exit 0
fi

sweep_layout="${TC_SWEEP_LAYOUT:-independent}"
if [[ "${SLURM_JOB_NUM_NODES:-0}" != "1" || \
      "${SLURM_CPUS_PER_TASK:-0}" != "1" || \
      -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    echo "Production requires one node, one CPU per task and an array task ID." >&2
    exit 6
fi

if [[ "$sweep_scope" == "pre_control" && \
      ( "$sweep_layout" != "independent" || \
        "${SLURM_ARRAY_TASK_ID:-0}" != "1" ) ]]; then
    echo 'pre_control requires independent array task 1.' >&2
    exit 6
fi

if [[ "$sweep_layout" == "independent" ]]; then
    if [[ "${SLURM_NTASKS:-0}" != "1" ]]; then
        echo "Independent layout requires one Slurm task." >&2
        exit 6
    fi
    unset TC_SWEEP_CASE_ID || true
    srun --kill-on-bad-exit=1 --cpu-bind=cores \
        env TC_SWEEP_LAYOUT=independent TC_SWEEP_SCOPE="$sweep_scope" \
        matlab -singleCompThread -batch \
        "addpath(fullfile(pwd,'$campaign_dir')); Run_Sweep_Task"
    exit 0
fi

if [[ "$sweep_scope" != "full" || "$sweep_layout" != "packed" || \
      "${SLURM_NTASKS:-0}" != "3" ]]; then
    echo "Packed layout requires exactly three Slurm tasks." >&2
    exit 6
fi

run_packed_smoke_rank() {
    local rank="$1"
    local smoke_dir="$2"
    srun --exclusive --nodes=1 --ntasks=1 --cpus-per-task=1 \
        --mem=8G --kill-on-bad-exit=0 --cpu-bind=cores \
        bash -c '
            set -Eeuo pipefail
            smoke_dir="$1"
            rank="$2"
            awk '\''/^Cpus_allowed_list:/ {print $2}'\'' /proc/self/status > \
                "$smoke_dir/cpu_$rank"
            : > "$smoke_dir/ready_$rank"
            for ((attempt=0; attempt<150; attempt++)); do
                shopt -s nullglob
                ready=("$smoke_dir"/ready_*)
                if (( ${#ready[@]} == 3 )); then exit 0; fi
                sleep 0.2
            done
            echo "Packed Slurm steps did not overlap." >&2
            exit 8
        ' bash "$smoke_dir" "$rank"
}

if [[ "${TC_PACKED_SMOKE:-0}" == "1" ]]; then
    smoke_dir="$scratch_dir/packed_smoke"
    mkdir -p "$smoke_dir"
    packed_smoke_pids=()
    for rank in 0 1 2; do
        run_packed_smoke_rank "$rank" "$smoke_dir" &
        packed_smoke_pids+=("$!")
    done
    packed_smoke_failed=0
    for packed_smoke_pid in "${packed_smoke_pids[@]}"; do
        if ! wait "$packed_smoke_pid"; then
            packed_smoke_failed=1
        fi
    done
    mapfile -t packed_cpu_masks < <(sort -u "$smoke_dir"/cpu_*)
    if (( packed_smoke_failed != 0 || ${#packed_cpu_masks[@]} != 3 )); then
        echo "Packed Slurm resource smoke failed." >&2
        exit 8
    fi
    printf 'PASS: three packed steps overlapped on distinct CPU masks: %s\n' \
        "${packed_cpu_masks[*]}"
elif [[ "${TC_PACKED_SMOKE:-0}" != "0" ]]; then
    echo "TC_PACKED_SMOKE must be 0 or 1." >&2
    exit 6
fi

run_packed_case() {
    local case_id="$1"
    local case_scratch="$scratch_dir/case_${case_id}"
    mkdir -p "$case_scratch/prefs" "$case_scratch/tmp"
    srun --exclusive --nodes=1 --ntasks=1 --cpus-per-task=1 \
        --mem=8G \
        --kill-on-bad-exit=0 --cpu-bind=cores \
        env TC_SWEEP_LAYOUT=packed TC_SWEEP_CASE_ID="$case_id" \
        MATLAB_PREFDIR="$case_scratch/prefs" TMPDIR="$case_scratch/tmp" \
        TMP="$case_scratch/tmp" TEMP="$case_scratch/tmp" \
        matlab -singleCompThread -batch \
        "addpath(fullfile(pwd,'$campaign_dir')); Run_Sweep_Task"
}

first_case=$(( (SLURM_ARRAY_TASK_ID - 1) * 3 + 1 ))
packed_pids=()
for case_id in "$first_case" "$((first_case + 1))" "$((first_case + 2))"; do
    run_packed_case "$case_id" &
    packed_pids+=("$!")
done
packed_failed=0
for packed_pid in "${packed_pids[@]}"; do
    if ! wait "$packed_pid"; then
        packed_failed=1
    fi
done
if (( packed_failed != 0 )); then
    echo "At least one packed trajectory failed; completed siblings remain resumable." >&2
    exit 7
fi
