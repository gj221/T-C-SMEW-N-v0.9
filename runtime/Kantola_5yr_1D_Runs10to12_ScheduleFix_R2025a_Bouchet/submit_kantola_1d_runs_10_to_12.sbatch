#!/usr/bin/env bash
#SBATCH --job-name=Kantola1Dfix
#SBATCH --partition=day
#SBATCH --constraint=cpumodel:9575f
#SBATCH --nodes=1
#SBATCH --ntasks=3
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=06:00:00
#SBATCH --requeue
#SBATCH --output=runtime/Kantola_5yr_1D_Runs10to12_ScheduleFix_R2025a_Bouchet/slurm_%j.out
#SBATCH --error=runtime/Kantola_5yr_1D_Runs10to12_ScheduleFix_R2025a_Bouchet/slurm_%j.err

set -euo pipefail
umask 027

repo_root="${SLURM_SUBMIT_DIR:?Submit this job from the repository root}"
if [[ ! -f "$repo_root/v2/run/Run_One_Subfield.m" ]]; then
    echo "Submit from the T&C-SMEW-N v0.9 repository root." >&2
    exit 2
fi
cd "$repo_root"

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
    'v2/T&C_Code/ERW_Application_Indices.m' \
    'v2/T&C_Code/ERW_Daily_Step_2D.m' \
    'v2/T&C_Code/MAIN_FRAME_BG_ERW_2D.m' \
    'v2/T&C_Code/Parameters/ERW_MOD_PARAM_excel.m' \
    'v2/run/Run_One_Subfield.m'

export TC_PARALLEL_WORKERS=0
export TC_CAMPAIGN_TAG='Kantola_5yr_1D_Runs10to12_ScheduleFix_R2025a_Bouchet'
export TC_NUMERICAL_IDENTITY='cpumodel=9575f;matlab_release=2025a;single_comp_thread=1;parallel_workers=0'
export TC_CKPT_MINUTES=30
unset TC_CKPT_FORCE_RESUME || true
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

scratch_parent="${TMPDIR:-/tmp/$USER}"
mkdir -p "$scratch_parent"
scratch_parent="$(cd "$scratch_parent" && pwd -P)"
scratch_root="$scratch_parent/tc_smewn_1d_${SLURM_JOB_ID}_${SLURM_RESTART_COUNT:-0}"
mkdir -p "$scratch_root"
cleanup_scratch() {
    case "$scratch_root" in
        "$scratch_parent"/tc_smewn_1d_*) rm -rf -- "$scratch_root" ;;
        *) echo "Refusing to remove unexpected scratch path: $scratch_root" >&2 ;;
    esac
}
trap cleanup_scratch EXIT

if [[ "${TC_REBUILD_EQUIVALENCE:-0}" == "1" ]]; then
    export MATLAB_PREFDIR="$scratch_root/prefs_baseline"
    export TMPDIR="$scratch_root/tmp_baseline"
    mkdir -p "$MATLAB_PREFDIR" "$TMPDIR"
    export TMP="$TMPDIR"
    export TEMP="$TMPDIR"
    export TC_CASE=5
    export TC_SCENARIO=EQUIV_1D_TMP
    export TC_NDAYS=10
    srun --cpu-bind=cores matlab -singleCompThread -batch \
        "assert(strcmp(version('-release'),'2025a'),'The equivalence baseline requires MATLAB R2025a.'); cd(fullfile(pwd,'v2','run')); ERW_run"
    exit 0
fi

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
        "addpath(fullfile(pwd,'runtime','Kantola_5yr_1D_Runs10to12_ScheduleFix_R2025a_Bouchet')); Validate_Bouchet_1D_Campaign"
    exit 0
fi

if [[ "${SLURM_JOB_NUM_NODES:-0}" != "1" || \
      "${SLURM_NTASKS:-0}" != "3" || \
      "${SLURM_CPUS_PER_TASK:-0}" != "1" ]]; then
    echo "Production requires one node, three tasks and one CPU per task." >&2
    exit 5
fi

campaign='Kantola_5yr_1D_Runs10to12_ScheduleFix_R2025a_Bouchet'
results_dir="$repo_root/v2/SIM/Subfield_12/$campaign"
logs_dir="$repo_root/v2/logs/Subfield_12/$campaign"
mkdir -p "$results_dir" "$logs_dir"

pids=()
launched_runs=()
# Whether an existing MAT may stand is a validation question, not a
# filename question, so the decision is left to the MATLAB wrapper: it
# compares the saved provenance with the stamp this job passed and exits
# cleanly only when the two agree.
for run_id in 10 11 12; do
    run_scratch="$scratch_root/run_${run_id}"
    mkdir -p "$run_scratch/prefs" "$run_scratch/tmp"
    (
        export TC_RUN_ID="$run_id"
        export MATLAB_PREFDIR="$run_scratch/prefs"
        export TMPDIR="$run_scratch/tmp"
        export TMP="$TMPDIR"
        export TEMP="$TMPDIR"
        srun --exclusive --nodes=1 --ntasks=1 --cpus-per-task=1 \
            --cpu-bind=cores \
            --output="$logs_dir/srun_run${run_id}_%j.out" \
            --error="$logs_dir/srun_run${run_id}_%j.err" \
            matlab -singleCompThread -batch \
            "addpath(fullfile(pwd,'runtime','$campaign')); Run_Bouchet_Kantola_1D_Task"
    ) &
    pids+=("$!")
    launched_runs+=("$run_id")
done

failed=0
for k in "${!pids[@]}"; do
    pid="${pids[$k]}"
    if ! wait "$pid"; then
        echo "Run ${launched_runs[$k]} failed." >&2
        failed=1
    fi
done
if [[ "$failed" -ne 0 ]]; then
    exit 6
fi
echo "All requested 1-D treatments completed."
