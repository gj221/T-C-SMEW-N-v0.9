#!/usr/bin/env bash
#SBATCH --job-name=Kantola2D
#SBATCH --array=1-9%9
#SBATCH --partition=day
#SBATCH --constraint=cpumodel:9575f
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=32G
#SBATCH --time=23:30:00
#SBATCH --requeue
#SBATCH --output=runtime/Kantola_5yr_2D_Runs01to09_DenitPH7p8_Alpha1p3_Bouchet/slurm_%A_%a.out
#SBATCH --error=runtime/Kantola_5yr_2D_Runs01to09_DenitPH7p8_Alpha1p3_Bouchet/slurm_%A_%a.err

set -euo pipefail
umask 027

repo_root="${SLURM_SUBMIT_DIR:?Submit this job from the repository root}"
if [[ ! -f "$repo_root/v2/run/Run_One_Subfield.m" ]]; then
    echo "Submit from the T&C-SMEW-N v0.9 repository root." >&2
    exit 2
fi

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

export TC_PARALLEL_WORKERS=6
export TC_CAMPAIGN_TAG='Kantola_5yr_2D_Runs01to09_DenitPH7p8_Alpha1p3_Bouchet'
export TC_NUMERICAL_IDENTITY='cpumodel=9575f;matlab_release=2025a;single_comp_thread=1;parallel_workers=6'
export TC_CKPT_MINUTES=10
unset TC_CKPT_FORCE_RESUME || true
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

scratch_parent="${TMPDIR:-/tmp/$USER}"
mkdir -p "$scratch_parent"
scratch_parent="$(cd "$scratch_parent" && pwd -P)"
scratch_dir="$scratch_parent/tc_smewn_2d_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}_${SLURM_RESTART_COUNT:-0}"
mkdir -p "$scratch_dir/prefs" "$scratch_dir/tmp"
cleanup_scratch() {
    case "$scratch_dir" in
        "$scratch_parent"/tc_smewn_2d_*) rm -rf -- "$scratch_dir" ;;
        *) echo "Refusing to remove unexpected scratch path: $scratch_dir" >&2 ;;
    esac
}
trap cleanup_scratch EXIT

export MATLAB_PREFDIR="$scratch_dir/prefs"
export TMPDIR="$scratch_dir/tmp"
export TMP="$scratch_dir/tmp"
export TEMP="$scratch_dir/tmp"

cd "$repo_root"
if [[ "${TC_REBUILD_EQUIVALENCE:-0}" == "1" ]]; then
    export TC_CASE=5
    export TC_SCENARIO=EQUIV_1D_TMP
    export TC_NDAYS=10
    unset TC_VALIDATE_ONLY || true
    srun --cpu-bind=cores matlab -singleCompThread -batch \
        "assert(strcmp(version('-release'),'2025a'),'The equivalence baseline requires MATLAB R2025a.'); cd(fullfile(pwd,'v2','run')); ERW_run"
    exit 0
fi

# An unqualified srun inherits the allocation, so a submission overriding
# --ntasks would start that many MATLAB processes for one array element,
# all writing the same status file, checkpoint and result path.
if [[ "${SLURM_JOB_NUM_NODES:-0}" != "1" || \
      "${SLURM_NTASKS:-0}" != "1" || \
      "${SLURM_CPUS_PER_TASK:-0}" -lt 6 ]]; then
    echo "Validation and production require one node, one task and at least six CPUs." >&2
    exit 5
fi

srun --nodes=1 --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" \
    --cpu-bind=cores matlab -singleCompThread -batch \
    "addpath(fullfile(pwd,'runtime','Kantola_5yr_2D_Runs01to09_DenitPH7p8_Alpha1p3_Bouchet')); Run_Bouchet_Kantola_2D_Array_Task"
