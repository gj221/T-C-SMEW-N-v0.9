#!/usr/bin/env bash
#SBATCH --job-name=Base2D
#SBATCH --array=1-9%9
#SBATCH --partition=week
#SBATCH --constraint=cpumodel:9575f
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=48G
#SBATCH --time=4-00:00:00
#SBATCH --requeue
#SBATCH --output=runtime/Base_2D/slurm_%A_%a.out
#SBATCH --error=runtime/Base_2D/slurm_%A_%a.err

set -euo pipefail
umask 027

repo_root="${SLURM_SUBMIT_DIR:?Submit this job from the repository root}"
if [[ ! -f "$repo_root/v2/run/Run_One_Subfield.m" ]]; then
    echo "Submit from the T&C-SMEW-N v0.9 repository root." >&2
    exit 2
fi

# Only the campaign controller may submit this file: it performs the
# preflight checks, pins the checkout and builds the dependency graph.
if [[ "${TC_ALL12_CONTROLLER_TOKEN:-}" != "CEC13_ALL12_V1" ]]; then
    echo "Submit through runtime/Base_1D/submit_all_12.sh, not sbatch." >&2
    exit 7
fi
# The equivalence baseline is a fingerprinted dependency shared with the
# finished CEC 13 sweep; rebuilding it would change every stamp.
if [[ "${TC_REBUILD_EQUIVALENCE:-0}" != "0" ]]; then
    echo "This campaign never rebuilds the equivalence baseline." >&2
    exit 7
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
    'v2/T&C_Code/Parameters/MOD_PARAM_Kantola_Maize_Soy.m' \
    'v2/T&C_Code/Inputs/ERW_MOD_PARAM_parameters.xlsx' \
    'v2/run/Run_One_Subfield.m' \
    'runtime/Base_1D/Campaign_Parameters.m'

export TC_PARALLEL_WORKERS=6
# A calcite-case campaign (TC_CALCITE_CASE) runs the calcite treatments only,
# under its own tag; the controller sets its array to runs 3, 6 and 9.
case_suffix=''
if [[ -n "${TC_CALCITE_CASE:-}" ]]; then
    case_suffix="_case${TC_CALCITE_CASE}"
fi
export TC_CAMPAIGN_TAG="Base_2D${case_suffix}"
export TC_NUMERICAL_IDENTITY='cpumodel=9575f;matlab_release=2025a;single_comp_thread=1;parallel_workers=6'
export TC_CKPT_MINUTES=10
unset TC_CKPT_FORCE_RESUME || true
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

scratch_parent="${TMPDIR:-/tmp/$USER}"
mkdir -p "$scratch_parent"
scratch_parent="$(cd "$scratch_parent" && pwd -P)"
scratch_dir="$scratch_parent/tc_smewn_cal2d_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}_${SLURM_RESTART_COUNT:-0}"
mkdir -p "$scratch_dir/prefs" "$scratch_dir/tmp"
cleanup_scratch() {
    case "$scratch_dir" in
        "$scratch_parent"/tc_smewn_cal2d_*) rm -rf -- "$scratch_dir" ;;
        *) echo "Refusing to remove unexpected scratch path: $scratch_dir" >&2 ;;
    esac
}
trap cleanup_scratch EXIT

export MATLAB_PREFDIR="$scratch_dir/prefs"
export TMPDIR="$scratch_dir/tmp"
export TMP="$scratch_dir/tmp"
export TEMP="$scratch_dir/tmp"

cd "$repo_root"

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
    "addpath(fullfile(pwd,'runtime','Base_2D')); Run_Base_2D_Task"
