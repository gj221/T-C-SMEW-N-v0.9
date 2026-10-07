#!/usr/bin/env bash
#SBATCH --job-name=Field2D
#SBATCH --array=1-9%9
#SBATCH --partition=week
#SBATCH --constraint=cpumodel:9575f
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=256G
#SBATCH --time=4-00:00:00
#SBATCH --requeue
#SBATCH --output=runtime/Field_2D/slurm_%A_%a.out
#SBATCH --error=runtime/Field_2D/slurm_%A_%a.err
#
# Runs 1-9 of the five-year Kantola campaign on a field grid saved by
# Build_Field_Grid_Topology (TC_FIELD_TOPOLOGY). The grid's file name enters
# the campaign tag, so every grid writes to its own tree.
#
# Memory and time scale with the grid's cell count, because the hourly
# T&C loop is serial over columns and every column carries its own state:
# the controller sets --mem and --time from TC_FIELD_MEM and TC_FIELD_TIME,
# so the directives above are only the defaults.
#
# Twelve workers for the daily ERW columns; see the task.

set -euo pipefail
umask 027

repo_root="${SLURM_SUBMIT_DIR:?Submit this job from the repository root}"
if [[ ! -f "$repo_root/v2/run/Run_One_Subfield.m" ]]; then
    echo "Submit from the T&C-SMEW-N v0.9 repository root." >&2
    exit 2
fi
cd "$repo_root"

# Only the campaign controller may submit this file: it runs the preflight
# and builds the dependency graph.
if [[ "${TC_FIELD_CONTROLLER_TOKEN:-}" != "FIELD2D_V1" ]]; then
    echo "Submit through runtime/Field_2D/submit_field_2d.sh, not sbatch." >&2
    exit 7
fi
# The equivalence baseline is a fingerprinted dependency; rebuilding it
# would change every stamp.
if [[ "${TC_REBUILD_EQUIVALENCE:-0}" != "0" ]]; then
    echo "This campaign never rebuilds the equivalence baseline." >&2
    exit 7
fi
: "${TC_FIELD_TOPOLOGY:?TC_FIELD_TOPOLOGY must name the grid file}"
if [[ ! -r "$TC_FIELD_TOPOLOGY" ]]; then
    echo "The grid file is missing or unreadable: $TC_FIELD_TOPOLOGY" >&2
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
    'v2/T&C_Code/MAIN_FRAME_BG_ERW_2D.m' \
    'v2/T&C_Code/Darcy_Lateral_Topology.m' \
    'v2/T&C_Code/Surface_Runoff_Routing.m' \
    'v2/run/Run_One_Subfield.m' \
    'runtime/Base_1D/Campaign_Parameters.m' \
    "$TC_FIELD_TOPOLOGY"

export TC_PARALLEL_WORKERS=12
topology_name="$(basename "$TC_FIELD_TOPOLOGY" .mat)"
export TC_CAMPAIGN_TAG="Field_2D__${topology_name}"
export TC_NUMERICAL_IDENTITY='cpumodel=9575f;matlab_release=2025a;single_comp_thread=1;parallel_workers=12'
export TC_CKPT_MINUTES=120
unset TC_CKPT_FORCE_RESUME || true
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

scratch_parent="${TMPDIR:-/tmp/$USER}"
mkdir -p "$scratch_parent"
scratch_parent="$(cd "$scratch_parent" && pwd -P)"
scratch_dir="$scratch_parent/tc_smewn_field2d_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}_${SLURM_RESTART_COUNT:-0}"
mkdir -p "$scratch_dir/prefs" "$scratch_dir/tmp"
cleanup_scratch() {
    case "$scratch_dir" in
        "$scratch_parent"/tc_smewn_field2d_*) rm -rf -- "$scratch_dir" ;;
        *) echo "Refusing to remove unexpected scratch path: $scratch_dir" >&2 ;;
    esac
}
trap cleanup_scratch EXIT

export MATLAB_PREFDIR="$scratch_dir/prefs"
export TMPDIR="$scratch_dir/tmp"
export TMP="$scratch_dir/tmp"
export TEMP="$scratch_dir/tmp"

# An unqualified srun inherits the allocation, so a submission overriding
# --ntasks would start that many MATLAB processes for one array element,
# all writing the same status file, checkpoint and result path.
if [[ "${SLURM_JOB_NUM_NODES:-0}" != "1" || \
      "${SLURM_NTASKS:-0}" != "1" || \
      "${SLURM_CPUS_PER_TASK:-0}" -lt 12 ]]; then
    echo "Validation and production require one node, one task and at least twelve CPUs." >&2
    exit 5
fi

srun --nodes=1 --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" \
    --cpu-bind=cores matlab -singleCompThread -batch \
    "addpath(fullfile(pwd,'runtime','Field_2D')); Run_Field_2D_Task"

echo "Field 2-D ${TC_CAMPAIGN_TAG} array element ${SLURM_ARRAY_TASK_ID} finished."
