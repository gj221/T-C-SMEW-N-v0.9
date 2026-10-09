#!/usr/bin/env bash
#SBATCH --job-name=FieldSSP2D
#SBATCH --array=1-9
#SBATCH --partition=day
#SBATCH --requeue
#SBATCH --signal=B:USR1@1800
#SBATCH --open-mode=append
#SBATCH --constraint=cpumodel:9575f
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=192G
#SBATCH --time=1-00:00:00
#SBATCH --output=runtime/Field_SSP_2D/slurm_%A_%a.out
#SBATCH --error=runtime/Field_SSP_2D/slurm_%A_%a.err
#
# Runs 1-9 of one 35-year SSP campaign on a field grid saved by
# Build_Field_Grid_Topology (TC_FIELD_TOPOLOGY): the scenario forcing
# (TC_SSP), rock every 10 November to 2050, optionally the 5 t/ha calcite
# case (TC_CALCITE_CASE). The grid's file name, the scenario and the case
# enter the campaign tag, so every grid and scenario writes to its own tree.
#
# Submit through runtime/Field_SSP_2D/submit_field_ssp_2d.sh, which sets
# --mem, --cpus-per-task, --partition, --time and the array from the
# campaign's choices; the directives above are only the defaults. The
# day partition stops a job at 24 hours: a run requeues itself half an
# hour before (requeue_near_limit.sh), at most TC_FIELD_RESTARTS times,
# and continues from its two-hourly checkpoint under the same job id and
# run lock. The validation element is never requeued.
#
# TC_FIELD_WORKERS workers (one CPU each) share out the squares' daily
# vegetation and biogeochemistry, hourly physics and daily ERW chemistry.

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
if [[ "${TC_FIELD_SSP_CONTROLLER_TOKEN:-}" != "FIELD_SSP2D_V1" ]]; then
    echo "Submit through runtime/Field_SSP_2D/submit_field_ssp_2d.sh, not sbatch." >&2
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
: "${TC_SSP:?TC_SSP must be ssp126, ssp245 or ssp585}"
case "$TC_SSP" in ssp126|ssp245|ssp585) ;; *) echo "Bad TC_SSP: $TC_SSP" >&2; exit 2 ;; esac
if [[ ! -f "v2/T&C_Code/Inputs/Data_Kantola_${TC_SSP}.mat" ]]; then
    echo "The scenario forcing v2/T&C_Code/Inputs/Data_Kantola_${TC_SSP}.mat is missing." >&2
    exit 6
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
    'v2/T&C_Code/ERW_Daily_Step_2D.m' \
    'v2/run/Run_One_Subfield.m' \
    'runtime/Base_1D/Campaign_Parameters.m' \
    "v2/T&C_Code/Inputs/Data_Kantola_${TC_SSP}.mat" \
    "$TC_FIELD_TOPOLOGY"

workers="${TC_FIELD_WORKERS:-32}"
if [[ ! "$workers" =~ ^[1-9][0-9]*$ ]] || (( workers < 2 || workers > 64 )); then
    echo "TC_FIELD_WORKERS must be a whole number from 2 to 64." >&2
    exit 2
fi
restarts="${TC_FIELD_RESTARTS:-12}"
if [[ ! "$restarts" =~ ^[0-9]+$ ]] || (( restarts > 30 )); then
    echo "TC_FIELD_RESTARTS must be a whole number from 0 to 30." >&2
    exit 2
fi
kernel="${TC_ERW_KERNEL:-matlab}"
case "$kernel" in
    matlab|fast) ;;
    *) echo "TC_ERW_KERNEL must be matlab or fast for a field run, not '$kernel'." >&2; exit 2 ;;
esac
export TC_ERW_KERNEL="$kernel"
case_suffix=''
if [[ -n "${TC_CALCITE_CASE:-}" ]]; then
    case_suffix="_case${TC_CALCITE_CASE}"
fi
export TC_FIELD_WORKERS="$workers"
export TC_PARALLEL_WORKERS="$workers"
topology_name="$(basename "$TC_FIELD_TOPOLOGY" .mat)"
export TC_CAMPAIGN_TAG="Field_SSP_2D_35yr_${TC_SSP}${case_suffix}__${topology_name}"
export TC_NUMERICAL_IDENTITY="cpumodel=9575f;matlab_release=2025a;single_comp_thread=1;parallel_workers=${workers}"
export TC_CKPT_MINUTES=120
unset TC_CKPT_FORCE_RESUME || true
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

# Each piece of each element keeps its own MATLAB preferences and temporary
# files, so the thread pools of elements sharing a node never share them.
scratch_parent="${TMPDIR:-/tmp/$USER}"
mkdir -p "$scratch_parent"
scratch_parent="$(cd "$scratch_parent" && pwd -P)"
scratch_dir="$scratch_parent/tc_smewn_fieldssp_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}_${SLURM_RESTART_COUNT:-0}"
mkdir -p "$scratch_dir/prefs" "$scratch_dir/tmp"
cleanup_scratch() {
    case "$scratch_dir" in
        "$scratch_parent"/tc_smewn_fieldssp_*) rm -rf -- "$scratch_dir" ;;
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
      "${SLURM_CPUS_PER_TASK:-0}" -lt "$workers" ]]; then
    echo "Validation and production require one node, one task and at least $workers CPUs." >&2
    exit 5
fi

task="addpath(fullfile(pwd,'runtime','Field_SSP_2D')); Run_Field_SSP_2D_Task"
if [[ "${TC_VALIDATE_ONLY:-}" == "1" ]]; then
    # The time-limit signal is for production pieces; a validation runs to
    # its end or to its limit.
    trap '' USR1
    srun --nodes=1 --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" \
        --cpu-bind=cores matlab -singleCompThread -batch "$task"
else
    source runtime/Uncertainty_1D/requeue_near_limit.sh
    max_restarts="$restarts"
    run_requeueing srun --nodes=1 --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" \
        --cpu-bind=cores matlab -singleCompThread -batch "$task"
fi

echo "Field SSP 2-D ${TC_CAMPAIGN_TAG} array element ${SLURM_ARRAY_TASK_ID} finished."
