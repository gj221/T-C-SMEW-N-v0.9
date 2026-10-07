#!/usr/bin/env bash
#SBATCH --job-name=WorkerScale
#SBATCH --partition=day
#SBATCH --constraint=cpumodel:9575f
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --array=1-3
#SBATCH --output=runtime/Worker_Scaling_Test/slurm_%A_%a.out
#SBATCH --error=runtime/Worker_Scaling_Test/slurm_%A_%a.err
#
# Three arms of a 365-day 2-D run, to size the 35-year campaigns:
#   1   6 workers, checkpoints every 10 min   (production settings)
#   2  12 workers, checkpoints every 10 min
#   3   6 workers, checkpoints disabled
#
# Every arm reserves 12 CPUs so that worker count is the only difference
# between arms 1 and 2. Arm 1 therefore leaves six cores idle; that is the
# point, not an oversight.
#
# Submit from the checkout root:
#   sbatch runtime/Worker_Scaling_Test/submit_worker_scaling_test.sbatch
#
# Afterwards:
#   sacct -j <jobid> --format=JobID%18,Elapsed,MaxRSS,State
# Arm 1 vs 2 gives the worker speedup, arm 1 vs 3 the checkpoint cost, and
# MaxRSS across all three (with the known five-year figures, 17.3 GB control
# and 24.2 GB basalt) gives a second point for sizing --mem at 35 years.
#
# Diagnostic only: no controller token, no validation stamp, no provenance
# gate. Delete this directory and v2/SIM/Worker_Scaling_Test once the
# numbers are recorded.

set -euo pipefail
umask 027

repo_root="${SLURM_SUBMIT_DIR:?Submit this job from the repository root}"
if [[ ! -f "$repo_root/v2/run/Run_One_Subfield.m" ]]; then
    echo "Submit from the T&C-SMEW-N v0.9 checkout root." >&2
    exit 2
fi
cd "$repo_root"

module purge
module load "${MATLAB_MODULE:-MATLAB/2025a}"

# The whole point is a timing comparable with production, so refuse a node
# that is not the one production pins.
node_record="$(scontrol show node "$SLURMD_NODENAME")"
if [[ "$node_record" != *"cpumodel:9575f"* ]]; then
    echo "Not running on a cpumodel:9575f node; the timing would not be comparable." >&2
    exit 3
fi
hostname
module list

# Worker count and checkpoint interval are set per arm inside MATLAB. Clear
# the environment fallbacks so nothing else can decide them.
unset TC_PARALLEL_WORKERS TC_CKPT_MINUTES TC_CKPT_FORCE_RESUME || true
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

# MATLAB deadlocks on a shared preferences directory, and an NFS prefdir is
# worse. Keep prefs and temp on node-local scratch, per arm.
scratch_parent="${TMPDIR:-/tmp/$USER}"
mkdir -p "$scratch_parent"
scratch_parent="$(cd "$scratch_parent" && pwd -P)"
scratch_root="$scratch_parent/tc_wscale_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$scratch_root/prefs" "$scratch_root/tmp"
cleanup_scratch() {
    case "$scratch_root" in
        "$scratch_parent"/tc_wscale_*) rm -rf -- "$scratch_root" ;;
        *) echo "Refusing to remove unexpected scratch path: $scratch_root" >&2 ;;
    esac
}
trap cleanup_scratch EXIT

export MATLAB_PREFDIR="$scratch_root/prefs"
export TMPDIR="$scratch_root/tmp"
export TMP="$TMPDIR"
export TEMP="$TMPDIR"

srun --cpu-bind=cores matlab -singleCompThread -batch \
    "addpath(fullfile(pwd,'runtime','Worker_Scaling_Test')); Run_Worker_Scaling_Test"

echo "Arm ${SLURM_ARRAY_TASK_ID} finished."
