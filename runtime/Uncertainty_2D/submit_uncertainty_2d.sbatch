#!/usr/bin/env bash
#SBATCH --job-name=Uncertainty2D
#SBATCH --requeue
#SBATCH --partition=day
#SBATCH --signal=B:USR1@1800
#SBATCH --open-mode=append
#SBATCH --constraint=cpumodel:9575f
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=256G
#SBATCH --time=1-00:00:00
#SBATCH --array=1-18
#SBATCH --output=runtime/Uncertainty_2D/slurm_%A_%a.out
#SBATCH --error=runtime/Uncertainty_2D/slurm_%A_%a.err
#
# One element of the Kantola 2-D uncertainty ensemble: the lower or upper
# set on one of runs 1-9. The directives above size a 35-year element as
# the SSP 2-D campaign sizes its memory and CPUs; submit through
# runtime/Uncertainty_1D/submit_uncertainty.sh, which sets memory, CPUs
# and the array for the horizon and calcite case. The workers are the
# CPUs reserved (6 for the 5-year ensemble, as Base_2D; 12 for the
# projections, as SSP_2D) and are declared in the numerical identity; the
# task refuses any other count. The day partition's per-user limits are
# far wider than the week partition's, and a 35-year run (about 43 hours)
# requeues itself before the 24-hour limit and continues from its
# checkpoint (requeue_near_limit.sh).

set -euo pipefail
umask 027

repo_root="${SLURM_SUBMIT_DIR:?Submit this job from the repository root}"
if [[ ! -f "$repo_root/v2/run/Run_One_Subfield.m" ]]; then
    echo "Submit from the T&C-SMEW-N v0.9 repository root." >&2
    exit 2
fi
cd "$repo_root"
if [[ "${TC_UC_CONTROLLER_TOKEN:-}" != "UNCERTAINTY_CAMPAIGN_V2" ]]; then
    echo "Submit through runtime/Uncertainty_1D/submit_uncertainty.sh, not sbatch." >&2
    exit 7
fi
if [[ ! -f v2/SIM/EQUIV_1D_TMP/TC_ERW_results_Kantola_crop_control.mat ]]; then
    echo "The equivalence baseline is missing; copy it from the original's checkout." >&2
    exit 6
fi
if [[ -n "${TC_SSP:-}" && ! -f "v2/T&C_Code/Inputs/Data_Kantola_${TC_SSP}.mat" ]]; then
    echo "The scenario forcing v2/T&C_Code/Inputs/Data_Kantola_${TC_SSP}.mat is missing." >&2
    exit 6
fi

module purge
module load "${MATLAB_MODULE:-MATLAB/2025a}"

node_record="$(scontrol show node "$SLURMD_NODENAME")"
if [[ "$node_record" != *"cpumodel:9575f"* ]]; then
    echo "The job is not running on the required cpumodel:9575f node." >&2
    exit 3
fi
hostname

workers="${SLURM_CPUS_PER_TASK:?The controller reserves one CPU per worker}"
export TC_PARALLEL_WORKERS="$workers"
export TC_NUMERICAL_IDENTITY="cpumodel=9575f;matlab_release=2025a;single_comp_thread=1;parallel_workers=$workers"
export TC_CKPT_MINUTES=30
unset TC_CKPT_FORCE_RESUME || true
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export TC_SSP="${TC_SSP:-}" TC_HORIZON="${TC_HORIZON:-}" TC_CALCITE_CASE="${TC_CALCITE_CASE:-}" \
    TC_UC_WITH_ORIGINAL="${TC_UC_WITH_ORIGINAL:-}"

# Each element keeps its own MATLAB preferences and temporary files, so the
# parallel pools of concurrent elements on one node do not share a job
# storage location.
scratch_parent="${TMPDIR:-/tmp/$USER}"
mkdir -p "$scratch_parent"
scratch_parent="$(cd "$scratch_parent" && pwd -P)"
scratch_root="$scratch_parent/tc_smewn_uc2d_${SLURM_JOB_ID}_${SLURM_RESTART_COUNT:-0}"
mkdir -p "$scratch_root/prefs" "$scratch_root/tmp"
cleanup_scratch() {
    case "$scratch_root" in
        "$scratch_parent"/tc_smewn_uc2d_*) rm -rf -- "$scratch_root" ;;
        *) echo "Refusing to remove unexpected scratch path: $scratch_root" >&2 ;;
    esac
}
trap cleanup_scratch EXIT
export MATLAB_PREFDIR="$scratch_root/prefs"
export TMPDIR="$scratch_root/tmp" TMP="$scratch_root/tmp" TEMP="$scratch_root/tmp"

source runtime/Uncertainty_1D/requeue_near_limit.sh
run_requeueing srun --cpu-bind=cores matlab -singleCompThread -batch \
    "addpath(fullfile(pwd,'runtime','Uncertainty_2D')); Run_Uncertainty_2D_Task(${SLURM_ARRAY_TASK_ID})"
