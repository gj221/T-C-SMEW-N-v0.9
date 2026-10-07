#!/usr/bin/env bash
#SBATCH --job-name=Uncertainty1D
#SBATCH --requeue
#SBATCH --partition=day
#SBATCH --signal=B:USR1@1800
#SBATCH --open-mode=append
#SBATCH --constraint=cpumodel:9575f
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=12:00:00
#SBATCH --array=1-6
#SBATCH --output=runtime/Uncertainty_1D/slurm_%A_%a.out
#SBATCH --error=runtime/Uncertainty_1D/slurm_%A_%a.err
#
# One element of the Kantola 1-D uncertainty ensemble: one CPU, single
# threaded numerics (-singleCompThread, as every Bouchet campaign). The
# directives above are the 5-year ensemble's (two sets x runs 10-12);
# submit through submit_uncertainty.sh, which compares the model with the
# original's checkout, sets time, memory and the array size for the
# horizon and calcite case, and exports TC_SSP, TC_HORIZON and
# TC_CALCITE_CASE. A run longer than its time limit requeues itself and
# continues from its checkpoint (requeue_near_limit.sh). An element whose
# result already carries its set's provenance exits at once, so a
# requeued or resubmitted array only runs what is missing.

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
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export TC_NUMERICAL_IDENTITY='cpumodel=9575f;matlab_release=2025a;single_comp_thread=1;parallel_workers=0'
export TC_CKPT_MINUTES=30
unset TC_CKPT_FORCE_RESUME || true
export TC_SSP="${TC_SSP:-}" TC_HORIZON="${TC_HORIZON:-}" TC_CALCITE_CASE="${TC_CALCITE_CASE:-}" \
    TC_UC_WITH_ORIGINAL="${TC_UC_WITH_ORIGINAL:-}"

# Each element keeps its own MATLAB preferences and temporary files: MATLAB
# can deadlock on a preferences directory it shares, and on NFS more so.
scratch_parent="${TMPDIR:-/tmp/$USER}"
mkdir -p "$scratch_parent"
scratch_parent="$(cd "$scratch_parent" && pwd -P)"
scratch_root="$scratch_parent/tc_smewn_uc1d_${SLURM_JOB_ID}_${SLURM_RESTART_COUNT:-0}"
mkdir -p "$scratch_root/prefs" "$scratch_root/tmp"
cleanup_scratch() {
    case "$scratch_root" in
        "$scratch_parent"/tc_smewn_uc1d_*) rm -rf -- "$scratch_root" ;;
        *) echo "Refusing to remove unexpected scratch path: $scratch_root" >&2 ;;
    esac
}
trap cleanup_scratch EXIT
export MATLAB_PREFDIR="$scratch_root/prefs"
export TMPDIR="$scratch_root/tmp" TMP="$scratch_root/tmp" TEMP="$scratch_root/tmp"

source runtime/Uncertainty_1D/requeue_near_limit.sh
run_requeueing matlab -singleCompThread -batch "addpath('runtime/Uncertainty_1D'); \
    Run_Uncertainty_Task(${SLURM_ARRAY_TASK_ID})"
