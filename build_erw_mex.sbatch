#!/usr/bin/env bash
#SBATCH --job-name=ERWmex
#SBATCH --partition=day
#SBATCH --constraint=cpumodel:9575f
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=build_erw_mex_%j.out
#
# The compiled ERW chemistry kernels for this root, and the proof that the
# restructured kernels reproduce the original code bit for bit.
#
#   1. Build_ERW_Chemistry_Mex compiles ERW_Appelo_Seed_Kernel and
#      ERW_Charge_Root_Kernel with MATLAB Coder into
#      v2/T&C_Code/ERW_*_Kernel_mex.mexa64 (+ a manifest of the source
#      hashes) and checks them against the MATLAB kernels;
#   2. test_erw_chemistry_kernels compares TC_ERW_KERNEL=fast and =mex with
#      the original code (3000 reduced Appelo solves over pore pH 3.5-9.5,
#      the 281-point charge scan, all seven ERW_function outputs on the
#      step-22053 fixture variants and three 48-step chains), every value
#      bit for bit;
#   3. erw_kernel_benchmark times one ERW_function step in each mode.
#
# A MEX runs only on the platform and MATLAB release it was built with, so
# every deployed root runs this once, before its validation jobs. The
# campaign fingerprints hash *.m sources, not the MEX files; the MEX files,
# their manifest and this job's .out are the only files written in the
# root. MEX functions do not run on thread-based workers, so a 2-D run whose
# daily ERW columns use parpool('Threads') (TC_PARALLEL_WORKERS or
# TC_FIELD_WORKERS > 1) takes TC_ERW_KERNEL=fast; =mex serves 1-D and
# serial 2-D runs.
#
# Submit from the package root:
#   sbatch --account=pi_np363 build_erw_mex.sbatch
# GCC_MODULE=<module> loads a compiler module when the node's default gcc
# is not one MATLAB accepts.

set -euo pipefail
umask 027

root="${SLURM_SUBMIT_DIR:?Submit from the package root}"
if [[ ! -f "$root/v2/T&C_Code/Build_ERW_Chemistry_Mex.m" ]]; then
    echo "No ERW kernel sources under $root; submit from the package root." >&2
    exit 2
fi
work="${ERW_MEX_DIR:-/nfs/roberts/scratch/pi_np363/$USER/erw_mex_${SLURM_JOB_ID}}"
mkdir -p "$work/prefs" "$work/tmp" "$work/build"

module purge
module load "${MATLAB_MODULE:-MATLAB/2025a}"
if [[ -n "${GCC_MODULE:-}" ]]; then
    module load "$GCC_MODULE"
fi
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
unset TC_ERW_KERNEL

echo "node ${SLURMD_NODENAME:-unknown}  root $root"
gcc --version | head -n 1 || echo "no gcc on PATH"

cd "$root"
MATLAB_PREFDIR="$work/prefs" TMPDIR="$work/tmp" matlab -singleCompThread -batch \
    "code_dir = fullfile(pwd,'v2','T&C_Code'); addpath(code_dir,fullfile(code_dir,'Tests')); Build_ERW_Chemistry_Mex('$work/build'); test_erw_chemistry_kernels(3000,48); erw_kernel_benchmark(25)"
echo "Coder working folder kept in $work"
