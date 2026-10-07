#!/usr/bin/env bash
#SBATCH --job-name=RegenEquiv1D
#SBATCH --partition=day
#SBATCH --constraint=cpumodel:9575f
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=regen_equiv_%j.out
#SBATCH --error=regen_equiv_%j.err
# The audit holds the 2-D driver to this baseline bit for bit, so it must be
# built under the identity every campaign runs with: a cpumodel:9575f node,
# R2025a and one computational thread. Built elsewhere it differs by one
# ulp within hours and the all-12 validation fails.

set -euo pipefail
module purge
module load "${MATLAB_MODULE:-MATLAB/2025a}"
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export TC_CKPT_MINUTES=0
export TC_CASE=5
export TC_SCENARIO=EQUIV_1D_TMP
export TC_NDAYS=10
srun --kill-on-bad-exit=1 --cpu-bind=cores matlab -singleCompThread -batch \
    "assert(strcmp(version('-release'),'2025a'),'The equivalence baseline requires MATLAB R2025a.'); cd(fullfile(pwd,'v2','run')); ERW_run"
