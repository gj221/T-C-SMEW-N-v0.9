#!/usr/bin/env bash
#SBATCH --job-name=CmpUV
#SBATCH --partition=day
#SBATCH --constraint=cpumodel:9575f
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=03:00:00
#SBATCH --output=compare_u_v_%j.out
#
# One year of 1-D basalt (run 11) from 5 November 2017, across the
# applications of 10 November 2017 and 2018, run in package u and in
# package v side by side, then every saved variable compared
# (compare_u_v.m). Two applications are what make v keep the particle
# histories in blocks, and a year is enough to drop hundreds of days.
# Submit from the v root with U_ROOT set to the u root:
#   sbatch --account=pi_np363 --export=ALL,U_ROOT=<u root> compare_u_v.sbatch
# The runs write to scratch (CMP_DIR); nothing is written inside either
# root except this job's output file beside it.

set -euo pipefail
v_root="${SLURM_SUBMIT_DIR:?Submit from the v root}"
u_root="${U_ROOT:?Set U_ROOT to the u root}"
out="${CMP_DIR:-/nfs/roberts/scratch/pi_np363/$USER/compare_u_v_${SLURM_JOB_ID}}"
for d in "$v_root" "$u_root"; do
    if [[ ! -f "$d/v2/run/Run_One_Subfield.m" ]]; then
        echo "No T&C-SMEW-N checkout at $d" >&2
        exit 2
    fi
done
mkdir -p "$out/u" "$out/v" "$out/prefs_u" "$out/prefs_v" "$out/tmp_u" "$out/tmp_v"

module purge
module load "${MATLAB_MODULE:-MATLAB/2025a}"
export TC_CKPT_MINUTES=0 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

run_one() {  # root label
    MATLAB_PREFDIR="$out/prefs_$2" TMPDIR="$out/tmp_$2" matlab -singleCompThread -batch \
        "cd(fullfile('$1','v2','run')); Run_One_Subfield(11,'StartDate','2017-11-05','NDays',375,'ResultDirectory','$out/$2')" \
        > "$out/run_$2.log" 2>&1
}
run_one "$u_root" u &
pid_u=$!
run_one "$v_root" v &
pid_v=$!
status=0
wait "$pid_u" || { echo "The u run failed; see $out/run_u.log" >&2; status=1; }
wait "$pid_v" || { echo "The v run failed; see $out/run_v.log" >&2; status=1; }
if (( status != 0 )); then exit 1; fi

cd "$v_root"
MATLAB_PREFDIR="$out/prefs_v" TMPDIR="$out/tmp_v" matlab -singleCompThread -batch \
    "compare_u_v('$out/u','$out/v')"
echo "Results kept in $out"
