#!/usr/bin/env bash
#SBATCH --job-name=CmpERWk
#SBATCH --partition=day
#SBATCH --constraint=cpumodel:9575f
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=06:00:00
#SBATCH --output=compare_erw_kernels_%j.out
#
# The same model runs under TC_ERW_KERNEL=matlab (the original code), fast
# and mex, every saved variable compared with compare_u_v.m:
#   - run 11 (basalt, one column) from 5 November 2017 for 375 days, across
#     the applications of 10 November 2017 and 2018, serial: matlab, fast
#     and, when this root holds the MEX files (build_erw_mex.sbatch), mex;
#   - run 2 (basalt on the whole 12-column plot) from 5 November 2017 for
#     40 days with 6 thread-pool workers: matlab and fast (MATLAB runs no MEX
#     function on a thread-based worker).
# The five runs share the node, so their wall times compare like for like;
# the PASS lines and the times are at the end of this job's .out.
# The runs write to scratch (CMP_DIR); nothing is written inside the root
# except this job's output file beside it.
# Submit from the package root:
#   sbatch --account=pi_np363 compare_erw_kernels.sbatch

set -euo pipefail
root="${SLURM_SUBMIT_DIR:?Submit from the package root}"
if [[ ! -f "$root/v2/T&C_Code/ERW_Chemistry_Kernel_Mode.m" ]]; then
    echo "No ERW kernel switch under $root; submit from the package root." >&2
    exit 2
fi
out="${CMP_DIR:-/nfs/roberts/scratch/pi_np363/$USER/compare_erw_kernels_${SLURM_JOB_ID}}"
mkdir -p "$out"

module purge
module load "${MATLAB_MODULE:-MATLAB/2025a}"
export TC_CKPT_MINUTES=0 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
unset TC_ERW_KERNEL TC_PARALLEL_WORKERS

have_mex=0
if compgen -G "$root/v2/T&C_Code/ERW_Appelo_Seed_Kernel_mex.mexa64" > /dev/null && \
   compgen -G "$root/v2/T&C_Code/ERW_Charge_Root_Kernel_mex.mexa64" > /dev/null; then
    have_mex=1
fi

run_one() {  # label kernel run-call
    mkdir -p "$out/$1" "$out/prefs_$1" "$out/tmp_$1"
    local start=$SECONDS
    TC_ERW_KERNEL="$2" MATLAB_PREFDIR="$out/prefs_$1" TMPDIR="$out/tmp_$1" \
        matlab -singleCompThread -batch \
        "cd(fullfile('$root','v2','run')); $3" > "$out/run_$1.log" 2>&1
    echo "$1 ($2) finished in $(( SECONDS - start )) s" >> "$out/times.txt"
}
call11="Run_One_Subfield(11,'StartDate','2017-11-05','NDays',375,'ResultDirectory'"
call2="Run_One_Subfield(2,'StartDate','2017-11-05','NDays',40,'ParallelWorkers',6,'ResultDirectory'"
pids=()
run_one r11_matlab matlab "${call11},'$out/r11_matlab')" & pids+=($!)
run_one r11_fast   fast   "${call11},'$out/r11_fast')"   & pids+=($!)
if (( have_mex )); then
    run_one r11_mex mex "${call11},'$out/r11_mex')" & pids+=($!)
else
    echo "No MEX files in this root: the mex run is skipped (build_erw_mex.sbatch builds them)."
fi
run_one r2_matlab matlab "${call2},'$out/r2_matlab')" & pids+=($!)
run_one r2_fast   fast   "${call2},'$out/r2_fast')"   & pids+=($!)
status=0
for p in "${pids[@]}"; do wait "$p" || status=1; done
cat "$out/times.txt"
if (( status != 0 )); then
    echo "A run failed; see $out/run_*.log" >&2
    exit 1
fi

cd "$root"
compare() {  # reference candidate
    echo "=== $2 against $1"
    MATLAB_PREFDIR="$out/prefs_r11_matlab" TMPDIR="$out/tmp_r11_matlab" \
        matlab -singleCompThread -batch "compare_u_v('$out/$1','$out/$2')"
}
compare r11_matlab r11_fast
if (( have_mex )); then compare r11_matlab r11_mex; fi
compare r2_matlab r2_fast
echo "Results kept in $out"
