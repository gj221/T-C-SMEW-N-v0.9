#!/usr/bin/env bash
#SBATCH --job-name=CmpERWk
#SBATCH --partition=day
#SBATCH --constraint=cpumodel:9575f
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=128G
#SBATCH --time=06:00:00
#SBATCH --output=compare_erw_kernels_%j.out
#
# The same model runs under TC_ERW_KERNEL=matlab (the original code), fast
# and mex, with the calibrated Base_1D set (Campaign_Parameters), every
# saved value compared (compare_erw_kernels.m):
#   - run 11 (basalt, one column, full records) from 5 November 2017 for
#     375 days, across the applications of 10 November 2017 and 2018,
#     serial: matlab, fast and, when this root holds the MEX files
#     (build_erw_mex.sbatch), mex;
#   - run 2 (basalt on the whole 12-column plot) from 5 November 2017 for
#     40 days with 6 thread-pool workers, kept as a multi-square run keeps
#     it (daily stores): matlab and fast (MATLAB runs no MEX function on a
#     thread-based worker).
# With X_ROOT set to the package-x root, x runs both as well (x is the
# original code) and each fast run is checked against x's with
# compare_x_y.m, every output against the same reduction of x's records.
# The runs share the node, so their wall times compare like for like; the
# PASS lines, the times and ALL PASS are at the end of this job's .out.
# The runs write to scratch (CMP_DIR); nothing is written inside the root
# except this job's output file beside it.
# Submit from the package root:
#   sbatch --account=pi_np363 compare_erw_kernels.sbatch
#   sbatch --account=pi_np363 --export=ALL,X_ROOT=<x root> compare_erw_kernels.sbatch

set -euo pipefail
umask 027
root="${SLURM_SUBMIT_DIR:?Submit from the package root}"
root="${root%/}"
if [[ ! -f "$root/v2/T&C_Code/ERW_Chemistry_Kernel_Mode.m" || \
      ! -f "$root/compare_erw_kernels.m" ]]; then
    echo "No ERW kernel switch under $root; submit from the package root." >&2
    exit 2
fi
x_root="${X_ROOT:-}"
x_root="${x_root%/}"
if [[ -n "$x_root" ]]; then
    if [[ ! -f "$x_root/v2/run/Run_One_Subfield.m" || \
          ! -f "$x_root/runtime/Base_1D/Campaign_Parameters.m" ]]; then
        echo "No T&C-SMEW-N v0.9 tree at X_ROOT=$x_root" >&2
        exit 2
    fi
    if grep -q "FULL_RECORDS" "$x_root/v2/T&C_Code/MAIN_FRAME_BG_ERW_2D.m"; then
        echo "X_ROOT is a package-y tree, not x: $x_root" >&2
        exit 2
    fi
fi
out="${CMP_DIR:-/nfs/roberts/scratch/pi_np363/${USER:-$(id -un)}/compare_erw_kernels_${SLURM_JOB_ID}}"
out="${out%/}"
case "$out/" in
    "$root"/*) echo "Keep CMP_DIR outside the root: $out" >&2; exit 2 ;;
esac
if [[ -n "$x_root" ]]; then
    case "$out/" in
        "$x_root"/*) echo "Keep CMP_DIR outside X_ROOT: $out" >&2; exit 2 ;;
    esac
fi
mkdir -p "$out"

module purge
module load "${MATLAB_MODULE:-MATLAB/2025a}"
export TC_CKPT_MINUTES=0 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
unset TC_ERW_KERNEL TC_CKPT_FORCE_RESUME TC_FULL_RECORDS TC_PARALLEL_WORKERS \
    TC_VALIDATION_FINGERPRINT TC_SWEEP_LAYOUT TC_SWEEP_CASE_ID ERW_DEBUG_CANDIDATES || true

have_mex=0
if compgen -G "$root/v2/T&C_Code/ERW_Appelo_Seed_Kernel_mex.mexa64" > /dev/null && \
   compgen -G "$root/v2/T&C_Code/ERW_Charge_Root_Kernel_mex.mexa64" > /dev/null; then
    have_mex=1
fi

# MATLAB deadlocks on a shared preferences directory, and an NFS one is
# worse: each MATLAB gets its own prefs and temp on node-local scratch.
scratch_parent="${TMPDIR:-/tmp/${USER:-$(id -un)}}"
mkdir -p "$scratch_parent"
scratch_parent="$(cd "$scratch_parent" && pwd -P)"
scratch_dir="$scratch_parent/tc_compare_erw_kernels_${SLURM_JOB_ID}"
cleanup() {
    case "$scratch_dir" in
        "$scratch_parent"/tc_compare_erw_kernels_*) rm -rf -- "$scratch_dir" ;;
    esac
}
trap cleanup EXIT

echo "root: $root"
if [[ -n "$x_root" ]]; then echo "x:    $x_root"; fi
echo "results in $out"
hostname

run_one() {  # name root kernel run_id days workers
    local name=$1 r=$2 kernel=$3 run_id=$4 days=$5 nw=$6
    local start=$SECONDS status=0
    mkdir -p "$out/$name" "$scratch_dir/$name/prefs" "$scratch_dir/$name/tmp"
    (
        export TC_ERW_KERNEL="$kernel"    # package x has no switch and ignores it
        export MATLAB_PREFDIR="$scratch_dir/$name/prefs" TMPDIR="$scratch_dir/$name/tmp"
        matlab -singleCompThread -batch \
            "addpath(fullfile('$r','runtime','Base_1D')); P = Campaign_Parameters; cd(fullfile('$r','v2','run')); Run_One_Subfield($run_id,'StartDate','2017-11-05','NDays',$days,'AppeloAlpha',P.AppeloAlpha,'AppeloLogKH',P.AppeloLogKH,'KNmaxd',P.KNmaxd,'Kbot',P.Kbot,'ParallelWorkers',$nw,'ResultDirectory','$out/$name')"
    ) > "$out/run_$name.log" 2>&1 || status=$?
    echo "$name ($kernel) finished in $(( SECONDS - start )) s, exit $status" >> "$out/times.txt"
    return $status
}

compare_kernels() {  # reference candidate reference_kernel candidate_kernel
    echo
    echo "=== $2 against $1"
    (
        export MATLAB_PREFDIR="$scratch_dir/compare/prefs" TMPDIR="$scratch_dir/compare/tmp"
        mkdir -p "$MATLAB_PREFDIR" "$TMPDIR"
        cd "$root"
        matlab -singleCompThread -batch \
            "addpath('$root'); if ~compare_erw_kernels('$out/$1','$out/$2','$3','$4'), exit(1); end"
    ) 2>&1 | tee "$out/compare_$2.log"
    return "${PIPESTATUS[0]}"
}

compare_x() {  # x_name candidate
    echo
    echo "=== $2 against package x ($1)"
    (
        export MATLAB_PREFDIR="$scratch_dir/compare/prefs" TMPDIR="$scratch_dir/compare/tmp"
        mkdir -p "$MATLAB_PREFDIR" "$TMPDIR"
        cd "$root"
        matlab -singleCompThread -batch \
            "addpath('$root',fullfile('$root','v2','T&C_Code')); if ~compare_x_y('$out/$1','$out/$2'), exit(1); end"
    ) 2>&1 | tee "$out/compare_$2_against_x.log"
    return "${PIPESTATUS[0]}"
}

: > "$out/times.txt"
pids=()
run_one r11_matlab "$root" matlab 11 375 0 & pids+=("$!:r11_matlab")
run_one r11_fast "$root" fast 11 375 0 & pids+=("$!:r11_fast")
if (( have_mex )); then
    run_one r11_mex "$root" mex 11 375 0 & pids+=("$!:r11_mex")
else
    echo "No MEX files in this root: the mex run is skipped (build_erw_mex.sbatch builds them)."
fi
run_one r2_matlab "$root" matlab 2 40 6 & pids+=("$!:r2_matlab")
run_one r2_fast "$root" fast 2 40 6 & pids+=("$!:r2_fast")
if [[ -n "$x_root" ]]; then
    run_one r11_x "$x_root" matlab 11 375 0 & pids+=("$!:r11_x")
    run_one r2_x "$x_root" matlab 2 40 6 & pids+=("$!:r2_x")
fi
status=0
for item in "${pids[@]}"; do
    if ! wait "${item%%:*}"; then
        echo "The ${item#*:} run failed; see $out/run_${item#*:}.log" >&2
        status=1
    fi
done
cat "$out/times.txt"
for f in "$out"/run_*.log; do
    echo "--- $(basename "$f"): $(grep -h 'COMPUTATIONAL TIME' "$f" | tail -1)"
done

compare_kernels r11_matlab r11_fast matlab fast || status=1
if (( have_mex )); then
    compare_kernels r11_matlab r11_mex matlab mex || status=1
fi
compare_kernels r2_matlab r2_fast matlab fast || status=1
if [[ -n "$x_root" ]]; then
    compare_x r11_x r11_fast || status=1
    compare_x r2_x r2_fast || status=1
fi
echo
what="the fast kernels"
if (( have_mex )); then what="$what and the MEX files"; fi
against=""
if [[ -n "$x_root" ]]; then against=", and the fast runs reproduce package x"; fi
if (( status == 0 )); then
    echo "ALL PASS: $what reproduce the original code bit for bit$against."
else
    echo "NOT ALL PASS: see the CHECK lines above and the run logs in $out" >&2
fi
echo "Results kept in $out"
exit $status
