#!/usr/bin/env bash
#SBATCH --job-name=CmpXY
#SBATCH --partition=day
#SBATCH --constraint=cpumodel:9575f
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=160G
#SBATCH --time=23:00:00
#SBATCH --output=compare_x_y_%j.out
#
# Package y against package x: the same runs in both, then every output y
# keeps compared with the same reduction of x's records (compare_x_y.m).
# Run 2 (basalt, whole field) with the calibrated Base_1D set, from
# CMP_START (2016-11-01):
#
#   plot        the 12-column plot, CMP_NDAYS days (40), across the rock
#               application of 10 November 2016; x and y (daily stores)
#   plot_full   the same, y with TC_FULL_RECORDS=1: whole records,
#               compared with x's array for array
#   year        the 12-column plot, CMP_YEAR_DAYS days (400), across the
#               applications of 2016 and 2017, so the particle histories
#               of two applications run in daily blocks; x and y
#   grid        the 10 m Kantola basalt DEM grid (340 squares), CMP_NDAYS
#               days; x and y
#
# The five plot runs go side by side first (12 workers each), then the two
# grid runs (CMP_WORKERS each, 32). Submit from the y root with X_ROOT set
# to the x root:
#   sbatch --account=pi_np363 --export=ALL,X_ROOT=<x root> compare_x_y.sbatch
# GRID_FILE defaults to kantola_basalt_dem_10m.mat beside the y root. The
# runs write to scratch (CMP_DIR); nothing is written inside either root
# except this job's output file beside it.
#
# Resources. A day of the plot takes x about 12 s (6 h for five years), so
# the year runs take about 1.5 h. x holds about 1.2 MB of records per
# square and day, so its 40-day grid run keeps ~16 GiB and the comparison
# loads it beside y's; 160G covers the runs and the comparison. A day of
# the grid takes x 2-6 minutes with 32 workers (field_profile), so the job
# ends in about 4-8 h; 23 h is the margin for a slow node.

set -euo pipefail
umask 027
y_root="${SLURM_SUBMIT_DIR:?Submit from the y root}"
y_root="${y_root%/}"
x_root="${X_ROOT:?Set X_ROOT to the x root}"
x_root="${x_root%/}"
grid="${GRID_FILE:-$(dirname "$y_root")/kantola_basalt_dem_10m.mat}"
ndays="${CMP_NDAYS:-40}"
year_days="${CMP_YEAR_DAYS:-400}"
start_date="${CMP_START:-2016-11-01}"
workers="${CMP_WORKERS:-32}"
out="${CMP_DIR:-/nfs/roberts/scratch/pi_np363/${USER:-$(id -un)}/compare_x_y_${SLURM_JOB_ID}}"
out="${out%/}"

for d in "$y_root" "$x_root"; do
    if [[ ! -f "$d/v2/run/Run_One_Subfield.m" || \
          ! -f "$d/runtime/Base_1D/Campaign_Parameters.m" ]]; then
        echo "No T&C-SMEW-N v0.9 tree at $d" >&2
        exit 2
    fi
done
if [[ ! -f "$y_root/compare_x_y.m" ]] || \
   ! grep -q "FULL_RECORDS" "$y_root/v2/T&C_Code/MAIN_FRAME_BG_ERW_2D.m"; then
    echo "Submit from the y root: $y_root is not package y." >&2
    exit 2
fi
if grep -q "FULL_RECORDS" "$x_root/v2/T&C_Code/MAIN_FRAME_BG_ERW_2D.m"; then
    echo "X_ROOT is a package-y tree, not x: $x_root" >&2
    exit 2
fi
if [[ ! -r "$grid" ]]; then
    echo "The grid file is missing or unreadable: $grid (set GRID_FILE)" >&2
    exit 2
fi
for v in "$ndays" "$year_days" "$workers"; do
    if [[ ! "$v" =~ ^[1-9][0-9]*$ ]]; then
        echo "CMP_NDAYS, CMP_YEAR_DAYS and CMP_WORKERS must be whole numbers." >&2
        exit 2
    fi
done
if [[ ! "$start_date" =~ ^[0-9]{4}-[0-9]{2}-[0-9]{2}$ ]]; then
    echo "CMP_START must be yyyy-mm-dd." >&2
    exit 2
fi
case "$out/" in
    "$y_root"/*|"$x_root"/*) echo "Keep CMP_DIR outside both trees: $out" >&2; exit 2 ;;
esac
mkdir -p "$out"

module purge
module load "${MATLAB_MODULE:-MATLAB/2025a}"
export TC_CKPT_MINUTES=0 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
unset TC_CKPT_FORCE_RESUME TC_FULL_RECORDS TC_PARALLEL_WORKERS TC_VALIDATION_FINGERPRINT \
    TC_SWEEP_LAYOUT TC_SWEEP_CASE_ID ERW_DEBUG_CANDIDATES || true

# MATLAB deadlocks on a shared preferences directory, and an NFS one is
# worse: each MATLAB gets its own prefs and temp on node-local scratch.
scratch_parent="${TMPDIR:-/tmp/${USER:-$(id -un)}}"
mkdir -p "$scratch_parent"
scratch_parent="$(cd "$scratch_parent" && pwd -P)"
scratch_dir="$scratch_parent/tc_compare_x_y_${SLURM_JOB_ID}"
cleanup() {
    case "$scratch_dir" in
        "$scratch_parent"/tc_compare_x_y_*) rm -rf -- "$scratch_dir" ;;
    esac
}
trap cleanup EXIT

echo "x:    $x_root"
echo "y:    $y_root"
echo "grid: $grid"
echo "run 2 from $start_date: $ndays days (plot, grid), $year_days days (year); results in $out"
hostname
sha256sum "$x_root/v2/T&C_Code/MAIN_FRAME_BG_ERW_2D.m" "$y_root/v2/T&C_Code/MAIN_FRAME_BG_ERW_2D.m" "$grid"

run_one() {  # name root workers full_records days [grid]
    local name=$1 root=$2 nw=$3 full=$4 days=$5 topo=${6:-}
    local grid_option=""
    if [[ -n "$topo" ]]; then
        grid_option=",'FieldTopologyFile','$topo'"
    fi
    mkdir -p "$out/$name" "$scratch_dir/$name/prefs" "$scratch_dir/$name/tmp"
    (
        if [[ "$full" == 1 ]]; then export TC_FULL_RECORDS=1; fi
        export MATLAB_PREFDIR="$scratch_dir/$name/prefs" TMPDIR="$scratch_dir/$name/tmp"
        matlab -singleCompThread -batch \
            "addpath(fullfile('$root','runtime','Base_1D')); P = Campaign_Parameters; cd(fullfile('$root','v2','run')); Run_One_Subfield(2,'StartDate','$start_date','NDays',$days,'AppeloAlpha',P.AppeloAlpha,'AppeloLogKH',P.AppeloLogKH,'KNmaxd',P.KNmaxd,'Kbot',P.Kbot,'ParallelWorkers',$nw$grid_option,'ResultDirectory','$out/$name')"
    ) > "$out/run_$name.log" 2>&1
}

compare_one() {  # x_dir y_dir label
    echo
    echo "=== $3 ==="
    (
        export MATLAB_PREFDIR="$scratch_dir/compare/prefs" TMPDIR="$scratch_dir/compare/tmp"
        mkdir -p "$MATLAB_PREFDIR" "$TMPDIR"
        cd "$y_root"
        matlab -singleCompThread -batch \
            "addpath(fullfile('$y_root','v2','T&C_Code')); if ~compare_x_y('$out/$1','$out/$2'), exit(1); end"
    ) 2>&1 | tee "$out/compare_$3.log"
    return "${PIPESTATUS[0]}"
}

wait_all() {  # pid:name ...
    local status=0 item
    for item in "$@"; do
        if ! wait "${item%%:*}"; then
            echo "The ${item#*:} run failed; see $out/run_${item#*:}.log" >&2
            status=1
        fi
    done
    return $status
}

status=0
run_one plot_x "$x_root" 12 0 "$ndays" & p1=$!
run_one plot_y "$y_root" 12 0 "$ndays" & p2=$!
run_one plot_yfull "$y_root" 12 1 "$ndays" & p3=$!
run_one year_x "$x_root" 12 0 "$year_days" & p4=$!
run_one year_y "$y_root" 12 0 "$year_days" & p5=$!
wait_all "$p1:plot_x" "$p2:plot_y" "$p3:plot_yfull" "$p4:year_x" "$p5:year_y" || status=1
run_one grid_x "$x_root" "$workers" 0 "$ndays" "$grid" & p6=$!
run_one grid_y "$y_root" "$workers" 0 "$ndays" "$grid" & p7=$!
wait_all "$p6:grid_x" "$p7:grid_y" || status=1

for f in "$out"/run_*.log; do
    echo "--- $(basename "$f"): $(grep -h 'COMPUTATIONAL TIME' "$f" | tail -1)"
done

compare_one plot_x plot_y plot || status=1
compare_one plot_x plot_yfull plot_full || status=1
compare_one year_x year_y year || status=1
compare_one grid_x grid_y grid || status=1
echo
if (( status == 0 )); then
    echo "ALL PASS: y reproduces x on the plot (40 days, and with full records), over a year of two applications, and on the 10 m grid."
else
    echo "NOT ALL PASS: see the CHECK lines above and the run logs in $out" >&2
fi
echo "Results kept in $out"
exit $status
