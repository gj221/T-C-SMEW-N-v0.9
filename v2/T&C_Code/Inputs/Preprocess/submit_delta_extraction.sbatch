#!/usr/bin/env bash
#SBATCH --job-name=DeltaExtract
#SBATCH --partition=day
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --requeue
#SBATCH --output=delta_extract_%j.out
#SBATCH --error=delta_extract_%j.err

# CMIP6 monthly change factors for the three SSPs.
#
# Reads the local single-cell NetCDFs staged by fetch_nexgddp_cell.py, so
# this is local file I/O, not network: one task and one CPU is right, and it
# runs in minutes. It was network-bound and asked for 12 hours before the
# fetch step existed. A short walltime backfills into gaps in the schedule
# instead of waiting for a large hole, which is worth real queue time here.
#
# Submit from the directory holding the .m files:
#   sbatch submit_delta_extraction.sbatch
#
# Resumable at scenario granularity. Run_Extract_Deltas_Local skips any
# scenario whose forcing/deltas_<scenario>.mat already exists, so a requeue
# or a timeout loses at most the scenario in flight.
#
# That skip is also why the cell-directory check below refuses to start on a
# partial fetch: a deltas file built from a fraction of the ensemble would be
# skipped on every later run rather than rebuilt, and the shortfall would
# carry silently into the forcing.

set -euo pipefail
umask 027

work_dir="${SLURM_SUBMIT_DIR:?Submit from the directory holding the .m files}"
cd "$work_dir"
for f in Run_Extract_Deltas_Local.m Delta_Forcing_Config.m Extract_CMIP6_Monthly_Deltas.m; do
    if [[ ! -f "$f" ]]; then
        echo "Missing $f in $work_dir" >&2
        exit 2
    fi
done

module purge
module load "${MATLAB_MODULE:-MATLAB/2025a}"

hostname
module list

# Fail now rather than after MATLAB has started. The old check here probed
# NCCS for OPeNDAP reachability, which this job no longer needs; the
# precondition that matters now is a complete cell directory.
cell_dir="${NEXGDDP_CELL_DIR:-$work_dir/nexgddp_cell}"
expected="${NEXGDDP_EXPECTED:-5500}"
if [[ ! -d "$cell_dir" ]]; then
    echo "Missing cell directory: $cell_dir" >&2
    echo "Run fetch_nexgddp_cell.py first." >&2
    exit 3
fi
have=$(find "$cell_dir" -maxdepth 1 -type f -name '*.nc' | wc -l)
if (( have < expected )); then
    echo "Cell directory holds $have of $expected files: the fetch is incomplete." >&2
    echo "Finish fetch_nexgddp_cell.py before extracting; it skips what is" >&2
    echo "already on disk, so resubmitting it only retrieves the gaps." >&2
    exit 3
fi
echo "Cell directory complete: $have files in $cell_dir"

# MATLAB deadlocks on a shared preferences directory, and an NFS prefdir is
# worse. Keep both prefs and temp on node-local scratch for this job only.
scratch_parent="${TMPDIR:-/tmp/$USER}"
mkdir -p "$scratch_parent"
scratch_root="$(mktemp -d "${scratch_parent}/tc_delta_${SLURM_JOB_ID}_XXXXXX")"
cleanup_scratch() {
    case "$scratch_root" in
        "$scratch_parent"/tc_delta_*) rm -rf -- "$scratch_root" ;;
        *) echo "Refusing to remove unexpected scratch path: $scratch_root" >&2 ;;
    esac
}
trap cleanup_scratch EXIT

export MATLAB_PREFDIR="$scratch_root/prefs"
export TMPDIR="$scratch_root/tmp"
export TMP="$TMPDIR"
export TEMP="$TMPDIR"
mkdir -p "$MATLAB_PREFDIR" "$TMPDIR"

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

srun --cpu-bind=cores matlab -singleCompThread -nodisplay -nosplash \
    -logfile "$work_dir/delta_extract_${SLURM_JOB_ID}.log" \
    -batch "addpath(pwd); Run_Extract_Deltas_Local"

echo "Delta extraction finished."
ls -la "$work_dir/forcing/"
