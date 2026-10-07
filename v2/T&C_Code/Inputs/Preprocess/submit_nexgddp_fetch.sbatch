#!/usr/bin/env bash
#SBATCH --job-name=NexGddpFetch
#SBATCH --partition=day
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --requeue
#SBATCH --output=nexgddp_fetch_%j.out
#SBATCH --error=nexgddp_fetch_%j.err

# Pull one NEX-GDDP grid cell per file into local NetCDFs.
#
# This runs under Slurm rather than on the login node for two reasons: the
# login node's watchdog kills a long multi-threaded job, and MATLAB's netCDF
# DAP client fails here (NC_EIO) so the retrieval has to happen in Python
# either way. Python's netCDF4 DAP client does work on compute nodes.
#
# Submit from the directory holding fetch_nexgddp_cell.py:
#   sbatch submit_nexgddp_fetch.sbatch
#
# Resumable: completed outputs are skipped, so a requeue or timeout costs
# only the files in flight. Re-submit until it reports failed=0.

set -euo pipefail
umask 027

work_dir="${SLURM_SUBMIT_DIR:?Submit from the directory holding fetch_nexgddp_cell.py}"
cd "$work_dir"
if [[ ! -f fetch_nexgddp_cell.py ]]; then
    echo "fetch_nexgddp_cell.py not found in $work_dir" >&2
    exit 2
fi

module purge
module load "${PYTHON_MODULE:-Python/3.12.3-GCCcore-13.3.0}"
source "${VENV:-$HOME/era5env}/bin/activate"

hostname
python -c "import xarray, netCDF4, requests; print('deps ok')"

# The work is server-side latency, not bandwidth or local CPU, so threads
# scale well past the core count. Kept moderate to stay a polite guest.
export NEXGDDP_WORKERS="${NEXGDDP_WORKERS:-24}"
echo "workers=$NEXGDDP_WORKERS"

# Partial writes from any previous killed run; the fetcher writes .part and
# renames, so a leftover .part is never a valid output.
rm -f nexgddp_cell/*.part 2>/dev/null || true

python -u fetch_nexgddp_cell.py

echo "Files fetched:"
ls nexgddp_cell | wc -l
