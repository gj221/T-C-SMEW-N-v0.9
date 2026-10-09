#!/usr/bin/env bash
#SBATCH --job-name=UCsummary
#SBATCH --partition=day
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=summarize_uncertainty_%j.out
#
# The uncertainty summary of one arm (Summarize_Uncertainty.m beside this
# file): the original, lower and upper sets of the runs UC_RUNS, as field
# means on the daily axis with yearly statistics, envelopes and the
# observations. Post-processing only: it reads results and writes the
# summary file and two CSV tables beside UC_OUTPUT; nothing in the root
# changes, so it may run while campaigns run.
#
#   sbatch --account=pi_np363 \
#     --export=ALL,UC_SOURCES="<original dir>|<ensemble dir>[;<original dir>|<ensemble dir>]",UC_RUNS="[1 2 3]",UC_OUTPUT=<file.mat> \
#     summarize_uncertainty.sbatch
#
# Each source pair is one campaign: the directory of the original set's
# results and the ensemble directory holding lower/ and upper/. A 35-year
# field-grid run reads about 2-5 GB of stores and day summaries per
# result; 64G covers one result at a time with room for the summary.

set -euo pipefail
umask 027
root="${SLURM_SUBMIT_DIR:?Submit from the package root}"
root="${root%/}"
if [[ ! -f "$root/Summarize_Uncertainty.m" ]]; then
    echo "Submit from the package root that holds Summarize_Uncertainty.m." >&2
    exit 2
fi
: "${UC_SOURCES:?Set UC_SOURCES to <original dir>|<ensemble dir> pairs separated by ;}"
: "${UC_RUNS:?Set UC_RUNS to the runs, such as [1 2 3] or 1:9}"
: "${UC_OUTPUT:?Set UC_OUTPUT to the summary file to write}"
runs_pattern='^[][0-9:, ]+$'
if [[ ! "$UC_RUNS" =~ $runs_pattern ]]; then
    echo "UC_RUNS must be a MATLAB run list such as [1 2 3] or 1:9." >&2
    exit 2
fi
if [[ "$UC_OUTPUT" != *.mat ]]; then
    echo "UC_OUTPUT must name a .mat file." >&2
    exit 2
fi

module purge
module load "${MATLAB_MODULE:-MATLAB/2025a}"
scratch_parent="${TMPDIR:-/tmp/$USER}"
mkdir -p "$scratch_parent"
scratch_dir="$(cd "$scratch_parent" && pwd -P)/tc_uc_summary_${SLURM_JOB_ID}"
mkdir -p "$scratch_dir/prefs" "$scratch_dir/tmp"
trap 'rm -rf -- "$scratch_dir"' EXIT
export MATLAB_PREFDIR="$scratch_dir/prefs" TMPDIR="$scratch_dir/tmp"

matlab -singleCompThread -batch "addpath('$root'); pairs = strsplit(getenv('UC_SOURCES'),';'); sources = cell(numel(pairs),2); for k = 1:numel(pairs), parts = strsplit(pairs{k},'|'); sources(k,:) = parts; end; Summarize_Uncertainty(sources,'Runs',str2num(getenv('UC_RUNS')),'Output',getenv('UC_OUTPUT'));"
