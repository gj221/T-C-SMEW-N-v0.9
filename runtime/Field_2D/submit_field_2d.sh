#!/usr/bin/env bash
# Submit the five-year Kantola runs 1 to 9 on a field grid (T&C-SMEW-N
# v1.0's polygon grid, its outer band the margin, with or without
# topography) on Bouchet: one validation element (the complete audit, the
# field-grid tests and a one-day run of the grid) followed by the
# production array. The grid is a file saved by Build_Field_Grid_Topology
# that lives beside the deployed root, never inside it.
#
# Usage, from the checkout root:
#   TC_FIELD_TOPOLOGY=/full/path/beside/the/root/<name>.mat \
#   TC_FIELD_MEM=320G TC_FIELD_TIME=4-00:00:00 \
#   MATLAB_MODULE=MATLAB/2025a TC_BOUCHET_ACCOUNT=<pi account> \
#   bash runtime/Field_2D/submit_field_2d.sh
# TC_FIELD_MEM and TC_FIELD_TIME have no default: they scale with the
# grid's cell count (see the README), so they are chosen per grid.
# Stages (TC_FIELD_STAGE): fresh (default; result and log trees must be
# empty) or resume (after TIMEOUT or node loss; refuses while any job of
# this campaign is queued and while a run lock remains).

set -Eeuo pipefail

repo_root="$(pwd -P)"
campaign_dir="runtime/Field_2D"
campaign_batch="$campaign_dir/submit_field_2d.sbatch"
pin_list="runtime/Base_1D/shared_pins.sha256"
equivalence_baseline="$repo_root/v2/SIM/EQUIV_1D_TMP/TC_ERW_results_Kantola_crop_control.mat"
minimum_free_gib=250
stage="${TC_FIELD_STAGE:-fresh}"

if [[ ! -f "$repo_root/v2/run/Run_One_Subfield.m" ]]; then
    echo "Run this controller from the T&C-SMEW-N v0.9 checkout root." >&2
    exit 2
fi
case "$stage" in
    fresh|resume) ;;
    *)
        echo "TC_FIELD_STAGE must be fresh or resume, not '$stage'. No jobs were submitted." >&2
        exit 2 ;;
esac
: "${TC_BOUCHET_ACCOUNT:?Set TC_BOUCHET_ACCOUNT to the pi account reported by slurm_checkup.sh}"
: "${TC_FIELD_TOPOLOGY:?Set TC_FIELD_TOPOLOGY to the full path of the grid file}"
: "${TC_FIELD_MEM:?Set TC_FIELD_MEM (for example 320G) from the grid's cell count; see the README}"
: "${TC_FIELD_TIME:?Set TC_FIELD_TIME (for example 4-00:00:00) from the grid's cell count; see the README}"
if [[ ! "$TC_FIELD_MEM" =~ ^[1-9][0-9]*G$ ]]; then
    echo "TC_FIELD_MEM must be a whole number of gigabytes such as 320G. No jobs were submitted." >&2
    exit 2
fi
if [[ ! "$TC_FIELD_TIME" =~ ^[0-9]+-[0-9]{2}:[0-9]{2}:[0-9]{2}$ ]]; then
    echo "TC_FIELD_TIME must read D-HH:MM:SS. No jobs were submitted." >&2
    exit 2
fi

# The grid: an existing <name>.mat outside the deployed root, its name of
# letters, digits and underscores (it names the campaign).
if [[ "$TC_FIELD_TOPOLOGY" != /* || ! -r "$TC_FIELD_TOPOLOGY" ]]; then
    echo "TC_FIELD_TOPOLOGY must be the full path of a readable grid file. No jobs were submitted." >&2
    exit 2
fi
topology_file="$(cd "$(dirname "$TC_FIELD_TOPOLOGY")" && pwd -P)/$(basename "$TC_FIELD_TOPOLOGY")"
case "$topology_file" in
    "$repo_root"/*)
        echo "Keep the grid file beside the deployed root, not inside it: $topology_file" >&2
        echo 'No jobs were submitted.' >&2
        exit 2 ;;
esac
topology_name="$(basename "$topology_file" .mat)"
if [[ "$(basename "$topology_file")" != "$topology_name.mat" || \
      ! "$topology_name" =~ ^[A-Za-z0-9_]{1,32}$ ]]; then
    echo "The grid file must be <name>.mat with a name of at most 32 letters, digits and underscores. No jobs were submitted." >&2
    exit 2
fi
export TC_FIELD_TOPOLOGY="$topology_file"
tag="Field_2D__${topology_name}"
results_dir="$repo_root/v2/SIM/Subfield_12/$tag"
logs_dir="$repo_root/v2/logs/Subfield_12/$tag"

required_files=(
    "$campaign_batch"
    "$campaign_dir/Run_Field_2D_Task.m"
    "$campaign_dir/Check_Lateral_Conservation.m"
    "$campaign_dir/README.md"
    "$pin_list"
    'runtime/Base_1D/Campaign_Parameters.m'
    'runtime/Base_1D/Assert_Base_Provenance.m'
    'runtime/Base_1D/Fingerprint_Digest.m'
    'v2/T&C_Code/MAIN_FRAME_BG_ERW_2D.m'
    'v2/T&C_Code/Bouchet_Campaign_Validation_Fingerprint.m'
    'v2/T&C_Code/Bouchet_Check_Existing_Result.m'
    'v2/T&C_Code/Inputs/Data_Kantola.mat'
    'v2/T&C_Code/Inputs/Data_Kantola_BG.mat'
    'v2/T&C_Code/Inputs/ERW_MOD_PARAM_parameters.xlsx'
    'v2/T&C_Code/Tests/run_post_legacy_equation_audit_tests.m'
    'v2/run/Run_Matrix_12.m'
    'v2/run/Run_Subfield_Worker.m'
    'v2/run/Run_One_Subfield.m'
)
for required_file in "${required_files[@]}"; do
    if [[ ! -r "$repo_root/$required_file" ]]; then
        echo "Required campaign file is missing or unreadable: $repo_root/$required_file" >&2
        echo 'No jobs were submitted.' >&2
        exit 2
    fi
done

if [[ ! -f "$equivalence_baseline" ]]; then
    echo "The equivalence baseline is missing: $equivalence_baseline" >&2
    echo 'Regenerate it for this deployment first (regen_equiv_1d.sbatch). No jobs were submitted.' >&2
    exit 6
fi
if ! (cd "$repo_root" && LC_ALL=C sha256sum -c --quiet --strict "$pin_list"); then
    echo "The checkout differs from the state pinned in $pin_list. No jobs were submitted." >&2
    exit 6
fi
echo "Checkout pins: $(wc -l < "$repo_root/$pin_list") files match $pin_list"
echo "Grid: $topology_file ($(sha256sum "$topology_file" | cut -c1-12))"

for state_path in "$results_dir" "$logs_dir"; do
    if [[ -L "$state_path" || ( -e "$state_path" && ! -d "$state_path" ) ]]; then
        echo "Campaign path is not a directory: $state_path" >&2
        echo 'No jobs were submitted.' >&2
        exit 4
    fi
    if [[ "$stage" == fresh && -d "$state_path" ]]; then
        if [[ -n "$(find "$state_path" -mindepth 1 -print -quit)" ]]; then
            echo "Existing campaign state found under $state_path." >&2
            echo 'A fresh submission needs empty campaign trees; use TC_FIELD_STAGE=resume to continue.' >&2
            echo 'No jobs were submitted.' >&2
            exit 4
        fi
    fi
done

unset TC_REBUILD_EQUIVALENCE TC_VALIDATE_ONLY TC_RUN_ID TC_CKPT_FORCE_RESUME
for submit_env in "${!SBATCH_@}"; do
    unset "$submit_env"
done
unset SLURM_CLUSTERS

export MATLAB_MODULE="${MATLAB_MODULE:-MATLAB/2025a}"
if ! module spider "$MATLAB_MODULE" >/dev/null 2>&1; then
    echo "Bouchet does not currently expose $MATLAB_MODULE. No jobs were submitted." >&2
    exit 3
fi
if ! week_features="$(sinfo -h -p week -o '%f' 2>/dev/null)" || \
        [[ "$week_features" != *cpumodel:9575f* ]]; then
    echo 'The week partition does not expose cpumodel:9575f nodes. No jobs were submitted.' >&2
    exit 3
fi

minimum_free_kib=$((minimum_free_gib * 1024 * 1024))
if ! df_record="$(LC_ALL=C df -Pk "$repo_root")"; then
    echo "Cannot read free space for $repo_root. No jobs were submitted." >&2
    exit 4
fi
available_kib="$(printf '%s\n' "$df_record" | awk 'NR == 2 {print $4}')"
if [[ ! "$available_kib" =~ ^[0-9]+$ ]] || (( available_kib < minimum_free_kib )); then
    echo "Fewer than $minimum_free_gib GiB are free on the repository filesystem. No jobs were submitted." >&2
    exit 4
fi
printf 'Storage preflight: %.1f GiB free on the repository filesystem.\n' \
    "$(awk -v kib="$available_kib" 'BEGIN {print kib/1024/1024}')"

mkdir -p "$logs_dir"
controller_lock="$logs_dir/.submission_controller.lock"
if ! mkdir "$controller_lock" 2>/dev/null; then
    echo "Submission lock already exists: $controller_lock" >&2
    echo 'No jobs were submitted. Check for an active controller or stale lock.' >&2
    exit 5
fi
release_controller_lock() {
    [[ -d "$controller_lock" ]] || return 0
    if ! rmdir "$controller_lock" 2>/dev/null; then
        echo "Warning: could not remove submission lock $controller_lock" >&2
    fi
}
trap release_controller_lock EXIT

job_names="Field2Dcheck_${topology_name},Field2D_${topology_name}"
if ! active_jobs="$(squeue -h -u "$USER" --name="$job_names" -o '%A|%j|%T')"; then
    echo 'Cannot check for an active campaign graph. No jobs were submitted.' >&2
    exit 5
fi
if [[ -n "$active_jobs" ]]; then
    echo 'An active graph for this grid already exists. No jobs were submitted.' >&2
    printf '%s\n' "$active_jobs" >&2
    exit 5
fi
stale_locks="$(find "$logs_dir" -type d -name 'run_*.lock' 2>/dev/null || true)"
if [[ -n "$stale_locks" ]]; then
    echo 'Stale run locks remain from killed jobs. No graph job is active, so remove them and rerun:' >&2
    while IFS= read -r lock_dir; do
        [[ -n "$lock_dir" ]] || continue
        printf '  rm -r -- %q\n' "$lock_dir" >&2
    done <<< "$stale_locks"
    echo 'No jobs were submitted.' >&2
    exit 5
fi

submitted_jobs=()
graph_complete=0
rollback() {
    local status="$1"
    trap - ERR HUP INT TERM
    if (( graph_complete == 0 && ${#submitted_jobs[@]} > 0 )); then
        echo 'Submission graph is incomplete. Cancelling submitted jobs.' >&2
        if ! scancel "${submitted_jobs[@]}"; then
            echo "Automatic cancellation failed. Manual cancellation required: scancel ${submitted_jobs[*]}" >&2
        fi
    fi
    release_controller_lock
    exit "$status"
}
trap 'rollback $?' ERR
trap 'rollback 129' HUP
trap 'rollback 130' INT
trap 'rollback 143' TERM

submit_job() {
    local result_var="$1"
    shift
    local submitted submitted_id
    submitted="$(sbatch --parsable --hold --account="$TC_BOUCHET_ACCOUNT" "$@")"
    submitted_id="${submitted%%;*}"
    printf -v "$result_var" '%s' "$submitted_id"
    submitted_jobs+=("$submitted_id")
}

export_common="ALL,TC_FIELD_CONTROLLER_TOKEN=FIELD2D_V1,TC_FIELD_TOPOLOGY=$topology_file"
submit_job validation_job \
    --partition=day --array=1 --ntasks=1 --cpus-per-task=12 --mem=96G --time=04:00:00 \
    --job-name="Field2Dcheck_${topology_name}" \
    --export="$export_common,TC_VALIDATE_ONLY=1" "$campaign_batch"
submit_job production_job \
    --dependency="afterok:${validation_job}" --kill-on-invalid-dep=yes \
    --array='1-9' --partition=week --mem="$TC_FIELD_MEM" --time="$TC_FIELD_TIME" \
    --job-name="Field2D_${topology_name}" \
    --export="$export_common" "$campaign_batch"

scontrol release "${submitted_jobs[@]}"
graph_complete=1
trap - ERR HUP INT TERM

printf 'Campaign: %s\n' "$tag"
printf 'Validation: %s (day, 4 h)\n' "$validation_job"
printf 'Production: runs 1-9 %s (week, %s, %s each), eligible after the validation passes\n' \
    "$production_job" "$TC_FIELD_MEM" "$TC_FIELD_TIME"
printf 'Monitor: sacct -j %s,%s -X --format=JobID%%16,JobName%%30,State%%12,Elapsed,MaxRSS,NodeList%%16\n' \
    "$validation_job" "$production_job"
