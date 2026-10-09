#!/usr/bin/env bash
# Submit one 35-year Kantola SSP campaign on a field grid on Bouchet: the
# runs of TC_FIELD_RUNS (default 1-9; 3,6,9 with TC_CALCITE_CASE) on the
# grid TC_FIELD_TOPOLOGY under the scenario TC_SSP, rock every 10 November
# to 2050. One validation element (the complete audit, the Linux checks
# and a one-day run of the grid under this scenario and case) is followed
# by the production array. The grid is a file saved by
# Build_Field_Grid_Topology that lives beside the deployed root, never
# inside it. One call per grid, scenario and case; each has its own tag,
# trees and job names, so calls run side by side.
#
# Usage, from the checkout root:
#   TC_FIELD_TOPOLOGY=/full/path/beside/the/root/<name>.mat \
#   TC_SSP=ssp245 TC_FIELD_RUNS=2 TC_FIELD_WORKERS=32 TC_FIELD_MEM=192G \
#   MATLAB_MODULE=MATLAB/2025a TC_BOUCHET_ACCOUNT=<pi account> \
#   bash runtime/Field_SSP_2D/submit_field_ssp_2d.sh
# TC_FIELD_MEM has no default: it scales with the grid's square count and
# is set from the pilot (README). TC_FIELD_WORKERS (default 32, 2 to 64)
# sets the pool and the CPUs of every element and enters the validation
# stamp. TC_FIELD_RUNS is the Slurm array of runs, such as 2 for the pilot,
# 1,4,7 or 2-9%4. TC_FIELD_PARTITION is day (default; 24 h pieces, each run
# requeues itself, at most TC_FIELD_RESTARTS times, default 12) or week
# (pieces of up to 7 days; the week partition holds 96 CPUs per user).
# TC_FIELD_TIME (default the partition's limit) is the length of a piece.
# The ERW chemistry runs package z's fast kernels (the original code's
# numbers, bit for bit); TC_ERW_KERNEL may be unset or fast, nothing else.
# Stages (TC_FIELD_STAGE): fresh (default; result and log trees must be
# empty) or resume (adds runs or continues them; refuses while a requested
# run is queued or locked, or while a validation of this campaign is
# queued). Runs of the campaign not requested may keep running, so the
# pilot's other runs are submitted with resume while the pilot runs.

set -Eeuo pipefail

repo_root="$(pwd -P)"
campaign_dir="runtime/Field_SSP_2D"
campaign_batch="$campaign_dir/submit_field_ssp_2d.sbatch"
pin_list="runtime/Base_1D/ssp_pins.sha256"
equivalence_baseline="$repo_root/v2/SIM/EQUIV_1D_TMP/TC_ERW_results_Kantola_crop_control.mat"
minimum_free_gib=300
stage="${TC_FIELD_STAGE:-fresh}"
token='FIELD_SSP2D_V1'

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
: "${TC_SSP:?Set TC_SSP to ssp126, ssp245 or ssp585}"
: "${TC_FIELD_MEM:?Set TC_FIELD_MEM (for example 192G) from the pilot; see the README}"
case "$TC_SSP" in
    ssp126|ssp245|ssp585) ;;
    *) echo "Bad TC_SSP: $TC_SSP. No jobs were submitted." >&2; exit 2 ;;
esac
if [[ ! "$TC_FIELD_MEM" =~ ^[1-9][0-9]*G$ ]]; then
    echo "TC_FIELD_MEM must be a whole number of gigabytes such as 192G. No jobs were submitted." >&2
    exit 2
fi
calcite_case="${TC_CALCITE_CASE:-}"
case_suffix=''
allowed_runs='1 2 3 4 5 6 7 8 9'
default_runs='1-9'
if [[ -n "$calcite_case" ]]; then
    if [[ ! "$calcite_case" =~ ^[1-9][0-9]*$ ]]; then
        echo "Bad TC_CALCITE_CASE: $calcite_case. No jobs were submitted." >&2
        exit 2
    fi
    case_suffix="_case${calcite_case}"
    allowed_runs='3 6 9'
    default_runs='3,6,9'
fi
workers="${TC_FIELD_WORKERS:-32}"
if [[ ! "$workers" =~ ^[1-9][0-9]*$ ]] || (( workers < 2 || workers > 64 )); then
    echo "TC_FIELD_WORKERS must be a whole number from 2 to 64. No jobs were submitted." >&2
    exit 2
fi
restarts="${TC_FIELD_RESTARTS:-12}"
if [[ ! "$restarts" =~ ^[0-9]+$ ]] || (( restarts > 30 )); then
    echo "TC_FIELD_RESTARTS must be a whole number from 0 to 30. No jobs were submitted." >&2
    exit 2
fi
# Package z's fast ERW chemistry kernels, the original code's numbers bit
# for bit, are the only kernels a field run runs (MATLAB runs no MEX
# function on the thread pool every field run has). Every element is told
# the kernel explicitly.
kernel="${TC_ERW_KERNEL:-fast}"
case "$kernel" in
    fast) ;;
    mex) echo "TC_ERW_KERNEL=mex cannot run on a field run's thread pool; fast computes the same numbers. No jobs were submitted." >&2; exit 2 ;;
    *) echo "Field runs use the fast ERW kernels: TC_ERW_KERNEL must be fast or unset, not '$kernel'. No jobs were submitted." >&2; exit 2 ;;
esac
partition="${TC_FIELD_PARTITION:-day}"
case "$partition" in
    day)  limit_s=86400;  default_time='1-00:00:00' ;;
    week) limit_s=604800; default_time='7-00:00:00' ;;
    *) echo "TC_FIELD_PARTITION must be day or week, not '$partition'. No jobs were submitted." >&2; exit 2 ;;
esac
piece_time="${TC_FIELD_TIME:-$default_time}"
if [[ ! "$piece_time" =~ ^([0-9]+)-([0-9]{2}):([0-9]{2}):([0-9]{2})$ ]]; then
    echo "TC_FIELD_TIME must read D-HH:MM:SS. No jobs were submitted." >&2
    exit 2
fi
piece_s=$(( 10#${BASH_REMATCH[1]}*86400 + 10#${BASH_REMATCH[2]}*3600 + 10#${BASH_REMATCH[3]}*60 + 10#${BASH_REMATCH[4]} ))
if (( piece_s < 7200 || piece_s > limit_s )); then
    echo "TC_FIELD_TIME must lie between two hours and the $partition partition's limit. No jobs were submitted." >&2
    exit 2
fi
runs="${TC_FIELD_RUNS:-$default_runs}"
if [[ ! "$runs" =~ ^[1-9](-[1-9])?(,[1-9](-[1-9])?)*(%[1-9])?$ ]]; then
    echo "TC_FIELD_RUNS must be a Slurm array of runs, such as 2, 1,4,7 or 1-9%4. No jobs were submitted." >&2
    exit 2
fi
# The runs the array names, each once and each a run of this campaign.
requested=' '
IFS=',' read -r -a run_parts <<< "${runs%%[%]*}"
for part in "${run_parts[@]}"; do
    first="${part%-*}"; last="${part#*-}"
    if (( first > last )); then
        echo "TC_FIELD_RUNS has a descending range ($part). No jobs were submitted." >&2
        exit 2
    fi
    for (( r = first; r <= last; r++ )); do
        if [[ " $allowed_runs " != *" $r "* || "$requested" == *" $r "* ]]; then
            echo "TC_FIELD_RUNS must name distinct runs of $allowed_runs; $r is not one. No jobs were submitted." >&2
            exit 2
        fi
        requested="$requested$r "
    done
done

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
tag="Field_SSP_2D_35yr_${TC_SSP}${case_suffix}__${topology_name}"
results_dir="$repo_root/v2/SIM/Subfield_12/$tag"
logs_dir="$repo_root/v2/logs/Subfield_12/$tag"

required_files=(
    "$campaign_batch"
    "$campaign_dir/Run_Field_SSP_2D_Task.m"
    "$campaign_dir/Field_Campaign_Grid.m"
    "$campaign_dir/Check_Lateral_Conservation.m"
    "$campaign_dir/README.md"
    "$pin_list"
    'runtime/Base_1D/Campaign_Parameters.m'
    'runtime/Base_1D/Calcite_Case.m'
    'runtime/Base_1D/Assert_Base_Provenance.m'
    'runtime/Base_1D/Fingerprint_Digest.m'
    'runtime/Uncertainty_1D/requeue_near_limit.sh'
    'v2/T&C_Code/MAIN_FRAME_BG_ERW_2D.m'
    'v2/T&C_Code/Bouchet_Campaign_Validation_Fingerprint.m'
    'v2/T&C_Code/Bouchet_Check_Existing_Result.m'
    "v2/T&C_Code/Inputs/Data_Kantola_${TC_SSP}.mat"
    'v2/T&C_Code/Inputs/Data_Kantola_BG.mat'
    'v2/T&C_Code/Inputs/ERW_MOD_PARAM_parameters.xlsx'
    'v2/T&C_Code/Tests/run_post_legacy_equation_audit_tests.m'
    'v2/run/Run_Matrix_12.m'
    'v2/run/Run_Subfield_Worker.m'
    'v2/run/Run_One_Subfield.m'
    'v2/run/Bounded_File_Tag.m'
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
            echo 'A fresh submission needs empty campaign trees; use TC_FIELD_STAGE=resume to add or continue runs.' >&2
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
for check_partition in day "$partition"; do
    if ! partition_features="$(sinfo -h -p "$check_partition" -o '%f' 2>/dev/null)" || \
            [[ "$partition_features" != *cpumodel:9575f* ]]; then
        echo "The $check_partition partition does not expose cpumodel:9575f nodes. No jobs were submitted." >&2
        exit 3
    fi
done

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

# Job names carry the scenario, case and grid, so campaigns do not see each
# other. Within this campaign a queued validation refuses the call, and so
# does a queued element of a requested run; elements of other runs may run
# on, which is how the runs after the pilot are added.
suffix="${TC_SSP}${case_suffix}_${topology_name}"
validation_name="FieldSSP2Dchk_${suffix}"
production_name="FieldSSP2D_${suffix}"
if ! queued_jobs="$(squeue -h -u "$USER" -r --name="$validation_name,$production_name" -o '%i|%j|%T')"; then
    echo 'Cannot check for an active campaign graph. No jobs were submitted.' >&2
    exit 5
fi
conflicts=()
while IFS='|' read -r job_id job_name job_state; do
    [[ -n "$job_id" ]] || continue
    if [[ "$job_name" == "$validation_name" ]]; then
        conflicts+=("$job_id $job_name $job_state (validation)")
        continue
    fi
    [[ "$job_name" == "$production_name" ]] || continue
    array_task="${job_id##*_}"
    if [[ "$array_task" == "$job_id" || ! "$array_task" =~ ^[0-9]+$ || \
          "$requested" == *" $array_task "* ]]; then
        conflicts+=("$job_id $job_name $job_state")
    fi
done <<< "$queued_jobs"
if (( ${#conflicts[@]} > 0 )); then
    echo 'Jobs of this campaign are queued for a requested run or for validation. No jobs were submitted.' >&2
    printf '  %s\n' "${conflicts[@]}" >&2
    exit 5
fi
stale_locks=()
for r in $requested; do
    lock_dir="$logs_dir/locks/$(printf 'run_%02d.lock' "$r")"
    if [[ -d "$lock_dir" ]]; then
        stale_locks+=("$lock_dir")
    fi
done
if (( ${#stale_locks[@]} > 0 )); then
    echo 'Stale run locks remain from killed jobs. No job of these runs is queued, so remove them and rerun:' >&2
    for lock_dir in "${stale_locks[@]}"; do
        printf '  rm -r -- %q\n' "$lock_dir" >&2
    done
    echo 'No jobs were submitted.' >&2
    exit 5
fi
if [[ "$stage" == resume ]]; then
    echo 'Resume: validation rewrites its stamp (identical when nothing changed); finished runs are kept, interrupted runs continue from their checkpoints.'
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

export_common="ALL,TC_FIELD_SSP_CONTROLLER_TOKEN=$token,TC_FIELD_TOPOLOGY=$topology_file,TC_SSP=$TC_SSP,TC_CALCITE_CASE=$calcite_case,TC_FIELD_WORKERS=$workers,TC_FIELD_RESTARTS=$restarts,TC_ERW_KERNEL=$kernel"
submit_job validation_job \
    --partition=day --array=1 --ntasks=1 --cpus-per-task="$workers" --mem=96G --time=06:00:00 \
    --no-requeue --job-name="$validation_name" \
    --export="$export_common,TC_VALIDATE_ONLY=1" "$campaign_batch"
submit_job production_job \
    --dependency="afterok:${validation_job}" --kill-on-invalid-dep=yes \
    --array="$runs" --partition="$partition" --cpus-per-task="$workers" \
    --mem="$TC_FIELD_MEM" --time="$piece_time" \
    --job-name="$production_name" \
    --export="$export_common" "$campaign_batch"

scontrol release "${submitted_jobs[@]}"
graph_complete=1
trap - ERR HUP INT TERM

printf 'Campaign: %s\n' "$tag"
printf 'Validation: %s (day, %d CPUs, 6 h)\n' "$validation_job" "$workers"
printf 'Production: runs %s %s (%s, %d CPUs, %s, pieces of %s, up to %d requeues, ERW kernel %s), eligible after the validation passes\n' \
    "$runs" "$production_job" "$partition" "$workers" "$TC_FIELD_MEM" "$piece_time" "$restarts" "$kernel"
printf 'Monitor: sacct -j %s,%s --format=JobID%%18,JobName%%34,State%%12,Elapsed,MaxRSS,NodeList%%16\n' \
    "$validation_job" "$production_job"
