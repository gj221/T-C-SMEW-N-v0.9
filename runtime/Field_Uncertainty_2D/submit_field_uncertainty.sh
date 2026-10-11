#!/usr/bin/env bash
# Submit the lower and upper sets of the Kantola uncertainty ensemble on a
# field grid for one horizon, scenario and calcite case: one validation
# element (the complete audit, the Linux checks and one day of the grid
# under each set) followed by an array over the sets and the runs of
# TC_FIELD_RUNS (member-major: every run under the lower set, then under
# the upper). The original set is not run: it is the central field
# campaign's result (runtime/Field_SSP_2D) in the checkout
# TC_UC_ORIGINAL_ROOT, whose model must be byte for byte this one's;
# "$PWD" when both run in this root. Without TC_SSP the sets run the
# five-year field record on the legacy forcing (rock to 2019, as
# runtime/Field_2D), with it the 35-year projection; TC_HORIZON may name
# the horizon (5yr or 35yr) or stay unset.
#
# Usage, from the checkout root:
#   TC_UC_ORIGINAL_ROOT="$PWD" \
#   TC_FIELD_TOPOLOGY=/full/path/beside/the/root/<name>.mat \
#   TC_SSP=ssp245 TC_FIELD_RUNS=2,3,5,6,8,9 TC_FIELD_WORKERS=32 TC_FIELD_MEM=96G \
#   MATLAB_MODULE=MATLAB/2025a TC_BOUCHET_ACCOUNT=<pi account> \
#   bash runtime/Field_Uncertainty_2D/submit_field_uncertainty.sh
# Add TC_CALCITE_CASE=8 for the 5 t/ha arm (runs 3, 6 and 9). TC_FIELD_RUNS
# (default every run of the arm) names runs as the field SSP controller
# does; a %N suffix holds N elements at once. TC_FIELD_MEM, TC_FIELD_WORKERS,
# TC_FIELD_PARTITION, TC_FIELD_TIME, TC_FIELD_RESTARTS and TC_ERW_KERNEL are
# the field SSP controller's. After a TIMEOUT or node loss, resubmit the same command once
# the call's jobs have left the queue: finished elements exit at once and
# the rest continue from their checkpoints.

set -Eeuo pipefail

repo_root="$(pwd -P)"
campaign_dir='runtime/Field_Uncertainty_2D'
campaign_batch="$campaign_dir/submit_field_uncertainty.sbatch"
pin_list='runtime/Base_1D/ssp_pins.sha256'
equivalence_baseline="$repo_root/v2/SIM/EQUIV_1D_TMP/TC_ERW_results_Kantola_crop_control.mat"
token='FIELD_UNCERTAINTY_V1'
run_members=2   # lower and upper (Uncertainty_Members.m); the original is reused
minimum_free_gib=600

if [[ ! -f "$repo_root/v2/run/Run_One_Subfield.m" || ! -f "$repo_root/$campaign_batch" ]]; then
    echo "Run this from the T&C-SMEW-N v0.9 checkout root." >&2
    exit 2
fi
: "${TC_BOUCHET_ACCOUNT:?Set TC_BOUCHET_ACCOUNT to the pi account reported by slurm_checkup.sh}"
: "${TC_UC_ORIGINAL_ROOT:?Set TC_UC_ORIGINAL_ROOT to the checkout whose field SSP campaigns hold the original set}"
: "${TC_FIELD_TOPOLOGY:?Set TC_FIELD_TOPOLOGY to the full path of the grid file}"
: "${TC_FIELD_MEM:?Set TC_FIELD_MEM (for example 96G) from the pilot; see the README}"
# The horizon, as the task reads it (Field_Campaign_Grid 'horizon'): the
# five-year field record without a scenario, the 35-year projection with one.
ssp="${TC_SSP:-}"
horizon="${TC_HORIZON:-}"
if [[ -z "$ssp" ]]; then
    if [[ -n "$horizon" && "$horizon" != 5yr ]]; then
        echo "Without TC_SSP the ensemble is the five-year field record; TC_HORIZON=$horizon is not it. No jobs were submitted." >&2
        exit 2
    fi
    horizon='5yr'; horizon_tag='5yr'; forcing='Data_Kantola.mat'
else
    case "$ssp" in
        ssp126|ssp245|ssp585) ;;
        *) echo "Bad TC_SSP: $ssp. No jobs were submitted." >&2; exit 2 ;;
    esac
    if [[ -n "$horizon" && "$horizon" != 35yr ]]; then
        echo "A scenario ensemble is the 35-year projection; TC_HORIZON=$horizon is not it. No jobs were submitted." >&2
        exit 2
    fi
    horizon='35yr'; horizon_tag="35yr_${ssp}"; forcing="Data_Kantola_${ssp}.mat"
fi
if [[ ! "$TC_FIELD_MEM" =~ ^[1-9][0-9]*G$ ]]; then
    echo "TC_FIELD_MEM must be a whole number of gigabytes such as 96G. No jobs were submitted." >&2
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
if [[ ! "$runs" =~ ^[1-9](-[1-9])?(,[1-9](-[1-9])?)*(%[1-9][0-9]*)?$ ]]; then
    echo "TC_FIELD_RUNS must name runs such as 2, 1,4,7 or 2-9%4. No jobs were submitted." >&2
    exit 2
fi
throttle=''
if [[ "$runs" == *%* ]]; then
    throttle="%${runs##*[%]}"
fi
# The runs, each once and each a run of this arm, as the comma list the
# task reads.
requested=' '
run_list=''
n_runs=0
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
        run_list="${run_list:+$run_list,}$r"
        n_runs=$(( n_runs + 1 ))
    done
done

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
tag="Field_Uncertainty_2D_${horizon_tag}${case_suffix}__${topology_name}"
results_dir="$repo_root/v2/SIM/Subfield_12/$tag"
logs_dir="$repo_root/v2/logs/Subfield_12/$tag"

required_files=(
    "$campaign_batch"
    "$campaign_dir/Run_Field_Uncertainty_Task.m"
    "$campaign_dir/Field_Uncertainty_Campaign.m"
    "$campaign_dir/README.md"
    "$pin_list"
    'runtime/Field_SSP_2D/Field_Campaign_Grid.m'
    'runtime/Field_SSP_2D/Check_Lateral_Conservation.m'
    'runtime/Uncertainty_1D/Uncertainty_Members.m'
    'runtime/Uncertainty_1D/Uncertainty_Knobs.m'
    'runtime/Uncertainty_1D/Uncertainty_Result_Matches.m'
    'runtime/Uncertainty_1D/Uncertainty_Result_State.m'
    'runtime/Uncertainty_1D/same_model_code.sh'
    'runtime/Uncertainty_1D/requeue_near_limit.sh'
    'runtime/Base_1D/Campaign_Parameters.m'
    'runtime/Base_1D/Calcite_Case.m'
    'runtime/Base_1D/Fingerprint_Digest.m'
    "v2/T&C_Code/Inputs/$forcing"
    'v2/T&C_Code/Tests/run_post_legacy_equation_audit_tests.m'
    'v2/run/Run_One_Subfield.m'
    'v2/run/Run_Subfield_Worker.m'
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
    echo 'Build it for this deployment first (regen_equiv_1d.sbatch). No jobs were submitted.' >&2
    exit 6
fi
if ! (cd "$repo_root" && LC_ALL=C sha256sum -c --quiet --strict "$pin_list"); then
    echo "The checkout differs from the state pinned in $pin_list. No jobs were submitted." >&2
    exit 6
fi
echo "Checkout pins: $(wc -l < "$repo_root/$pin_list") files match $pin_list"
echo "Grid: $topology_file ($(sha256sum "$topology_file" | cut -c1-12))"
bash runtime/Uncertainty_1D/same_model_code.sh "$TC_UC_ORIGINAL_ROOT" || {
    echo 'No jobs were submitted.' >&2; exit 6; }

for state_path in "$results_dir" "$logs_dir"; do
    if [[ -L "$state_path" || ( -e "$state_path" && ! -d "$state_path" ) ]]; then
        echo "Campaign path is not a directory: $state_path" >&2
        echo 'No jobs were submitted.' >&2
        exit 4
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

if ! available_kib="$(LC_ALL=C df -Pk "$repo_root" | awk 'NR == 2 {print $4}')" || \
        [[ ! "$available_kib" =~ ^[0-9]+$ ]]; then
    echo "Cannot read the free space of $repo_root. No jobs were submitted." >&2
    exit 4
fi
if (( available_kib < minimum_free_gib * 1024 * 1024 )); then
    echo "Only $(( available_kib / 1024 / 1024 )) GiB free; this call needs $minimum_free_gib GiB. No jobs were submitted." >&2
    exit 4
fi

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

suffix="${ssp:-5yr}${case_suffix}_${topology_name}"
validation_name="FieldUC2Dchk_${suffix}"
production_name="FieldUC2D_${suffix}"
if ! active_jobs="$(squeue -h -u "$USER" --name="$validation_name,$production_name" -o '%A|%j|%T')"; then
    echo 'Cannot check for active ensemble jobs. No jobs were submitted.' >&2
    exit 5
fi
if [[ -n "$active_jobs" ]]; then
    echo "Ensemble jobs for $tag are still queued or running. No jobs were submitted." >&2
    printf '%s\n' "$active_jobs" >&2
    exit 5
fi
stale_locks="$(find "$logs_dir" -type d -name 'run_*.lock' 2>/dev/null || true)"
if [[ -n "$stale_locks" ]]; then
    echo 'Stale run locks remain from killed jobs. No job of this ensemble is queued, so remove them and rerun:' >&2
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

n_elements=$(( run_members * n_runs ))
# The run list holds commas, which separate the variables of --export, so
# it reaches the jobs through the environment --export=ALL passes on.
export TC_FIELD_RUNS="$run_list"
export_common="ALL,TC_FIELD_UC_CONTROLLER_TOKEN=$token,TC_FIELD_TOPOLOGY=$topology_file,TC_SSP=$ssp,TC_HORIZON=$horizon,TC_CALCITE_CASE=$calcite_case,TC_FIELD_WORKERS=$workers,TC_FIELD_RESTARTS=$restarts,TC_ERW_KERNEL=$kernel"
submit_job validation_job \
    --partition=day --array=1 --ntasks=1 --cpus-per-task="$workers" --mem=96G --time=06:00:00 \
    --no-requeue --job-name="$validation_name" \
    --export="$export_common,TC_VALIDATE_ONLY=1" "$campaign_batch"
submit_job production_job \
    --dependency="afterok:${validation_job}" --kill-on-invalid-dep=yes \
    --array="1-${n_elements}${throttle}" --partition="$partition" --cpus-per-task="$workers" \
    --mem="$TC_FIELD_MEM" --time="$piece_time" \
    --job-name="$production_name" \
    --export="$export_common" "$campaign_batch"

scontrol release "${submitted_jobs[@]}"
graph_complete=1
trap - ERR HUP INT TERM

printf 'Uncertainty ensemble %s (original %s in %s)\n' "$tag" \
    "Field_SSP_2D_${horizon_tag}${case_suffix}__${topology_name}" "$TC_UC_ORIGINAL_ROOT"
printf 'Validation: %s (day, %d CPUs, 6 h)\n' "$validation_job" "$workers"
printf 'Production: %s, %d elements (lower then upper over runs %s; %s, %d CPUs, %s, pieces of %s, up to %d requeues, ERW kernel %s), eligible after the validation passes\n' \
    "$production_job" "$n_elements" "$run_list" "$partition" "$workers" "$TC_FIELD_MEM" \
    "$piece_time" "$restarts" "$kernel"
printf 'Monitor: sacct -j %s,%s --format=JobID%%18,JobName%%40,State%%12,Elapsed,MaxRSS,NodeList%%16\n' \
    "$validation_job" "$production_job"
