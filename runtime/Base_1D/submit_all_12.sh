#!/usr/bin/env bash
# Submit the calibrated CEC 13 all-12 Kantola campaign on Bouchet:
# validation (1-D audit, pre-application identity, smoke; 2-D wrapper)
# followed by production of runs 10-12 (one node) and runs 1-9 (array).
# The parameter set is alpha 3.4 (as the sweep forms it), KNmaxd 0.0048,
# free drainage (Kbot NaN) with CEC 13 and the Appelo exchange scheme; see
# Campaign_Parameters.m and the campaign README.
#
# Usage, from the checkout root:
#   MATLAB_MODULE=MATLAB/2025a TC_BOUCHET_ACCOUNT=pi_np363 \
#   bash runtime/Base_1D/submit_all_12.sh
# Add TC_CALCITE_CASE=8 to run the calcite treatments (1-D run 12 and 2-D
# runs 3, 6 and 9) from ERW workbook case 8, the agronomic 5 t/ha rate, as a
# campaign of its own (Base_1D_case8, Base_2D_case8) beside the standard one.
# Stages (TC_ALL12_STAGE): fresh (default; both campaign result and log
# trees must be empty) or resume (after TIMEOUT or node loss; refuses while
# any campaign job is queued and while a run lock remains).
# Dimensions (TC_ALL12_DIMS): all (default) or 1d, which submits the 1-D
# validation and runs 10-12 (run 12 for a calcite case) alone. The 2-D half
# then follows from the same command with TC_ALL12_STAGE=resume and no
# TC_ALL12_DIMS: the finished 1-D runs are kept.

set -Eeuo pipefail

repo_root="$(pwd -P)"
campaign_name_1d='Base_1D'
campaign_name_2d='Base_2D'
campaign_1d="runtime/$campaign_name_1d/submit_base_1d.sbatch"
campaign_2d="runtime/$campaign_name_2d/submit_base_2d.sbatch"
calcite_case="${TC_CALCITE_CASE:-}"
case_suffix=''
prod_1d_ntasks=3; prod_1d_runs='10-12'
prod_2d_array='1-9%9'; prod_2d_runs='1-9'
if [[ -n "$calcite_case" ]]; then
    if [[ ! "$calcite_case" =~ ^[1-9][0-9]*$ ]]; then
        echo "TC_CALCITE_CASE must be a positive integer, not '$calcite_case'. No jobs were submitted." >&2
        exit 2
    fi
    case_suffix="_case${calcite_case}"
    prod_1d_ntasks=1; prod_1d_runs='12'
    prod_2d_array='3,6,9'; prod_2d_runs='3,6,9'
fi
tag_1d="Base_1D${case_suffix}"
tag_2d="Base_2D${case_suffix}"
results_1d="$repo_root/v2/SIM/Subfield_12/$tag_1d"
results_2d="$repo_root/v2/SIM/Subfield_12/$tag_2d"
logs_1d="$repo_root/v2/logs/Subfield_12/$tag_1d"
logs_2d="$repo_root/v2/logs/Subfield_12/$tag_2d"
equivalence_baseline="$repo_root/v2/SIM/EQUIV_1D_TMP/TC_ERW_results_Kantola_crop_control.mat"
pin_list="runtime/$campaign_name_1d/shared_pins.sha256"
minimum_free_gib=40
stage="${TC_ALL12_STAGE:-fresh}"
dims="${TC_ALL12_DIMS:-all}"

if [[ ! -f "$repo_root/v2/run/Run_One_Subfield.m" ]]; then
    echo "Run this controller from the T&C-SMEW-N v0.9 checkout root." >&2
    exit 2
fi
case "$stage" in
    fresh|resume) ;;
    *)
        echo "TC_ALL12_STAGE must be fresh or resume, not '$stage'. No jobs were submitted." >&2
        exit 2 ;;
esac
case "$dims" in
    all|1d) ;;
    *)
        echo "TC_ALL12_DIMS must be all or 1d, not '$dims'. No jobs were submitted." >&2
        exit 2 ;;
esac

required_files=(
    "$campaign_1d"
    "$campaign_2d"
    "$pin_list"
    "runtime/$campaign_name_1d/Campaign_Parameters.m"
    "runtime/$campaign_name_1d/Calcite_Case.m"
    "runtime/$campaign_name_1d/Assert_Base_Provenance.m"
    "runtime/$campaign_name_1d/Fingerprint_Digest.m"
    "runtime/$campaign_name_1d/Run_Base_1D_Task.m"
    "runtime/$campaign_name_1d/Validate_Base_1D.m"
    "runtime/$campaign_name_1d/Validate_Preapplication_Identity.m"
    "runtime/$campaign_name_1d/Validate_All12_Smoke.m"
    "runtime/$campaign_name_1d/test_all_12_configuration.m"
    "runtime/$campaign_name_1d/README.md"
    "runtime/$campaign_name_2d/Run_Base_2D_Task.m"
    "runtime/$campaign_name_2d/Check_Lateral_Conservation.m"
    "runtime/$campaign_name_2d/README.md"
    "runtime/Calibration_Sweep_1D/Sweep_Cases.m"
    'v2/T&C_Code/ERW_function.m'
    'v2/T&C_Code/Bouchet_Campaign_Validation_Fingerprint.m'
    'v2/T&C_Code/Bouchet_Check_Existing_Result.m'
    'v2/T&C_Code/Inputs/Data_Kantola.mat'
    'v2/T&C_Code/Inputs/Data_Kantola_BG.mat'
    'v2/T&C_Code/Inputs/ERW_MOD_PARAM_parameters.xlsx'
    'v2/T&C_Code/Parameters/ERW_MOD_PARAM_excel.m'
    'v2/T&C_Code/Parameters/MOD_PARAM_Kantola_Maize_Soy.m'
    'v2/T&C_Code/Tests/run_post_legacy_equation_audit_tests.m'
    'v2/T&C_Code/Tests/Fixtures/no3_immobilization_20190906.mat'
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
: "${TC_BOUCHET_ACCOUNT:?Set TC_BOUCHET_ACCOUNT to the pi account reported by slurm_checkup.sh}"

# The equivalence baseline is hashed into every campaign stamp, including
# the finished CEC 13 sweep's. It is reused, never rebuilt.
if [[ ! -f "$equivalence_baseline" ]]; then
    echo "The equivalence baseline is missing: $equivalence_baseline" >&2
    echo 'This campaign does not rebuild it (a rebuild would orphan the CEC 13 sweep results). No jobs were submitted.' >&2
    exit 6
fi
echo "Equivalence baseline: reused ($equivalence_baseline)"

# The checkout must be the validated CEC 13 state: the shared code the
# audit suite pins, the old all-12 campaign files that suite reads, and the
# workbook and parameter script that supply CEC 13 and the Appelo scheme.
if ! (cd "$repo_root" && LC_ALL=C sha256sum -c --quiet --strict "$pin_list"); then
    echo "The checkout differs from the validated CEC 13 state pinned in $pin_list." >&2
    echo 'No jobs were submitted.' >&2
    exit 6
fi
echo "Checkout pins: $(wc -l < "$repo_root/$pin_list") files match $pin_list"

for state_path in "$results_1d" "$results_2d" "$logs_1d" "$logs_2d"; do
    if [[ -L "$state_path" || ( -e "$state_path" && ! -d "$state_path" ) ]]; then
        echo "Campaign path is not a directory: $state_path" >&2
        echo 'No jobs were submitted.' >&2
        exit 4
    fi
    if [[ "$stage" == fresh && -d "$state_path" ]]; then
        existing_state="$(find "$state_path" -mindepth 1 -print -quit)"
        if [[ -n "$existing_state" ]]; then
            echo "Existing campaign state found under $state_path." >&2
            echo 'A fresh submission needs empty campaign trees; use TC_ALL12_STAGE=resume to continue an interrupted campaign.' >&2
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
# The 2-D production runs on the week partition under the same node
# constraint as everything else; refuse before submitting anything if that
# partition no longer carries such nodes.
if [[ "$dims" == all ]] && { ! week_features="$(sinfo -h -p week -o '%f' 2>/dev/null)" || [[ "$week_features" != *cpumodel:9575f* ]]; }; then
    echo 'The week partition does not expose cpumodel:9575f nodes. No jobs were submitted.' >&2
    exit 3
fi

minimum_free_kib=$((minimum_free_gib * 1024 * 1024))
if ! df_record="$(LC_ALL=C df -Pk "$repo_root")"; then
    echo "Cannot read free space for $repo_root. No jobs were submitted." >&2
    exit 4
fi
available_kib="$(printf '%s\n' "$df_record" | awk 'NR == 2 {print $4}')"
if [[ ! "$available_kib" =~ ^[0-9]+$ ]]; then
    echo "Cannot parse free space from df -Pk. No jobs were submitted." >&2
    printf '%s\n' "$df_record" >&2
    exit 4
fi
if (( available_kib < minimum_free_kib )); then
    printf 'Only %.1f GiB is free on the repository filesystem; %d GiB is required.\n' \
        "$(awk -v kib="$available_kib" 'BEGIN {print kib/1024/1024}')" \
        "$minimum_free_gib" >&2
    echo 'No jobs were submitted.' >&2
    exit 4
fi
printf 'Storage preflight: %.1f GiB free on the repository filesystem.\n' \
    "$(awk -v kib="$available_kib" 'BEGIN {print kib/1024/1024}')"

mkdir -p "$logs_1d" "$logs_2d"
controller_lock="$logs_1d/.submission_controller.lock"
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

job_names="Base1Dcheck${case_suffix},Base2Dcheck${case_suffix},Base1D${case_suffix},Base2D${case_suffix}"
if ! active_jobs="$(squeue -h -u "$USER" --name="$job_names" -o '%A|%j|%T')"; then
    echo 'Cannot check for an active campaign graph. No jobs were submitted.' >&2
    exit 5
fi
if [[ -n "$active_jobs" ]]; then
    echo 'An active calibrated all-12 graph already exists. No jobs were submitted.' >&2
    printf '%s\n' "$active_jobs" >&2
    exit 5
fi

# No graph job is active, so any per-run lock left under the campaign logs
# belongs to a job that was killed before MATLAB could release it. A
# resubmitted element would fail with Bouchet:RunLocked, so refuse until the
# stale locks are removed by hand.
stale_locks="$(find "$logs_1d" "$logs_2d" -type d -name 'run_*.lock' 2>/dev/null || true)"
if [[ -n "$stale_locks" ]]; then
    echo 'Stale run locks remain from killed jobs. No graph job is active, so remove them and rerun:' >&2
    while IFS= read -r lock_dir; do
        [[ -n "$lock_dir" ]] || continue
        printf '  rm -r -- %q\n' "$lock_dir" >&2
    done <<< "$stale_locks"
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
        echo "Submitted job IDs: ${submitted_jobs[*]}" >&2
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

export_common="ALL,TC_ALL12_CONTROLLER_TOKEN=CEC13_ALL12_V1,TC_CALCITE_CASE=$calcite_case"
submit_job validation_1d_job \
    --partition=day --ntasks=1 --cpus-per-task=6 --mem=40G --time=04:00:00 \
    --job-name="Base1Dcheck${case_suffix}" \
    --export="$export_common,TC_VALIDATE_ONLY=1" "$campaign_1d"
if [[ "$dims" == all ]]; then
    submit_job validation_2d_job \
        --partition=day --array=1 --ntasks=1 --cpus-per-task=6 --mem=32G --time=04:00:00 \
        --job-name="Base2Dcheck${case_suffix}" \
        --export="$export_common,TC_VALIDATE_ONLY=1" "$campaign_2d"
    production_dependency="afterok:${validation_1d_job}:${validation_2d_job}"
else
    # 1-D first: the 1-D validation holds the audit, the pre-application
    # identity and both smoke runs.
    production_dependency="afterok:${validation_1d_job}"
fi
submit_job production_1d_job \
    --dependency="$production_dependency" \
    --kill-on-invalid-dep=yes --ntasks="$prod_1d_ntasks" \
    --job-name="Base1D${case_suffix}" \
    --export="$export_common,TC_VALIDATE_ONLY=0" "$campaign_1d"
if [[ "$dims" == all ]]; then
    submit_job production_2d_job \
        --dependency="$production_dependency" \
        --kill-on-invalid-dep=yes --array="$prod_2d_array" \
        --job-name="Base2D${case_suffix}" \
        --export="$export_common" "$campaign_2d"
fi

scontrol release "${submitted_jobs[@]}"
graph_complete=1
trap - ERR HUP INT TERM

if [[ "$dims" == all ]]; then
    printf 'Campaign: %s and %s\n' "$tag_1d" "$tag_2d"
    printf 'Validations: 1-D %s, 2-D %s\n' "$validation_1d_job" "$validation_2d_job"
    printf 'Production: runs %s %s (day, 6 h), runs %s %s (week, array %s, 4 days each)\n' \
        "$prod_1d_runs" "$production_1d_job" "$prod_2d_runs" "$production_2d_job" "$prod_2d_array"
    printf 'Both production jobs become eligible after both validations pass. Slurm controls their start times.\n'
    printf 'Monitor: sacct -j %s,%s,%s,%s -X --format=JobID%%16,JobName%%18,State%%12,Elapsed,NodeList%%16\n' \
        "$validation_1d_job" "$validation_2d_job" "$production_1d_job" "$production_2d_job"
else
    printf 'Campaign: %s (1-D only; %s follows with TC_ALL12_STAGE=resume)\n' "$tag_1d" "$tag_2d"
    printf 'Validation: 1-D %s\n' "$validation_1d_job"
    printf 'Production: runs %s %s (day, 6 h), eligible after the validation passes.\n' \
        "$prod_1d_runs" "$production_1d_job"
    printf 'Monitor: sacct -j %s,%s -X --format=JobID%%16,JobName%%18,State%%12,Elapsed,NodeList%%16\n' \
        "$validation_1d_job" "$production_1d_job"
fi
