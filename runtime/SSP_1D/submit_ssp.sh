#!/usr/bin/env bash
# Submit ONE Kantola SSP campaign on Bouchet: validation (1-D audit and 2-D
# wrapper) followed by production of runs 10-12 and runs 1-9.
#
# One campaign per invocation, deliberately. The four campaigns -- the 5yr
# validation and the three 35yr projections -- are independent: different
# scenarios, different horizons, different result and log trees, no shared
# writable state. Submitting them separately lets Slurm run all four at once,
# so three projections cost the wall-clock of one rather than three. Nothing
# here serialises them; only the validation-before-production dependency
# WITHIN a campaign is enforced.
#
# Usage, from the checkout root:
#   MATLAB_MODULE=MATLAB/2025a TC_BOUCHET_ACCOUNT=pi_np363 \
#   TC_SSP=ssp245 TC_HORIZON=5yr \
#   bash runtime/SSP_1D/submit_ssp.sh
#
# Add TC_CALCITE_CASE=8 to run the calcite treatments (2-D runs 3, 6, 9 and
# 1-D run 12) from ERW workbook case 8, the agronomic 5 t/ha rate, as a
# campaign of its own. The 35yr horizon applies rock every 10 November to
# 2050; 5yr keeps the field record, November 2016-2019.
#
# Stages (TC_SSP_STAGE): fresh (default; both campaign trees must be empty)
# or resume (after TIMEOUT or node loss; refuses while any campaign job is
# queued and while a run lock remains).

set -Eeuo pipefail

repo_root="$(pwd -P)"
campaign_name_1d='SSP_1D'
campaign_name_2d='SSP_2D'
calibrated_1d='Base_1D'
campaign_1d="runtime/$campaign_name_1d/submit_ssp_1d.sbatch"
campaign_2d="runtime/$campaign_name_2d/submit_ssp_2d.sbatch"
pin_list="runtime/$calibrated_1d/ssp_pins.sha256"
equivalence_baseline="$repo_root/v2/SIM/EQUIV_1D_TMP/TC_ERW_results_Kantola_crop_control.mat"
minimum_free_gib=80
stage="${TC_SSP_STAGE:-fresh}"
token='SSP_CAMPAIGN_V1'

if [[ ! -f "$repo_root/v2/run/Run_One_Subfield.m" ]]; then
    echo "Run this controller from the T&C-SMEW-N v0.9 checkout root." >&2
    exit 2
fi
case "$stage" in
    fresh|resume) ;;
    *) echo "TC_SSP_STAGE must be fresh or resume, not '$stage'. No jobs were submitted." >&2; exit 2 ;;
esac

: "${TC_SSP:?Set TC_SSP to ssp126, ssp245 or ssp585}"
: "${TC_HORIZON:?Set TC_HORIZON to 5yr or 35yr}"
case "$TC_SSP" in ssp126|ssp245|ssp585) ;; *) echo "Bad TC_SSP: $TC_SSP" >&2; exit 2 ;; esac
export TC_SSP
calcite_case="${TC_CALCITE_CASE:-}"
case_suffix=''
prod_1d_ntasks=3; prod_1d_runs='10-12'
prod_2d_array='1-9'
if [[ -n "$calcite_case" ]]; then
    if [[ ! "$calcite_case" =~ ^[1-9][0-9]*$ ]]; then
        echo "Bad TC_CALCITE_CASE: $calcite_case" >&2; exit 2
    fi
    case_suffix="_case${calcite_case}"
    prod_1d_ntasks=1; prod_1d_runs='12'
    prod_2d_array='3,6,9'
fi

# Resources differ only by horizon. 35yr is roughly seven times the span, so
# it needs the week partition and more memory. The week partition lets one
# user hold 1.5 TB at once, so the 2-D figure sets how many runs go
# together: six at 256 G. Package u's 35-year whole-field basalt run had
# reached 205 GB by half way, most of it particle-class histories kept at
# every ten-minute row; from package v those keep the last three days and
# each day's last row (ERW_Particle_History), which puts that run near
# 150-200 GB. The controls need about 105 GB and can be lowered with
# scontrol update before they start; read a finished run's MaxRSS before
# lowering anything else.
#
# The 1-D job-level figure must exceed the sum of its three srun steps, which
# request run_mem each inside submit_ssp_1d.sbatch: 3 x 12 G at 5yr
# and 3 x 56 G at 35yr. 48 G and 192 G clear those with slack.
case "$TC_HORIZON" in
    5yr)
        prod_1d_partition=day;  prod_1d_time=06:00:00;     prod_1d_mem=48G
        prod_2d_partition=day;  prod_2d_time=12:00:00;     prod_2d_mem=64G
        ;;
    35yr)
        prod_1d_partition=week; prod_1d_time=1-12:00:00;   prod_1d_mem=192G
        prod_2d_partition=week; prod_2d_time=4-00:00:00;   prod_2d_mem=256G
        ;;
    *) echo "Bad TC_HORIZON: $TC_HORIZON" >&2; exit 2 ;;
esac
export TC_HORIZON

tag_1d="SSP_1D_${TC_HORIZON}_${TC_SSP}${case_suffix}"
tag_2d="SSP_2D_${TC_HORIZON}_${TC_SSP}${case_suffix}"
results_1d="$repo_root/v2/SIM/Subfield_12/$tag_1d"
results_2d="$repo_root/v2/SIM/Subfield_12/$tag_2d"
logs_1d="$repo_root/v2/logs/Subfield_12/$tag_1d"
logs_2d="$repo_root/v2/logs/Subfield_12/$tag_2d"

required_files=(
    "$campaign_1d"
    "$campaign_2d"
    "$pin_list"
    "runtime/$campaign_name_1d/Run_SSP_1D_Task.m"
    "runtime/$campaign_name_1d/Validate_SSP_1D.m"
    "runtime/$campaign_name_2d/Run_SSP_2D_Task.m"
    "runtime/$campaign_name_2d/Check_Lateral_Conservation.m"
    "runtime/$calibrated_1d/Campaign_Parameters.m"
    "runtime/$calibrated_1d/Assert_Base_Provenance.m"
    "runtime/$calibrated_1d/Fingerprint_Digest.m"
    "runtime/$calibrated_1d/test_all_12_configuration.m"
    'v2/T&C_Code/ERW_function.m'
    'v2/T&C_Code/Bouchet_Campaign_Validation_Fingerprint.m'
    'v2/T&C_Code/Bouchet_Check_Existing_Result.m'
    "v2/T&C_Code/Inputs/Data_Kantola_${TC_SSP}.mat"
    'v2/T&C_Code/Inputs/Data_Kantola_BG.mat'
    'v2/T&C_Code/Inputs/ERW_MOD_PARAM_parameters.xlsx'
    'v2/T&C_Code/Parameters/ERW_MOD_PARAM_excel.m'
    'v2/T&C_Code/Parameters/MOD_PARAM_Kantola_Maize_Soy.m'
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

if [[ ! -f "$equivalence_baseline" ]]; then
    echo "The equivalence baseline is missing: $equivalence_baseline" >&2
    echo 'Copy it from a root that has it; this campaign never rebuilds it.' >&2
    exit 6
fi
echo "Equivalence baseline: reused ($equivalence_baseline)"

# The checkout must be the state this campaign was pinned against. The pin
# list has to be regenerated whenever the deployed chemistry moves -- it is
# NOT the CEC13 list, which describes a different (earlier) stack.
if ! (cd "$repo_root" && LC_ALL=C sha256sum -c --quiet --strict "$pin_list"); then
    echo "The checkout differs from the state pinned in $pin_list." >&2
    echo 'No jobs were submitted.' >&2
    exit 6
fi
echo "Checkout pins: $(wc -l < "$repo_root/$pin_list") files match $pin_list"

for state_path in "$results_1d" "$results_2d" "$logs_1d" "$logs_2d"; do
    if [[ -L "$state_path" || ( -e "$state_path" && ! -d "$state_path" ) ]]; then
        echo "Campaign path is not a directory: $state_path" >&2
        exit 4
    fi
    if [[ "$stage" == fresh && -d "$state_path" ]]; then
        if [[ -n "$(find "$state_path" -mindepth 1 -print -quit)" ]]; then
            echo "Existing campaign state found under $state_path." >&2
            echo 'A fresh submission needs empty campaign trees; use TC_SSP_STAGE=resume to continue.' >&2
            exit 4
        fi
    fi
done

unset TC_REBUILD_EQUIVALENCE TC_VALIDATE_ONLY TC_RUN_ID TC_CKPT_FORCE_RESUME
for submit_env in "${!SBATCH_@}"; do unset "$submit_env"; done
unset SLURM_CLUSTERS

export MATLAB_MODULE="${MATLAB_MODULE:-MATLAB/2025a}"
if ! module spider "$MATLAB_MODULE" >/dev/null 2>&1; then
    echo "Bouchet does not currently expose $MATLAB_MODULE. No jobs were submitted." >&2
    exit 3
fi
if ! week_features="$(sinfo -h -p "$prod_2d_partition" -o '%f' 2>/dev/null)" || \
   [[ "$week_features" != *cpumodel:9575f* ]]; then
    echo "The $prod_2d_partition partition does not expose cpumodel:9575f nodes." >&2
    exit 3
fi

minimum_free_kib=$((minimum_free_gib * 1024 * 1024))
if ! df_record="$(LC_ALL=C df -Pk "$repo_root")"; then
    echo "Cannot read free space for $repo_root." >&2; exit 4
fi
available_kib="$(printf '%s\n' "$df_record" | awk 'NR == 2 {print $4}')"
if [[ ! "$available_kib" =~ ^[0-9]+$ ]]; then
    echo "Cannot parse free space from df -Pk." >&2; exit 4
fi
if (( available_kib < minimum_free_kib )); then
    printf 'Only %.1f GiB free; %d GiB required.\n' \
        "$(awk -v k="$available_kib" 'BEGIN {print k/1024/1024}')" "$minimum_free_gib" >&2
    exit 4
fi
printf 'Storage preflight: %.1f GiB free.\n' \
    "$(awk -v k="$available_kib" 'BEGIN {print k/1024/1024}')"

mkdir -p "$logs_1d" "$logs_2d"
controller_lock="$logs_1d/.submission_controller.lock"
if ! mkdir "$controller_lock" 2>/dev/null; then
    echo "Submission lock already exists: $controller_lock" >&2
    echo 'No jobs were submitted. Check for an active controller or a stale lock.' >&2
    exit 5
fi
release_controller_lock() {
    [[ -d "$controller_lock" ]] || return 0
    rmdir "$controller_lock" 2>/dev/null || \
        echo "Warning: could not remove submission lock $controller_lock" >&2
}
trap release_controller_lock EXIT

# Job names carry the scenario and horizon so that concurrent campaigns do
# not see each other as an active graph.
suffix="${TC_SSP}_${TC_HORIZON}${case_suffix}"
job_names="SSP1Dchk_$suffix,SSP2Dchk_$suffix,SSP1D_$suffix,SSP2D_$suffix"
if ! active_jobs="$(squeue -h -u "$USER" --name="$job_names" -o '%A|%j|%T')"; then
    echo 'Cannot check for an active campaign graph.' >&2; exit 5
fi
if [[ -n "$active_jobs" ]]; then
    echo "An active graph already exists for $suffix. No jobs were submitted." >&2
    printf '%s\n' "$active_jobs" >&2
    exit 5
fi

stale_locks="$(find "$logs_1d" "$logs_2d" -type d -name 'run_*.lock' 2>/dev/null || true)"
if [[ -n "$stale_locks" ]]; then
    echo 'Stale run locks remain from killed jobs. Remove them and rerun:' >&2
    while IFS= read -r lock_dir; do
        [[ -n "$lock_dir" ]] || continue
        printf '  rm -r -- %q\n' "$lock_dir" >&2
    done <<< "$stale_locks"
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
        scancel "${submitted_jobs[@]}" || \
            echo "Automatic cancellation failed. Run: scancel ${submitted_jobs[*]}" >&2
    fi
    release_controller_lock
    exit "$status"
}
trap 'rollback $?' ERR
trap 'rollback 129' HUP
trap 'rollback 130' INT
trap 'rollback 143' TERM

submit_job() {
    local result_var="$1"; shift
    local submitted submitted_id
    submitted="$(sbatch --parsable --hold --account="$TC_BOUCHET_ACCOUNT" "$@")"
    submitted_id="${submitted%%;*}"
    printf -v "$result_var" '%s' "$submitted_id"
    submitted_jobs+=("$submitted_id")
}

# Validation needs 12 CPUs on both sides: the 2-D task resolves twelve
# workers before it reaches the validate-only branch, and the 1-D validator
# asserts at least six.
export_common="ALL,TC_SSP_CONTROLLER_TOKEN=$token,TC_SSP=$TC_SSP,TC_HORIZON=$TC_HORIZON,TC_CALCITE_CASE=$calcite_case"
submit_job validation_1d_job \
    --partition=day --ntasks=1 --cpus-per-task=12 --mem=48G --time=04:00:00 \
    --job-name="SSP1Dchk_$suffix" \
    --export="$export_common,TC_VALIDATE_ONLY=1" "$campaign_1d"
submit_job validation_2d_job \
    --partition=day --array=1 --ntasks=1 --cpus-per-task=12 --mem=48G --time=04:00:00 \
    --job-name="SSP2Dchk_$suffix" \
    --export="$export_common,TC_VALIDATE_ONLY=1" "$campaign_2d"

production_dependency="afterok:${validation_1d_job}:${validation_2d_job}"
submit_job production_1d_job \
    --dependency="$production_dependency" --kill-on-invalid-dep=yes \
    --partition="$prod_1d_partition" --time="$prod_1d_time" --mem="$prod_1d_mem" \
    --nodes=1 --ntasks="$prod_1d_ntasks" --cpus-per-task=1 \
    --job-name="SSP1D_$suffix" \
    --export="$export_common,TC_VALIDATE_ONLY=0" "$campaign_1d"
submit_job production_2d_job \
    --dependency="$production_dependency" --kill-on-invalid-dep=yes \
    --partition="$prod_2d_partition" --time="$prod_2d_time" --mem="$prod_2d_mem" \
    --array="$prod_2d_array" --nodes=1 --ntasks=1 --cpus-per-task=12 \
    --job-name="SSP2D_$suffix" \
    --export="$export_common" "$campaign_2d"

scontrol release "${submitted_jobs[@]}"
graph_complete=1
trap - ERR HUP INT TERM

printf '\nCampaign: %s %s%s\n' "$TC_SSP" "$TC_HORIZON" "$case_suffix"
printf 'Validations: 1-D %s, 2-D %s\n' "$validation_1d_job" "$validation_2d_job"
printf 'Production:  runs %s %s (%s, %s, %s)\n' \
    "$prod_1d_runs" "$production_1d_job" "$prod_1d_partition" "$prod_1d_time" "$prod_1d_mem"
printf '             runs %s %s (%s, array %s, %s, %s, 12 workers)\n' \
    "$prod_2d_array" "$production_2d_job" "$prod_2d_partition" "$prod_2d_array" "$prod_2d_time" "$prod_2d_mem"
printf 'Both production jobs become eligible once both validations pass.\n'
printf 'This campaign is independent of the others: run this script again\n'
printf 'with a different TC_SSP/TC_HORIZON to add a concurrent campaign.\n'
printf 'Monitor: sacct -j %s,%s,%s,%s -X --format=JobID%%16,JobName%%18,State%%12,Elapsed,MaxRSS\n' \
    "$validation_1d_job" "$validation_2d_job" "$production_1d_job" "$production_2d_job"
