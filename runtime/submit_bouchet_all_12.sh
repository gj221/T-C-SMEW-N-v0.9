#!/usr/bin/env bash

set -Eeuo pipefail

repo_root="$(pwd -P)"
campaign_name_1d='Kantola_5yr_1D_Runs10to12_ScheduleFix_R2025a_Bouchet'
campaign_name_2d='Kantola_5yr_2D_Runs01to09_DenitPH7p8_Alpha1p3_Bouchet'
campaign_1d="runtime/$campaign_name_1d/submit_kantola_1d_runs_10_to_12.sbatch"
campaign_2d="runtime/$campaign_name_2d/submit_kantola_2d_runs_01_to_09.sbatch"
results_1d="$repo_root/v2/SIM/Subfield_12/$campaign_name_1d"
results_2d="$repo_root/v2/SIM/Subfield_12/$campaign_name_2d"
logs_1d="$repo_root/v2/logs/Subfield_12/$campaign_name_1d"
logs_2d="$repo_root/v2/logs/Subfield_12/$campaign_name_2d"

required_files=(
    "$campaign_1d"
    "$campaign_2d"
    'v2/T&C_Code/ERW_function.m'
    'v2/T&C_Code/Inputs/Data_Kantola.mat'
    'v2/T&C_Code/Inputs/Data_Kantola_BG.mat'
    'v2/T&C_Code/Parameters/ERW_MOD_PARAM_excel.m'
    'v2/T&C_Code/Tests/run_post_legacy_equation_audit_tests.m'
    'v2/T&C_Code/Tests/Fixtures/no3_immobilization_20190906.mat'
    'v2/run/Run_Matrix_12.m'
)
for required_file in "${required_files[@]}"; do
    if [[ ! -r "$repo_root/$required_file" ]]; then
        echo "Required campaign file is missing or unreadable: $repo_root/$required_file" >&2
        echo 'No jobs were submitted.' >&2
        exit 2
    fi
done
: "${TC_BOUCHET_ACCOUNT:?Set TC_BOUCHET_ACCOUNT to the pi account reported by slurm_checkup.sh}"

for state_path in "$results_1d" "$results_2d" "$logs_1d" "$logs_2d"; do
    if [[ -L "$state_path" || ( -e "$state_path" && ! -d "$state_path" ) ]]; then
        echo "Campaign path is not a fresh directory: $state_path" >&2
        echo 'Use a fresh deployment root. No jobs were submitted.' >&2
        exit 4
    fi
    if [[ -d "$state_path" ]]; then
        existing_state="$(find "$state_path" -mindepth 1 -print -quit)"
        if [[ -n "$existing_state" ]]; then
            echo "Existing campaign state found under $state_path." >&2
            echo 'Use a fresh deployment root so old results and checkpoints remain archived.' >&2
            echo 'No jobs were submitted.' >&2
            exit 4
        fi
    fi
done

unset TC_REBUILD_EQUIVALENCE TC_VALIDATE_ONLY
for submit_env in "${!SBATCH_@}"; do
    unset "$submit_env"
done
unset SLURM_CLUSTERS

export MATLAB_MODULE="${MATLAB_MODULE:-MATLAB/2025a}"
if ! module spider "$MATLAB_MODULE" >/dev/null 2>&1; then
    echo "Bouchet does not currently expose $MATLAB_MODULE. No jobs were submitted." >&2
    exit 3
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

submit_job baseline_job \
    --partition=day --ntasks=1 --cpus-per-task=1 --mem=16G --time=01:00:00 \
    --job-name=KantolaBase \
    --export=ALL,TC_REBUILD_EQUIVALENCE=1,TC_VALIDATE_ONLY=0 "$campaign_1d"

validation_dependency="afterok:${baseline_job}"
submit_job validation_1d_job \
    --dependency="$validation_dependency" \
    --kill-on-invalid-dep=yes \
    --partition=day --ntasks=1 --cpus-per-task=6 --mem=40G --time=04:00:00 \
    --job-name=Kantola1Dcheck \
    --export=ALL,TC_REBUILD_EQUIVALENCE=0,TC_VALIDATE_ONLY=1 "$campaign_1d"
submit_job validation_2d_job \
    --dependency="$validation_dependency" \
    --kill-on-invalid-dep=yes \
    --partition=day --array=1 --ntasks=1 --cpus-per-task=6 --mem=32G --time=04:00:00 \
    --job-name=Kantola2Dcheck \
    --export=ALL,TC_REBUILD_EQUIVALENCE=0,TC_VALIDATE_ONLY=1 "$campaign_2d"

production_dependency="afterok:${validation_1d_job}:${validation_2d_job}"
submit_job production_1d_job \
    --dependency="$production_dependency" \
    --kill-on-invalid-dep=yes \
    --export=ALL,TC_REBUILD_EQUIVALENCE=0,TC_VALIDATE_ONLY=0 "$campaign_1d"
submit_job production_2d_job \
    --dependency="$production_dependency" \
    --kill-on-invalid-dep=yes \
    --export=ALL,TC_REBUILD_EQUIVALENCE=0 "$campaign_2d"

scontrol release "${submitted_jobs[@]}"
graph_complete=1
trap - ERR HUP INT TERM

printf 'Baseline: %s\n' "$baseline_job"
printf 'Validations: 1-D %s, 2-D %s\n' "$validation_1d_job" "$validation_2d_job"
printf 'Production: runs 10-12 %s, runs 1-9 %s\n' "$production_1d_job" "$production_2d_job"
printf 'Both production jobs become eligible after both validations pass. Slurm controls their start times.\n'
