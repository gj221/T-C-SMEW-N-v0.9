#!/usr/bin/env bash

set -Eeuo pipefail

repo_root="$(pwd -P)"
campaign_name='Calibration_Sweep_1D'
campaign_dir="runtime/$campaign_name"
batch_file="$campaign_dir/submit_sweep.sbatch"
required_files=(
    "$batch_file"
    "$campaign_dir/Sweep_Cases.m"
    "$campaign_dir/Pre_Control_pH_Score.m"
    "$campaign_dir/Run_Sweep_Task.m"
    "$campaign_dir/Validate_Sweep.m"
    "$campaign_dir/Summarize_Sweep.m"
    "$campaign_dir/Load_Sweep_Observations.m"
    "$campaign_dir/Sweep_Fingerprint.m"
    "$campaign_dir/Fingerprint_Digest.m"
    "$campaign_dir/Assert_Sweep_Provenance.m"
    "$campaign_dir/README.md"
    'v2/T&C_Code/BIOGEOCHEMISTRY_DYNAMIC3.m'
    'v2/T&C_Code/Bouchet_Acquire_Run_Lock.m'
    'v2/T&C_Code/Bouchet_Campaign_Validation_Fingerprint.m'
    'v2/T&C_Code/Bouchet_Check_Existing_Result.m'
    'v2/T&C_Code/Bouchet_Single_Thread_Runtime.m'
    'v2/T&C_Code/Bouchet_Write_Validation_Stamp.m'
    'v2/T&C_Code/ERW_Checkpoint_Fingerprint.m'
    'v2/T&C_Code/ERW_Equivalence_Fingerprint.m'
    'v2/T&C_Code/Inputs/ERW_MOD_PARAM_parameters.xlsx'
    'v2/T&C_Code/Inputs/Data_Kantola.mat'
    'v2/T&C_Code/Inputs/Data_Kantola_BG.mat'
    'v2/T&C_Code/Inputs/Ca_Data.mat'
    'v2/T&C_Code/Inputs/All_deposition_data.mat'
    'v2/T&C_Code/Inputs/blueridge_CDF_data.json'
    'v2/T&C_Code/Tests/run_post_legacy_equation_audit_tests.m'
    'v2/T&C_Code/Tests/Kantola_N_Magnitude_Audit.m'
    'v2/T&C_Code/Tests/test_bouchet_1d_joint_sweep_configuration.m'
    'v2/T&C_Code/Tests/test_bouchet_linux_portability.m'
    'v2/Experimental observations/Web plot digitizer/Data_Beerling.mat'
    'v2/Experimental observations/Web plot digitizer/figure3a_n2o_digitised.csv'
    'v2/Experimental observations/Web plot digitizer/Kantola_maize_extracted_data.mat'
    'v2/Experimental observations/Web plot digitizer/Kantola_soybean_extracted_data.mat'
    'v2/run/ERW_run.m'
    'v2/run/Run_One_Subfield.m'
    'v2/run/Run_All_12_Parallel.m'
    'v2/run/Run_Subfield_Worker.m'
)
for required_file in "${required_files[@]}"; do
    if [[ ! -r "$repo_root/$required_file" ]]; then
        echo "Required campaign file is missing or unreadable: $repo_root/$required_file" >&2
        echo 'No jobs were submitted.' >&2
        exit 2
    fi
done

campaign_logs="$repo_root/v2/logs/Subfield_12/$campaign_name"
stamp_file="$campaign_logs/validation_stamp.txt"
pre_status="$campaign_logs/pre_control_gate/pre_control_pH_discrepancy.status"
gate_results_glob="$repo_root/v2/SIM/Subfield_12/$campaign_name/pre_control_gate/*/TC_2D_run10_*.mat"
production_results_glob="$repo_root/v2/SIM/Subfield_12/$campaign_name/*/TC_2D_run*_*.mat"

# Files that score, validate, submit or document results but never run
# inside a trajectory. A rescore accepts a validation stamp that differs
# from the current one only in these files. The campaign fingerprint on
# Bouchet also hashes this controller (runtime/<campaign>/*.sh). Keep
# identical to allowed_relative_files in Accepted_Stamp_Digests.m;
# the configuration test compares them.
rescore_allowed_files=(
    'v2/T&C_Code/Tests/Kantola_N_Magnitude_Audit.m'
    'v2/T&C_Code/Tests/test_kantola_n_magnitude_audit.m'
    'v2/T&C_Code/Tests/test_bouchet_1d_joint_sweep_configuration.m'
    "runtime/$campaign_name/Summarize_Sweep.m"
    "runtime/$campaign_name/Validate_Sweep.m"
    "runtime/$campaign_name/Assert_Sweep_Provenance.m"
    "runtime/$campaign_name/Accepted_Stamp_Digests.m"
    "runtime/$campaign_name/README.md"
    "runtime/$campaign_name/submit_sweep.sh"
)
rescore_mode=0
lineage_changed=()

recovery_hint() {
    local finished
    finished="$(count_production_results)"
    if (( finished > 0 )); then
        echo "$finished finished production results exist under an earlier validation stamp." >&2
        echo 'To score them after a scoring-side change, run TC_SWEEP_STAGE=rescore TC_ACCEPT_PRIOR_STAMP=YES (README.md, Rescoring a complete sweep after a scoring fix).' >&2
        echo 'Recomputing any trajectory under the current stamp is a new campaign: move the campaign result and log directories aside first.' >&2
        echo 'No jobs were submitted.' >&2
        return 0
    fi
    echo 'Move the stale pre_control_gate result and logs aside (README.md, Recovering a stale one-year gate), rerun the exploratory stage, review the new one-year result, then release the full stage:' >&2
    printf '  stale=$(date -u +%%Y%%m%%dT%%H%%M%%SZ)\n' >&2
    printf '  mkdir -p %q %q\n' "$repo_root/v2/SIM/Subfield_12/stale_gates" "$repo_root/v2/logs/Subfield_12/stale_gates" >&2
    printf '  mv %q "%s/pre_control_gate_$stale"\n' "$repo_root/v2/SIM/Subfield_12/$campaign_name/pre_control_gate" "$repo_root/v2/SIM/Subfield_12/stale_gates" >&2
    printf '  mv %q "%s/pre_control_gate_$stale"\n' "$campaign_logs/pre_control_gate" "$repo_root/v2/logs/Subfield_12/stale_gates" >&2
    echo 'No jobs were submitted.' >&2
}

# Fingerprint_Digest hashes strtrim(fileread(stamp)); strip leading
# and trailing whitespace the same way before hashing.
read_stamp() {
    stamp_text="$(<"$stamp_file")"
    stamp_text="${stamp_text#"${stamp_text%%[![:space:]]*}"}"
    stamp_text="${stamp_text%"${stamp_text##*[![:space:]]}"}"
    current_digest="$(printf '%s' "$stamp_text" | sha256sum | awk '{print $1}')"
    archive_stamp
}

# Every stamp read by the controller is kept under its own digest, before
# any job can overwrite validation_stamp.txt, so that a rescore can later
# prove which files changed between the stamp a result carries and the
# stamp now on disk (Accepted_Stamp_Digests).
archive_stamp() {
    local archive_dir="$campaign_logs/validation_stamps"
    local archive="$archive_dir/$current_digest.txt" archived_text
    mkdir -p "$archive_dir"
    if [[ ! -f "$archive" ]]; then
        cp -p -- "$stamp_file" "$archive"
    fi
    archived_text="$(<"$archive")"
    archived_text="${archived_text#"${archived_text%%[![:space:]]*}"}"
    archived_text="${archived_text%"${archived_text##*[![:space:]]}"}"
    if [[ "$(printf '%s' "$archived_text" | sha256sum | awk '{print $1}')" != "$current_digest" ]]; then
        echo "The archived validation stamp $archive does not digest to its name. No jobs were submitted." >&2
        exit 6
    fi
}

# Digests of every archived stamp other than the current one, for the
# summary to judge (Accepted_Stamp_Digests accepts a prior stamp only
# when it differs from the current one in scoring-side files alone). The
# list is ':'-separated because sbatch splits its --export value on commas.
archived_prior_digests() {
    local archive name list=''
    for archive in "$campaign_logs"/validation_stamps/*.txt; do
        [[ -f "$archive" ]] || continue
        name="${archive##*/}"
        name="${name%.txt}"
        [[ "$name" =~ ^[0-9a-f]{64}$ ]] || continue
        [[ "$name" == "$current_digest" ]] && continue
        list+="${list:+:}$name"
    done
    printf '%s' "$list"
}

count_production_results() {
    local result count=0
    while IFS= read -r result; do
        [[ -n "$result" ]] && count=$((count + 1))
    done < <(compgen -G "$production_results_glob" || true)
    printf '%d' "$count"
}

# The stamp lists every fingerprinted file as <path>:<sha256>. Validation
# recomputes the stamp from the tree, so a listed file that changed or
# disappeared, or a file that appeared under one of the fingerprinted globs,
# means the next validation writes a different stamp and every result
# carrying the old digest becomes stale. Check all three here so that is
# refused before any job. The stamp's own path prefix is used for the glob
# scan so that a symlinked checkout path cannot make listed files look new.
# In rescore mode a divergence in a scoring-side file is recorded instead of
# refused; every other divergence is refused in every mode.
stamp_divergence() {
    local kind="$1" path="$2" rel allowed
    if (( rescore_mode )) && [[ -n "$stamp_root" && "$path" == "$stamp_root"/* ]]; then
        rel="${path#"$stamp_root"/}"
        for allowed in "${rescore_allowed_files[@]}"; do
            if [[ "$rel" == "$allowed" ]]; then
                lineage_changed+=("$rel")
                return 0
            fi
        done
    fi
    echo "A fingerprinted file $kind after the validation stamp was written: $path" >&2
    return 1
}

assert_stamp_current() {
    local part path expected actual glob file
    local -a parts
    local -A listed=()
    IFS='|' read -r -a parts <<< "$stamp_text"
    stamp_root=''
    for part in "${parts[@]}"; do
        case "$part" in NN=*|runtime=*) continue ;; esac
        path="${part%:*}"
        if [[ "$path" == */runtime/"$campaign_name"/README.md ]]; then
            stamp_root="${path%/runtime/"$campaign_name"/README.md}"
        fi
    done
    if [[ -z "$stamp_root" ]]; then
        echo 'The validation stamp does not list the campaign README.md.' >&2
        return 1
    fi
    for part in "${parts[@]}"; do
        case "$part" in NN=*|runtime=*) continue ;; esac
        path="${part%:*}"
        expected="${part##*:}"
        if [[ ! "$expected" =~ ^[0-9a-f]{64}$ ]]; then
            echo "The validation stamp has an unreadable entry: $part" >&2
            return 1
        fi
        listed["$path"]=1
        if [[ ! -f "$path" ]]; then
            stamp_divergence disappeared "$path" || return 1
            continue
        fi
        if [[ ! -r "$path" ]]; then
            echo "A fingerprinted file is not readable: $path" >&2
            return 1
        fi
        if ! actual="$(sha256sum -- "$path" | awk '{print $1}')" || \
                [[ ! "$actual" =~ ^[0-9a-f]{64}$ ]]; then
            echo "Cannot hash the fingerprinted file: $path" >&2
            return 1
        fi
        if [[ "$actual" != "$expected" ]]; then
            stamp_divergence changed "$path" || return 1
        fi
    done
    # Bouchet_Campaign_Validation_Fingerprint expands these globs at
    # validation time; a file that matches one of them but is not in the
    # stamp would enter the next stamp.
    for glob in \
            "$stamp_root/v2/T&C_Code/*.m" \
            "$stamp_root/v2/T&C_Code/Parameters/*.m" \
            "$stamp_root/v2/T&C_Code/Tests/*.m" \
            "$stamp_root/v2/T&C_Code/Tests/Fixtures/*.mat" \
            "$stamp_root/v2/run/*.m" \
            "$stamp_root/v2/audit/*.m" \
            "$stamp_root/v2/plotting/*.m" \
            "$stamp_root/runtime/$campaign_name/*.m" \
            "$stamp_root/runtime/$campaign_name/*.sbatch" \
            "$stamp_root/runtime/$campaign_name/*.sh"; do
        while IFS= read -r file; do
            [[ -n "$file" && -f "$file" ]] || continue
            if [[ -z "${listed[$file]+x}" ]]; then
                stamp_divergence appeared "$file" || return 1
            fi
        done < <(compgen -G "$glob" || true)
    done
    return 0
}

sweep_stage="${TC_SWEEP_STAGE:-exploratory}"
case "$sweep_stage" in
    exploratory)
        minimum_free_gib=10
        ;;
    full)
        minimum_free_gib=85
        if [[ "${TC_EXPLORATORY_REVIEWED:-}" != "YES" ]]; then
            echo 'Set TC_EXPLORATORY_REVIEWED=YES only after reviewing the one-year pH score.' >&2
            echo 'No jobs were submitted.' >&2
            exit 6
        fi
        if [[ ! -f "$pre_status" ]] || ! grep -Fxq 'state=complete' "$pre_status"; then
            echo 'The one-year control pH gate is not complete.' >&2
            if [[ -f "$pre_status" ]]; then
                sed -n 's/^error_identifier=/gate error identifier: /p' "$pre_status" >&2
            fi
            echo 'Run the exploratory stage: it recovers a current one-year result without recomputing it, or refuses with the move-aside instructions when the result is stale.' >&2
            echo 'No jobs were submitted.' >&2
            exit 6
        fi
        gate_result_count=0
        while IFS= read -r gate_result; do
            [[ -n "$gate_result" ]] && gate_result_count=$((gate_result_count + 1))
        done < <(compgen -G "$gate_results_glob" || true)
        if (( gate_result_count != 1 )); then
            printf 'The reviewed one-year control result should be exactly one MAT file under v2/SIM/Subfield_12/%s/pre_control_gate; found %d.\n' "$campaign_name" "$gate_result_count" >&2
            echo 'Move the pre_control_gate logs aside as well (README.md, Recovering a stale one-year gate) and run the exploratory stage. No jobs were submitted.' >&2
            exit 6
        fi
        # The reviewed one-year result records the SHA-256 digest of the
        # validation stamp it ran under. MATLAB rechecks this later, but a
        # mismatch is decided here so the refusal costs seconds, not a
        # validation job.
        if [[ ! -f "$stamp_file" ]]; then
            echo 'The campaign validation stamp is missing. Run the exploratory stage first. No jobs were submitted.' >&2
            exit 6
        fi
        read_stamp
        reviewed_digest="$(sed -n 's/^validation_fingerprint=//p' "$pre_status" | head -n 1)"
        if [[ ! "$reviewed_digest" =~ ^[0-9a-f]{64}$ || "$reviewed_digest" != "$current_digest" ]]; then
            echo 'The reviewed one-year control result was produced under a different validation stamp than the one now on disk.' >&2
            printf 'reviewed digest: %s\ncurrent digest:  %s\n' "$reviewed_digest" "$current_digest" >&2
            recovery_hint
            exit 6
        fi
        if ! assert_stamp_current; then
            echo 'The next validation would write a new stamp and the reviewed one-year result would no longer match it.' >&2
            recovery_hint
            exit 6
        fi
        ;;
    rescore)
        # Score a complete sweep again after a change to scoring-side code.
        # Nothing is recomputed: validation writes the new stamp, and the
        # summary accepts results carrying the previous stamp after proving
        # that the two stamps differ only in scoring-side files.
        minimum_free_gib=10
        if [[ "${TC_ACCEPT_PRIOR_STAMP:-}" != "YES" ]]; then
            echo 'Set TC_ACCEPT_PRIOR_STAMP=YES to score results produced under the previous validation stamp.' >&2
            echo 'No jobs were submitted.' >&2
            exit 6
        fi
        if [[ ! -f "$stamp_file" ]]; then
            echo 'The campaign validation stamp is missing; there is nothing to rescore. No jobs were submitted.' >&2
            exit 6
        fi
        read_stamp
        prior_digest="$current_digest"
        rescore_mode=1
        if ! assert_stamp_current; then
            echo 'A file outside the scoring-side allow-list changed after the validation stamp was written, so the finished trajectories cannot be scored by this tree.' >&2
            echo 'Restore that file to its stamped content, or treat the sweep as a new campaign. No jobs were submitted.' >&2
            exit 6
        fi
        production_result_count="$(count_production_results)"
        complete_status_count=0
        for status_path in "$campaign_logs"/*/run_*.status; do
            [[ -f "$status_path" ]] || continue
            if grep -Fxq 'state=complete' "$status_path"; then
                complete_status_count=$((complete_status_count + 1))
            fi
        done
        if (( production_result_count != 336 || complete_status_count != 336 )); then
            printf 'A rescore needs the complete sweep: found %d result files and %d complete status files, expected 336 of each.\n' \
                "$production_result_count" "$complete_status_count" >&2
            echo 'Finish the sweep with the full stage first. No jobs were submitted.' >&2
            exit 6
        fi
        ;;
    *)
        echo 'TC_SWEEP_STAGE must be exploratory, full or rescore. No jobs were submitted.' >&2
        exit 6
        ;;
esac

# The R2025a equivalence baseline is a dependency of the validation
# fingerprint and, like every MAT file, embeds its creation timestamp. A
# rebuilt baseline therefore changes the fingerprint even when its numbers are
# identical, and every result validated under the previous stamp reads as
# stale. So the baseline is kept unless it is missing or a rebuild is asked
# for explicitly, and the full stage never rebuilds it.
equivalence_baseline="$repo_root/v2/SIM/EQUIV_1D_TMP/TC_ERW_results_Kantola_crop_control.mat"
baseline_policy="${TC_EQUIVALENCE_BASELINE:-keep}"
case "$baseline_policy" in
    keep)
        if [[ -f "$equivalence_baseline" ]]; then
            rebuild_baseline=0
        elif [[ "$sweep_stage" == "exploratory" ]]; then
            rebuild_baseline=1
        else
            echo "The R2025a equivalence baseline is missing: $equivalence_baseline" >&2
            echo 'Run the exploratory stage first. No jobs were submitted.' >&2
            exit 6
        fi
        ;;
    rebuild)
        if [[ "$sweep_stage" != "exploratory" ]]; then
            echo 'Only the exploratory stage rebuilds the equivalence baseline: the rebuilt MAT file changes the validation fingerprint and orphans every result validated under the current stamp.' >&2
            echo 'Rebuild it with TC_SWEEP_STAGE=exploratory TC_EQUIVALENCE_BASELINE=rebuild, review the new one-year result, then release the full stage. No jobs were submitted.' >&2
            exit 6
        fi
        rebuild_baseline=1
        ;;
    *)
        echo 'TC_EQUIVALENCE_BASELINE must be keep or rebuild. No jobs were submitted.' >&2
        exit 6
        ;;
esac

if [[ "$sweep_stage" == "exploratory" ]] && compgen -G "$gate_results_glob" > /dev/null; then
    # A one-year result already exists. The control task recovers it only if
    # the stamp it carries is the stamp validation will write again: no
    # baseline rebuild, a complete status naming the current stamp digest and
    # no fingerprinted file changed since. Otherwise the task fails with
    # Bouchet:StaleResult after validation, so refuse here instead.
    gate_ok=1
    gate_state=''
    gate_result_count=0
    while IFS= read -r gate_result; do
        [[ -n "$gate_result" ]] && gate_result_count=$((gate_result_count + 1))
    done < <(compgen -G "$gate_results_glob" || true)
    if (( gate_result_count != 1 )); then
        printf 'Found %d one-year control result MAT files; the control task accepts exactly one.\n' "$gate_result_count" >&2
        gate_ok=0
    elif (( rebuild_baseline )); then
        echo 'A one-year control result exists and a baseline rebuild would change the validation stamp it carries.' >&2
        gate_ok=0
    elif [[ ! -f "$stamp_file" ]]; then
        echo 'A one-year control result exists but the campaign validation stamp is missing.' >&2
        gate_ok=0
    else
        read_stamp
        if ! assert_stamp_current; then
            echo 'A one-year control result exists but a fingerprinted file changed after its validation stamp was written.' >&2
            gate_ok=0
        else
            gate_state=''
            gate_error=''
            if [[ -f "$pre_status" ]]; then
                gate_state="$(sed -n 's/^state=//p' "$pre_status" | head -n 1)"
                gate_error="$(sed -n 's/^error_identifier=//p' "$pre_status" | head -n 1)"
            fi
            if [[ "$gate_state" == complete ]]; then
                reviewed_digest="$(sed -n 's/^validation_fingerprint=//p' "$pre_status" | head -n 1)"
                if [[ "$reviewed_digest" != "$current_digest" ]]; then
                    echo 'A one-year control result exists but its status names a different validation stamp than the one now on disk.' >&2
                    gate_ok=0
                fi
            else
                case "$gate_error" in
                    Bouchet_Joint_Sweep:WrongProvenance|Bouchet_Joint_Sweep:StaleValidation|Bouchet:StaleResult|Bouchet_Joint_Sweep:IncompletePreControl)
                        echo "The last summary rejected the existing one-year control result ($gate_error)." >&2
                        gate_ok=0
                        ;;
                    *)
                        # The status file is missing, still running or failed for
                        # a reason unrelated to the stamp, so the controller
                        # cannot tell whether the result carries the current
                        # digest; only the control task can, by reading the MAT.
                        # No fingerprinted file changed, so validation rewrites
                        # the same stamp and the task decides in seconds.
                        echo "Existing one-year control result: status is ${gate_state:-missing}${gate_error:+ ($gate_error)}, so its stamp digest cannot be confirmed here."
                        echo 'The control task recovers it if its recorded digest matches the current stamp, and otherwise fails with Bouchet:StaleResult; in that case follow README.md, Recovering a stale one-year gate.'
                        ;;
                esac
            fi
        fi
    fi
    if (( gate_ok == 0 )); then
        echo 'The control task would fail with Bouchet:StaleResult after validation.' >&2
        recovery_hint
        exit 6
    fi
    if [[ "$gate_state" == complete ]]; then
        echo 'Existing one-year control result is current; the control task will recover it without recomputation.'
    fi
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
if (( rebuild_baseline )); then
    if [[ -f "$equivalence_baseline" ]]; then
        echo 'Equivalence baseline: rebuild requested; the validation stamp will change and any earlier one-year result becomes stale.'
    else
        echo 'Equivalence baseline: missing, so a baseline job will build it.'
    fi
else
    echo "Equivalence baseline: reused $equivalence_baseline"
fi

: "${TC_BOUCHET_ACCOUNT:?Set TC_BOUCHET_ACCOUNT to the pi account reported by slurm_checkup.sh}"
if [[ "$sweep_stage" == "exploratory" ]]; then
    if [[ "${TC_SWEEP_LAYOUT:-independent}" != "independent" ]]; then
        echo 'The one-year control gate requires the independent layout.' >&2
        echo 'No jobs were submitted.' >&2
        exit 6
    fi
    sweep_layout=independent
    production_concurrency=1
elif [[ "$sweep_stage" == "rescore" ]]; then
    # The rescore submits no array; a layout left in the environment from
    # the full stage is irrelevant and ignored.
    sweep_layout=independent
    production_concurrency=1
else
    production_concurrency="${TC_SWEEP_MAX_CONCURRENT:-48}"
    if [[ ! "$production_concurrency" =~ ^[0-9]+$ ]] || \
            (( production_concurrency < 1 || production_concurrency > 48 )); then
        echo 'TC_SWEEP_MAX_CONCURRENT must be an integer from 1 through 48.' >&2
        echo 'No jobs were submitted.' >&2
        exit 6
    fi
    sweep_layout="${TC_SWEEP_LAYOUT:-independent}"
    case "$sweep_layout" in
        independent)
            pilot_array=336
            production_array="1-336%${production_concurrency}"
            production_tasks=1
            production_mem=8G
            ;;
        packed)
            if (( production_concurrency % 3 != 0 )); then
                echo 'Packed layout requires TC_SWEEP_MAX_CONCURRENT to be divisible by 3.' >&2
                echo 'No jobs were submitted.' >&2
                exit 6
            fi
            packed_concurrency=$((production_concurrency / 3))
            pilot_array=112
            production_array="1-112%${packed_concurrency}"
            production_tasks=3
            production_mem=24G
            ;;
        *)
            echo 'TC_SWEEP_LAYOUT must be independent or packed.' >&2
            echo 'No jobs were submitted.' >&2
            exit 6
            ;;
    esac
fi

unset TC_REBUILD_EQUIVALENCE TC_VALIDATE_ONLY TC_SUMMARIZE_ONLY \
    TC_SWEEP_CASE_ID TC_PACKED_SMOKE TC_SWEEP_SCOPE \
    TC_SWEEP_SUMMARY_SCOPE TC_SWEEP_CONTROLLER_TOKEN TC_EQUIVALENCE_BASELINE \
    TC_ACCEPTED_RESULT_DIGESTS TC_ACCEPT_PRIOR_STAMP
for submit_env in "${!SBATCH_@}"; do
    unset "$submit_env"
done
unset SLURM_CLUSTERS
export MATLAB_MODULE="${MATLAB_MODULE:-MATLAB/2025a}"
if ! module spider "$MATLAB_MODULE" >/dev/null 2>&1; then
    echo "Bouchet does not currently expose $MATLAB_MODULE. No jobs were submitted." >&2
    exit 3
fi

controller_lock_root="$repo_root/v2/logs/Subfield_12/$campaign_name"
mkdir -p "$controller_lock_root"
controller_lock="$controller_lock_root/.submission_controller.lock"
if ! mkdir "$controller_lock" 2>/dev/null; then
    echo "Submission lock already exists: $controller_lock" >&2
    echo 'No jobs were submitted. Check for an active controller or stale lock.' >&2
    exit 5
fi
release_controller_lock() {
    if ! rmdir "$controller_lock" 2>/dev/null; then
        echo "Warning: could not remove submission lock $controller_lock" >&2
    fi
}
trap release_controller_lock EXIT

job_names='SweepBase,SweepCheck,SweepExplore,SweepExploreSum,SweepPilot,Sweep1D,SweepSummary'
if ! active_jobs="$(squeue -h -u "$USER" --name="$job_names" -o '%A|%j|%T')"; then
    echo 'Cannot check for an active joint-sweep graph. No jobs were submitted.' >&2
    exit 5
fi
if [[ -n "$active_jobs" ]]; then
    echo 'An active joint-sweep graph already exists. No jobs were submitted.' >&2
    printf '%s\n' "$active_jobs" >&2
    exit 5
fi

# No graph job is active, so any per-case run lock left under the campaign
# logs belongs to a job that was killed before MATLAB could release it. A
# resubmitted element would fail with Bouchet:RunLocked, so refuse until the
# stale locks are removed by hand.
stale_locks="$(find "$controller_lock_root" -type d -name 'run_*.lock' 2>/dev/null || true)"
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
            echo "Manual cancellation required: scancel ${submitted_jobs[*]}" >&2
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
    if ! submitted="$(sbatch --parsable --account="$TC_BOUCHET_ACCOUNT" "$@")"; then
        echo 'sbatch submission failed.' >&2
        return 1
    fi
    submitted_id="${submitted%%;*}"
    if [[ ! "$submitted_id" =~ ^[0-9]+$ ]]; then
        echo "sbatch returned an invalid job id: $submitted" >&2
        return 1
    fi
    printf -v "$result_var" '%s' "$submitted_id"
    submitted_jobs+=("$submitted_id")
}

graph_jobs=()
if [[ "$sweep_stage" == "rescore" ]]; then
    # Validation may rewrite the stamp before the summary runs, so the
    # summary is offered every archived stamp (the one on disk now included
    # when it changes); it accepts each one only after proving that it
    # differs from the stamp it validates against in scoring-side files.
    accepted_digests="$(archived_prior_digests)"
    if (( ${#lineage_changed[@]} > 0 )); then
        accepted_digests="${accepted_digests:+$accepted_digests:}$prior_digest"
    fi
    if (( ${#lineage_changed[@]} > 0 )); then
        echo "Rescore: scoring-side files changed since stamp $prior_digest:"
        printf '  %s\n' "${lineage_changed[@]}"
        echo "Validation will write a new stamp; the summary is offered the archived stamps $accepted_digests."
        submit_job validation_job \
            --array=1 --partition=day --ntasks=1 \
            --cpus-per-task=6 --mem=40G --time=08:00:00 \
            --job-name=SweepCheck \
            --export=ALL,TC_SWEEP_CONTROLLER_TOKEN=CEC13_STAGED_V4,TC_REBUILD_EQUIVALENCE=0,TC_VALIDATE_ONLY=1,TC_SUMMARIZE_ONLY=0,TC_SWEEP_SCOPE=full,TC_SWEEP_SUMMARY_SCOPE=full "$batch_file"
        summary_dependency=(--dependency="afterok:${validation_job}" --kill-on-invalid-dep=yes)
        graph_jobs+=("$validation_job")
    else
        echo "Rescore: no fingerprinted file changed since stamp $prior_digest; the summary runs under the current stamp."
        if [[ -n "$accepted_digests" ]]; then
            echo "The summary is offered the archived prior stamps $accepted_digests for results that carry them."
        fi
        validation_job='not needed'
        summary_dependency=()
    fi
    submit_job summary_job \
        --array=1 ${summary_dependency[@]+"${summary_dependency[@]}"} \
        --partition=day --ntasks=1 \
        --cpus-per-task=1 --mem=12G --time=02:00:00 \
        --job-name=SweepSummary \
        --export=ALL,TC_SWEEP_CONTROLLER_TOKEN=CEC13_STAGED_V4,TC_REBUILD_EQUIVALENCE=0,TC_VALIDATE_ONLY=0,TC_SUMMARIZE_ONLY=1,TC_SWEEP_SCOPE=full,TC_SWEEP_SUMMARY_SCOPE=full,TC_ACCEPTED_RESULT_DIGESTS=$accepted_digests "$batch_file"
    graph_jobs+=("$summary_job")
    graph_complete=1
    trap - ERR HUP INT TERM
    printf 'Validation: %s\nSummary: %s\n' "$validation_job" "$summary_job"
    for job_id in "${graph_jobs[@]}"; do
        if ! scontrol show jobid="$job_id" -dd | \
                grep -E 'JobId=|JobState=|Dependency=|Constraint=|Partition=|SubmitLine='; then
            printf 'Warning: could not inspect submitted job %s. The job remains submitted.\n' \
                "$job_id" >&2
        fi
    done
    printf '\nRESCORE: no trajectory is recomputed; the summary re-scores the 336 finished results.\n'
    exit 0
fi

if (( rebuild_baseline )); then
    submit_job baseline_job \
        --array=1 --partition=day --ntasks=1 --cpus-per-task=1 \
        --mem=16G --time=02:00:00 --job-name=SweepBase \
        --export=ALL,TC_SWEEP_CONTROLLER_TOKEN=CEC13_STAGED_V4,TC_REBUILD_EQUIVALENCE=1,TC_VALIDATE_ONLY=0,TC_SUMMARIZE_ONLY=0,TC_SWEEP_SCOPE=full,TC_SWEEP_SUMMARY_SCOPE=full "$batch_file"
    validation_dependency=(--dependency="afterok:${baseline_job}" --kill-on-invalid-dep=yes)
    graph_jobs+=("$baseline_job")
else
    baseline_job='reused'
    validation_dependency=()
fi

submit_job validation_job \
    --array=1 ${validation_dependency[@]+"${validation_dependency[@]}"} \
    --partition=day --ntasks=1 \
    --cpus-per-task=6 --mem=40G --time=08:00:00 \
    --job-name=SweepCheck \
    --export=ALL,TC_SWEEP_CONTROLLER_TOKEN=CEC13_STAGED_V4,TC_REBUILD_EQUIVALENCE=0,TC_VALIDATE_ONLY=1,TC_SUMMARIZE_ONLY=0,TC_SWEEP_SCOPE=full,TC_SWEEP_SUMMARY_SCOPE=full "$batch_file"

if [[ "$sweep_stage" == "exploratory" ]]; then
    submit_job exploratory_job \
        --array=1 --dependency="afterok:${validation_job}" \
        --kill-on-invalid-dep=yes --partition=day --ntasks=1 \
        --cpus-per-task=1 --mem=8G --time=02:00:00 \
        --job-name=SweepExplore \
        --export=ALL,TC_SWEEP_CONTROLLER_TOKEN=CEC13_STAGED_V4,TC_REBUILD_EQUIVALENCE=0,TC_VALIDATE_ONLY=0,TC_SUMMARIZE_ONLY=0,TC_SWEEP_LAYOUT=independent,TC_SWEEP_SCOPE=pre_control,TC_SWEEP_SUMMARY_SCOPE=full,TC_PACKED_SMOKE=0 "$batch_file"

    submit_job exploratory_summary_job \
        --array=1 --dependency="afterany:${exploratory_job}" \
        --kill-on-invalid-dep=yes --partition=day --ntasks=1 \
        --cpus-per-task=1 --mem=4G --time=00:30:00 \
        --job-name=SweepExploreSum \
        --export=ALL,TC_SWEEP_CONTROLLER_TOKEN=CEC13_STAGED_V4,TC_REBUILD_EQUIVALENCE=0,TC_VALIDATE_ONLY=0,TC_SUMMARIZE_ONLY=1,TC_SWEEP_SCOPE=full,TC_SWEEP_SUMMARY_SCOPE=pre_control,TC_PRE_CONTROL_RECHECK=0 "$batch_file"

    graph_complete=1
    trap - ERR HUP INT TERM
    printf 'Baseline: %s\nValidation: %s\nExploratory control: %s\nExploratory summary: %s\n' \
        "$baseline_job" "$validation_job" "$exploratory_job" \
        "$exploratory_summary_job"
    graph_jobs+=("$validation_job" "$exploratory_job" "$exploratory_summary_job")
    for job_id in "${graph_jobs[@]}"; do
        if ! scontrol show jobid="$job_id" -dd | \
                grep -E 'JobId=|JobState=|Dependency=|Constraint=|Partition=|SubmitLine='; then
            printf 'Warning: could not inspect submitted job %s. The job remains submitted.\n' \
                "$job_id" >&2
        fi
    done
    printf '\nSTOP: no pilot or production sweep was submitted.\n'
    exit 0
fi

submit_job exploratory_review_job \
    --array=1 --dependency="afterok:${validation_job}" \
    --kill-on-invalid-dep=yes --partition=day --ntasks=1 \
    --cpus-per-task=1 --mem=4G --time=00:30:00 \
    --job-name=SweepExploreSum \
    --export=ALL,TC_SWEEP_CONTROLLER_TOKEN=CEC13_STAGED_V4,TC_REBUILD_EQUIVALENCE=0,TC_VALIDATE_ONLY=0,TC_SUMMARIZE_ONLY=1,TC_SWEEP_SCOPE=full,TC_SWEEP_SUMMARY_SCOPE=pre_control,TC_PRE_CONTROL_RECHECK=1 "$batch_file"

submit_job pilot_job \
    --array="$pilot_array" --dependency="afterok:${exploratory_review_job}" \
    --kill-on-invalid-dep=yes --partition=day --ntasks="$production_tasks" \
    --cpus-per-task=1 --mem="$production_mem" --time=05:00:00 \
    --job-name=SweepPilot \
    --export=ALL,TC_SWEEP_CONTROLLER_TOKEN=CEC13_STAGED_V4,TC_REBUILD_EQUIVALENCE=0,TC_VALIDATE_ONLY=0,TC_SUMMARIZE_ONLY=0,TC_SWEEP_LAYOUT="$sweep_layout",TC_SWEEP_SCOPE=full,TC_SWEEP_SUMMARY_SCOPE=full,TC_PACKED_SMOKE=1 "$batch_file"

submit_job production_job \
    --array="$production_array" --dependency="afterok:${pilot_job}" \
    --kill-on-invalid-dep=yes --partition=day --ntasks="$production_tasks" \
    --cpus-per-task=1 --mem="$production_mem" --time=05:00:00 \
    --job-name=Sweep1D \
    --export=ALL,TC_SWEEP_CONTROLLER_TOKEN=CEC13_STAGED_V4,TC_REBUILD_EQUIVALENCE=0,TC_VALIDATE_ONLY=0,TC_SUMMARIZE_ONLY=0,TC_SWEEP_LAYOUT="$sweep_layout",TC_SWEEP_SCOPE=full,TC_SWEEP_SUMMARY_SCOPE=full,TC_PACKED_SMOKE=0 "$batch_file"

submit_job summary_job \
    --array=1 --dependency="afterany:${production_job}" \
    --kill-on-invalid-dep=yes --partition=day --ntasks=1 \
    --cpus-per-task=1 --mem=12G --time=02:00:00 \
    --job-name=SweepSummary \
    --export=ALL,TC_SWEEP_CONTROLLER_TOKEN=CEC13_STAGED_V4,TC_REBUILD_EQUIVALENCE=0,TC_VALIDATE_ONLY=0,TC_SUMMARIZE_ONLY=1,TC_SWEEP_SCOPE=full,TC_SWEEP_SUMMARY_SCOPE=full "$batch_file"

graph_complete=1
trap - ERR HUP INT TERM
printf 'Baseline: %s\nValidation: %s\nExploratory recheck: %s\nPilot: %s\nProduction: %s\nSummary: %s\n' \
    "$baseline_job" "$validation_job" "$exploratory_review_job" \
    "$pilot_job" "$production_job" "$summary_job"
printf 'Layout: %s, maximum concurrent MATLAB cases: %s\n' \
    "$sweep_layout" "$production_concurrency"
graph_jobs+=("$validation_job" "$exploratory_review_job" \
    "$pilot_job" "$production_job" "$summary_job")
for job_id in "${graph_jobs[@]}"; do
    if ! scontrol show jobid="$job_id" -dd | \
            grep -E 'JobId=|JobState=|Dependency=|Constraint=|Partition=|SubmitLine='; then
        printf 'Warning: could not inspect submitted job %s. The job remains submitted.\n' \
            "$job_id" >&2
    fi
done

printf '\nThe summary waits for every production element, including failures.\n'
printf 'This controller does not submit the final all-12 campaign.\n'
