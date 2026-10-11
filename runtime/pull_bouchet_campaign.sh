#!/usr/bin/env bash
#
# Pull a finished Bouchet twelve-run campaign into the OneDrive archive.
#
# Runs on the workstation in Git Bash, not on Bouchet. It refuses to archive a
# campaign that is not complete, so a half-finished transfer never becomes the
# comparison set. One ssh control connection is opened, so Duo is answered once.
#
#   bash runtime/pull_bouchet_campaign.sh \
#       --netid gj299 \
#       --remote-root '/home/gj299/project/tc-smewn-v09-p-return' \
#       --label p_return_nresponse_v1
#
# Options
#   --netid            Yale netid used for ssh                       (required)
#   --remote-root      deployment root on Bouchet                    (required)
#   --label            short archive name, no spaces                 (required)
#   --host             login host              [bouchet.ycrc.yale.edu]
#   --archive-base     where the archive folder is created  [v2/SIM/Subfield_12]
#   --stage            staging folder outside OneDrive     [~/bouchet_staging]
#   --logs-only        skip the result MATs, take logs and slurm output only
#   --yes              do not ask before transferring

set -Eeuo pipefail

host='bouchet.ycrc.yale.edu'
netid=''
remote_root=''
label=''
repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
archive_base="$repo_root/v2/SIM/Subfield_12"
stage_base="$HOME/bouchet_staging"
logs_only=0
assume_yes=0

campaign_1d='Kantola_5yr_1D_Runs10to12_ScheduleFix_R2025a_Bouchet'
campaign_2d='Kantola_5yr_2D_Runs01to09_DenitPH7p8_Alpha1p3_Bouchet'

while (( $# > 0 )); do
    case "$1" in
        --netid)        netid="$2";        shift 2 ;;
        --remote-root)  remote_root="$2";  shift 2 ;;
        --label)        label="$2";        shift 2 ;;
        --host)         host="$2";         shift 2 ;;
        --archive-base) archive_base="$2"; shift 2 ;;
        --stage)        stage_base="$2";   shift 2 ;;
        --campaign-1d)  campaign_1d="$2";  shift 2 ;;
        --campaign-2d)  campaign_2d="$2";  shift 2 ;;
        --logs-only)    logs_only=1;       shift   ;;
        --yes|-y)       assume_yes=1;      shift   ;;
        *) echo "Unknown option: $1" >&2; exit 2 ;;
    esac
done

for required in netid remote_root label; do
    if [[ -z "${!required}" ]]; then
        echo "--${required//_/-} is required." >&2
        exit 2
    fi
done
if [[ "$label" == *[[:space:]]* || "$label" == */* ]]; then
    echo 'The label must be a single path-safe word.' >&2
    exit 2
fi

stamp="$(date -u +%Y%m%dT%H%M%SZ)"
archive_dir="$archive_base/archive_${stamp}_${label}"
stage_dir="$stage_base/$stamp"
if [[ -e "$archive_dir" ]]; then
    echo "Archive path already exists: $archive_dir" >&2
    exit 4
fi

control_path="$HOME/.ssh/cm-bouchet-$$"
cleanup() {
    ssh -O exit -o "ControlPath=$control_path" "$netid@$host" >/dev/null 2>&1 || true
}
trap cleanup EXIT

echo "Opening one ssh connection to $host. Answer Duo once."
ssh -fN -o ControlMaster=yes -o "ControlPath=$control_path" -o ControlPersist=30m "$netid@$host"

remote() {
    ssh -o "ControlPath=$control_path" "$netid@$host" "$@"
}

## The campaign is complete only when every expected 1826-day MAT exists and
## every status file says so. Anything else stops here, before a partial set is
## archived and later mistaken for the comparison basis.
echo 'Checking the campaign on Bouchet.'
read -r -d '' precheck_script <<REMOTE_PRECHECK || true
set -uo pipefail
root='${remote_root}'
sim1d="\$root/v2/SIM/Subfield_12/${campaign_1d}"
sim2d="\$root/v2/SIM/Subfield_12/${campaign_2d}"
log1d="\$root/v2/logs/Subfield_12/${campaign_1d}"
log2d="\$root/v2/logs/Subfield_12/${campaign_2d}"
missing=''
for run_id in 10 11 12; do
    ls "\$sim1d"/TC_2D_run\${run_id}_*.mat >/dev/null 2>&1 || missing="\$missing \$run_id"
done
for run_id in 01 02 03 04 05 06 07 08 09; do
    ls "\$sim2d"/TC_2D_run\${run_id}_*.mat >/dev/null 2>&1 || missing="\$missing \$run_id"
done
incomplete=''
for status_file in "\$log1d"/run_*.status "\$log2d"/run_*.status; do
    [ -e "\$status_file" ] || continue
    grep -q 'state=complete' "\$status_file" || incomplete="\$incomplete \$(basename "\$status_file")"
done
echo "MISSING_RUNS:\$missing"
echo "INCOMPLETE_STATUS:\$incomplete"
echo "LIVE_JOBS:\$(squeue -u "\$USER" -h -o '%i %j %T' | grep Kantola | tr '\n' ';')"
echo "BYTES:\$(du -sb --exclude=checkpoints --exclude=locks "\$sim1d" "\$sim2d" "\$log1d" "\$log2d" 2>/dev/null | awk '{s+=\$1} END {print s+0}')"
REMOTE_PRECHECK

precheck="$(remote bash -s <<<"$precheck_script")"
printf '%s\n' "$precheck"

missing_runs="$(sed -n 's/^MISSING_RUNS://p' <<<"$precheck")"
incomplete_status="$(sed -n 's/^INCOMPLETE_STATUS://p' <<<"$precheck")"
live_jobs="$(sed -n 's/^LIVE_JOBS://p' <<<"$precheck")"
byte_total="$(sed -n 's/^BYTES://p' <<<"$precheck")"

if [[ -n "${live_jobs//[[:space:];]/}" ]]; then
    echo 'A Kantola job is still queued or running. Archive after it finishes.' >&2
    exit 5
fi
if [[ -n "${missing_runs// /}" && $logs_only -eq 0 ]]; then
    echo "Missing final results for runs:$missing_runs" >&2
    echo 'Nothing was transferred.' >&2
    exit 5
fi
if [[ -n "${incomplete_status// /}" ]]; then
    echo "Status files that do not say complete:$incomplete_status" >&2
    echo 'Nothing was transferred.' >&2
    exit 5
fi

human_size="$(awk -v b="${byte_total:-0}" 'BEGIN {printf "%.1f GB", b/1073741824}')"
echo "Campaign size on Bouchet, checkpoints excluded: $human_size"
if (( assume_yes == 0 )); then
    read -r -p "Transfer into $archive_dir ? [y/N] " reply
    if [[ "$reply" != [yY] ]]; then
        echo 'Stopped. Nothing was transferred.'
        exit 0
    fi
fi

## The checksums are taken on Bouchet before the bytes move, then verified here
## after they land, so a truncated stream is caught rather than archived.
mkdir -p "$stage_dir/slurm"
echo 'Recording checksums on Bouchet.'
remote "cd '$remote_root' && find \
    'v2/SIM/Subfield_12/$campaign_1d' 'v2/SIM/Subfield_12/$campaign_2d' \
    'v2/logs/Subfield_12/$campaign_1d' 'v2/logs/Subfield_12/$campaign_2d' \
    -type f -not -path '*/checkpoints/*' -not -path '*/locks/*' -print0 \
    | sort -z | xargs -0 sha256sum" > "$stage_dir/MANIFEST.sha256"

echo 'Recording the Slurm accounting record.'
remote "sacct -u '$netid' --starttime now-14days \
    --name=KantolaBase,Kantola1Dcheck,Kantola2Dcheck,Kantola1Dfix,Kantola2D \
    --format=JobID%20,JobName%18,State,ExitCode,Submit,Start,End,Elapsed,MaxRSS,NodeList%20" \
    > "$stage_dir/slurm/sacct.txt" || true

tar_excludes=(--exclude='checkpoints' --exclude='locks')
if (( logs_only == 1 )); then
    tar_excludes+=(--exclude='*.mat')
fi

echo 'Transferring. Checkpoints and locks are left on Bouchet.'
remote "cd '$remote_root' && tar -cf - --exclude=checkpoints --exclude=locks \
    'v2/SIM/Subfield_12/$campaign_1d' 'v2/SIM/Subfield_12/$campaign_2d' \
    'v2/logs/Subfield_12/$campaign_1d' 'v2/logs/Subfield_12/$campaign_2d'" \
    | tar -xf - -C "$stage_dir" "${tar_excludes[@]}"

remote "cd '$remote_root' && tar -cf - \
    \$(find 'runtime/$campaign_1d' 'runtime/$campaign_2d' \
        \\( -name 'slurm_*.out' -o -name 'slurm_*.err' \\) -type f)" \
    | tar -xf - -C "$stage_dir" || echo 'No Slurm startup logs were found.'

if (( logs_only == 0 )); then
    echo 'Verifying checksums locally.'
    ( cd "$stage_dir" && sha256sum -c --quiet MANIFEST.sha256 )
    echo 'All transferred files match their Bouchet checksums.'
fi

## Plot_1D_Then_All_12 reads one ResultDirectory, so the twelve MATs are
## flattened into results/. The two campaign log trees stay separate.
echo 'Laying out the archive.'
mkdir -p "$stage_dir/results" "$stage_dir/logs" "$stage_dir/figures"
find "$stage_dir/v2/SIM/Subfield_12" -name 'TC_2D_run*.mat' -exec mv -n {} "$stage_dir/results/" ';'

## Run_One_Subfield names a result from every option that was set, so the
## Bouchet files carry __grid2x2, __alpha1p3 and __from20160101 as well as
## __1826d. Plot_1D_Then_All_12 builds only "<treatment>__<layout>__1826d", so
## the extra tags are trimmed here, after the checksums were verified against
## the delivered names. RENAMED.txt keeps the mapping.
: > "$stage_dir/RENAMED.txt"
for delivered in "$stage_dir/results"/TC_2D_run*.mat; do
    [[ -e "$delivered" ]] || continue
    delivered_name="$(basename "$delivered")"
    plot_name="$(sed -E 's/__grid[0-9]+x[0-9]+//; s/__alpha[0-9p]+//; s/__logkh[0-9p_]+//; s/__from[0-9]{8}//' <<<"$delivered_name")"
    if [[ "$plot_name" != "$delivered_name" ]]; then
        mv -n "$delivered" "$stage_dir/results/$plot_name"
        printf '%s -> %s\n' "$delivered_name" "$plot_name" >> "$stage_dir/RENAMED.txt"
    fi
done
mv "$stage_dir/v2/logs/Subfield_12/$campaign_1d" "$stage_dir/logs/1D_runs10to12"
mv "$stage_dir/v2/logs/Subfield_12/$campaign_2d" "$stage_dir/logs/2D_runs01to09"
if [[ -d "$stage_dir/runtime" ]]; then
    find "$stage_dir/runtime" -type f -name 'slurm_*' -exec mv -n {} "$stage_dir/slurm/" ';'
fi
rm -rf "$stage_dir/v2" "$stage_dir/runtime"

{
    printf 'Display name: archive_%s_%s\n' "$stamp" "$label"
    printf 'Campaign tag 1-D: %s\n' "$campaign_1d"
    printf 'Campaign tag 2-D: %s\n' "$campaign_2d"
    if (( logs_only == 1 )); then
        printf 'Status: LOGS ONLY, result MATs left on Bouchet\n'
    else
        printf 'Status: COMPLETE\n'
    fi
    printf 'Runs: 1-3 whole_field; 4-6 margin_crop; 7-9 margin_prairie; 10-12 vertical_1d\n'
    printf 'Treatments: control, basalt, calcite in each of the four layouts\n'
    printf 'Duration: 1826 days from 2016-01-01\n'
    printf 'Appelo alpha: 1.3\n'
    printf 'Execution: Bouchet MATLAB/2025a, %s:%s\n' "$host" "$remote_root"
    printf 'Imported: %s\n' "$(date -u +%Y-%m-%d)"
    printf 'Checksums: MANIFEST.sha256, verified after transfer\n'
    printf 'Slurm accounting: slurm/sacct.txt\n'
    printf '\nResult files\n'
    ( cd "$stage_dir/results" && ls -1 ) || true
} > "$stage_dir/RUN_INFO.txt"

## OneDrive only sees the archive once it is complete and verified, so a
## partial folder is never uploaded and never syncs to the other machines.
mkdir -p "$archive_base"
mv "$stage_dir" "$archive_dir"
rmdir "$stage_base" 2>/dev/null || true

echo
echo "Archived to: $archive_dir"
echo 'Compare the twelve runs with:'
echo "  matlab -batch \"cd('v2/plotting'); NDays=1826; ResultDirectory='$archive_dir/results'; FigureDirectory='$archive_dir/figures'; Plot_1D_Then_All_12\""
