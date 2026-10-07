# Sourced by the ensemble batch files. A day-partition job stops at 24
# hours, so a run that needs longer requeues itself when Slurm signals that
# the limit is near (--signal=B:USR1@1800, one checkpoint interval ahead)
# and continues from its last checkpoint under the same job id, which also
# keeps its 2-D run lock. Five restarts cover a 35-year 2-D run with room
# to spare; a run that needs more is not progressing and is left to time
# out. The batch files open their logs in append mode, so each piece adds
# to the same log instead of replacing it.

max_restarts=5

requeue_self() {
    signalled=1
    local restarts="${SLURM_RESTART_COUNT:-0}"
    # An array task is named <array>_<task>: the first task's own job id is
    # also the array's, and requeueing that would restart every task.
    local target="${SLURM_ARRAY_JOB_ID:+${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}}"
    if (( restarts < max_restarts )); then
        echo "Time limit near: requeueing ${target:-$SLURM_JOB_ID} (restart $(( restarts + 1 )) of $max_restarts)."
        # A busy controller can refuse a call for a moment; the signal comes
        # half an hour early, so there is time to try again. If Slurm still
        # refuses, the run keeps checkpointing until the limit, and a
        # resubmission continues it from there.
        local try
        for try in 1 2 3; do
            scontrol requeue "${target:-$SLURM_JOB_ID}" && return 0
            sleep 20
        done
        echo "The requeue of ${target:-$SLURM_JOB_ID} failed; it will stop at its time limit." >&2
    else
        echo "Time limit near after $max_restarts restarts; ${target:-$SLURM_JOB_ID} is not requeued." >&2
    fi
}

# The step runs in the background because bash runs a trap only between
# commands: a foreground step would hold the signal until it had ended.
run_requeueing() {
    trap requeue_self USR1
    "$@" &
    local step=$! status
    while true; do
        signalled=0
        wait "$step" && status=0 || status=$?
        # A trapped signal ends the wait early; waiting again returns the
        # step's own status, at once if it ended while the trap ran.
        if (( signalled && status > 128 )); then continue; fi
        return "$status"
    done
}
