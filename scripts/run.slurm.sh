#!/usr/bin/env bash
# Submit with site-owned sbatch flags; no lab defaults are embedded.
set -Eeuo pipefail
: "${SLURM_JOB_ID:?Submit this driver with sbatch}"
root=${SV_WORKFLOW_ROOT:-${SLURM_SUBMIT_DIR:?Submit from the repository root or set SV_WORKFLOW_ROOT}}
root=$(cd "$root" && pwd -P)
[[ -f "$root/scripts/workflow.py" && -f "$root/main.nf" ]] || {
    echo 'Repository not found: submit from its root, or export SV_WORKFLOW_ROOT=/absolute/repository/path' >&2
    exit 2
}
cd "$root"
exec "${SV_PYTHON:-python3}" scripts/workflow.py "$@" --profile slurm
