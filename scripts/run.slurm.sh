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
if [[ -n ${SV_PYTHON:-} ]]; then
    launcher_python=$SV_PYTHON
elif [[ -x "$root/.venv/bin/python" ]]; then
    launcher_python="$root/.venv/bin/python"
else
    launcher_python=python3
fi
if [[ -d "$root/.venv/bin" ]]; then
    export PATH="$root/.venv/bin:$PATH"
fi
# workflow.py loads YAML modules before invoking Nextflow.
exec "$launcher_python" scripts/workflow.py "$@" --profile slurm
