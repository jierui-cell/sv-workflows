#!/usr/bin/env bash
# Submit with site-owned sbatch flags; no lab defaults are embedded.
set -Eeuo pipefail
: "${SLURM_JOB_ID:?Submit this driver with sbatch}"
root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$root"
exec python3 scripts/workflow.py "$@" --profile slurm
