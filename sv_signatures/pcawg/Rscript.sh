#!/usr/bin/env bash
set -euo pipefail
pipeline_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
if [[ -n ${SV_R_MODULE:-} ]]; then module load "$SV_R_MODULE"; fi
library_file="$pipeline_root/environments/clustersv-r43/library-path.txt"
[[ -s $library_file ]] || { echo 'Run Rscript sv_signatures/pcawg/setup_env.R PROJECT_ROOT first' >&2; exit 2; }
export R_LIBS_USER="$(< "$library_file")"
exec Rscript --vanilla "$@"
