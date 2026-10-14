#!/bin/bash
# Invoke the pinned BRASS image with its validated resources and mandatory overlay.
# Usage (inside a Slurm job): bash scripts/brass_with_bundle.sh BUNDLE [brass.pl options]
set -Eeuo pipefail
: "${SV_JOB_ID:?Launch through scripts/workflow.py}"
root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
: "${1:?Provide a validated BRASS bundle directory}"
bundle=$(realpath "$1")
shift
test -f "$bundle/COMPATIBILITY_VALIDATED"
test -s "$bundle/native-validation.json"
test -s "$bundle/get_rg_cns.R"
# Resource arguments are owned by this adapter, not silently overridable.
for option in "$@"; do
    option=${option%%=*}
    if [[ $option == --* ]]; then option=${option#-}; fi
    case "$option" in
        -g|-genome|-d|-depth|-gc|-g_cache|-vi|-viral|-mi|-microbe|-b|-gcbins|-cb|-cytoband|-ct|-centtel|-f|-filter)
            echo "Resource override $option is not supported by this validated bundle adapter" >&2
            exit 2 ;;
    esac
done
python3 - "$bundle/native-validation.json" <<'PY'
import json
import sys
with open(sys.argv[1]) as handle:
    report = json.load(handle)
if report.get('status') != 'passed':
    raise SystemExit('BRASS native resource validation did not pass')
PY
overlay=/opt/wtsi-cgp/lib/perl5/auto/share/module/Sanger-CGP-Brass-Implement/Rscripts/get_rg_cns.R
# Keep native code in per-job provenance; never modify the pinned image.
patch_dir="${BRASS_RUN_DIR:?}/sample-name-adapter"
mkdir -p "$patch_dir"
native_module=/opt/wtsi-cgp/lib/perl5/Sanger/CGP/Brass/Implement.pm
if [[ -f $patch_dir/adapter.sha256 ]]; then
    (cd "$patch_dir"; sha256sum -c adapter.sha256)
    cmp "$root/scripts/patch_brass_sample_name.py" "$patch_dir/patch_brass_sample_name.py"
else
    # Native source is read-only. Never overwrite it on subsequent stage calls.
    if [[ ! -e $patch_dir/Implement.original.pm ]]; then
        singularity exec --cleanenv --bind "$root:$root" \
            "${BRASS_IMAGE:?}" cp "$native_module" "$patch_dir/Implement.original.pm"
    fi
    native_sha=$(singularity exec --cleanenv "${BRASS_IMAGE:?}" sha256sum "$native_module")
    original_sha=$(sha256sum "$patch_dir/Implement.original.pm")
    [[ ${native_sha%% *} == "${original_sha%% *}" ]]
    python3 "$root/scripts/patch_brass_sample_name.py" "$patch_dir/Implement.original.pm" "$patch_dir/Implement.pm"
    cp "$root/scripts/patch_brass_sample_name.py" "$patch_dir/"
    (cd "$patch_dir"; sha256sum Implement.original.pm Implement.pm Implement.patch.json patch_brass_sample_name.py > adapter.sha256)
fi
implementation="$patch_dir/Implement.pm"
if [[ ${BRASS_STATISTICS_STATE:-} == MIXED ]]; then
    python3 "$root/scripts/patch_brass_mixed_statistics.py" "$implementation" "$patch_dir/Implement.mixed.pm"
    cp "$root/scripts/patch_brass_mixed_statistics.py" "$patch_dir/"
    sha256sum "$patch_dir/Implement.mixed.pm" > "$patch_dir/mixed.sha256"
    implementation="$patch_dir/Implement.mixed.pm"
fi
exec singularity exec --cleanenv \
    --bind "$root:$root,$bundle:$bundle:ro,${BRASS_EXTRA_BINDS:?},$bundle/get_rg_cns.R:$overlay:ro,$implementation:$native_module:ro" \
    "${BRASS_IMAGE:?}" brass.pl \
    -g "${BRASS_REFERENCE:?}" \
    -d "$bundle/HiDepth.bed.gz" -gc "$bundle/vagrent.cache.gz" \
    -vi "$bundle/viral.genomic.fa.2bit" -mi "$bundle/all_ncbi_bacteria" \
    -b "$bundle/500bp_windows.gc.bed.gz" -cb "$bundle/cytoband.txt" \
    -ct "$bundle/CentTelo.tsv" "$@"
