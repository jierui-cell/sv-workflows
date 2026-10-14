#!/bin/bash -l
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
set -Eeuo pipefail
: "${SV_JOB_ID:?Launch through scripts/workflow.py}"
cd "${SV_PROJECT_DIR:?}"
root=$PWD
cfg=${1:?Expected resolved task configuration}
run=${BRASS_RUN_DIR:-"$root/logs/brass-$SV_JOB_ID"}
mkdir -p "$run/preflight"
exec > >(tee -a "$run/driver.log") 2>&1
step=initialize
trap 'rc=$?; echo "FAILED step=$step line=$LINENO exit=$rc time=$(date -Is)"; exit "$rc"' ERR
echo "START job=$SV_JOB_ID host=$(hostname) time=$(date -Is)"
cp "$cfg" "$run/requested-config.json"
cp scripts/benchmark_brass.sh scripts/brass_with_bundle.sh scripts/measure_command.py "$run/"
cp "${BRASS_MANIFEST:?}" "$run/samplesheet.tsv"
cfg="$run/requested-config.json"
pair=$(jq -er .pair_id "$cfg")
[[ $pair =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*$ ]]
bundle=$(realpath "$(jq -er .bundle "$cfg")")
stats=$(realpath "$(jq -er .bam_stats "$cfg")")
out=${BRASS_OUTDIR:-"$root/results/benchmark/$pair/brass-$SV_JOB_ID"}
input_dir=${BRASS_INPUTS_DIR:-"$out/inputs"}
mkdir -p "$out" "$input_dir"
ref=${BRASS_REFERENCE:?}
test -f "$bundle/COMPATIBILITY_VALIDATED"
python3 - "$bundle/bundle-receipt.json" "$ref" "${BRASS_IMAGE:?}" <<'PY'
import hashlib,json,sys
def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for block in iter(lambda:f.read(1048576),b''): h.update(block)
    return h.hexdigest()
r=json.load(open(sys.argv[1]))
if r.get('status')!='passed' or r.get('reference_sha256')!=sha(sys.argv[2]) or r.get('image_sha256')!=sha(sys.argv[3]):
    raise ValueError('BRASS bundle is not validated against this reference and image; rebuild/revalidate it')
PY
step=sample_statistics
python3 - "$run/samplesheet.tsv" "$pair" "$cfg" "$out" "$stats" "$SV_CPUS" "$input_dir" <<'PY'
import csv, json, math, pathlib, sys
sys.path.insert(0, 'scripts')
from brass_statistics import resolve_statistics
manifest, pair, cfg, dest, stats, cpus, input_dir = sys.argv[1:]
with open(manifest) as handle:
    rows = [r for r in csv.DictReader(handle, delimiter='\t') if r['patient_id'] == pair]
if len(rows) != 1: raise ValueError('Expected one manifest row')
r = rows[0]
with open(cfg) as handle: settings = json.load(handle)
resolved, statistics_text = resolve_statistics(r, settings.get('missing_statistics_policy', 'native_defaults'))
settings.update(resolved)
if type(settings['min_group']) is not int or settings['min_group'] < 1 or not 0 < settings['min_cn_change']:
    raise ValueError('Invalid BRASS scientific thresholds')
out = pathlib.Path(dest)
for role in ('tumor', 'normal'):
    bam = pathlib.Path(r[role+'_bam']).resolve()
    bai = pathlib.Path(r[role+'_bai']).resolve()
    bas = pathlib.Path(stats)/(role+'.bam.bas')
    for f in (bam, bai, bas):
        if not f.is_file() or not f.stat().st_size: raise ValueError('Missing input: '+str(f))
    staged = pathlib.Path(input_dir)/bam.name
    staged.symlink_to(bam)
    pathlib.Path(str(staged)+'.bai').symlink_to(bai)
    pathlib.Path(str(staged)+'.bas').symlink_to(bas)
    settings[role+'_bam'] = str(staged)
(out/'sample-statistics.txt').write_text(statistics_text)
settings.update(sex=r['sex'],
    statistics_source=manifest, sex_chromosome_policy='derived from supplied sex metadata, not measured Y copy number',
    cpus=int(cpus), input_cpus=2, normal_panel='disabled', scientific_defaults='min_group=2; min_cn_change=0.3; remaining pinned native defaults')
(out/'effective-settings.json').write_text(json.dumps(settings, indent=2)+'\n')
PY
cp scripts/brass_statistics.py "$run/"
step=verify_native_statistics
export BRASS_STATISTICS_STATE=$(jq -er .statistics_state "$out/effective-settings.json")
stats_binds="$root:$root"
if [[ $BRASS_STATISTICS_STATE == MIXED ]]; then
    native=/opt/wtsi-cgp/lib/perl5/Sanger/CGP/Brass/Implement.pm
    singularity exec --cleanenv --bind "$root:$root" "${BRASS_IMAGE:?}" cp "$native" "$run/Statistics.original.pm"
    python3 scripts/patch_brass_mixed_statistics.py "$run/Statistics.original.pm" "$run/Statistics.mixed.pm"
    stats_binds="$stats_binds,$run/Statistics.mixed.pm:$native:ro"
    cp scripts/patch_brass_mixed_statistics.py "$run/"
fi
singularity exec --cleanenv --bind "$stats_binds" "${BRASS_IMAGE:?}" \
    perl -MSanger::CGP::Brass::Implement -MJSON::PP -e \
    'my %o=(ascat_summary=>$ARGV[0]); Sanger::CGP::Brass::Implement::get_ascat_summary(\%o); print encode_json(\%o);' \
    "$out/sample-statistics.txt" > "$run/native-statistics.json"
python3 - "$out/effective-settings.json" "$run/native-statistics.json" <<'PY'
import json,sys
settings,native=[json.load(open(p)) for p in sys.argv[1:]]
assert float(native['Acf']) == settings['purity']
assert float(native['Ploidy']) == settings['ploidy']
assert native['PloidyAcfState'] == settings['statistics_state']
assert native['GenderChr'] == settings['gender_chr']
assert native['GenderChrFound'] == settings['gender_chr_found']
print('Native statistics verified:', settings['statistics_state'])
PY
tumor=$(jq -er .tumor_bam "$out/effective-settings.json")
normal=$(jq -er .normal_bam "$out/effective-settings.json")
step=preflight
(cd "$run/preflight"; bash "$root/scripts/validate_pair.sh" "$pair" "$ref.fai" "$tumor" "$normal" "${SV_ALLOW_LOCAL:-0}") > "$run/preflight.log" 2>&1
sha256sum "${BRASS_IMAGE:?}" "$ref" "$bundle/adaptation.json" "$bundle/native-validation.json" "$bundle/get_rg_cns.R" "$stats"/*.bas scripts/benchmark_brass.sh scripts/brass_with_bundle.sh scripts/measure_command.py > "$run/artifacts.sha256"
common=(bash -x "$root/scripts/brass_with_bundle.sh" "$bundle"
  -t "$tumor" -n "$normal" -o "$out/calls" -ss "$out/sample-statistics.txt"
  -s "$(jq -er .species "$cfg")" -as "$(jq -er .assembly "$cfg")"
  -pr "$(jq -er .protocol "$cfg")" -pl "$(jq -er .platform "$cfg")"
  -j "$(jq -er .min_group "$cfg")" -cn "$(jq -er .min_cn_change "$cfg")")
# Keep native stage order; input preparation uses at most two workers.
for step in input cover merge group isize normcn filter split assemble grass tabix; do
    cpus=$SV_CPUS
    [[ $step == input ]] && cpus=2
    cmd=("${common[@]}" -c "$cpus" -p "$step")
    printf '%q ' "${cmd[@]}" >> "$run/commands.sh"
    printf '\n' >> "$run/commands.sh"
    echo "STAGE $step time=$(date -Is)"
    if python3 "$run/measure_command.py" "$run/$step.time.txt" "${cmd[@]}" \
        2>&1 | tee "$run/$step.command.log"; then
        :
    elif [[ $step == tabix ]] && \
         grep -Fq 'tar: intermediates: file changed as we read it' "$run/tabix.command.log" && \
         grep -Fq 'unexpectedly returned exit value 1' "$run/tabix.command.log"; then
        echo 'Native BRASS archive warning: independently validating and recovering the completed callset'
        python3 "$root/scripts/recover_brass_archive.py" "$out/calls" "$run/archive-recovery.json"
    else
        echo "BRASS stage $step failed; no automatic recovery applies" >&2
        exit 1
    fi
done
step=output_qc
shopt -s nullglob
vcfs=("$out/calls/"*.vcf.gz)
(( ${#vcfs[@]} > 0 ))
for vcf in "${vcfs[@]}"; do
    bcftools view -Ov -o /dev/null "$vcf"
    bcftools stats "$vcf" > "$run/$(basename "$vcf").stats.txt"
    records=$(bcftools view -H "$vcf" | awk 'END {print NR+0}')
    printf '%s\t%s\n' "$vcf" "$records" >> "$run/record-counts.tsv"
    if (( records == 0 )); then
        echo "REVIEW_REQUIRED: empty BRASS VCF; not evidence of biological absence" | tee "$run/REVIEW_REQUIRED"
    fi
    sha256sum "$vcf" >> "$run/outputs.sha256"
done
if [[ -f $run/REVIEW_REQUIRED ]]; then
    echo "COMPUTATION_COMPLETE_REVIEW_REQUIRED output=$out time=$(date -Is)"
else
    echo "COMPLETE output=$out time=$(date -Is)"
fi
