#!/usr/bin/env bash
# Generic safe pair IDs; historical launcher remains unchanged for active pilots.
# Execute the existing, validated standalone caller inside one Nextflow Slurm task.
set -Eeuo pipefail
export SV_JOB_ID=${SLURM_JOB_ID:-${PBS_JOBID:-${LSB_JOBID:-${JOB_ID:-local-$$}}}}
export SV_CPUS=${SLURM_CPUS_PER_TASK:-1}
spec=${1:?Expected task.json}
task_dir=$PWD
root=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["project_dir"])' "$spec")
tool=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["tool"])' "$spec")
pair=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["id"])' "$spec")
manifest=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["manifest"])' "$spec")
case "$tool" in brass|lumpy) ;; *) echo "Unsupported caller: $tool" >&2; exit 2;; esac
[[ $pair =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*$ ]] || { echo "Unsafe pair ID: $pair" >&2; exit 2; }
[[ -f $manifest && -f $root/main.nf ]] || { echo 'Missing staged manifest or project root' >&2; exit 2; }
manifest=$(realpath -e "$manifest")
mkdir -p "$task_dir/output"
cp "$spec" "$task_dir/output/effective-task.json"
# Historical wrappers resolve code and the pinned reference from SLURM_SUBMIT_DIR.
# Their output/log overrides below keep every new artifact in this task.
export SV_PROJECT_DIR=$root
export SV_CPUS=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["settings"]["cpus"])' "$spec")
if [[ $tool == brass ]]; then
    export BRASS_IMAGE=$(python3 -c 'import json,sys,pathlib; s=json.load(open(sys.argv[1])); print((pathlib.Path(s["project_dir"])/s["settings"]["image"]).resolve())' "$spec")
    export BRASS_REFERENCE=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["settings"]["reference"])' "$spec")
    export BRASS_EXTRA_BINDS=$(python3 - "$spec" "$task_dir" <<'PY'
import csv,json,pathlib,sys
s=json.load(open(sys.argv[1]))
rows=[r for r in csv.DictReader(open(s['manifest']),delimiter='\t') if r['patient_id']==s['id']]
if len(rows)!=1: raise ValueError('Expected exactly one matched pair')
parents={str(pathlib.Path(rows[0][k]).resolve().parent) for k in ('tumor_bam','normal_bam','tumor_bai','normal_bai')}
parents.add(str(pathlib.Path(s['settings']['reference']).resolve().parent))
print(','.join(p+':'+p+':ro' for p in sorted(parents))+','+sys.argv[2]+':'+sys.argv[2])
PY
)
    export BRASS_MANIFEST=$manifest
    export BRASS_STATS_RUN_DIR=$task_dir/output/stats-log
    export BRASS_STATS_OUTDIR=$task_dir/output/stats
    export BRASS_RUN_DIR=$task_dir/output/brass-log
    export BRASS_OUTDIR=$task_dir/output/callset
    export BRASS_QC_OUTDIR=$task_dir/output/qc-log
    export BRASS_INPUTS_DIR=$task_dir/staged-inputs
    python3 - "$spec" "$task_dir/output/brass-config.json" "$BRASS_STATS_OUTDIR" <<'PY'
import json, pathlib, sys
spec = json.load(open(sys.argv[1]))
s = spec['settings']
root = pathlib.Path(spec['project_dir'])
bundle = pathlib.Path(s['bundle'])
if not bundle.is_absolute(): bundle = root / bundle
if not (bundle / 'COMPATIBILITY_VALIDATED').is_file():
    raise ValueError(f'BRASS bundle is not validated: {bundle}')
required = ('missing_statistics_policy','min_group','min_cn_change','species','assembly','protocol','platform')
missing = [k for k in required if k not in s]
if missing: raise ValueError(f'Missing BRASS settings: {missing}')
c = {k: s[k] for k in required}
c.update(pair_id=spec['id'],bundle=str(bundle.resolve()),bam_stats=sys.argv[3])
pathlib.Path(sys.argv[2]).write_text(json.dumps(c,indent=2)+'\n')
PY
    SV_CPUS=2 bash "$root/scripts/prepare_brass_stats_v2.sh" "$pair"
    bash "$root/scripts/benchmark_brass.sh" "$task_dir/output/brass-config.json"
    bash "$root/scripts/qc_brass_first_pair.sh" "$BRASS_OUTDIR/calls"
    cp "$BRASS_QC_OUTDIR/summary.json" "$BRASS_OUTDIR/technical-qc.json"
else
    echo 'LUMPY is discontinued in this release' >&2
    exit 2
    export LUMPY_PAIR=$pair
    export LUMPY_MANIFEST=$manifest
    export LUMPY_RUN_DIR=$task_dir/output/lumpy-log
    export LUMPY_OUTDIR=$task_dir/output/callset
    export LUMPY_SHARED_TMP=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["settings"]["shared_tmp"])' "$spec")
    export LUMPY_GENOTYPE_WORKERS=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["settings"]["genotype_workers"])' "$spec")
    policy=$(python3 -c 'import json,sys,pathlib; s=json.load(open(sys.argv[1])); p=pathlib.Path(s["settings"]["filter_policy"]); print(p if p.is_absolute() else pathlib.Path(s["project_dir"])/p)' "$spec")
    [[ -f $policy ]] || { echo "Missing LUMPY filter policy: $policy" >&2; exit 2; }
    bash "$root/scripts/benchmark_lumpy_v2.sh"
    python3 "$root/scripts/filter_lumpy.py" \
      --vcf "$LUMPY_OUTDIR/joint.genotyped.vcf" \
      --tumor "$(cat "$LUMPY_RUN_DIR/tumor.samples.txt")" \
      --normal "$(cat "$LUMPY_RUN_DIR/normal.samples.txt")" \
      --config "$policy" --out "$task_dir/output/research-filter"
fi
printf '%s\n' "COMPLETE caller=$tool pair=$pair job=$SV_JOB_ID"
