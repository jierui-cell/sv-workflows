#!/bin/bash -l
#SBATCH --time=00:15:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
set -Eeuo pipefail
: "${SV_JOB_ID:?Launch through scripts/workflow.py}"
cd "${SV_PROJECT_DIR:?}"
out=${BRASS_QC_OUTDIR:-"$PWD/logs/brass-output-qc-$SV_JOB_ID"}
mkdir -p "$out"
calls=${1:?Expected BRASS calls directory}
vcfs=("$calls/"*.annot.vcf.gz)
[[ ${#vcfs[@]} == 1 ]]
bcftools view -Ov -o /dev/null "${vcfs[0]}"
bcftools index -n "${vcfs[0]}" > "$out/index-records.txt"
python3 - "$calls" "$out" <<'PY'
import collections, gzip, json, pathlib, re, subprocess, sys
calls,out=map(pathlib.Path,sys.argv[1:])
settings=json.loads((calls.parent/'effective-settings.json').read_text())
expected=[]
for role in ('normal','tumor'):
    header=subprocess.check_output(['samtools','view','-H',settings[role+'_bam']],universal_newlines=True)
    names={field[3:] for row in header.splitlines() if row.startswith('@RG\t') for field in row.split('\t') if field.startswith('SM:')}
    assert len(names)==1
    expected.append(next(iter(names)))
records={}; filters=collections.Counter(); types=collections.Counter(); sample_mapping={}
with gzip.open(next(calls.glob('*.annot.vcf.gz')),'rt') as f:
    for line in f:
        if line.startswith('##SAMPLE='):
            label=re.search(r'ID=([^,>]+)',line).group(1)
            sample_mapping[label]=re.search(r'SampleName=([^,>]+)',line).group(1)
        if line.startswith('#CHROM'):
            samples=line.rstrip().split('\t')[9:]
        if line.startswith('#'): continue
        row=line.rstrip().split('\t'); info=dict(x.split('=',1) if '=' in x else (x,True) for x in row[7].split(';'))
        assert row[2] not in records
        records[row[2]]=(row,info); filters[row[6]]+=1; types[info.get('SVTYPE','sequence-resolved')]+=1
is_control = settings.get('pair_id', '').startswith('NN_')
assert records or is_control, 'Unexpected empty tumour callset'
assert len(samples)==2
# BRASS may sanitize output sample labels; require the exact original or safe names.
safe=[re.sub(r'[^a-zA-Z0-9_.-]','_',s) for s in expected]
if samples == ['NORMAL','TUMOUR']:
    assert [sample_mapping[x] for x in samples] == expected, sample_mapping
else:
    assert set(samples) in (set(expected),set(safe)), (samples,expected,safe)
bad=[]; pairs=set()
for key,(row,info) in records.items():
    if 'MATEID' in info:
        mate=info['MATEID']
        if mate not in records or records[mate][1].get('MATEID') != key: bad.append(key)
        pairs.add(tuple(sorted((key,mate))))
assert not bad, bad[:10]
summary={'status':'technical_qc_passed','records':len(records),'filters':dict(filters),'types':dict(types),'sample_labels':samples,'sample_mapping':sample_mapping,'reciprocal_mate_pairs':len(pairs),'missing_mates':bad,'limitation':'Not a truth benchmark, somatic-specificity validation or manual read review'}
if is_control:
    summary['control_assessment'] = 'no_records_observed' if not records else 'nonzero_control_review_required'
(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary))
PY
