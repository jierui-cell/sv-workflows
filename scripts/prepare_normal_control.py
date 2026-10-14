#!/usr/bin/env python3
"""Make one isolated SM-relabeled normal BAM for a self-normal plumbing control."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys


def relabel(header, name):
    lines, samples = [], set()
    for line in header.splitlines():
        fields = line.split('\t')
        if fields[0] == '@RG':
            sm = [f for f in fields if f.startswith('SM:')]
            if len(sm) != 1:
                raise ValueError('Each read group requires exactly one SM field')
            samples.add(sm[0][3:])
            fields = ['SM:'+name if f.startswith('SM:') else f for f in fields]
        lines.append('\t'.join(fields))
    if len(samples) != 1 or name in samples:
        raise ValueError('Expected one original sample, distinct from control label')
    return '\n'.join(lines)+'\n'


def main():
    if not os.environ.get('SLURM_JOB_ID'):
        raise RuntimeError('Control preparation requires a compute-node Slurm job')
    p = argparse.ArgumentParser()
    p.add_argument('--input', required=True)
    p.add_argument('--pair', required=True)
    p.add_argument('--out', required=True)
    args = p.parse_args()
    with open(args.input) as handle:
        rows = [row for row in csv.DictReader(handle, delimiter='\t') if row['patient_id'] == args.pair]
    if len(rows) != 1:
        raise ValueError('Expected exactly one source pair')
    row = rows[0]
    source = Path(row['normal_bam']).resolve()
    out = Path(args.out).resolve()
    out.mkdir(parents=True, exist_ok=False)
    if shutil.disk_usage(str(out)).free < source.stat().st_size * 1.2:
        raise ValueError('Insufficient space for isolated control BAM and index')
    pair = 'NN_' + args.pair
    bam = out/(pair+'.pseudo_tumor.bam')
    step = 'header'

    def run(command, binary_output=None):
        with (out/'commands.jsonl').open('a') as handle:
            handle.write(json.dumps({'step': step, 'argv': command})+'\n')
        if binary_output:
            with binary_output.open('xb') as handle:
                subprocess.run(command, stdout=handle, check=True)
            return b''
        return subprocess.check_output(command)

    try:
        header = run(['samtools', 'view', '--no-PG', '-H', str(source)]).decode()
        expected = relabel(header, pair+'_PSEUDO_TUMOR')
        new_header = out/'pseudo_tumor.header.sam'
        new_header.write_text(expected)
        step = 'reheader_only'
        run([sys.executable, str(Path(__file__).with_name('measure_command.py')), str(out/'reheader.time.txt'), 'samtools',
             'reheader', '-P', str(new_header), str(source)], bam)
        step = 'index_and_qc'
        run(['samtools', 'index', '-@', os.environ.get('SLURM_CPUS_PER_TASK', '1'), str(bam)])
        run(['samtools', 'quickcheck', '-v', str(bam)])
        if run(['samtools', 'view', '--no-PG', '-H', str(bam)]).decode() != expected:
            raise ValueError('Unexpected header changes')
        if run(['samtools', 'idxstats', str(source)]) != run(['samtools', 'idxstats', str(bam)]):
            raise ValueError('Alignment-count/index mismatch')
        checks = {}
        for region in ['chr1:1000000-1001000', 'chr11:70000000-70001000', 'chrX:1000000-1001000']:
            before = run(['samtools', 'view', str(source), region])
            after = run(['samtools', 'view', str(bam), region])
            if before != after:
                raise ValueError('Alignment mismatch at '+region)
            checks[region] = hashlib.sha256(before).hexdigest()
        columns = ['patient_id', 'tumor_bam', 'tumor_bai', 'normal_bam', 'normal_bai',
                   'sex', 'tumor_phenotype', 'purity', 'ploidy']
        with (out/'samplesheet.tsv').open('x') as handle:
            writer = csv.DictWriter(handle, fieldnames=columns, delimiter='\t', lineterminator='\n')
            writer.writeheader()
            writer.writerow(dict(zip(columns, [pair, str(bam), str(bam)+'.bai', str(source),
                row['normal_bai'], row['sex'], 'NORMAL_SELF_CONTROL_NOT_TUMOR', '', ''])))
        (out/'SUCCESS.json').write_text(json.dumps({'source_normal': str(source),
            'source_bytes': source.stat().st_size, 'source_header_sha256': hashlib.sha256(header.encode()).hexdigest(),
            'operation': 'SM-only reheader; original BAM unchanged', 'qc': 'full index counts, header, EOF and three sampled read intervals; not a full payload checksum',
            'sampled_alignment_sha256': checks, 'interpretation': 'identical reads, not an independent normal replicate'}, indent=2)+'\n')
    except Exception as exc:
        (out/'FAILED.json').write_text(json.dumps({'step': step, 'error': str(exc)})+'\n')
        raise


if __name__ == '__main__':
    main()
