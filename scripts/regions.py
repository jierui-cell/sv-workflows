#!/usr/bin/env python3
"""Validate optional BED exclusions; derive Manta's inclusion-list complement."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import subprocess
import shutil
from reference_policy import read_fai

def read_exclusions(path, contigs):
    opener = gzip.open if str(path).endswith('.gz') else open
    intervals = {}
    with opener(path, 'rt') as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith(('#','track ','browser ')):
                continue
            fields = line.split()
            if len(fields) < 3:
                raise ValueError(f'BED line {number}: expected chromosome, start, end')
            chrom = fields[0]
            try:
                start, end = map(int, fields[1:3])
            except ValueError:
                raise ValueError(f'BED line {number}: coordinates must be integers')
            if chrom not in contigs:
                raise ValueError(f'BED line {number}: chromosome {chrom} absent from alignment reference; aliases are not inferred')
            if not 0 <= start < end <= contigs[chrom]:
                raise ValueError(f'BED line {number}: interval outside reference bounds')
            intervals.setdefault(chrom, []).append((start,end))
    merged = {}
    for chrom, values in intervals.items():
        result = []
        for start,end in sorted(values):
            if result and start <= result[-1][1]:
                result[-1] = (result[-1][0], max(end,result[-1][1]))
            else:
                result.append((start,end))
        merged[chrom] = result
    if not merged:
        raise ValueError('Exclusion BED is empty; use null to explicitly disable exclusions')
    return merged

def complement(contigs, intervals):
    result = []
    for chrom,length in contigs.items():
        cursor = 0
        for start,end in intervals.get(chrom, []):
            if cursor < start:
                result.append((chrom,cursor,start))
            cursor = end
        if cursor < length:
            result.append((chrom,cursor,length))
    if not result:
        raise ValueError('Exclusion BED removes the entire reference')
    return result

def prepare(path, fai, output, manta=False):
    path = Path(path).resolve()
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    contigs = read_fai(fai)
    intervals = read_exclusions(path, contigs)
    normalized = output/'exclude.normalized.bed'
    normalized.write_text(''.join(f'{chrom}\t{start}\t{end}\n'
        for chrom in contigs for start,end in intervals.get(chrom, [])))
    audit = {'policy':'explicit_user_supplied_exclusions', 'source':str(path),
        'source_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
        'normalized_sha256':hashlib.sha256(normalized.read_bytes()).hexdigest(),
        'excluded_bases':sum(end-start for values in intervals.values() for start,end in values),
        'coordinate_convention':'0-based half-open BED; no inferred contig aliases'}
    if manta:
        if not all(shutil.which(tool) for tool in ('bgzip','tabix')):
            raise ValueError('Manta exclusions require bgzip and tabix on the compute-node PATH')
        included = complement(contigs, intervals)
        bed = output/'manta.call-regions.bed'
        # Tabix requires chromosome-grouped, coordinate-sorted records.
        bed.write_text(''.join(f'{chrom}\t{start}\t{end}\n' for chrom,start,end in included))
        compressed = Path(str(bed)+'.gz')
        with compressed.open('wb') as handle:
            subprocess.run(['bgzip','-c',str(bed)], stdout=handle, check=True)
        subprocess.run(['tabix','-p','bed',str(compressed)], check=True)
        audit['manta_inclusion_sha256'] = hashlib.sha256(compressed.read_bytes()).hexdigest()
        audit['included_bases'] = sum(end-start for _,start,end in included)
        audit['compression_tools'] = {tool:subprocess.check_output([tool,'--version'],stderr=subprocess.STDOUT,text=True).splitlines()[0]
                                      for tool in ('bgzip','tabix')}
    (output/'regions.audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    return normalized

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exclude', required=True, type=Path)
    parser.add_argument('--fai', required=True, type=Path)
    parser.add_argument('--out', required=True, type=Path)
    parser.add_argument('--manta', action='store_true')
    args = parser.parse_args()
    prepare(args.exclude,args.fai,args.out,args.manta)

if __name__ == '__main__': main()
