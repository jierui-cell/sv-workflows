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

def subtract_targets(targets, exclusions):
    """Subtract exclusion intervals from a normalized capture design."""
    result = []
    for chrom, values in targets.items():
        masks = exclusions.get(chrom, [])
        index = 0
        for start, end in values:
            cursor = start
            while index < len(masks) and masks[index][1] <= start:
                index += 1
            scan = index
            while scan < len(masks) and masks[scan][0] < end:
                left, right = masks[scan]
                if cursor < left:
                    result.append((chrom, cursor, min(left,end)))
                cursor = max(cursor, right)
                if right >= end:
                    break
                scan += 1
            index = scan
            if cursor < end:
                result.append((chrom,cursor,end))
    if not result:
        raise ValueError('No capture regions remain after exclusions')
    return result

def prepare(path, fai, output, manta=False, targets=None):
    path = Path(path).resolve() if path else None
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    contigs = read_fai(fai)
    intervals = read_exclusions(path, contigs) if path else {}
    capture = read_exclusions(targets, contigs) if targets else None
    normalized = output/'exclude.normalized.bed'
    normalized.write_text(''.join(f'{chrom}\t{start}\t{end}\n'
        for chrom in contigs for start,end in intervals.get(chrom, [])))
    audit = {'policy':'explicit_capture_targets_and_exclusions' if capture else 'explicit_user_supplied_exclusions', 'source':str(path) if path else None,
        'source_sha256':hashlib.sha256(path.read_bytes()).hexdigest() if path else None,
        'normalized_sha256':hashlib.sha256(normalized.read_bytes()).hexdigest(),
        'excluded_bases':sum(end-start for values in intervals.values() for start,end in values),
        'coordinate_convention':'0-based half-open BED; no inferred contig aliases'}
    if capture:
        audit['target_source'] = str(Path(targets).resolve())
        audit['target_sha256'] = hashlib.sha256(Path(targets).read_bytes()).hexdigest()
        audit['target_bases'] = sum(end-start for values in capture.values() for start,end in values)
    if manta:
        if not all(shutil.which(tool) for tool in ('bgzip','tabix')):
            raise ValueError('Manta region lists require bgzip and tabix on the compute-node PATH')
        included = subtract_targets(capture, intervals) if capture else complement(contigs, intervals)
        order = {chrom:index for index,chrom in enumerate(contigs)}
        included.sort(key=lambda row: (order[row[0]], row[1]))
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
    parser.add_argument('--exclude', type=Path)
    parser.add_argument('--targets', type=Path)
    parser.add_argument('--fai', required=True, type=Path)
    parser.add_argument('--out', required=True, type=Path)
    parser.add_argument('--manta', action='store_true')
    args = parser.parse_args()
    if not args.exclude and not args.targets:
        parser.error('Supply --exclude and/or --targets')
    prepare(args.exclude,args.fai,args.out,args.manta,args.targets)

if __name__ == '__main__': main()
