#!/usr/bin/env python3
"""Descriptive Manta SV QC; BND mate pairs, not complex rearrangements, are collapsed."""
import argparse
import collections
import gzip
import json
import re
import statistics
from pathlib import Path


def summarize(path):
    records = {}
    with gzip.open(str(path), 'rt') as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            fields = line.rstrip().split('\t')
            info = dict(item.split('=', 1) if '=' in item else (item, True)
                        for item in fields[7].split(';'))
            if fields[2] in records:
                raise ValueError('Duplicate VCF ID')
            records[fields[2]] = (fields, info)
    passed = {k: v for k, v in records.items() if v[0][6] == 'PASS'}
    counts = collections.Counter()
    chrom_pairs = collections.Counter()
    windows = collections.Counter()
    lengths = collections.defaultdict(list)
    seen = set()
    warnings = []
    examples = []
    for key, (f, info) in passed.items():
        mate = info.get('MATEID')
        event = tuple(sorted([key, mate])) if mate else (key,)
        if event in seen:
            continue
        seen.add(event)
        kind = info['SVTYPE']
        endpoints = [(f[0], int(f[1]))]
        if kind == 'BND':
            match = re.search(r'[\[\]]([^\[\]]+):(\d+)[\[\]]', f[4])
            if not match:
                raise ValueError('Unparsed BND ALT: ' + f[4])
            endpoints.append((match.group(1), int(match.group(2))))
            if mate not in passed or passed[mate][1].get('MATEID') != key:
                warnings.append('Missing/nonreciprocal PASS mate: ' + key)
            kind = 'BND_interchromosomal' if endpoints[0][0] != endpoints[1][0] else 'BND_intrachromosomal'
            chrom_pairs['--'.join(sorted([x[0] for x in endpoints]))] += 1
        elif 'END' in info:
            endpoints.append((f[0], int(info['END'])))
            lengths[kind].append(abs(int(info.get('SVLEN', int(info['END']) - int(f[1])))))
        counts[kind] += 1
        for chrom, window in {(c, (pos - 1) // 10000000) for c, pos in endpoints}:
            windows['{}:{}-{}'.format(chrom, window * 10000000 + 1, (window + 1) * 10000000)] += 1
        examples.append({'id': key, 'type': kind, 'endpoints': endpoints,
                         'somatic_score': int(info.get('SOMATICSCORE', 0)),
                         'imprecise': 'IMPRECISE' in info})
    return {'vcf_records': len(records), 'pass_records': len(passed),
            'pass_adjacencies_or_symbolic_events': len(seen), 'counts': dict(counts),
            'mate_warnings': warnings,
            'size_bp': {k: {'min': min(v), 'median': statistics.median(v), 'max': max(v),
                           'at_least_1Mb': sum(x >= 1000000 for x in v)} for k, v in lengths.items()},
            'top_chromosome_pairs': chrom_pairs.most_common(10),
            'top_10Mb_windows_by_events_touching_window': windows.most_common(10),
            'high_score_examples': sorted(examples, key=lambda e: e['somatic_score'], reverse=True)[:10],
            'limitations': ['Mate collapsing does not reconstruct complex events.',
                            'Unnormalized window counts are descriptive, not enrichment tests.',
                            'No gene annotation, copy-number integration, or orthogonal validation.']}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('vcf', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    args.output.write_text(json.dumps(summarize(args.vcf), indent=2) + '\n')
