#!/usr/bin/env python3
"""Export validated, paired SV caller performance tables from one Nextflow outdir."""
import argparse
import csv
import json
from pathlib import Path
import statistics
from datetime import datetime

CALLERS = ('manta', 'delly', 'svaba', 'gridss', 'brass', 'lumpy', 'cue')
FIELDS = ('dataset_id', 'pair_id', 'caller', 'output_class', 'elapsed_seconds',
          'elapsed_hours', 'peak_sampled_resident_kib', 'peak_sampled_resident_gib',
          'memory_method_id', 'allocated_cpus', 'requested_memory_gb',
          'start_utc', 'end_utc', 'overlap_with_discovered_tasks',
          'benchmark_path')


def parse_utc(value):
    if not value.endswith('+00:00'):
        raise ValueError('Benchmark timestamp is not UTC: ' + value)
    plain = value[:-6]
    for pattern in ('%Y-%m-%dT%H:%M:%S.%f', '%Y-%m-%dT%H:%M:%S'):
        try:
            return datetime.strptime(plain, pattern)
        except ValueError:
            pass
    raise ValueError('Malformed benchmark timestamp: ' + value)


def read_report(path, caller, pair):
    record = json.loads(path.read_text())
    if record.get('schema') != 'sv-caller-benchmark/v1':
        raise ValueError(f'Wrong benchmark schema: {path}')
    if record.get('caller') != caller or record.get('pair_id') != pair:
        raise ValueError(f'Caller/pair identity mismatch: {path}')
    if record.get('exit_code') != 0:
        raise ValueError(f'Caller did not complete successfully: {path}')
    if record.get('memory_status') != 'complete' or not record.get('peak_sampled_resident_kib'):
        raise ValueError(f'Memory measurement incomplete: {path}')
    if record.get('elapsed_seconds', 0) <= 0 or record.get('allocated_cpus', 0) <= 0:
        raise ValueError(f'Invalid runtime or CPU allocation: {path}')
    if not record.get('dataset_id') or not record.get('reference'):
        raise ValueError(f'Missing dataset/reference identity: {path}')
    if not record.get('settings_sha256'):
        raise ValueError(f'Missing caller-settings fingerprint: {path}')
    if not record.get('memory_method_id'):
        raise ValueError(f'Missing memory accounting method: {path}')
    return record


def paths_for(root, caller):
    for pair_dir in sorted((root / caller).glob('*')):
        if not pair_dir.is_dir():
            continue
        benchmark = pair_dir / ('benchmark.json' if caller == 'manta' else 'output/benchmark.json')
        if benchmark.is_file():
            yield pair_dir.name, benchmark


def write_tsv(path, fields, rows):
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def summarize(root, destination, callers, samplesheet=None, allow_incomplete=False):
    roots = [root] if isinstance(root, Path) else list(root)
    reports, issues = [], []
    for caller in callers:
        for one_root in roots:
            for pair, path in paths_for(one_root, caller):
                try:
                    reports.append((path, read_report(path, caller, pair)))
                except (OSError, ValueError, KeyError, json.JSONDecodeError) as exc:
                    issues.append(str(exc))
    if not reports:
        issues.append('No completed benchmark reports were found')
    identities = {(r['dataset_id'], r['reference'], r['boundary'], r.get('memory_method_id')) for _, r in reports}
    if len(identities) > 1:
        issues.append('Mixed dataset, reference, benchmark boundary, or memory method: refusing pooled comparison')
    pairs = set()
    if samplesheet:
        with samplesheet.open(newline='') as handle:
            pairs = {r['patient_id'] for r in csv.DictReader(handle, delimiter='\t')}
        extra = {r['pair_id'] for _, r in reports} - pairs
        if extra:
            issues.append('Reports contain pairs absent from samplesheet: ' + ', '.join(sorted(extra)))
    else:
        pairs = {r['pair_id'] for _, r in reports}
    by_key = {}
    for path, record in reports:
        key = (record['pair_id'], record['caller'])
        if key in by_key:
            issues.append(f'Duplicate benchmark for {key}; use a separate outdir per experiment')
        by_key[key] = (path, record)
    for caller in callers:
        config_ids = {r['settings_sha256'] for _, r in reports if r['caller'] == caller}
        output_classes = {r.get('output_class') for _, r in reports if r['caller'] == caller}
        if len(config_ids) > 1 or len(output_classes) > 1:
            issues.append('Mixed settings or output classes for caller ' + caller)
    missing = [(pair, caller) for pair in sorted(pairs) for caller in callers
               if (pair, caller) not in by_key]
    if missing and not allow_incomplete:
        issues.append(f'Incomplete caller-by-pair matrix: {len(missing)} missing runs')
    balanced = [pair for pair in sorted(pairs)
                if all((pair, caller) in by_key for caller in callers)]
    rows = []
    for path, r in reports:
        begin, end = parse_utc(r['start_utc']), parse_utc(r['end_utc'])
        overlaps = sum(1 for other_path, other in reports if other_path != path and
                       other['dataset_id'] == r['dataset_id'] and
                       parse_utc(other['start_utc']) < end and
                       begin < parse_utc(other['end_utc']))
        rows.append(dict(dataset_id=r['dataset_id'], pair_id=r['pair_id'], caller=r['caller'],
                         output_class=r.get('output_class', ''),
                         elapsed_seconds=r['elapsed_seconds'],
                         elapsed_hours=r['elapsed_seconds'] / 3600,
                         peak_sampled_resident_kib=r['peak_sampled_resident_kib'],
                         peak_sampled_resident_gib=r['peak_sampled_resident_kib'] / 1048576,
                         memory_method_id=r['memory_method_id'],
                         allocated_cpus=r['allocated_cpus'],
                         requested_memory_gb=r['requested_memory_gb'],
                         start_utc=r['start_utc'], end_utc=r['end_utc'],
                         overlap_with_discovered_tasks=overlaps,
                         benchmark_path=str(path.resolve())))
    rows.sort(key=lambda x: (x['pair_id'], x['caller']))
    destination.mkdir(parents=True, exist_ok=True)
    write_tsv(destination / 'all_runs.tsv', FIELDS, rows)
    write_tsv(destination / 'balanced_runs.tsv', FIELDS,
              [row for row in rows if row['pair_id'] in balanced])
    write_tsv(destination / 'missing_runs.tsv', ('pair_id', 'caller'),
              [dict(pair_id=p, caller=c) for p, c in missing])
    summary = []
    for caller in callers:
        selected = [row for row in rows if row['caller'] == caller and row['pair_id'] in balanced]
        summary.append(dict(caller=caller, balanced_pairs=len(selected),
                            median_elapsed_hours=statistics.median(row['elapsed_hours'] for row in selected) if selected else '',
                            median_peak_resident_gib=statistics.median(row['peak_sampled_resident_gib'] for row in selected) if selected else '',
                            allocated_cpus=','.join(map(str, sorted({row['allocated_cpus'] for row in selected}))),
                            requested_memory_gb=','.join(map(str, sorted({row['requested_memory_gb'] for row in selected})))))
    write_tsv(destination / 'balanced_summary.tsv',
              ('caller', 'balanced_pairs', 'median_elapsed_hours', 'median_peak_resident_gib',
               'allocated_cpus', 'requested_memory_gb'), summary)
    status = 'review_required' if issues else ('partial' if missing else 'ready')
    manifest = {'status': status,
                'callers': callers, 'dataset_reference_boundary': list(identities),
                'expected_pairs': len(pairs), 'balanced_pairs': balanced,
                'missing_runs': len(missing), 'issues': issues,
                'limitations': ['Cgroup resident memory is sampled once per second, not an exact peak; shared mapped pages can be charged elsewhere',
                                'overlap count sees only benchmark tasks in this outdir; other cluster/storage load is unmeasured',
                                'Different CPU allocations and caller-specific algorithms affect runtime',
                                'Cue and LUMPY outputs are research candidates, not validated somatic calls']}
    (destination / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir', type=Path, required=True, action='append',
                        help='Nextflow output directory; repeat to combine distinct pairs')
    parser.add_argument('--output', type=Path, required=True, help='New summary directory')
    parser.add_argument('--samplesheet', type=Path)
    parser.add_argument('--callers', default=','.join(CALLERS))
    parser.add_argument('--allow-incomplete', action='store_true')
    args = parser.parse_args()
    callers = tuple(args.callers.split(','))
    if not callers or any(c not in CALLERS for c in callers) or len(set(callers)) != len(callers):
        parser.error('Invalid or duplicate caller name')
    result = summarize(args.outdir, args.output, callers, args.samplesheet, args.allow_incomplete)
    print(json.dumps(result, indent=2))
    raise SystemExit(0 if result['status'] in ('ready', 'partial') else 2)


if __name__ == '__main__':
    main()
