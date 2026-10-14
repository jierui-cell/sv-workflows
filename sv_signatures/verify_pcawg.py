#!/usr/bin/env python3
"""Verify ClusterSV membership, optional 44-feature matrices and optional NMF exports."""
import argparse
import csv
import json
import math
from pathlib import Path

from run_clustersv import CLUSTER_COLUMNS, REQUIRED, table
from clustersv_features import CHANNELS, PROFILE, channel


def verify(root, expected):
    sources = table(root / 'sources.tsv', {'pair_id', 'retained_events', 'sha256'})
    summary = table(root / 'cluster_summary.tsv', {'pair_id', 'junctions', 'clustered_junctions', 'singleton_junctions'})
    if len(sources) != len(expected) or {r['pair_id'] for r in sources} != expected:
        raise ValueError('Incorrect source sample set')
    if len(summary) != len(expected) or {r['pair_id'] for r in summary} != expected:
        raise ValueError('Incorrect cluster summary sample set')
    summaries = {r['pair_id']: r for r in summary}
    for source in sources:
        pair = source['pair_id']
        clustered = table(root / 'clusters' / (pair + '.clusters.tsv'), set(CLUSTER_COLUMNS))
        counts = {}
        for row in clustered:
            counts[row['cluster_id']] = counts.get(row['cluster_id'], 0) + 1
        observed = sum(n for n in counts.values() if n > 1)
        if (len(clustered) != int(source['retained_events']) or
                int(summaries[pair]['junctions']) != len(clustered) or
                int(summaries[pair]['clustered_junctions']) != observed or
                int(summaries[pair]['singleton_junctions']) != len(clustered) - observed or
                any(int(row['cluster_size']) != counts[row['cluster_id']] for row in clustered)):
            raise ValueError('Membership/count mismatch for ' + pair)
    qc = json.loads((root / 'clustering_qc.json').read_text())
    if qc['status'] != 'completed' or qc['samples'] != len(expected) or qc['junctions'] != sum(int(r['junctions']) for r in summary):
        raise ValueError('Clustering QC mismatch')
    status = json.loads((root / 'signature_status.json').read_text())
    feature_status = 'not_requested'
    cnfree = (root/'ClusterSV38.matrix.tsv').exists()
    channel_count = 38 if cnfree else 44
    prefix = 'ClusterSV.NMF' if cnfree else 'PCAWG.NMF'
    feature_path = root / ('ClusterSV38.matrix.tsv' if cnfree else 'PCAWG.features.tsv')
    if feature_path.exists() or (root / 'feature_status.json').exists() or status['status'] == 'completed':
        with feature_path.open() as handle:
            feature_reader = csv.DictReader(handle, delimiter='\t')
            columns = feature_reader.fieldnames
            features = list(feature_reader)
        if (not columns or columns[0] != 'Feature' or len(features) != channel_count
                or set(columns[1:]) != expected or len(columns[1:]) != len(expected)
                or len({r.get('Feature') for r in features}) != channel_count
                or any(not r.get('Feature') or None in r or any(v is None for v in r.values()) for r in features)):
            raise ValueError('Incorrect feature matrix')
        values = [float(r[s]) for r in features for s in columns[1:]]
        if any(not math.isfinite(v) or v < 0 for v in values):
            raise ValueError('Invalid 44-feature count')
        feature_qc = json.loads((root / ('feature_status.json' if cnfree else 'feature_qc.json')).read_text())
        if (feature_qc.get('profile') != (PROFILE if cnfree else 'panprostate_2026_core') or feature_qc.get('channels') != channel_count
                or feature_qc.get('samples') != len(expected)
                or not math.isclose(feature_qc['feature_total'], math.fsum(values), rel_tol=1e-9, abs_tol=1e-6)):
            raise ValueError('Feature QC mismatch')
        if (root / 'feature_status.json').exists():
            metadata = json.loads((root / 'feature_status.json').read_text())
            if (metadata.get('status') != 'completed' or metadata.get('channels') != channel_count
                    or metadata.get('samples') != len(expected) or metadata.get('profile') != (PROFILE if cnfree else 'panprostate_2026_core')
                    or metadata.get('matrix') != feature_path.name):
                raise ValueError('Feature status mismatch')
        feature_status = 'completed'
        if cnfree:
            from collections import Counter
            if [r['Feature'] for r in features] != list(CHANNELS):
                raise ValueError('Incorrect CN-free channel order')
            for pair in expected:
                observed = Counter(channel(r) for r in table(root/'clusters'/f'{pair}.clusters.tsv',set(CLUSTER_COLUMNS)))
                if any(float(row[pair]) != observed[row['Feature']] for row in features):
                    raise ValueError('CN-free matrix does not reproduce source junctions: ' + pair)
        if status['status'] == 'not_requested' and list(root.glob(prefix + '.*')):
            raise ValueError('Unexpected NMF exports when decomposition was not requested')
    if status['status'] == 'completed':
        signatures = table(root / (prefix+'.signatures.tsv'), {'Feature'})
        activity = table(root / (prefix+'.activities.tsv'), {'Samples'})
        fractions = table(root / (prefix+'.fractions.tsv'), {'Samples'})
        if ({r['Samples'] for r in activity} != expected or len(activity) != len(expected)
                or {r['Samples'] for r in fractions} != expected or len(fractions) != len(expected)):
            raise ValueError('Incorrect activity sample set')
        names = [k for k in signatures[0] if k != 'Feature']
        if (len(names) != status['rank'] or len(signatures) != channel_count
                or [r['Feature'] for r in signatures] != [r['Feature'] for r in features]):
            raise ValueError('NMF rank/channel mismatch')
        for name in names:
            values = [float(r[name]) for r in signatures]
            if any(not math.isfinite(v) or v < 0 for v in values) or abs(sum(values) - 1) > 1e-6:
                raise ValueError('Invalid normalized NMF signature')
        for row in activity + fractions:
            if any(not math.isfinite(float(row[n])) or float(row[n]) < 0 for n in names):
                raise ValueError('Invalid NMF activity')
        if any(abs(sum(float(row[n]) for n in names) - 1) > 1e-6
               and sum(float(row[n]) for n in names) != 0 for row in fractions):
            raise ValueError('Invalid NMF activity fractions')
    elif status['status'] != 'not_requested':
        raise ValueError('Unexpected signature status')
    return dict(samples=len(expected), junctions=qc['junctions'], clustered_junctions=qc['clustered_junctions'],
                features=feature_status, signatures=status['status'], output=str(root.resolve()))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--group-output',type=Path,help='Verify a single completed task output')
    p.add_argument('--expected-count',type=int)
    p.add_argument('--manifest', type=Path, default=Path('sv_signatures/input_manifest.tsv'))
    p.add_argument('--outdir', type=Path, default=Path('results/sv_signatures_pcawg'))
    p.add_argument('--groups', default='all')
    p.add_argument('--pair-ids', default='')
    p.add_argument('--pair-id', default='', help='Single-pair smoke layout')
    p.add_argument('--output', type=Path)
    args = p.parse_args()
    if args.group_output:
        sources = table(args.group_output/'sources.tsv',{'pair_id'})
        pairs = {row['pair_id'] for row in sources}
        if args.expected_count is None or len(pairs) != args.expected_count:
            p.error('Single-group verification requires matching --expected-count')
        payload = dict(status='verified',**verify(args.group_output,pairs))
        text = json.dumps(payload,indent=2)+'\n'
        if args.output:
            args.output.write_text(text)
        print(text,end='')
        return
    rows = table(args.manifest, REQUIRED)
    groups = {}
    if args.pair_id and args.pair_ids:
        p.error('Use pair-id or pair-ids, not both')
    ids = {args.pair_id} if args.pair_id else (set(args.pair_ids.split(',')) if args.pair_ids else None)
    wanted = set(args.groups.split(',')) if args.groups != 'all' else None
    for row in rows:
        name = row['cohort'] + ':' + row['caller']
        if (not wanted or name in wanted) and (not ids or row['pair_id'] in ids):
            groups.setdefault((row['cohort'], row['caller']), set()).add(row['pair_id'])
    if not groups:
        p.error('No selected groups')
    payload = {'status': 'verified', 'method': 'pcawg_clustersv', 'groups': []}
    for (cohort, caller), samples in sorted(groups.items()):
        location = (args.outdir / 'smoke' / cohort / caller / args.pair_id / 'group_output'
                    if args.pair_id else args.outdir / ('pilot' if ids else '') / cohort / caller / 'group_output')
        payload['groups'].append(dict(cohort=cohort, caller=caller, **verify(location, samples)))
    text = json.dumps(payload, indent=2) + '\n'
    if args.output:
        args.output.write_text(text)
    print(text, end='')


if __name__ == '__main__':
    main()
