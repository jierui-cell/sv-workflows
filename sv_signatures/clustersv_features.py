#!/usr/bin/env python3
"""One junction, one channel: CN-free ClusterSV38 counts for de novo discovery."""
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

PROFILE = 'clustersv_cnfree38_v1'
SIZE_BINS = ('lt1kb','1kb_10kb','10kb_100kb','100kb_1Mb','1Mb_10Mb','ge10Mb')
STATES = ('singleton','clustered')
KINDS = ('DEL','DUP','INV')
CHANNELS = tuple(f'{state}:{kind}:{size}' for state in STATES
                 for kind in KINDS for size in SIZE_BINS) + tuple(f'{state}:TRA' for state in STATES)

def size_bin(span):
    if span < 1:
        raise ValueError('Intrachromosomal junction must have a positive span')
    for threshold, label in zip((1000,10000,100000,1000000,10000000), SIZE_BINS):
        if span < threshold:
            return label
    return SIZE_BINS[-1]

def channel(row):
    if int(row['cluster_size']) < 1:
        raise ValueError('Cluster size must be positive')
    state = 'clustered' if int(row['cluster_size']) > 1 else 'singleton'
    strands = row['strand1'] + row['strand2']
    if strands not in ('++','--','+-','-+'):
        raise ValueError('Invalid junction orientation')
    if row['chrom1'] != row['chrom2']:
        return state + ':TRA'
    left, right = int(row['start1']), int(row['start2'])
    if right < left:
        strands = strands[::-1]
    kind = {'+-':'DEL','-+':'DUP','++':'INV','--':'INV'}[strands]
    return f'{state}:{kind}:{size_bin(abs(right-left))}'

def build(root, pairs):
    root = Path(root)
    pairs = sorted(pairs)
    if not pairs or len(set(pairs)) != len(pairs):
        raise ValueError('Expected a unique, nonempty sample set')
    counts, hashes = {}, {}
    with (root/'feature_events.tsv').open('w',newline='') as handle:
        writer = csv.writer(handle,delimiter='\t',lineterminator='\n')
        writer.writerow(('sample','junction_id','channel','cluster_id','cluster_size'))
        for pair in pairs:
            source = root/'clusters'/f'{pair}.clusters.tsv'
            hashes[pair] = hashlib.sha256(source.read_bytes()).hexdigest()
            with source.open(newline='') as reader_handle:
                rows = list(csv.DictReader(reader_handle,delimiter='\t'))
            memberships = Counter(row['cluster_id'] for row in rows)
            ids = [row['name'] for row in rows]
            if len(set(ids)) != len(ids) or any(not ident for ident in ids):
                raise ValueError('Duplicate/missing junction ID: ' + pair)
            counts[pair] = Counter()
            for row in rows:
                if int(row['cluster_size']) != memberships[row['cluster_id']]:
                    raise ValueError('Cluster membership mismatch: ' + pair)
                label = channel(row)
                counts[pair][label] += 1
                writer.writerow((pair,row['name'],label,row['cluster_id'],row['cluster_size']))
    with (root/'ClusterSV38.matrix.tsv').open('w',newline='') as handle:
        writer = csv.writer(handle,delimiter='\t',lineterminator='\n')
        writer.writerow(['Feature'] + pairs)
        for label in CHANNELS:
            writer.writerow([label] + [counts[pair][label] for pair in pairs])
    status = dict(status='completed',profile=PROFILE,channels=len(CHANNELS),samples=len(pairs),
                  matrix='ClusterSV38.matrix.tsv',feature_total=sum(sum(c.values()) for c in counts.values()),
                  sample_totals={pair:sum(counts[pair].values()) for pair in pairs},
                  cluster_inputs_sha256=hashes,channel_order=list(CHANNELS),
                  size_boundary_policy='lower-inclusive, upper-exclusive; <1kb retained',
                  counting_unit='one retained junction; reciprocal BND mates collapsed upstream',
                  timing_policy='endpoint timing remains in audit tables; not an NMF channel',
                  interpretation='Local CN-free feature scheme, not the published PCAWG 44-feature matrix')
    (root/'feature_status.json').write_text(json.dumps(status,indent=2)+'\n')
    return status
