#!/usr/bin/env python3
"""Independently verify reference hashes, grid, labels and direct-base signal means."""
import argparse
import json
import math
from contextlib import ExitStack
from pathlib import Path
import numpy as np
import pyBigWig
from replication_timing import categories, digest, probabilities, write_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    metadata = json.loads((args.reference / 'reference.json').read_text())
    for filename, expected in metadata['outputs_sha256'].items():
        if digest(args.reference / filename) != expected:
            raise ValueError('Output hash mismatch: ' + filename)
    checked, measured, unknown, values = 0, 0, 0, []
    generator = np.random.default_rng(1234)
    with ExitStack() as stack:
        tracks = []
        for source in metadata['sources']:
            path = args.reference / 'source' / source['filename']
            if digest(path,'md5') != source['md5'] or digest(path) != source['sha256']:
                raise ValueError('Original track hash mismatch: ' + str(path))
            tracks.append(stack.enter_context(pyBigWig.open(str(path))))
        bw = stack.enter_context(pyBigWig.open(str(args.reference / 'replication_timing.bw')))
        arrays = stack.enter_context(np.load(args.reference / 'replication_timing.npz', allow_pickle=False))
        for chrom, length in metadata['chromosome_sizes'].items():
            scores, labels = arrays[chrom+'_score'], arrays[chrom+'_label']
            n = math.ceil(length / metadata['bin_bp'])
            if len(scores) != n or len(labels) != n:
                raise ValueError('Bin count mismatch for ' + chrom)
            valid = np.flatnonzero(np.isfinite(scores)); missing = np.flatnonzero(~np.isfinite(scores))
            measured += len(valid); unknown += len(missing); values.append(scores[valid])
            expected_labels = categories(scores, metadata['label_policy'], metadata['thresholds'])
            np.testing.assert_array_equal(labels, expected_labels)
            indices = set([0, n-1] + generator.integers(0,n,size=5).tolist() + valid[:2].tolist() + missing[:2].tolist())
            for index in indices:
                start = index * metadata['bin_bp']; end = min(start+metadata['bin_bp'],length)
                cell_means = []
                for track in tracks:
                    # Direct per-base values, not the builder's stats() implementation.
                    key = 'chr'+chrom if 'chr'+chrom in track.chroms() else chrom
                    if track.chroms()[key] != length:
                        raise ValueError('Source chromosome length mismatch')
                    raw = np.asarray(track.values(key,start,end), dtype=np.float64)
                    finite = np.isfinite(raw)
                    cell_means.append(float(raw[finite].mean()) if finite.mean() >= metadata['min_input_coverage']-1e-12 else np.nan)
                expected = float(np.mean(cell_means))
                np.testing.assert_allclose(scores[index], expected, rtol=1e-10, atol=1e-10, equal_nan=True)
                # BigWig is float32; annotation arrays retain full float64 precision.
                observed_bw = bw.values('chr'+chrom,start,start+1)[0]
                np.testing.assert_allclose(observed_bw, expected, rtol=1e-6, atol=1e-6, equal_nan=True)
                checked += 1
    if measured != metadata['measured_bins'] or unknown != metadata['unknown_bins']:
        raise ValueError('Coverage totals mismatch')
    quartiles = np.quantile(np.concatenate(values),[.25,.75],method='linear')
    np.testing.assert_allclose(quartiles,metadata['diagnostic_quartiles'],rtol=0,atol=0)
    if metadata['thresholds'] is not None:
        cutoffs = np.quantile(np.concatenate(values),probabilities(metadata['label_policy']),method='linear')
        np.testing.assert_allclose(cutoffs,metadata['thresholds'],rtol=0,atol=0)
    result = dict(status='verified', reference=str(args.reference.resolve()),
                  reference_sha256=digest(args.reference/'reference.json'),
                  total_bins=measured+unknown, measured_bins=measured, unknown_bins=unknown,
                  independently_checked_bins=checked, label_policy=metadata['label_policy'],
                  thresholds=metadata['thresholds'])
    write_json(args.output,result)
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    main()
