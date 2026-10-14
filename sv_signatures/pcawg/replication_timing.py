#!/usr/bin/env python3
"""Checksum-pinned PCAWG timing reference construction and breakpoint annotation."""
import argparse
import csv
import gzip
import hashlib
import json
import math
import platform
import shutil
import subprocess
import sys
import time
from contextlib import ExitStack
from pathlib import Path

PRIMARY = [str(i) for i in range(1, 23)] + ['X']
LABELS = ['unknown', 'late', 'mid', 'early']
POLICIES = ['continuous', 'panprostate_quartiles', 'tertiles']


def probabilities(policy):
    if policy == 'tertiles':
        return [1/3, 2/3]
    if policy == 'panprostate_quartiles':
        return [.25, .75]
    if policy == 'continuous':
        return None
    raise ValueError('Unknown timing label policy: ' + policy)


def policy_note(policy):
    if policy == 'tertiles':
        return 'User-selected genome-wide tertiles on the PCAWG three-cell-line mean; not verified published PCAWG cutoffs'
    if policy == 'continuous':
        return 'Continuous reference only; no early/mid/late assignments'
    return '2026 prostate quartile rule applied to the PCAWG three-cell-line mean: a documented adaptation'


def digest(path, algorithm='sha256'):
    h = hashlib.new(algorithm)
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def load_sources(path):
    lock = json.loads(Path(path).read_text())
    if [s['cell_line'] for s in lock['sources']] != ['NHEK', 'GM12878', 'IMR90']:
        raise ValueError('The PCAWG profile requires exactly NHEK, GM12878 and IMR90 in that order')
    for source in lock['sources']:
        if Path(source['filename']).name != source['filename'] or len(source['md5']) != 32:
            raise ValueError('Invalid pinned source filename/checksum')
    return lock


def download(args):
    lock = load_sources(args.sources)
    source = next(s for s in lock['sources'] if s['cell_line'] == args.cell_line)
    target = args.outdir / source['filename']
    args.outdir.mkdir(parents=True, exist_ok=True)
    if target.exists():
        if digest(target, 'md5') != source['md5']:
            raise ValueError('Existing reference has the wrong MD5: ' + str(target))
        return
    partial = target.with_suffix(target.suffix + '.partial')
    # curl preserves TLS verification, fails on HTTP errors, and retries transient failures.
    argv = ['curl', '--fail', '--location', '--silent', '--show-error', '--retry', '3',
            '--connect-timeout', '30', '--max-time', '900', '--output', str(partial),
            lock['download_base'] + source['filename']]
    print(json.dumps({'argv': argv}), flush=True)
    subprocess.run(argv, check=True)
    if digest(partial, 'md5') != source['md5']:
        raise ValueError('Downloaded reference failed official MD5 verification: ' + str(partial))
    partial.replace(target)


def chromosome_sizes(path):
    sizes = {}
    for line in Path(path).read_text().splitlines():
        fields = line.split()
        if len(fields) < 2:
            raise ValueError('Malformed chromosome-size line')
        chrom = fields[0].removeprefix('chr')
        if chrom in PRIMARY:
            if chrom in sizes or int(fields[1]) <= 0:
                raise ValueError('Duplicate/invalid primary chromosome length')
            sizes[chrom] = int(fields[1])
    if set(sizes) != set(PRIMARY):
        raise ValueError('Chromosome sizes must include autosomes and X')
    return {chrom: sizes[chrom] for chrom in PRIMARY}


def binned_track(track, chrom, length, bin_bp, chunk_bins, min_coverage):
    import numpy as np
    key = 'chr' + chrom if 'chr' + chrom in track.chroms() else chrom
    if track.chroms().get(key) != length:
        raise ValueError('Reference assembly/chromosome length mismatch: ' + chrom)
    n = math.ceil(length / bin_bp)
    scores = np.full(n, np.nan, dtype=np.float64)
    coverage = np.zeros(n, dtype=np.float64)
    # exact=True prevents BigWig zoom summaries from approximating the 1-kb means.
    full_bins = length // bin_bp
    for left in range(0, full_bins, chunk_bins):
        right = min(left + chunk_bins, full_bins)
        scores[left:right] = track.stats(key, left * bin_bp, right * bin_bp,
                                        nBins=right-left, type='mean', exact=True)
        coverage[left:right] = track.stats(key, left * bin_bp, right * bin_bp,
                                          nBins=right-left, type='coverage', exact=True)
    if n > full_bins:
        scores[-1] = track.stats(key, full_bins * bin_bp, length, type='mean', exact=True)[0]
        coverage[-1] = track.stats(key, full_bins * bin_bp, length, type='coverage', exact=True)[0]
    coverage = np.nan_to_num(coverage, nan=0.0)
    scores[coverage < min_coverage - 1e-12] = np.nan
    return scores, coverage


def categories(scores, policy, cutoffs=None):
    import numpy as np
    result = np.zeros(scores.shape, dtype=np.uint8)
    valid = np.isfinite(scores)
    if policy == 'continuous':
        return result
    if policy not in ['panprostate_quartiles', 'tertiles'] or cutoffs is None:
        raise ValueError('Invalid category policy/cutoffs')
    result[valid] = 2
    result[valid & (scores < cutoffs[0])] = 1
    result[valid & (scores > cutoffs[1])] = 3
    return result


def build(args):
    import numpy as np
    import pyBigWig
    started = time.time()
    if args.outdir.exists() and any(args.outdir.iterdir()):
        raise ValueError('Refusing nonempty reference output; choose another outdir')
    if args.bin_bp < 1 or args.chunk_bins < 1 or not 0 < args.min_coverage <= 1:
        raise ValueError('Invalid bin size, chunk size or coverage requirement')
    lock = load_sources(args.sources)
    sizes = chromosome_sizes(args.chrom_sizes)
    paths = {p.name: p for p in args.tracks}
    if set(paths) != {s['filename'] for s in lock['sources']}:
        raise ValueError('Track filenames must match the pinned three sources exactly')
    provenance = []
    for s in lock['sources']:
        path = paths[s['filename']]
        if digest(path, 'md5') != s['md5']:
            raise ValueError('Reference MD5 mismatch: ' + str(path))
        provenance.append(dict(s, sha256=digest(path), bytes=path.stat().st_size,
                               url=lock['download_base'] + s['filename']))
    args.outdir.mkdir(parents=True, exist_ok=True)
    scores, coverage = {}, {}
    with ExitStack() as stack:
        tracks = [stack.enter_context(pyBigWig.open(str(paths[s['filename']]))) for s in lock['sources']]
        for chrom, length in sizes.items():
            pieces = [binned_track(t, chrom, length, args.bin_bp, args.chunk_bins,
                                   args.min_coverage) for t in tracks]
            # Equal cell-line weights; never average only the available subset.
            scores[chrom] = np.mean(np.stack([p[0] for p in pieces]), axis=0)
            coverage[chrom] = np.min(np.stack([p[1] for p in pieces]), axis=0)
            print('Binned chromosome', chrom, 'complete bins', int(np.isfinite(scores[chrom]).sum()), flush=True)
    valid_scores = np.concatenate([v[np.isfinite(v)] for v in scores.values()])
    if not len(valid_scores):
        raise ValueError('No bins have valid measurements from all three cell lines')
    # numpy linear quantiles match R quantile(type=7); compute once genome-wide,
    # never on each sample, caller, chromosome or observed breakpoint subset.
    quantiles = np.quantile(valid_scores, [0.25, 0.75], method='linear').tolist()
    probs = probabilities(args.label_policy)
    cutoff = np.quantile(valid_scores, probs, method='linear').tolist() if probs else None
    archive = {}
    qc_rows = []
    with gzip.open(args.outdir / 'replication_timing.bed.tsv.gz', 'wt') as bed, \
            pyBigWig.open(str(args.outdir / 'replication_timing.bw'), 'w') as bw:
        bw.addHeader([('chr' + c, n) for c, n in sizes.items()])
        writer = csv.writer(bed, delimiter='\t', lineterminator='\n')
        writer.writerow(['chrom', 'start', 'end', 'mean_signal', 'minimum_input_coverage', 'label'])
        for chrom, length in sizes.items():
            value = scores[chrom]
            code = categories(value, args.label_policy, cutoff)
            archive[chrom + '_score'] = value
            archive[chrom + '_label'] = code
            for i, score in enumerate(value):
                label = LABELS[code[i]] if cutoff else ('unclassified' if np.isfinite(score) else 'unknown')
                writer.writerow(['chr' + chrom, i * args.bin_bp, min((i+1)*args.bin_bp, length),
                                 format(score, '.17g') if np.isfinite(score) else 'NA',
                                 format(coverage[chrom][i], '.9g'), label])
            for left in range(0, len(value), args.chunk_bins):
                idx = np.arange(left, min(left + args.chunk_bins, len(value)))
                idx = idx[np.isfinite(value[idx])]
                if len(idx):
                    bw.addEntries(['chr' + chrom] * len(idx), (idx * args.bin_bp).tolist(),
                                  ends=np.minimum((idx+1)*args.bin_bp, length).tolist(), values=value[idx].tolist())
            qc_rows.append({'chrom': chrom, 'bins': len(value), 'measured_bins': int(np.isfinite(value).sum()),
                            'unknown_bins': int((~np.isfinite(value)).sum()),
                            'late': int((code == 1).sum()), 'mid': int((code == 2).sum()), 'early': int((code == 3).sum())})
    np.savez_compressed(args.outdir / 'replication_timing.npz', **archive)
    with (args.outdir / 'chromosome_qc.tsv').open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=qc_rows[0].keys(), delimiter='\t', lineterminator='\n')
        writer.writeheader(); writer.writerows(qc_rows)
    metadata = dict(profile=lock['profile'], genome_build=lock['genome_build'],
        coordinate_system='0-based half-open', chromosome_sizes=sizes, bin_bp=args.bin_bp,
        total_bins=sum(r['bins'] for r in qc_rows), measured_bins=len(valid_scores),
        unknown_bins=sum(r['unknown_bins'] for r in qc_rows), sources=provenance,
        source_lock_sha256=digest(args.sources), chromosome_sizes_sha256=digest(args.chrom_sizes),
        code_sha256=digest(__file__), label_policy=args.label_policy, thresholds=cutoff,
        diagnostic_quartiles=quantiles, label_probabilities=probs,
        quantile_method='R type 7 / numpy linear; equal weight per valid bin',
        min_input_coverage=args.min_coverage, missing_policy='unknown; no interpolation, no zero fill',
        original_pcawg_numeric_category_cutoffs=None,
        categorization_note=policy_note(args.label_policy),
        binning_note='Exact covered-base means per cell line, then equal three-cell-line mean; clipped terminal bins retained',
        elapsed_seconds=time.time()-started, python=platform.python_version(), numpy=np.__version__, pyBigWig=pyBigWig.__version__,
        outputs_sha256={p.name: digest(p) for p in args.outdir.iterdir() if p.is_file()})
    write_json(args.outdir / 'reference.json', metadata)
    print(json.dumps({k: metadata[k] for k in ['total_bins', 'measured_bins', 'unknown_bins', 'label_policy', 'thresholds']}), flush=True)


def relabel(args):
    """Reuse measured bins without downloading or repeatedly re-binning BigWigs."""
    import numpy as np
    started = time.time()
    original = json.loads((args.reference / 'reference.json').read_text())
    if args.outdir.exists() and any(args.outdir.iterdir()):
        raise ValueError('Refusing nonempty reference output')
    for filename, expected in original['outputs_sha256'].items():
        if digest(args.reference / filename) != expected:
            raise ValueError('Original reference hash mismatch: ' + filename)
    probs = probabilities(args.label_policy)
    with np.load(args.reference / 'replication_timing.npz', allow_pickle=False) as archive:
        scores = {c: archive[c+'_score'].copy() for c in original['chromosome_sizes']}
    valid = np.concatenate([v[np.isfinite(v)] for v in scores.values()])
    if len(valid) != original['measured_bins']:
        raise ValueError('Original coverage totals disagree with arrays')
    cutoff = np.quantile(valid, probs, method='linear').tolist() if probs else None
    arrays, rows = {}, []
    for chrom, score in scores.items():
        label = categories(score, args.label_policy, cutoff)
        arrays[chrom+'_score'] = score; arrays[chrom+'_label'] = label
        rows.append(dict(chrom=chrom, bins=len(score), measured_bins=int(np.isfinite(score).sum()),
                         unknown_bins=int((~np.isfinite(score)).sum()), late=int((label==1).sum()),
                         mid=int((label==2).sum()), early=int((label==3).sum())))
    args.outdir.mkdir(parents=True, exist_ok=True)
    (args.outdir / 'source').mkdir()
    for source in original['sources']:
        path = args.reference / 'source' / source['filename']
        if digest(path) != source['sha256']:
            raise ValueError('Original source hash mismatch: ' + str(path))
        shutil.copyfile(path, args.outdir / 'source' / source['filename'])
    shutil.copyfile(args.reference/'replication_timing.bw', args.outdir/'replication_timing.bw')
    np.savez_compressed(args.outdir/'replication_timing.npz', **arrays)
    with gzip.open(args.reference/'replication_timing.bed.tsv.gz','rt') as source, \
            gzip.open(args.outdir/'replication_timing.bed.tsv.gz','wt') as output:
        reader = csv.DictReader(source,delimiter='\t')
        writer = csv.DictWriter(output,fieldnames=reader.fieldnames,delimiter='\t',lineterminator='\n')
        writer.writeheader()
        total = 0
        for row in reader:
            chrom = row['chrom'].removeprefix('chr')
            i = int(row['start']) // original['bin_bp']
            code = arrays[chrom+'_label'][i]
            row['label'] = LABELS[code] if cutoff else ('unclassified' if np.isfinite(scores[chrom][i]) else 'unknown')
            writer.writerow(row); total += 1
    if total != original['total_bins']:
        raise ValueError('Original BED row count mismatch')
    with (args.outdir/'chromosome_qc.tsv').open('w') as handle:
        writer = csv.DictWriter(handle,fieldnames=rows[0].keys(),delimiter='\t',lineterminator='\n')
        writer.writeheader(); writer.writerows(rows)
    metadata = dict(original)
    metadata.update(label_policy=args.label_policy, label_probabilities=probs, thresholds=cutoff,
                    categorization_note=policy_note(args.label_policy), code_sha256=digest(__file__),
                    derived_from_sha256=digest(args.reference/'reference.json'),
                    original_construction_elapsed_seconds=original['elapsed_seconds'],
                    elapsed_seconds=time.time()-started,
                    outputs_sha256={p.name:digest(p) for p in args.outdir.iterdir() if p.is_file()})
    write_json(args.outdir/'reference.json',metadata)
    print(json.dumps(dict(label_policy=args.label_policy,thresholds=cutoff,
                          category_bins={k:sum(r[k] for r in rows) for k in ['late','mid','early']})),flush=True)


def input_junctions(args):
    """Yield common endpoint tuples from classified 1-based TSV or 0-based BEDPE."""
    if getattr(args, 'bedpe', None):
        with args.bedpe.open() as handle:
            for row in csv.reader(handle,delimiter='\t'):
                if len(row) != 10 or int(row[2])-int(row[1]) != 1 or int(row[5])-int(row[4]) != 1:
                    raise ValueError('Timing annotation requires ten-column, one-base-endpoint BEDPE')
                yield row[6], [(row[0],int(row[1])), (row[3],int(row[4]))]
    else:
        with args.classified_sv.open() as handle:
            reader = csv.DictReader(handle,delimiter='\t')
            required = {'id','V7','V1','V4','V11','V12'}
            if not required.issubset(reader.fieldnames or []):
                raise ValueError('Classified SV table missing columns: ' + ','.join(sorted(required-set(reader.fieldnames or []))))
            for row in reader:
                if row['V7'] != row['id']:
                    raise ValueError('Classification id must match V7')
                yield row['id'], [(row['V1'],int(row['V11'])-1), (row['V4'],int(row['V12'])-1)]


def annotate(args):
    import numpy as np
    metadata = json.loads((args.reference / 'reference.json').read_text())
    if metadata['label_policy'] not in ['panprostate_quartiles', 'tertiles']:
        raise ValueError('Reference has continuous scores only; explicitly select a documented category policy first')
    archive_path = args.reference / 'replication_timing.npz'
    if digest(archive_path) != metadata['outputs_sha256']['replication_timing.npz']:
        raise ValueError('Reference array checksum mismatch')
    if args.output.exists():
        raise ValueError('Refusing existing annotation output')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    unknown, seen = [], set()
    with np.load(archive_path, allow_pickle=False) as arrays, args.output.open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=['id','name_regionA1','dist_regionA1','name_regionA2','dist_regionA2',
                                                   'signalA1','signalA2','statusA1','statusA2'], delimiter='\t', lineterminator='\n')
        writer.writeheader()
        for ident, endpoints in input_junctions(args):
            if not ident or ident in seen:
                raise ValueError('Duplicate/missing junction id: ' + ident)
            seen.add(ident)
            out = {'id': ident}
            for endpoint, (chrom, position) in enumerate(endpoints,1):
                chrom = chrom.removeprefix('chr')
                if chrom not in metadata['chromosome_sizes'] or not 0 <= position < metadata['chromosome_sizes'][chrom]:
                    raise ValueError('Invalid GRCh37 primary-chromosome coordinate for ' + ident)
                index = position // metadata['bin_bp']
                code = int(arrays[chrom + '_label'][index]); score = arrays[chrom + '_score'][index]
                out[f'name_regionA{endpoint}'] = LABELS[code]
                out[f'dist_regionA{endpoint}'] = 0 if code else 'NA'
                out[f'signalA{endpoint}'] = format(score, '.17g') if np.isfinite(score) else 'NA'
                out[f'statusA{endpoint}'] = 'measured' if code else 'missing_reference'
                if not code: unknown.append({'id': ident, 'endpoint': endpoint, 'chrom': chrom, 'position': position+1})
            writer.writerow(out)
    write_json(args.output.with_suffix(args.output.suffix + '.qc.json'),
               dict(junctions=len(seen), unknown_endpoints=unknown, reference=str(args.reference.resolve()),
                    reference_sha256=digest(args.reference / 'reference.json'),
                    input_sha256=digest(args.bedpe if getattr(args,'bedpe',None) else args.classified_sv),
                    input_format='0-based BEDPE' if getattr(args,'bedpe',None) else '1-based classified SV TSV',
                    label_policy=metadata['label_policy'], thresholds=metadata['thresholds'],
                    status='needs_missing_annotation_review' if unknown else 'completed'))
    if unknown:
        raise ValueError(f'{len(unknown)} endpoints lack measured timing; audit written, no labels fabricated')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    d = commands.add_parser('download')
    d.add_argument('--sources', type=Path, required=True)
    d.add_argument('--cell-line', choices=['NHEK','GM12878','IMR90'], required=True)
    d.add_argument('--outdir', type=Path, default=Path('.'))
    d.set_defaults(func=download)
    b = commands.add_parser('build')
    b.add_argument('--sources', type=Path, required=True)
    b.add_argument('--tracks', type=Path, nargs=3, required=True)
    b.add_argument('--chrom-sizes', type=Path, required=True)
    b.add_argument('--outdir', type=Path, required=True)
    b.add_argument('--bin-bp', type=int, default=1000)
    b.add_argument('--chunk-bins', type=int, default=5000)
    b.add_argument('--min-coverage', type=float, default=1.0)
    b.add_argument('--label-policy', choices=POLICIES, default='continuous')
    b.set_defaults(func=build)
    r = commands.add_parser('relabel')
    r.add_argument('--reference',type=Path,required=True)
    r.add_argument('--outdir',type=Path,required=True)
    r.add_argument('--label-policy',choices=POLICIES,default='tertiles')
    r.set_defaults(func=relabel)
    a = commands.add_parser('annotate')
    a.add_argument('--reference', type=Path, required=True)
    source = a.add_mutually_exclusive_group(required=True)
    source.add_argument('--classified-sv', type=Path)
    source.add_argument('--bedpe', type=Path)
    a.add_argument('--output', type=Path, required=True)
    a.set_defaults(func=annotate)
    args = parser.parse_args()
    args.func(args)


if __name__ == '__main__':
    main()
