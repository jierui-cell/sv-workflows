#!/usr/bin/env python3
"""Run pinned ClusterSV, optionally build annotated cSV features and extract NMF signatures."""
import argparse
import csv
import json
import shutil
import subprocess
import sys
import time
from collections import Counter
from pathlib import Path

from run_group import sha256
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from reference_policy import check_vcf
from vcf_to_bedpe import CLUSTERSV_FIELDS, convert
from clustersv_features import build as build_cnfree, PROFILE

CLUSTER_COLUMNS = CLUSTERSV_FIELDS + ('cluster_id', 'cluster_size', 'footprint1', 'footprint2',
                                    'footprint1_bounds', 'footprint2_bounds', 'last_agglomeration_value')
REQUIRED = {'cohort', 'pair_id', 'caller', 'vcf', 'source_status', 'source_sha256'}
CLASS_REQUIRED = {'cohort', 'pair_id', 'caller', 'classified_sv', 'replication_timing',
                  'fragile_sites', 'refined_bedpe', 'clusters_sha256', 'genome_build', 'source_sha256'}
PRIMARY = [str(i) for i in range(1, 23)] + ['X']


def table(path, required):
    with Path(path).open(newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError('{} missing columns: {}'.format(path, ', '.join(sorted(missing))))
        data = list(reader)
    if any(None in row or any(row.get(k) is None for k in required) for row in data):
        raise ValueError('Malformed tabular row in ' + str(path))
    return data


def write_table(path, fields, rows):
    with Path(path).open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def references(source, output):
    """Use upstream hg19 masks, not its malformed hs37d5 mask rows."""
    sizes = dict(line.split() for line in (source / 'references/hs37d5.chrom_sizes').read_text().splitlines())
    mask_file = source / 'references/hg19_centromere_and_telomere_coords.txt'
    masks = {}
    for line in mask_file.read_text().splitlines()[1:]:
        fields = line.split('#', 1)[0].split()
        if not fields:
            continue
        chrom = fields[0]
        values = [int(x) for x in fields[1:]]
        if len(values) != 5 or values != sorted(values) or not 0 <= values[0] < values[-1] <= int(sizes[chrom]):
            raise ValueError('Invalid centromere/telomere mask for ' + chrom)
        masks[chrom] = values
    if set(masks) != set(PRIMARY):
        raise ValueError('Expected GRCh37 autosomes + X masks')
    sizes_out = output / 'primary.chrom_sizes'
    sizes_out.write_text(''.join('{}\t{}\n'.format(c, sizes[c]) for c in PRIMARY))
    return sizes_out, mask_file, masks, sizes


def mask_bedpe(path, masks, sizes, excluded_path, policy='local_cen_telo_exclusion'):
    if policy not in ('local_cen_telo_exclusion','primary_only'):
        raise ValueError('Unknown mask policy: ' + policy)
    rows = []
    excluded = []
    names = set()
    for line in path.read_text().splitlines():
        fields = line.split('\t')
        if len(fields) != 10:
            raise ValueError('Expected ten-column headerless ClusterSV BEDPE')
        if fields[6] in names or not fields[6] or fields[8] not in ('+', '-') or fields[9] not in ('+', '-'):
            raise ValueError('Duplicate/missing junction ID or invalid orientation: ' + line)
        names.add(fields[6])
        reasons = []
        for c, start, end in ((fields[0], int(fields[1]), int(fields[2])),
                              (fields[3], int(fields[4]), int(fields[5]))):
            if c not in masks:
                reasons.append('outside_autosomes_X')
                continue
            if not 0 <= start < end <= int(sizes[c]):
                raise ValueError('Out-of-reference breakpoint: ' + line)
            ptel, cen_start, _, cen_end, qtel = masks[c]
            if policy == 'local_cen_telo_exclusion' and (end <= ptel or end >= qtel):
                reasons.append('telomere_mask')
            if policy == 'local_cen_telo_exclusion' and cen_start <= end <= cen_end:
                reasons.append('centromere_mask')
        if reasons:
            excluded.append(dict(record_id=fields[6], decision=','.join(sorted(set(reasons)))))
        else:
            rows.append(line)
    path.write_text(''.join(row + '\n' for row in rows))
    write_table(excluded_path, ['record_id', 'decision'], excluded)
    return len(rows), len(excluded)


def check_clusters(native, expected, output):
    with native.open() as handle:
        raw = list(csv.reader(handle, delimiter='\t'))
    if len(raw) != expected or any(len(row) != len(CLUSTER_COLUMNS) for row in raw):
        raise ValueError('ClusterSV output count/schema mismatch: ' + str(native))
    counts = Counter(row[10] for row in raw)
    names = [row[6] for row in raw]
    if len(set(names)) != len(names):
        raise ValueError('Duplicate junction IDs in ClusterSV result')
    if any(int(row[11]) != counts[row[10]] for row in raw):
        raise ValueError('ClusterSV membership/size mismatch')
    write_table(output, CLUSTER_COLUMNS, [dict(zip(CLUSTER_COLUMNS, row)) for row in raw])
    return {'junctions': len(raw), 'clusters': len(counts),
            'multi_junction_clusters': sum(n > 1 for n in counts.values()),
            'clustered_junctions': sum(n for n in counts.values() if n > 1),
            'singleton_junctions': sum(n == 1 for n in counts.values())}


def check_revision(source, expected):
    revision = subprocess.check_output(['git', '-C', str(source), 'rev-parse', 'HEAD'], text=True).strip()
    if revision != expected:
        raise ValueError('Unreviewed upstream revision: ' + revision)
    changed = subprocess.check_output(['git', '-C', str(source), 'status', '--porcelain'], text=True)
    if changed.strip():
        raise ValueError('Upstream checkout contains local modifications: ' + str(source))


def run(args):
    out = args.output.resolve()
    if out.exists() and any(out.iterdir()):
        raise ValueError('Refusing nonempty output: ' + str(out))
    out.mkdir(parents=True, exist_ok=True)
    upstream = json.loads(args.upstream_lock.read_text())
    check_revision(args.clustersv_source, upstream['clustersv']['commit'])
    rows = [r for r in table(args.manifest, REQUIRED)
            if r['cohort'] == args.cohort and r['caller'] == args.caller
            and (not args.pair_ids or r['pair_id'] in args.pair_ids.split(','))]
    if len(rows) != args.expected_pairs or len({r['pair_id'] for r in rows}) != len(rows):
        raise ValueError('Incorrect or duplicate selected sample set')
    timing_metadata = None
    if args.timing_reference:
        timing_metadata = json.loads((args.timing_reference / 'reference.json').read_text())
        if timing_metadata['genome_build'] != 'GRCh37' or timing_metadata['label_policy'] != 'tertiles':
            raise ValueError('This workflow requires the verified GRCh37 tertile timing reference')
        verification = json.loads((args.timing_reference / 'verification.json').read_text())
        if verification.get('status') != 'verified' or verification.get('reference_sha256') != sha256(args.timing_reference / 'reference.json'):
            raise ValueError('Timing reference verification is absent or stale')
        if sha256(args.timing_reference / 'replication_timing.npz') != timing_metadata['outputs_sha256']['replication_timing.npz']:
            raise ValueError('Timing reference arrays changed after verification')
    class_rows = []
    if args.stage in ('features', 'signatures') and getattr(args,'classified_manifest',None):
        if not args.classified_manifest:
            raise ValueError('PCAWG 44-feature construction requires --classified-manifest: cSVc classes + annotations + source checksums')
        check_revision(args.panprostate_source, upstream['panprostate']['commit'])
        selected = {r['pair_id'] for r in rows}
        class_rows = [r for r in table(args.classified_manifest, CLASS_REQUIRED)
                      if r['cohort'] == args.cohort and r['caller'] == args.caller and r['pair_id'] in selected]
        if {r['pair_id'] for r in class_rows} != selected or len(class_rows) != len(selected):
            raise ValueError('One classification manifest row is required for every selected pair')
        for r in class_rows:
            if r['genome_build'] != 'GRCh37':
                raise ValueError('Classification genome_build must be GRCh37')
            fields = ['classified_sv', 'fragile_sites', 'refined_bedpe']
            if not args.timing_reference:
                fields.append('replication_timing')
            for field in fields:
                p = Path(r[field])
                if not p.is_absolute() or not p.is_file():
                    raise ValueError('Classification input must be an existing absolute file: ' + str(p))
    sizes_path, mask_path, masks, sizes = references(args.clustersv_source, out)
    for dirname in ('bedpe', 'audit', 'clusters'):
        (out / dirname).mkdir()
    function_source = args.clustersv_source / 'R/clustering_index.R'
    effective = {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()}
    effective.update({'method': 'pcawg_clustersv', 'upstream': upstream,
                      'feature_profile': None if args.stage == 'clustering' else ('panprostate_2026_core' if class_rows else PROFILE),
                      'coordinate_system': '0-based half-open BEDPE, 1bp intervals',
                      'chromosomes': PRIMARY, 'mask_file': str(mask_path),
                      'clustering_function_sha256': sha256(function_source),
                      'r_environment_lock_sha256': sha256(args.rscript.parent / 'renv.lock'),
                      'mask_sha256': sha256(mask_path), 'chrom_sizes_sha256': sha256(sizes_path)})
    effective['mask_policy'] = getattr(args,'mask_policy','local_cen_telo_exclusion')
    effective['mask_policy_note'] = 'Local endpoint filtering is separate from upstream footprint boundary references.'
    if timing_metadata:
        effective.update(timing_source='reference_tertiles', timing_label_policy='tertiles',
                         timing_thresholds=timing_metadata['thresholds'],
                         timing_reference_sha256=sha256(args.timing_reference / 'reference.json'))
    (out / 'effective-settings.json').write_text(json.dumps(effective, indent=2) + '\n')
    source_rows, sample_qc = [], []
    commands = out / 'commands.jsonl'
    def command(argv, log, memory_log):
        with commands.open('a') as handle:
            handle.write(json.dumps({'argv': argv, 'log': str(log), 'started_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}) + '\n')
        with log.open('w') as handle:
            subprocess.run([sys.executable, str(args.memory_monitor), str(memory_log)] + argv,
                           stdout=handle, stderr=subprocess.STDOUT, check=True)
    for row in sorted(rows, key=lambda r: r['pair_id']):
        pair, source = row['pair_id'], Path(row['vcf'])
        digest = sha256(source)
        if getattr(args,'genome_build',None):
            check_vcf(source,args.genome_build)
        if row['source_sha256'] and row['source_sha256'] != digest:
            raise ValueError('Source VCF changed: ' + str(source))
        if class_rows:
            classification = next(r for r in class_rows if r['pair_id'] == pair)
            if classification['source_sha256'] != digest:
                raise ValueError('Classifications do not match source VCF SHA-256 for ' + pair)
        bedpe = out / 'bedpe' / (pair + '.bedpe')
        if class_rows:
            shutil.copyfile(classification['refined_bedpe'], bedpe)
            before_mask = len(bedpe.read_text().splitlines())
        else:
            before_mask = convert(source, args.caller, pair, bedpe, out / 'audit' / (pair + '.events.tsv'),
                                  out / 'audit' / (pair + '.conversion.json'), args.min_size_bp, 'clustersv')
        n, excluded = mask_bedpe(bedpe, masks, sizes, out / 'audit' / (pair + '.masks.tsv'),
                                getattr(args,'mask_policy','local_cen_telo_exclusion'))
        prefix = out / 'clusters' / pair
        command(['bash', str(args.rscript), str(args.cluster_wrapper), str(bedpe), str(sizes_path),
                 str(mask_path), str(prefix), str(function_source), str(args.cpus),
                 str(args.fdr_first), str(args.fdr_final), str(args.footprint_fdr)],
                out / 'clusters' / (pair + '.log'), out / 'clusters' / (pair + '.memory.json'))
        qc = check_clusters(Path(str(prefix) + '.sv_clusters_and_footprints.tsv'), n,
                            out / 'clusters' / (pair + '.clusters.tsv'))
        if class_rows and sha256(Path(str(prefix) + '.sv_clusters_and_footprints.tsv')) != classification['clusters_sha256']:
            raise ValueError('Imported classifications used different ClusterSV output/settings for ' + pair)
        if args.timing_reference:
            timing_output = out / 'audit' / (pair + '.replication.tsv')
            input_args = (['--classified-sv', classification['classified_sv']] if class_rows else ['--bedpe', str(bedpe)])
            command([sys.executable, str(args.timing_wrapper), 'annotate', '--reference', str(args.timing_reference),
                     '--output', str(timing_output)] + input_args,
                    out / 'audit' / (pair + '.replication.log'), out / 'audit' / (pair + '.replication.memory.json'))
            if class_rows:
                classification['replication_timing'] = str(timing_output)
        qc.update(pair_id=pair, pre_mask_junctions=before_mask, masked_junctions=excluded)
        sample_qc.append(qc)
        source_rows.append(dict(pair_id=pair, vcf=str(source.resolve()), sha256=digest,
                                source_status=row['source_status'], retained_events=n))
    write_table(out / 'sources.tsv', source_rows[0].keys(), source_rows)
    write_table(out / 'cluster_summary.tsv', sample_qc[0].keys(), sample_qc)
    (out / 'clustering_qc.json').write_text(json.dumps({'status': 'completed', 'samples': len(rows),
        'junctions': sum(r['junctions'] for r in sample_qc),
        'clustered_junctions': sum(r['clustered_junctions'] for r in sample_qc),
        'note': 'Multi-junction membership is not a complex-event mechanism classification'}, indent=2) + '\n')
    status = {'status': 'not_requested', 'reason': 'Clustering stage only'}
    if args.stage in ('features', 'signatures') and not class_rows:
        build_cnfree(out, [r['pair_id'] for r in rows])
        status = {'status':'not_requested','reason':'CN-free features stage only'}
    if args.stage in ('features', 'signatures') and class_rows:
        selected_manifest = out / 'classification_inputs.tsv'
        write_table(selected_manifest, sorted(CLASS_REQUIRED), [{k: r[k] for k in sorted(CLASS_REQUIRED)} for r in class_rows])
        input_hashes = {r[k]: sha256(Path(r[k])) for r in class_rows for k in ('classified_sv', 'replication_timing', 'fragile_sites', 'refined_bedpe')}
        (out / 'classification_sha256.json').write_text(json.dumps(input_hashes, indent=2) + '\n')
        upstream_features = args.panprostate_source / upstream['panprostate']['feature_script']
        command(['bash', str(args.rscript), str(args.features_wrapper), str(selected_manifest),
                 str(upstream_features), str(out)], out / 'features.log', out / 'features.memory.json')
        (out / 'feature_status.json').write_text(json.dumps({
            'status': 'completed', 'profile': 'panprostate_2026_core',
            'channels': 44, 'samples': len(rows), 'matrix': 'PCAWG.features.tsv'
        }, indent=2) + '\n')
        status = {'status': 'not_requested', 'reason': '44-feature matrix only; decomposition explicitly deferred'}
    if args.stage == 'signatures':
        matrix = 'PCAWG.features.tsv' if class_rows else 'ClusterSV38.matrix.tsv'
        command(['bash', str(args.rscript), str(args.nmf_wrapper), str(out / matrix),
                 str(out), str(args.rank_min), str(args.rank_max), str(args.rank), str(args.nrun),
                 str(args.seed), str(args.cpus), str(args.nmf_epsilon), str(args.nmf_max_iter)] +
                ([] if class_rows else [PROFILE, 'ClusterSV.NMF']),
                out / 'nmf.log', out / 'nmf.memory.json')
        status = json.loads((out / 'signature_status.json').read_text())
    (out / 'signature_status.json').write_text(json.dumps(status, indent=2) + '\n')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('manifest', 'output', 'clustersv_source', 'upstream_lock',
                 'rscript', 'cluster_wrapper', 'nmf_wrapper', 'memory_monitor'):
        p.add_argument('--' + name.replace('_', '-'), type=Path, required=True)
    p.set_defaults(panprostate_source=None, features_wrapper=None)
    p.set_defaults(classified_manifest=None)  # Public route constructs features directly from ClusterSV.
    p.add_argument('--timing-reference', type=Path, help='Verified GRCh37 genome-wide tertile reference; annotate automatically')
    p.add_argument('--timing-wrapper', type=Path, default=Path(__file__).parent / 'pcawg/replication_timing.py')
    p.add_argument('--cohort', required=True)
    p.add_argument('--caller', required=True)
    p.add_argument('--expected-pairs', type=int, required=True)
    p.add_argument('--pair-ids', default='')
    p.add_argument('--stage', choices=('clustering','features','signatures'), default='signatures',help='Stop after clustering/features or extract CN-free de novo signatures')
    p.add_argument('--genome-build', choices=('GRCh37',), default='GRCh37')
    p.add_argument('--mask-policy', choices=('local_cen_telo_exclusion','primary_only'), default='local_cen_telo_exclusion')
    p.add_argument('--min-size-bp', type=int, default=1)
    p.add_argument('--cpus', type=int, default=1)
    p.add_argument('--nmf-epsilon', type=float, default=2.220446049250313e-16)
    p.add_argument('--nmf-max-iter', type=int, default=2000)
    for name, default in (('fdr-first', .01), ('fdr-final', .05), ('footprint-fdr', .01)):
        p.add_argument('--' + name, type=float, default=default)
    for name, default in (('rank-min', 2), ('rank-max', 10), ('rank', 0), ('nrun', 100), ('seed', 1234)):
        p.add_argument('--' + name, type=int, default=default)
    args = p.parse_args()
    if (args.cpus < 1 or args.expected_pairs < 1 or args.min_size_bp < 1 or args.nrun < 1 or
            not 0 < args.fdr_first <= args.fdr_final < 1 or not 0 < args.footprint_fdr < 1 or
            args.rank_min < 2 or args.rank_max < args.rank_min or args.rank < 0 or args.rank == 1 or
            args.nmf_epsilon <= 0 or args.nmf_max_iter < 1):
        p.error('Invalid positive resources, FDR cutoffs, or rank range')
    run(args)


if __name__ == '__main__':
    main()
