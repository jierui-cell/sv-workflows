"""Compare Cue single-sample VCFs; never label unmatched calls as somatic."""
import argparse
import collections
import csv
import hashlib
import json
from pathlib import Path

import yaml

PRIMARY = {'chr' + str(i) for i in range(1, 23)} | {'chrX', 'chrY'}
FIELDS = ('chrom', 'pos', 'end', 'svtype', 'qual', 'gt', 'normal_matches_strict',
          'nearest_same_type_normal_start_bp', 'nearest_same_type_normal_end_bp',
          'no_normal_match_strict')


def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def validate_cue_provenance(tumor_vcf, normal_vcf, pair, self_control):
    expected_id = pair[3:] if self_control and pair.startswith('NN_') else pair
    metadata = []
    for vcf, expected_role in ((tumor_vcf, 'normal' if self_control else 'tumor'),
                               (normal_vcf, 'normal')):
        base = vcf.parent
        effective = json.loads((base / 'effective-settings.json').read_text())
        success = json.loads((base / 'SUCCESS.json').read_text())
        if effective.get('id') != expected_id or effective.get('role') != expected_role:
            raise ValueError('Cue output pair ID or tumor/normal role mismatch: ' + str(vcf))
        if success.get('stage') != 'single_sample_not_somatic':
            raise ValueError('Input is not a completed Cue single-sample call: ' + str(vcf))
        metadata.append((effective, success))
    a, b = metadata
    if a[0]['fai'] != b[0]['fai'] or a[1]['source_commit'] != b[1]['source_commit']:
        raise ValueError('Cue input reference or source version mismatch')
    if a[1]['scope'] != b[1]['scope']:
        raise ValueError('Cue chromosome scopes differ')
    if not self_control and a[0]['bam'] == b[0]['bam']:
        raise ValueError('Tumor and normal Cue outputs point to the same BAM')
    return dict(tumor_bam=a[0]['bam'], normal_bam=b[0]['bam'],
                fai=a[0]['fai'], cue_source_commit=a[1]['source_commit'])


def parse_vcf(path, primary_only):
    calls = []
    with path.open() as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 8 or f[6] != 'PASS':
                continue
            info = dict(item.split('=', 1) for item in f[7].split(';') if '=' in item)
            if 'END' not in info or 'SVTYPE' not in info:
                continue
            if primary_only and f[0] not in PRIMARY:
                continue
            pos, end = int(f[1]), int(info['END'])
            if end <= pos:
                continue
            gt = f[9].split(':')[0] if len(f) > 9 else '.'
            calls.append(dict(chrom=f[0], pos=pos, end=end,
                              svtype=info['SVTYPE'], qual=f[5], gt=gt,
                              raw_fields=f))
    return calls


def reciprocal_overlap(a, b):
    overlap = max(0, min(a['end'], b['end']) - max(a['pos'], b['pos']))
    return min(overlap / (a['end'] - a['pos']), overlap / (b['end'] - b['pos']))


def compatible(a, b, tolerance, min_overlap):
    if a['chrom'] != b['chrom'] or a['svtype'] != b['svtype']:
        return False
    return (abs(a['pos'] - b['pos']) <= tolerance
            and abs(a['end'] - b['end']) <= tolerance
            and reciprocal_overlap(a, b) >= min_overlap)


def compare(tumor, normal, settings):
    tolerances = settings['exploratory_endpoint_bp']
    strict = settings['strict_endpoint_bp']
    min_overlap = settings['min_reciprocal_overlap']
    by_locus = collections.defaultdict(list)
    for n in normal:
        by_locus[(n['chrom'], n['svtype'])].append(n)
    rows = []
    sensitivity = {str(t): 0 for t in tolerances}
    for t in tumor:
        same = by_locus[(t['chrom'], t['svtype'])]
        matches = [n for n in same if compatible(t, n, strict, min_overlap)]
        for tol in tolerances:
            sensitivity[str(tol)] += any(compatible(t, n, tol, min_overlap) for n in same)
        nearest = min(same, key=lambda n: max(abs(t['pos'] - n['pos']),
                                               abs(t['end'] - n['end']))) if same else None
        rows.append(dict(t, normal_matches_strict=len(matches),
                         nearest_same_type_normal_start_bp=(abs(t['pos'] - nearest['pos']) if nearest else ''),
                         nearest_same_type_normal_end_bp=(abs(t['end'] - nearest['end']) if nearest else ''),
                         no_normal_match_strict=int(not matches)))
    return rows, sensitivity


def write_candidate_vcf(source, destination, rows, pair):
    """Write tumor-only review candidates; deliberately emit no PASS somatic records."""
    with source.open() as src, destination.open('w') as out:
        found_header = False
        for line in src:
            if line.startswith('##'):
                out.write(line)
            elif line.startswith('#CHROM'):
                out.write('##source=Cue_pair_comparison_v1\n')
                out.write('##FILTER=<ID=UNVALIDATED,Description="Tumor Cue prediction without a matched normal Cue prediction under the configured rule; somatic origin unverified">\n')
                out.write('##INFO=<ID=CUE_NMATCH,Number=1,Type=Integer,Description="Number of compatible normal Cue predictions under the strict matching rule">\n')
                columns = line.rstrip('\n').split('\t')
                if len(columns) > 9:
                    columns[9] = pair + '_TUMOR'
                out.write('\t'.join(columns) + '\n')
                found_header = True
                break
        if not found_header:
            raise ValueError('Tumor VCF is missing #CHROM header')
        for row in rows:
            if not row['no_normal_match_strict']:
                continue
            fields = row['raw_fields'].copy()
            fields[6] = 'UNVALIDATED'
            fields[7] += ';CUE_NMATCH=0'
            out.write('\t'.join(fields) + '\n')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--tumor', type=Path, required=True)
    ap.add_argument('--normal', type=Path, required=True)
    ap.add_argument('--config', type=Path, required=True)
    ap.add_argument('--outdir', type=Path, required=True)
    ap.add_argument('--pair', required=True)
    ap.add_argument('--self-control', action='store_true')
    args = ap.parse_args()
    settings = yaml.safe_load(args.config.read_text())
    strict = settings['strict_endpoint_bp']
    tolerances = settings['exploratory_endpoint_bp']
    overlap = settings['min_reciprocal_overlap']
    if not isinstance(tolerances, list) or not tolerances or any(not isinstance(x, int) or x < 0 for x in tolerances):
        ap.error('exploratory_endpoint_bp must be a list of nonnegative integers')
    if not isinstance(strict, int) or strict < 0 or strict not in tolerances:
        ap.error('strict_endpoint_bp must be a nonnegative member of exploratory_endpoint_bp')
    if not isinstance(overlap, (float, int)) or not 0 < overlap <= 1:
        ap.error('min_reciprocal_overlap must be within (0, 1]')
    same = args.tumor.resolve() == args.normal.resolve()
    if same != args.self_control:
        ap.error('identical VCF paths require --self-control; self-control requires identical paths')
    provenance = validate_cue_provenance(args.tumor, args.normal, args.pair, args.self_control)
    tumor = parse_vcf(args.tumor, settings['primary_chromosomes_only'])
    normal = parse_vcf(args.normal, settings['primary_chromosomes_only'])
    if not tumor or not normal:
        ap.error('both VCFs must contain PASS calls with SVTYPE and END')
    rows, sensitivity = compare(tumor, normal, settings)
    if args.self_control and any(r['no_normal_match_strict'] for r in rows):
        raise AssertionError('self-comparison must match every call')
    args.outdir.mkdir(parents=True, exist_ok=False)
    with (args.outdir / 'tumor_calls.tsv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)
    write_candidate_vcf(args.tumor, args.outdir / 'cue_somatic_candidates.unvalidated.vcf', rows, args.pair)
    summary = dict(pair=args.pair, status='exploratory_comparison_complete',
                   self_control=args.self_control, tumor_records=len(tumor),
                   normal_records=len(normal), normal_matched_by_tolerance=sensitivity,
                   unmatched_tumor_strict=sum(r['no_normal_match_strict'] for r in rows),
                   settings=settings, tumor_vcf=str(args.tumor.resolve()),
                   normal_vcf=str(args.normal.resolve()),
                   tumor_sha256=sha256(args.tumor), normal_sha256=sha256(args.normal),
                   candidate_vcf='cue_somatic_candidates.unvalidated.vcf',
                   input_provenance=provenance,
                   caveat='No normal Cue prediction is not evidence of somatic origin. '
                          'This is a VCF-level comparison only; read support, low-coverage normals, '
                          'mapping artifacts and CNV breakpoints have not been adjudicated.')
    (args.outdir / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps({k: v for k, v in summary.items() if not k.endswith('sha256')}, indent=2))


if __name__ == '__main__':
    main()
