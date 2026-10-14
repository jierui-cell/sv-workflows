"""BAM-aware screen of Cue tumor-only predictions; emits no PASS somatic calls."""
import argparse
import collections
import csv
import hashlib
import json
import subprocess
from pathlib import Path

import yaml
from summarize_cue_review import analyze


def checksum(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def validate_settings(settings):
    integer_keys = ('breakpoint_tolerance_bp', 'fetch_radius_bp',
                    'depth_window_radius_bp', 'min_mapping_quality',
                    'min_base_quality', 'min_tumor_junction_links',
                    'max_normal_junction_links')
    for key in integer_keys:
        value = settings.get(key)
        if type(value) is not int or value < 0:
            raise ValueError('{} must be a nonnegative integer'.format(key))
    if settings['fetch_radius_bp'] < settings['breakpoint_tolerance_bp']:
        raise ValueError('fetch radius must cover breakpoint tolerance')
    depth = settings.get('min_normal_breakpoint_depth')
    if type(depth) not in (int, float) or depth < 0:
        raise ValueError('min_normal_breakpoint_depth must be nonnegative')
    kinds = settings.get('screened_svtypes')
    if not isinstance(kinds, list) or not kinds or any(x not in ('DEL', 'DUP') for x in kinds) or len(set(kinds)) != len(kinds):
        raise ValueError('screened_svtypes must contain unique DEL and/or DUP values')


def classify(tumor, normal, settings, svtype):
    if svtype not in settings['screened_svtypes']:
        return ['UNASSESSED_TYPE']
    flags = []
    if min(normal['depth_mean_start'], normal['depth_mean_end']) < settings['min_normal_breakpoint_depth']:
        flags.append('LOW_NORMAL_DEPTH')
    if normal['junction_link_unique_reads'] > settings['max_normal_junction_links']:
        flags.append('NORMAL_READ_LINK')
    if tumor['junction_link_unique_reads'] < settings['min_tumor_junction_links']:
        flags.append('LOW_TUMOR_READ_LINK')
    return flags or ['UNVALIDATED']


def read_candidate_vcf(path):
    header, records = [], []
    with path.open() as handle:
        for line in handle:
            if line.startswith('#'):
                header.append(line)
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 8 or fields[6] != 'UNVALIDATED':
                raise ValueError('Expected UNVALIDATED Cue candidate VCF records')
            info = dict(x.split('=', 1) for x in fields[7].split(';') if '=' in x)
            if info.get('CUE_NMATCH') != '0' or 'END' not in info:
                raise ValueError('Candidate missing CUE_NMATCH=0 or END')
            records.append((fields, int(info['END'])))
    if not any(x.startswith('#CHROM') for x in header):
        raise ValueError('VCF missing #CHROM header')
    return header, records


def annotate(fields, tumor, normal, flags):
    result = fields.copy()
    result[6] = ';'.join(flags)
    result[7] += ';CUE_TLINK={};CUE_NLINK={};CUE_NDP1={:.2f};CUE_NDP2={:.2f}'.format(
        tumor['junction_link_unique_reads'], normal['junction_link_unique_reads'],
        normal['depth_mean_start'], normal['depth_mean_end'])
    return '\t'.join(result) + '\n'


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--pair-summary', type=Path, required=True)
    ap.add_argument('--config', type=Path, required=True)
    ap.add_argument('--outdir', type=Path, required=True)
    args = ap.parse_args()
    pair_summary = json.loads(args.pair_summary.read_text())
    if pair_summary.get('status') != 'exploratory_comparison_complete' or pair_summary.get('self_control'):
        ap.error('Expected completed non-self Cue pair comparison')
    settings = yaml.safe_load(args.config.read_text())
    validate_settings(settings)
    vcf = args.pair_summary.parent / pair_summary['candidate_vcf']
    source_header, records = read_candidate_vcf(vcf)
    if len(records) != pair_summary['unmatched_tumor_strict']:
        raise ValueError('Candidate VCF count differs from pair summary')
    provenance = pair_summary['input_provenance']
    bams = {role: Path(provenance[role + '_bam']) for role in ('tumor', 'normal')}
    for role, bam in bams.items():
        if not bam.exists() or not Path(str(bam) + '.bai').exists():
            raise FileNotFoundError('Missing {} BAM or index: {}'.format(role, bam))
        subprocess.run(['samtools', 'quickcheck', '-v', str(bam)], check=True)
    rows, annotated, provisional = [], [], []
    for fields, end in records:
        chrom, pos = fields[0], int(fields[1])
        info = dict(x.split('=', 1) for x in fields[7].split(';') if '=' in x)
        evidence = {}
        for role in ('tumor', 'normal'):
            evidence[role] = analyze(
                bams[role], chrom, pos, end, indexed_regions=True,
                tolerance=settings['breakpoint_tolerance_bp'],
                fetch_radius=settings['fetch_radius_bp'],
                depth_radius=settings['depth_window_radius_bp'],
                min_mapq=settings['min_mapping_quality'],
                min_baseq=settings['min_base_quality'])
        flags = classify(evidence['tumor'], evidence['normal'], settings, info['SVTYPE'])
        line = annotate(fields, evidence['tumor'], evidence['normal'], flags)
        annotated.append(line)
        if flags == ['UNVALIDATED']:
            provisional.append(line)
        rows.append(dict(chrom=chrom, pos=pos, end=end, svtype=info['SVTYPE'],
                         cue_qual=fields[5], filter=';'.join(flags),
                         tumor_links=evidence['tumor']['junction_link_unique_reads'],
                         normal_links=evidence['normal']['junction_link_unique_reads'],
                         tumor_depth_start=evidence['tumor']['depth_mean_start'],
                         tumor_depth_end=evidence['tumor']['depth_mean_end'],
                         normal_depth_start=evidence['normal']['depth_mean_start'],
                         normal_depth_end=evidence['normal']['depth_mean_end']))
    args.outdir.mkdir(parents=True, exist_ok=False)
    header = []
    for line in source_header:
        if line.startswith('#CHROM'):
            for name, description in (
                    ('NORMAL_READ_LINK', 'Normal BAM has a breakpoint-to-breakpoint SA or mate link'),
                    ('LOW_NORMAL_DEPTH', 'Normal breakpoint depth is below the configured minimum'),
                    ('LOW_TUMOR_READ_LINK', 'Tumor BAM has too few breakpoint-to-breakpoint links')):
                header.append('##FILTER=<ID={},Description="{}">\n'.format(name, description))
            header.append('##FILTER=<ID=UNASSESSED_TYPE,Description="This SV type is not assessed by the research BAM screen">\n')
            for name, kind, description in (
                    ('CUE_TLINK', 'Integer', 'Unique tumor reads with SA or mate linkage'),
                    ('CUE_NLINK', 'Integer', 'Unique normal reads with SA or mate linkage'),
                    ('CUE_NDP1', 'Float', 'Mean normal depth at first breakpoint window'),
                    ('CUE_NDP2', 'Float', 'Mean normal depth at second breakpoint window')):
                header.append('##INFO=<ID={},Number=1,Type={},Description="{}">\n'.format(name, kind, description))
        header.append(line)
    for filename, calls in (('cue_candidates.with_bam_evidence.vcf', annotated),
                            ('cue_provisional_somatic.unvalidated.vcf', provisional)):
        (args.outdir / filename).write_text(''.join(header + calls))
    with (args.outdir / 'read_evidence.tsv').open('w', newline='') as handle:
        fields = list(rows[0]) if rows else ('chrom', 'pos', 'end', 'svtype', 'cue_qual', 'filter',
                                           'tumor_links', 'normal_links', 'tumor_depth_start',
                                           'tumor_depth_end', 'normal_depth_start', 'normal_depth_end')
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    summary = dict(status='research_bam_screen_complete', pair=pair_summary['pair'],
                   input_candidate_vcf=str(vcf.resolve()), input_vcf_sha256=checksum(vcf),
                   input_pair_summary=str(args.pair_summary.resolve()),
                   tumor_bam=str(bams['tumor'].resolve()), normal_bam=str(bams['normal'].resolve()),
                   settings=settings, input_candidates=len(records),
                   provisional_candidates=len(provisional),
                   filter_counts=dict(collections.Counter(r['filter'] for r in rows)),
                   caveat='UNVALIDATED records are provisional, not PASS somatic calls. Links ignore orientation/type compatibility and may include mapping artifacts. Thresholds are uncalibrated; no truth benchmark or complete manual IGV review.')
    (args.outdir / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
