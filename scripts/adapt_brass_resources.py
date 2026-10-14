#!/usr/bin/env python3
"""Adapt verified GRCh37d5 resources without modifying coordinates or originals."""
import argparse
from collections import Counter
import gzip
import hashlib
import json
from pathlib import Path
import re
import subprocess

ABSENT = {'hs37d5', 'NC_007605'}
CHR = re.compile(r"('_chr'\s*=>\s*')([^']+)(')")


def annotation_only(name):
    # Explicit Ensembl 75 patch/haplotype classes, inventoried from this cache.
    return bool(re.fullmatch(r'HG[0-9A-Z_]+_(PATCH|NOVEL_TEST)|HSCHR[0-9A-Z_]+', name))


def dictionary(path):
    result = {}
    for line in Path(path).read_text().splitlines():
        if line.startswith('@SQ\t'):
            fields = dict(x.split(':', 1) for x in line.split('\t')[1:])
            result[fields['SN']] = (int(fields['LN']), fields['M5'])
    if not result:
        raise ValueError('Empty sequence dictionary')
    return result


def mapping(source, target):
    aliases = {}
    for name, identity in source.items():
        new = 'chr' + name if name in {str(i) for i in range(1, 23)} | {'X', 'Y', 'MT'} else name
        if new not in target:
            if name not in ABSENT:
                raise ValueError('Unexpected source-only sequence: ' + name)
            continue
        if target[new] != identity:
            raise ValueError('Sequence length/MD5 mismatch: ' + name)
        aliases[name] = new
    if set(aliases.values()) != set(target):
        raise ValueError('Not every target sequence has an identical source sequence')
    return aliases


def change_cache(fields, new):
    old = fields[0]
    embedded = CHR.findall(fields[6])
    if not embedded or any(value != old for _, value, _ in embedded):
        raise ValueError('Unexpected embedded exon chromosome for ' + fields[3])
    fields[6] = CHR.sub(lambda m: m.group(1) + new + m.group(3), fields[6])
    fields[0] = new
    return fields


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--source', required=True)
    p.add_argument('--source-dict', required=True)
    p.add_argument('--target-dict', required=True)
    p.add_argument('--out', required=True)
    args = p.parse_args()
    out = Path(args.out)
    out.mkdir(exist_ok=False, parents=True)
    report = {'status': 'building', 'coordinate_changes': 0, 'resources': {},
              'limitations': ['Optional normal-panel groups are NOT adapted or enabled.',
                  'External masks and annotations are inherited, not newly validated biologically.',
                  'BRASS valid_seqs uses CentTelo: calling scope is chr1-22,chrX,chrY.']}
    try:
        target = dictionary(args.target_dict)
        aliases = mapping(dictionary(args.source_dict), target)
        report['sequence_verified_aliases'] = aliases
        report['source_dictionary_sha256'] = sha(args.source_dict)
        report['target_dictionary_sha256'] = sha(args.target_dict)
        (out/'contig-map.tsv').write_text('source\ttarget\tlength\tsequence_md5\n' + ''.join(
            '{}\t{}\t{}\t{}\n'.format(a, b, *target[b]) for a, b in aliases.items()))
        source = Path(args.source)
        brass = source/'CNV_SV_ref/brass'
        vagrent = source/'VAGrENT_ref_GRCh37d5_ensembl_75/vagrent'

        def transform(src, kind):
            dest = out/src.name
            compressed = src.suffix == '.gz'
            raw = out/(src.name[:-3] if compressed else src.name)
            kept, removed, outside_scope, annotation_excluded = Counter(), Counter(), Counter(), Counter()
            previous = {}
            opener = gzip.open if compressed else open
            with opener(str(src), 'rt') as inp, raw.open('x') as result:
                for number, line in enumerate(inp, 1):
                    if line.startswith('#') or (kind == 'centtel' and number == 1):
                        result.write(line)
                        continue
                    fields = line.rstrip('\n').split('\t')
                    old = fields[0]
                    if kind == 'cytoband' and old not in {str(i) for i in range(1, 23)} | {'X', 'Y'}:
                        # This publisher file includes UCSC haplotypes and chrM
                        # (16571 bases, not our 16569-base MT). None are used by
                        # the supplied CentTelo calling scope. Do not lift/clip.
                        if not re.fullmatch(r'M|[0-9XY]+_gl[0-9]+_random|Un_gl[0-9]+|[0-9]+_[a-z0-9]+_hap[0-9]+', old):
                            raise ValueError('Unexpected non-primary cytoband name: ' + old)
                        outside_scope[old] += 1
                        continue
                    if old not in aliases:
                        if kind == 'cache' and annotation_only(old):
                            annotation_excluded[old] += 1
                            continue
                        if old not in ABSENT:
                            raise ValueError('{}:{} unknown contig {}'.format(src, number, old))
                        removed[old] += 1
                        continue
                    new = aliases[old]
                    length = target[new][0]
                    if kind == 'centtel':
                        coords = [int(x) for x in fields[1:5]]
                        if not (0 <= coords[0] <= coords[1] <= coords[2] <= coords[3] <= length):
                            raise ValueError('Invalid centromere/telomere coordinates: ' + line)
                    else:
                        start, end = int(fields[1]), int(fields[2])
                        if not (0 <= start < end <= length):
                            raise ValueError('Out-of-bounds interval {}:{}: {}'.format(src.name, number, line[:120]))
                        if start < previous.get(new, -1):
                            raise ValueError('Unsorted resource: ' + src.name)
                        previous[new] = start
                    if kind == 'cache':
                        fields = change_cache(fields, new)
                    else:
                        fields[0] = new
                    result.write('\t'.join(fields) + '\n')
                    kept[new] += 1
            if not kept:
                raise ValueError('Empty adapted resource: ' + src.name)
            if compressed:
                subprocess.run(['bgzip', str(raw)], check=True)
                subprocess.run(['tabix', '-p', 'bed', str(dest)], check=True)
            report['resources'][src.name] = {'source_sha256': sha(src), 'adapted_sha256': sha(dest),
                'kept_rows_by_contig': dict(kept), 'removed_rows_by_absent_contig': dict(removed),
                'excluded_annotation_only_contigs': dict(annotation_excluded),
                'excluded_noncalling_cytoband_rows': dict(outside_scope)}

        for name, kind in [('CentTelo.tsv', 'centtel'), ('cytoband.txt', 'cytoband'),
                           ('HiDepth.bed.gz', 'bed'), ('500bp_windows.gc.bed.gz', 'bed')]:
            transform(brass/name, kind)
        transform(vagrent/'vagrent.cache.gz', 'cache')
        # Transcript FASTA identifiers are ENST accessions, not genome contigs.
        # Microbial sequence names are unrelated to the human contig mapping.
        for src in [vagrent/'vagrent.fa', vagrent/'vagrent.fa.fai'] + sorted(brass.glob('*.2bit')):
            (out/src.name).symlink_to(src.resolve())
            report['resources'][src.name] = {'action': 'unchanged_symlink', 'sha256': sha(src)}
        report['status'] = 'adapted_requires_native_reader_validation'
    except Exception as exc:
        report['status'] = 'failed'
        report['error'] = str(exc)
        raise
    finally:
        (out/'adaptation.json').write_text(json.dumps(report, indent=2) + '\n')


if __name__ == '__main__':
    main()
