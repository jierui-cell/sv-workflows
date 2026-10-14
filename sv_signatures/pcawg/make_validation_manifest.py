#!/usr/bin/env python3
"""Synthetic input contract for the annotated Nextflow branch; no real biological inference."""
import argparse
import csv
import subprocess
import sys
from pathlib import Path


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--pipeline-root', type=Path, required=True)
    p.add_argument('--outdir', type=Path, required=True)
    args = p.parse_args()
    root = args.pipeline_root.resolve()
    out = args.outdir.resolve()
    out.mkdir(parents=True, exist_ok=False)
    sys.path.insert(0, str(root / 'sv_signatures'))
    from run_clustersv import CLASS_REQUIRED, REQUIRED, references, write_table
    from run_group import sha256
    sizes, masks, _, _ = references(root / 'resources/clustersv/source', out)
    manifest, classified = [], []
    for i in range(7):
        sample = 'SYNTH_{:03d}'.format(i + 1)
        cls, repl, fragile, lines, vcf_lines = [], [], [], [], []
        for j, (kind, start, span, strands) in enumerate([
                ('del', 10_000_000, [10_000, 50_000, 500_000][i % 3], '+-'),
                ('td', 50_000_000, [25_000, 250_000, 2_000_000][i % 3], '-+')]):
            if i == 6:  # Explicit supplied zero-event sample.
                continue
            end = start + span
            ident = kind
            lines.append('\t'.join(map(str, ['1', start-1, start, '1', end-1, end, ident, 0, *strands])))
            label = 'DEL' if kind == 'del' else 'DUP'
            vcf_lines.append('1\t{}\t{}\tN\t<{}>\t.\tPASS\tSVTYPE={};END={}'.format(start, ident, label, label, end))
            row = {'V' + str(k): '0' for k in range(1, 19)}
            row.update(V1='1', V2=start-1, V3=start, V4='1', V5=end-1, V6=end,
                       V7=ident, V9=strands[0], V10=strands[1], V11=start, V12=end,
                       V13=str(j+1), V18=kind, id=ident, tumour=sample)
            cls.append(row)
            timing = ['early', 'mid', 'late'][(i+j) % 3]
            repl.append(dict(id=ident, name_regionA1=timing, name_regionA2=timing, dist_regionA1=0, dist_regionA2=0))
            fragile.append(dict(id=ident, name_regionA1='none', name_regionA2='none', dist_regionA1=10000000, dist_regionA2=10000000))
        vcf = out / (sample + '.vcf')
        vcf.write_text('##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' + ''.join(x + '\n' for x in vcf_lines))
        bedpe = out / (sample + '.bedpe')
        bedpe.write_text(''.join(x + '\n' for x in lines))
        prefix = out / sample
        with (out / (sample + '.cluster.log')).open('w') as log:
            subprocess.run(['bash', str(root / 'sv_signatures/pcawg/Rscript.sh'),
                            str(root / 'sv_signatures/pcawg/cluster.R'), str(bedpe), str(sizes),
                            str(masks), str(prefix), str(root / 'resources/clustersv/source/R/clustering_index.R'),
                            '1', '.01', '.05', '.01'], stdout=log, stderr=subprocess.STDOUT, check=True)
        clustered = list(csv.reader(Path(str(prefix) + '.sv_clusters_and_footprints.tsv').open(), delimiter='\t'))
        cluster_ids = {r[6]: r[10] for r in clustered}
        for r in cls:
            r['V13'] = cluster_ids[r['V7']]
        paths = {}
        for field, rows in [('classified_sv', cls), ('replication_timing', repl), ('fragile_sites', fragile)]:
            path = out / (sample + '.' + field + '.tsv')
            fields = (['V'+str(k) for k in range(1,19)] + ['id','tumour'] if field == 'classified_sv'
                      else ['id','name_regionA1','name_regionA2','dist_regionA1','dist_regionA2'])
            write_table(path, fields, rows)
            paths[field] = str(path)
        digest = sha256(vcf)
        manifest.append(dict(cohort='synthetic', pair_id=sample, caller='manta', vcf=str(vcf),
                             source_status='synthetic_control', source_sha256=digest))
        classified.append(dict(cohort='synthetic', pair_id=sample, caller='manta', **paths,
                               refined_bedpe=str(bedpe), genome_build='GRCh37', source_sha256=digest,
                               clusters_sha256=sha256(Path(str(prefix) + '.sv_clusters_and_footprints.tsv'))))
    write_table(out / 'input_manifest.tsv', ['cohort','pair_id','caller','vcf','source_status','source_sha256'], manifest)
    write_table(out / 'classified_manifest.tsv', sorted(CLASS_REQUIRED), classified)
    print('Synthetic annotated manifests ready:', out)


if __name__ == '__main__':
    main()
