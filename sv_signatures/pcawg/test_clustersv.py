#!/usr/bin/env python3
"""Compute-node positive controls: official ClusterSV equivalence and cSV feature/NMF fixtures."""
import argparse
import csv
from contextlib import nullcontext
import difflib
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pipeline-root', type=Path, required=True)
    parser.add_argument('--clustering-only', action='store_true')
    parser.add_argument('--rscript',help='Use an explicit Rscript executable for distribution-level CI controls')
    parser.add_argument('--clustersv-source',type=Path,help='Reuse an existing pinned checkout without modifying it')
    parser.add_argument('--workdir', type=Path, help='Keep synthetic fixtures and control outputs for review')
    args = parser.parse_args()
    root = args.pipeline_root.resolve()
    sys.path.insert(0, str(root / 'sv_signatures'))
    from run_clustersv import check_clusters, references, write_table
    rscript = root / 'sv_signatures/pcawg/Rscript.sh'
    rcommand = [shutil.which(args.rscript) or args.rscript,'--vanilla'] if args.rscript else ['bash',str(rscript)]
    source = args.clustersv_source or root / 'resources/clustersv/source'
    from run_clustersv import check_revision
    check_revision(source, json.loads((root/'sv_signatures/pcawg/upstream.json').read_text())['clustersv']['commit'])
    if args.workdir:
        args.workdir.mkdir(parents=True, exist_ok=False)
    with (nullcontext(str(args.workdir)) if args.workdir else tempfile.TemporaryDirectory()) as directory:
        tmp = Path(directory)
        sizes, masks, _, _ = references(source, tmp)
        bedpe = tmp / 'positive.bedpe'
        bedpe.write_text('1\t19999999\t20000000\t1\t21999999\t22000000\tclose1\t0\t+\t-\n'
                         '1\t20000099\t20000100\t1\t22000099\t22000100\tclose2\t0\t+\t-\n'
                         '1\t69999999\t70000000\t1\t70399999\t70400000\tfar\t0\t+\t-\n'
                         '2\t119999999\t120000000\t2\t124999999\t125000000\tinv\t0\t+\t+\n'
                         '3\t49999999\t50000000\t4\t69999999\t70000000\ttra\t0\t-\t+\n')
        def cluster(input_path, prefix):
            subprocess.run(rcommand + [str(root / 'sv_signatures/pcawg/cluster.R'),
                            str(input_path), str(sizes), str(masks), str(prefix),
                            str(source / 'R/clustering_index.R'), '1', '.01', '.05', '.01'], check=True)
        prefix = tmp / 'wrapped'
        cluster(bedpe, prefix)
        native = Path(str(prefix) + '.sv_clusters_and_footprints.tsv')
        qc = check_clusters(native, 5, tmp / 'clusters.tsv')
        assert qc['clustered_junctions'] >= 2 and qc['singleton_junctions'] >= 1, qc
        rows = {r['name']: r for r in csv.DictReader((tmp / 'clusters.tsv').open(), delimiter='\t')}
        assert rows['close1']['cluster_id'] == rows['close2']['cluster_id']
        stock = tmp / 'stock'
        subprocess.run(rcommand + [str(source / 'R/run_cluster_sv.R'),
                        '-bedpe', str(bedpe), '-chr', str(sizes), '-cen_telo', str(masks),
                        '-out', str(stock), '-n', '1'], check=True)
        original = Path(str(stock) + '.sv_clusters_and_footprints.tsv').read_text()
        assert native.read_text() == original, 'Wrapper differs from upstream defaults:\n' + ''.join(difflib.unified_diff(original.splitlines(True), native.read_text().splitlines(True)))
        assert Path(str(prefix) + '.sv_distance_pvals').read_text() == Path(str(stock) + '.sv_distance_pvals').read_text()
        singleton = tmp / 'singleton.bedpe'
        singleton.write_text(bedpe.read_text().splitlines()[0] + '\n')
        cluster(singleton, tmp / 'single')
        assert check_clusters(tmp / 'single.sv_clusters_and_footprints.tsv', 1, tmp / 'single.tsv')['singleton_junctions'] == 1
        empty = tmp / 'empty.bedpe'; empty.touch()
        cluster(empty, tmp / 'empty')
        assert check_clusters(tmp / 'empty.sv_clusters_and_footprints.tsv', 0, tmp / 'empty.tsv')['junctions'] == 0
        print('ClusterSV positive control and upstream equivalence passed:', qc, flush=True)
        if args.clustering_only:
            return
        # Include absent template classes, size boundaries, overlapping annotations,
        # and repeated cluster IDs across patients, which break the original script.
        fixtures = {'A': [('del', 100000, 150000, '1'), ('td', 200000, 249999, '2')],
                    'B': [('del', 300000, 309999, '1')],
                    'C': [('inv', 300000, 350000, '7')],
                    'D': [('chromoplexy_cycle:2', 400000, 450000, '7')] * 2,
                    'E': [('chromoplexy_cycle:2', 400000, 450000, '7')] * 2,
                    'F': [('direct_inversion', 600000, 650000, '7')] * 2,
                    'G': [('shard_cycle:2', 1001000, 2000000, '9.1'),
                          ('shard_cycle:2', 2001000, 1000000, '9.1'),
                          ('shard_cycle:2', 3001000, 4000000, '9x1'),
                          ('shard_cycle:2', 4001000, 3000000, '9x1')], 'H': [],
                    'I': [('chromoplexy_cycle:20', 700000, 750000, '10')] * 2,
                    'J': [('shard_cycle:12', 10_000_000 + i * 1_000_000 + 1000,
                           10_000_000 + ((i + 1) % 12) * 1_000_000, '12') for i in range(12)]}
        manifests = []
        for sample, events in fixtures.items():
            d, r, f = [], [], []
            for i, (kind, start, end, cluster_id) in enumerate(events):
                ident = str(i)
                row = {'V' + str(j): '0' for j in range(1, 19)}
                row.update(V1='1', V4='1', V7=ident, V9='+', V10='-', V11=start, V12=end,
                           V13=cluster_id, V18=kind, id=ident, tumour=sample)
                if sample == 'G':
                    row.update(V1='1' if i % 2 == 0 else '2', V4='2' if i % 2 == 0 else '1')
                d.append(row)
                r.append(dict(id=ident, name_regionA1='early' if sample == 'B' else 'mid',
                              name_regionA2='late' if sample == 'B' else 'mid', dist_regionA1=0, dist_regionA2=0))
                f.append(dict(id=ident, name_regionA1='fragile_site' if sample == 'B' else 'none',
                              name_regionA2='none', dist_regionA1=0, dist_regionA2=1000000))
            paths = {}
            for field, content in (('classified_sv', d), ('replication_timing', r), ('fragile_sites', f)):
                path = tmp / (sample + '.' + field + '.tsv')
                fallback = {'classified_sv': ['V' + str(j) for j in range(1, 19)] + ['id', 'tumour'],
                            'replication_timing': ['id', 'name_regionA1', 'name_regionA2', 'dist_regionA1', 'dist_regionA2'],
                            'fragile_sites': ['id', 'name_regionA1', 'name_regionA2', 'dist_regionA1', 'dist_regionA2']}
                write_table(path, content[0].keys() if content else fallback[field], content)
                paths[field] = str(path)
            manifests.append(dict(pair_id=sample, **paths))
        manifest = tmp / 'classified.tsv'
        write_table(manifest, manifests[0].keys(), manifests)
        subprocess.run(['bash', str(rscript), str(root / 'sv_signatures/pcawg/features.R'), str(manifest),
                        str(root / 'resources/pcawg_signatures/panprostate/complex_sv_signatures/get_complex_nmf_sigs.R'),
                        str(tmp)], check=True)
        features = tmp / 'PCAWG.features.tsv'
        counts = {r['Feature']: r for r in csv.DictReader(features.open(), delimiter='\t')}
        assert len(counts) == 44
        for feature, sample in [('del_less_500kb_mid', 'A'), ('td_less_50kb_mid', 'A'),
                                ('del_less_50kb_early', 'B'), ('del_less_50kb_late', 'B'),
                                ('del_fragileSite', 'B'), ('trans_reciprocal', 'D'),
                                ('trans_reciprocal', 'E'), ('inv_reci_less_100kb', 'F')]:
            assert int(counts[feature][sample]) == 1, (feature, sample, counts[feature])
        assert int(counts['tmpIns_less_5kb_cycle']['A']) == 0
        assert int(counts['tmpIns_less_5kb_cycle']['G']) == 4
        assert sum(int(r['H']) for r in counts.values()) == 0
        assert int(counts['chromoplexy_cycle']['I']) == 1
        assert int(counts['tmpIns_less_5kb_cycle']['J']) == 12
        print('44-feature control passed: size boundaries, absent classes, overlap and patient namespaces', flush=True)
        for name in ('nmf1', 'nmf2'):
            directory = tmp / name; directory.mkdir()
            subprocess.run(['bash', str(rscript), str(root / 'sv_signatures/pcawg/nmf.R'), str(features),
                            str(directory), '2', '4', '2', '5', '1234', '1', '2.220446049250313e-16', '2000'], check=True)
            assert json.loads((directory / 'signature_status.json').read_text())['status'] == 'completed'
            activities = list(csv.DictReader((directory / 'PCAWG.NMF.activities.tsv').open(), delimiter='\t'))
            zero = next(r for r in activities if r['Samples'] == 'H')
            assert all(float(value) == 0 for key, value in zero.items() if key != 'Samples')
        for name in ('PCAWG.NMF.signatures.tsv', 'PCAWG.NMF.activities.tsv'):
            assert (tmp / 'nmf1' / name).read_text() == (tmp / 'nmf2' / name).read_text(), 'NMF seed repeat differs'
        print('KL-NMF positive control and deterministic seed repeat passed', flush=True)
        auto = tmp / 'auto'; auto.mkdir()
        with (auto / 'log.txt').open('w') as handle:
            result = subprocess.run(['bash', str(rscript), str(root / 'sv_signatures/pcawg/nmf.R'),
                                     str(features), str(auto), '2', '4', '0', '5', '1234', '1', '2.220446049250313e-16', '2000'],
                                    stdout=handle, stderr=subprocess.STDOUT)
        if result.returncode:
            assert 'No finite cophenetic decline' in (auto / 'log.txt').read_text(), (auto / 'log.txt').read_text()
            print('Automatic rank control correctly declined to infer a rank from a flat survey', flush=True)
        else:
            survey = list(csv.DictReader((auto / 'NMF.rank_survey.tsv').open(), delimiter='\t'))
            declines = [(float(survey[i+1]['cophenetic']) - float(survey[i]['cophenetic']), int(survey[i]['rank']))
                        for i in range(len(survey)-1)]
            declines = [x for x in declines if x[0] < 0]
            drop = min(x[0] for x in declines)
            expected = max(rank for change, rank in declines if change == drop)
            assert json.loads((auto / 'signature_status.json').read_text())['rank'] == expected
            print('Automatic rank survey control passed', flush=True)


if __name__ == '__main__':
    main()
