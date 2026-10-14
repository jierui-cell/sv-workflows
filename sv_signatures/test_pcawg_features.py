"""Lightweight orchestration/QC tests; no R, NMF or biological inference is run."""
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from run_clustersv import CLASS_REQUIRED, run, write_table
from run_group import sha256
from verify_pcawg import verify


class MatrixOnlyTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        root = Path(self.tmp.name)
        (root / 'R').mkdir()
        (root / 'R/clustering_index.R').write_text('fixture\n')
        for name in ('source.vcf', 'classes.tsv', 'timing.tsv', 'fragile.tsv',
                     'renv.lock', 'sizes.tsv', 'masks.tsv', 'clustering_index.R'):
            (root / name).write_text('fixture\n')
        self.fields = ['1', '9999999', '10000000', '1', '10009999', '10010000', 'del', '0', '+', '-']
        native = '\t'.join(self.fields + ['1', '1', '1', '1', '10000000-10000000', '10010000-10010000', '0']) + '\n'
        self.native = native
        (root / 'refined.bedpe').write_text('\t'.join(self.fields) + '\n')
        digest = sha256(root / 'source.vcf')
        (root / 'native.tsv').write_text(native)
        write_table(root / 'manifest.tsv', ['cohort', 'pair_id', 'caller', 'vcf', 'source_status', 'source_sha256'], [dict(
            cohort='synthetic', pair_id='PAIR', caller='manta', vcf=str(root / 'source.vcf'),
            source_status='synthetic_control', source_sha256=digest)])
        write_table(root / 'classified.tsv', sorted(CLASS_REQUIRED), [dict(
            cohort='synthetic', pair_id='PAIR', caller='manta', classified_sv=str(root / 'classes.tsv'),
            replication_timing=str(root / 'timing.tsv'), fragile_sites=str(root / 'fragile.tsv'),
            refined_bedpe=str(root / 'refined.bedpe'), clusters_sha256=sha256(root / 'native.tsv'),
            genome_build='GRCh37', source_sha256=digest)])
        (root / 'upstream.json').write_text(json.dumps({
            'clustersv': {'commit': 'fixture'}, 'panprostate': {'commit': 'fixture', 'feature_script': 'upstream.R'}}))
        self.args = SimpleNamespace(output=root / 'output', manifest=root / 'manifest.tsv',
            classified_manifest=root / 'classified.tsv', upstream_lock=root / 'upstream.json',
            clustersv_source=root, panprostate_source=root, rscript=root / 'Rscript.sh',
            cluster_wrapper=root / 'cluster.R', features_wrapper=root / 'features.R', nmf_wrapper=root / 'nmf.R',
            memory_monitor=root / 'monitor.py', timing_reference=None, stage='features', cohort='synthetic',
            caller='manta', pair_ids='', expected_pairs=1, min_size_bp=1, cpus=1,
            fdr_first=.01, fdr_final=.05, footprint_fdr=.01, rank_min=3, rank_max=10, rank=2,
            nrun=5, seed=1234, nmf_epsilon=2.220446049250313e-16, nmf_max_iter=2000)
        self.commands = []

    def fake_process(self, argv, **kwargs):
        self.commands.append(argv)
        if str(self.args.cluster_wrapper) in argv:
            prefix = Path(argv[argv.index(str(self.args.cluster_wrapper)) + 4])
            Path(str(prefix) + '.sv_clusters_and_footprints.tsv').write_text(self.native)
        elif str(self.args.features_wrapper) in argv:
            out = self.args.output
            write_table(out / 'PCAWG.features.tsv', ['Feature', 'PAIR'],
                        [dict(Feature='feature_{:02d}'.format(i), PAIR=int(i == 0)) for i in range(44)])
            (out / 'feature_qc.json').write_text(json.dumps(dict(
                profile='panprostate_2026_core', channels=44, samples=1, feature_total=1)))
        elif str(self.args.nmf_wrapper) in argv:
            (self.args.output / 'signature_status.json').write_text(json.dumps(dict(status='completed', rank=2)))
        else:
            raise AssertionError('Unexpected command: ' + str(argv))

    def execute(self):
        root = self.args.output.parent
        with patch('run_clustersv.check_revision'), patch('run_clustersv.references', return_value=(
                root / 'sizes.tsv', root / 'masks.tsv', {'1': [10000, 120000000, 125000000, 130000000, 149000000]},
                {'1': '150000000'})), patch('run_clustersv.subprocess.run', side_effect=self.fake_process):
            run(self.args)

    def test_matrix_only_never_runs_nmf(self):
        self.execute()
        self.assertEqual(len(self.commands), 2)
        self.assertFalse(any(str(self.args.nmf_wrapper) in argv for argv in self.commands))
        status = json.loads((self.args.output / 'signature_status.json').read_text())
        self.assertEqual(status['status'], 'not_requested')
        self.assertIn('deferred', status['reason'])
        result = verify(self.args.output, {'PAIR'})
        self.assertEqual(result['features'], 'completed')
        self.assertEqual(result['signatures'], 'not_requested')

    def test_nmf_still_requires_explicit_signature_stage(self):
        self.args.stage = 'signatures'
        self.execute()
        self.assertTrue(any(str(self.args.nmf_wrapper) in argv for argv in self.commands))

    def test_cnfree_features_do_not_require_classification_manifest(self):
        self.args.classified_manifest = None
        with patch('run_clustersv.convert', side_effect=lambda source, caller, pair, target, *rest:
                   (target.write_text('\t'.join(self.fields) + '\n'), 1)[1]):
            self.execute()
        self.assertEqual(len(self.commands),1)
        self.assertEqual(verify(self.args.output, {'PAIR'})['features'],'completed')
        self.assertTrue((self.args.output/'ClusterSV38.matrix.tsv').is_file())

    def test_clustering_does_not_require_classifications(self):
        self.args.stage = 'clustering'
        self.args.classified_manifest = None
        with patch('run_clustersv.convert', side_effect=lambda source, caller, pair, target, *rest:
                   (target.write_text('\t'.join(self.fields) + '\n'), 1)[1]):
            self.execute()
        self.assertEqual(verify(self.args.output, {'PAIR'})['features'], 'not_requested')
        self.assertFalse((self.args.output / 'PCAWG.features.tsv').exists())

    def test_nonfinite_counts_are_rejected(self):
        self.execute()
        path = self.args.output / 'PCAWG.features.tsv'
        path.write_text(path.read_text().replace('feature_00\t1', 'feature_00\tnan'))
        with self.assertRaisesRegex(ValueError, 'Invalid 44-feature count'):
            verify(self.args.output, {'PAIR'})

    def test_duplicate_features_are_rejected(self):
        self.execute()
        path = self.args.output / 'PCAWG.features.tsv'
        path.write_text(path.read_text().replace('feature_01', 'feature_00'))
        with self.assertRaisesRegex(ValueError, 'Incorrect feature matrix'):
            verify(self.args.output, {'PAIR'})

    def test_missing_matrix_is_rejected(self):
        self.execute()
        (self.args.output / 'PCAWG.features.tsv').unlink()
        with self.assertRaises(FileNotFoundError):
            verify(self.args.output, {'PAIR'})

    def test_negative_counts_are_rejected(self):
        self.execute()
        path = self.args.output / 'PCAWG.features.tsv'
        path.write_text(path.read_text().replace('feature_00\t1', 'feature_00\t-1'))
        with self.assertRaisesRegex(ValueError, 'Invalid 44-feature count'):
            verify(self.args.output, {'PAIR'})

    def test_feature_qc_total_must_match(self):
        self.execute()
        path = self.args.output / 'feature_qc.json'
        qc = json.loads(path.read_text())
        qc['feature_total'] = 100
        path.write_text(json.dumps(qc))
        with self.assertRaisesRegex(ValueError, 'Feature QC mismatch'):
            verify(self.args.output, {'PAIR'})

    def test_unrequested_nmf_files_are_rejected(self):
        self.execute()
        (self.args.output / 'PCAWG.NMF.activities.tsv').write_text('stale\n')
        with self.assertRaisesRegex(ValueError, 'Unexpected NMF exports'):
            verify(self.args.output, {'PAIR'})


if __name__ == '__main__':
    unittest.main()
