"""Synthetic caller-to-signature handoff; no patient data, caller execution or fitting."""
from contextlib import redirect_stderr, redirect_stdout
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import MagicMock, patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import workflow
from check_resources import check_files, SAMPLE_FIELDS
from signature_inputs import discover, write_auto_manifest


class SignatureInputTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        for name in ('config/defaults.yaml', 'config/signatures.yaml', 'main.nf', 'nextflow.config'):
            self.write(name, (ROOT / name).read_text())
        self.caller = {'input': 'inputs/samplesheet.tsv', 'outdir': 'results/cohort_run01',
                       'dataset_id': 'synthetic', 'callers': 'manta,delly,gridss,svaba,brass'}
        self.save_caller()
        self.config = workflow.read_yaml(ROOT / 'config/signatures.yaml')
        self.config['python'] = sys.executable
        self.write('config/signatures.yaml', workflow.yaml.safe_dump(self.config))
        with self.write('inputs/samplesheet.tsv', '').open('w', newline='') as handle:
            writer = csv.writer(handle, delimiter='\t')
            writer.writerow(SAMPLE_FIELDS)
            for pair in ('PAIR_001', 'PAIR_002'):
                # BAM paths deliberately do not exist: signatures consume VCFs, not BAMs.
                writer.writerow([pair, 'inputs/bams/t.bam', 'inputs/bams/t.bam.bai',
                                 'inputs/bams/n.bam', 'inputs/bams/n.bam.bai', '', '', '', ''])

    def write(self, name, content):
        path = self.root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content)
        return path

    def save_caller(self):
        self.write('config/private-cohort.yaml', workflow.yaml.safe_dump(self.caller))

    def complete(self, caller, pair='PAIR_001', compact=False, legacy_manta=False):
        base = Path(self.caller['outdir']) / caller / pair
        task = base if caller == 'manta' or compact else base / 'output'
        vcf = {
            'manta': Path('somaticSV.vcf.gz') if compact else Path('manta/results/variants/somaticSV.vcf.gz'),
            'delly': Path('somaticSV.vcf.gz'),
            'svaba': Path(pair + '.svaba.somatic.sv.vcf'),
            'gridss': Path('somatic/somatic.high_confidence.vcf.gz'),
            'brass': Path('callset/calls/' + pair + '.annot.vcf.gz'),
        }[caller]
        path = self.write(task / vcf, '')
        content = ('##fileformat=VCFv4.2\n##contig=<ID=chr1,length=249250621>\n'
                   '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
                   'chr1\t1000\ttoy\tN\t<DEL>\t100\tPASS\tSVTYPE=DEL;END=3000\n')
        if path.suffix == '.gz':
            with gzip.open(path, 'wt') as handle:
                handle.write(content)
        else:
            path.write_text(content)
        reference = self.root / 'resources/reference/reference.fa'
        benchmark = {'schema': 'sv-caller-benchmark/v1', 'exit_code': 0,
                     'dataset_id': self.caller['dataset_id'], 'caller': caller, 'pair_id': pair,
                     'reference': str(reference), 'output_class': 'somatic',
                     'settings': {'genome_build': 'GRCh37', 'sequencing_type': 'WGS'}}
        benchmark_dir = task / 'output' if caller == 'manta' and not compact and not legacy_manta else task
        self.write(benchmark_dir / 'benchmark.json', json.dumps(benchmark))
        if caller == 'brass':
            self.write(task / 'callset/technical-qc.json', '{"status":"technical_qc_passed"}')
        elif caller != 'manta':
            stage = 'somatic_filtered_no_pon' if caller == 'gridss' else 'somatic'
            self.write(task / 'SUCCESS.json', json.dumps({'stage': stage, 'vcfs': [str(vcf)]}))
            if caller == 'gridss':
                self.write(task / 'somatic/SUCCESS.json', json.dumps({'stage': stage}))
        return path

    def complete_cohort(self, compact=False):
        for caller in self.caller['callers'].split(','):
            for pair in ('PAIR_001', 'PAIR_002'):
                self.complete(caller, pair, compact=compact and caller != 'brass')

    def test_discovers_all_five_callers_and_never_selects_gridss_joint_or_full(self):
        self.complete_cohort()
        for name in ('joint.vcf.gz', 'somatic/somatic.full.vcf.gz'):
            self.write('results/cohort_run01/gridss/PAIR_001/output/' + name, 'not the input')
        rows, source = discover(self.config, self.root)
        self.assertEqual(len(rows), 10)
        self.assertEqual({row['caller'] for row in rows}, {'manta', 'delly', 'svaba', 'gridss', 'brass'})
        self.assertTrue(all(row['vcf'].endswith('somatic.high_confidence.vcf.gz') for row in rows if row['caller'] == 'gridss'))
        self.assertTrue(all(row['source_sha256'] == '' for row in rows))
        self.assertEqual(source['cohort'], 'synthetic')
        self.assertFalse((self.root / 'inputs/bams').exists())

    def test_compact_publishing_layouts(self):
        self.complete_cohort(compact=True)
        self.assertEqual(len(discover(self.config, self.root)[0]), 10)

    def test_historical_manta_benchmarks_still_select_native_vcfs(self):
        self.caller.update(callers='manta', pair_ids='PAIR_001')
        self.save_caller()
        path = self.complete('manta', legacy_manta=True)
        rows, _ = discover(self.config, self.root)
        self.assertEqual([row['vcf'] for row in rows], [str(path)])

    def test_group_and_pair_selections_do_not_require_unselected_outputs(self):
        self.complete('delly', 'PAIR_002')
        self.config.update(groups='synthetic:delly', pair_ids='PAIR_002')
        rows, _ = discover(self.config, self.root)
        self.assertEqual([(row['caller'], row['pair_id']) for row in rows], [('delly', 'PAIR_002')])
        for overrides in ({'groups': 'wrong:delly'}, {'groups': 'synthetic:delly,synthetic:delly'},
                          {'pair_ids': 'UNKNOWN'}, {'pair_id': 'PAIR_001', 'pair_ids': 'PAIR_002'}):
            with self.subTest(overrides=overrides), self.assertRaises(ValueError):
                discover(dict(self.config, **overrides), self.root)

    def test_caller_pair_selection_and_custom_output_directory_are_respected(self):
        self.caller.update(outdir='results/experiment02', callers='manta', pair_ids='PAIR_001')
        self.save_caller()
        path = self.complete('manta')
        rows, _ = discover(self.config, self.root)
        self.assertEqual([row['vcf'] for row in rows], [str(path)])

    def test_saved_caller_settings_use_the_saved_samplesheet_when_cache_is_removed(self):
        self.caller.update(callers='manta', pair_ids='PAIR_001')
        self.save_caller()
        self.complete('manta')
        saved = dict(self.caller, input='.nextflow/input-manifests/removed.tsv')
        self.write('logs/callers-previous/resolved.yaml', workflow.yaml.safe_dump(saved))
        self.write('logs/callers-previous/resolved-samplesheet.tsv',
                   (self.root / 'inputs/samplesheet.tsv').read_text())
        self.config['caller_config'] = 'logs/callers-previous/resolved.yaml'
        rows, source = discover(self.config, self.root)
        self.assertEqual(len(rows), 1)
        self.assertTrue(source['samplesheet'].endswith('resolved-samplesheet.tsv'))

    def test_build_and_assay_must_match_the_previous_caller_run(self):
        for overrides in ({'genome_build': 'GRCh38'}, {'sequencing_type': 'WES'}):
            self.caller.update(overrides)
            self.save_caller()
            with self.subTest(overrides=overrides), self.assertRaisesRegex(ValueError, 'differs from signature settings'):
                discover(self.config, self.root)

    def test_missing_failed_stub_or_mismatched_tasks_are_not_silently_included(self):
        self.caller.update(callers='manta', pair_ids='PAIR_001')
        self.save_caller()
        vcf = self.complete('manta')
        benchmark = vcf.parents[3] / 'output/benchmark.json'
        original = benchmark.read_text()
        for overrides in ({'exit_code': 1}, {'exit_code': None}, {'exit_code': False},
                          {'schema': 'stub_not_benchmarked'}, {'dataset_id': 'another'},
                          {'pair_id': 'PAIR_OTHER'}, {'caller': 'delly'}, {'reference': '/other/reference.fa'}):
            benchmark.write_text(json.dumps(dict(json.loads(original), **overrides)))
            with self.subTest(overrides=overrides), self.assertRaisesRegex(ValueError, 'not ready'):
                discover(self.config, self.root)
        benchmark.write_text(original)
        vcf.unlink()
        with self.assertRaisesRegex(ValueError, 'final SV VCF'):
            discover(self.config, self.root)

    def test_native_completion_markers_are_required(self):
        for caller, marker in (('delly', 'SUCCESS.json'), ('gridss', 'somatic/SUCCESS.json'),
                               ('brass', 'callset/technical-qc.json')):
            self.caller.update(callers=caller, pair_ids='PAIR_001')
            self.save_caller()
            self.complete(caller)
            (self.root / self.caller['outdir'] / caller / 'PAIR_001/output' / marker).unlink()
            with self.subTest(caller=caller), self.assertRaisesRegex(ValueError, 'not ready'):
                discover(self.config, self.root)

    def test_raw_gridss_configuration_is_not_a_signature_input(self):
        self.caller.update(callers='gridss', gridss={'somatic_filter': False})
        self.save_caller()
        with self.assertRaisesRegex(ValueError, 'somatic_filter=true'):
            discover(self.config, self.root)

    def test_ambiguous_final_outputs_are_reported(self):
        self.caller.update(callers='manta', pair_ids='PAIR_001')
        self.save_caller()
        self.complete('manta')
        self.complete('manta', compact=True)
        with self.assertRaisesRegex(ValueError, 'expected one file, found 2'):
            discover(self.config, self.root)

    def test_check_files_is_read_only_and_does_not_hash_vcfs(self):
        self.complete_cohort()
        before = {str(p): p.stat().st_mtime_ns for p in self.root.rglob('*') if p.is_file()}
        original = Path.open
        def guard(path, *args, **kwargs):
            if '.vcf' in path.name:
                raise AssertionError('Read-only inventory check must not read/hash VCF contents')
            return original(path, *args, **kwargs)
        with patch.object(Path, 'open', guard):
            report = check_files(self.config, 'signatures', self.root)
        self.assertEqual(report['status'], 'passed', report['issues'])
        self.assertEqual(report['selected_rows'], 10)
        self.assertEqual(before, {str(p): p.stat().st_mtime_ns for p in self.root.rglob('*') if p.is_file()})

    def test_auto_manifest_hashes_inputs_and_keeps_resume_path_and_mtime(self):
        self.complete_cohort()
        first_log, next_log = self.root / 'first', self.root / 'second'
        first_log.mkdir(); next_log.mkdir()
        first = write_auto_manifest(self.config, self.root, first_log)
        mtime = first.stat().st_mtime_ns
        second = write_auto_manifest(self.config, self.root, next_log)
        self.assertEqual(first, second)
        self.assertEqual(first.stat().st_mtime_ns, mtime)
        self.assertEqual(first.read_bytes(), (first_log / 'input_manifest.tsv').read_bytes())
        with first.open() as handle:
            rows = list(csv.DictReader(handle, delimiter='\t'))
        for row in rows:
            self.assertEqual(row['source_sha256'], hashlib.sha256(Path(row['vcf']).read_bytes()).hexdigest())
        Path(rows[0]['vcf']).write_bytes(Path(rows[0]['vcf']).read_bytes() + b'changed')
        self.assertNotEqual(write_auto_manifest(self.config, self.root, next_log), first)

    def test_launcher_creates_and_passes_the_manifest_without_manual_inventory(self):
        self.complete_cohort()
        process = MagicMock(); process.stdout = iter(['synthetic completion\n']); process.wait.return_value = 0
        with patch.object(workflow, 'ROOT', self.root), patch.object(sys, 'argv', ['workflow.py', 'signatures']), \
             patch.object(workflow.shutil, 'which', return_value='/synthetic/nextflow'), \
             patch.object(workflow.subprocess, 'check_output', return_value='nextflow version 25.10.7'), \
             patch.object(workflow.subprocess, 'Popen', return_value=process) as submit, redirect_stdout(io.StringIO()):
            with self.assertRaises(SystemExit) as stopped:
                workflow.main()
        self.assertEqual(stopped.exception.code, 0)
        log = next((self.root / 'logs').glob('signatures-*'))
        resolved = workflow.read_yaml(log / 'resolved.yaml')
        self.assertEqual(Path(resolved['manifest']).read_bytes(), (log / 'input_manifest.tsv').read_bytes())
        self.assertEqual(json.loads((log / 'signature-inputs.json').read_text())['callsets'], 10)
        self.assertTrue(submit.called)
        self.assertFalse((self.root / 'inputs/callsets.tsv').exists())

    def test_dry_run_and_missing_output_checks_submit_nothing(self):
        for option in ('--dry-run', '--check-files', None):
            with patch.object(workflow, 'ROOT', self.root), \
                 patch.object(sys, 'argv', ['workflow.py', 'signatures'] + ([option] if option else [])), \
                 patch.object(workflow.subprocess, 'Popen') as submit, \
                 redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()):
                if option == '--dry-run':
                    workflow.main()
                else:
                    with self.assertRaises(SystemExit) as stopped:
                        workflow.main()
                    self.assertEqual(stopped.exception.code, 2)
                submit.assert_not_called()
        self.assertFalse((self.root / 'logs').exists())
        self.assertFalse((self.root / '.nextflow').exists())

    def test_external_manifest_does_not_require_a_caller_configuration(self):
        vcf = self.write('inputs/external.vcf', (ROOT / 'examples/toy.somatic.vcf').read_text())
        manifest = self.write('manifests/external.tsv',
            'cohort\tpair_id\tcaller\tvcf\tsource_status\tsource_sha256\n' +
            f'external\tTOY_PAIR\tmanta\t{vcf}\texternal_reviewed\t{hashlib.sha256(vcf.read_bytes()).hexdigest()}\n')
        (self.root / 'config/private-cohort.yaml').unlink()
        self.config['manifest'] = str(manifest)
        with patch('signature_inputs.discover', side_effect=AssertionError('external route must not discover caller outputs')):
            report = check_files(self.config, 'signatures', self.root)
        self.assertEqual(report['status'], 'passed', report['issues'])
        self.assertEqual(report['selected_rows'], 1)


if __name__ == '__main__':
    unittest.main()
