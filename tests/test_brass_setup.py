import csv
from io import BytesIO
import json
from pathlib import Path
import subprocess
import sys
import tarfile
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import setup_brass as setup
from adapt_brass_resources import dictionary
from workflow import read_yaml


def archive(path, members):
    with tarfile.open(path, 'w:gz') as stream:
        for name, content in members:
            member = tarfile.TarInfo(name)
            if content is None:
                member.type = tarfile.SYMTYPE
                member.linkname = '/outside'
                stream.addfile(member)
            else:
                member.size = len(content)
                stream.addfile(member, BytesIO(content))


class BrassSetupTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def test_runtime_detection_and_explicit_choices(self):
        for installed, preferred, explicit, expected in (
            ({'apptainer', 'singularity'}, 'apptainer', False, 'apptainer'),
            ({'singularity'}, 'apptainer', False, 'singularity'),
            ({'apptainer'}, 'apptainer', False, 'apptainer'),
            ({'apptainer', 'singularity'}, 'singularity', True, 'singularity'),
            ({'apptainer', 'singularity'}, 'apptainer', True, 'apptainer'),
        ):
            with self.subTest(installed=installed, explicit=explicit), \
                 patch.object(setup.shutil, 'which', side_effect=lambda x: '/bin/' + x if x in installed else None):
                self.assertEqual(setup.select_engine(preferred, explicit), expected)
        with patch.object(setup.shutil, 'which', side_effect=lambda x: '/bin/singularity' if x == 'singularity' else None):
            with self.assertRaisesRegex(ValueError, 'Requested container runtime is unavailable: apptainer'):
                setup.select_engine('apptainer', explicit=True)
        with patch.object(setup.shutil, 'which', return_value=None):
            with self.assertRaisesRegex(ValueError, 'No container runtime found'):
                setup.select_engine('apptainer')
        with self.assertRaisesRegex(ValueError, 'container_engine must be'):
            setup.select_engine('docker')

    def test_cli_runtime_precedence(self):
        overlay = self.root / 'cohort.yaml'
        overlay.write_text('container_engine: singularity\n')
        paths_only = self.root / 'paths.yaml'
        paths_only.write_text('reference: resources/reference/reference.fa\n')
        cases = [([], 'apptainer', False),
                 (['--config', str(paths_only)], 'apptainer', False),
                 (['--config', str(overlay)], 'singularity', True),
                 (['--config', str(overlay), '--config', str(paths_only)], 'singularity', True),
                 (['--config', str(overlay), '--set', 'container_engine=apptainer'], 'apptainer', True),
                 (['--set', 'container_engine=singularity'], 'singularity', True),
                 (['--set', 'container_engine=apptainer', '--engine', 'singularity'], 'singularity', True)]
        for argv, engine, explicit in cases:
            with self.subTest(argv=argv), patch.object(sys, 'argv', ['setup_brass.py'] + argv), \
                 patch.object(setup, 'setup') as run_setup:
                setup.main()
            self.assertEqual(run_setup.call_args.args[0]['container_engine'], engine)
            self.assertEqual(run_setup.call_args.kwargs['engine_explicit'], explicit)

    def test_auto_detected_singularity_is_recorded_for_setup(self):
        config = self.configured_root()
        (self.root / 'resources/brass/bundle').mkdir(parents=True)
        with patch.object(setup.shutil, 'which', side_effect=lambda x: None if x == 'apptainer' else '/bin/' + x), \
             patch.object(setup, 'check_existing_bundle'), \
             patch.object(setup.subprocess, 'run', side_effect=AssertionError('not needed')):
            result = setup.setup(config, self.root, engine_explicit=False)
        self.assertEqual(result['container_engine'], 'singularity')
        self.assertEqual(config['container_engine'], 'apptainer')  # User YAML/defaults stay unchanged.
        status = next((self.root / 'logs').glob('brass-setup-*/status.json'))
        self.assertEqual(json.loads(status.read_text())['container_engine'], 'singularity')

    def test_versioned_source_dictionary_is_public_metadata(self):
        path = ROOT / 'config/brass-GRCh37d5.dict'
        identities = dictionary(path)
        self.assertEqual(len(identities), 86)
        self.assertEqual(identities['1'], (249250621, '1b22b98cdeb4a9304cb5d48026a85128'))
        self.assertEqual(identities['MT'], (16569, 'c68f52674c9fb33aef52dcf399755519'))
        text = path.read_text()
        self.assertNotIn('file://', text)
        self.assertIn('https://ftp.1000genomes.ebi.ac.uk/', text)

    def test_manifest_pins_only_used_archives(self):
        with (ROOT / 'config/brass-resources.tsv').open() as handle:
            rows = list(csv.DictReader(handle, delimiter='\t'))
        self.assertEqual({row['file'] for row in rows}, set(setup.ARCHIVES))
        for row in rows:
            self.assertRegex(row['md5'], r'^[a-f0-9]{32}$')
            self.assertRegex(row['sha256'], r'^[a-f0-9]{64}$')
            self.assertTrue(row['url'].startswith('https://ftp.sanger.ac.uk/'))
            self.assertTrue(row['url'].endswith('/' + row['file']))

    def test_download_resume_verifies_and_promotes(self):
        target = self.root / 'download.tar.gz'
        partial = target.with_name(target.name + '.partial')
        partial.write_bytes(b'part')
        content = b'complete archive'
        row = {'file': target.name, 'md5': setup.hashlib.md5(content).hexdigest(),
               'sha256': setup.hashlib.sha256(content).hexdigest(),
               'url': 'https://example.org/download.tar.gz'}
        commands = []
        def run(argv):
            commands.append(argv)
            partial.write_bytes(content)
        record = setup.fetch(row, self.root, run)
        self.assertEqual(target.read_bytes(), content)
        self.assertFalse(partial.exists())
        self.assertIn('--continue-at', commands[0])
        self.assertEqual(record['sha256'], setup.digest(target))
        self.assertEqual(setup.fetch(row, self.root, self.fail), record)

    def test_sha256_mismatch_is_not_promoted_even_if_md5_matches(self):
        content = b'fixture'
        row = {'file': 'bad.tar.gz', 'md5': setup.hashlib.md5(content).hexdigest(),
               'sha256': '0' * 64, 'url': 'https://example.org/bad'}
        def run(argv):
            (self.root / 'bad.tar.gz.partial').write_bytes(content)
        with self.assertRaisesRegex(ValueError, 'SHA-256 checksum mismatch'):
            setup.fetch(row, self.root, run)
        self.assertFalse((self.root / 'bad.tar.gz').exists())

    def test_bad_download_is_not_promoted(self):
        row = {'file': 'bad.tar.gz', 'md5': '0' * 32, 'url': 'https://example.org/bad'}
        def run(argv):
            (self.root / 'bad.tar.gz.partial').write_bytes(b'bad')
        with self.assertRaisesRegex(ValueError, 'checksum mismatch'):
            setup.fetch(row, self.root, run)
        self.assertFalse((self.root / 'bad.tar.gz').exists())
        self.assertTrue((self.root / 'bad.tar.gz.partial').exists())

    def test_corrupt_cached_archive_is_not_replaced(self):
        target = self.root / 'bad.tar.gz'
        target.write_bytes(b'bad')
        with self.assertRaisesRegex(ValueError, 'Cached archive checksum mismatch'):
            setup.fetch({'file': target.name, 'md5': '0' * 32}, self.root, self.fail)
        self.assertEqual(target.read_bytes(), b'bad')

    def test_extract_only_required_regular_files(self):
        path = self.root / 'archive.tar.gz'
        archive(path, [('a/keep', b'needed'), ('unrelated/skip', b'not needed')])
        out = self.root / 'out'
        setup.extract_selected(path, out, ('a/keep',))
        self.assertEqual((out / 'a/keep').read_bytes(), b'needed')
        self.assertFalse((out / 'unrelated').exists())

    def test_reject_unsafe_paths_links_duplicates_and_missing_members(self):
        cases = [([('../escape', b'bad')], 'Unsafe archive path'),
                 ([('/absolute', b'bad')], 'Unsafe archive path'),
                 ([('wanted', None)], 'regular archive member'),
                 ([('wanted', b'a'), ('wanted', b'b')], 'regular archive member'),
                 ([('other', b'a')], 'Missing archive members')]
        for i, (members, message) in enumerate(cases):
            with self.subTest(message=message, i=i):
                path = self.root / f'archive{i}.tar.gz'
                archive(path, members)
                with self.assertRaisesRegex(ValueError, message):
                    setup.extract_selected(path, self.root / str(i), ('wanted',))

    def source_fixture(self):
        downloads = self.root / 'downloads'
        downloads.mkdir()
        for filename, members in setup.ARCHIVES.items():
            archive(downloads / filename, [(name, b'fixture') for name in members])
        source_dict = self.root / 'reference.dict'
        source_dict.write_text('@SQ\tSN:1\tLN:10\tM5:abc\n')
        return downloads, source_dict

    def test_source_reuse_checks_file_hashes(self):
        downloads, source_dict = self.source_fixture()
        source = self.root / 'source'
        records = [{'file': name} for name in setup.ARCHIVES]
        setup.prepare_source(downloads, source, records, source_dict)
        setup.prepare_source(downloads, source, records, source_dict)
        changed = source / setup.BRASS / 'CentTelo.tsv'
        changed.write_text('changed')
        with self.assertRaisesRegex(ValueError, 'Extracted source changed'):
            setup.prepare_source(downloads, source, records, source_dict)

    def test_incomplete_extraction_does_not_publish_source(self):
        downloads, source_dict = self.source_fixture()
        archive(downloads / next(iter(setup.ARCHIVES)), [('unrelated', b'no')])
        source = self.root / 'source'
        with self.assertRaisesRegex(ValueError, 'Missing archive members'):
            setup.prepare_source(downloads, source, [], source_dict)
        self.assertFalse(source.exists())
        self.assertEqual(list(self.root.glob('.brass-source-*')), [])

    def configured_root(self):
        (self.root / 'config').mkdir()
        (self.root / 'config/brass-GRCh37d5.dict').write_text('@SQ\tSN:1\tLN:10\tM5:abc\n')
        reference = self.root / 'resources/reference/reference.fa'
        reference.parent.mkdir(parents=True)
        reference.write_text('>chr1\nACGT\n')
        Path(str(reference) + '.fai').write_text('chr1\t10\n')
        reference.with_suffix('.dict').write_text('@SQ\tSN:chr1\tLN:10\tM5:abc\n')
        image = self.root / 'containers/images/brass-v6.3.4.sif'
        image.parent.mkdir(parents=True)
        image.write_bytes(b'fake image')
        return read_yaml(ROOT / 'config/defaults.yaml')

    def test_incompatible_reference_fails_before_tools_or_downloads(self):
        config = self.configured_root()
        path = self.root / 'resources/reference/reference.dict'
        path.write_text('@SQ\tSN:chr1\tLN:10\tM5:wrong\n')
        with patch.object(setup.shutil, 'which', side_effect=AssertionError('too late')):
            with self.assertRaisesRegex(ValueError, 'Sequence length/MD5 mismatch'):
                setup.setup(config, self.root)
        self.assertFalse((self.root / 'resources/brass').exists())

    def test_grch38_rejected_before_download(self):
        config = read_yaml(ROOT / 'config/defaults.yaml')
        config['genome_build'] = 'GRCh38'
        with self.assertRaisesRegex(ValueError, 'bundle is for GRCh37'):
            setup.setup(config, self.root)

    def test_completed_bundle_reused_without_download_or_preparation(self):
        config = self.configured_root()
        bundle = self.root / 'resources/brass/bundle'
        bundle.mkdir(parents=True)
        with patch.object(setup.shutil, 'which', return_value='/tool'), \
             patch.object(setup, 'check_existing_bundle') as check, \
             patch.object(setup.subprocess, 'run', side_effect=AssertionError('not needed')):
            result = setup.setup(config, self.root)
        self.assertEqual(result['outcome'], 'reused_verified_bundle')
        check.assert_called_once()
        status = next((self.root / 'logs').glob('brass-setup-*/status.json'))
        self.assertEqual(json.loads(status.read_text())['status'], 'passed')

    def complete_bundle(self):
        self.configured_root()
        reference = self.root / 'resources/reference/reference.fa'
        image = self.root / 'containers/images/brass-v6.3.4.sif'
        bundle = self.root / 'resources/brass/bundle'
        bundle.mkdir(parents=True)
        for name in setup.BRASS_FILES:
            (bundle / name).write_text('fixture')
        (bundle / 'COMPATIBILITY_VALIDATED').touch()
        (bundle / 'native-validation.json').write_text('{"status":"passed"}')
        (bundle / 'adaptation.json').write_text(json.dumps({
            'status': 'adapted_requires_native_reader_validation', 'coordinate_changes': 0,
            'target_dictionary_sha256': setup.digest(reference.with_suffix('.dict')),
            'resources': {'HiDepth.bed.gz': {'adapted_sha256': setup.digest(bundle / 'HiDepth.bed.gz')}}}))
        inputs = {'reference_sha256': reference, 'reference_fai_sha256': Path(str(reference) + '.fai'),
                  'target_dictionary_sha256': reference.with_suffix('.dict'), 'image_sha256': image,
                  'overlay_sha256': bundle / 'get_rg_cns.R'}
        receipt = dict(status='passed', assembly='GRCh37', **{
            key: setup.digest(path) for key, path in inputs.items()})
        (bundle / 'bundle-receipt.json').write_text(json.dumps(receipt))
        return bundle, reference, image

    def test_bundle_reuse_verifies_hashes_and_native_status(self):
        bundle, reference, image = self.complete_bundle()
        setup.check_existing_bundle(bundle, reference, image)
        (bundle / 'HiDepth.bed.gz').write_text('modified')
        with self.assertRaisesRegex(ValueError, 'resource changed'):
            setup.check_existing_bundle(bundle, reference, image)
        (bundle / 'native-validation.json').write_text('{"status":"failed"}')
        with self.assertRaisesRegex(ValueError, 'validation did not pass'):
            setup.check_existing_bundle(bundle, reference, image)

    def test_bundle_reuse_rejects_changed_shared_reference(self):
        bundle, reference, image = self.complete_bundle()
        reference.write_text('changed')
        with self.assertRaisesRegex(ValueError, 'reference_sha256'):
            setup.check_existing_bundle(bundle, reference, image)

    def test_failed_preparation_is_logged_without_publishing_bundle(self):
        config = self.configured_root()
        with (ROOT / 'config/brass-resources.tsv').open() as handle:
            rows = list(csv.DictReader(handle, delimiter='\t'))
        (self.root / 'config/brass-resources.tsv').write_text(
            (ROOT / 'config/brass-resources.tsv').read_text())
        with patch.object(setup.shutil, 'which', return_value='/tool'), \
             patch.object(setup, 'fetch', side_effect=rows), \
             patch.object(setup, 'prepare_source'), \
             patch.object(setup.subprocess, 'run', side_effect=subprocess.CalledProcessError(1, 'prepare')):
            with self.assertRaisesRegex(ValueError, 'adaptation_and_native_validation'):
                setup.setup(config, self.root)
        self.assertFalse((self.root / 'resources/brass/bundle').exists())
        status = next((self.root / 'logs').glob('brass-setup-*/status.json'))
        self.assertEqual(json.loads(status.read_text())['step'], 'adaptation_and_native_validation')

    def test_successful_preparation_publishes_bundle_atomically(self):
        config = self.configured_root()
        (self.root / 'config/brass-resources.tsv').write_text(
            (ROOT / 'config/brass-resources.tsv').read_text())
        final = self.root / 'resources/brass/bundle'
        def run(argv, **kwargs):
            self.assertFalse(final.exists())
            out = Path(argv[argv.index('--output') + 1])
            self.assertNotEqual(out, final)
            out.mkdir()
            (out / 'test-complete').touch()
        with patch.object(setup.shutil, 'which', return_value='/tool'), \
             patch.object(setup, 'fetch', return_value={}), \
             patch.object(setup, 'prepare_source'), \
             patch.object(setup, 'check_existing_bundle'), \
             patch.object(setup.subprocess, 'run', side_effect=run):
            result = setup.setup(config, self.root)
        self.assertEqual(result['outcome'], 'prepared_and_validated')
        self.assertTrue((final / 'test-complete').is_file())
        self.assertEqual(list(final.parent.glob('.brass-bundle-*')), [])


if __name__ == '__main__':
    unittest.main()
