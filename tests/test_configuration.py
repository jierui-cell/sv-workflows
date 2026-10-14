import configparser
import gzip
import sys
import os
import json
import subprocess
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from configure_manta import merge
from summarize_manta import summarize
from run_caller import validate, prepare_reference_cache
from prepare_normal_control import relabel


class ConfigurationTests(unittest.TestCase):
    def test_measurement_without_gnu_time_preserves_exit(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp)/'usage.txt'
            script = Path(__file__).resolve().parents[1]/'scripts/measure_command.py'
            env = dict(os.environ, PATH='')
            result = subprocess.run([sys.executable, str(script), str(output),
                sys.executable, '-c', 'raise SystemExit(7)'], env=env)
            self.assertEqual(result.returncode, 7)
            report = json.loads(output.read_text())
            self.assertEqual(report['exit_code'], 7)
            self.assertGreater(report['max_rss_kib'], 0)

    def test_control_reheader_changes_only_sm(self):
        original = '@HD\tSO:coordinate\n@SQ\tSN:chr1\tLN:100\n@RG\tID:rg1\tSM:NORMAL\tLB:lib\n'
        changed = relabel(original, 'CONTROL_T')
        self.assertEqual(changed, original.replace('SM:NORMAL', 'SM:CONTROL_T'))
        with self.assertRaises(ValueError):
            relabel(original, 'NORMAL')

    def test_reference_cache_reuse_and_conflicts(self):
        with tempfile.TemporaryDirectory() as tmp:
            ref = Path(tmp)/'reference.fa'
            for suffix in ['', '.fai', '.amb', '.ann', '.bwt', '.pac', '.sa']:
                Path(str(ref)+suffix).write_text('fixture')
            ref.with_suffix('.dict').write_text('fixture')
            cache = Path(tmp)/'cache'
            prepare_reference_cache(cache, ref)
            prepare_reference_cache(cache, ref)
            self.assertEqual((cache/'ref.fa').resolve(), ref.resolve())
            (cache/'ref.fa').unlink()
            (cache/'ref.fa').symlink_to(ref.with_suffix('.dict'))
            with self.assertRaisesRegex(ValueError, 'Conflicting'):
                prepare_reference_cache(cache, ref)
            (cache/'ref.fa').unlink()
            (cache/'ref.fa').write_text('do not overwrite')
            with self.assertRaisesRegex(ValueError, 'Refusing to overwrite'):
                prepare_reference_cache(cache, ref)
            self.assertEqual((cache/'ref.fa').read_text(), 'do not overwrite')

    def test_caller_defaults_and_overrides(self):
        options, filters = validate('delly', {'cpus': 8, 'options': {'map_qual': 2}})
        self.assertEqual(options, {'map_qual': 2, 'minclip': 25, 'mad_cutoff': 9})
        self.assertEqual(filters['quality'], 300)
        options, _ = validate('svaba', {'cpus': 8})
        self.assertEqual(options['chunk_size'], 25000)
        options, _ = validate('gridss', {'cpus': 8, 'memory_gb': 96,
            'jvmheap_gb': 40, 'otherjvmheap_gb': 4,
            'picard_validation_stringency': 'LENIENT'})
        self.assertEqual(options, {})

    def test_caller_rejects_invalid_settings(self):
        for cfg in ({'cpus': 8, 'options': {'typo': 1}},
                    {'cpus': 8, 'options': {'map_qual': True}},
                    {'cpus': 8, 'filter_options': {'altaf': 2}},
                    {'cpus': 0}):
            with self.assertRaises(ValueError):
                validate('delly', cfg)
        with self.assertRaises(ValueError):
            validate('gridss', {'cpus': 8, 'memory_gb': 40,
                               'jvmheap_gb': 40, 'otherjvmheap_gb': 4})
        with self.assertRaises(ValueError):
            validate('gridss', {'cpus': 8, 'memory_gb': 96,
                'jvmheap_gb': 40, 'otherjvmheap_gb': 4,
                'picard_validation_stringency': 'SILENT'})

    def test_defaults_override_and_reference(self):
        with tempfile.TemporaryDirectory() as tmp:
            base, out = Path(tmp) / 'base.ini', Path(tmp) / 'out.ini'
            base.write_text('[manta]\nminSomaticScore=10\nminPassSomaticScore=30\n')
            merge(str(base), {'minPassSomaticScore': 40}, str(out), 'ref.fa')
            cfg = configparser.ConfigParser()
            cfg.read(str(out))
            self.assertEqual(cfg.getint('manta', 'minPassSomaticScore'), 40)
            self.assertEqual(cfg.getint('manta', 'minSomaticScore'), 10)
            self.assertEqual(cfg.get('manta', 'referenceFasta'), 'ref.fa')

    def test_reject_bad_options(self):
        with tempfile.TemporaryDirectory() as tmp:
            base, out = Path(tmp) / 'base.ini', Path(tmp) / 'out.ini'
            base.write_text('[manta]\nminSomaticScore=10\nflag=0\n')
            for overrides in ({'typo': 3}, {'referenceFasta': 'x'},
                              {'minSomaticScore': -1}, {'minSomaticScore': 1.5},
                              {'flag': 2}, {'flag': True}):
                with self.assertRaises(ValueError):
                    merge(str(base), overrides, str(out))

    def test_mates_count_once_and_filters_excluded(self):
        with tempfile.TemporaryDirectory() as tmp:
            vcf = Path(tmp) / 'input.vcf.gz'
            with gzip.open(str(vcf), 'wt') as out:
                out.write('chr1\t100\ta\tA\tA]chr2:200]\t.\tPASS\tSVTYPE=BND;MATEID=b\n')
                out.write('chr2\t200\tb\tA\tA]chr1:100]\t.\tPASS\tSVTYPE=BND;MATEID=a\n')
                out.write('chr1\t300\tc\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=400;SVLEN=-100\n')
                out.write('chr1\t500\td\tA\t<DEL>\t.\tLowQual\tSVTYPE=DEL;END=600\n')
            summary = summarize(vcf)
            self.assertEqual(summary['pass_records'], 3)
            self.assertEqual(summary['pass_adjacencies_or_symbolic_events'], 2)
            self.assertEqual(summary['mate_warnings'], [])
            self.assertEqual(summary['counts']['BND_interchromosomal'], 1)


if __name__ == '__main__':
    unittest.main()
