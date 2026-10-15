import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from summarize_benchmarks import CALLERS, paths_for, summarize


ROOT = Path(__file__).resolve().parents[1]


class BenchmarkReportingTests(unittest.TestCase):
    @unittest.skipUnless(Path('/proc/self/stat').is_file(),'Linux /proc memory recorder')
    def test_recorder_measures_child_and_preserves_failure(self):
        with tempfile.TemporaryDirectory() as temp:
            base = Path(temp)
            meta = base / 'meta.json'
            report = base / 'output/benchmark.json'
            meta.write_text(json.dumps(dict(caller='test', pair_id='PAIR1', dataset_id='test',
                                            reference='ref.fa', allocated_cpus=1,
                                            requested_memory_gb=1, settings={})))
            command = [sys.executable, str(ROOT / 'scripts/measure_benchmark.py'),
                       str(report), str(meta), '--', sys.executable, '-c',
                       'import time; data=bytearray(20_000_000); time.sleep(1.2); raise SystemExit(7)']
            result = subprocess.run(command, check=False)
            data = json.loads(report.read_text())
            self.assertEqual(result.returncode, 7)
            self.assertEqual(data['schema'], 'sv-caller-benchmark/v1')
            self.assertEqual(data['exit_code'], 7)
            self.assertEqual(data['memory_status'], 'complete')
            self.assertGreater(data['peak_sampled_resident_kib'], 10000)
            self.assertGreater(data['peak_sampled_pss_kib'], 10000)

    def test_balanced_summary_never_substitutes_missing_caller(self):
        with tempfile.TemporaryDirectory() as temp:
            base = Path(temp)
            root, out = base / 'calls', base / 'summary'
            for caller in ('manta', 'delly'):
                p = root / caller / 'PAIR1/output'
                p.mkdir(parents=True)
                (p / 'benchmark.json').write_text(json.dumps(dict(
                    schema='sv-caller-benchmark/v1', caller=caller, pair_id='PAIR1',
                    dataset_id='DATA1', reference='ref.fa', boundary='common boundary',
                    output_class='somatic', settings_sha256='same-config',
                    exit_code=0, memory_status='complete',
                    elapsed_seconds=3600, peak_sampled_pss_kib=1048576,
                    peak_sampled_resident_kib=1048576,
                    memory_method_id='cgroup_v1_total_rss_plus_total_mapped_file',
                    peak_sampled_rss_sum_kib=2000000, allocated_cpus=2,
                    requested_memory_gb=8, start_utc='2026-09-23T12:00:00+00:00',
                    end_utc='2026-09-23T13:00:00+00:00')))
            ready = summarize(root, out, ('manta', 'delly'))
            self.assertEqual(ready['status'], 'ready')
            self.assertEqual(ready['balanced_pairs'], ['PAIR1'])
            with (out / 'balanced_runs.tsv').open() as handle:
                self.assertEqual(len(list(csv.DictReader(handle, delimiter='\t'))), 2)
            partial = summarize(root, out, ('manta', 'delly', 'svaba'), allow_incomplete=True)
            self.assertEqual(partial['status'], 'partial')
            self.assertEqual(partial['balanced_pairs'], [])
            strict = summarize(root, out, ('manta', 'delly', 'svaba'))
            self.assertEqual(strict['status'], 'review_required')
            duplicate = summarize([root, root], out, ('manta', 'delly'))
            self.assertEqual(duplicate['status'], 'review_required')

            # Do not silently select one of two benchmark files for a pair.
            current = root / 'manta/PAIR1/output/benchmark.json'
            (current.parent.parent / 'benchmark.json').write_text(current.read_text())
            conflicting = summarize(root, out, ('manta', 'delly'))
            self.assertEqual(conflicting['status'], 'review_required')
            self.assertTrue(any('Duplicate benchmark' in issue for issue in conflicting['issues']))

    def test_discovers_standard_historical_and_compact_benchmark_paths(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            for caller in CALLERS:
                for pair, relative in (('CURRENT', 'output/benchmark.json'),
                                       ('LEGACY_OR_COMPACT', 'benchmark.json')):
                    path = root / caller / pair / relative
                    path.parent.mkdir(parents=True)
                    path.write_text('{}')
                found = dict(paths_for(root, caller))
                self.assertEqual(set(found), {'CURRENT', 'LEGACY_OR_COMPACT'})
                self.assertEqual(found['CURRENT'], root / caller / 'CURRENT/output/benchmark.json')
                self.assertEqual(found['LEGACY_OR_COMPACT'], root / caller / 'LEGACY_OR_COMPACT/benchmark.json')

    def test_manta_declares_and_writes_the_standard_benchmark_path(self):
        process = (ROOT / 'main.nf').read_text().split('process MANTA {', 1)[1].split('\nworkflow {', 1)[0]
        self.assertIn("path 'output/benchmark.json'", process)
        self.assertIn("python3 '${benchmarker}' output/benchmark.json benchmark.meta.json", process)
        self.assertIn("mkdir -p output\n    echo '{\"status\":\"stub_not_benchmarked\"}' > output/benchmark.json", process)


if __name__ == '__main__':
    unittest.main()
