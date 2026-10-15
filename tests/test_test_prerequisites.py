"""Installation failures are reported before synthetic tasks or output creation."""
from contextlib import redirect_stderr, redirect_stdout
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import test_main_stub as stub
import workflow


class TestPrerequisiteTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def test_stub_lists_missing_tools_without_creating_output(self):
        out = self.root / 'check'
        with patch.object(sys, 'argv', ['test_main_stub.py', '--out', str(out)]), \
             patch.object(stub.shutil, 'which', return_value=None), \
             patch.object(stub.subprocess, 'run') as task, redirect_stderr(io.StringIO()) as error:
            with self.assertRaises(SystemExit) as stopped:
                stub.main()
        self.assertEqual(stopped.exception.code, 2)
        for name in ('nextflow', 'java', 'samtools'):
            self.assertIn(name, error.getvalue())
        self.assertIn('docs/setup.md', error.getvalue())
        self.assertFalse(out.exists())
        task.assert_not_called()

    def test_stub_preserves_existing_output_and_suggests_a_fresh_path(self):
        out = self.root / 'check'
        out.mkdir()
        marker = out / 'keep.txt'; marker.write_text('previous attempt')
        with patch.object(sys, 'argv', ['test_main_stub.py', '--out', str(out)]), \
             patch.object(stub.shutil, 'which', return_value='/tool'), \
             patch.object(stub.subprocess, 'check_output') as version, redirect_stderr(io.StringIO()) as error:
            with self.assertRaises(SystemExit):
                stub.main()
        self.assertIn('choose a fresh --out', error.getvalue())
        self.assertEqual(marker.read_text(), 'previous attempt')
        version.assert_not_called()

    def test_stub_reports_nextflow_startup_or_version_errors_before_creating_output(self):
        for value in ('nextflow version 24.10.0', subprocess.CalledProcessError(1, ['nextflow'], output='Java 8 is unsupported')):
            out = self.root / 'check'
            kwargs = {'side_effect': value} if isinstance(value, Exception) else {'return_value': value}
            with self.subTest(value=str(value)), \
                 patch.object(sys, 'argv', ['test_main_stub.py', '--out', str(out)]), \
                 patch.object(stub.shutil, 'which', return_value='/tool'), \
                 patch.object(stub.subprocess, 'check_output', **kwargs), redirect_stderr(io.StringIO()) as error:
                with self.assertRaises(SystemExit):
                    stub.main()
            self.assertIn('25.10.7', error.getvalue())
            self.assertFalse(out.exists())

    def test_stub_checks_optional_manta_region_tools(self):
        image = self.root / 'manta.sif'; image.write_text('synthetic path only')
        out = self.root / 'check'
        with patch.object(sys, 'argv', ['test_main_stub.py', '--out', str(out), '--manta-image', str(image),
                                       '--container-engine', 'singularity', '--with-exclusions']), \
             patch.object(stub.shutil, 'which', side_effect=lambda name: None if name in ('singularity','bgzip','tabix') else '/tool'), \
             redirect_stderr(io.StringIO()) as error:
            with self.assertRaises(SystemExit):
                stub.main()
        for name in ('singularity', 'bgzip', 'tabix'):
            self.assertIn(name, error.getvalue())
        self.assertFalse(out.exists())

    def test_launcher_records_java_nextflow_startup_failure_without_submitting_tasks(self):
        for name in ('config/test.yaml', 'main.nf', 'nextflow.config'):
            path = self.root / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text((ROOT / name).read_text())
        with patch.object(workflow, 'ROOT', self.root), \
             patch.object(sys, 'argv', ['workflow.py', 'test', '--profile', 'test']), \
             patch.object(workflow.shutil, 'which', return_value='/tool/nextflow'), \
             patch.object(workflow.subprocess, 'check_output', side_effect=subprocess.CalledProcessError(
                 1, ['nextflow', '-version'], output='Java 8 is unsupported')), \
             patch.object(workflow.subprocess, 'Popen') as task, redirect_stderr(io.StringIO()) as error:
            with self.assertRaises(SystemExit) as stopped:
                workflow.main()
        self.assertEqual(stopped.exception.code, 2)
        self.assertIn('Java 17+', error.getvalue())
        task.assert_not_called()
        log = next((self.root / 'logs').glob('test-*'))
        self.assertEqual((log / 'startup.log').read_text(), 'Java 8 is unsupported')
        self.assertEqual(json.loads((log / 'status.json').read_text())['stage'], 'nextflow_startup')


if __name__ == '__main__':
    unittest.main()
