import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / 'scripts/measure_process_tree.py'
METER = Path(__file__).resolve().parents[1] / 'scripts/measure_command.py'


class ProcessTreeMemoryTests(unittest.TestCase):
    def test_records_process_tree_and_exit_code(self):
        with tempfile.TemporaryDirectory() as temp:
            report = Path(temp) / 'usage.json'
            standard = Path(temp) / 'call.time.txt'
            command = [sys.executable, str(METER), str(standard), sys.executable,
                       str(SCRIPT), str(report), sys.executable,
                       '-c', 'import time; time.sleep(2.2); raise SystemExit(7)']
            result = subprocess.run(command, check=False, env={**os.environ, 'PATH': ''})
            data = json.loads(report.read_text())
            standard_data = json.loads(standard.read_text())
            self.assertEqual(result.returncode, 7)
            self.assertEqual(data['exit_code'], 7)
            self.assertEqual(standard_data['exit_code'], 7)
            self.assertGreater(standard_data['max_rss_kib'], 0)
            self.assertGreater(data['max_concurrent_tree_rss_kib'], 0)
            self.assertGreater(data['max_single_process_rss_kib'], 0)
            self.assertGreaterEqual(data['samples'], 1)


if __name__ == '__main__':
    unittest.main()
