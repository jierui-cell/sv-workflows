import contextlib,io,sys,unittest
from pathlib import Path
from types import SimpleNamespace
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from cue_read_policy import compatible_reads,patch_io

class CueReadPolicyTest(unittest.TestCase):
    def test_narrow_exclusion_and_log(self):
        mapped=SimpleNamespace(query_sequence='ACT',is_unmapped=False)
        unmapped=SimpleNamespace(query_sequence='ACT',is_unmapped=True)
        missing=SimpleNamespace(query_sequence=None,is_unmapped=True)
        log=io.StringIO()
        with contextlib.redirect_stdout(log):
            self.assertEqual(list(compatible_reads([mapped,missing,unmapped])),[mapped,unmapped])
        self.assertIn('"skipped_unmapped_missing_sequence": 1',log.getvalue())

    def test_mapped_missing_fails(self):
        with contextlib.redirect_stdout(io.StringIO()):
            with self.assertRaisesRegex(ValueError,'mapped alignment'):
                list(compatible_reads([SimpleNamespace(query_sequence=None,is_unmapped=False)]))

    def test_patch_guard_and_iterators(self):
        source='for i, read in enumerate(input_bam.fetch(chr_name)):\n    pass\nfor read in input_bam.fetch(chr_name, start, end):\n    pass\n'
        patched=patch_io(source);compile(patched,'adapted_io','exec')
        self.assertEqual(patched.count('compatible_reads(input_bam.fetch'),2)
        with self.assertRaises(ValueError):patch_io(patched)
