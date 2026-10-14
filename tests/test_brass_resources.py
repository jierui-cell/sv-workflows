import sys
import unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from adapt_brass_resources import mapping, change_cache, annotation_only


class BrassResourceTests(unittest.TestCase):
    def test_annotation_only_classes(self):
        for name in ['HG1032_PATCH', 'HG142_HG150_NOVEL_TEST', 'HSCHR6_MHC_APD']:
            self.assertTrue(annotation_only(name))
        for name in ['1', 'GL000191.1', 'MT', 'unknown', 'chr1']:
            self.assertFalse(annotation_only(name))

    def test_sequence_verified_mapping(self):
        self.assertEqual(mapping({'1': (10, 'abc'), 'hs37d5': (20, 'def')},
                                 {'chr1': (10, 'abc')}), {'1': 'chr1'})

    def test_reject_sequence_mismatch_and_missing_contig(self):
        for source, target in [({'1': (10, 'abc')}, {'chr1': (10, 'def')}),
                               ({'1': (10, 'abc')}, {'chr1': (10, 'abc'), 'chr2': (20, 'def')})]:
            with self.assertRaises(ValueError):
                mapping(source, target)

    def test_only_chromosome_fields_change(self):
        original = ['1', '10', '20', 'ENST1', 'GENE1', '10',
                    "{'_chr' => '1', '_acc' => 'ENST1', '_minpos' => 11}"]
        changed = change_cache(original[:], 'chr1')
        self.assertEqual(changed[1:6], original[1:6])
        self.assertEqual(changed[6], "{'_chr' => 'chr1', '_acc' => 'ENST1', '_minpos' => 11}")
        with self.assertRaises(ValueError):
            change_cache(['2'] + original[1:], 'chr2')


if __name__ == '__main__':
    unittest.main()
