import importlib.util
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]


def load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / 'scripts' / (name + '.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class SampleNamePatchTests(unittest.TestCase):
    def test_lumpy_wrapper_owned_cleanup(self):
        patch = load('patch_lumpy_sample_names').patch
        source = ('id:${DISC_SAMPLE}, id:${SPL_SAMPLE}, ' * 2) + '\ntrap cleanup EXIT\n'
        self.assertIn('trap cleanup EXIT', patch(source))
        self.assertNotIn('trap cleanup EXIT', patch(source, keep_temp=True))
        with self.assertRaises(ValueError):
            patch(source.replace('trap cleanup EXIT', ''), keep_temp=True)

    def test_brass_changes_lookup_only(self):
        patch = load('patch_brass_sample_name').patch
        source = "$options->{'safe_tumour_name'},\n                                      $groups,"
        unrelated = "filename=$options->{'safe_tumour_name'};"
        self.assertEqual(patch(source + unrelated), source.replace('safe_tumour_name', 'tumour_name') + unrelated)
        with self.assertRaises(ValueError):
            patch(source + source)
        with self.assertRaises(ValueError):
            patch('unexpected source')

    def test_lumpy_patches_all_four_option_ids_only(self):
        patch = load('patch_lumpy_sample_names').patch
        source = 'id:${DISC_SAMPLE},mean:100 id:${SPL_SAMPLE},weight:1 ' * 2
        self.assertEqual(patch(source), 'id:LUMPY_SAMPLE_${i},mean:100 id:LUMPY_SAMPLE_${i},weight:1 ' * 2)
        with self.assertRaises(ValueError):
            patch(source + 'id:${DISC_SAMPLE},')
        with self.assertRaises(ValueError):
            patch('unexpected source')


if __name__ == '__main__':
    unittest.main()
