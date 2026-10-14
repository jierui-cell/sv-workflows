import sys
from pathlib import Path
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from filter_gridss import validate, patch_r_compatibility


class GridssFilterTests(unittest.TestCase):
    def test_defaults_and_override(self):
        self.assertEqual(validate({})['min_qual'],350)
        self.assertEqual(validate({'min_normal_depth':12})['min_normal_depth'],12)

    def test_reject_invalid(self):
        for options in [{'typo':1},{'min_qual':-1},{'min_af':2},{'min_qual':float('nan')},{'min_qual':True}]:
            with self.assertRaises(ValueError): validate(options)

    def test_guarded_compatibility_patch(self):
        src='  isbp = str_detect(VariantAnnotation::fixed(vcf)$ALT, "pattern")'
        self.assertIn('elementNROWS',patch_r_compatibility(src))
        self.assertIn('as.character(unlist',patch_r_compatibility(src))
        with self.assertRaises(ValueError): patch_r_compatibility('unexpected')

    def test_pipeline_configuration(self):
        from run_caller import validate as pipeline_validate
        cfg={'cpus':8,'memory_gb':64,'jvmheap_gb':40,'otherjvmheap_gb':4,
             'somatic_filter':True,'somatic_filter_options':{'min_normal_depth':12}}
        pipeline_validate('gridss',cfg)
        with self.assertRaises(ValueError): pipeline_validate('gridss',dict(cfg,somatic_filter='yes'))
