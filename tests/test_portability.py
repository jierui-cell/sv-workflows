import csv
import importlib.util
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from workflow import merge, override, read_yaml, validate
from build_sv_manifest import build

class PortabilityTests(unittest.TestCase):
    def test_copyable_caller_template(self):
        cfg=merge(read_yaml(ROOT/'config/defaults.yaml'),read_yaml(ROOT/'examples/cohort.yaml'))
        validate(cfg,'callers')
        self.assertEqual(cfg['callers'],'manta,delly')
        self.assertEqual(cfg['manta_options']['minPassSomaticScore'],30)
        self.assertEqual(cfg['brass']['missing_statistics_policy'],'native_defaults')

    def test_copyable_signature_template(self):
        cfg=merge(read_yaml(ROOT/'config/signatures.yaml'),read_yaml(ROOT/'examples/signatures.yaml'))
        validate(cfg,'signatures')
        cfg['method']='pcawg_clustersv'
        validate(cfg,'signatures')
        self.assertEqual(cfg['pcawg_stage'],'clustering')

    def test_copyable_timing_template(self):
        cfg=merge(read_yaml(ROOT/'config/replication.yaml'),read_yaml(ROOT/'examples/timing.yaml'))
        self.assertEqual(cfg['timing_label_policy'],'tertiles')
        self.assertEqual(cfg['timing_bin_bp'],1000)
        self.assertEqual(cfg['timing_min_coverage'],1.0)

    def test_documented_toy_output(self):
        sys.path.insert(0,str(ROOT/'sv_signatures'))
        from vcf_to_bedpe import convert
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            n=convert(ROOT/'examples/toy.somatic.vcf','manta','TOY_PAIR',root/'bedpe.tsv',root/'audit.tsv',root/'qc.json')
            self.assertEqual(n,3)
            with (root/'bedpe.tsv').open() as actual,(ROOT/'examples/toy_expected/toy.bedpe.tsv').open() as expected:
                self.assertEqual(list(csv.DictReader(actual,delimiter='\t')),list(csv.DictReader(expected,delimiter='\t')))

    def test_nested_override_preserves_defaults(self):
        base=read_yaml(ROOT/'config/defaults.yaml')
        cfg=merge(base,{'delly':{'options':{'map_qual':2}}})
        self.assertEqual(cfg['delly']['options']['minclip'],25)
        self.assertEqual(cfg['delly']['options']['map_qual'],2)
        override(cfg,'max_parallel_callers=8')
        self.assertEqual(cfg['max_parallel_callers'],8)
        self.assertEqual(base['max_parallel_callers'],3)

    def test_typos_and_duplicate_yaml_rejected(self):
        with self.assertRaisesRegex(ValueError,'Unknown'):
            merge({'delly':{'options':{'map_qual':1}}},{'delly':{'options':{'mapqual':2}}})
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'bad.yaml'; p.write_text('callers: manta\ncallers: delly\n')
            with self.assertRaisesRegex(ValueError,'Duplicate'):
                read_yaml(p)

    def test_missing_cn_and_retired_caller_fail_closed(self):
        cfg=read_yaml(ROOT/'config/signatures.yaml')
        cfg.update(method='pcawg_clustersv',pcawg_stage='features')
        with self.assertRaisesRegex(ValueError,'44 features'):
            validate(cfg,'signatures')
        cfg=read_yaml(ROOT/'config/defaults.yaml'); cfg['callers']='lumpy'
        with self.assertRaisesRegex(ValueError,'discontinued'):
            validate(cfg,'callers')

    def test_manifest_hash_and_duplicate_guard(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); inv=root/'inventory.tsv'; out=root/'manifest.tsv'
            row='toy\tTOY_PAIR\tmanta\t'+str(ROOT/'examples/toy.somatic.vcf')+'\tsynthetic_control\n'
            inv.write_text('cohort\tpair_id\tcaller\tvcf\tsource_status\n'+row)
            self.assertEqual(build(inv,out),1)
            with out.open() as handle:
                self.assertEqual(len(next(csv.DictReader(handle,delimiter='\t'))['source_sha256']),64)
            inv.write_text(inv.read_text()+row)
            with self.assertRaisesRegex(ValueError,'Duplicate'):
                build(inv,root/'duplicate.tsv')

    def test_scheduler_dry_runs_submit_nothing(self):
        for profile in ('slurm','pbs','lsf','sge','local'):
            result=subprocess.run([sys.executable,str(ROOT/'scripts/workflow.py'),'callers',
                '--dry-run','--profile',profile,'--allow-local'],capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stderr)
            self.assertIn('-profile '+profile,result.stdout)

if __name__=='__main__': unittest.main()
