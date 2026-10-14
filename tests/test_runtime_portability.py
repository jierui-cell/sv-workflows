"""Regression controls for launchers, recording failures and reference policies."""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
import measure_benchmark as benchmark
from filter_gridss import completed_exit_code
from reference_policy import LENGTHS, check_dictionary, check_vcf
from regions import complement, prepare, read_exclusions
from workflow import read_yaml, validate
from prepare_brass_bundle import validate_adapted

class RuntimePortabilityTests(unittest.TestCase):
    def test_reused_brass_resources_require_dictionary_and_content_identity(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); dictionary=root/'ref.dict'; dictionary.write_text('verified dictionary')
            resource=root/'resource.bed'; resource.write_text('verified resource')
            fasta=root/'resource.fa'; fasta.write_text('unchanged public sequence')
            report={'status':'adapted_requires_native_reader_validation','coordinate_changes':0,
                    'target_dictionary_sha256':hashlib.sha256(dictionary.read_bytes()).hexdigest(),
                    'resources':{'resource.bed':{'adapted_sha256':hashlib.sha256(resource.read_bytes()).hexdigest()},
                                 'resource.fa':{'action':'unchanged_symlink','sha256':hashlib.sha256(fasta.read_bytes()).hexdigest()}}}
            (root/'adaptation.json').write_text(json.dumps(report))
            self.assertEqual(validate_adapted(root,dictionary),report)
            resource.write_text('changed')
            with self.assertRaisesRegex(ValueError,'changed'):
                validate_adapted(root,dictionary)
            dictionary.write_text('different dictionary')
            with self.assertRaisesRegex(ValueError,'another reference'):
                validate_adapted(root,dictionary)
    def test_spooled_slurm_launcher_uses_submit_directory(self):
        with tempfile.TemporaryDirectory() as tmp:
            spool = Path(tmp)/'slurm_script'
            spool.write_text((ROOT/'scripts/run.slurm.sh').read_text())
            env = dict(os.environ, SLURM_JOB_ID='synthetic_launcher_test',
                       SLURM_SUBMIT_DIR=str(ROOT), SV_PYTHON=sys.executable)
            env.pop('SV_WORKFLOW_ROOT', None)
            result = subprocess.run(['bash',str(spool),'callers','--dry-run'],
                                    env=env, capture_output=True, text=True)
            self.assertEqual(result.returncode,0,result.stderr)
            self.assertIn(str(ROOT/'main.nf'),result.stdout)

    def test_null_module_directives_are_safe(self):
        for path in [ROOT/'main.nf',ROOT/'modules/call_sv.nf',ROOT/'modules/integrated_caller.nf']:
            for line in path.read_text().splitlines():
                if 'module { params.modules_' in line:
                    self.assertIn("?: ''",line,str(path))

    def test_reference_dictionary_rejects_other_assembly(self):
        names = [str(i) for i in range(1,23)]+['X']
        contigs = dict(zip(names,LENGTHS['GRCh37']))
        self.assertEqual(len(check_dictionary(contigs,'GRCh37')),23)
        with self.assertRaisesRegex(ValueError,'mismatch'):
            check_dictionary(contigs,'GRCh38')
        with self.assertRaisesRegex(ValueError,'missing primary'):
            check_dictionary({'1':contigs['1']},'GRCh37')
        with self.assertRaisesRegex(ValueError,'aliases'):
            check_dictionary(dict(contigs,chr1=contigs['1']),'GRCh37')

    def test_vcf_header_build_check(self):
        with tempfile.TemporaryDirectory() as tmp:
            vcf = Path(tmp)/'input.vcf'
            vcf.write_text('##fileformat=VCFv4.2\n##contig=<ID=chr1,length=248956422>\n#CHROM\n')
            with self.assertRaisesRegex(ValueError,'mismatch'):
                check_vcf(vcf,'GRCh37')
            self.assertEqual(check_vcf(vcf,'GRCh38'),['1'])
            vcf.write_text('##fileformat=VCFv4.2\n#CHROM\n')
            with self.assertRaisesRegex(ValueError,'No primary'):
                check_vcf(vcf,'GRCh37')

    def test_public_stage_and_exclusion_guards(self):
        signatures = read_yaml(ROOT/'config/signatures.yaml')
        signatures.update(method='pcawg_clustersv',pcawg_stage='features',pcawg_classified_manifest='supplied.tsv')
        with self.assertRaisesRegex(ValueError,'outside'):
            validate(signatures,'signatures')
        signatures.update(pcawg_stage='clustering',genome_build='GRCh38')
        with self.assertRaisesRegex(ValueError,'GRCh37'):
            validate(signatures,'signatures')
        callers = read_yaml(ROOT/'config/defaults.yaml')
        callers.update(callers='svaba',exclude_bed='mask.bed')
        with self.assertRaisesRegex(ValueError,'supported'):
            validate(callers,'callers')

    def test_exclusions_merge_and_manta_complement(self):
        with tempfile.TemporaryDirectory() as tmp:
            bed = Path(tmp)/'exclude.bed'
            bed.write_text('1\t10\t20\n1\t15\t30\n1\t30\t40\n2\t0\t50\n')
            contigs = {'1':100,'2':50}
            regions = read_exclusions(bed,contigs)
            self.assertEqual(regions,{'1':[(10,40)],'2':[(0,50)]})
            self.assertEqual(complement(contigs,regions),[('1',0,10),('1',40,100)])
            fai = Path(tmp)/'ref.fa.fai'; fai.write_text('1\t100\n2\t50\n')
            out = Path(tmp)/'out'
            self.assertEqual(prepare(bed,fai,out).read_text(),'1\t10\t40\n2\t0\t50\n')
            self.assertEqual(json.loads((out/'regions.audit.json').read_text())['excluded_bases'],80)
            for content in ('chr1\t0\t10\n','1\t-1\t10\n','1\t0\t101\n',''):
                bed.write_text(content)
                with self.assertRaises(ValueError): read_exclusions(bed,contigs)
            with self.assertRaisesRegex(ValueError,'entire reference'):
                complement(contigs,{'1':[(0,100)],'2':[(0,50)]})

    def test_gnu_and_json_timer_resume_formats(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'time.txt'
            path.write_text('\tElapsed (wall clock) time: 0:01\n\tExit status: 0\n')
            self.assertEqual(completed_exit_code(path),0)
            path.write_text(json.dumps({'exit_code':7})+'\n')
            self.assertEqual(completed_exit_code(path),7)
            sidecar = Path(str(path)+'.status.json')
            sidecar.write_text(json.dumps({'exit_code':0,'timing_sha256':hashlib.sha256(path.read_bytes()).hexdigest()}))
            self.assertEqual(completed_exit_code(path),0)
            path.write_text('modified')
            with self.assertRaisesRegex(ValueError,'changed'):
                completed_exit_code(path)

    def test_cgroup_failure_does_not_prevent_caller_execution(self):
        meta = dict(caller='test',pair_id='SYNTHETIC',dataset_id='synthetic',
                    reference='ref.fa',allocated_cpus=1,requested_memory_gb=1,settings={})
        with (tempfile.TemporaryDirectory() as tmp, patch.dict(os.environ,SLURM_JOB_ID='123'),
              patch.object(benchmark,'memory_cgroup',side_effect=PermissionError('fixture')),
              patch.object(benchmark,'sample',return_value=(100,200,1,[]))):
            report = Path(tmp)/'benchmark.json'
            code = benchmark.measure(report,meta,[sys.executable,'-c','raise SystemExit(7)'])
            data = json.loads(report.read_text())
            self.assertEqual(code,7)
            self.assertEqual(data['exit_code'],7)
            self.assertEqual(data['memory_method_id'],'process_tree_sampled_rss_sum')
            self.assertEqual(data['peak_sampled_resident_kib'],200)
            self.assertTrue(data['measurement_warnings'])

    def test_midrun_cgroup_failure_does_not_mix_memory_methods(self):
        meta = dict(caller='test',pair_id='SYNTHETIC',dataset_id='synthetic',
                    reference='ref.fa',allocated_cpus=1,requested_memory_gb=1,settings={})
        with tempfile.TemporaryDirectory() as tmp, patch.dict(os.environ,SLURM_JOB_ID='123'), \
             patch.object(benchmark,'memory_cgroup',return_value=('fixture',Path('/fixture/job_123'))), \
             patch.object(benchmark,'cgroup_resident_kib',side_effect=[900000,PermissionError('fixture')]), \
             patch.object(benchmark,'sample',return_value=(100,200,1,[])):
            report = Path(tmp)/'benchmark.json'
            self.assertEqual(benchmark.measure(report,meta,[sys.executable,'-c','pass']),0)
            data = json.loads(report.read_text())
            self.assertEqual(data['memory_method_id'],'process_tree_sampled_rss_sum')
            self.assertEqual(data['peak_sampled_resident_kib'],200)

if __name__ == '__main__': unittest.main()
