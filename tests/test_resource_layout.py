"""Selected-caller setup checks use synthetic files, never patient data or inference."""
from contextlib import redirect_stderr, redirect_stdout
import csv
import io
import json
import os
import re
from pathlib import Path
import sys
import tempfile
from unittest.mock import patch, MagicMock
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from check_resources import (BRASS_FILES, BWA_SUFFIXES, SAMPLE_FIELDS, MANIFEST_FIELDS,
                             absolute, check_files, normalize_paths, write_resolved_samplesheet)
import workflow

class ResourceLayoutTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.config = workflow.read_yaml(ROOT/'config/defaults.yaml')
        self.config['callers'] = 'manta,delly'
        for name in ('resources/reference/reference.fa','resources/reference/reference.fa.fai',
                     self.config['manta_image'],self.config['delly']['image']):
            self.touch(name)
        self.sheet(['PAIR_001','PAIR_002'])

    def touch(self,name,content='synthetic fixture\n'):
        path = self.root/name
        path.parent.mkdir(parents=True,exist_ok=True)
        path.write_text(content)
        return path

    def sheet(self,pairs):
        path = self.root/'inputs/samplesheet.tsv'; path.parent.mkdir(parents=True,exist_ok=True)
        with path.open('w',newline='') as handle:
            writer = csv.writer(handle,delimiter='\t')
            writer.writerow(SAMPLE_FIELDS)
            for pair in pairs:
                files = [f'inputs/bams/{pair}.{role}.{extension}'
                         for role in ('tumor','normal') for extension in ('bam','bam.bai')]
                for name in files: self.touch(name)
                writer.writerow([pair]+files+['unknown','','',''])

    def report(self): return check_files(self.config,'callers',self.root)

    def test_defaults_match_downloaded_images(self):
        with (ROOT/'config/containers.tsv').open() as handle:
            rows = list(csv.DictReader(handle,delimiter='\t'))
        for row in rows:
            if row['tool'] == 'bwa': continue
            configured = self.config['manta_image'] if row['tool'] == 'manta' else self.config[row['tool']]['image']
            self.assertEqual(configured,f"containers/images/{row['tool']}-{row['tag']}.sif")
        self.assertEqual(self.config['input'],'inputs/samplesheet.tsv')
        self.assertNotIn('/path/to/',(ROOT/'config/defaults.yaml').read_text())
        source = (ROOT/'lib/Samplesheet.groovy').read_text()
        schema = re.search(r'static final List REQUIRED = \[(.*?)\]',source,re.S).group(1)
        self.assertEqual(tuple(re.findall(r"'([^']+)'",schema)),SAMPLE_FIELDS)

    def test_selected_callers_do_not_require_unused_resources(self):
        result = self.report()
        self.assertEqual(result['status'],'passed',result['issues'])
        self.assertEqual(result['selected_rows'],2)
        self.config['callers'] = 'svaba'
        result = self.report()
        self.assertEqual(result['status'],'failed')
        self.assertTrue(any('svaba image' in error for error in result['issues']))
        self.assertTrue(any('BWA index .bwt' in error for error in result['issues']))
        for suffix in BWA_SUFFIXES: self.touch(self.config['reference']+suffix)
        self.touch(self.config['svaba']['image'])
        self.assertEqual(self.report()['status'],'passed')

    def test_missing_files_empty_images_and_broken_links(self):
        image = self.root/self.config['manta_image']; image.write_text('')
        index = self.root/'inputs/bams/PAIR_001.tumor.bam.bai'; index.unlink()
        index.symlink_to(self.root/'missing.bai')
        result = self.report()
        self.assertTrue(any('manta image' in error and 'empty' in error for error in result['issues']))
        self.assertTrue(any('PAIR_001 tumor_bai' in error and 'missing' in error for error in result['issues']))

    def test_pair_selection_and_header_errors(self):
        (self.root/'inputs/bams/PAIR_002.normal.bam').unlink()
        self.config['pair_ids'] = 'PAIR_001'
        self.assertEqual(self.report()['status'],'passed')
        self.config['pair_ids'] = 'UNKNOWN'
        self.assertTrue(any('unknown pair_ids' in error for error in self.report()['issues']))
        self.config['pair_ids'] = ''
        path = self.root/'inputs/samplesheet.tsv'
        path.write_text(path.read_text().replace('tumor_bam','tumour_bam',1))
        self.assertTrue(any('tumor_bam' in error and 'Missing columns' in error for error in self.report()['issues']))

    def test_linked_reference_source_indices_match_native_adapter(self):
        source = self.touch('external/alignment.fa')
        alias = self.root/self.config['reference']; alias.unlink(); alias.symlink_to(source)
        self.assertTrue(any('linked reference source FASTA index' in error for error in self.report()['issues']))
        self.touch(str(source)+'.fai')
        self.assertEqual(self.report()['status'],'passed')

    def test_unreadable_files_are_reported(self):
        with patch('check_resources.os.access',return_value=False):
            report = self.report()
        self.assertEqual(report['status'],'failed')
        self.assertTrue(any('not readable' in error for error in report['issues']))

    def test_relative_paths_are_repository_relative_and_venv_links_survive(self):
        cfg = normalize_paths(self.config,'callers',self.root)
        self.assertEqual(cfg['reference'],str(self.root/'resources/reference/reference.fa'))
        self.assertEqual(self.config['reference'],'resources/reference/reference.fa')
        interpreter = self.root/'sv_signatures/.venv/bin/python'
        interpreter.parent.mkdir(parents=True); interpreter.symlink_to(sys.executable)
        cfg = normalize_paths(workflow.read_yaml(ROOT/'config/signatures.yaml'),'signatures',self.root)
        self.assertEqual(cfg['python'],str(interpreter))
        self.assertNotEqual(cfg['python'],str(interpreter.resolve()))
        self.assertEqual(absolute(str(self.root/'outside.fa'),ROOT),self.root/'outside.fa')

    def test_normalized_sheet_is_stable_for_resume_and_original_is_unchanged(self):
        original = (self.root/'inputs/samplesheet.tsv').read_bytes()
        first = write_resolved_samplesheet(self.config,self.root,self.root/'first.tsv')
        modified = first.stat().st_mtime_ns
        second = write_resolved_samplesheet(self.config,self.root,self.root/'second.tsv')
        self.assertEqual(first,second)
        self.assertEqual(first.stat().st_mtime_ns,modified)
        self.assertEqual((self.root/'inputs/samplesheet.tsv').read_bytes(),original)
        with first.open() as handle:
            row = next(csv.DictReader(handle,delimiter='\t'))
        self.assertTrue(Path(row['tumor_bam']).is_absolute())
        self.config['pair_ids'] = 'PAIR_001'
        selected = write_resolved_samplesheet(self.config,self.root,self.root/'selected.tsv')
        self.assertNotEqual(selected,first)

    def test_gridss_dictionary_and_brass_native_bundle_checks(self):
        self.config['callers'] = 'gridss,brass'
        for caller in ('gridss','brass'): self.touch(self.config[caller]['image'])
        for suffix in BWA_SUFFIXES: self.touch(self.config['reference']+suffix)
        self.touch('resources/reference/reference.dict')
        bundle = self.root/self.config['brass']['bundle']
        for name in BRASS_FILES: self.touch(str(bundle/name))
        self.touch(str(bundle/'COMPATIBILITY_VALIDATED'),'')
        self.touch(str(bundle/'native-validation.json'),json.dumps({'status':'passed'}))
        self.touch(str(bundle/'bundle-receipt.json'),json.dumps({'status':'passed','assembly':'GRCh37'}))
        self.assertEqual(self.report()['status'],'passed')
        self.touch(str(bundle/'native-validation.json'),json.dumps({'status':'failed'}))
        self.assertTrue(any('not passed' in error for error in self.report()['issues']))
        (bundle/'all_ncbi_bacteria.4.fa.2bit').unlink()
        self.assertTrue(any('all_ncbi_bacteria.4' in error for error in self.report()['issues']))

    def test_signature_environment_and_manifest_inputs(self):
        cfg = workflow.read_yaml(ROOT/'config/signatures.yaml')
        cfg['python'] = sys.executable
        path = self.root/cfg['manifest']; path.parent.mkdir()
        vcf = self.touch('inputs/synthetic.vcf')
        with path.open('w') as handle:
            writer = csv.writer(handle,delimiter='\t')
            writer.writerow(MANIFEST_FIELDS)
            writer.writerow(('synthetic','PAIR_001','manta',str(vcf),'synthetic_control','0'*64))
        self.assertEqual(check_files(cfg,'signatures',self.root)['status'],'passed')
        vcf.unlink()
        self.assertEqual(check_files(cfg,'signatures',self.root)['status'],'failed')

    def test_checker_cli_is_read_only_and_automatic_failures_submit_nothing(self):
        config = self.root/'config/defaults.yaml'; config.parent.mkdir()
        config.write_text((ROOT/'config/defaults.yaml').read_text())
        for name in ('main.nf','nextflow.config'): self.touch(name)
        for option in ('--check-files',None):
            argv = ['workflow.py','callers','--set','callers=manta,delly'] + ([option] if option else [])
            (self.root/'inputs/bams/PAIR_001.normal.bam').unlink(missing_ok=True)
            with patch.object(workflow,'ROOT',self.root),patch.object(sys,'argv',argv),\
                    patch.object(workflow.subprocess,'Popen') as submit,\
                    redirect_stdout(io.StringIO()) as output,redirect_stderr(io.StringIO()):
                with self.assertRaises(SystemExit) as stopped: workflow.main()
                self.assertEqual(stopped.exception.code,2)
                submit.assert_not_called()
                if option: self.assertEqual(json.loads(output.getvalue())['status'],'failed')
            self.assertFalse((self.root/'logs').exists())
            self.assertFalse((self.root/'.nextflow').exists())

    def test_real_launcher_uses_checked_absolute_input_and_keeps_snapshot(self):
        path = self.root/'config/defaults.yaml'; path.parent.mkdir()
        path.write_text((ROOT/'config/defaults.yaml').read_text())
        for name in ('main.nf','nextflow.config'): self.touch(name)
        process = MagicMock(); process.stdout = iter(['synthetic completion\n']); process.wait.return_value = 0
        with patch.object(workflow,'ROOT',self.root),patch.object(sys,'argv',['workflow.py','callers']),\
                patch.object(workflow.shutil,'which',return_value='/synthetic/nextflow'),\
                patch.object(workflow.subprocess,'check_output',return_value='nextflow version 25.10.7'),\
                patch.object(workflow.subprocess,'Popen',return_value=process),redirect_stdout(io.StringIO()):
            with self.assertRaises(SystemExit) as stopped: workflow.main()
            self.assertEqual(stopped.exception.code,0)
        logs = list((self.root/'logs').iterdir()); self.assertEqual(len(logs),1)
        resolved = workflow.read_yaml(logs[0]/'resolved.yaml')
        self.assertTrue(Path(resolved['input']).is_file())
        self.assertEqual((logs[0]/'setup-check.json').is_file(),True)
        self.assertEqual((logs[0]/'resolved-samplesheet.tsv').read_text(),Path(resolved['input']).read_text())

if __name__ == '__main__': unittest.main()
