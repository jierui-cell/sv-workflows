"""Synthetic module managers exercise environment transfer and driver bootstrapping."""
from contextlib import redirect_stderr, redirect_stdout
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from hpc_environment import load_modules, module_names
import workflow


class HpcEnvironmentTests(unittest.TestCase):
    def test_module_list_validation_and_empty_installations(self):
        self.assertEqual(module_names(['java/20','samtools/1.20','java/20']),['java/20','samtools/1.20'])
        with patch('hpc_environment.subprocess.run') as command:
            self.assertEqual(load_modules([]),{'modules':[]})
            command.assert_not_called()
        for value in ('java/20', ['--force'], ['x;touch-file'], ['x y'], ['$(id)'], ['one:two'], [None]):
            with self.subTest(value=value), self.assertRaises(ValueError):
                module_names(value)

    def test_module_function_environment_is_preserved_without_polluting_stdout(self):
        function = '() { export SV_MODULE_TEST="loaded $*"; export LOADEDMODULES="java/test:samtools/test"; echo module-banner; }'
        with patch.dict(os.environ,{'BASH_FUNC_module%%':function,'BASH_ENV':'/dev/null','SV_KEEP_TEST':'keep'}), \
             redirect_stdout(io.StringIO()) as out, redirect_stderr(io.StringIO()):
            report = load_modules(['java/test','samtools/test'])
            self.assertEqual(os.environ['SV_MODULE_TEST'],'loaded load java/test samtools/test')
            self.assertEqual(os.environ['SV_KEEP_TEST'],'keep')
        self.assertEqual(out.getvalue(),'')
        self.assertEqual(report['loaded_modules'],['java/test','samtools/test'])

    def test_module_failure_keeps_the_original_environment(self):
        function = '() { echo synthetic-module-error >&2; return 9; }'
        with patch.dict(os.environ,{'BASH_FUNC_module%%':function,'BASH_ENV':'/dev/null','SV_KEEP_TEST':'keep'}):
            before = dict(os.environ)
            with self.assertRaisesRegex(ValueError,'synthetic-module-error'):
                load_modules(['not-installed/1'])
            self.assertEqual(dict(os.environ),before)

    def test_signatures_and_tests_inherit_only_the_cohort_modules(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root/'config').mkdir()
            (root/'config/defaults.yaml').write_text((ROOT/'config/defaults.yaml').read_text())
            (root/'config/private-cohort.yaml').write_text('modules: [java/test, samtools/test]\n')
            for mode in ('signatures','test'):
                cfg = workflow.read_yaml(ROOT/'config'/workflow.DEFAULTS[mode])
                self.assertEqual(workflow.effective_modules(cfg,root),['java/test','samtools/test'])
                cfg['modules'] = []
                self.assertEqual(workflow.effective_modules(cfg,root),[])

    def test_slurm_driver_selects_repository_python_without_activation(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root/'scripts').mkdir(); (root/'.venv/bin').mkdir(parents=True)
            (root/'main.nf').touch()
            (root/'scripts/workflow.py').write_text('import json,os,sys; print(json.dumps({"python":sys.executable,"path":os.environ["PATH"],"args":sys.argv[1:]}))\n')
            python = root/'.venv/bin/python'; python.symlink_to(sys.executable)
            spool = root/'spooled-script'
            spool.write_text((ROOT/'scripts/run.slurm.sh').read_text())
            env = dict(os.environ,SLURM_JOB_ID='synthetic',SLURM_SUBMIT_DIR=str(root))
            env.pop('SV_PYTHON',None); env.pop('SV_WORKFLOW_ROOT',None)
            result = subprocess.run(['bash',str(spool),'callers','--dry-run'],env=env,text=True,capture_output=True,check=True)
            report = json.loads(result.stdout)
            self.assertEqual(report['python'],str(python))
            self.assertEqual(report['path'].split(os.pathsep)[0],str(python.parent))
            self.assertEqual(report['args'],['callers','--dry-run','--profile','slurm'])

    def test_launcher_loads_modules_before_invoking_nextflow(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root/'config').mkdir(); (root/'tools').mkdir()
            for name in ('main.nf','nextflow.config'):
                (root/name).write_text((ROOT/name).read_text())
            (root/'config/test.yaml').write_text('outdir: results/toy\nmodules: [java/test]\ncaller_config: null\n')
            executable = root/'tools/nextflow'
            executable.write_text('#!/usr/bin/env bash\nset -eu\ntest "$SV_MODULE_TEST" = loaded\nif [[ "$1" == -version ]]; then echo "nextflow version 25.10.7"; else echo "synthetic task completed"; fi\n')
            executable.chmod(0o755)
            function = '() { export SV_MODULE_TEST=loaded; export PATH="$SV_TEST_TOOLS:$PATH"; }'
            with patch.dict(os.environ,{'BASH_FUNC_module%%':function,'BASH_ENV':'/dev/null','SV_TEST_TOOLS':str(root/'tools')}), \
                 patch.object(workflow,'ROOT',root), \
                 patch.object(sys,'argv',['workflow.py','test','--profile','test']), \
                 redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()):
                with self.assertRaises(SystemExit) as stopped:
                    workflow.main()
            self.assertEqual(stopped.exception.code,0)
            log = next((root/'logs').glob('test-*'))
            self.assertEqual(json.loads((log/'environment.json').read_text())['modules'],['java/test'])
            self.assertIn('synthetic task completed',(log/'driver.log').read_text())

    def test_failed_module_setup_does_not_start_nextflow(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); (root/'config').mkdir()
            (root/'config/test.yaml').write_text('outdir: results/toy\nmodules: [missing/1]\n')
            with patch.object(workflow,'ROOT',root), patch.object(sys,'argv',['workflow.py','test','--profile','test']), \
                 patch.object(workflow,'load_modules',side_effect=ValueError('missing module')), \
                 patch.object(workflow.subprocess,'Popen') as submit, \
                 redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()):
                with self.assertRaises(SystemExit):
                    workflow.main()
            submit.assert_not_called()
            log = next((root/'logs').glob('test-*'))
            self.assertEqual(json.loads((log/'status.json').read_text())['stage'],'module_setup')


if __name__ == '__main__':
    unittest.main()
