import sys,unittest,tempfile
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from run_cue import configs,interpreter_path
class CueConfigTest(unittest.TestCase):
 def test_virtualenv_interpreter_symlink_is_preserved(self):
  with tempfile.TemporaryDirectory() as tmp:
   root=Path(tmp); target=root/'system-python';target.touch()
   link=root/'venv'/'bin'/'python';link.parent.mkdir(parents=True);link.symlink_to(target)
   self.assertEqual(interpreter_path(root,'venv/bin/python'),str(link))
   self.assertNotEqual(interpreter_path(root,'venv/bin/python'),str(link.resolve()))
 def test_config(self):
  s=dict(root='/project',role='tumor',bam='/tumor.bam',fai='/ref.fai',settings=dict(cpus=2,model='cue.pt',data_options=dict(chr_names=['chr21']),model_options=dict(gpu_ids=[],n_cpus=2)))
  d,m=configs(s);self.assertEqual(d['bam'],'/tumor.bam');self.assertEqual(m['model_path'],'/project/cue.pt')
  s['settings']['model_options']['gpu_ids']=[0]
  with self.assertRaises(ValueError):configs(s)
