import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from gridss_identity import validate_identity

class IdentityTest(unittest.TestCase):
    def test_identity_and_swaps(self):
        row={'normal_bam':'/normal.bam','tumor_bam':'/tumor.bam'}
        spec={'normal':'normal.bam','tumor':'tumor.bam'}
        argv=['gridss','--labels','NORMAL,TUMOR','/normal.bam','/tumor.bam']
        validate_identity(row,spec,[{'step':'raw_joint_call','argv':argv}],['NORMAL','TUMOR'])
        for bad in [argv[:-2]+argv[-2:][::-1],argv[:2]+['TUMOR,NORMAL']+argv[3:]]:
            with self.assertRaises(ValueError):
                validate_identity(row,spec,[{'step':'raw_joint_call','argv':bad}],['NORMAL','TUMOR'])
