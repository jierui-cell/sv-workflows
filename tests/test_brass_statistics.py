import sys
from pathlib import Path
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from brass_statistics import resolve_statistics


class StatisticsTests(unittest.TestCase):
    def test_defaults_explicit_and_marked(self):
        row={'purity':'','ploidy':'','sex':'female'}
        with self.assertRaises(ValueError): resolve_statistics(row,'require')
        resolved,text=resolve_statistics(row,'native_defaults')
        self.assertEqual(resolved['statistics_state'],'DEFAULT')
        self.assertEqual(resolved['purity'],0.75)
        self.assertEqual(resolved['ploidy'],2)
        self.assertIn('Ploidy ?\n',text)
        self.assertIn('GenderChr chrY\nGenderChrFound N\n',text)

    def test_measured_values_preserved(self):
        resolved,text=resolve_statistics({'purity':'0.953','ploidy':'1.842','sex':'female'},'native_defaults')
        self.assertEqual(resolved['statistics_state'],'REAL')
        self.assertEqual(resolved['purity'],0.953)

    def test_invalid_and_partial_values_rejected(self):
        for p,l in [('nan','2'),('1.1','2'),('bad','2'),('0','2')]:
            with self.assertRaises(ValueError):
                resolve_statistics({'purity':p,'ploidy':l,'sex':'male'},'native_defaults')

    def test_partial_preserves_measured_value(self):
        for p,l,expected in [('0.9','',(0.9,2.0)),('','3.1',(0.75,3.1))]:
            r,t=resolve_statistics({'purity':p,'ploidy':l,'sex':'female'})
            self.assertEqual((r['purity'],r['ploidy']),expected)
            self.assertEqual(r['statistics_state'],'MIXED')
            self.assertIn('# BRASS_PIPELINE_STATE MIXED',t)
