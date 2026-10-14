import csv
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
sys.path.insert(0,str(ROOT/'sv_signatures'))
from workflow import read_yaml, validate, merge
from regions import subtract_targets, prepare
from clustersv_features import channel, size_bin, build, CHANNELS
from run_clustersv import CLUSTER_COLUMNS

class CaptureAndNMFTTests(unittest.TestCase):
    def test_capture_template_and_policies(self):
        cfg = merge(read_yaml(ROOT/'config/defaults.yaml'),read_yaml(ROOT/'examples/capture.yaml'))
        validate(cfg,'callers')
        cfg['sequencing_type'] = 'targeted'
        with self.assertRaisesRegex(ValueError,'target_bed'): validate(cfg,'callers')
        cfg['target_bed'] = 'panel.bed'
        validate(cfg,'callers')
        cfg['callers'] = 'manta,delly'
        with self.assertRaisesRegex(ValueError,'capture|Capture'): validate(cfg,'callers')
        cfg['callers'] = 'manta'; cfg['sequencing_type'] = 'WGS'
        with self.assertRaisesRegex(ValueError,'target_bed'): validate(cfg,'callers')

    def test_target_subtraction(self):
        targets = {'1':[(10,100),(200,250)]}
        exclusions = {'1':[(0,20),(30,40),(80,220)]}
        self.assertEqual(subtract_targets(targets,exclusions),[('1',20,30),('1',40,80),('1',220,250)])
        self.assertEqual(subtract_targets(targets,{}),[('1',10,100),('1',200,250)])
        with self.assertRaisesRegex(ValueError,'No capture'):
            subtract_targets(targets,{'1':[(0,300)]})
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root/'ref.fai').write_text('1\t300\n')
            (root/'targets.bed').write_text('1\t20\t30\n1\t10\t20\n')
            prepare(None,root/'ref.fai',root/'out',targets=root/'targets.bed')
            qc = json.loads((root/'out/regions.audit.json').read_text())
            self.assertEqual(qc['target_bases'],20)
            self.assertEqual(qc['excluded_bases'],0)

    def test_size_bins_and_orientation(self):
        expected = ['lt1kb','1kb_10kb','10kb_100kb','100kb_1Mb','1Mb_10Mb','ge10Mb']
        self.assertEqual([size_bin(i) for i in (999,1000,10000,100000,1000000,10000000)],expected)
        row = dict(chrom1='1',start1='100',chrom2='1',start2='1100',strand1='+',strand2='-',cluster_size='2')
        self.assertEqual(channel(row),'clustered:DEL:1kb_10kb')
        row.update(start1='1100',start2='100')
        self.assertEqual(channel(row),'clustered:DUP:1kb_10kb')
        row.update(chrom2='2')
        self.assertEqual(channel(row),'clustered:TRA')
        self.assertEqual(len(set(CHANNELS)),38)

    def test_count_conservation_zero_samples_and_membership(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); (root/'clusters').mkdir()
            row = dict(zip(CLUSTER_COLUMNS,['1','100','101','1','1100','1101','sv1','0','+','-','1','1','1','1','x','x','0']))
            for pair, rows in [('A',[row]),('Z',[])]:
                with (root/'clusters'/f'{pair}.clusters.tsv').open('w') as handle:
                    writer = csv.DictWriter(handle,fieldnames=CLUSTER_COLUMNS,delimiter='\t')
                    writer.writeheader(); writer.writerows(rows)
            qc = build(root,['A','Z'])
            self.assertEqual(qc['sample_totals'],{'A':1,'Z':0})
            self.assertEqual(qc['feature_total'],1)
            with (root/'ClusterSV38.matrix.tsv').open() as handle:
                matrix = list(csv.DictReader(handle,delimiter='\t'))
            self.assertEqual(len(matrix),38)
            self.assertEqual(sum(int(r['A']) for r in matrix),1)
            self.assertEqual(sum(int(r['Z']) for r in matrix),0)
            row['cluster_size'] = '2'
            with (root/'clusters/A.clusters.tsv').open('w') as handle:
                writer = csv.DictWriter(handle,fieldnames=CLUSTER_COLUMNS,delimiter='\t')
                writer.writeheader(); writer.writerow(row)
            with self.assertRaisesRegex(ValueError,'membership'): build(root,['A','Z'])

    def test_nmf_settings_and_capture_signature_guard(self):
        cfg = read_yaml(ROOT/'config/signatures.yaml'); cfg['method'] = 'pcawg_clustersv'
        for stage in ('clustering','features','signatures'):
            cfg['pcawg_stage'] = stage; validate(cfg,'signatures')
        for key,value in [('pcawg_rank',1),('pcawg_nrun',0),('pcawg_rank_min',1),('pcawg_nmf_epsilon',0)]:
            broken = dict(cfg,**{key:value})
            with self.assertRaises(ValueError): validate(broken,'signatures')
        cfg['sequencing_type'] = 'WES'
        with self.assertRaisesRegex(ValueError,'WGS'): validate(cfg,'signatures')

if __name__ == '__main__': unittest.main()
