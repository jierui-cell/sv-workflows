#!/usr/bin/env python3
"""Execute sparse/zero-sample, repeatability and rank-survey NMF controls."""
import argparse
import csv
import json
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from clustersv_features import CHANNELS, PROFILE

def rows(path):
    with path.open() as handle:
        return list(csv.DictReader(handle,delimiter='\t'))

def run(rscript):
    with tempfile.TemporaryDirectory(prefix='sv-nmf-control-') as tmp:
        root = Path(tmp); matrix = root/'counts.tsv'
        with matrix.open('w') as handle:
            writer = csv.writer(handle,delimiter='\t')
            writer.writerow(['Feature'] + [f'S{i}' for i in range(6)] + ['ZERO'])
            for index,label in enumerate(CHANNELS):
                # Two independently varying components, four positive channels.
                vals = ([80,100,120,2,4,6] if index in (0,1) else
                        [2,4,6,80,100,120] if index in (6,7) else [0]*6)
                writer.writerow([label] + vals + [0])
        snapshots = []
        for name, rank in [('fixed1',2),('fixed2',2),('survey',0)]:
            out = root/name; out.mkdir()
            subprocess.run([rscript,str(ROOT/'pcawg/nmf.R'),str(matrix),str(out),
                            '2','3',str(rank),'3','1234','1','2.220446049250313e-16','500',
                            PROFILE,'ClusterSV.NMF'],check=True)
            status = json.loads((out/'signature_status.json').read_text())
            assert status['status'] == 'completed' and status['feature_profile'] == PROFILE
            assert status['channels'] == 38 and status['positive_samples'] == 6
            signatures = rows(out/'ClusterSV.NMF.signatures.tsv')
            names = [key for key in signatures[0] if key != 'Feature']
            for label in names:
                assert abs(sum(float(row[label]) for row in signatures)-1) < 1e-6
            activity = rows(out/'ClusterSV.NMF.activities.tsv')
            zero = next(row for row in activity if row['Samples'] == 'ZERO')
            assert all(float(zero[label]) == 0 for label in names)
            qc = rows(out/'NMF.reconstruction_qc.tsv')
            assert all(float(row['cosine']) > .99 for row in qc if row['Samples'] != 'ZERO')
            if rank:
                snapshots.append((out/'ClusterSV.NMF.signatures.tsv').read_bytes())
            else:
                assert (out/'NMF.rank_survey.tsv').is_file()
                assert 2 <= status['rank'] <= 3
        assert snapshots[0] == snapshots[1], 'Seeded single-thread profiles changed across repeat runs'
    print('PASS: sparse CN-free NMF, zero samples, repeated seeded fits, rank survey and exports')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rscript',default='Rscript')
    run(parser.parse_args().rscript)
