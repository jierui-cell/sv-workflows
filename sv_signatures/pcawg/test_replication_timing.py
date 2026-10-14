#!/usr/bin/env python3
"""Synthetic controls for measured means, missingness, bins and category boundaries."""
import argparse
import csv
import gzip
import json
import shutil
import tempfile
import unittest
from pathlib import Path
import numpy as np
import pyBigWig
from replication_timing import PRIMARY, annotate, build, categories, digest, probabilities, relabel


class TimingTests(unittest.TestCase):
    def test_category_boundaries(self):
        values = np.array([np.nan, 0, 1, 2, 3, 4.])
        self.assertEqual(categories(values, 'panprostate_quartiles', [1.,3.]).tolist(), [0,1,2,2,2,3])
        self.assertEqual(categories(values, 'continuous').tolist(), [0]*6)
        self.assertEqual(categories(values, 'tertiles', [1.,3.]).tolist(), [0,1,2,2,2,3])
        self.assertEqual(probabilities('tertiles'),[1/3,2/3])

    def test_reference_and_annotations(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            # No true patient data: 23 toy chromosomes, three different cell lines.
            sizes = {c: 5500 if c == '1' else 5000 for c in PRIMARY}
            size_file = root / 'sizes.tsv'
            size_file.write_text(''.join(f'{c}\t{n}\n' for c,n in sizes.items()))
            locked = json.loads((Path(__file__).parent / 'replication_sources.json').read_text())
            paths = []
            for offset, source in enumerate(locked['sources']):
                path = root / source['filename']; paths.append(path)
                with pyBigWig.open(str(path), 'w') as bw:
                    bw.addHeader([('chr'+c,n) for c,n in sizes.items()])
                    for chrom,length in sizes.items():
                        # Half a bin missing in the third track on chr1; must be unknown.
                        starts = [0,1000,2000,3000,4000]
                        ends = [1000,2000,3000,4000,5000]
                        if offset == 2 and chrom == '1': starts[2] = 2500
                        bw.addEntries(['chr'+chrom]*5, starts, ends=ends,
                                      values=[float(10*i+offset) for i in range(5)])
                        if length > 5000:
                            bw.addEntries(['chr'+chrom], [5000], ends=[length], values=[51.+offset])
                source['md5'] = digest(path, 'md5')
            lock_path = root / 'sources.json'; lock_path.write_text(json.dumps(locked))
            settings = dict(sources=lock_path, tracks=paths, chrom_sizes=size_file, bin_bp=1000,
                            chunk_bins=2, min_coverage=1., label_policy='panprostate_quartiles')
            output = root / 'reference'
            build(argparse.Namespace(**settings, outdir=output))
            metadata = json.loads((output / 'reference.json').read_text())
            self.assertEqual(metadata['total_bins'], 116)
            self.assertEqual(metadata['unknown_bins'], 1)
            with np.load(output / 'replication_timing.npz') as archive:
                np.testing.assert_allclose(archive['1_score'], [1,11,np.nan,31,41,52], equal_nan=True)
                self.assertEqual(archive['1_label'][2],0)
            with gzip.open(output / 'replication_timing.bed.tsv.gz','rt') as handle:
                rows=list(csv.DictReader(handle, delimiter='\t'))
            self.assertEqual(rows[5]['end'],'5500')
            # Test 1-based positions at both sides of a 1-kb boundary.
            classes = root / 'classes.tsv'
            classes.write_text('id\tV7\tV1\tV4\tV11\tV12\nJ\tJ\t1\tchr1\t1000\t1001\n')
            annotations = root / 'annotation.tsv'
            annotate(argparse.Namespace(reference=output, classified_sv=classes, output=annotations))
            with annotations.open() as handle: row=next(csv.DictReader(handle,delimiter='\t'))
            self.assertEqual((row['signalA1'],row['signalA2']),('1','11'))
            (output/'source').mkdir()
            for path in paths: shutil.copyfile(path,output/'source'/path.name)
            relabel(argparse.Namespace(reference=output,outdir=root/'tertiles',label_policy='tertiles'))
            tertile_meta=json.loads((root/'tertiles'/'reference.json').read_text())
            self.assertEqual(tertile_meta['label_probabilities'],[1/3,2/3])
            self.assertEqual(tertile_meta['unknown_bins'],1)
            with np.load(output/'replication_timing.npz') as a, np.load(root/'tertiles'/'replication_timing.npz') as b:
                for chrom in PRIMARY: np.testing.assert_array_equal(a[chrom+'_score'],b[chrom+'_score'])
            bedpe=root/'positions.bedpe'
            bedpe.write_text('1\t999\t1000\t1\t1000\t1001\tJ\t0\t+\t-\n')
            annotate(argparse.Namespace(reference=root/'tertiles',bedpe=bedpe,output=root/'bedpe.annotation.tsv'))
            with (root/'bedpe.annotation.tsv').open() as handle: row=next(csv.DictReader(handle,delimiter='\t'))
            self.assertEqual((row['signalA1'],row['signalA2']),('1','11'))
            classes.write_text('id\tV7\tV1\tV4\tV11\tV12\nJ\tJ\t1\t1\t2001\t1001\n')
            with self.assertRaisesRegex(ValueError,'lack measured timing'):
                annotate(argparse.Namespace(reference=output, classified_sv=classes, output=root/'unknown.tsv'))
            # Smaller chunks must not change reference values or labels.
            build(argparse.Namespace(**dict(settings,chunk_bins=1), outdir=root/'repeat'))
            with np.load(output/'replication_timing.npz') as a, np.load(root/'repeat'/'replication_timing.npz') as b:
                for key in a.files: np.testing.assert_array_equal(a[key],b[key])
            # Continuous-only profiles cannot silently supply category labels.
            build(argparse.Namespace(**dict(settings,label_policy='continuous'), outdir=root/'continuous'))
            with self.assertRaisesRegex(ValueError,'continuous scores only'):
                annotate(argparse.Namespace(reference=root/'continuous',classified_sv=classes,output=root/'blocked.tsv'))
            with self.assertRaisesRegex(ValueError,'nonempty'):
                build(argparse.Namespace(**settings,outdir=output))
            paths[0].write_bytes(b'corrupt')
            with self.assertRaisesRegex(ValueError,'MD5 mismatch'):
                build(argparse.Namespace(**settings,outdir=root/'bad'))


if __name__ == '__main__':
    unittest.main()
