import csv
import tempfile
import unittest
from pathlib import Path

from run_clustersv import mask_bedpe, references, table
from vcf_to_bedpe import convert, oriented_strands


class ClusterSVAdapterTests(unittest.TestCase):
    def test_ambiguous_inversion_is_not_invented(self):
        record = ['1', '100', 'inv', 'N', '<INV>', '.', 'PASS', 'SVTYPE=INV;END=2000']
        self.assertIsNone(oriented_strands(record))
        record[7] += ';CT=3to3'
        self.assertEqual(oriented_strands(record), '++')
        record[7] += ';STRANDS=--:5'
        self.assertIsNone(oriented_strands(record))

    def test_headerless_and_short_adjacency(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            vcf = root / 'test.vcf'
            vcf.write_text('##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
                           '1\t100\ta\tN\tN[1:150[\t.\tPASS\tSVTYPE=BND;MATEID=b\n'
                           '1\t150\tb\tN\t]1:100]N\t.\tPASS\tSVTYPE=BND;MATEID=a\n'
                           '2\t100\ti\tN\t<INV>\t.\tPASS\tSVTYPE=INV;END=2000\n')
            n = convert(vcf, 'manta', 'P1', root / 'input.bedpe', root / 'audit.tsv', root / 'qc.json', 1, 'clustersv')
            self.assertEqual(n, 1)
            fields = (root / 'input.bedpe').read_text().strip().split('\t')
            self.assertEqual(len(fields), 10)
            self.assertEqual(fields, ['1','99','100','1','149','150','a','0','+','-'])

    def test_delly_fallback_keeps_strands_with_their_endpoints(self):
        record = ['chr2', '1000', 'bnd', 'N', '<BND>', '.', 'PASS',
                  'SVTYPE=BND;CHR2=chr1;POS2=3000;CT=3to5']
        self.assertEqual(oriented_strands(record), '+-')
        record[4] = 'N[chr1:3000['
        self.assertEqual(oriented_strands(record), '-+')  # ALT endpoint is canonically sorted.

    def test_reference_validation_and_explicit_masks(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / 'source'
            refs = source / 'references'; refs.mkdir(parents=True)
            chroms = [str(i) for i in range(1,23)] + ['X']
            (refs / 'hs37d5.chrom_sizes').write_text(''.join(c+'\t150000000\n' for c in chroms))
            (refs / 'hg19_centromere_and_telomere_coords.txt').write_text('chrom\ta\tb\tc\td\te\n'+
                ''.join(c+'\t1000\t60000000\t61000000\t149000000\t150000000\n' for c in chroms))
            _, _, masks, sizes = references(source, root)
            self.assertGreater(masks['17'][-1], 80_000_000)
            self.assertNotIn('Y', masks)
            bedpe = root / 'test.bedpe'
            bedpe.write_text('1\t99\t100\t1\t199\t200\ttel\t0\t+\t-\n'
                             '1\t9999999\t10000000\t1\t10000999\t10001000\tkeep\t0\t+\t-\n'
                             'Y\t99999\t100000\t1\t10000999\t10001000\ty\t0\t+\t-\n')
            self.assertEqual(mask_bedpe(bedpe, masks, sizes, root / 'excluded.tsv'), (1, 2))

    def test_missing_manifest_header_names_the_column(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'bad.tsv'
            path.write_text('cohort\ttumour_vcf\nmelanoma\tfoo\n')
            with self.assertRaisesRegex(ValueError, 'missing columns: vcf'):
                table(path, {'cohort', 'vcf'})

    def test_mask_policy_does_not_claim_upstream_exclusion(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            original = '1\t120\t121\t1\t179\t180\tcentromere\t0\t+\t-\n'
            masks = {'1':[10,100,150,200,990]}
            sizes = {'1':'1000'}
            bed = root/'input.bedpe'; bed.write_text(original)
            self.assertEqual(mask_bedpe(bed,masks,sizes,root/'local.tsv'),(0,1))
            bed.write_text(original)
            self.assertEqual(mask_bedpe(bed,masks,sizes,root/'primary.tsv','primary_only'),(1,0))
            self.assertEqual(bed.read_text(),original)


if __name__ == '__main__':
    unittest.main()
