import csv
import gzip
import json
import tempfile
import unittest
from pathlib import Path

from vcf_to_bedpe import bnd_endpoint, classify, convert


class ConversionTests(unittest.TestCase):
    def test_orientations(self):
        for alt, expected in (("A[chr1:2000[", "deletion"),
                              ("]chr1:2000]A", "tandem-duplication"),
                              ("[chr1:2000[A", "inversion"),
                              ("A]chr1:2000]", "inversion")):
            fields = ["chr1", "1000", "x", "A", alt, ".", "PASS", "SVTYPE=BND"]
            self.assertEqual(classify(fields)[0][4], expected)

    def test_recipient_mate(self):
        first = bnd_endpoint("chr1", 1000, "A[chr1:2000[")
        second = bnd_endpoint("chr1", 2000, "]chr1:1000]A")
        self.assertEqual(first[:4], second[:4])

    def test_collapse_filter_and_size(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            vcf = root / "test.vcf.bgz"
            with gzip.open(vcf, "wt") as handle:
                handle.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
                handle.write("chr1\t1000\ta\tA\tA[chr1:2000[\t.\tPASS\tSVTYPE=BND;MATEID=b\n")
                handle.write("chr1\t2000\tb\tA\t]chr1:1000]A\t.\tPASS\tSVTYPE=BND;MATEID=a\n")
                handle.write("chr2\t100\tc\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=500\n")
                handle.write("chr3\t100\td\tA\t<DUP>\t.\tLowQual\tSVTYPE=DUP;END=3000\n")
            n = convert(vcf, "manta", "P1", root / "P1.bedpe.tsv", root / "audit.tsv", root / "summary.json")
            self.assertEqual(n, 1)
            summary = json.loads((root / "summary.json").read_text())
            self.assertEqual(summary["decisions"], {"retained": 1, "below_min_size": 1})
            self.assertEqual(summary["selected_records"], 3)
            with (root / "P1.bedpe.tsv").open() as handle:
                event = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual((event["chrom1"], event["chrom2"]), ("1", "1"))

    def test_brass_dot_filter(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            vcf = root / "test.vcf"
            vcf.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                           "chr1\t1000\ta\tA\tA[chr1:2000[\t4\t.\tSVTYPE=BND;MATEID=b;SVCLASS=deletion\n"
                           "chr1\t2000\tb\tA\t]chr1:1000]A\t4\t.\tSVTYPE=BND;MATEID=a;SVCLASS=deletion\n")
            n = convert(vcf, "brass", "P1", root / "P1.bedpe.tsv", root / "audit.tsv", root / "summary.json")
            self.assertEqual(n, 1)


if __name__ == "__main__":
    unittest.main()
