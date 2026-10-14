import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from run_jabba_coverage import bind_directories, validate_options


class JaBbACoverageTests(unittest.TestCase):
    def test_documented_defaults(self):
        self.assertEqual(validate_options({}), (1000, 20, 1e-5))

    def test_invalid_settings_rejected(self):
        for options in ({"window": 0}, {"minmapq": 100}, {"cnsignif": 0},
                        {"made_up": 1}):
            with self.assertRaises(ValueError):
                validate_options({"coverage_options": options})

    def test_bam_parent_outside_project_is_mounted_read_only(self):
        project = Path("/scratch/project")
        bams = {"tumor": Path("/scratch/data/t.bam"),
                "normal": Path("/scratch/data/n.bam")}
        self.assertEqual(bind_directories(project, bams),
                         "/scratch/data:/scratch/data:ro,/scratch/project:/scratch/project")


if __name__ == "__main__":
    unittest.main()
