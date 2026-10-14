#!/usr/bin/env python3
"""Offline installation checks. Optional numerical controls require their pinned environments."""
from pathlib import Path
import subprocess
import sys
import unittest
ROOT = Path(__file__).resolve().parents[1]
loader=unittest.TestLoader()
suite=loader.discover(str(ROOT/'tests'))
sys.path.insert(0,str(ROOT/'sv_signatures'))
for name in ('test_vcf_to_bedpe','test_clustersv_adapter','test_pcawg_features'):
    suite.addTests(loader.loadTestsFromName(name))
result=unittest.TextTestRunner(verbosity=2).run(suite)
raise SystemExit(0 if result.wasSuccessful() else 1)
