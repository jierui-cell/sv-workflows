#!/usr/bin/env python3
"""Positive-control the pinned SigProfiler SV clustering implementation."""

import tempfile

import pandas as pd

from run_group import install_sigprofiler_sv_compatibility


def synthetic_events():
    rows = []
    for i in range(50):
        start = 10_000_000 + i * 5_000
        rows.append(dict(sample="CONTROL", chrom1="1", start1=start, end1=start + 1,
                         chrom2="1", start2=start + 2_000, end2=start + 2_001,
                         svclass="deletion"))
    for i in range(30):
        start = 30_000_000 + i * 4_000_000
        rows.append(dict(sample="CONTROL", chrom1="1", start1=start, end1=start + 1,
                         chrom2="1", start2=start + 50_000, end2=start + 50_001,
                         svclass="deletion"))
    return pd.DataFrame(rows)


def main():
    install_sigprofiler_sv_compatibility()
    from SigProfilerMatrixGenerator.scripts import SVMatrixGenerator as sv

    result = sv.annotateBedpe(synthetic_events())
    annotated = result["sv_bedpe"]
    clustered = int(annotated["is_clustered"].sum())
    if len(result["kat_regions"]) < 1 or clustered < 30 or clustered >= len(annotated):
        raise AssertionError("Synthetic SV clusters were not detected correctly")
    with tempfile.TemporaryDirectory() as directory:
        matrix = sv.tsv2matrix([annotated], "control", directory)
        clustered_count = int(matrix.loc[matrix.index.str.startswith("clustered_"), "CONTROL"].sum())
        nonclustered_count = int(matrix.loc[matrix.index.str.startswith("non-clustered_"), "CONTROL"].sum())
    if clustered_count != clustered or clustered_count + nonclustered_count != len(annotated):
        raise AssertionError("SV32 matrix did not retain the clustered/non-clustered split")
    print("Clustering positive control passed: {} clustered, {} non-clustered events".format(
        clustered_count, nonclustered_count))


if __name__ == "__main__":
    main()
