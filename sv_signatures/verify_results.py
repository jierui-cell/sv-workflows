#!/usr/bin/env python3
"""Verify all published cohort/caller SV32 matrices and COSMIC exposures."""

import argparse
import csv
import json
import re
from collections import defaultdict
from pathlib import Path


def read_tsv(path):
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return reader.fieldnames, list(reader)


def verify(manifest, outdir):
    _, manifest_rows = read_tsv(manifest)
    groups = defaultdict(set)
    manifest_by_key = {}
    for row in manifest_rows:
        key = (row["cohort"], row["caller"], row["pair_id"])
        if key in manifest_by_key:
            raise ValueError("Duplicate manifest row: " + repr(key))
        manifest_by_key[key] = row
        groups[key[:2]].add(key[2])
    if len(groups) != 10:
        raise ValueError("Expected ten cohort/caller groups, found {}".format(len(groups)))
    report = []
    for (cohort, caller), expected in sorted(groups.items()):
        root = outdir / cohort / caller / "group_output"
        qc = json.loads((root / "matrix_qc.json").read_text())
        fit = json.loads((root / "fit_status.json").read_text())
        _, sources = read_tsv(root / "sources.tsv")
        matrix_header, matrix_rows = read_tsv(root / "SV32.matrix.tsv")
        activity_header, activities = read_tsv(root / "SV32.exposures.tsv")
        exposure_qc = json.loads((root / "exposure_qc.json").read_text())
        samples = matrix_header[1:]
        source_ids = {row["pair_id"] for row in sources}
        activity_ids = {row["Samples"] for row in activities}
        if (source_ids != expected or set(samples) != expected or activity_ids != expected or
                len(samples) != len(expected) or len(sources) != len(expected) or
                len(activities) != len(expected)):
            raise ValueError("Sample-set mismatch in {}:{}".format(cohort, caller))
        for source in sources:
            original = manifest_by_key[(cohort, caller, source["pair_id"])]
            if (source["vcf"] != original["vcf"] or
                    not re.fullmatch(r"[0-9a-f]{64}", source["sha256"]) or
                    (original["source_sha256"] and original["source_sha256"] != source["sha256"])):
                raise ValueError("Source VCF provenance mismatch: {} {} {}".format(cohort, caller, source["pair_id"]))
        if len(matrix_rows) != 32 or len({row[matrix_header[0]] for row in matrix_rows}) != 32:
            raise ValueError("Invalid SV32 channels in {}:{}".format(cohort, caller))
        if len(activity_header) != 13 or len(set(activity_header[1:])) != 12:
            raise ValueError("Expected twelve COSMIC SV signatures in {}:{}".format(cohort, caller))
        counts = {sample: sum(int(row[sample]) for row in matrix_rows) for sample in samples}
        clustered = sum(sum(int(row[sample]) for sample in samples) for row in matrix_rows
                        if row[matrix_header[0]].startswith("clustered_"))
        nonclustered = sum(sum(int(row[sample]) for sample in samples) for row in matrix_rows
                           if row[matrix_header[0]].startswith("non-clustered_"))
        retained = {row["pair_id"]: int(row["retained_events"]) for row in sources}
        if counts != retained:
            raise ValueError("Matrix and adapter counts differ in {}:{}".format(cohort, caller))
        if (qc["channels"] != 32 or qc["samples"] != len(expected) or
                qc["total_events"] != sum(counts.values()) or
                qc["clustered_events"] != clustered or qc["nonclustered_events"] != nonclustered or
                clustered + nonclustered != sum(counts.values()) or fit["status"] != "completed"):
            raise ValueError("QC or fit status mismatch in {}:{}".format(cohort, caller))
        if cohort == "melanoma" and clustered == 0:
            raise ValueError("No clustered melanoma events in {}:{}; investigate SV clustering".format(cohort, caller))
        if (exposure_qc["samples"] != len(expected) or
                exposure_qc["zero_event_samples"] != sorted(sample for sample, n in counts.items() if n == 0)):
            raise ValueError("Exposure QC mismatch in {}:{}".format(cohort, caller))
        exposure_totals = {row["Samples"]: sum(float(row[name]) for name in activity_header[1:])
                           for row in activities}
        if any(abs(exposure_totals[sample] - counts[sample]) > 1 for sample in expected):
            raise ValueError("Exposure totals differ from event counts in {}:{}".format(cohort, caller))
        report.append({"cohort": cohort, "caller": caller, "samples": len(expected),
                       "events": sum(counts.values()), "clustered_events": clustered,
                       "zero_event_samples": sum(n == 0 for n in counts.values()),
                       "output": str(root.resolve())})
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=Path("sv_signatures/input_manifest.tsv"))
    parser.add_argument("--outdir", type=Path, default=Path("results/sv_signatures_corrected"))
    parser.add_argument("--output", type=Path, help="Optionally save the verification report as JSON")
    args = parser.parse_args()
    report = verify(args.manifest, args.outdir)
    payload = json.dumps({"status": "verified", "groups": report}, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload)
    print(payload, end="")


if __name__ == "__main__":
    main()
