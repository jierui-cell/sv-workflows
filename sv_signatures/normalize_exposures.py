#!/usr/bin/env python3
"""Make a complete SV32 COSMIC activity table, adding zero-event samples."""

import argparse
import csv
import json
from pathlib import Path


def rows(path):
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return reader.fieldnames, list(reader)


def normalize(group_output):
    _, sources = rows(group_output / "sources.tsv")
    names, native = rows(group_output / "cosmic_fit/Assignment_Solution/Activities/Assignment_Solution_Activities.txt")
    if names[0] != "Samples" or len(names) != 13 or len(set(names)) != 13:
        raise ValueError("Expected sample plus twelve COSMIC SV signatures")
    counts = {row["pair_id"]: int(row["retained_events"]) for row in sources}
    if len(counts) != len(sources):
        raise ValueError("Duplicate pair in sources.tsv")
    native_by_sample = {row["Samples"]: row for row in native}
    if len(native_by_sample) != len(native):
        raise ValueError("Duplicate sample in native activities")
    positive = {sample for sample, n in counts.items() if n > 0}
    if not positive <= native_by_sample.keys() or not native_by_sample.keys() <= counts.keys():
        raise ValueError("Native activity sample set does not match positive-event samples")
    normalized = []
    for sample in sorted(counts):
        row = native_by_sample.get(sample, {"Samples": sample, **{name: "0" for name in names[1:]}})
        exposure_sum = sum(float(row[name]) for name in names[1:])
        if abs(exposure_sum - counts[sample]) > 1:
            raise ValueError("Exposure count mismatch: " + sample)
        normalized.append(row)
    out = group_output / "SV32.exposures.tsv"
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=names, delimiter="\t")
        writer.writeheader()
        writer.writerows(normalized)
    qc = {"samples": len(counts), "signatures": names[1:], "total_events": sum(counts.values()),
          "zero_event_samples": sorted(sample for sample, n in counts.items() if n == 0),
          "native_activity_rows": len(native), "standardized_activity_rows": len(counts)}
    (group_output / "exposure_qc.json").write_text(json.dumps(qc, indent=2) + "\n")
    return qc


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("group_output", type=Path)
    args = parser.parse_args()
    print(json.dumps(normalize(args.group_output), indent=2))


if __name__ == "__main__":
    main()
