#!/usr/bin/env python3
"""Build deterministic tumour-normal manifests from EGA sample/file metadata."""

import argparse
import csv
import json
import re
from collections import defaultdict
from pathlib import Path
from typing import Dict, List


UUID_BAM = re.compile(
    r"([0-9a-f]{8}(?:-[0-9a-f]{4}){3}-[0-9a-f]{12}\.bam)$", re.IGNORECASE
)


def read_csv(path: Path) -> List[Dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def donor_id(sample: Dict[str, str]) -> str:
    attributes = json.loads(sample["extra_attributes"])
    values = [item.get("value", "") for item in attributes if item.get("tag") == "donor_id"]
    if len(values) != 1 or not values[0]:
        raise ValueError(f"Expected one donor_id for {sample['accession_id']}")
    return values[0]


def sample_role(phenotype: str) -> str:
    value = phenotype.lower()
    if "normal" in value:
        return "normal"
    if "tumour" in value or "tumor" in value:
        return "tumor"
    raise ValueError(f"Cannot classify phenotype: {phenotype!r}")


def write_tsv(path: Path, rows: List[Dict[str, str]], fields: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--metadata-dir", required=True, type=Path)
    parser.add_argument("--bam-dir", required=True, type=Path)
    parser.add_argument("--dataset-id", required=True)
    parser.add_argument("--purity-table", type=Path,
                        help="Optional CSV keyed by sample_alias; missing values stay blank")
    parser.add_argument("--pairs", required=True, type=Path)
    parser.add_argument("--runnable-pairs", type=Path,
                        help="Optional sheet containing only pairs with both BAM indexes")
    parser.add_argument("--unpaired", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    args = parser.parse_args()

    samples = read_csv(args.metadata_dir / "samples.csv")
    sample_files = read_csv(args.metadata_dir / "sample_file.csv")
    purity_rows = read_csv(args.purity_table) if args.purity_table else []

    samples_by_accession = {row["accession_id"]: row for row in samples}
    if len(samples_by_accession) != len(samples):
        raise ValueError("Duplicate sample accession in metadata")
    purity_by_alias = {row["sample_alias"]: row for row in purity_rows}
    if len(purity_by_alias) != len(purity_rows):
        raise ValueError("Duplicate sample alias in purity table")
    file_records = {}

    for row in sample_files:
        match = UUID_BAM.search(row["file_name"])
        if not match:
            raise ValueError(f"Could not extract BAM UUID from {row['file_name']!r}")
        local_name = match.group(1)
        if local_name in file_records:
            raise ValueError(f"Duplicate file mapping for {local_name}")
        file_records[local_name] = row

    grouped = defaultdict(dict)
    unmapped = []

    for bam in sorted(args.bam_dir.glob("*.bam")):
        file_record = file_records.get(bam.name)
        if file_record is None:
            unmapped.append(bam.name)
            continue
        sample = samples_by_accession.get(file_record["sample_accession_id"])
        if sample is None:
            raise ValueError(f"File mapping has unknown sample accession: {file_record['sample_accession_id']}")
        role = sample_role(sample["phenotype"])
        patient = donor_id(sample)
        if role in grouped[patient]:
            raise ValueError(f"Multiple local {role} samples for {patient}")
        grouped[patient][role] = {
            "bam": str(bam.resolve()),
            "bai": str(Path(f"{bam.resolve()}.bai")),
            "sample_accession": sample["accession_id"],
            "sample_alias": sample["alias"],
            "file_accession": file_record["file_accession_id"],
            "phenotype": sample["phenotype"],
            "sex": sample["biological_sex"],
            "subject_id": sample["subject_id"],
        }

    if unmapped:
        raise ValueError(f"Local BAMs absent from metadata: {', '.join(unmapped)}")

    pair_rows = []
    unpaired_rows = []
    missing_indexes = []

    for patient, records in sorted(grouped.items()):
        for record in records.values():
            if not Path(record["bai"]).is_file():
                missing_indexes.append(record["bai"])

        if set(records) != {"tumor", "normal"}:
            for role, record in sorted(records.items()):
                unpaired_rows.append(
                    {
                        "patient_id": patient,
                        "role": role,
                        "phenotype": record["phenotype"],
                        "sample_accession": record["sample_accession"],
                        "sample_alias": record["sample_alias"],
                        "file_accession": record["file_accession"],
                        "bam": record["bam"],
                        "bai": record["bai"],
                    }
                )
            continue

        tumor = records["tumor"]
        normal = records["normal"]
        if tumor["subject_id"] != normal["subject_id"]:
            raise ValueError(f"Subject mismatch within {patient}")
        if tumor["sex"] != normal["sex"]:
            raise ValueError(f"Sex mismatch within {patient}")
        purity = purity_by_alias.get(tumor["sample_alias"], {})
        pair_rows.append(
            {
                "patient_id": patient,
                "subject_id": tumor["subject_id"],
                "sex": tumor["sex"],
                "tumor_phenotype": tumor["phenotype"],
                "tumor_sample_accession": tumor["sample_accession"],
                "tumor_sample_alias": tumor["sample_alias"],
                "tumor_file_accession": tumor["file_accession"],
                "tumor_bam": tumor["bam"],
                "tumor_bai": tumor["bai"],
                "normal_sample_accession": normal["sample_accession"],
                "normal_sample_alias": normal["sample_alias"],
                "normal_file_accession": normal["file_accession"],
                "normal_bam": normal["bam"],
                "normal_bai": normal["bai"],
                "purity": purity.get("purity", ""),
                "ploidy": purity.get("ploidy", ""),
                "purity_conf_mad": purity.get("purity_conf_mad", ""),
                "wgd_status": purity.get("wgd_status", ""),
                "wgd_uncertain": purity.get("wgd_uncertain", ""),
            }
        )

    pair_fields = [
        "patient_id", "subject_id", "sex", "tumor_phenotype",
        "tumor_sample_accession", "tumor_sample_alias", "tumor_file_accession",
        "tumor_bam", "tumor_bai", "normal_sample_accession", "normal_sample_alias",
        "normal_file_accession", "normal_bam", "normal_bai", "purity", "ploidy",
        "purity_conf_mad", "wgd_status", "wgd_uncertain",
    ]
    unpaired_fields = [
        "patient_id", "role", "phenotype", "sample_accession", "sample_alias",
        "file_accession", "bam", "bai",
    ]
    write_tsv(args.pairs, pair_rows, pair_fields)
    runnable_rows = [row for row in pair_rows
                     if Path(row["tumor_bai"]).is_file() and Path(row["normal_bai"]).is_file()]
    if args.runnable_pairs:
        write_tsv(args.runnable_pairs, runnable_rows, pair_fields)
    write_tsv(args.unpaired, unpaired_rows, unpaired_fields)

    summary = {
        "dataset": args.dataset_id,
        "metadata_samples": len(samples),
        "local_bams": sum(len(records) for records in grouped.values()),
        "local_patients": len(grouped),
        "complete_pairs": len(pair_rows),
        "runnable_pairs_with_indexes": len(runnable_rows),
        "pairs_missing_indexes": sorted({row["patient_id"] for row in pair_rows
                                         if row not in runnable_rows}),
        "unpaired_samples": len(unpaired_rows),
        "pairs_with_purity": sum(bool(row["purity"]) for row in pair_rows),
        "missing_indexes": missing_indexes,
    }
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    args.summary.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))

    if missing_indexes and not args.runnable_pairs:
        raise SystemExit("Manifest written, but one or more BAM indexes are missing")


if __name__ == "__main__":
    main()
