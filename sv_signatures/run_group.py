#!/usr/bin/env python3
"""Build one cohort/caller SV32 matrix, then refit COSMIC signatures."""

import argparse
import csv
import hashlib
import json
import sys
import time
from importlib.metadata import version
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))

from vcf_to_bedpe import convert
from normalize_exposures import normalize


def install_sigprofiler_sv_compatibility():
    """Guard tested pandas semantics and bridge two removed upstream APIs."""
    import numpy as np
    import pandas as pd

    if version("SigProfilerMatrixGenerator") != "1.3.6":
        raise RuntimeError("SV compatibility layer is reviewed only for SigProfilerMatrixGenerator 1.3.6")
    if pd.__version__ != "2.2.3" or pd.options.mode.copy_on_write:
        raise RuntimeError("SV clustering is validated only with pandas 2.2.3 and Copy-on-Write disabled")
    shims = []
    if not hasattr(np, "in1d"):
        # np.isin has the same membership semantics/defaults used upstream.
        np.in1d = np.isin
        shims.append("numpy.in1d=numpy.isin")
    if not hasattr(pd.DataFrame, "append"):
        # Upstream invokes DataFrame.append(other) with default options only.
        def append_default(self, other):
            return pd.concat([self, other], ignore_index=False, sort=False)
        pd.DataFrame.append = append_default
        shims.append("pandas.DataFrame.append=pandas.concat(defaults)")
    return shims


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def matrix_candidate(root):
    candidates = []
    for path in root.rglob("*"):
        if not path.is_file() or "SV32" not in path.name or path.stat().st_size == 0:
            continue
        try:
            with path.open(newline="") as handle:
                reader = csv.reader(handle, delimiter="\t")
                header = next(reader)
                rows = list(reader)
            if len(rows) == 32 and len(header) >= 2 and all(len(row) == len(header) for row in rows):
                candidates.append((path, header, rows))
        except (UnicodeError, csv.Error, StopIteration):
            continue
    if len(candidates) != 1:
        raise ValueError("Expected one native SV32 matrix, found: " + repr([str(x[0]) for x in candidates]))
    return candidates[0]


def normalize_matrix(source, header, rows, samples, output, counts):
    if len(set(header[1:])) != len(header) - 1:
        raise ValueError("Native matrix has duplicate sample columns")
    extra = set(header[1:]) - set(samples)
    if extra:
        raise ValueError("Native matrix has unexpected samples: " + repr(sorted(extra)))
    lookup = {name: index for index, name in enumerate(header)}
    for name in samples:
        observed = sum(int(row[lookup[name]]) for row in rows) if name in lookup else 0
        if observed != counts[name]:
            raise ValueError("SV32 count mismatch for {}: generator {}, converter {}".format(name, observed, counts[name]))
    with output.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([header[0]] + samples)
        for row in rows:
            writer.writerow([row[0]] + [row[lookup[name]] if name in lookup else "0" for name in samples])
    clustered = sum(sum(int(row[lookup[name]]) for name in header[1:])
                    for row in rows if row[0].startswith("clustered_"))
    nonclustered = sum(sum(int(row[lookup[name]]) for name in header[1:])
                       for row in rows if row[0].startswith("non-clustered_"))
    if clustered + nonclustered != sum(counts.values()):
        raise ValueError("SV32 clustered/non-clustered totals do not match source events")
    return {"native_matrix": str(source.relative_to(output.parent)), "samples": len(samples),
            "zero_event_samples": [name for name in samples if counts[name] == 0],
            "channels": len(rows), "total_events": sum(counts.values()),
            "clustered_events": clustered, "nonclustered_events": nonclustered}


def run(args):
    out = args.output.resolve()
    if out.exists() and any(out.iterdir()):
        raise ValueError("Output directory is not empty: " + str(out))
    out.mkdir(parents=True, exist_ok=True)
    with args.manifest.open(newline="") as handle:
        all_rows = list(csv.DictReader(handle, delimiter="\t"))
    rows = [r for r in all_rows if r["cohort"] == args.cohort and r["caller"] == args.caller
            and (not args.pair_id or r["pair_id"] == args.pair_id)
            and (not args.pair_ids or r["pair_id"] in args.pair_ids.split(','))]
    if len(rows) != args.expected_pairs or len({r["pair_id"] for r in rows}) != len(rows):
        raise ValueError("Expected {} unique {} {} pairs, found {}".format(args.expected_pairs, args.cohort, args.caller, len(rows)))
    inputs = out / "bedpe"
    audits = out / "audit"
    inputs.mkdir()
    audits.mkdir()
    source_rows, counts = [], {}
    for row in sorted(rows, key=lambda r: r["pair_id"]):
        pair, source = row["pair_id"], Path(row["vcf"])
        if not source.is_file():
            raise ValueError("Missing source VCF: " + str(source))
        digest = sha256(source)
        if row["source_sha256"] and row["source_sha256"] != digest:
            raise ValueError("Source VCF SHA256 changed: " + str(source))
        if getattr(args,'genome_build',None):
            from reference_policy import check_vcf
            check_vcf(source,args.genome_build)
        counts[pair] = convert(source, args.caller, pair, inputs / (pair + ".bedpe"),
                               audits / (pair + ".audit.tsv"), audits / (pair + ".summary.json"), args.min_size_bp)
        source_rows.append({"pair_id": pair, "vcf": str(source.resolve()),
                            "sha256": digest, "source_status": row["source_status"], "retained_events": counts[pair]})
    with (out / "sources.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=source_rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(source_rows)
    start = time.time()
    compatibility_shims = install_sigprofiler_sv_compatibility()
    from SigProfilerMatrixGenerator.scripts import SVMatrixGenerator as sv
    native = out / "sigprofiler_native"
    native.mkdir()
    sv.generateSVMatrix(str(inputs), "{}_{}".format(args.cohort, args.caller), str(native))
    source, header, matrix_rows = matrix_candidate(native)
    matrix = out / "SV32.matrix.tsv"
    matrix_info = normalize_matrix(source, header, matrix_rows, sorted(counts), matrix, counts)
    matrix_info.update({"cohort": args.cohort, "caller": args.caller,
                        "generator_seconds": time.time() - start,
                        "allocated_cpus": args.cpus,
                        "compatibility_shims": compatibility_shims,
                        "min_size_bp": args.min_size_bp,
                        "source_statuses": sorted({r["source_status"] for r in source_rows})})
    (out / "matrix_qc.json").write_text(json.dumps(matrix_info, indent=2) + "\n")
    if args.skip_fit:
        return
    if not sum(counts.values()):
        (out / "fit_status.json").write_text(json.dumps({"status": "no_events", "cosmic_version": args.cosmic_version}, indent=2) + "\n")
        return
    from SigProfilerAssignment import Analyzer as analyze
    fit = out / "cosmic_fit"
    start = time.time()
    analyze.cosmic_fit(str(matrix), str(fit), input_type="matrix", cosmic_version=args.cosmic_version,
                       genome_build="GRCh37", collapse_to_SBS96=False,
                       make_plots=False, verbose=False, cpu=args.cpus)
    exposure_qc = normalize(out)
    (out / "fit_status.json").write_text(json.dumps({"status": "completed", "cosmic_version": args.cosmic_version,
                                                     "assigned_cpus": args.cpus,
                                                     "standardized_samples": exposure_qc["standardized_activity_rows"],
                                                     "elapsed_seconds": time.time() - start}, indent=2) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--cohort", required=True)
    parser.add_argument("--caller", choices=("manta", "delly", "gridss", "svaba", "brass"), required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--expected-pairs", type=int, required=True)
    parser.add_argument("--min-size-bp", type=int, default=1000)
    parser.add_argument("--cosmic-version", type=float, default=3.6)
    parser.add_argument("--genome-build",choices=('GRCh37',),default='GRCh37')
    parser.add_argument("--cpus", type=int, default=1)
    parser.add_argument("--pair-id", default="", help="One-sample technical smoke test only")
    parser.add_argument("--pair-ids", default="", help="Comma-separated pilot sample IDs")
    parser.add_argument("--skip-fit", action="store_true")
    args = parser.parse_args()
    if args.expected_pairs < 1 or args.min_size_bp < 1 or args.cpus < 1:
        parser.error("expected-pairs, min-size-bp, and cpus must be positive")
    run(args)


if __name__ == "__main__":
    main()
