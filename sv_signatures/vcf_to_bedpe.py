#!/usr/bin/env python3
"""Convert one caller's somatic VCF to auditable SigProfiler SV BEDPE input.

Each retained row represents one two-ended adjacency or symbolic SV record.
Reciprocal BND mates count once. This does not reconstruct multi-junction events.
"""

import argparse
import csv
import gzip
import json
import re
from collections import Counter
from pathlib import Path

ALT_BND = re.compile(r"[\[\]]([^\[\]]+):(\d+)[\[\]]")
CHROMS = {str(i) for i in range(1, 23)} | {"X", "Y"}
BEDPE_FIELDS = ("sample", "chrom1", "start1", "end1", "chrom2", "start2", "end2", "svclass")
AUDIT_FIELDS = ("record_id", "mate_id", "chrom1", "pos1", "chrom2", "pos2", "native_svtype", "svclass", "decision")
CLUSTERSV_FIELDS = ("chrom1", "start1", "end1", "chrom2", "start2", "end2", "name", "score", "strand1", "strand2")
CLASSES = {"DEL": "deletion", "DUP": "tandem-duplication", "INV": "inversion", "TRA": "translocation",
           "deletion": "deletion", "tandem-duplication": "tandem-duplication", "inversion": "inversion", "translocation": "translocation"}


def chrom_name(name):
    name = name.removeprefix("chr")
    return name if name in CHROMS else None


def info_map(value):
    return dict(item.split("=", 1) if "=" in item else (item, True) for item in value.split(";"))


def bnd_endpoint(chrom, pos, alt):
    match = ALT_BND.search(alt)
    if not match:
        return None
    other, other_pos = match.group(1), int(match.group(2))
    orientation = ("-" if alt[0] in "[]" else "+") + ("+" if match.group(0)[0] == "]" else "-")
    if (other, other_pos) < (chrom, pos):
        return other, other_pos, chrom, pos, orientation[::-1]
    return chrom, pos, other, other_pos, orientation


def mate_position_agrees(record, mate, tolerance=400):
    match = ALT_BND.search(record[4])
    return bool(match and match.group(1) == mate[0]
                and abs(int(match.group(2)) - int(mate[1])) <= tolerance)


def classify(record):
    chrom, pos, ident, _ref, alt, _qual, _filter, raw_info = record[:8]
    pos = int(pos)
    info = info_map(raw_info)
    native = info.get("SVTYPE", "")
    if native == "BND":
        endpoint = bnd_endpoint(chrom, pos, alt)
        if endpoint is None and info.get("CHR2") and info.get("POS2"):
            endpoint = (chrom, pos, info["CHR2"], int(info["POS2"]), None)
        if endpoint is None:
            return None, "unresolved_breakend"
        c1, p1, c2, p2, orientation = endpoint
        if c1 != c2:
            kind = "translocation"
        elif orientation in ("++", "--"):
            kind = "inversion"
        elif orientation == "+-":
            kind = "deletion"
        elif orientation == "-+":
            kind = "tandem-duplication"
        else:
            return None, "unknown_breakend_orientation"
    elif native in ("DEL", "DUP", "INV", "TRA") and info.get("END"):
        c1, p1, c2, p2 = chrom, pos, info.get("CHR2", chrom), int(info["END"])
        if (c2, p2) < (c1, p1):
            c1, p1, c2, p2 = c2, p2, c1, p1
        kind = CLASSES[native]
    else:
        return None, "unsupported_svtype_or_endpoint"
    c1, c2 = chrom_name(c1), chrom_name(c2)
    if c1 is None or c2 is None:
        return None, "nonstandard_contig"
    if p1 < 1 or p2 < 1:
        return None, "invalid_coordinate"
    if c1 == c2 and p1 == p2:
        return None, "zero_span"
    label = info.get("SVCLASS")
    if label and label in CLASSES and CLASSES[label] != kind:
        return None, "brass_class_orientation_disagreement"
    if label and label not in CLASSES:
        return None, "unsupported_brass_class"
    return (c1, p1, c2, p2, kind, native), "retained"


def oriented_strands(record):
    """Return ordered endpoint strands; never invent an inversion orientation."""
    info = info_map(record[7])
    if info.get("SVTYPE") == "BND":
        endpoint = bnd_endpoint(record[0], int(record[1]), record[4])
        if endpoint:
            return endpoint[4]
    if info.get("SVTYPE") == "DEL":
        strands = "+-"
    elif info.get("SVTYPE") == "DUP":
        strands = "-+"
    else:
        options = set(re.findall(r"([+-]{2})(?=:|$|,)", str(info.get("STRANDS", ""))))
        ct = {"3to3": "++", "5to5": "--", "3to5": "+-", "5to3": "-+"}.get(info.get("CT"))
        if ct:
            options.add(ct)
        if info.get("INV3"):
            options.add("++")
        if info.get("INV5"):
            options.add("--")
        if len(options) != 1:
            return None
        strands = options.pop()
    # DELLY's CHR2/POS2 fallback preserves raw endpoint order in classify().
    if info.get("SVTYPE") == "BND":
        return strands
    second = (info.get("CHR2", record[0]), int(info.get("END", info.get("POS2", record[1]))))
    return strands[::-1] if second < (record[0], int(record[1])) else strands


def convert(vcf, caller, sample, bedpe, audit, summary, min_size=1000, output_format="sigprofiler"):
    if output_format not in ("sigprofiler", "clustersv"):
        raise ValueError("Unknown BEDPE output format")
    if caller not in ("manta", "delly", "gridss", "svaba", "brass"):
        raise ValueError("Unknown caller")
    with vcf.open("rb") as handle:
        compressed = handle.read(2) == b"\x1f\x8b"
    opener = gzip.open if compressed else open
    records = {}
    with opener(vcf, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8 or not fields[2] or fields[2] == "." or fields[2] in records:
                raise ValueError("Malformed, missing, or duplicate VCF ID in " + str(vcf))
            records[fields[2]] = fields
    passed = {key: fields for key, fields in records.items()
              if fields[6] == ("." if caller == "brass" else "PASS")}
    output, checks, decisions = [], [], Counter()
    seen = set()
    for key, fields in passed.items():
        if key in seen:
            continue
        info = info_map(fields[7])
        mate_id = info.get("MATEID", "")
        native = info.get("SVTYPE", "")
        if native == "BND" and mate_id:
            mate = passed.get(mate_id)
            if not mate or info_map(mate[7]).get("MATEID") != key:
                event, decision = None, "missing_or_nonreciprocal_pass_mate"
            else:
                if not mate_position_agrees(fields, mate) or not mate_position_agrees(mate, fields):
                    event, decision = None, "mate_endpoint_disagreement"
                else:
                    event, decision = classify(fields)
                    if event:
                        # POS from the two records is authoritative. ALT remote
                        # coordinates may differ slightly within microhomology.
                        a, b = sorted(((chrom_name(fields[0]), int(fields[1])),
                                       (chrom_name(mate[0]), int(mate[1]))))
                        event = (a[0], a[1], b[0], b[1], event[4], event[5])
                seen.add(mate_id)
        elif native == "BND" and caller != "delly":
            event, decision = None, "unpaired_breakend"
        else:
            event, decision = classify(fields)
        seen.add(key)
        c1 = p1 = c2 = p2 = kind = ""
        if event:
            c1, p1, c2, p2, kind, _native = event
            if c1 == c2 and abs(p2 - p1) < min_size:
                decision = "below_min_size"
            else:
                if output_format == "clustersv":
                    strands = oriented_strands(fields)
                    if strands is None:
                        decision = "unresolved_orientation"
                    else:
                        output.append(dict(chrom1=c1, start1=p1 - 1, end1=p1,
                                           chrom2=c2, start2=p2 - 1, end2=p2, name=key,
                                           score=0, strand1=strands[0], strand2=strands[1]))
                else:
                    output.append(dict(sample=sample, chrom1=c1, start1=p1 - 1, end1=p1,
                                       chrom2=c2, start2=p2 - 1, end2=p2, svclass=kind))
        decisions[decision] += 1
        checks.append(dict(record_id=key, mate_id=mate_id, chrom1=c1, pos1=p1,
                           chrom2=c2, pos2=p2, native_svtype=native, svclass=kind, decision=decision))
    bedpe.parent.mkdir(parents=True, exist_ok=True)
    with bedpe.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CLUSTERSV_FIELDS if output_format == "clustersv" else BEDPE_FIELDS, delimiter="\t")
        if output_format == "sigprofiler":
            writer.writeheader()
        writer.writerows(output)
    with audit.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=AUDIT_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(checks)
    summary.write_text(json.dumps({"sample": sample, "caller": caller, "source_vcf": str(vcf.resolve()),
                                   "total_records": len(records), "selected_records": len(passed),
                                   "collapsed_selected_events": len(checks), "retained_events": len(output),
                                   "min_intrachromosomal_size_bp": min_size, "output_format": output_format,
                                   "decisions": dict(decisions)},
                                  indent=2) + "\n")
    return len(output)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vcf", required=True, type=Path)
    parser.add_argument("--caller", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--min-size-bp", type=int, default=1000)
    args = parser.parse_args()
    if args.min_size_bp < 1:
        parser.error("--min-size-bp must be positive")
    args.outdir.mkdir(parents=True, exist_ok=True)
    count = convert(args.vcf, args.caller, args.sample, args.outdir / (args.sample + ".bedpe.tsv"),
                    args.outdir / "event_audit.tsv", args.outdir / "summary.json", args.min_size_bp)
    print("Retained {} two-ended events for {}".format(count, args.sample))


if __name__ == "__main__":
    main()
