#!/usr/bin/env python3
"""Recover BRASS's final archive warning only after checking the finished callset.

This is deliberately narrow: it does not rerun or change SV calling, and refuses
to recover any task without one readable, indexed native VCF and BEDPE.
"""

import hashlib
import json
import subprocess
import sys
import tarfile
from pathlib import Path


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def inventory(source):
    result = {}
    for path in sorted(source.rglob("*")):
        if path.is_symlink():
            raise ValueError("Symlink in BRASS intermediates: " + str(path))
        if path.is_file():
            result[str(path.relative_to(source.parent))] = (path.stat().st_size, digest(path))
        elif not path.is_dir():
            raise ValueError("Unexpected BRASS intermediate: " + str(path))
    if not result:
        raise ValueError("No BRASS intermediates to archive")
    return result


def recover(calls, report):
    calls = Path(calls).resolve()
    vcfs = list(calls.glob("*.annot.vcf.gz"))
    bedpes = list(calls.glob("*.annot.bedpe.gz"))
    if len(vcfs) != 1 or len(bedpes) != 1:
        raise ValueError("Expected exactly one annotated VCF and BEDPE")
    protected = vcfs + bedpes + [Path(str(p) + ".tbi") for p in vcfs + bedpes]
    if any(not p.is_file() or p.stat().st_size == 0 for p in protected):
        raise ValueError("Missing or empty BRASS callset/index")
    before = {str(p): digest(p) for p in protected}
    subprocess.run(["bcftools", "view", "-Ov", "-o", "/dev/null", str(vcfs[0])], check=True)
    count = int(subprocess.check_output(["bcftools", "index", "-n", str(vcfs[0])],
                                        universal_newlines=True))
    if count < 1:
        raise ValueError("Empty BRASS VCF; manual review required")
    image = Path(__file__).resolve().parents[1] / "containers/images/brass-v6.3.4.sif"
    subprocess.run(["singularity", "exec", "--cleanenv", "--bind",
                    str(calls) + ":" + str(calls) + ":ro", str(image),
                    "tabix", "-l", str(bedpes[0])],
                   stdout=subprocess.DEVNULL, check=True)
    source = calls / "intermediates"
    if not source.is_dir():
        raise ValueError("Missing BRASS intermediates")
    original = inventory(source)
    archive = calls / "intermediates.recovered.tar.gz"
    if archive.exists():
        raise ValueError("Refusing to overwrite existing recovery archive")
    temporary = calls / "intermediates.recovered.tar.gz.partial"
    if temporary.exists():
        raise ValueError("Refusing to overwrite partial recovery archive")
    try:
        subprocess.run(["tar", "-C", str(calls), "-zcf", str(temporary), "intermediates"], check=True)
        if inventory(source) != original:
            raise ValueError("BRASS intermediates changed while archiving")
        observed = {}
        with tarfile.open(temporary, "r:gz") as handle:
            for member in handle:
                if member.isdir():
                    continue
                if not member.isfile() or member.name in observed:
                    raise ValueError("Unexpected or duplicate archive member")
                h = hashlib.sha256()
                with handle.extractfile(member) as stream:
                    for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b""):
                        h.update(chunk)
                observed[member.name] = (member.size, h.hexdigest())
        if observed != original:
            raise ValueError("Recovered archive differs from source intermediates")
        if {str(p): digest(p) for p in protected} != before:
            raise ValueError("BRASS callset changed during archive recovery")
        temporary.rename(archive)
    except Exception:
        # Retain the partial archive for diagnosis; never claim success.
        raise
    report = Path(report)
    report.write_text(json.dumps({"status": "archive_recovered_and_verified",
                                  "reason": "native tar reported intermediates changed during read",
                                  "records": count, "archive": str(archive),
                                  "archive_sha256": digest(archive),
                                  "callset_sha256": before,
                                  "intermediate_file_count": len(original)}, indent=2) + "\n")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit("Usage: recover_brass_archive.py CALLS_DIR REPORT_JSON")
    recover(sys.argv[1], sys.argv[2])
