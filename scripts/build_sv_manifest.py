#!/usr/bin/env python3
"""Build a caller-independent signature manifest from a private callset inventory."""
import argparse
import csv
import hashlib
from pathlib import Path
import re

FIELDS = ('cohort','pair_id','caller','vcf','source_status','source_sha256')
CALLERS = {'manta','delly','gridss','svaba','brass'}

def validate_rows(rows, hash_sources=True):
    """Shared validation for discovered pipeline outputs and external inventories."""
    if not rows:
        raise ValueError('Empty inventory')
    seen = set()
    output = []
    for row in rows:
        key = tuple(row.get(k) for k in FIELDS[:3])
        if key in seen:
            raise ValueError('Duplicate cohort/pair/caller: ' + repr(key))
        seen.add(key)
        if any(not re.fullmatch(r'[A-Za-z0-9_-]+', row[k] or '') for k in ('cohort','pair_id')) or row['caller'] not in CALLERS:
            raise ValueError('Invalid cohort, pair or caller')
        if not row['source_status']:
            raise ValueError('QC status is required; hashing is not biological validation')
        path = Path(row['vcf']).resolve(strict=True)
        if not path.is_file() or not path.stat().st_size:
            raise ValueError('Missing or empty VCF: ' + str(path))
        checksum = ''
        if hash_sources:
            h = hashlib.sha256()
            with path.open('rb') as data:
                for chunk in iter(lambda: data.read(1024*1024), b''):
                    h.update(chunk)
            checksum = h.hexdigest()
        output.append(dict(zip(FIELDS, key + (str(path),row['source_status'],checksum))))
    return output


def build(source, destination):
    with Path(source).open(newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        required = set(FIELDS[:-1])
        absent = required - set(reader.fieldnames or [])
        if absent:
            raise ValueError('Missing inventory columns: ' + ','.join(sorted(absent)))
        rows = list(reader)
    output = validate_rows(rows)
    with Path(destination).open('x', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter='\t')
        writer.writeheader(); writer.writerows(output)
    return len(output)

if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--inventory', required=True, type=Path)
    ap.add_argument('--output', required=True, type=Path)
    args = ap.parse_args()
    print('Wrote', build(args.inventory,args.output), 'callsets')
