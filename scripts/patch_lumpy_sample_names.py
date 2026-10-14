#!/usr/bin/env python3
"""Use per-input safe IDs in LUMPY's colon-delimited option strings.

Caller output must be renamed back to the original BAM SM before SVTyper.
No evidence extraction or scientific parameters are changed.
"""
from pathlib import Path
import hashlib
import json
import sys
import argparse


def patch(source, keep_temp=False):
    for variable in ('DISC_SAMPLE', 'SPL_SAMPLE'):
        old = 'id:${' + variable + '},'
        if source.count(old) != 2:
            raise ValueError('Expected two ID arguments for ' + variable)
        source = source.replace(old, 'id:LUMPY_SAMPLE_${i},')
    if keep_temp:
        old = 'trap cleanup EXIT'
        if source.count(old) != 1:
            raise ValueError('Expected one native cleanup trap')
        source = source.replace(old, '# Cleanup owned by wrapper after output validation')
    return source


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('source', type=Path)
    parser.add_argument('destination', type=Path)
    parser.add_argument('--keep-temp', action='store_true')
    args = parser.parse_args()
    source, destination = args.source, args.destination
    original = source.read_text()
    updated = patch(original, keep_temp=args.keep_temp)
    destination.write_text(updated)
    destination.chmod(0o755)
    destination.with_suffix('.patch.json').write_text(json.dumps({
        'change': 'safe per-input IDs; restore original SM before genotyping',
        'original_sha256': hashlib.sha256(original.encode()).hexdigest(),
        'patched_sha256': hashlib.sha256(updated.encode()).hexdigest(),
        'scientific_threshold_changes': False,
        'wrapper_owns_cleanup': args.keep_temp,
    }, indent=2) + '\n')
