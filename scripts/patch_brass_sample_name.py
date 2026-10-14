#!/usr/bin/env python3
"""Patch only the initial BRASS filter's sample lookup, not output filenames.

The groups header preserves BAM SM verbatim. The native filename-safe alias
must not be used to look up that header. Fail closed on unexpected source.
"""
import hashlib
import json
from pathlib import Path
import sys


def patch(source):
    old = "$options->{'safe_tumour_name'},\n                                      $groups,"
    new = "$options->{'tumour_name'},\n                                      $groups,"
    if source.count(old) != 1:
        raise ValueError('Expected exactly one initial BRASS filter sample argument')
    return source.replace(old, new)


if __name__ == '__main__':
    source, destination = map(Path, sys.argv[1:])
    original = source.read_text()
    updated = patch(original)
    destination.write_text(updated)
    destination.with_suffix('.patch.json').write_text(json.dumps({
        'change': 'initial filter: match original BAM SM against groups header',
        'original_sha256': hashlib.sha256(original.encode()).hexdigest(),
        'patched_sha256': hashlib.sha256(updated.encode()).hexdigest(),
        'scientific_threshold_changes': False,
    }, indent=2) + '\n')
