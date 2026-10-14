#!/usr/bin/env python3
"""Fetch immutable third-party source revisions without altering existing checkouts."""
import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[2]
lock = json.loads((root / 'sv_signatures/pcawg/upstream.json').read_text())
for name, destination in [('clustersv', root / 'resources/clustersv/source'),
                          ('panprostate', root / 'resources/pcawg_signatures/panprostate')]:
    info = lock[name]
    if destination.exists():
        head = subprocess.check_output(['git', '-C', str(destination), 'rev-parse', 'HEAD'], text=True).strip()
        dirty = subprocess.check_output(['git', '-C', str(destination), 'status', '--porcelain'], text=True).strip()
        if head != info['commit'] or dirty:
            raise SystemExit('Existing source requires review; refusing to modify ' + str(destination))
    else:
        destination.parent.mkdir(parents=True, exist_ok=True)
        subprocess.run(['git', 'clone', '--no-checkout', info['url'] + '.git', str(destination)], check=True)
        subprocess.run(['git', '-C', str(destination), 'checkout', '--detach', info['commit']], check=True)
    print(name, info['commit'])
