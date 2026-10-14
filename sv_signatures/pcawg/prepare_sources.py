#!/usr/bin/env python3
"""Fetch immutable third-party source revisions without altering existing checkouts."""
import json
import argparse
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[2]
lock = json.loads((root / 'sv_signatures/pcawg/upstream.json').read_text())
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--clustersv-only',action='store_true',help='Fetch only the supported clustering source')
args = parser.parse_args()
sources = [('clustersv', root / 'resources/clustersv/source')]
if not args.clustersv_only:
    sources.append(('panprostate', root / 'resources/pcawg_signatures/panprostate'))
for name, destination in sources:
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
