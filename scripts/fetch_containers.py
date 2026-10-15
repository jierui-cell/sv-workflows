#!/usr/bin/env python3
"""Pull digest-pinned SIF images and record their local hashes; never overwrite an image."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import yaml
from workflow import resolve_config
from hpc_environment import load_modules

ROOT = Path(__file__).resolve().parents[1]
if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--engine', choices=('apptainer','singularity'))
    ap.add_argument('--config', type=Path, action='append', default=[])
    ap.add_argument('--tools', default='manta,delly,svaba,gridss,brass')
    args = ap.parse_args()
    default_config = ROOT/'config/private-cohort.yaml'
    paths = args.config or ([default_config] if default_config.is_file() else [])
    try:
        config = resolve_config('callers',paths,root=ROOT)
        load_modules(config.get('modules'))
    except (ValueError,OSError,yaml.YAMLError) as exc:
        ap.error(str(exc))
    args.engine = args.engine or config['container_engine']
    rows = list(csv.DictReader((ROOT/'config/containers.tsv').open(),delimiter='\t'))
    wanted = args.tools.split(',')
    if set(wanted)-{r['tool'] for r in rows}: ap.error('Unknown tool')
    folder = ROOT/'containers/images'; folder.mkdir(parents=True,exist_ok=True)
    records = []
    for row in rows:
        if row['tool'] not in wanted: continue
        if '@sha256:' not in row['uri']: raise ValueError('Unpinned container URI')
        image = folder/(row['tool']+'-'+row['tag']+'.sif')
        receipt = Path(str(image)+'.json')
        if image.exists() and not receipt.exists():
            raise ValueError('Existing image has no receipt; review manually: '+str(image))
        if not image.exists(): subprocess.run([args.engine,'pull',str(image),row['uri']],check=True)
        h=hashlib.sha256()
        with image.open('rb') as data:
            for block in iter(lambda:data.read(1048576),b''): h.update(block)
        record=dict(row,path=str(image),sha256=h.hexdigest())
        if receipt.exists() and json.loads(receipt.read_text()) != record:
            raise ValueError('Image changed since download: '+str(image))
        subprocess.run([args.engine,'inspect',str(image)],check=True,stdout=subprocess.DEVNULL)
        receipt.write_text(json.dumps(record,indent=2)+'\n')
        records.append(record)
    print(json.dumps(records,indent=2))
