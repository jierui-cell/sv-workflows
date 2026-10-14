#!/usr/bin/env python3
"""Review candidate source files for accidental private data; not a secret-scanner guarantee."""
import argparse
import json
from pathlib import Path
import re
import subprocess
ROOT=Path(__file__).resolve().parents[1]
EXCLUDE={'.git','.venv','__pycache__','inputs','results','work','logs','resources','environments','containers','site','.nextflow','.nextflow-home'}
SENSITIVE_SUFFIXES={'.bam','.bai','.cram','.crai','.sif','.bin','.bigWig','.pt','.pem','.key'}
def candidates():
    if (ROOT/'.git').exists():
        for name in subprocess.check_output(['git','-C',str(ROOT),'ls-files'],text=True).splitlines():
            yield ROOT/name
    else:
        for p in ROOT.rglob('*'):
            rel=p.relative_to(ROOT)
            if p.is_file() and not any(x in EXCLUDE for x in rel.parts) and not p.name.startswith('private'):
                yield p
def audit():
    issues=[]; count=0
    patterns=[r'/(?:scratch|data\d+)/[^/\s]+/[^/\s]+/(?:EGA_data|data/EGA)',r'\b(?:MELA|PNET)_\d{4}\b',
              r'gh[pousr]_[A-Za-z0-9]{20,}',r'github_pat_[A-Za-z0-9_]{20,}',r'AKIA[A-Z0-9]{16}',
              r'-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----']
    for path in candidates():
        rel=str(path.relative_to(ROOT)); count+=1
        if any(x in EXCLUDE for x in path.relative_to(ROOT).parts) or path.name.startswith('private') or path.suffix in SENSITIVE_SUFFIXES:
            issues.append(rel+': excluded runtime/private artifact')
        if '.vcf' in path.name and rel!='examples/toy.somatic.vcf':
            issues.append(rel+': non-allowlisted genomic file')
        try: text=path.read_text()
        except (UnicodeError,OSError):
            issues.append(rel+': binary/unreadable file needs review'); continue
        if any(re.search(p,text) for p in patterns): issues.append(rel+': private path, patient identifier or credential pattern')
    return {'status':'passed' if not issues else 'failed','files_checked':count,'issues':issues,
            'limitation':'Pattern scan supplements manual review; it cannot guarantee de-identification or absence of secrets'}
if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__); ap.add_argument('--require-license',action='store_true')
    args=ap.parse_args(); result=audit()
    if args.require_license and not (ROOT/'LICENSE').is_file():
        result['issues'].append('Owner/lab license approval and LICENSE file are required before publication'); result['status']='failed'
    print(json.dumps(result,indent=2)); raise SystemExit(0 if result['status']=='passed' else 1)
