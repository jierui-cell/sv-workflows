#!/usr/bin/env python3
"""Run the actual ClusterSV signature entry point end to end on synthetic VCFs."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import yaml
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'sv_signatures'))
from verify_pcawg import verify
from workflow import read_yaml

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--clustersv-source',type=Path,default=ROOT/'resources/clustersv/source')
    parser.add_argument('--rscript',default='Rscript')
    args = parser.parse_args()
    out = args.out.resolve(); out.mkdir(parents=True,exist_ok=False)
    source = args.clustersv_source.resolve()
    manifest = out/'manifest.tsv'
    with manifest.open('w') as handle:
        writer = csv.writer(handle,delimiter='\t')
        writer.writerow(('cohort','pair_id','caller','vcf','source_status','source_sha256'))
        for index in range(3):
            vcf = out/f'S{index}.vcf'
            records = [f'1\t20000000\ta\tN\t<DEL>\t100\tPASS\tSVTYPE=DEL;END=22000000',
                       f'1\t20000100\tb\tN\t<DEL>\t100\tPASS\tSVTYPE=DEL;END=22000100',
                       f'1\t70000000\tc\tN\t<DUP>\t100\tPASS\tSVTYPE=DUP;END={70400000+index*10000}',
                       f'1\t90000000\td\tN\t<INV>\t100\tPASS\tSVTYPE=INV;END=95000000;CT=3to3']
            vcf.write_text('##fileformat=VCFv4.2\n##contig=<ID=1,length=249250621>\n'
                           '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'+'\n'.join(records)+'\n')
            writer.writerow(('synthetic',f'S{index}','manta',vcf,'synthetic_control',digest(vcf)))
    # Synthetic, fully measured chr1 timing fixture; never a biological reference.
    timing = out/'timing'; timing.mkdir()
    size = 249250621; bins = (size+999)//1000
    arrays = timing/'replication_timing.npz'
    np.savez_compressed(arrays,**{'1_label':np.full(bins,3,dtype=np.uint8),
                                  '1_score':np.full(bins,.8,dtype=np.float64)})
    reference = timing/'reference.json'
    reference.write_text(json.dumps(dict(genome_build='GRCh37',label_policy='tertiles',bin_bp=1000,
        chromosome_sizes={'1':size},thresholds=[.333,.667],synthetic_control=True,
        outputs_sha256={'replication_timing.npz':digest(arrays)})))
    (timing/'verification.json').write_text(json.dumps(dict(status='verified',reference_sha256=digest(reference))))
    rwrapper = out/'Rscript.sh'
    executable = shutil.which(args.rscript) or args.rscript
    rwrapper.write_text('#!/bin/bash\nset -euo pipefail\nexec '+json.dumps(executable)+' --vanilla "$@"\n')
    rwrapper.chmod(0o755)
    shutil.copyfile(ROOT/'sv_signatures/pcawg/renv.lock',out/'renv.lock')
    cfg = read_yaml(ROOT/'config/signatures.yaml')
    cfg.update(manifest=str(manifest),method='pcawg_clustersv',pcawg_stage='signatures',
               python=sys.executable,pcawg_rscript=str(rwrapper),clustersv_source=str(source),
               pcawg_timing_reference=str(timing),outdir=str(out/'results'),group_cpus=1,
               group_memory_gb=8,max_parallel_groups=1,pcawg_rank=2,pcawg_nrun=3,
               pcawg_nmf_max_iter=500)
    params = out/'params.yaml'; params.write_text(yaml.safe_dump(cfg))
    command = ['nextflow','-log',str(out/'nextflow.log'),'run',str(ROOT/'sv_signatures/main.nf'),
               '-profile','local','-params-file',str(params),'-work-dir',str(out/'work')]
    subprocess.run(command,check=True,cwd=out,env=dict(os.environ,NXF_VER='25.10.7'))
    group = out/'results/synthetic/manta/group_output'
    result = verify(group,{'S0','S1','S2'})
    assert result['junctions'] == 12 and result['features'] == result['signatures'] == 'completed'
    assert json.loads((group/'verification.json').read_text())['status'] == 'verified'
    print('PASS actual signature main.nf: VCF → clustering/timing → 38 features → NMF → verification')

if __name__ == '__main__': main()
