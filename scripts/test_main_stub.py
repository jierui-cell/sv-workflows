#!/usr/bin/env python3
"""Exercise main.nf with real synthetic BAM preflight, caller stubs and resume.

No caller inference is performed. CI tests four caller branches without containers;
--manta-image additionally tests configuration inside the real pinned Manta image.
"""
import argparse
import csv
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import yaml

ROOT = Path(__file__).resolve().parents[1]

def run(command, **kwargs):
    print('+ ' + ' '.join(map(str,command)), flush=True)
    return subprocess.run(list(map(str,command)),check=True,**kwargs)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',required=True,type=Path)
    parser.add_argument('--manta-image',type=Path)
    parser.add_argument('--container-engine',choices=('apptainer','singularity'),help='Runtime for the optional Manta container check')
    parser.add_argument('--debug-cache',action='store_true',help='Include detailed hash/cache diagnostics in test logs')
    parser.add_argument('--with-exclusions',action='store_true',help='Exercise the optional DELLY/GRIDSS/Manta region policy')
    args = parser.parse_args()
    base = args.out.resolve()
    base.mkdir(parents=True,exist_ok=False)
    data = base/'data'; data.mkdir()
    aliases = base/'aliases'; aliases.mkdir()
    ref = data/'reference.fa'
    ref.write_text('>1\n'+('ACGT'*2500)+'\n')
    run(['samtools','faidx',ref])
    for suffix in ('.amb','.ann','.bwt','.pac','.sa'):
        Path(str(ref)+suffix).write_text('STUB ONLY: no BWA index\n')
    ref.with_suffix('.dict').write_text('@HD\tVN:1.6\n@SQ\tSN:1\tLN:10000\n')
    files = {}
    for role in ('tumor','normal'):
        sam = data/(role+'.sam')
        sam.write_text('@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:1\tLN:10000\n'
                       '@RG\tID:rg\tSM:'+role.upper()+'\n'+
                       role+'_read\t0\t1\t100\t60\t20M\t*\t0\t0\tACGTACGTACGTACGTACGT\t'+
                       'I'*20+'\tRG:Z:rg\n')
        bam = data/(role+'.bam')
        run(['samtools','view','-b','-o',bam,sam])
        run(['samtools','index',bam])
        for path in (bam,Path(str(bam)+'.bai')):
            alias = aliases/path.name
            alias.symlink_to(path)
        files[role] = aliases/bam.name
    sheet = base/'samples.tsv'
    with sheet.open('w',newline='') as handle:
        writer = csv.writer(handle,delimiter='\t')
        writer.writerow(['patient_id','tumor_bam','tumor_bai','normal_bam','normal_bai',
                         'sex','tumor_phenotype','purity','ploidy'])
        writer.writerow(['SYNTHETIC',files['tumor'],str(files['tumor'])+'.bai',
                         files['normal'],str(files['normal'])+'.bai','unknown','synthetic','',''])
    cfg = yaml.safe_load((ROOT/'config/defaults.yaml').read_text())
    engine = args.container_engine or os.environ.get('SV_CONTAINER_ENGINE') or ('apptainer' if shutil.which('apptainer') else 'singularity')
    cfg['container_engine'] = engine
    callers = ['delly','gridss'] if args.with_exclusions else ['delly','svaba','gridss','brass']
    image = data/'STUB_ONLY.sif'; image.write_text('Not a container; caller stubs only\n')
    cfg.update(input=str(sheet),reference=str(ref),dataset_id='synthetic_stub',
               genome_build='synthetic',allow_local_execution=True,
               outdir=str(base/'results'),max_parallel_callers=1,max_parallel_preflight=1)
    for caller in callers:
        cfg[caller].update(image=str(image),cpus=1,memory_gb=4)
    cfg['gridss'].update(jvmheap_gb=1,otherjvmheap_gb=1)
    if args.with_exclusions:
        mask = data/'exclude.bed'; mask.write_text('1\t500\t600\n')
        cfg['exclude_bed'] = str(mask)
    if args.manta_image:
        callers.append('manta')
        cfg.update(manta_image=str(args.manta_image.resolve(strict=True)),manta_cpus=1,
                   manta_memory_gb=4,manta_workflow_memory_gb=2)
    cfg['callers'] = ','.join(callers)
    params = base/'params.yaml'; params.write_text(yaml.safe_dump(cfg))
    env = dict(os.environ,NXF_VER='25.10.7')
    env['SV_CONTAINER_ENGINE'] = engine
    env['PATH'] = str(ROOT/'bin')+os.pathsep+env.get('PATH','')
    # Exercise scheduler recognition and explicit local acknowledgement, without
    # pretending these fixtures constitute execution on other schedulers.
    clean = {k:v for k,v in env.items() if k not in ('SLURM_JOB_ID','PBS_JOBID','LSB_JOBID','JOB_ID')}
    preflight = ['bash',ROOT/'scripts/validate_pair.sh','SYNTHETIC',str(ref)+'.fai',
                 files['tumor'],files['normal']]
    for key in (None,'PBS_JOBID','LSB_JOBID','JOB_ID'):
        context = base/('preflight-'+(key or 'local')); context.mkdir()
        context_env = dict(clean)
        if key: context_env[key] = 'synthetic_context_test'
        run(preflight+['1' if key is None else '0'],cwd=context,env=context_env)
    rejected = base/'preflight-unacknowledged'; rejected.mkdir()
    result = subprocess.run(list(map(str,preflight+['0'])),cwd=rejected,env=clean)
    if result.returncode == 0:
        raise RuntimeError('Unacknowledged local preflight must be rejected')
    command = ['run',ROOT/'main.nf','-stub-run','-profile','local',
               '-params-file',params,'-work-dir',base/'work']
    diagnostic = ['nextflow']
    if args.debug_cache:
        diagnostic += ['-trace','nextflow.processor.TaskProcessor,nextflow.cache.CacheDB']
        command += ['-dump-hashes','json']
    run(diagnostic+['-log',base/'first.nextflow.log']+command+['-with-trace',base/'first.trace.tsv'],cwd=base,env=env)
    run(diagnostic+['-log',base/'resume.nextflow.log']+command+['-resume','-with-trace',base/'resume.trace.tsv'],cwd=base,env=env)
    with (base/'resume.trace.tsv').open() as handle:
        tasks = list(csv.DictReader(handle,delimiter='\t'))
    if len(tasks) != len(callers)+1 or any(row['status']!='CACHED' for row in tasks):
        raise RuntimeError('Expected preflight and every requested caller stub to resume from cache: '+repr(tasks))
    for role in ('tumor','normal'):
        count = base/'results/preflight/SYNTHETIC'/(role+'.interval_read_count.txt')
        if int(count.read_text().strip()) != 1:
            raise RuntimeError('Indexed query on numbered chromosome did not retrieve the synthetic read')
    summary = {'status':'passed','callers':callers,'cached_tasks':len(tasks),
               'scope':'real synthetic BAM preflight, task wiring and resume; no caller inference',
               'manta_container_configuration_tested':bool(args.manta_image),
               'exclusions_tested':args.with_exclusions}
    if args.with_exclusions:
        for caller in ('delly','gridss'):
            audit = json.loads((base/'results'/caller/'SYNTHETIC/output/regions.audit.json').read_text())
            if audit['excluded_bases'] != 100:
                raise RuntimeError('Incorrect excluded-base count')
        if args.manta_image:
            matches = list((base/'results/manta/SYNTHETIC').rglob('regions.audit.json'))
            if len(matches)!=1 or json.loads(matches[0].read_text())['included_bases'] != 9900:
                raise RuntimeError('Incorrect Manta inclusion-list complement')
    (base/'validation.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))

if __name__ == '__main__': main()
