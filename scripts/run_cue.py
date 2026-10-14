"""Run pinned Cue independently per sample; never claim somatic classification."""
import hashlib,json,os,shutil,subprocess,sys
from pathlib import Path
from cue_read_policy import POLICY,patch_io

def interpreter_path(root, configured):
    # Preserve the virtualenv symlink: resolving it bypasses pyvenv.cfg.
    return os.path.abspath(str(Path(root)/configured))

def configs(spec):
    c=spec['settings'];d=c['data_options'];m=c['model_options']
    if set(d)-{'chr_names','logging_level','min_sv_len','min_qual_score','refine_disable'}:raise ValueError('Unsupported data options')
    if set(m)-{'gpu_ids','n_jobs_per_gpu','n_cpus','batch_size','report_interval'}:raise ValueError('Unsupported model options')
    if m['gpu_ids']:raise ValueError('This initial workflow is CPU-only; GPU execution needs an explicit Slurm profile')
    if m['n_cpus']!=c['cpus'] or m['n_cpus']<1:raise ValueError('CPU configuration mismatch')
    if spec['role'] not in ['normal','tumor']:raise ValueError('Invalid role')
    root=Path(spec['root'])
    return dict(d,bam=spec['bam'],fai=spec['fai']),dict(m,model_path=str((root/c['model']).resolve()))

def main():
    spec=json.load(open(sys.argv[1]));data,model=configs(spec);root=Path(spec['root']);c=spec['settings']
    subprocess.check_call(['sha256sum','-c',str(root/'resources/cue/READY.sha256')],cwd=str(root))
    source=(root/c['source']).resolve()
    commit=subprocess.check_output(['git','-C',str(source),'rev-parse','HEAD'],universal_newlines=True).strip()
    if commit!=c['source_commit']:raise ValueError('Source commit differs')
    subprocess.check_call(['git','-C',str(source),'diff','--exit-code','HEAD'])
    out=Path('output');out.mkdir();native=Path('native');native.mkdir()
    # Task-local source adapter and BAM aliases keep upstream and input files intact.
    adapted=Path('cue-source').absolute()
    shutil.copytree(str(source),str(adapted),ignore=shutil.ignore_patterns('.git','__pycache__','*.pyc'))
    io_path=adapted/'seq/io.py'
    original=io_path.read_text();io_path.write_text(patch_io(original))
    shutil.copy2(str(root/'scripts/cue_read_policy.py'),str(adapted/'cue_read_policy.py'))
    shutil.copy2(str(adapted/'cue_read_policy.py'),str(out/'cue_read_policy.py'))
    shutil.copy2(str(io_path),str(out/'cue-adapted-io.py'))
    (out/'compatibility.json').write_text(json.dumps(dict(policy=POLICY,
        upstream_io_sha256=hashlib.sha256(original.encode()).hexdigest(),
        adapted_io_sha256=hashlib.sha256(io_path.read_bytes()).hexdigest(),
        mapped_missing_sequence='fail', cache='task_local_bam_alias'),indent=2)+'\n')
    bam=Path(spec['bam']).resolve();alias=native/'input.bam'
    alias.symlink_to(bam)
    indexes=[Path(str(bam)+'.bai'),bam.with_suffix('.bai')]
    index=next((p for p in indexes if p.is_file()),None)
    if index is None:raise ValueError('Missing BAM index for Cue')
    (native/'input.bam.bai').symlink_to(index)
    data['bam']=str(alias.absolute())
    # JSON is valid YAML and avoids dependence on host PyYAML.
    for name,obj in [('data',data),('model',model)]:
        (native/(name+'.yaml')).write_text(json.dumps(obj,indent=2)+'\n')
        shutil.copy2(str(native/(name+'.yaml')),str(out/(name+'.yaml')))
    (out/'effective-settings.json').write_text(json.dumps(dict(spec,classification='single_sample_not_somatic',scope=data['chr_names']),indent=2)+'\n')
    cmd=['singularity','exec','--cleanenv','--bind',str(root.parent)+':'+str(root.parent)+','+str(bam.parent)+':'+str(bam.parent)+':ro',str((root/c['image']).resolve()),
         'env','PYTHONPATH='+str(adapted),'LD_LIBRARY_PATH='+str(root/'resources/cue/runtime-libs/root/usr/lib/x86_64-linux-gnu')+':'+str(root/'resources/cue/runtime-libs/root/lib/x86_64-linux-gnu'),'OMP_NUM_THREADS=1','MKL_NUM_THREADS=1',interpreter_path(root,c['python']),str(adapted/'engine/call.py'),
         '--data_config',str((native/'data.yaml').resolve()),'--model_config',str((native/'model.yaml').resolve())]
    (out/'command.json').write_text(json.dumps(cmd,indent=2)+'\n')
    with (out/'driver.log').open('w') as log:
        subprocess.check_call(['python3',str(root/'scripts/measure_command.py'),str(out/'time.txt')]+cmd,stdout=log,stderr=subprocess.STDOUT)
    vcf=native/'reports/svs.vcf'
    if not vcf.is_file():raise ValueError('Missing native VCF')
    shutil.copy2(str(vcf),str(out/'single_sample.svs.vcf'))
    subprocess.check_call(['bcftools','view','-Ov','-o','/dev/null',str(out/'single_sample.svs.vcf')])
    shutil.copy2(str(source/'engine/config_utils.py'),str(out/'upstream-defaults.py'))
    shutil.copy2(str(root/'resources/cue/READY.sha256'),str(out/'artifacts.sha256'))
    if (native/'logs/main.log').exists():shutil.copy2(str(native/'logs/main.log'),str(out/'native-config-and-run.log'))
    (out/'SUCCESS.json').write_text(json.dumps(dict(stage='single_sample_not_somatic',scope=data['chr_names'],source_commit=commit))+'\n')
if __name__=='__main__':main()
