#!/usr/bin/env python3
"""Typed pilot caller adapter. Commands/settings are recorded once; inputs stay read-only."""
import argparse
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

SCHEMAS = {
    'delly': {'map_qual': (1, '-q'), 'minclip': (25, '-c'), 'mad_cutoff': (9, '-s')},
    'svaba': {'mate_lookup_min': (3, '-L'), 'chunk_size': (25000, '-c'),
              'max_reads': (50000, '-x'), 'max_coverage': (500, '-C')},
    'gridss': {},
}
FILTER_DEFAULTS = {'quality': 300, 'altaf': 0.03, 'coverage': 10, 'controlcont': 0}


def prepare_reference_cache(cache, ref):
    """Reuse only links to exactly the requested reference; never overwrite data."""
    cache = Path(cache)
    cache.mkdir(exist_ok=True)
    expected = {'ref.fa'+ext: Path(str(ref)+ext) for ext in
                ['', '.fai', '.amb', '.ann', '.bwt', '.pac', '.sa']}
    expected['ref.dict'] = Path(ref).with_suffix('.dict')
    for name, target in expected.items():
        if not target.is_file():
            raise ValueError('Missing reference cache source: ' + str(target))
        link = cache/name
        if link.is_symlink():
            if link.resolve() != target.resolve():
                raise ValueError('Conflicting reference cache link: ' + str(link))
        elif link.exists():
            raise ValueError('Refusing to overwrite reference cache file: ' + str(link))
        else:
            link.symlink_to(target.resolve())
    return cache.resolve()


def validate(tool, cfg):
    allowed = {'cpus', 'memory_gb', 'time', 'image', 'options'}
    if tool == 'delly': allowed.add('filter_options')
    if tool == 'gridss': allowed.update(['jvmheap_gb', 'otherjvmheap_gb', 'somatic_filter', 'somatic_filter_options',
                                        'picard_validation_stringency'])
    unknown = set(cfg) - allowed
    if unknown: raise ValueError('Unknown caller configuration keys: ' + str(unknown))
    opts = {k: v[0] for k, v in SCHEMAS[tool].items()}
    supplied = cfg.get('options', {})
    if set(supplied) - set(opts): raise ValueError('Unknown scientific option for ' + tool)
    opts.update(supplied)
    for k, v in opts.items():
        if isinstance(v, bool) or not isinstance(v, int) or v < 0:
            raise ValueError('Expected nonnegative integer option: ' + k)
    if int(cfg['cpus']) < 1: raise ValueError('CPUs must be positive')
    filters = dict(FILTER_DEFAULTS)
    if tool == 'delly':
        if set(cfg.get('filter_options', {})) - set(filters): raise ValueError('Unknown DELLY filter option')
        filters.update(cfg.get('filter_options', {}))
        for k, v in filters.items():
            if isinstance(v, bool) or not isinstance(v, (int, float)) or v < 0:
                raise ValueError('Invalid DELLY filter value: ' + k)
        if filters['altaf'] > 1 or filters['controlcont'] > 1: raise ValueError('Allele fractions must be <=1')
    elif cfg.get('filter_options'): raise ValueError('Filtering not implemented for this caller')
    if tool == 'gridss' and int(cfg['jvmheap_gb']) + int(cfg['otherjvmheap_gb']) >= int(cfg['memory_gb']):
        raise ValueError('GRIDSS Java heaps require memory headroom')
    if tool == 'gridss':
        from filter_gridss import validate as validate_filter
        validate_filter(cfg.get('somatic_filter_options', {}))
        if type(cfg.get('somatic_filter', True)) is not bool:
            raise ValueError('GRIDSS somatic_filter must be boolean')
        if cfg.get('picard_validation_stringency', 'STRICT') not in ('STRICT', 'LENIENT'):
            raise ValueError('GRIDSS picard_validation_stringency must be STRICT or LENIENT')
    return opts, filters


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('spec')
    parser.add_argument('--stub', action='store_true')
    args = parser.parse_args()
    spec = json.loads(Path(args.spec).read_text())
    if spec.get('project_dir'): sys.path.insert(0,str(Path(spec['project_dir'])/'scripts'))
    spec['runner_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    tool, cfg = spec['tool'], spec['settings']
    opts, filters = validate(tool, cfg)
    out = Path('output').resolve()
    out.mkdir(exist_ok=True)
    cfg['options'] = opts
    if tool == 'gridss': cfg['picard_validation_stringency'] = cfg.get('picard_validation_stringency', 'STRICT')
    if tool == 'delly': cfg['filter_options'] = filters
    spec['stage'] = 'somatic' if tool in ('delly', 'svaba') else 'raw_joint_not_somatic'
    spec['optional_external_resources'] = 'none; technical pilot baseline'
    exclude = None
    if spec.get('exclude_bed'):
        if tool not in ('delly','gridss'):
            raise ValueError('Exclusions are not implemented for ' + tool)
        from regions import prepare
        exclude = str(prepare(spec['exclude_bed'], spec['reference']+'.fai', out).resolve())
        spec['optional_external_resources'] = 'explicit validated exclusion BED; see regions.audit.json'
    (out / 'effective-settings.json').write_text(json.dumps(spec, indent=2) + '\n')
    if args.stub:
        (out / 'STUB_ONLY.txt').write_text('Validated settings only; no variants called.\n')
        return
    image = str(Path(spec['image']).resolve())
    ref = str(Path(spec['reference']).resolve())
    tumor, normal = [str(Path(spec[k]).resolve()) for k in ('tumor', 'normal')]
    mounts = {str(Path.cwd()), spec['project_dir']}
    binds = [p + ':' + p for p in sorted(mounts)]
    binds += [p + ':' + p + ':ro' for p in sorted({str(Path(x).parent) for x in [tumor, normal, ref, image]})]
    prefix = ['singularity', 'exec', '--cleanenv', '--bind', ','.join(binds), image]
    step = 'initialize'

    def run(command, target=None, help_only=False):
        if not help_only and command[:len(prefix)] == prefix:
            command = [sys.executable, str(Path(spec['project_dir'])/'scripts/measure_command.py'),
                       str(out/(step+'.time.txt'))] + command
        with (out / 'commands.jsonl').open('a') as handle:
            handle.write(json.dumps({'step': step, 'argv': command}) + '\n')
        if target:
            with Path(target).open('w') as handle:
                result = subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT)
        else: result = subprocess.run(command)
        if not help_only: result.check_returncode()

    try:
        step = 'artifact_checksums'
        with (out / 'artifacts.sha256').open('w') as handle:
            for name in [image, ref, ref + '.fai']:
                digest = hashlib.sha256()
                with open(name, 'rb') as data:
                    for chunk in iter(lambda: data.read(1024*1024), b''): digest.update(chunk)
                handle.write(digest.hexdigest() + '  ' + name + '\n')
        flags = [v for k in opts for v in (SCHEMAS[tool][k][1], str(opts[k]))]
        step = 'capture_cli_defaults'
        help_cmd = {'delly': ['delly', 'sr', '--help'], 'svaba': ['svaba','run','--help'], 'gridss':['gridss','--help']}[tool]
        run(prefix + help_cmd, out/'native-help.txt', help_only=True)
        if tool == 'delly':
            run(prefix + ['delly','filter','--help'], out/'native-filter-help.txt', help_only=True)
            step = 'sample_roles'
            names = []
            for bam in [tumor,normal]:
                header = subprocess.check_output(['samtools','view','-H',bam]).decode()
                samples = {v[3:] for line in header.splitlines() if line.startswith('@RG\t') for v in line.split('\t') if v.startswith('SM:')}
                if len(samples) != 1: raise ValueError('Expected exactly one BAM sample name')
                names.append(next(iter(samples)))
            if names[0] == names[1]: raise ValueError('Identical tumor/normal names')
            roles = out/'samples.tsv'
            roles.write_text(names[0]+'\ttumor\n'+names[1]+'\tcontrol\n')
            step = 'discovery'
            run(prefix+['delly','sr','--threads',str(cfg['cpus']),'-g',ref,'-o',str(out/'discovery.bcf')]+flags+(['-x',exclude] if exclude else [])+[tumor,normal])
            observed = subprocess.check_output(['bcftools','query','-l',str(out/'discovery.bcf')]).decode().splitlines()
            if sorted(observed) != sorted(names): raise ValueError('VCF sample mapping mismatch')
            step = 'somatic_filter'
            filter_flags = [v for k in filters for v in ('--'+k,str(filters[k]))]
            run(prefix+['delly','filter','-f','somatic','-s',str(roles),'-o',str(out/'somatic.bcf')]+filter_flags+[str(out/'discovery.bcf')])
            run(['bcftools','view','-Oz','-o',str(out/'somaticSV.vcf.gz'),str(out/'somatic.bcf')])
            run(['bcftools','index',str(out/'somaticSV.vcf.gz')])
            vcfs = [out/'somaticSV.vcf.gz']
        elif tool == 'svaba':
            for ext in ['amb','ann','bwt','pac','sa']:
                if not Path(ref+'.'+ext).is_file(): raise ValueError('Missing BWA index '+ext)
            step = 'call'
            run(prefix+['svaba','run','-t',tumor,'-n',normal,'-G',ref,'-a',str(out/spec['id']),'-p',str(cfg['cpus'])]+flags)
            vcfs = list(out.glob('*.vcf'))
            if not list(out.glob('*.svaba.somatic.sv.vcf')): raise ValueError('Missing somatic SV output')
        else:
            step = 'reference_cache'
            cache = prepare_reference_cache('reference_cache', ref)
            step = 'raw_joint_call'
            picard_options = (['--picardoptions', 'VALIDATION_STRINGENCY=LENIENT']
                              if cfg['picard_validation_stringency'] == 'LENIENT' else [])
            run(prefix+['gridss','--reference',str(cache/'ref.fa'),'--output',str(out/'joint.vcf.gz'),
                '--assembly',str(Path('assembly.bam').resolve()),'--workingdir',str(Path('gridss_work').resolve()),
                '--threads',str(cfg['cpus']),'--jvmheap',str(cfg['jvmheap_gb'])+'g',
                '--otherjvmheap',str(cfg['otherjvmheap_gb'])+'g','--labels','NORMAL,TUMOR']
                +picard_options+(['--blacklist',exclude] if exclude else [])+[normal,tumor])
            observed = subprocess.check_output(['bcftools','query','-l',str(out/'joint.vcf.gz')]).decode().splitlines()
            if observed != ['NORMAL','TUMOR']: raise ValueError('Unexpected GRIDSS sample order')
            vcfs = [out/'joint.vcf.gz']
            if cfg.get('somatic_filter', True):
                step = 'somatic_filter'
                options_path=out/'somatic-filter-options.json'
                options_path.write_text(json.dumps(cfg.get('somatic_filter_options',{}))+'\n')
                run([sys.executable,str(Path(spec['project_dir'])/'scripts/filter_gridss.py'),
                     '--input',str(out/'joint.vcf.gz'),'--out',str(out/'somatic'),
                     '--config',str(options_path),'--image',image])
                spec['stage']='somatic_filtered_no_pon'
                spec['optional_external_resources']='PON and BSgenome omitted; reference-dependent filtering/equivalence checks disabled'+('; explicit exclusions recorded in regions.audit.json' if exclude else '')
                (out/'effective-settings.json').write_text(json.dumps(spec,indent=2)+'\n')
                vcfs += [out/'somatic/somatic.high_confidence.vcf.gz',out/'somatic/somatic.full.vcf.gz']
        step = 'output_qc'
        for vcf in vcfs:
            run(['bcftools','view','-Ov','-o','/dev/null',str(vcf)])
            run(['bcftools','stats',str(vcf)], out/(vcf.name+'.stats.txt'))
        (out/'SUCCESS.json').write_text(json.dumps({'stage':spec['stage'],'vcfs':[str(p.relative_to(out)) for p in vcfs]},indent=2)+'\n')
    except Exception as exc:
        (out/'FAILED.json').write_text(json.dumps({'step':step,'error':str(exc)},indent=2)+'\n')
        raise


if __name__ == '__main__': main()
