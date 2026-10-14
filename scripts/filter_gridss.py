#!/usr/bin/env python3
"""Pinned native GRIDSS somatic filter; preserves raw joint VCF unchanged."""
import argparse
import gzip
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path

DEFAULTS={'allowable_normal_contamination':0.03,'min_normal_depth':8,
          'min_direct_read_support':1,'min_af':0.005,'min_qual':350,
          'max_homology_length':50,'max_inexact_homology_length':50}


def validate(options):
    if set(options)-set(DEFAULTS): raise ValueError('Unknown GRIDSS somatic filter option')
    resolved=dict(DEFAULTS); resolved.update(options)
    for k,v in resolved.items():
        import math
        if isinstance(v,bool) or not isinstance(v,(int,float)) or not math.isfinite(v) or v<0:
            raise ValueError('Invalid GRIDSS filter value: '+k)
        if k in ('allowable_normal_contamination','min_af') and v>1:
            raise ValueError('Fraction must be <=1: '+k)
    return resolved


def digest(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for block in iter(lambda:f.read(1048576),b''): h.update(block)
    return h.hexdigest()


def completed_exit_code(path):
    """Read a verified structured sidecar, or either historical timing format."""
    path = Path(path)
    sidecar = Path(str(path) + '.status.json')
    if sidecar.exists():
        status = json.loads(sidecar.read_text())
        if status['timing_sha256'] != digest(path):
            raise ValueError('Timing file changed since command completion')
        return int(status['exit_code'])
    text = path.read_text()
    try:
        return int(json.loads(text.splitlines()[-1])['exit_code'])
    except (ValueError, KeyError, IndexError):
        exits = re.findall(r'^\s*Exit status:\s*(\d+)\s*$', text, re.M)
        if not exits:
            raise ValueError('No completed exit status in timing file')
        return int(exits[-1])


def patch_r_compatibility(source):
    old='  isbp = str_detect(VariantAnnotation::fixed(vcf)$ALT,'
    if source.count(old)!=1: raise ValueError('Unexpected pinned GRIDSS ALT conversion source')
    return source.replace(old,
        '  stopifnot(all(elementNROWS(VariantAnnotation::fixed(vcf)$ALT) == 1L))\n'
        '  isbp = str_detect(as.character(unlist(VariantAnnotation::fixed(vcf)$ALT)),')


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--input',required=True,type=Path)
    p.add_argument('--out',required=True,type=Path)
    p.add_argument('--config',type=Path)
    p.add_argument('--image',type=Path)
    p.add_argument('--resume-qc',action='store_true')
    args=p.parse_args()
    root=Path(__file__).resolve().parents[1]
    image=(args.image or root/'containers/images/gridss-2.13.2--h96c455f_6.sif').resolve()
    source=args.input.resolve(); out=args.out.resolve()
    out.mkdir(parents=True,exist_ok=args.resume_qc)
    options=validate(json.loads(args.config.read_text()) if args.config else {})
    # Mount the read-only input parent before its writable output child when
    # invoked inside a Nextflow task (output/joint.vcf.gz -> output/somatic/).
    prefix=['singularity','exec','--cleanenv','--bind',str(root)+':'+str(root)+','+str(source.parent)+':'+str(source.parent)+':ro,'+str(out)+':'+str(out),str(image)]
    step='initialize'
    def run(command,target=None):
        with (out/'commands.jsonl').open('a') as f: f.write(json.dumps({'step':step,'argv':command})+'\n')
        if target:
            with open(target,'w') as f: subprocess.check_call(command,stdout=f,stderr=subprocess.STDOUT)
        else: subprocess.check_call(command)
    try:
        names=subprocess.check_output(['bcftools','query','-l',str(source)],universal_newlines=True).splitlines()
        if names!=['NORMAL','TUMOR']: raise ValueError('Unexpected GRIDSS sample order: '+str(names))
        if args.resume_qc:
            provenance=json.loads((out/'effective-settings.json').read_text())
            if completed_exit_code(out/'filter.time.txt')!=0 or provenance['input_sha256']!=digest(source):
                raise ValueError('Cannot resume QC without completed filter and identical input')
            if provenance['options'] != options or provenance['image_sha256'] != digest(image):
                raise ValueError('Cannot resume QC with changed filter options or image')
            output_qc(out,names,run,provenance)
            return
        native='/usr/local/share/gridss-2.13.2-6'
        step='native_sources'
        for name in ['gridss_somatic_filter','gridss.config.R','libgridss.R']:
            run(prefix+['cp',native+'/'+name,str(out/('native-'+name))])
        library=out/'libgridss.compat.R'
        library.write_text(patch_r_compatibility((out/'native-libgridss.R').read_text()))
        filter_prefix=list(prefix)
        filter_prefix[filter_prefix.index('--bind')+1]+=','+str(library)+':'+native+'/libgridss.R:ro'
        run(filter_prefix+['Rscript','-e','library(IRanges); library(stringr); x <- CharacterList("N[chr1:20[", "N."); stopifnot(all(elementNROWS(x)==1L)); stopifnot(identical(str_detect(as.character(unlist(x)), fixed("[")), c(TRUE,FALSE)))'],out/'compatibility-test.txt')
        text=(out/'native-gridss.config.R').read_text()
        text+='\n# Explicit pilot settings; remaining settings retain pinned native defaults\n'
        for k,v in options.items(): text+='gridss.'+k+' = '+str(v)+'\n'
        (out/'gridss.config.R').write_text(text)
        run(prefix+['gridss_somatic_filter','--help'],out/'native-help.txt')
        provenance={'stage':'somatic_filtered_no_pon','input':str(source),'input_sha256':digest(source),
                    'runner_sha256':digest(__file__),
                    'image_sha256':digest(image),'options':options,'normalordinal':1,'tumourordinal':2,
                    'panel_of_normals':'not supplied; PON filtering disabled',
                    'bsgenome':'not supplied; optional reference-dependent filter disabled',
                    'native_config_sha256':digest(out/'native-gridss.config.R'),
                    'effective_config_sha256':digest(out/'gridss.config.R')}
        provenance['r_compatibility_patch']={'original_sha256':digest(out/'native-libgridss.R'),'patched_sha256':digest(library),'change':'Explicit single-ALT character vector conversion for stringr; scientific thresholds unchanged'}
        (out/'effective-settings.json').write_text(json.dumps(provenance,indent=2)+'\n')
        step='somatic_filter'
        command=filter_prefix+['gridss_somatic_filter','--input',str(source),
                '--output',str(out/'somatic.high_confidence.vcf.gz'),
                '--fulloutput',str(out/'somatic.full.vcf.gz'),
                '--normalordinal','1','--tumourordinal','2','--configdir',str(out)]
        run([sys.executable,str(root/'scripts/measure_command.py'),str(out/'filter.time.txt')]+command)
        step='output_qc'
        output_qc(out,names,run,provenance)
    except Exception as exc:
        (out/'FAILED.json').write_text(json.dumps({'step':step,'error':str(exc)},indent=2)+'\n')
        raise


def output_qc(out,names,run,provenance):
    summaries={}
    for name in ['somatic.high_confidence.vcf.gz','somatic.full.vcf.gz']:
        path=out/name
        # VariantAnnotation's writeVcf(index=TRUE) emits .bgz in this image.
        # Preserve native files; expose stable .gz names via relative symlinks.
        native=path.with_suffix('.bgz')
        if not path.exists() and native.is_file():
            path.symlink_to(native.name)
        native_index=Path(str(native)+'.tbi')
        index=Path(str(path)+'.tbi')
        if not index.exists() and native_index.is_file():
            index.symlink_to(native_index.name)
        run(['bcftools','view','-Ov','-o','/dev/null',str(path)])
        if not Path(str(path)+'.tbi').exists(): run(['bcftools','index','-t',str(path)])
        run(['bcftools','index','-n',str(path)],out/(name+'.index-count.txt'))
        observed=subprocess.check_output(['bcftools','query','-l',str(path)],universal_newlines=True).splitlines()
        if observed!=names: raise ValueError('Filtered VCF sample order changed')
        filters={}; records={}
        with gzip.open(path,'rt') as f:
            for line in f:
                if line.startswith('#'): continue
                fields=line.rstrip().split('\t'); flag=fields[6]; filters[flag]=filters.get(flag,0)+1
                if fields[2] in records: raise ValueError('Duplicate filtered VCF ID')
                records[fields[2]]=dict(x.split('=',1) if '=' in x else (x,True) for x in fields[7].split(';'))
        pairs=set(); missing=[]
        for key,info in records.items():
            mate=info.get('MATEID')
            if mate:
                if mate not in records or records[mate].get('MATEID')!=key: missing.append(key)
                pairs.add(tuple(sorted((key,mate))))
        if missing: raise ValueError('Missing/nonreciprocal filtered mates: '+str(missing[:5]))
        summaries[name]={'records':len(records),'filters':filters,'paired_adjacencies':len(pairs),
                         'unpaired_records':sum('MATEID' not in info for info in records.values()),'sha256':digest(path)}
    (out/'SUCCESS.json').write_text(json.dumps({'stage':provenance['stage'],'outputs':summaries,'limitation':'Native somatic filtering completed; biological accuracy/control validation not established'},indent=2)+'\n')


if __name__=='__main__': main()
