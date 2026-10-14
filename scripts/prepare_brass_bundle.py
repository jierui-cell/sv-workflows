#!/usr/bin/env python3
"""Adapt trusted publisher resources to sequence-verified chr-prefixed GRCh37 and validate native readers."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess

ROOT=Path(__file__).resolve().parents[1]
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda:handle.read(1048576),b''): h.update(chunk)
    return h.hexdigest()

def prepare(args):
    for key in ('source','source_dict','reference','image','output'):
        setattr(args,key,getattr(args,key).resolve())
    if args.output.exists(): raise ValueError('Refusing existing output; choose a new bundle directory')
    args.output.parent.mkdir(parents=True,exist_ok=True)
    target_dict=args.reference.with_suffix('.dict')
    for path in (args.source_dict,target_dict,args.image,Path(str(args.reference)+'.fai')):
        if not path.is_file(): raise ValueError('Missing prerequisite: '+str(path))
    binds={ROOT,args.output.parent,args.source,args.source_dict.parent,args.reference.parent,args.image.parent}
    prefix=[args.engine,'exec','--cleanenv','--bind',','.join(str(p)+':'+str(p) for p in sorted(binds)),str(args.image)]
    commands=[]
    def run(argv,output=None):
        commands.append(argv)
        if output:
            with output.open('w') as handle: subprocess.run(argv,check=True,stdout=handle)
        else: subprocess.run(argv,check=True)
    try:
        run(prefix+['python3',str(ROOT/'scripts/adapt_brass_resources.py'),
            '--source',str(args.source),'--source-dict',str(args.source_dict),
            '--target-dict',str(target_dict),'--out',str(args.output)])
        run(prefix+['perl',str(ROOT/'scripts/validate_brass_cache.pl'),str(args.output),str(args.reference)],
            args.output/'native-validation.json')
        if json.loads((args.output/'native-validation.json').read_text()).get('status')!='passed':
            raise ValueError('Native BRASS resource readers did not validate')
        original=args.output/'get_rg_cns.R.original'
        native='/opt/wtsi-cgp/lib/perl5/auto/share/module/Sanger-CGP-Brass-Implement/Rscripts/get_rg_cns.R'
        run(prefix+['cp',native,str(original)])
        run(['patch','--fuzz=0','--output='+str(args.output/'get_rg_cns.R'),str(original),
             str(ROOT/'config/brass-gender-chromosome.patch')])
        run(prefix+['Rscript','-e','p<-commandArgs(TRUE)[1]; x<-readLines(p); parse(p); stopifnot(sum(grepl("if (seg_chr[i] == gender_chr)",x,fixed=TRUE))==1)',
             str(args.output/'get_rg_cns.R')])
        receipt={'status':'passed','assembly':'GRCh37','reference_sha256':sha(args.reference),
            'reference_fai_sha256':sha(str(args.reference)+'.fai'),'target_dictionary_sha256':sha(target_dict),
            'image_sha256':sha(args.image),'overlay_sha256':sha(args.output/'get_rg_cns.R'),
            'commands':commands,'limitations':['No normal panel; inherited external masks not independently biologically validated']}
        (args.output/'bundle-receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
        (args.output/'COMPATIBILITY_VALIDATED').touch()
    except Exception as exc:
        if args.output.is_dir():
            (args.output/'FAILED.json').write_text(json.dumps({'error':str(exc),'commands':commands},indent=2)+'\n')
        raise

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    for key in ('source','source-dict','reference','image','output'):
        ap.add_argument('--'+key,type=Path,required=True)
    ap.add_argument('--engine',choices=('apptainer','singularity'),default='apptainer')
    prepare(ap.parse_args())
