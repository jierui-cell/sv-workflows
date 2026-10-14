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

def validate_adapted(source, target_dict):
    """Revalidate an existing adapted resource identity, never just trust a marker."""
    report=json.loads((source/'adaptation.json').read_text())
    if report.get('status')!='adapted_requires_native_reader_validation' or report.get('coordinate_changes')!=0:
        raise ValueError('Unexpected prior resource-adaptation state')
    if report.get('target_dictionary_sha256')!=sha(target_dict):
        raise ValueError('Existing BRASS resources were adapted to another reference dictionary')
    if not report.get('resources'):
        raise ValueError('Existing resource checksum inventory is absent')
    for name, info in report['resources'].items():
        expected=info.get('adapted_sha256') or (info.get('sha256') if info.get('action')=='unchanged_symlink' else None)
        if Path(name).name!=name or expected!=sha(source/name):
            raise ValueError('Existing adapted resource changed: '+name)
    return report

def prepare(args):
    for key in ('source','reference','image','output'):
        setattr(args,key,getattr(args,key).resolve())
    if args.source_dict: args.source_dict=args.source_dict.resolve()
    if args.output.exists(): raise ValueError('Refusing existing output; choose a new bundle directory')
    args.output.parent.mkdir(parents=True,exist_ok=True)
    target_dict=args.reference.with_suffix('.dict')
    required=[target_dict,args.image,Path(str(args.reference)+'.fai')]
    if not args.reuse_adapted:
        if not args.source_dict: raise ValueError('--source-dict is required for publisher-resource adaptation')
        required.append(args.source_dict)
    for path in required:
        if not path.is_file(): raise ValueError('Missing prerequisite: '+str(path))
    binds={ROOT,args.output.parent,args.source,args.reference.parent,args.image.parent}
    if args.reuse_adapted:
        binds.update(p.resolve().parent for p in args.source.iterdir() if p.is_file())
    if args.source_dict: binds.add(args.source_dict.parent)
    prefix=[args.engine,'exec','--cleanenv','--bind',','.join(str(p)+':'+str(p) for p in sorted(binds)),str(args.image)]
    commands=[]
    def run(argv,output=None):
        commands.append(argv)
        if output:
            with output.open('w') as handle: subprocess.run(argv,check=True,stdout=handle)
        else: subprocess.run(argv,check=True)
    try:
        if args.reuse_adapted:
            validate_adapted(args.source,target_dict)
            args.output.mkdir()
            skipped={'COMPATIBILITY_VALIDATED','native-validation.json','bundle-receipt.json',
                     'get_rg_cns.R','get_rg_cns.R.original','FAILED.json'}
            for path in args.source.iterdir():
                if path.is_file() and path.name not in skipped:
                    (args.output/path.name).symlink_to(path.resolve())
        else:
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
            'preparation':'revalidated_existing_adaptation' if args.reuse_adapted else 'publisher_adaptation',
            'source_adaptation_sha256':sha(args.source/'adaptation.json') if args.reuse_adapted else None,
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
    for key in ('source','reference','image','output'):
        ap.add_argument('--'+key,type=Path,required=True)
    ap.add_argument('--source-dict',type=Path)
    ap.add_argument('--reuse-adapted',action='store_true',help='Verify existing adapted resource hashes/dictionary, then rerun native readers and regenerate the overlay/receipt in a new directory')
    ap.add_argument('--engine',choices=('apptainer','singularity'),default='apptainer')
    prepare(ap.parse_args())
