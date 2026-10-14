#!/usr/bin/env python3
"""Fast, read-only resource checks; no caller execution or large-file hashing."""
from copy import deepcopy
import csv
import hashlib
import io
import json
import os
from pathlib import Path
import shutil
import tempfile

SAMPLE_FIELDS = ('patient_id','tumor_bam','tumor_bai','normal_bam','normal_bai',
                 'sex','tumor_phenotype','purity','ploidy')
MANIFEST_FIELDS = ('cohort','pair_id','caller','vcf','source_status','source_sha256')
BWA_SUFFIXES = ('.amb','.ann','.bwt','.pac','.sa')
BRASS_FILES = ('native-validation.json','bundle-receipt.json','get_rg_cns.R',
               'HiDepth.bed.gz','HiDepth.bed.gz.tbi','500bp_windows.gc.bed.gz',
               '500bp_windows.gc.bed.gz.tbi','vagrent.cache.gz','vagrent.cache.gz.tbi',
               'cytoband.txt','CentTelo.tsv','viral.genomic.fa.2bit') + tuple(
                   f'all_ncbi_bacteria.{i}.fa.2bit' for i in range(1,5))

def absolute(value, root):
    """Keep symlinks intact, especially a virtualenv's Python executable."""
    if not isinstance(value,str) or not value.strip():
        raise ValueError('Expected a nonempty file/directory path')
    path = Path(value).expanduser()
    return Path(os.path.abspath(path if path.is_absolute() else Path(root)/path))

def normalize_paths(config, mode, root):
    config = deepcopy(config)
    keys = {'callers':('input','reference','outdir','manta_image','exclude_bed','target_bed'),
            'signatures':('manifest','outdir','python','clustersv_source','pcawg_rscript','pcawg_timing_reference'),
            'timing':('timing_outdir','timing_python','timing_input_reference'),
            'test':('outdir',)}[mode]
    for key in keys:
        if config.get(key) not in (None,'') and not (key == 'python' and config[key] == 'python3'):
            config[key] = str(absolute(config[key],root))
    if mode == 'callers':
        for caller in ('delly','svaba','gridss','brass'):
            for key in ('image','bundle'):
                if config.get(caller,{}).get(key):
                    config[caller][key] = str(absolute(config[caller][key],root))
    return config

def read_table(path, required):
    with Path(path).open(newline='',encoding='utf-8-sig') as handle:
        rows = list(csv.reader(handle,delimiter='\t',quoting=csv.QUOTE_NONE))
    if not rows:
        raise ValueError('Table is empty: '+str(path))
    header = rows[0]
    if any(not name or name != name.strip() for name in header) or len(set(header)) != len(header):
        raise ValueError('Empty, whitespace-padded or duplicate column names: '+str(path))
    missing = set(required)-set(header)
    if missing:
        raise ValueError('Missing columns in '+str(path)+': '+', '.join(sorted(missing)))
    if len(rows) == 1:
        raise ValueError('Table has no data rows: '+str(path))
    data = []
    for number,cells in enumerate(rows[1:],2):
        if len(cells) != len(header):
            raise ValueError(f'{path} line {number}: expected {len(header)} columns, found {len(cells)}')
        data.append(dict(zip(header,cells)))
    return header,data

def selected_samples(config,root):
    header,rows = read_table(absolute(config['input'],root),SAMPLE_FIELDS)
    seen = set()
    for row in rows:
        for key in SAMPLE_FIELDS[:5]:
            if not row[key].strip():
                raise ValueError('Samplesheet has empty required value '+key)
        if row['patient_id'] in seen:
            raise ValueError('Duplicate patient_id: '+row['patient_id'])
        seen.add(row['patient_id'])
    if not isinstance(config.get('pair_ids',''),str):
        raise ValueError('pair_ids must be a comma-separated string')
    ids = config.get('pair_ids','').split(',') if config.get('pair_ids') else []
    if len(set(ids)) != len(ids) or set(ids)-seen:
        raise ValueError('Duplicate or unknown pair_ids selection')
    selected = [row for row in rows if not ids or row['patient_id'] in ids]
    for row in selected:
        for key in SAMPLE_FIELDS[1:5]:
            row[key] = str(absolute(row[key],root))
    return header,selected

def write_resolved_samplesheet(config,root,destination):
    """Use a stable content-addressed input path so resumed BRASS tasks stay cached."""
    header,rows = selected_samples(config,root)
    handle = io.StringIO()
    writer = csv.DictWriter(handle,fieldnames=header,delimiter='\t',lineterminator='\n')
    writer.writeheader(); writer.writerows(rows)
    content = handle.getvalue()
    digest = hashlib.sha256(content.encode()).hexdigest()
    cached = Path(root)/'.nextflow/input-manifests'/(digest+'.tsv')
    cached.parent.mkdir(parents=True,exist_ok=True)
    if not cached.exists():
        with tempfile.NamedTemporaryFile(mode='w',dir=cached.parent,delete=False) as temp:
            temp.write(content)
            temporary = Path(temp.name)
        try:
            try:
                os.link(temporary,cached)  # Atomic create; never replace an existing cache file.
            except FileExistsError:
                pass
        finally:
            temporary.unlink()
    if cached.read_text() != content:
        raise ValueError('Conflicting resolved input cache: '+str(cached))
    with Path(destination).open('x') as handle:
        handle.write(content)
    return cached

def check_files(config,mode,root):
    if mode not in ('callers','signatures'):
        raise ValueError('--check-files supports callers and signatures')
    issues,checked = [],set()
    count = 0
    def require(value,label,allow_empty=False):
        path = absolute(str(value),root)
        if str(path) in checked:
            return path.is_file()
        checked.add(str(path))
        try:
            if not path.is_file():
                raise ValueError('missing file (or broken symlink)')
            if not os.access(path,os.R_OK):
                raise ValueError('file is not readable')
            if not allow_empty and path.stat().st_size == 0:
                raise ValueError('file is empty')
        except (OSError,ValueError) as exc:
            issues.append(f'{label}: {path} — {exc}')
            return False
        return True
    if mode == 'callers':
        reference = absolute(config['reference'],root)
        require(str(reference),'reference')
        require(str(reference)+'.fai','reference FASTA index')
        callers = config['callers'].split(',')
        if set(callers)&{'svaba','gridss'}:
            for suffix in BWA_SUFFIXES:
                require(str(reference)+suffix,'BWA index '+suffix)
        if set(callers)&{'gridss','brass'}:
            require(str(reference.with_suffix('.dict')),'reference sequence dictionary')
        native_reference = reference.resolve()
        if native_reference != reference and set(callers)&{'delly','svaba','gridss'}:
            # Those adapters canonicalize FASTA symlinks before native execution.
            require(str(native_reference)+'.fai','linked reference source FASTA index')
            if set(callers)&{'svaba','gridss'}:
                for suffix in BWA_SUFFIXES:
                    require(str(native_reference)+suffix,'linked reference source BWA index '+suffix)
            if 'gridss' in callers:
                require(str(native_reference.with_suffix('.dict')),'linked reference source dictionary')
        for caller in callers:
            if caller != 'preflight':
                require(config['manta_image'] if caller == 'manta' else config[caller]['image'],caller+' image')
        for key in ('exclude_bed','target_bed'):
            if config.get(key):
                require(config[key],key)
        if require(config['input'],'samplesheet'):
            try:
                _,rows = selected_samples(config,root)
                count = len(rows)
                for row in rows:
                    for key in SAMPLE_FIELDS[1:5]:
                        require(row[key],row['patient_id']+' '+key)
            except (OSError,ValueError) as exc:
                issues.append(str(exc))
        if 'brass' in callers:
            bundle = absolute(config['brass']['bundle'],root)
            require(str(bundle/'COMPATIBILITY_VALIDATED'),'BRASS validation marker',allow_empty=True)
            for name in BRASS_FILES:
                require(str(bundle/name),'BRASS bundle '+name)
            for name in ('native-validation.json','bundle-receipt.json'):
                path = bundle/name
                if path.is_file() and os.access(path,os.R_OK):
                    try:
                        record = json.loads(path.read_text())
                        if not isinstance(record,dict) or record.get('status') != 'passed':
                            raise ValueError('validation status is not passed')
                        if name == 'bundle-receipt.json' and record.get('assembly') != config['genome_build']:
                            raise ValueError('bundle assembly differs from genome_build')
                    except (OSError,ValueError) as exc:
                        issues.append('BRASS '+str(path)+' — '+str(exc))
    else:
        python = config['python']
        require((shutil.which(python) or python) if python == 'python3' else python,'signature Python environment')
        if require(config['manifest'],'signature manifest'):
            try:
                header,rows = read_table(absolute(config['manifest'],root),MANIFEST_FIELDS)
                if tuple(header) != MANIFEST_FIELDS:
                    raise ValueError('Signature manifest must use the documented six-column order')
                groups = config.get('groups','all')
                if not isinstance(groups,str) or not isinstance(config.get('pair_ids',''),str):
                    raise ValueError('groups and pair_ids must be comma-separated strings')
                ids = [config['pair_id']] if config.get('pair_id') else config.get('pair_ids','').split(',')
                rows = [row for row in rows if (groups == 'all' or row['cohort']+':'+row['caller'] in groups.split(','))
                        and (ids == [''] or row['pair_id'] in ids)]
                if not rows:
                    raise ValueError('Signature selection matches no manifest rows')
                count = len(rows)
                for row in rows:
                    if not Path(row['vcf']).is_absolute():
                        raise ValueError('Signature VCF paths must be absolute; build the manifest with scripts/build_sv_manifest.py')
                    require(row['vcf'],row['cohort']+':'+row['caller']+':'+row['pair_id']+' VCF')
            except (OSError,ValueError) as exc:
                issues.append(str(exc))
    return dict(status='passed' if not issues else 'failed',workflow=mode,
                input_source=str(absolute(config['input'] if mode == 'callers' else config['manifest'],root)),
                files_checked=len(checked),selected_rows=count,issues=issues,
                scope='file presence/readability only; BAM/reference and native tool checks run on compute nodes')
