#!/usr/bin/env python3
"""One launcher for portable calling, signature analysis, reference setup and smoke tests."""
import argparse
from copy import deepcopy
from datetime import datetime, timezone
import hashlib
import json
import os
import re
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import yaml
from check_resources import check_files, normalize_paths, write_resolved_samplesheet

ROOT = Path(__file__).resolve().parents[1]
DEFAULTS = {'callers': 'defaults.yaml', 'signatures': 'signatures.yaml',
            'timing': 'replication.yaml', 'test': 'test.yaml'}
ENTRY = {'callers': 'main.nf', 'signatures': 'sv_signatures/main.nf',
         'timing': 'sv_signatures/pcawg/replication_timing.nf', 'test': 'tests/smoke.nf'}

class UniqueLoader(yaml.SafeLoader):
    pass

def unique_mapping(loader, node, deep=False):
    result = {}
    for key_node, value_node in node.value:
        key = loader.construct_object(key_node, deep=deep)
        if key in result:
            raise ValueError(f'Duplicate YAML key: {key}')
        result[key] = loader.construct_object(value_node, deep=deep)
    return result

UniqueLoader.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, unique_mapping)

def read_yaml(path):
    result = yaml.load(Path(path).read_text(), Loader=UniqueLoader)
    if not isinstance(result, dict):
        raise ValueError(f'Expected a YAML mapping: {path}')
    return result

def merge(base, overlay, prefix=''):
    """Deep override, but never silently accept a mistyped parameter name."""
    result = deepcopy(base)
    for key, value in overlay.items():
        if key not in base:
            raise ValueError(f'Unknown configuration key: {prefix}{key}')
        if isinstance(base[key], dict):
            if not isinstance(value, dict):
                raise ValueError(f'Expected a mapping: {prefix}{key}')
            result[key] = merge(base[key], value, prefix + key + '.')
        else:
            result[key] = value
    return result

def override(config, text):
    key, sep, raw = text.partition('=')
    if not sep:
        raise ValueError('--set expects dotted.key=YAML_VALUE')
    parts = key.split('.')
    node = config
    for part in parts[:-1]:
        if part not in node or not isinstance(node[part], dict):
            raise ValueError(f'Unknown configuration key: {key}')
        node = node[part]
    if parts[-1] not in node:
        raise ValueError(f'Unknown configuration key: {key}')
    node[parts[-1]] = yaml.load(raw, Loader=UniqueLoader)

def validate(config, mode):
    for key in ('max_parallel_callers', 'max_parallel_preflight', 'max_parallel_groups',
                'group_cpus', 'group_memory_gb'):
        if key in config and (type(config[key]) is not int or config[key] < 1):
            raise ValueError(f'{key} must be a positive integer')
    if config.get('container_engine', 'apptainer') not in ('apptainer', 'singularity'):
        raise ValueError('container_engine must be apptainer or singularity (SIF images)')
    if mode == 'callers':
        if config['genome_build'] not in ('GRCh37','GRCh38'):
            raise ValueError('genome_build must be GRCh37 or GRCh38')
        callers = config['callers'].split(',')
        supported = {'preflight', 'manta', 'delly', 'svaba', 'gridss', 'brass'}
        if not callers or set(callers) - supported:
            raise ValueError('Unknown/discontinued caller; supported: ' + ','.join(sorted(supported)))
        if config['sequencing_type'] not in ('WGS', 'WES', 'targeted'):
            raise ValueError('sequencing_type must be WGS, WES or targeted')
        if config['sequencing_type'] != 'WGS' and set(callers) - {'preflight','manta'}:
            raise ValueError('The implemented capture-assay route uses callers=manta; other adapters currently use WGS settings')
        if config['sequencing_type'] == 'targeted' and not config['target_bed']:
            raise ValueError('targeted sequencing requires target_bed')
        if config['target_bed'] and (config['sequencing_type'] == 'WGS' or set(callers)-{'preflight','manta'}):
            raise ValueError('target_bed is implemented for WES/targeted Manta runs')
        if config.get('exclude_bed') and set(callers)-{'preflight','manta','delly','gridss'}:
            raise ValueError('exclude_bed is supported for Manta, DELLY and GRIDSS only; run other callers separately')
        if 'brass' in callers and config['genome_build'] != config['brass']['assembly']:
            raise ValueError('BRASS assembly and genome_build must agree')
        for caller in ('delly', 'svaba', 'gridss'):
            if caller in callers:
                sys.path.insert(0, str(ROOT / 'scripts'))
                from run_caller import validate as validate_caller
                validate_caller(caller, config[caller])
        if config['manta_workflow_memory_gb'] >= config['manta_memory_gb']:
            raise ValueError('Manta workflow memory must be below the allocated memory')
    if mode == 'signatures':
        if config['method'] == 'pcawg_clustersv':
            raise ValueError('PCAWG/ClusterSV signature analysis is on hold; use method=cosmic_sv32')
        if config['genome_build'] != 'GRCh37':
            raise ValueError('The supported signature routes require genome_build=GRCh37')
        if config['method'] != 'cosmic_sv32':
            raise ValueError('Unknown signature method')
        for key in ('groups','pair_id','pair_ids'):
            if not isinstance(config[key],str):
                raise ValueError(key + ' must be a string')
        source_key = 'manifest' if config.get('manifest') else 'caller_config'
        if not isinstance(config.get(source_key),str) or not config[source_key].strip():
            raise ValueError('Set caller_config to your caller YAML, or manifest for external VCFs')
        if config['pcawg_stage'] not in ('clustering', 'features', 'signatures'):
            raise ValueError('Unknown pcawg_stage')
        if config['sequencing_type'] != 'WGS':
            raise ValueError('Genome-wide signature methods require sequencing_type=WGS; capture callsets have a different ascertainment/null model')
        for key in ('pcawg_rank_min','pcawg_rank_max','pcawg_nrun','pcawg_seed','pcawg_nmf_max_iter'):
            if type(config[key]) is not int or config[key] < (0 if key == 'pcawg_seed' else 1):
                raise ValueError(key + ' must be an appropriate nonnegative/positive integer')
        if config['pcawg_rank_min'] < 2 or config['pcawg_rank_max'] < config['pcawg_rank_min']:
            raise ValueError('NMF ranks must satisfy 2 <= rank_min <= rank_max')
        if type(config['pcawg_rank']) is not int or config['pcawg_rank'] < 0 or config['pcawg_rank'] == 1:
            raise ValueError('pcawg_rank must be 0 (automatic) or an integer >=2')
        if type(config['pcawg_nmf_epsilon']) not in (int,float) or not 0 < config['pcawg_nmf_epsilon'] < 1:
            raise ValueError('pcawg_nmf_epsilon must be between 0 and 1')
        if config['pcawg_timing_source'] not in ('none','reference'):
            raise ValueError('pcawg_timing_source must be none or reference')
        if config.get('pcawg_mask_policy') not in ('local_cen_telo_exclusion','primary_only'):
            raise ValueError('Unknown ClusterSV mask policy')

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('workflow', choices=DEFAULTS)
    ap.add_argument('--config', type=Path, action='append', default=[], help='Repeatable YAML overlays, in order')
    ap.add_argument('--set', action='append', default=[], metavar='KEY=VALUE', help='Typed dotted-key override')
    ap.add_argument('--profile', choices=('local','slurm','pbs','lsf','sge','test'), default='slurm')
    ap.add_argument('--site-config', type=Path, help='Additional Nextflow scheduler configuration')
    ap.add_argument('--resume', nargs='?', const='latest')
    inspection = ap.add_mutually_exclusive_group()
    inspection.add_argument('--dry-run', action='store_true', help='Validate YAML and show command; submit nothing')
    inspection.add_argument('--check-files', action='store_true', help='Check selected input/resource files; submit nothing')
    ap.add_argument('--allow-local', action='store_true', help='Confirm a workstation/allocated node, never a login node')
    args = ap.parse_args()
    try:
        config = read_yaml(ROOT / 'config' / DEFAULTS[args.workflow])
        for path in args.config:
            config = merge(config, read_yaml(path))
        for item in args.set:
            override(config, item)
        validate(config, args.workflow)
        config = normalize_paths(config,args.workflow,ROOT)
        if args.check_files:
            report = check_files(config,args.workflow,ROOT)
            if args.site_config and not args.site_config.is_file():
                report['issues'].append('Missing --site-config: '+str(args.site_config))
                report['status'] = 'failed'
            print(json.dumps(report,indent=2))
            raise SystemExit(0 if report['status'] == 'passed' else 2)
        if args.profile == 'test' and args.workflow != 'test':
            raise ValueError('The test profile is reserved for synthetic installation tests')
        if args.profile == 'local' and args.workflow != 'test' and not args.allow_local:
            raise ValueError('Use --allow-local only on a workstation/allocated compute node')
        if args.workflow == 'callers':
            config['allow_local_execution'] = args.profile == 'local' and args.allow_local
        stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S.%fZ')
        log_dir = ROOT / 'logs' / (args.workflow + '-' + stamp)
        command = ['nextflow', '-log', str(log_dir / 'nextflow.log'), 'run',
                   str(ROOT / ENTRY[args.workflow]), '-profile', args.profile,
                   '-params-file', str(log_dir / 'resolved.yaml'),
                   '-work-dir', str(ROOT / 'work' / args.workflow),
                   '-with-trace', str(log_dir / 'trace.tsv'),
                   '-with-report', str(log_dir / 'report.html'),
                   '-with-timeline', str(log_dir / 'timeline.html')]
        if args.site_config:
            command.extend(['-c', str(args.site_config.resolve())])
        if args.resume:
            command.extend(['-resume'] + ([] if args.resume == 'latest' else [args.resume]))
        if args.dry_run:
            print(yaml.safe_dump(config, sort_keys=False))
            print(shlex.join(command))
            return
        resource_report = None
        if args.workflow in ('callers','signatures'):
            resource_report = check_files(config,args.workflow,ROOT)
            if resource_report['status'] != 'passed':
                raise ValueError('Setup files are missing or invalid:\n'+'\n'.join(resource_report['issues'])+
                                 '\nUse --check-files to inspect setup without submitting jobs.')
        if args.site_config and not args.site_config.is_file():
            raise ValueError('Missing --site-config: '+str(args.site_config))
        if not shutil.which('nextflow'):
            raise ValueError('Install Nextflow 25.10.7 and Java 17+ or load the corresponding HPC modules')
        log_dir.mkdir(parents=True)
        if resource_report:
            (log_dir/'setup-check.json').write_text(json.dumps(resource_report,indent=2)+'\n')
        if args.workflow == 'callers':
            config['input'] = str(write_resolved_samplesheet(config,ROOT,log_dir/'resolved-samplesheet.tsv'))
        if args.workflow == 'signatures' and not config.get('manifest'):
            from signature_inputs import write_auto_manifest
            config['manifest'] = str(write_auto_manifest(config,ROOT,log_dir))
        (log_dir / 'resolved.yaml').write_text(yaml.safe_dump(config, sort_keys=False))
        hashes = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                  for folder in ('scripts','modules','lib','sv_signatures','config')
                  for p in (ROOT / folder).rglob('*')
                  if p.is_file() and p.suffix in ('.py','.sh','.pl','.nf','.groovy','.yaml','.R','.json','.lock','.txt','.tsv','.patch','.config')}
        for name in ('main.nf','nextflow.config'):
            hashes[name] = hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
        env = dict(os.environ, NXF_VER='25.10.7',
                   SV_CONTAINER_ENGINE=config.get('container_engine','apptainer'))
        nf_version = subprocess.check_output(['nextflow','-version'], text=True,env=env)
        if not re.search(r'\bversion\s+25\.10\.7\b',nf_version):
            raise ValueError('Expected Nextflow 25.10.7; load/install the pinned version before running')
        (log_dir / 'provenance.json').write_text(json.dumps({'argv': command, 'code_sha256': hashes,
            'nextflow_version': nf_version,
            'overlays': [str(p.resolve()) for p in args.config], 'overrides': args.set}, indent=2) + '\n')
        env['PATH'] = str(ROOT / 'bin') + os.pathsep + env.get('PATH','')
        with (log_dir / 'driver.log').open('w') as log:
            process = subprocess.Popen(command, cwd=ROOT, env=env, stdout=subprocess.PIPE,
                                       stderr=subprocess.STDOUT, text=True)
            for line in process.stdout:
                sys.stdout.write(line); log.write(line); log.flush()
            rc = process.wait()
        (log_dir / 'status.json').write_text(json.dumps({'exit_code': rc, 'command': command}, indent=2)+'\n')
        raise SystemExit(rc)
    except (ValueError, OSError, yaml.YAMLError) as exc:
        ap.error(str(exc))

if __name__ == '__main__':
    main()
