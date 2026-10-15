#!/usr/bin/env python3
"""Connect completed caller outputs to signatures without a hand-written VCF inventory."""
import csv
import hashlib
import io
import json
import os
from pathlib import Path
import re
import tempfile

from build_sv_manifest import FIELDS, CALLERS, validate_rows
from check_resources import absolute, normalize_paths, selected_samples


def one_file(paths, label):
    found = sorted({path.resolve() for path in paths if path.is_file()})
    if len(found) != 1:
        raise ValueError(f'{label}: expected one file, found {len(found)}; checked ' +
                         ', '.join(str(p) for p in paths))
    if not found[0].stat().st_size:
        raise ValueError(label + ': empty file ' + str(found[0]))
    return found[0]


def record(path):
    result = json.loads(path.read_text())
    if not isinstance(result, dict):
        raise ValueError('Expected JSON object: ' + str(path))
    return result


def discover(config, root):
    """Read only: select exact final outputs and verify task completion, never hash VCFs."""
    # Reuse the launcher's strict YAML loader/overlay rules at runtime.
    from workflow import merge, read_yaml
    source = absolute(config.get('caller_config'), root)
    defaults = Path(root) / 'config/defaults.yaml'
    caller = normalize_paths(merge(read_yaml(defaults), read_yaml(source)), 'callers', root)
    if caller['genome_build'] != config['genome_build'] or caller['sequencing_type'] != config['sequencing_type']:
        raise ValueError('Caller genome_build/sequencing_type differs from signature settings')
    if not isinstance(caller['callers'],str):
        raise ValueError('callers must be a comma-separated string in ' + str(source))
    requested = caller['callers'].split(',')
    if len(requested) != len(set(requested)) or set(requested) - (CALLERS | {'preflight'}):
        raise ValueError('Invalid or duplicate callers in ' + str(source))
    requested = sorted(set(requested) & CALLERS)
    if not requested:
        raise ValueError('Caller configuration selects no SV callers: ' + str(source))
    if not Path(caller['input']).is_file() and source.name == 'resolved.yaml':
        snapshot = source.parent / 'resolved-samplesheet.tsv'
        if snapshot.is_file():
            caller['input'] = str(snapshot)
    _, samples = selected_samples(caller, root)  # Metadata only; no BAMs/indexes are opened.
    cohort = caller['dataset_id']
    if not isinstance(cohort,str) or not re.fullmatch(r'[A-Za-z0-9_-]+',cohort):
        raise ValueError('Signature cohort identifiers must contain only letters, digits, underscores or hyphens')
    if any(not re.fullmatch(r'[A-Za-z0-9_-]+',row['patient_id']) for row in samples):
        raise ValueError('Signature pair identifiers must contain only letters, digits, underscores or hyphens')
    all_groups = {cohort + ':' + name for name in requested}
    groups = config.get('groups', 'all')
    if not isinstance(groups, str):
        raise ValueError('groups must be all or comma-separated cohort:caller names')
    groups = sorted(all_groups) if groups == 'all' else groups.split(',')
    if not groups or len(set(groups)) != len(groups) or set(groups) - all_groups:
        raise ValueError('Unknown or duplicate signature groups; available: ' + ', '.join(sorted(all_groups)))
    if any(not isinstance(config.get(key,''),str) for key in ('pair_id','pair_ids')):
        raise ValueError('pair_id and pair_ids must be strings')
    if config.get('pair_id') and config.get('pair_ids'):
        raise ValueError('Use pair_id or pair_ids, not both')
    ids = [config['pair_id']] if config.get('pair_id') else (config.get('pair_ids', '').split(',') if config.get('pair_ids') else [])
    if len(set(ids)) != len(ids) or set(ids) - {row['patient_id'] for row in samples}:
        raise ValueError('Unknown or duplicate signature pair selection')
    samples = [row for row in samples if not ids or row['patient_id'] in ids]
    rows, issues = [], []
    for name in requested:
        if cohort + ':' + name not in groups:
            continue
        if name == 'gridss' and caller['gridss']['somatic_filter'] is not True:
            raise ValueError('GRIDSS signature inputs require somatic_filter=true; joint VCFs are not somatic inputs')
        for sample in sorted(samples, key=lambda row: row['patient_id']):
            pair = sample['patient_id']
            folder = Path(caller['outdir']) / name / pair
            label = f'{cohort}:{name}:{pair}'
            try:
                benchmark_path = one_file([folder / 'benchmark.json', folder / 'output/benchmark.json'], label + ' benchmark')
                benchmark = record(benchmark_path)
                if benchmark.get('schema') != 'sv-caller-benchmark/v1' or type(benchmark.get('exit_code')) is not int or benchmark['exit_code'] != 0:
                    raise ValueError('caller task is incomplete, failed or a stub: ' + str(benchmark_path))
                for key, expected in (('dataset_id', cohort), ('caller', name), ('pair_id', pair)):
                    if benchmark.get(key) != expected:
                        raise ValueError('benchmark ' + key + ' does not match the caller configuration')
                if Path(benchmark.get('reference', '')).resolve() != Path(caller['reference']).resolve():
                    raise ValueError('benchmark reference differs from the caller configuration')
                settings = benchmark.get('settings', {})
                if not isinstance(settings,dict):
                    raise ValueError('benchmark settings must be a JSON object')
                if settings.get('genome_build', caller['genome_build']) != config['genome_build'] or settings.get('sequencing_type', caller['sequencing_type']) != config['sequencing_type']:
                    raise ValueError('benchmark build/assay differs from signature settings')
                # Manta keeps its native VCF tree at the pair root; its benchmark
                # is now in output/, with historical root-level records supported.
                task = folder if name == 'manta' else benchmark_path.parent
                patterns = {
                    'manta': ('manta/results/variants/somaticSV.vcf.gz', 'somaticSV.vcf.gz'),
                    'delly': ('somaticSV.vcf.gz',),
                    'svaba': (pair + '.svaba.somatic.sv.vcf',),
                    'gridss': ('somatic/somatic.high_confidence.vcf.gz',),
                    'brass': ('callset/calls/*.annot.vcf.gz',),
                }[name]
                candidates = [p for pattern in patterns for p in (list(task.glob(pattern)) or [task / pattern])]
                vcf = one_file(candidates, label + ' final SV VCF')
                if name == 'manta':
                    if benchmark.get('output_class') != 'somatic':
                        raise ValueError('Manta benchmark does not describe somatic outputs')
                elif name == 'brass':
                    qc = record(task / 'callset/technical-qc.json')
                    if qc.get('status') != 'technical_qc_passed':
                        raise ValueError('BRASS technical output QC has not passed')
                else:
                    success = record(task / 'SUCCESS.json')
                    expected_stage = 'somatic_filtered_no_pon' if name == 'gridss' else 'somatic'
                    if success.get('stage') != expected_stage or not isinstance(success.get('vcfs'),list) or str(vcf.relative_to(task.resolve())) not in success['vcfs']:
                        raise ValueError('final VCF is not registered as a completed somatic output')
                    if name == 'gridss' and record(task / 'somatic/SUCCESS.json').get('stage') != expected_stage:
                        raise ValueError('GRIDSS native somatic filtering has not completed')
                rows.append(dict(cohort=cohort, pair_id=pair, caller=name, vcf=str(vcf),
                                 source_status='pipeline_completed_technical_qc'))
            except (OSError, ValueError, KeyError, TypeError) as exc:
                issues.append(label + ': ' + str(exc))
    if issues:
        raise ValueError('Caller outputs are not ready for signatures:\n' + '\n'.join(issues))
    rows = validate_rows(rows, hash_sources=False)
    provenance = {'mode': 'caller_outputs', 'caller_config': str(source),
                  'caller_config_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
                  'caller_defaults_sha256': hashlib.sha256(defaults.read_bytes()).hexdigest(),
                  'caller_outdir': caller['outdir'], 'cohort': cohort,
                  'samplesheet': caller['input'], 'samplesheet_sha256': hashlib.sha256(Path(caller['input']).read_bytes()).hexdigest(),
                  'genome_build': caller['genome_build'], 'sequencing_type': caller['sequencing_type'],
                  'groups': groups, 'callsets': len(rows), 'gridss_input': 'somatic.high_confidence.vcf.gz'}
    return rows, provenance


def write_auto_manifest(config, root, log_dir):
    """Hash final VCFs and save a stable input path plus a per-run snapshot."""
    rows, provenance = discover(config, root)
    rows = validate_rows(rows, hash_sources=True)
    handle = io.StringIO()
    writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter='\t', lineterminator='\n')
    writer.writeheader(); writer.writerows(rows)
    content = handle.getvalue()
    checksum = hashlib.sha256(content.encode()).hexdigest()
    cache = Path(root) / '.nextflow/signature-manifests' / (checksum + '.tsv')
    cache.parent.mkdir(parents=True, exist_ok=True)
    if not cache.exists():
        with tempfile.NamedTemporaryFile(mode='w', dir=cache.parent, delete=False) as temp:
            temp.write(content)
            temporary = Path(temp.name)
        try:
            try:
                os.link(temporary, cache)
            except FileExistsError:
                pass
        finally:
            temporary.unlink()
    if cache.read_text() != content:
        raise ValueError('Conflicting signature manifest cache: ' + str(cache))
    (Path(log_dir) / 'input_manifest.tsv').write_text(content)
    provenance['manifest_sha256'] = checksum
    (Path(log_dir) / 'signature-inputs.json').write_text(json.dumps(provenance, indent=2) + '\n')
    return cache
