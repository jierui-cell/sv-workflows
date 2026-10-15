#!/usr/bin/env python3
"""Download, prepare and validate BRASS's extra resources against the shared FASTA."""
import argparse
import csv
from datetime import datetime, timezone
import fcntl
import hashlib
import json
from pathlib import Path, PurePosixPath
import shutil
import subprocess
import sys
import tarfile
import tempfile
import yaml

from adapt_brass_resources import dictionary, mapping
from check_resources import BRASS_FILES, normalize_paths
from prepare_brass_bundle import validate_adapted
from workflow import merge, override, read_yaml

ROOT = Path(__file__).resolve().parents[1]
BRASS = 'CNV_SV_ref/brass/'
VAGRENT = 'VAGrENT_ref_GRCh37d5_ensembl_75/vagrent/'
ARCHIVES = {
    'CNV_SV_ref_GRCh37d5_brass6+.tar.gz': tuple(BRASS + name for name in (
        'CentTelo.tsv', 'cytoband.txt', 'HiDepth.bed.gz', '500bp_windows.gc.bed.gz',
        'viral.genomic.fa.2bit', *(f'all_ncbi_bacteria.{i}.fa.2bit' for i in range(1, 5)))),
    'VAGrENT_ref_GRCh37d5_ensembl_75.tar.gz': tuple(VAGRENT + name for name in (
        'vagrent.cache.gz', 'vagrent.fa', 'vagrent.fa.fai')),
}


def digest(path, algorithm='sha256'):
    h = hashlib.new(algorithm)
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def fetch(row, folder, run):
    """Resume partial downloads; accept a completed file only after checksum verification."""
    target = folder / row['file']
    cached = target.exists()
    candidate = target if cached else target.with_name(target.name + '.partial')
    if not cached:
        run(['curl', '--fail', '--location', '--retry', '3', '--connect-timeout', '30',
             '--continue-at', '-', '--output', str(candidate), row['url']])
    if digest(candidate, 'md5') != row['md5']:
        label = 'Cached archive' if cached else 'Download'
        raise ValueError(label + ' checksum mismatch; inspect retained file: ' + str(candidate))
    checksum = digest(candidate)
    if checksum != row['sha256']:
        raise ValueError('Archive SHA-256 checksum mismatch; inspect retained file: ' + str(candidate))
    if not cached:
        candidate.rename(target)
    return {'file': row['file'], 'url': row['url'], 'md5': row['md5'],
            'sha256': checksum, 'bytes': target.stat().st_size}


def extract_selected(archive, destination, required):
    """Extract only needed regular files, never tar links, devices or traversing paths."""
    required, found = set(required), set()
    with tarfile.open(archive, 'r|gz') as stream:
        for member in stream:
            name = PurePosixPath(member.name)
            if name.is_absolute() or '..' in name.parts:
                raise ValueError('Unsafe archive path: ' + member.name)
            key = str(name)
            if key not in required:
                continue
            if not member.isfile() or key in found:
                raise ValueError('Expected one regular archive member: ' + key)
            target = destination / key
            target.parent.mkdir(parents=True, exist_ok=True)
            with stream.extractfile(member) as source, target.open('xb') as output:
                shutil.copyfileobj(source, output)
            found.add(key)
    if found != required:
        raise ValueError('Missing archive members: ' + ', '.join(sorted(required - found)))


def prepare_source(downloads, source, records, source_dict):
    expected = {name for names in ARCHIVES.values() for name in names} | {'reference.dict'}
    if source.exists():
        receipt = json.loads((source / 'source-receipt.json').read_text())
        if receipt.get('archives') != records or set(receipt.get('files', {})) != expected:
            raise ValueError('Existing source inventory differs: ' + str(source))
        for name, checksum in receipt['files'].items():
            if digest(source / name) != checksum:
                raise ValueError('Extracted source changed: ' + str(source / name))
        if digest(source / 'reference.dict') != digest(source_dict):
            raise ValueError('Existing source dictionary differs from the versioned dictionary')
        return
    temporary = Path(tempfile.mkdtemp(prefix='.brass-source-', dir=source.parent))
    try:
        for filename, names in ARCHIVES.items():
            extract_selected(downloads / filename, temporary, names)
        shutil.copyfile(source_dict, temporary / 'reference.dict')
        receipt = {'archives': records, 'files': {
            name: digest(temporary / name) for name in sorted(expected)}}
        (temporary / 'source-receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
        temporary.rename(source)
    finally:
        if temporary.exists():
            shutil.rmtree(temporary)


def check_existing_bundle(bundle, reference, image):
    """Reuse a validated bundle only when its resources, FASTA, image and overlay still match."""
    if not (bundle / 'COMPATIBILITY_VALIDATED').is_file():
        raise ValueError('Incomplete bundle at ' + str(bundle) + '; choose a new brass.bundle path')
    for name in BRASS_FILES:
        if not (bundle / name).is_file() or not (bundle / name).stat().st_size:
            raise ValueError('Missing/empty bundle resource: ' + str(bundle / name))
    receipt = json.loads((bundle / 'bundle-receipt.json').read_text())
    native = json.loads((bundle / 'native-validation.json').read_text())
    if receipt.get('status') != 'passed' or receipt.get('assembly') != 'GRCh37' or native.get('status') != 'passed':
        raise ValueError('Bundle validation did not pass')
    checks = {'reference_sha256': reference, 'reference_fai_sha256': Path(str(reference) + '.fai'),
              'target_dictionary_sha256': reference.with_suffix('.dict'), 'image_sha256': image,
              'overlay_sha256': bundle / 'get_rg_cns.R'}
    for key, path in checks.items():
        if receipt.get(key) != digest(path):
            raise ValueError('Existing bundle differs from current input: ' + key)
    validate_adapted(bundle, reference.with_suffix('.dict'))


def select_engine(preferred, explicit=False):
    """Detect an installed SIF runtime without overriding an explicit user choice."""
    if preferred not in ('apptainer', 'singularity'):
        raise ValueError('container_engine must be apptainer or singularity')
    candidates = (preferred,) if explicit else (preferred, 'singularity' if preferred == 'apptainer' else 'apptainer')
    for engine in candidates:
        if shutil.which(engine):
            return engine
    if explicit:
        raise ValueError(f'Requested container runtime is unavailable: {preferred}. '
                         'Load its module or choose --engine singularity / --engine apptainer.')
    raise ValueError('No container runtime found. Load an Apptainer or Singularity module, '
                     'then rerun; --engine singularity and --engine apptainer are supported.')


def setup(config, root=ROOT, engine_explicit=True):
    config = normalize_paths(config, 'callers', root)
    reference = Path(config['reference']).resolve()
    image = Path(config['brass']['image'])
    bundle = Path(config['brass']['bundle'])
    source_dict = root / 'config/brass-GRCh37d5.dict'
    if config['genome_build'] != 'GRCh37' or config['brass']['assembly'] != 'GRCh37':
        raise ValueError('The downloadable BRASS bundle is for GRCh37; use a GRCh37 configuration')
    for path in (reference, Path(str(reference) + '.fai'), reference.with_suffix('.dict'), source_dict):
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError('Missing reference prerequisite: ' + str(path))
    # Refuse incompatible builds before any large download. Native checks also
    # verify the supported chr1-22/chrX/chrY calling scope during preparation.
    mapping(dictionary(source_dict), dictionary(reference.with_suffix('.dict')))
    engine = select_engine(config['container_engine'], explicit=engine_explicit)
    print('Container runtime: ' + engine, flush=True)
    for command in ('curl', 'patch'):
        if not shutil.which(command):
            raise ValueError('Load/install required setup tool: ' + command)
    folder = root / 'resources/brass'
    folder.mkdir(parents=True, exist_ok=True)
    with (folder / '.setup.lock').open('a') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise ValueError('Another BRASS setup is active in this checkout') from None
        log_dir = root / 'logs' / ('brass-setup-' + datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S.%fZ'))
        log_dir.mkdir(parents=True)
        step, commands = 'start', []
        print('Setup log: ' + str(log_dir), flush=True)
        with (log_dir / 'setup.log').open('w') as log:
            def run(argv):
                commands.append(argv)
                print('Running: ' + ' '.join(argv), flush=True)
                log.write(json.dumps(argv) + '\n'); log.flush()
                subprocess.run(argv, check=True, stdout=log, stderr=subprocess.STDOUT)
            try:
                step = 'caller_image'
                if not image.is_file():
                    default_image = root / 'containers/images/brass-v6.3.4.sif'
                    if image != default_image:
                        raise ValueError('Configured custom BRASS image is absent: ' + str(image))
                    run([sys.executable, str(root / 'scripts/fetch_containers.py'),
                         '--engine', engine, '--tools', 'brass'])
                step = 'existing_bundle'
                if bundle.exists():
                    check_existing_bundle(bundle, reference, image)
                    outcome = 'reused_verified_bundle'
                else:
                    step = 'downloads'
                    downloads = folder / 'downloads'
                    downloads.mkdir(exist_ok=True)
                    with (root / 'config/brass-resources.tsv').open() as handle:
                        rows = list(csv.DictReader(handle, delimiter='\t'))
                    if {row['file'] for row in rows} != set(ARCHIVES) or len(rows) != len(ARCHIVES):
                        raise ValueError('Unexpected BRASS archive manifest')
                    records = [fetch(row, downloads, run) for row in rows]
                    step = 'extraction'
                    source = folder / 'source'
                    prepare_source(downloads, source, records, source_dict)
                    step = 'adaptation_and_native_validation'
                    bundle.parent.mkdir(parents=True, exist_ok=True)
                    staging = Path(tempfile.mkdtemp(prefix='.brass-bundle-', dir=bundle.parent))
                    prepared = staging / 'bundle'
                    # Publish only a complete bundle. A failed preparation retains
                    # its staged files/logs; the same setup command can be retried.
                    run([sys.executable, str(root / 'scripts/prepare_brass_bundle.py'),
                         '--source', str(source), '--source-dict', str(source / 'reference.dict'),
                         '--reference', str(reference), '--image', str(image),
                         '--output', str(prepared), '--engine', engine])
                    check_existing_bundle(prepared, reference, image)
                    prepared.rename(bundle)
                    staging.rmdir()
                    outcome = 'prepared_and_validated'
                report = {'status': 'passed', 'outcome': outcome, 'container_engine': engine, 'reference': str(reference),
                          'bundle': str(bundle), 'commands': commands}
            except Exception as exc:
                report = {'status': 'failed', 'step': step, 'container_engine': engine, 'error': str(exc), 'commands': commands}
                (log_dir / 'status.json').write_text(json.dumps(report, indent=2) + '\n')
                raise ValueError(f'BRASS setup failed at {step}: {exc}; log: {log_dir}/setup.log') from exc
            (log_dir / 'status.json').write_text(json.dumps(report, indent=2) + '\n')
    print('BRASS resources ready: ' + str(bundle))
    print('Shared human reference (unchanged): ' + str(reference))
    return report


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--config', type=Path, action='append', default=[], help='Same YAML overlays as the caller launcher')
    ap.add_argument('--set', action='append', default=[], metavar='KEY=VALUE')
    ap.add_argument('--engine', choices=('apptainer', 'singularity'),
                    help='Choose runtime; otherwise honor YAML/--set, or detect Apptainer then Singularity')
    args = ap.parse_args()
    try:
        config = read_yaml(ROOT / 'config/defaults.yaml')
        engine_explicit = False
        for path in args.config:
            overlay = read_yaml(path)
            engine_explicit = engine_explicit or 'container_engine' in overlay
            config = merge(config, overlay)
        for item in args.set:
            override(config, item)
            engine_explicit = engine_explicit or item.partition('=')[0] == 'container_engine'
        if args.engine:
            config['container_engine'] = args.engine
            engine_explicit = True
        setup(config, engine_explicit=engine_explicit)
    except (ValueError, OSError, KeyError, yaml.YAMLError) as exc:
        ap.error(str(exc))


if __name__ == '__main__':
    main()
