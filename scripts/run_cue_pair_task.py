#!/usr/bin/env python3
"""Sequential Cue tumor+normal inference and paired research screen in one task."""
import json
from pathlib import Path
import shutil
import subprocess
import sys


def main():
    spec = json.loads(Path(sys.argv[1]).read_text())
    root = Path(spec['root'])
    pair = spec['id']
    settings = spec['settings']
    if settings.get('pair_execution') != 'sequential':
        raise ValueError('Cue pair_execution must be sequential for this benchmark profile')
    output = Path('output')
    output.mkdir(exist_ok=True)
    (output / 'effective-task.json').write_text(json.dumps(spec, indent=2) + '\n')
    for role in ('tumor', 'normal'):
        work = Path(role + '-work')
        work.mkdir()
        role_spec = dict(id=pair, role=role, bam=spec[role + '_bam'],
                         fai=spec['fai'], settings=settings, root=str(root))
        task_json = work / 'task.json'
        task_json.write_text(json.dumps(role_spec, indent=2) + '\n')
        with (output / (role + '.driver.log')).open('w') as log:
            subprocess.run([sys.executable, str(root / 'scripts/run_cue.py'),
                            str(task_json.resolve())], cwd=work, stdout=log,
                           stderr=subprocess.STDOUT, check=True)
        source = work / 'output'
        if not (source / 'SUCCESS.json').is_file():
            raise ValueError('Cue ' + role + ' did not produce SUCCESS.json')
        shutil.move(str(source), str(output / role))
    paired = output / 'paired'
    with (output / 'paired.driver.log').open('w') as log:
        subprocess.run([sys.executable, str(root / 'scripts/run_cue_paired_screen.py'),
                        '--tumor', str((output / 'tumor/single_sample.svs.vcf').resolve()),
                        '--normal', str((output / 'normal/single_sample.svs.vcf').resolve()),
                        '--pair', pair, '--outdir', str(paired.resolve()),
                        '--comparison-config', str((root / settings['comparison_config']).resolve()),
                        '--bam-config', str((root / settings['bam_screen_config']).resolve())],
                       stdout=log, stderr=subprocess.STDOUT, check=True)
    success = json.loads((paired / 'SUCCESS.json').read_text())
    (output / 'SUCCESS.json').write_text(json.dumps({
        'stage': 'cue_pair_research_only_unvalidated', 'pair': pair,
        'tumor_and_normal_inference': 'complete', 'paired_screen': success}, indent=2) + '\n')


if __name__ == '__main__':
    main()
