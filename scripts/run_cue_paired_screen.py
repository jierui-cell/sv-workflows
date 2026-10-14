"""One-command, research-only Cue paired VCF comparison and BAM screen."""
import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path


def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def main():
    project = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--tumor', type=Path, required=True, help='Cue tumor single_sample.svs.vcf')
    ap.add_argument('--normal', type=Path, required=True, help='Cue normal single_sample.svs.vcf')
    ap.add_argument('--pair', required=True)
    ap.add_argument('--outdir', type=Path, required=True)
    ap.add_argument('--comparison-config', type=Path,
                    default=project / 'config/cue-pair-comparison.yaml')
    ap.add_argument('--bam-config', type=Path,
                    default=project / 'config/cue-bam-screen.yaml')
    args = ap.parse_args()
    for path in (args.tumor, args.normal, args.comparison_config, args.bam_config):
        if not path.is_file():
            ap.error('Missing input: {}'.format(path))
    args.outdir.mkdir(parents=True, exist_ok=False)
    scripts = [project / 'scripts/cue_pair_candidates.py',
               project / 'scripts/cue_bam_screen.py',
               project / 'scripts/summarize_cue_review.py', Path(__file__).resolve()]
    manifest = dict(pair=args.pair, tumor_vcf=str(args.tumor.resolve()),
                    normal_vcf=str(args.normal.resolve()),
                    code_sha256={str(p): sha256(p) for p in scripts},
                    config_sha256={str(p.resolve()): sha256(p) for p in
                                   (args.comparison_config, args.bam_config)})
    (args.outdir / 'invocation.json').write_text(json.dumps(manifest, indent=2) + '\n')
    comparison = args.outdir / 'comparison'
    screen = args.outdir / 'bam_screen'
    with (args.outdir / 'comparison.stdout.json').open('w') as log:
        subprocess.run([sys.executable, str(scripts[0]), '--tumor', str(args.tumor),
                        '--normal', str(args.normal), '--pair', args.pair,
                        '--config', str(args.comparison_config), '--outdir', str(comparison)],
                       check=True, stdout=log)
    with (args.outdir / 'bam_screen.stdout.json').open('w') as log:
        subprocess.run([sys.executable, str(scripts[1]), '--pair-summary',
                        str(comparison / 'summary.json'), '--config', str(args.bam_config),
                        '--outdir', str(screen)], check=True, stdout=log)
    summary = json.loads((screen / 'summary.json').read_text())
    if summary['status'] != 'research_bam_screen_complete':
        raise ValueError('BAM screen did not complete')
    for name in ('cue_candidates.with_bam_evidence.vcf',
                 'cue_provisional_somatic.unvalidated.vcf'):
        subprocess.run(['bcftools', 'view', '-Ov', '-o', '/dev/null', str(screen / name)],
                       check=True)
    (args.outdir / 'SUCCESS.json').write_text(json.dumps({
        'stage': 'research_only_unvalidated_somatic_candidates',
        'pair': args.pair, 'input_candidates': summary['input_candidates'],
        'provisional_candidates': summary['provisional_candidates']}, indent=2) + '\n')


if __name__ == '__main__':
    main()
