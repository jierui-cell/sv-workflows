"""Exploratory breakpoint-read metrics on small, indexed IGV BAM slices."""
import csv
import json
import re
import subprocess
import sys
from pathlib import Path


def ref_end(start, cigar):
    consumed = sum(int(n) for n, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar)
                   if op in 'MDN=X')
    return start + consumed - 1


def near_endpoint(start, end, position, tolerance):
    return min(abs(start - position), abs(end - position)) <= tolerance


def sa_bridge(fields, breakpoints, tolerance, min_sa_mapq=20):
    tags = [x[5:] for x in fields[11:] if x.startswith('SA:Z:')]
    if not tags:
        return False
    chrom = fields[2]
    start = int(fields[3])
    end = ref_end(start, fields[5])
    for entry in tags[0].rstrip(';').split(';'):
        parts = entry.split(',')
        if len(parts) < 6 or parts[0] != chrom or int(parts[4]) < min_sa_mapq:
            continue
        sa_start = int(parts[1])
        sa_end = ref_end(sa_start, parts[3])
        a, b = breakpoints
        if ((near_endpoint(start, end, a, tolerance) and near_endpoint(sa_start, sa_end, b, tolerance))
                or (near_endpoint(start, end, b, tolerance) and near_endpoint(sa_start, sa_end, a, tolerance))):
            return True
    return False


def mean_depth(bam, chrom, position, radius=500, min_mapq=20, min_baseq=20):
    left = max(1, position - radius)
    right = position + radius
    result = subprocess.run(['samtools', 'depth', '-aa', '-q', str(min_baseq), '-Q', str(min_mapq),
                             '-r', '{}:{}-{}'.format(chrom, left, right), str(bam)],
                            check=True, universal_newlines=True,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    values = [int(row.split('\t')[2]) for row in result.stdout.splitlines()]
    if len(values) != right - left + 1:
        raise ValueError('Depth window incomplete for {}:{} in {}'.format(chrom, position, bam))
    return round(sum(values) / len(values), 2)


def analyze(bam, chrom, pos, end, tolerance=500, indexed_regions=False,
            fetch_radius=2000, depth_radius=500, min_mapq=20, min_baseq=20):
    command = ['samtools', 'view', '-F', '0xF04', '-q', str(min_mapq)]
    if indexed_regions:
        command.append('-M')
    command.append(str(bam))
    if indexed_regions:
        for point in (pos, end):
            command.append('{}:{}-{}'.format(chrom, max(1, point - fetch_radius), point + fetch_radius))
    result = subprocess.run(command,
                            check=True, universal_newlines=True,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    sa, discordant, clipped = set(), set(), set()
    passing = 0
    for line in result.stdout.splitlines():
        f = line.split('\t')
        if f[2] != chrom:
            continue
        passing += 1
        name = f[0]
        start = int(f[3])
        stop = ref_end(start, f[5])
        at_pos = near_endpoint(start, stop, pos, tolerance)
        at_end = near_endpoint(start, stop, end, tolerance)
        if not (at_pos or at_end):
            continue
        if sa_bridge(f, (pos, end), tolerance, min_mapq):
            sa.add(name)
        if (int(f[1]) & 0x1 and f[6] in ('=', chrom)
                and ((at_pos and abs(int(f[7]) - end) <= tolerance)
                     or (at_end and abs(int(f[7]) - pos) <= tolerance))):
            discordant.add(name)
        if re.search(r'\d+S', f[5]):
            clipped.add(name)
    return dict(primary_high_mapq_reads=passing,
                sa_bridge_unique_reads=len(sa),
                breakpoint_mate_link_unique_reads=len(discordant),
                junction_link_unique_reads=len(sa | discordant),
                softclip_near_breakpoint_unique_reads=len(clipped),
                depth_mean_start=mean_depth(bam, chrom, pos, depth_radius, min_mapq, min_baseq),
                depth_mean_end=mean_depth(bam, chrom, end, depth_radius, min_mapq, min_baseq))


def main():
    review = Path(sys.argv[1])
    out = Path(sys.argv[2])
    if out.exists():
        raise FileExistsError(out)
    rows = list(csv.DictReader((review / 'review_regions.tsv').open(), delimiter='\t'))
    report = []
    for row in rows:
        event = row['event']
        details = dict(row)
        for role in ('tumor', 'normal'):
            bam = review / event / (role + '.bam')
            if not Path(str(bam) + '.bai').exists():
                raise FileNotFoundError(str(bam) + '.bai')
            details[role] = analyze(bam, row['chrom'], int(row['pos']), int(row['end']))
        report.append(details)
    out.mkdir()
    result = dict(status='descriptive_read_review_complete',
                  policy='primary MAPQ>=20, exclude unmapped/secondary/QC-fail/duplicate/supplementary; SA MAPQ>=20; breakpoint tolerance 500 bp; depth baseQ>=20 MAPQ>=20',
                  scope='User-selected breakpoint-region BAM slices',
                  events=report,
                  caveat='Counts are exploratory read-level features, not variant-support probabilities. SA and mate links do not establish SV type or somatic origin. BAM slices cannot assess genome-wide context or full-interval copy number.')
    (out / 'summary.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
