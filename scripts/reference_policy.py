#!/usr/bin/env python3
"""Assembly checks for human references. Contig lengths are not sequence identity."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import re

# Primary chromosome lengths from GRCh37 and GRCh38 reference dictionaries.
LENGTHS = {
    'GRCh37': (249250621,243199373,198022430,191154276,180915260,171115067,
               159138663,146364022,141213431,135534747,135006516,133851895,
               115169878,107349540,102531392,90354753,81195210,78077248,
               59128983,63025520,48129895,51304566,155270560),
    'GRCh38': (248956422,242193529,198295559,190214555,181538259,170805979,
               159345973,145138636,138394717,133797422,135086622,133275309,
               114364328,107043718,101991189,90338345,83257441,80373285,
               58617616,64444167,46709983,50818468,156040895),
}

def primary_name(name):
    return name[3:] if name.startswith('chr') else name

def check_dictionary(contigs, build, require_primary=True):
    if build not in LENGTHS:
        raise ValueError('Supported genome_build values are GRCh37 and GRCh38')
    expected = dict(zip([str(i) for i in range(1,23)] + ['X'], LENGTHS[build]))
    seen = set()
    for name, length in contigs.items():
        key = primary_name(name)
        if key in expected:
            if key in seen:
                raise ValueError('Mixed/duplicate primary chromosome aliases: ' + key)
            seen.add(key)
            if int(length) != expected[key]:
                raise ValueError(f'{build} mismatch: {name} has length {length}, expected {expected[key]}')
    if require_primary and seen != set(expected):
        raise ValueError('Reference missing primary chromosomes: ' + ','.join(sorted(set(expected)-seen)))
    if not seen:
        raise ValueError('No primary chromosome lengths available to verify the claimed assembly')
    return sorted(seen)

def read_fai(path):
    contigs = {}
    for line in Path(path).read_text().splitlines():
        fields = line.split('\t')
        if len(fields) < 2 or fields[0] in contigs or int(fields[1]) <= 0:
            raise ValueError('Invalid/duplicate FASTA index entry: ' + line)
        contigs[fields[0]] = int(fields[1])
    if not contigs:
        raise ValueError('Empty FASTA index')
    return contigs

def check_vcf(path, build):
    opener = gzip.open if str(path).endswith(('.gz','.bgz')) else open
    contigs = {}
    with opener(path, 'rt') as handle:
        for line in handle:
            if not line.startswith('##'):
                break
            if line.startswith('##contig=<'):
                name = re.search(r'(?:<|,)ID=([^,>]+)', line)
                length = re.search(r'(?:<|,)length=(\d+)', line, re.I)
                if name and length:
                    if name[1] in contigs:
                        raise ValueError('Duplicate VCF contig: ' + name[1])
                    contigs[name[1]] = int(length[1])
    return check_dictionary(contigs, build, require_primary=False)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fai', required=True, type=Path)
    parser.add_argument('--build', required=True)
    parser.add_argument('--allow-synthetic', action='store_true')
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    contigs = read_fai(args.fai)
    synthetic = args.build == 'synthetic' and args.allow_synthetic
    checked = [] if synthetic else check_dictionary(contigs, args.build)
    args.output.write_text(json.dumps({'genome_build':args.build,
        'status':'synthetic_stub_only' if synthetic else 'dictionary_lengths_verified',
        'primary_chromosomes_checked':checked,
        'fai_sha256':hashlib.sha256(args.fai.read_bytes()).hexdigest(),
        'limitation':'Length checks cannot distinguish references with identical dictionaries; supply the exact alignment reference.'}, indent=2)+'\n')

if __name__ == '__main__':
    main()
