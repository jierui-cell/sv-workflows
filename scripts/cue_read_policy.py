"""Narrow Cue compatibility policy; never invent missing sequence."""
import json

POLICY = 'skip_unmapped_missing_sequence_v1'

def compatible_reads(reads):
    seen = skipped = 0
    try:
        for read in reads:
            seen += 1
            if read.query_sequence is None:
                if not read.is_unmapped:
                    raise ValueError('Cue: mapped alignment has no sequence; review required')
                skipped += 1
                continue
            yield read
    finally:
        print('CUE_READ_POLICY ' + json.dumps(dict(policy=POLICY, seen=seen,
                                                   skipped_unmapped_missing_sequence=skipped)), flush=True)

def patch_io(text):
    # Fail closed if the pinned upstream iterator implementation changes.
    for old in ['enumerate(input_bam.fetch(chr_name))',
                'for read in input_bam.fetch(chr_name, start, end):']:
        if text.count(old) != 1:
            raise ValueError('Unexpected Cue iterator source: ' + old)
    text = text.replace('enumerate(input_bam.fetch(chr_name))',
                        'enumerate(compatible_reads(input_bam.fetch(chr_name)))')
    text = text.replace('for read in input_bam.fetch(chr_name, start, end):',
                        'for read in compatible_reads(input_bam.fetch(chr_name, start, end)):')
    return 'from cue_read_policy import compatible_reads\n' + text
