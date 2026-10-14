"""Validate generic GRIDSS labels against the actual executed input order."""
from pathlib import Path

def validate_identity(row, spec, commands, labels):
    expected=[str(Path(row[k+'_bam']).resolve()) for k in ('normal','tumor')]
    if expected[0]==expected[1]: raise ValueError('Identical input paths')
    for role,path in zip(('normal','tumor'),expected):
        if Path(spec[role]).name != Path(path).name:
            raise ValueError('Resolved settings disagree with manifest: '+role)
    calls=[x['argv'] for x in commands if x.get('step')=='raw_joint_call']
    if len(calls)!=1: raise ValueError('Expected exactly one raw GRIDSS command')
    argv=calls[0]
    if argv[-2:]!=expected: raise ValueError('Executed normal/tumour order mismatch')
    if argv[argv.index('--labels')+1]!='NORMAL,TUMOR' or labels!=['NORMAL','TUMOR']:
        raise ValueError('Unexpected label mapping')
    return dict(zip(labels,expected))
