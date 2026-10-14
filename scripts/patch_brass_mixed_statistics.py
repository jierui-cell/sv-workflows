#!/usr/bin/env python3
"""Preserve numeric partial estimates; label mixed-origin statistics honestly."""
from pathlib import Path
import sys


def patch(source):
    anchor = "    $options->{'GenderChrFound'} = 'Y' unless(defined $options->{'GenderChrFound'});\n  }\n  return $ret_val;"
    if source.count(anchor) != 1:
        raise ValueError('Expected unique native statistics return block')
    extra = """
  if(defined $options->{'ascat_summary'}) {
    open my $META, '<', $options->{'ascat_summary'} or die $!;
    while(my $line = <$META>) {
      if($line =~ /^# BRASS_PIPELINE_STATE MIXED\\s*$/) {
        $options->{'PloidyAcfState'} = 'MIXED';
      }
    }
    close $META;
  }
"""
    return source.replace(anchor, anchor.replace('  return $ret_val;',extra+'  return $ret_val;'))


if __name__ == '__main__':
    Path(sys.argv[2]).write_text(patch(Path(sys.argv[1]).read_text()))
