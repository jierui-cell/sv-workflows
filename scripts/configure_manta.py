#!/usr/bin/env python
"""Merge validated overrides into the defaults shipped in the pinned image (Python 2/3)."""
from __future__ import print_function
import argparse
import json
import shutil
from distutils.spawn import find_executable
try:
    import configparser
except ImportError:
    import ConfigParser as configparser


def merge(defaults, overrides, output, reference=None):
    config = configparser.RawConfigParser()
    config.optionxform = str
    if not config.read(defaults):
        raise ValueError('Missing Manta default configuration')
    for key, value in overrides.items():
        if key == 'referenceFasta' or not config.has_option('manta', key):
            raise ValueError('Unknown or reserved Manta setting: ' + key)
        if isinstance(value, bool) or not isinstance(value, int) or value < 0:
            raise ValueError('Manta numeric overrides must be nonnegative integers: ' + key)
        if config.get('manta', key) in ('0', '1') and value not in (0, 1):
            raise ValueError('Manta boolean settings require 0 or 1: ' + key)
        config.set('manta', key, str(value))
    if reference is not None:
        config.set('manta', 'referenceFasta', reference)
    with open(output, 'w') as handle:
        config.write(handle)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('overrides')
    parser.add_argument('--output', default='manta.effective.ini')
    parser.add_argument('--reference')
    args = parser.parse_args()
    defaults = find_executable('configManta.py') + '.ini'
    shutil.copyfile(defaults, 'manta.defaults.ini')
    with open(args.overrides) as handle:
        merge(defaults, json.load(handle), args.output, args.reference)
