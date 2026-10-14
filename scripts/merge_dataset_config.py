#!/usr/bin/env python3
"""Merge shared caller settings with small dataset/run YAML overlays."""

import argparse
from pathlib import Path

import yaml


def merge(base, override):
    result = dict(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(result.get(key), dict):
            result[key] = merge(result[key], value)
        else:
            result[key] = value
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, required=True)
    parser.add_argument("--overlay", type=Path, action="append", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    with args.base.open(encoding="utf-8") as handle:
        config = yaml.safe_load(handle)
    if not isinstance(config, dict):
        raise ValueError("Base configuration must be a YAML mapping")
    for path in args.overlay:
        with path.open(encoding="utf-8") as handle:
            overlay = yaml.safe_load(handle)
        if not isinstance(overlay, dict):
            raise ValueError(f"Overlay must be a YAML mapping: {path}")
        config = merge(config, overlay)
    for key in ("input", "dataset_id", "outdir", "callers", "reference"):
        if not config.get(key):
            raise ValueError(f"Missing required resolved setting: {key}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(yaml.safe_dump(config, default_flow_style=False), encoding="utf-8")


if __name__ == "__main__":
    main()
