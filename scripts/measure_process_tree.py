#!/usr/bin/env python3
"""Sample concurrent RSS of a caller's process tree, without counting page cache.

RSS sums shared pages more than once; this is an observational upper-biased
estimate, not Slurm cgroup MaxRSS, PSS, or a guaranteed instantaneous maximum.
"""
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time


def process_table():
    result = {}
    with os.scandir('/proc') as entries:
        for entry in entries:
            if not entry.name.isdigit():
                continue
            try:
                pid = int(entry.name)
                with open(f'/proc/{pid}/status') as handle:
                    fields = {}
                    for line in handle:
                        if line.startswith(('PPid:', 'VmRSS:')):
                            key, value = line.split(':', 1)
                            fields[key] = int(value.split()[0])
                result[pid] = (fields.get('PPid', 0), fields.get('VmRSS', 0))
            except (FileNotFoundError, ProcessLookupError, PermissionError, ValueError):
                continue
    return result


def descendant_rss(table, root):
    children = {}
    for pid, (parent, _) in table.items():
        children.setdefault(parent, []).append(pid)
    seen, pending = set(), [root]
    while pending:
        pid = pending.pop()
        if pid in seen:
            continue
        seen.add(pid)
        pending.extend(children.get(pid, ()))
    rss = [table[pid][1] for pid in seen if pid in table]
    return sum(rss), max(rss, default=0), len(rss)


def main():
    if len(sys.argv) < 3:
        raise SystemExit('Usage: measure_process_tree.py OUTPUT COMMAND [ARG ...]')
    output, command = Path(sys.argv[1]), sys.argv[2:]
    started = time.monotonic()
    child = subprocess.Popen(command)
    def forward(signum, _frame):
        if child.poll() is None:
            child.send_signal(signum)
    signal.signal(signal.SIGTERM, forward)
    signal.signal(signal.SIGINT, forward)
    peak_sum = peak_single = peak_processes = samples = 0
    peak_elapsed = 0.0
    try:
        while True:
            table = process_table()
            total, single, count = descendant_rss(table, child.pid)
            samples += 1
            if total > peak_sum:
                peak_sum, peak_processes = total, count
                peak_elapsed = time.monotonic() - started
            peak_single = max(peak_single, single)
            if child.poll() is not None:
                break
            if peak_sum == 0:
                # The first snapshot can land before exec has populated VmRSS.
                time.sleep(0.05)
                continue
            time.sleep(2)
    finally:
        code = child.wait()
        report = {
            'measurement': 'sampled /proc VmRSS of command process tree',
            'command': command,
            'elapsed_seconds': time.monotonic() - started,
            'sampling_interval_seconds': 2,
            'samples': samples,
            'max_concurrent_tree_rss_kib': peak_sum,
            'max_single_process_rss_kib': peak_single,
            'processes_at_tree_peak': peak_processes,
            'tree_peak_elapsed_seconds': peak_elapsed,
            'exit_code': code,
            'limitations': 'RSS can double-count shared pages; 2-second sampling can miss short peaks; descendants that detach/reparent are not included; this is not cgroup MaxRSS or PSS',
        }
        output.parent.mkdir(parents=True, exist_ok=True)
        temp = output.with_name(output.name + '.tmp')
        temp.write_text(json.dumps(report, indent=2) + '\n')
        os.replace(temp, output)
    raise SystemExit(code)


if __name__ == '__main__':
    main()
