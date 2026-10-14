#!/usr/bin/env python3
"""Uniform post-preflight SV caller benchmark: wall clock and sampled memory.

The primary memory metric is cgroup charged resident memory, excluding unmapped
streaming page cache. Process PSS is supplementary because setuid container
helpers can make their smaps inaccessible. Neither is an exact high-water mark.
"""
import hashlib
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time
from datetime import datetime, timezone

SCHEMA = 'sv-caller-benchmark/v1'
INTERVAL_SECONDS = 1.0
BOUNDARY = ('after shared pair preflight; includes caller-specific preparation, '
            'tumor/normal inference where applicable, native/research filtering, '
            'final output QC; excludes Slurm queue wait and common pair preflight')


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def session_pids(session_id):
    found = []
    with os.scandir('/proc') as entries:
        for entry in entries:
            if not entry.name.isdigit():
                continue
            try:
                stat = (Path(entry.path) / 'stat').read_text()
                fields = stat[stat.rfind(')') + 2:].split()
                if int(fields[3]) == session_id and fields[0] not in ('Z', 'X', 'x'):
                    found.append(int(entry.name))
            except (FileNotFoundError, ProcessLookupError, PermissionError,
                    ValueError, IndexError):
                continue
    return found


def rollup_kib(pid):
    values = {}
    try:
        with open(f'/proc/{pid}/smaps_rollup') as handle:
            for line in handle:
                if line.startswith(('Rss:', 'Pss:')):
                    key, number = line.split(':', 1)
                    values[key] = int(number.split()[0])
    except (FileNotFoundError, ProcessLookupError):
        return None  # Process exited between discovery and measurement.
    return values


def sample(session_id):
    total_pss = total_rss = count = 0
    unavailable = []
    for pid in session_pids(session_id):
        try:
            values = rollup_kib(pid)
        except PermissionError as exc:
            unavailable.append('{}:{}'.format(pid, type(exc).__name__))
            continue
        if values is None:
            continue
        if 'Pss' not in values or 'Rss' not in values:
            # A process may become a zombie after its stat was inspected.
            try:
                stat = Path('/proc/{}/stat'.format(pid)).read_text()
                state = stat[stat.rfind(')') + 2:].split()[0]
            except (FileNotFoundError, ProcessLookupError):
                continue
            if state not in ('Z', 'X', 'x'):
                unavailable.append('{}:missing_rollup'.format(pid))
            continue
        total_pss += values['Pss']
        total_rss += values['Rss']
        count += 1
    return total_pss, total_rss, count, unavailable


def memory_cgroup():
    v2 = None
    for line in Path('/proc/self/cgroup').read_text().splitlines():
        parts = line.split(':', 2)
        if len(parts) != 3:
            continue
        if 'memory' in parts[1].split(','):
            return ('cgroup_v1_total_rss_plus_total_mapped_file',
                    Path('/sys/fs/cgroup/memory') / parts[2].lstrip('/'))
        if parts[0] == '0' and parts[1] == '':
            v2 = ('cgroup_v2_anon_plus_file_mapped',
                  Path('/sys/fs/cgroup') / parts[2].lstrip('/'))
    if v2:
        return v2
    raise ValueError('No memory cgroup found for this Slurm task')


def cgroup_resident_kib(method, location):
    values = {}
    for line in (location / 'memory.stat').read_text().splitlines():
        key, number = line.split()
        values[key] = int(number)
    if method == 'cgroup_v1_total_rss_plus_total_mapped_file':
        keys = ('total_rss', 'total_mapped_file')
    else:
        keys = ('anon', 'file_mapped')
    if any(key not in values for key in keys):
        raise ValueError('Missing cgroup resident counters: ' + ', '.join(keys))
    return sum(values[key] for key in keys) // 1024


def measure(output, metadata, command):
    required = {'caller', 'pair_id', 'dataset_id', 'reference', 'settings',
                'allocated_cpus', 'requested_memory_gb'}
    if required - set(metadata):
        raise ValueError('Missing benchmark metadata: ' + ', '.join(sorted(required - set(metadata))))
    if not command:
        raise ValueError('Missing benchmark command')
    if not isinstance(metadata['settings'], dict):
        raise ValueError('Benchmark settings must be a JSON object')
    output.parent.mkdir(parents=True, exist_ok=True)
    slurm_id = os.environ.get('SLURM_JOB_ID')
    warnings = []
    def warn(message):
        warnings.append(message)
        print('BENCHMARK WARNING: ' + message, file=sys.stderr, flush=True)
    cgroup_method, cgroup_path = 'process_tree_sampled_rss_sum', None
    baseline_resident = 0
    if slurm_id:
        try:
            cgroup_method, cgroup_path = memory_cgroup()
            if ('job_' + slurm_id) not in cgroup_path.parts:
                raise ValueError('Memory cgroup is not scoped to this Slurm job')
            baseline_resident = cgroup_resident_kib(cgroup_method, cgroup_path)
        except (OSError, ValueError) as exc:
            warn('Cgroup unavailable; using isolated process-tree RSS: ' + str(exc))
            cgroup_method, cgroup_path = 'process_tree_sampled_rss_sum', None
    started_utc = utc_now()
    started = time.monotonic()
    child = subprocess.Popen(command, start_new_session=True)

    def forward(signum, _frame):
        if child.poll() is None:
            try:
                os.killpg(child.pid, signum)
            except ProcessLookupError:
                pass

    signal.signal(signal.SIGTERM, forward)
    signal.signal(signal.SIGINT, forward)
    peak_pss = peak_rss = peak_resident = samples = unavailable_samples = 0
    unavailable_examples = []
    peak_elapsed = 0.0
    peak_processes = 0
    try:
        while True:
            try:
                pss, rss, count, unavailable = sample(child.pid)
            except OSError as exc:
                pss, rss, count, unavailable = 0, 0, 0, [type(exc).__name__]
                if not unavailable_samples:
                    warn('Process memory unavailable; caller will continue: ' + str(exc))
            if cgroup_path:
                try:
                    resident = cgroup_resident_kib(cgroup_method, cgroup_path)
                except (OSError, ValueError) as exc:
                    warn('Cgroup sampling failed; using process-tree RSS for the entire run: ' + str(exc))
                    cgroup_method, cgroup_path = 'process_tree_sampled_rss_sum', None
                    peak_resident = peak_rss
                    resident = rss
            else:
                resident = rss
            peak_resident = max(peak_resident, resident)
            samples += 1
            if unavailable:
                unavailable_samples += 1
                unavailable_examples.extend(unavailable[:3])
            if pss > peak_pss and not unavailable:
                peak_pss, peak_elapsed, peak_processes = pss, time.monotonic() - started, count
            peak_rss = max(peak_rss, rss)
            try:
                child.wait(timeout=INTERVAL_SECONDS)
                break
            except subprocess.TimeoutExpired:
                pass
    finally:
        code = child.wait()
        report = dict(metadata)
        report.update({
            'schema': SCHEMA,
            'boundary': BOUNDARY,
            'start_utc': started_utc,
            'end_utc': utc_now(),
            'elapsed_seconds': time.monotonic() - started,
            'exit_code': code,
            'slurm_job_id': slurm_id,
            'hostname': os.uname().nodename,
            'sampling_interval_seconds': INTERVAL_SECONDS,
            'memory_method': ('sampled job-cgroup resident memory, excluding unmapped page cache' if cgroup_path
                              else 'sampled command-session process-tree RSS sum; shared pages may be double counted'),
            'memory_method_id': cgroup_method,
            'memory_cgroup': str(cgroup_path) if cgroup_path else None,
            'baseline_cgroup_resident_kib': baseline_resident,
            'peak_sampled_resident_kib': peak_resident if peak_resident else None,
            'supplementary_process_pss_method': 'sum of /proc/PID/smaps_rollup Pss for isolated command session',
            'peak_sampled_pss_kib': peak_pss if peak_pss and not unavailable_samples else None,
            'peak_sampled_rss_sum_kib': peak_rss if peak_rss else None,
            'peak_pss_elapsed_seconds': peak_elapsed if peak_pss else None,
            'processes_at_pss_peak': peak_processes,
            'samples': samples,
            'samples_with_unavailable_process_memory': unavailable_samples,
            'unavailable_process_examples': unavailable_examples[:12],
            'memory_status': ('unavailable' if not peak_resident else
                              'partial' if not cgroup_path and unavailable_samples else 'complete'),
            'measurement_warnings': warnings,
            'limitations': 'One-second sampling may miss brief peaks. Job-cgroup resident includes small recorder overhead; process-tree RSS may double-count shared pages. PSS can be incomplete when setuid helpers block smaps access. Never pool different memory_method_id values or substitute page-cache-inclusive Slurm MaxRSS.',
            'command': command,
            'recorder_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'settings_sha256': hashlib.sha256(json.dumps(metadata.get('settings'),
                                            sort_keys=True, separators=(',', ':')).encode()).hexdigest(),
        })
        temporary = output.with_name(output.name + '.tmp')
        temporary.write_text(json.dumps(report, indent=2) + '\n')
        os.replace(temporary, output)
    return code


def main():
    if len(sys.argv) < 5 or sys.argv[3] != '--':
        raise SystemExit('Usage: measure_benchmark.py OUTPUT METADATA.json -- COMMAND [ARG ...]')
    metadata = json.loads(Path(sys.argv[2]).read_text())
    raise SystemExit(measure(Path(sys.argv[1]), metadata, sys.argv[4:]))


if __name__ == '__main__':
    main()
