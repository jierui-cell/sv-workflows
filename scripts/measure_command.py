#!/usr/bin/env python3
"""GNU time if available; otherwise Linux wait4 command resource accounting."""
import json
import hashlib
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time


def main():
    output, command = sys.argv[1], sys.argv[2:]
    if not command:
        raise SystemExit('Usage: measure_command.py OUTPUT COMMAND [ARG ...]')
    timer = shutil.which('time')
    if timer:
        probe = subprocess.run([timer, '--version'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        if b'GNU' in probe.stdout:
            started = time.monotonic()
            code = subprocess.call([timer, '-v', '-a', '-o', output] + command)
            write_status(output, command, code, time.monotonic()-started, 'GNU time -v')
            raise SystemExit(code)
    started = time.monotonic()
    process = subprocess.Popen(command)
    _, status, usage = os.wait4(process.pid, 0)
    code = os.WEXITSTATUS(status) if os.WIFEXITED(status) else 128 + os.WTERMSIG(status)
    process.returncode = code
    report = {'measurement': 'Linux wait4 rusage; GNU time unavailable', 'argv': command,
        'elapsed_seconds': time.monotonic()-started, 'user_seconds': usage.ru_utime,
        'system_seconds': usage.ru_stime, 'max_rss_kib': usage.ru_maxrss, 'exit_code': code,
        'limitation': 'command/awaited descendants; not summed concurrent process memory'}
    with Path(output).open('a') as handle:
        handle.write(json.dumps(report)+'\n')
    write_status(output, command, code, report['elapsed_seconds'], report['measurement'])
    raise SystemExit(code)


def write_status(output, command, code, elapsed, method):
    timing = Path(output)
    target = Path(str(output) + '.status.json')
    temporary = target.with_name(target.name + '.tmp')
    temporary.write_text(json.dumps({'argv':command, 'exit_code':code,
        'elapsed_seconds':elapsed, 'method':method,
        'timing_sha256':hashlib.sha256(timing.read_bytes()).hexdigest()})+'\n')
    os.replace(temporary, target)


if __name__ == '__main__':
    main()
