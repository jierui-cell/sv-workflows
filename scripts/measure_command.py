#!/usr/bin/env python3
"""GNU time if available; otherwise Linux wait4 command resource accounting."""
import json
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
            os.execv(timer, [timer, '-v', '-a', '-o', output] + command)
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
    raise SystemExit(code)


if __name__ == '__main__':
    main()
