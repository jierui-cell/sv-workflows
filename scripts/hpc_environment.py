"""Load configured Environment Modules before launching host tools."""
import os
import re
import subprocess
import sys


def module_names(value):
    if value is None:
        return []
    if not isinstance(value,list) or any(not isinstance(name,str) or
            not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_./+@-]*',name) for name in value):
        raise ValueError('modules must be a YAML list of module names, for example [java/20.0.1, samtools/1.20]')
    return list(dict.fromkeys(value))


def load_modules(value):
    names = module_names(value)
    if not names:
        return {'modules': []}
    script = r'''
set -e
if ! type module >/dev/null 2>&1; then
    for init in /etc/profile.d/lmod.sh /etc/profile.d/modules.sh /usr/share/Modules/init/bash /usr/share/lmod/lmod/init/bash; do
        if [[ -r "$init" ]]; then
            source "$init" >&2
            if type module >/dev/null 2>&1; then break; fi
        fi
    done
fi
if ! type module >/dev/null 2>&1; then
    echo 'Environment Modules is unavailable; use modules: [] for installed tools.' >&2
    exit 2
fi
module load "$@" >&2
env -0
'''
    result = subprocess.run(['bash','-c',script,'sv-workflows-modules',*names],
                            stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    if result.returncode:
        raise ValueError('Module setup failed for '+', '.join(names)+':\n'+os.fsdecode(result.stderr).strip())
    environment = {}
    for entry in result.stdout.split(b'\0'):
        if not entry:
            continue
        key,sep,value = entry.partition(b'=')
        if not sep:
            raise ValueError('Unexpected output from module initialization')
        environment[os.fsdecode(key)] = os.fsdecode(value)
    if not environment.get('PATH'):
        raise ValueError('Module initialization returned no PATH')
    os.environ.clear()
    os.environ.update(environment)
    print('Loaded modules: '+', '.join(names),file=sys.stderr,flush=True)
    return {'modules':names,'loaded_modules':environment.get('LOADEDMODULES','').split(':') if environment.get('LOADEDMODULES') else []}


def main():
    import argparse
    from pathlib import Path
    import yaml
    from workflow import resolve_config
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description='Run a setup command using the modules in your caller YAML.')
    parser.add_argument('--config',type=Path,action='append',default=[])
    parser.add_argument('command',nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command[:1] == ['--'] else args.command
    if not command:
        parser.error('Supply a command after --, for example -- samtools --version')
    default = root/'config/private-cohort.yaml'
    paths = args.config or ([default] if default.is_file() else [])
    try:
        config = resolve_config('callers',paths,root=root)
        load_modules(config.get('modules'))
        raise SystemExit(subprocess.run(command).returncode)
    except (ValueError,OSError,yaml.YAMLError) as exc:
        parser.error(str(exc))


if __name__ == '__main__':
    main()
