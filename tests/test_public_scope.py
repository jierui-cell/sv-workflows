"""Keep the public guides and configuration focused on the workflows in use."""
from pathlib import Path
import re
import sys
import unittest
import yaml

ROOT = Path(__file__).resolve().parents[1]
CALLERS = {'manta', 'delly', 'gridss', 'svaba', 'brass'}

class PublicScopeTests(unittest.TestCase):
    def test_public_documents_do_not_advertise_unused_workflows(self):
        sources = [ROOT/'README.md', ROOT/'CONFIG_PARAMETERS.md', ROOT/'THIRD_PARTY.md']
        sources += list((ROOT/'docs').glob('*.md'))
        sources += list((ROOT/'examples').rglob('*.md'))
        forbidden = re.compile(r'\b(?:cue|jabba|lumpy|svtyper|nmf|cplex|gurobi|csvc|panprostate)\b', re.I)
        for path in sources:
            self.assertIsNone(forbidden.search(path.read_text()), str(path.relative_to(ROOT)))

    def test_public_configuration_contains_only_active_caller_settings(self):
        config = yaml.safe_load((ROOT/'config/defaults.yaml').read_text())
        maps = {key for key, value in config.items() if isinstance(value, dict)}
        self.assertEqual(maps, {'delly','svaba','gridss','brass','manta_options'})
        with (ROOT/'config/containers.tsv').open() as handle:
            tools = {line.split('\t')[0] for line in handle.read().splitlines()[1:]}
        self.assertEqual(tools, CALLERS | {'bwa'})

    def test_performance_export_defaults_match_active_callers(self):
        sys.path.insert(0, str(ROOT/'scripts'))
        from summarize_benchmarks import CALLERS as defaults
        self.assertEqual(set(defaults), CALLERS)

if __name__ == '__main__': unittest.main()
