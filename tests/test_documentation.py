import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]

@unittest.skipUnless(importlib.util.find_spec('markdown'), 'Optional docs renderer is not installed')
class DocumentationTests(unittest.TestCase):
    def builder(self):
        sys.path.insert(0,str(ROOT/'scripts'))
        import build_docs
        return build_docs

    def test_markdown_and_source_links(self):
        docs=self.builder()
        destinations={(ROOT/'README.md').resolve():'quickstart.html',
                      (ROOT/'CONFIG_PARAMETERS.md').resolve():'configuration.html'}
        source=ROOT/'docs/index.md'
        self.assertEqual(docs.rewrite_link('../README.md#usage-instructions',source,ROOT,destinations),
                         'quickstart.html#usage-instructions')
        self.assertEqual(docs.rewrite_link('../CONFIG_PARAMETERS.md',source,ROOT,destinations),'configuration.html')
        self.assertEqual(docs.rewrite_link('examples/cohort.yaml',ROOT/'README.md',ROOT,destinations),
                         'https://github.com/jierui-cell/sv-workflows/blob/main/examples/cohort.yaml')

    def test_missing_and_external_local_paths_rejected(self):
        docs=self.builder()
        for link in ('absent.md','../../outside.md'):
            with self.assertRaises(ValueError):
                docs.rewrite_link(link,ROOT/'README.md',ROOT,{})

    def test_build_navigation_and_escaping(self):
        docs=self.builder()
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); (root/'docs').mkdir()
            (root/'docs/index.md').write_text('# Home\n\n[Start](../README.md)\n')
            (root/'README.md').write_text('# *Quick start*\n\n[Parameters](CONFIG_PARAMETERS.md)\n')
            (root/'CONFIG_PARAMETERS.md').write_text('# Parameters\n')
            self.assertEqual(docs.build(root),3)
            page=(root/'site/quickstart.html').read_text()
            self.assertIn('href="configuration.html"',page)
            self.assertIn('href="quickstart.html" aria-current="page"',page)
            self.assertIn('<title>Quick start — SV Workflows</title>',page)

    def test_long_page_contents_and_keyboard_table(self):
        docs=self.builder()
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); (root/'docs').mkdir()
            (root/'docs/index.md').write_text('# Home\n')
            (root/'README.md').write_text('# Guide\n\n## Setup\n\n## Run\n\n## Results\n\n| Item | Path |\n| --- | --- |\n| Toy | Example |\n')
            (root/'CONFIG_PARAMETERS.md').write_text('# Parameters\n')
            docs.build(root)
            page=(root/'site/quickstart.html').read_text()
            self.assertIn('<summary>On this page</summary>',page)
            self.assertLess(page.index('<h1'),page.index('<details class="page-toc">'))
            self.assertIn('href="#setup"',page)
            self.assertIn('class="table-scroll" role="region"',page)
            self.assertIn('aria-label="Scrollable table" tabindex="0"',page)

    def test_missing_heading_anchor_is_rejected(self):
        docs=self.builder()
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); (root/'docs').mkdir()
            (root/'docs/index.md').write_text('# Home\n\n[Wrong heading](../README.md#missing)\n')
            (root/'README.md').write_text('# Guide\n')
            (root/'CONFIG_PARAMETERS.md').write_text('# Parameters\n')
            with self.assertRaisesRegex(ValueError,'Broken generated page/anchor'):
                docs.build(root)

if __name__=='__main__': unittest.main()
