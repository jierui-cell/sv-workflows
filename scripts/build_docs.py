#!/usr/bin/env python3
"""Build the Markdown manual as a navigable static site; never read runtime inputs."""
import html
from html.parser import HTMLParser
from pathlib import Path
import re
from urllib.parse import urlsplit
import markdown

ROOT = Path(__file__).resolve().parents[1]
ORDER = [('index', 'Home'), ('quickstart', 'Quick start'), ('configuration', 'Parameters'),
         ('setup', 'Installation'), ('hpc', 'HPC settings'), ('inputs-outputs', 'Inputs and results'),
         ('signatures', 'SV signatures'), ('validation', 'Validation'),
         ('workflow-status', 'Readiness'), ('performance', 'Performance'), ('development', 'Development')]
STYLE = """
*{box-sizing:border-box}body{margin:0;color:#1f2328;background:#fff;font:16px/1.65 system-ui,sans-serif}
header{border-bottom:1px solid #d1d9e0;background:#f6f8fa;padding:1rem 2rem}
header a{color:#1f2328;font-weight:650;text-decoration:none}
.layout{max-width:1320px;margin:auto;display:grid;grid-template-columns:210px minmax(0,1fr)}
nav{padding:1.5rem 1rem;position:sticky;top:0;align-self:start}
nav a{display:block;padding:.4rem .8rem;border-radius:6px;text-decoration:none;color:#59636e}
nav a[aria-current=page]{background:#eaf3ff;color:#0969da;font-weight:600}
main{min-width:0;max-width:1000px;padding:1.5rem 2.5rem 3rem;border-left:1px solid #d1d9e0}
a{color:#0969da}a:hover{text-decoration:underline}a:focus-visible{outline:2px solid #0969da}
pre{background:#f6f8fa;border-radius:6px;padding:1rem;overflow:auto;font-size:14px}
code{font-size:.9em}table{border-collapse:collapse;width:100%}
.table-scroll{overflow-x:auto;max-width:100%;margin:1rem 0}.table-scroll:focus-visible{outline:2px solid #0969da;outline-offset:2px}
td,th{padding:.55rem .75rem;border:1px solid #d1d9e0;text-align:left}
th{background:#f6f8fa}tbody tr:nth-child(even){background:#fafbfc}h1,h2,h3{line-height:1.3}h1,h2{padding-bottom:.4rem;border-bottom:1px solid #d1d9e0}
main>p,main>ul,main>ol{max-width:80ch}
.page-toc{max-width:80ch;border:1px solid #d1d9e0;border-radius:6px;padding:.6rem 1rem;margin:1rem 0;background:#f6f8fa}
.page-toc summary{cursor:pointer;font-weight:600}.page-toc ul{padding-left:1.2rem;margin:.4rem 0;line-height:1.5}
footer{margin-top:3rem;border-top:1px solid #d1d9e0;padding-top:1rem;color:#59636e;font-size:.85rem}
.skip{position:absolute;left:-10000px}.skip:focus{left:1rem;top:1rem;background:white;padding:.5rem}
img{max-width:100%}@media(max-width:800px){.layout{display:block}nav{position:static;display:flex;flex-wrap:wrap;padding:.75rem}
nav a{padding:.25rem .6rem}main{border-left:0;padding:1rem}header{padding:1rem}}
"""

def rewrite_link(target, source, root, destinations):
    parts = urlsplit(target)
    if parts.scheme or parts.netloc or not parts.path:
        return target
    resolved = (source.parent / parts.path).resolve()
    if resolved in destinations:
        return destinations[resolved] + ('#' + parts.fragment if parts.fragment else '')
    try:
        rel = resolved.relative_to(root)
    except ValueError:
        raise ValueError('Documentation link points outside the source repository: ' + target)
    if not resolved.exists():
        raise ValueError('Broken local documentation link: ' + str(rel))
    return 'https://github.com/jierui-cell/sv-workflows/blob/main/' + rel.as_posix() + (
        '#' + parts.fragment if parts.fragment else '')

class PageLinks(HTMLParser):
    def __init__(self):
        super().__init__()
        self.ids, self.links = set(), []

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if 'id' in attrs:
            self.ids.add(attrs['id'])
        if tag == 'a' and 'href' in attrs:
            self.links.append(attrs['href'])

def verify_generated_links(out):
    """Check local pages and heading anchors without making network requests."""
    pages = {}
    for path in out.glob('*.html'):
        parser = PageLinks()
        parser.feed(path.read_text())
        pages[path.resolve()] = parser
    for source, parser in pages.items():
        for link in parser.links:
            parts = urlsplit(link)
            if parts.scheme or parts.netloc:
                continue
            target = (source.parent / parts.path).resolve() if parts.path else source
            if target not in pages or (parts.fragment and parts.fragment not in pages[target].ids):
                raise ValueError('Broken generated page/anchor link: ' + source.name + ' -> ' + link)

def build(root=ROOT):
    root = Path(root).resolve()
    out = root / 'site'
    out.mkdir(exist_ok=True)
    sources = {p.stem: p for p in (root / 'docs').glob('*.md')}
    sources.update(quickstart=root / 'README.md', configuration=root / 'CONFIG_PARAMETERS.md')
    destinations = {p.resolve(): slug + '.html' for slug, p in sources.items()}
    order = [(slug, label) for slug, label in ORDER if slug in sources]
    order += [(s, s.replace('-', ' ').title()) for s in sorted(sources) if s not in dict(order)]
    for slug, source in sources.items():
        text = source.read_text()
        renderer = markdown.Markdown(extensions=['fenced_code', 'tables', 'toc'],
                                     extension_configs={'toc': {'toc_depth': '2-3'}})
        body = renderer.convert(text)
        if len(re.findall(r'^## ', text, flags=re.MULTILINE)) >= 3:
            contents = ('<details class="page-toc"><summary>On this page</summary>' +
                        renderer.toc + '</details>')
            body = re.sub(r'</h1>', lambda m: m.group() + contents, body, count=1)
        body = re.sub(r'<table>(.*?)</table>',
                      r'<div class="table-scroll" role="region" aria-label="Scrollable table" tabindex="0"><table>\1</table></div>',
                      body, flags=re.DOTALL)
        body = re.sub(r'href="([^"]+)"', lambda m: 'href="' + html.escape(
            rewrite_link(html.unescape(m[1]), source, root, destinations), quote=True) + '"', body)
        nav = '<nav aria-label="Documentation">' + ''.join(
            '<a href="' + s + '.html"' + (' aria-current="page"' if s == slug else '') +
            '>' + html.escape(label) + '</a>' for s, label in order) + '</nav>'
        title = text.splitlines()[0].lstrip('# ').replace('*', '')
        page = ('<!doctype html><html lang="en"><head><meta charset="utf-8">'
                '<meta name="viewport" content="width=device-width, initial-scale=1"><title>' +
                html.escape(title) + ' — SV Workflows</title><style>' + STYLE + '</style></head><body>'
                '<a class="skip" href="#content">Skip to content</a><header>'
                '<a href="index.html">SV Workflows</a> · paired SV calling and signatures</header>'
                '<div class="layout">' + nav + '<main id="content">' + body +
                '<footer>Research use only · Synthetic examples only · '
                '<a href="https://github.com/jierui-cell/sv-workflows">Source repository</a>'
                '</footer></main></div></body></html>')
        (out / (slug + '.html')).write_text(page)
    (out / '.nojekyll').touch()
    verify_generated_links(out)
    print('Built', len(sources), 'documentation pages in', out)
    return len(sources)

if __name__ == '__main__':
    build()
