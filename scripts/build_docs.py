#!/usr/bin/env python3
"""Build the Markdown manual as a navigable static site; never read runtime inputs."""
import html
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
code{font-size:.9em}table{border-collapse:collapse;width:100%;display:block;overflow:auto}
td,th{padding:.55rem .75rem;border:1px solid #d1d9e0;text-align:left}
th{background:#f6f8fa}h1,h2,h3{line-height:1.3}h1,h2{padding-bottom:.4rem;border-bottom:1px solid #d1d9e0}
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
        body = markdown.markdown(text, extensions=['fenced_code', 'tables', 'toc'])
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
    print('Built', len(sources), 'documentation pages in', out)
    return len(sources)

if __name__ == '__main__':
    build()
