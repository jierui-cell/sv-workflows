#!/usr/bin/env python3
"""Build a static GitHub Pages manual, without any patient/runtime inputs."""
import html
from pathlib import Path
import re
import markdown
ROOT=Path(__file__).resolve().parents[1]
out=ROOT/'site'; out.mkdir(exist_ok=True)
pages=sorted((ROOT/'docs').glob('*.md'))
nav='<nav>'+''.join('<a href="'+p.stem+'.html">'+html.escape(p.stem.replace('-',' ').title())+'</a>' for p in pages)+'</nav>'
style='body{margin:0;color:#173042;background:#f4f7f9;font:17px/1.65 system-ui,sans-serif}header{background:#173042;color:white;padding:1.5rem}main{max-width:980px;margin:2rem auto;padding:0 1.5rem}nav{display:flex;gap:1rem;flex-wrap:wrap;margin:1rem 0}a{color:#006e82}pre{background:#e6edf1;padding:1rem;overflow:auto}code{font-size:.9em}table{border-collapse:collapse;width:100%;display:block;overflow:auto}td,th{padding:.6rem;border-bottom:1px solid #bbc9d2;text-align:left}h1,h2,h3{line-height:1.25}footer{margin-top:3rem;border-top:1px solid #bbc9d2;padding:1rem 0;font-size:.9rem}'
for page in pages:
    body=markdown.markdown(page.read_text(),extensions=['fenced_code','tables','toc'])
    body=re.sub(r'href="([A-Za-z0-9_-]+)\.md([#][^"]*)?"',lambda m:'href="'+m[1]+'.html'+(m[2] or '')+'"',body)
    title=page.read_text().splitlines()[0].lstrip('# ').strip()
    content='<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>'+html.escape(title)+' — SV Workflows</title><style>'+style+'</style></head><body><header>SV Workflows · reproducible HPC research</header><main>'+nav+body+'<footer>Research use only. Public examples are synthetic. <a href="https://github.com/jierui-cell/sv-workflows">Source and releases</a></footer></main></body></html>'
    (out/(page.stem+'.html')).write_text(content)
(out/'.nojekyll').touch()
print('Built',len(pages),'documentation pages in',out)
