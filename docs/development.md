# Development and Publication

## Structure

~~~text
main.nf + modules/       # paired calling
sv_signatures/main.nf   # active COSMIC SV32 analysis
sv_signatures/pcawg/    # retained PCAWG/ClusterSV work, currently on hold
config/                 # defaults and site templates
scripts/workflow.py     # shared YAML launcher and provenance
tests/ + examples/      # synthetic regression controls
docs/                   # GitHub Pages sources
~~~

Add a caller through a module/adapter with typed defaults, pinned tool identity, benchmark wrapping and output QC. Reuse the paired schema and launcher across cohorts. Keep scientific settings separate from scheduler settings.

## Test an Update

Run `python scripts/test_installation.py` and the actual caller entry-point stub/resume test. CI uses synthetic inputs and verifies that the paused PCAWG/ClusterSV route is not launched. COSMIC's numerical clustering control runs within its signature environment. Real-pair comparisons and biological benchmarks provide additional release evidence; record their scope with the version.

OCI digests, source revisions, package locks and resolved commands make changes traceable. Give scientific experiments separate output directories.

## Publish Source and Documentation

Run `python scripts/audit_publication.py` before publication. Public content includes code, guides and synthetic fixtures; datasets, metadata, credentials, binaries and runtime files stay in private storage.

Independent code is MIT-licensed; [upstream notices](../THIRD_PARTY.md) preserve tool licenses.

For GitHub Pages, select **Settings → Pages → Source: GitHub Actions**, then run the manually triggered `pages` workflow. It builds this Markdown documentation as a navigable website. [GitHub Pages setup](https://docs.github.com/en/pages/getting-started-with-github-pages/using-custom-workflows-with-github-pages).

Pull requests run synthetic CI without private cohort inputs or credentials.
