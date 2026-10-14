# Development and publication

## Structure

```text
main.nf + modules/        # paired calling and research extensions
sv_signatures/main.nf    # independently selectable signature workflow
sv_signatures/pcawg/     # ClusterSV, timing reference, advanced feature adapters
config/                  # versioned scientific defaults and site templates
scripts/workflow.py      # shared strict-YAML launcher and provenance
tests/ + examples/       # synthetic installation/regression tests
docs/                    # GitHub Pages documentation source
```

Add a caller as a module/adapter with typed option defaults, exact tool identity, task benchmark wrapping, output QC and control tests. Reuse the paired schema and inventory; never add patient IDs or historical job paths to a generic resolver. Separate scientific changes from scheduler changes.

## Releasing a version

Run offline tests, Nextflow smoke, signature numerical controls where relevant, and a real-pair equivalence pilot. Review resolved commands and normalized VCF differences. Tag a version only after documenting these validation results. OCI digests, upstream Git revisions and package/R lockfiles record software identities; lockfiles do not by themselves guarantee identical compiled libraries or biology.

Before public upload, run `python scripts/audit_publication.py`; review every tracked file. Keep BAMs/VCFs from patients, samplesheets, clinical metadata, solver installers/licenses, private paths, credentials, model weights and runtime logs outside Git. Synthetic VCFs under `examples/` are the only genomic examples included.

The independent workflow code is MIT-licensed with owner approval. Upstream tools retain their own licenses; source/resources are fetched separately. See `THIRD_PARTY.md`.

## GitHub and GitHub Pages

Create an empty `jierui-cell/sv-workflows` repository and grant the connected GitHub app access. Upload the reviewed source only. For Pages, select **Settings → Pages → Source: GitHub Actions**, then run the manually triggered `pages` workflow. It builds the Markdown docs and publishes the site. See [GitHub's official Pages workflow guide](https://docs.github.com/en/pages/getting-started-with-github-pages/using-custom-workflows-with-github-pages).

CI runs offline tests on pushes/pull requests; optional real-call and numerical validation should run on properly provisioned compute infrastructure, never untrusted PRs with protected inputs or solver credentials. No constant monitoring is installed.
