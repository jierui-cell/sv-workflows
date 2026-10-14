# SV Signature Analysis

**COSMIC SV32 is the active signature-analysis route.** It builds a 32-channel matrix from caller-specific WGS VCFs and fits COSMIC SV1–SV12 reference signatures.

**PCAWG/ClusterSV signature analysis is on hold.** Clustering, feature construction, de novo NMF and related analysis are not being run. The implementation and existing outputs are retained for a later review. Selecting `method=pcawg_clustersv` reports the hold before submitting analysis tasks.

## Supply Your Callsets

Copy the five-column [template](../examples/callsets.tsv) to `inputs/callsets.tsv` and edit its rows to reference your caller outputs, then:

~~~bash
mkdir -p inputs manifests
python scripts/build_sv_manifest.py --inventory inputs/callsets.tsv \
  --output manifests/sv_manifest.tsv
cp examples/signatures.yaml config/private-signatures.yaml
~~~

The manifest (`manifests/sv_manifest.tsv`) and Python environment (`sv_signatures/.venv/bin/python`) already match the defaults. Set your account/partition and output directory. Use GRCh37 WGS callsets. The pipeline verifies source hashes and VCF primary-contig lengths before conversion. Groups run independently by cohort and caller.

## Install Once

~~~bash
python3.11 -m venv sv_signatures/.venv
sv_signatures/.venv/bin/python -m pip install -r sv_signatures/requirements.lock.txt
~~~

COSMIC does not require R, ClusterSV, copy-number inputs or replication-timing reference construction.

You can verify input files and the environment executable with `python scripts/workflow.py signatures --config config/private-signatures.yaml --check-files`. The same check runs automatically at launch.

## Launch COSMIC SV32

From the repository root with your site's driver allocation:

~~~bash
mkdir -p logs
sbatch --account=YOUR_ACCOUNT --partition=YOUR_PARTITION \
  --cpus-per-task=1 --mem=4G --time=24:00:00 --output=logs/cosmic-%j.log \
  scripts/run.slurm.sh signatures --config config/private-signatures.yaml --set method=cosmic_sv32
~~~

~~~text
Somatic VCFs → BEDPE → SV clustering → 32-channel counts → COSMIC SV1–SV12 fit
~~~

The workflow converts retained VCF events to BEDPE, applies SigProfiler clustering, builds SV32 counts and refits COSMIC v3.6 reference signatures. The intrachromosomal minimum is 1 kb. The locked pandas semantics and numerical clustering control are checked automatically.

## Retrieve Results

Results are grouped under `OUTDIR/COHORT/CALLER/group_output/`:

| File | Contents |
| --- | --- |
| `SV32.matrix.tsv` | 32-channel feature counts per sample |
| `SV32.exposures.tsv` | Fitted reference-signature exposures |
| Conversion/source audits | Retained events, source hashes and input decisions |
| Reconstruction QC | Fit quality and residuals |

Set `skip_fit: true` for counts only. Scientific/resource settings and overrides are listed in [the parameter guide](../CONFIG_PARAMETERS.md#signature-settings).

## PCAWG Status

The `pcawg_clustersv` route, its CN-free feature/NMF extension and associated numerical workflows are paused. They are excluded from active launch instructions and CI execution. Source, settings, tests and earlier outputs remain available; no files are deleted. We will revisit the method design before enabling it again.
