# Selectable SV signature analysis

Both methods consume caller-specific final VCFs through the same hashed manifest. BND mates are collapsed to one adjacency; caller-specific filter policies and excluded records are audited. These are five separate callsets, not five independent votes or a consensus truth set.

## COSMIC SV32

Uses SigProfilerMatrixGenerator 1.3.6 and the tested pandas 2.2.3 semantics. The current fitting adapter explicitly uses **GRCh37**; do not supply another assembly without an implementation change and revalidation.

Install the locked environment separately. Run the numerical control only in an allocated compute session (not on an HPC login node):

```bash
python3.11 -m venv sv_signatures/.venv
sv_signatures/.venv/bin/python -m pip install -r sv_signatures/requirements.lock.txt
# Numerical computation: allocated compute node only.
sv_signatures/.venv/bin/python sv_signatures/test_clustering.py
```

Run the numerical control on an allocated node. It checks actual clustered/nonclustered assignment; configuration parsing alone is insufficient. The compatibility layer fails closed if tested versions/semantics change. The default minimum intrachromosomal SV size is 1 kb. The 32-feature matrix is refitted to COSMIC v3.6 SV1–SV12, not de-novo extraction. Check reconstruction errors and exposure normalization before interpretation. Very small/empty callsets can be uninformative.

Prepare a private `config/private-signatures.yaml` containing the manifest, output directory, environment path and scheduler settings. `python` should be an absolute environment executable path.

```bash
mkdir -p logs
sbatch --account=YOUR_ACCOUNT --partition=YOUR_PARTITION \
  --cpus-per-task=1 --mem=4G --time=24:00:00 --output=logs/cosmic-%j.log \
  scripts/run.slurm.sh signatures --config config/private-signatures.yaml --set method=cosmic_sv32
```

## Original ClusterSV, without copy number

The pinned [ClusterSV](https://github.com/cancerit/ClusterSV) method reproduces the proximity clustering/footprint part of the [PCAWG structural-variation study](https://www.nature.com/articles/s41586-019-1913-9). It is not the newer cSVc classifier.

```bash
python sv_signatures/pcawg/prepare_sources.py
# In an allocated compute session with R 4.3 and build dependencies:
Rscript sv_signatures/pcawg/setup_env.R "$PWD"
python sv_signatures/pcawg/test_clustersv.py --pipeline-root "$PWD" --clustering-only
```

**Before submitting ClusterSV analysis, complete [replication-timing reference construction](#replication-timing) below.** The default workflow requires both `reference.json` and `verification.json`. If you choose a nondefault timing output location, set `pcawg_timing_reference` to that verified directory in your signature YAML.

```bash
# Submit a separate driver for analysis; use a fresh output directory.
mkdir -p logs
sbatch --account=YOUR_ACCOUNT --partition=YOUR_PARTITION \
  --cpus-per-task=1 --mem=4G --time=24:00:00 --output=logs/clustersv-%j.log \
  scripts/run.slurm.sh signatures --config config/private-signatures.yaml \
  --set method=pcawg_clustersv --set pcawg_stage=clustering --set outdir=results/clustersv_run01
```

Uses GRCh37 autosomes+X, upstream hg19 centromere/telomere masks, a 1-bp junction minimum and explicit orientation. Ambiguous inversion orientation is excluded rather than invented. Raw distance matrices are quadratic in junction count; high-burden samples need more memory/time. Inspect conversion/mask audits to understand differences from SV32 counts.

Without copy-number information, clusters do **not** establish balanced rearrangements, templated insertions, chromoplexy, chromothripsis or foldback mechanisms. The CN-free release provides clusters/footprints/timing; a supported observable-feature summary is a next-stage extension, not the complete original 44-feature schema. The original PCAWG 16-signature numerical reference is not bundled and reference decomposition is deferred.

## Replication timing

Public ENCODE/UW GRCh37 NHEK, GM12878 and IMR90 tracks are averaged in 1-kb bins. A separate reference-construction workflow creates fixed genome-wide tertiles: lower scores late, middle mid, higher early. Tertiles are a chosen policy, **not verified original-paper early/mid/late thresholds**. Unmeasured endpoint bins are audited and cause annotation failure rather than imputation.

```bash
python3.11 -m venv environments/replication-timing-python
environments/replication-timing-python/bin/python -m pip install \
  -r sv_signatures/pcawg/replication_requirements.lock.txt
cp examples/timing.yaml config/private-timing.yaml
```

Edit the copied timing YAML for your account/partition and output location. ClusterSV's pinned source must already be prepared, because it supplies chromosome sizes. Submit the reference construction to compute nodes:

```bash
mkdir -p logs
sbatch --account=YOUR_ACCOUNT --partition=YOUR_PARTITION \
  --cpus-per-task=1 --mem=4G --time=4:00:00 --output=logs/timing-%j.log \
  scripts/run.slurm.sh timing --config config/private-timing.yaml
```

Source URLs/checksums and output hashes are recorded. A verified reference may be reused across cohorts aligned to this assembly; do not recalculate category cutoffs per tumour. Other assemblies require an explicitly validated reference.

## Advanced 44-feature and optional NMF paths

Retained separately for compatibility with supplied, validated mechanistic classifications and annotations. `pcawg_stage=features` constructs the pinned **panprostate_2026_core** 44-channel profile, not a claim of complete CN-free PCAWG reproduction. Several features overlap or count clusters rather than junctions; totals need not equal SV count. This stage rejects missing classifications and annotation mismatches. `pcawg_stage=signatures` explicitly requests de-novo NMF, not fitting to the original PCAWG 16 signatures. Neither is automatically run in the CN-free workflow.
