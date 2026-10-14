# SV Signature Analysis

Two routes answer different questions: **COSMIC SV32 fits known reference signatures**; **ClusterSV + NMF discovers patterns in your cohort**. Both start from a hashed inventory of caller-specific WGS VCFs and run independently by cohort/caller.

## Supply Your Callsets

Create a five-column inventory from [the template](../examples/callsets.tsv), then:

~~~bash
mkdir -p manifests
python scripts/build_sv_manifest.py --inventory /path/to/private/callsets.tsv \
  --output manifests/sv_manifest.tsv
cp examples/signatures.yaml config/private-signatures.yaml
~~~

Set the manifest, output directory, Python path and account/partition. Use GRCh37 callsets. The pipeline verifies source hashes and VCF primary-contig lengths before conversion.

## COSMIC SV32

Install the locked Python environment once:

~~~bash
python3.11 -m venv sv_signatures/.venv
sv_signatures/.venv/bin/python -m pip install -r sv_signatures/requirements.lock.txt
~~~

Launch from the repository root with your site's driver allocation:

~~~bash
mkdir -p logs
sbatch --account=YOUR_ACCOUNT --partition=YOUR_PARTITION \
  --cpus-per-task=1 --mem=4G --time=24:00:00 --output=logs/cosmic-%j.log \
  scripts/run.slurm.sh signatures --config config/private-signatures.yaml --set method=cosmic_sv32
~~~

The workflow converts retained VCF events to BEDPE, applies SigProfiler clustering, builds 32-channel counts and refits COSMIC v3.6 SV1–SV12. The intrachromosomal minimum is 1 kb. The locked pandas semantics and numerical clustering control run automatically.

Retrieve `SV32.matrix.tsv`, `SV32.exposures.tsv`, conversion audits and reconstruction QC. Set `skip_fit: true` for counts only.

## ClusterSV Setup

The pinned [ClusterSV](https://github.com/cancerit/ClusterSV) implementation performs the original [PCAWG proximity clustering and footprint calculation](https://www.nature.com/articles/s41586-019-1913-9).

For ClusterSV alone, use the lean Python environment below; an existing COSMIC environment can be reused. Fetch source, then restore the R 4.3 environment in a compute allocation:

~~~bash
python3.11 -m venv sv_signatures/.venv
sv_signatures/.venv/bin/python -m pip install -r requirements.txt \
  -r sv_signatures/requirements-clustersv.lock.txt
python sv_signatures/pcawg/prepare_sources.py --clustersv-only
Rscript sv_signatures/pcawg/setup_env.R "$PWD"
~~~

The same R environment contains NMF. You can launch directly from your VCF manifest; timing annotation is optional:

~~~bash
sbatch --account=YOUR_ACCOUNT --partition=YOUR_PARTITION \
  --cpus-per-task=1 --mem=4G --time=24:00:00 --output=logs/clustersv-%j.log \
  scripts/run.slurm.sh signatures --config config/private-signatures.yaml \
  --set method=pcawg_clustersv --set outdir=results/clustersv_run01
~~~

The default `pcawg_stage: signatures` runs the full sequence:

~~~text
VCFs → Oriented BEDPE → ClusterSV → CN-free feature matrix → Seeded KL-NMF
                                 └→ Footprints + optional endpoint timing
~~~

## CN-Free Features and De Novo NMF

No copy-number files or classified-event manifest are needed. Our local **ClusterSV38** scheme assigns each retained junction to exactly one channel:

| Axis | Categories |
| --- | --- |
| Cluster membership | Singleton; clustered (cluster size >1) |
| Intrachromosomal orientation class | DEL; DUP; INV |
| Span | <1 kb; 1–10 kb; 10–100 kb; 100 kb–1 Mb; 1–10 Mb; ≥10 Mb |
| Interchromosomal junction | TRA, with no span bin |

This gives `2 × (3 × 6 + 1) = 38` channels. Bin lower bounds are inclusive and upper bounds exclusive. Reciprocal BND mates are counted once. Endpoint replication timing is exported separately, not included in these NMF channels.

The 38-channel scheme is a **CN-free extension**, not the published PCAWG 44-feature matrix. Its orientation labels describe individual junctions, not reconstructed complex-event mechanisms.

NMF factorizes feature counts into signature profiles × sample activities. It uses repeated seeded KL-NMF fits with per-iteration factor flooring for sparse counts. Default settings are 100 fits, seed 1234 and ranks 2–10, limited by nonzero dimensions. Automatic selection takes the rank before the steepest adjacent cophenetic decline; if there is no decline, it chooses the smallest rank at maximal finite cophenetic correlation, or the smallest admissible rank if correlations are undefined. The rank survey and selection rule are saved.

~~~yaml
method: pcawg_clustersv
pcawg_stage: signatures  # clustering | features | signatures
pcawg_rank: 0            # Automatic; set >=2 for a fixed rank
pcawg_rank_min: 2
pcawg_rank_max: 10
pcawg_nrun: 100
pcawg_seed: 1234
~~~

At least two nonzero samples and two nonzero channels are needed. Zero-count samples are retained with zero activities. New signatures are labeled `cSV1…` within each cohort/caller group; labels are not matched across independent runs. Compare profiles before comparing activities across groups.

Main files:

| File | Contents |
| --- | --- |
| `ClusterSV38.matrix.tsv` | Feature × sample integer counts |
| `feature_events.tsv` | Junction-to-channel assignments |
| `ClusterSV.NMF.signatures.tsv` | 38-channel component profiles; each sums to one |
| `ClusterSV.NMF.activities.tsv` | Sample × component fitted activities |
| `ClusterSV.NMF.fractions.tsv` | Within-sample normalized activities |
| `NMF.rank_survey.tsv` | Automatic rank survey; absent for fixed rank |
| `NMF.reconstruction_qc.tsv` | Cosine similarity and residuals per sample |
| `signature_status.json` | Rank, seed, solver and counting profile |

Activities are fitted contributions, not necessarily exact observed SV totals. Inspect reconstruction and rank stability when interpreting cohort patterns.

## ClusterSV Annotations and Policies

Readable `.clusters.tsv` files contain one row per junction with cluster ID/size, footprint IDs/bounds and `last_agglomeration_value`. A cluster is a group of statistically nearby junctions; a footprint describes the genomic intervals it covers.

The last-agglomeration and `.sv_distance_pvals` values are upstream comparison-count-scaled statistics, not ordinary 0–1 p-values or somatic confidence scores.

The default `local_cen_telo_exclusion` additionally removes endpoints in the upstream boundary file's centromere/telomere intervals. This is our filtering policy; upstream ClusterSV uses those coordinates for footprints. The chromosome 1 centromere interval is broad (~28.7 Mb). Select `primary_only` to retain those endpoints on autosomes+X. When timing is enabled, the annotator reports any endpoint without a measured reference value.

The historical method key `pcawg_clustersv` is retained for compatibility. De novo NMF is a local downstream analysis, not reconstruction of the published PCAWG 16 reference signatures. High-junction-count samples can need more memory because ClusterSV builds quadratic distance matrices.

## Replication Timing

This optional annotation is independent of the CN-free features/NMF. Leave `pcawg_timing_source: none` (default) for direct VCF-to-signature analysis. Set it to `reference` after preparing the reusable reference below.

The reference averages public ENCODE/UW GRCh37 NHEK, GM12878 and IMR90 tracks in 1-kb bins, then fixes genome-wide tertiles: low scores late, middle mid, high early. Tertiles are our category policy.

~~~bash
python3.11 -m venv environments/replication-timing-python
environments/replication-timing-python/bin/python -m pip install \
  -r sv_signatures/pcawg/replication_requirements.lock.txt
cp examples/timing.yaml config/private-timing.yaml
~~~

Set account/partition and output location, then build on compute nodes:

~~~bash
sbatch --account=YOUR_ACCOUNT --partition=YOUR_PARTITION \
  --cpus-per-task=1 --mem=4G --time=4:00:00 --output=logs/timing-%j.log \
  scripts/run.slurm.sh timing --config config/private-timing.yaml
~~~

Set `pcawg_timing_reference` to the verified directory and `pcawg_timing_source: reference` in signature YAML. Reference construction records input/output hashes and verification; reuse it across GRCh37 cohorts without recalculating cutoffs per sample.
