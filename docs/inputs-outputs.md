# Inputs, Results and Provenance

## Paired Samplesheet

Copy [examples/samplesheet.tsv](../examples/samplesheet.tsv) to `inputs/samplesheet.tsv`. Place or link BAMs/indexes in `inputs/bams/`, or override their paths. `inputs/` is Git-ignored. The same schema supports WGS, WES and targeted runs; set assay type in YAML.

All nine TSV headers are retained: `patient_id`, `tumor_bam`, `tumor_bai`, `normal_bam`, `normal_bai`, `sex`, `tumor_phenotype`, `purity`, `ploidy`.

The first five values are required. Metadata can be blank. BAM/index paths can be repository-root-relative (for example, `inputs/bams/TOY_PAIR.tumor.bam`) or absolute. Use unique IDs and tabs; keep trailing empty metadata fields. Provided purity is >0 and ≤1, and ploidy is >0. Schema and selected IDs are checked before tasks are created.

File checks run before caller-task submission. The launcher leaves the original sheet untouched, saves `resolved-samplesheet.tsv` with absolute paths in the run logs, and uses a stable hashed copy under `.nextflow/input-manifests/` so repeated launches can resume cached tasks.

Missing BRASS purity/ploidy uses `native_defaults`, recorded in effective settings as defaults rather than measured values.

## Signature Manifest

The exact TSV header is `cohort pair_id caller vcf source_status source_sha256`. One row identifies one caller's output for one pair. The manifest builder takes `inputs/callsets.tsv`, validates identifiers, resolves VCF paths to absolute paths, computes source hashes and writes `manifests/sv_manifest.tsv`.

Use final somatic VCFs for signature input. GRIDSS joint calls remain a separate output before somatic filtering.

## Output Layout

~~~text
results/RUN/
  preflight/PAIR/
  manta/PAIR/
  delly/PAIR/output/
  svaba/PAIR/output/
  gridss/PAIR/output/
  brass/PAIR/output/
  provenance/RUN_NAME/
logs/WORKFLOW-TIMESTAMP/   # settings, trace, driver log and status
work/WORKFLOW/            # task commands and resume cache
~~~

Native caller names/hierarchies are retained. [Quick-start results table](../README.md#retrieve-your-results) identifies the main VCFs.

Signature results use `OUTDIR/COHORT/CALLER/group_output/`. Optional pair selections retain the existing smoke/subset grouping layout.

## Signature Outputs

COSMIC exports `SV32.matrix.tsv` and `SV32.exposures.tsv`.

PCAWG/ClusterSV signature analysis is on hold. Earlier clustering, feature and NMF outputs are retained; no new analysis is being run through that route. [Signature status](signatures.md#pcawg-status).

## Provenance

Every run saves resolved YAML, code/container/reference identities, exact commands, trace and QC. Keep these with the results; huge BAM inputs are referenced by path rather than copied into provenance. Private paths, metadata and runtime outputs stay outside Git.
