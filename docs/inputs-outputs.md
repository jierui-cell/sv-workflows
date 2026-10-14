# Inputs, Results and Provenance

## Paired Samplesheet

Copy [examples/samplesheet.tsv](../examples/samplesheet.tsv) and edit it in private storage. The same schema supports WGS, WES and targeted runs; set assay type in YAML.

All nine TSV headers are retained: `patient_id`, `tumor_bam`, `tumor_bai`, `normal_bam`, `normal_bai`, `sex`, `tumor_phenotype`, `purity`, `ploidy`.

The first five values are required. Metadata can be blank. Use absolute paths, unique IDs and tabs; keep trailing empty metadata fields. Provided purity is >0 and ≤1, and ploidy is >0. Schema and selected IDs are checked before tasks are created.

Missing BRASS purity/ploidy uses `native_defaults`, recorded in effective settings as defaults rather than measured values.

## Signature Manifest

The exact TSV header is `cohort pair_id caller vcf source_status source_sha256`. One row identifies one caller's output for one pair. The manifest builder takes the five-column inventory template, validates identifiers and computes source hashes.

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

ClusterSV exports per-junction cluster/footprint tables and optional timing annotations. The `features` stage adds `ClusterSV38.matrix.tsv`, `feature_events.tsv` and feature QC. The `signatures` stage adds `ClusterSV.NMF.signatures.tsv`, `activities.tsv` and `fractions.tsv` with the same prefix, plus rank/reconstruction QC and the fitted R object.

NMF profiles are features × components; activities are samples × components. Zero-count samples retain zero activity. [Feature definitions](signatures.md#cn-free-features-and-de-novo-nmf).

## Provenance

Every run saves resolved YAML, code/container/reference identities, exact commands, trace and QC. Keep these with the results; huge BAM inputs are referenced by path rather than copied into provenance. Private paths, metadata and runtime outputs stay outside Git.
