# Inputs, outputs and provenance

## Paired WGS samplesheet

Copy the already tab-separated [samplesheet template](../examples/samplesheet.tsv) to an existing private folder outside the checkout, then edit it:

```bash
cp examples/samplesheet.tsv /path/to/private/paired_samples.tsv
```

Keep all nine headers and the tabs between columns. The required headers are `patient_id`, `tumor_bam`, `tumor_bai`, `normal_bam`, `normal_bai`, `sex`, `tumor_phenotype`, `purity` and `ploidy`. Do not remove trailing empty metadata fields. The parser accepts unquoted TSV, not space-delimited text or CSV.

The first five values must be present. Metadata values may be blank and generate warnings; no purity/ploidy is inferred from missingness. Purity must be >0 and ≤1 when provided; ploidy must be >0. Headers, row widths, duplicate/unsafe pair IDs, same-file pairings, selected IDs and metadata ranges are validated before creating task channels. BAM biological identity still needs separate confirmation. Use absolute paths.

BRASS's missing-statistics policy is `native_defaults`: absent values use native baseline values and are labeled as defaults, not measurements. Mixed measured/default cases use a narrowly tested adapter; retain `effective-settings.json`. These defaults do not validate BRASS's copy-number interpretation.

## Signature manifest

The exact header is `cohort pair_id caller vcf source_status source_sha256` (tab-separated). One row per cohort/pair/caller. `scripts/build_sv_manifest.py` takes a five-column private callset inventory, validates IDs/duplicates, and calculates hashes. It does not certify somatic validity. The manifest references final caller-specific VCFs, not a concatenation of all VCFs.

## Output layout

```text
results/RUN/
  preflight/PAIR/
  manta/PAIR/
  delly/PAIR/output/
  svaba/PAIR/output/
  gridss/PAIR/output/       # raw joint + somatic subsets kept separate
  brass/PAIR/output/       # callset, statistics, stage logs, QC
  provenance/RUN_NAME/
logs/WORKFLOW-TIMESTAMP/   # resolved YAML, invocation/hashes, trace, driver status
work/WORKFLOW/            # exact task commands and resume cache
```

Manta's native result hierarchy and other callers' native files are deliberately retained. Output roles are consistent even where native filename conventions differ. The public release does not flatten files that may share names.

Signature outputs group by cohort/caller and include conversion audits, source hashes, per-sample totals, clustering controls, feature matrices and verification reports. Failed or research-only outputs remain explicitly labeled. `SUCCESS` means technical completion under the stated policy, not clinical validity.

Keep the resolved YAML, code/container/reference hashes, task command logs, trace and key QC. Do not duplicate huge inputs into provenance. Runtime logs and manifests contain private paths/IDs and must never be pushed to GitHub.
