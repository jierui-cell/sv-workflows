# Validation: setup is not accuracy

## Fast offline tests

`python scripts/test_installation.py` runs unit tests for configuration/schema errors, paired-role handling, scientific option defaults, reference-cache reuse, VCF orientations/mate collapse, filtering adapters, BRASS statistics, optimizer convergence guards, and memory/exit-code recording. Fixtures are synthetic and require no EGA files or licensed software.

`python scripts/workflow.py test --profile test` then executes a real Nextflow task on a synthetic VCF. It checks three retained events: deletion, duplication and one interchromosomal adjacency represented by two reciprocal BND records. This test does not run a structural variant caller. On an HPC, launch it within a small compute allocation.

## Container/setup checks

Pull selected digest-pinned images, check the caller executable inside each with `apptainer exec IMAGE command --help`, and validate the exact alignment reference/indexes. BRASS also requires native bundle validation. ClusterSV and SigProfiler have genuine numerical controls that need their pinned R/Python environments; they are not part of the lightweight offline suite.

## Before a new cohort or release

1. Run the unchanged original and consolidated workflow on the same matched pair, comparing normalized biological events, sample roles, filters and resolved commands. Differences require review, not automatic acceptance.
2. Normal-self negative control for each caller, using two read-only aliases with distinct sample labels and **distinct read names** where the caller requires them. A low count is expected but is not a guaranteed zero; investigate every retained PASS candidate. Never interpret an input-label error as successful control validation.
3. Positive benchmark with a suitable truth set and confident regions. GIAB germline SV truth measures discovery, not tumour-specific somatic precision; do not treat a normal-self GIAB run as a positive somatic benchmark. Use validated somatic spike-in/truth or report the scope honestly.
4. Tumour/normal swap and repeated-run checks. Compare event-normalized VCFs, not date headers alone. Document assembly/threading variation.
5. Review emitted commands after parameter/version changes, plus representative calls in IGV. Inspect reference gaps, centromeric noise, adapter artifacts, depth and matched-normal evidence.

Only then expand to three concurrent pairs, then measured cohort concurrency. CI passing does not waive biological validation. Reference/tool/filter changes create a new configuration requiring revalidation.
