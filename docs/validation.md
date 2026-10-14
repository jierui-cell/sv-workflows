# Validation: setup is not accuracy

## Fast offline tests

`python scripts/test_installation.py` runs unit tests for configuration/schema errors, paired-role handling, scientific option defaults, reference-cache reuse, VCF orientations/mate collapse, filtering adapters, BRASS statistics, optimizer convergence guards, and memory/exit-code recording. Fixtures are synthetic and require no EGA files or licensed software.

`python scripts/workflow.py test --profile test` then executes a real Nextflow task on a synthetic VCF. It checks three retained events: deletion, duplication and one interchromosomal adjacency represented by two reciprocal BND records. This test does not run a structural variant caller. On an HPC, launch it within a small compute allocation.

## Actual entry-point and numerical controls

`python scripts/test_main_stub.py --out results/main_stub_check01` executes `main.nf -stub-run` using a tiny reference and real indexed synthetic tumour/normal BAMs. Preflight is not stubbed. It checks numbered-contig queries, symlinked inputs, four caller branches, cached resume, scheduler-environment recognition, and rejection of unacknowledged local preflight. Add `--manta-image /absolute/path/to/manta.sif` for the fifth branch and real Manta container configuration. Stub tasks do not measure performance or establish biological accuracy. The synthetic assembly exception exists only in stub runs.

GitHub CI runs that four-branch control and the pinned original ClusterSV positive control, including byte-for-byte equivalence against upstream defaults and empty/singleton input cases. CI uses Ubuntu distribution tools and R packages; the production R environment remains separately locked. This is not execution of every caller or validation of non-Slurm schedulers.

`--with-exclusions` runs a separate region-policy stub check for DELLY/GRIDSS, and Manta when its image is supplied. The Manta check verifies the inclusion-list complement and requires `bgzip`/`tabix`; it does not establish biological equivalence of caller region handling.

Specify `--container-engine singularity` or `apptainer` when needed for the optional Manta check. `--debug-cache` adds detailed task-hash diagnostics when troubleshooting resume; normal test logs stay compact.

## Container/setup checks

Pull selected digest-pinned images, check the caller executable inside each with `apptainer exec IMAGE command --help`, and validate the exact alignment reference/indexes. BRASS also requires native bundle validation. ClusterSV and SigProfiler have genuine numerical controls that need their pinned R/Python environments; they are not part of the lightweight offline suite.

## Before a new cohort or release

1. Run the unchanged original and consolidated workflow on the same matched pair, comparing normalized biological events, sample roles, filters and resolved commands. Differences require review, not automatic acceptance.
2. Normal-self negative control for each caller, using two read-only aliases with distinct sample labels and **distinct read names** where the caller requires them. A low count is expected but is not a guaranteed zero; investigate every retained PASS candidate. Never interpret an input-label error as successful control validation.
3. Positive benchmark with a suitable truth set and confident regions. GIAB germline SV truth measures discovery, not tumour-specific somatic precision; do not treat a normal-self GIAB run as a positive somatic benchmark. Use validated somatic spike-in/truth or report the scope honestly.
4. Tumour/normal swap and repeated-run checks. Compare event-normalized VCFs, not date headers alone. Document assembly/threading variation.
5. Review emitted commands after parameter/version changes, plus representative calls in IGV. Inspect reference gaps, centromeric noise, adapter artifacts, depth and matched-normal evidence.

Only then expand to three concurrent pairs, then measured cohort concurrency. CI passing does not waive biological validation. Reference/tool/filter changes create a new configuration requiring revalidation.
