# Workflow Readiness and Limitations

The core workflows were used on paired tumour–normal WGS. The public refactor has passed synthetic tests but still needs a real-pair equivalence pilot before cohort expansion. Technical completion is not clinical validation.

| Workflow | Output | Current boundary |
| --- | --- | --- |
| Manta, DELLY, SvABA | Paired somatic SV VCFs and QC | Core workflows; validate refactor equivalence and controls at your site |
| GRIDSS | Joint VCF plus native somatic subsets | No panel of normals or BSgenome in the baseline; raw joint calls are not somatic calls |
| BRASS | Annotated paired SV VCF/BEDPE | Validated chr-prefixed GRCh37 bundle; missing purity/ploidy explicitly use native defaults |
| COSMIC SV32 | 32-channel matrix and SV1–SV12 exposures | Pinned numerical clustering control required; caller-specific, not consensus truth |
| CN-free ClusterSV | Clusters, footprints, timing annotations | Original proximity method; no copy-number-dependent mechanism inference |

## Backend Validation

Slurm is production-tested. PBS, LSF and SGE are site templates, not verified portability claims. Local execution is for a workstation or allocated compute node, never an HPC login node. Each new site needs a setup test and pilot.

## Scientific Controls

Use normal-self controls, suitable positive benchmarks, role-swap/determinism tests, resolved-command review and representative read-level review. See [validation](validation.md). Retained normal-self candidates require investigation; `SUCCESS` does not establish somatic correctness.

Replication-timing tertiles are a recorded reconstruction policy, not verified original-paper category thresholds.
