# Workflow Readiness and Limitations

The core workflows were used on paired tumour–normal WGS. The public refactor has passed synthetic tests but still needs a real-pair equivalence pilot before cohort expansion. Technical completion is not clinical validation.

| Workflow | Output | Current boundary |
| --- | --- | --- |
| Manta, DELLY, SvABA | Paired somatic SV VCFs and QC | Core workflows; validate refactor equivalence and controls at your site |
| GRIDSS | Joint VCF plus native somatic subsets | No panel of normals or BSgenome in the baseline; raw joint calls are not somatic calls |
| BRASS | Annotated paired SV VCF/BEDPE | Validated chr-prefixed GRCh37 bundle; missing purity/ploidy explicitly use native defaults |
| Cue | Predictions and provisional paired candidates | Research-only; resource provisioning is not turnkey |
| JaBbA | Copy-number/junction graph | Licensed solver; toy passed, real pilot did not converge; not an independent BAM-to-VCF caller |
| LUMPY | Historical helpers/tests only | Discontinued; launcher rejects it |
| COSMIC SV32 | 32-channel matrix and SV1–SV12 exposures | Pinned numerical clustering control required; caller-specific, not consensus truth |
| CN-free ClusterSV | Clusters, footprints, timing annotations | Original proximity method; not CN-dependent classification or the complete 44 features |
| Advanced 44 features | Matrix from supplied validated classifications | Separate newer-paper profile; overlapping counts; no fabricated classes |

## Backend Validation

Slurm is production-tested. PBS, LSF and SGE are site templates, not verified portability claims. Local execution is for a workstation or allocated compute node, never an HPC login node. Each new site needs a setup test and pilot.

## Scientific Controls

Use normal-self controls, suitable positive benchmarks, role-swap/determinism tests, resolved-command review and representative read-level review. See [validation](validation.md). Retained normal-self candidates require investigation; `SUCCESS` does not establish somatic correctness.

The original PCAWG 16-signature reference is not included; reference decomposition is deferred. Replication-timing tertiles are a recorded reconstruction policy, not verified original-paper category thresholds.
