# SV Workflows

**Launch a cohort. Compare five callers. Discover SV signatures.**

A reusable Nextflow framework for matched tumour–normal structural-variant analysis. Configure once in YAML, run across your cohort, and retrieve caller-specific VCFs, QC and performance tables.

## Start Here

1. [Install tools and prepare resources](setup.md).
2. Copy your samplesheet and configuration, then [launch the full cohort](../README.md#usage-instructions).
3. Choose [COSMIC reference fitting or ClusterSV + de novo NMF](signatures.md).

Input validation is built into every run. There is no separate pilot requirement.

## What You Can Run

| Analysis | Output |
| --- | --- |
| WGS: Manta, DELLY, GRIDSS, SvABA, BRASS | Paired SV calls, QC, runtime and memory |
| WES / targeted: Manta | Capture-aware paired SV calls with the same launcher |
| COSMIC SV32 | 32-channel counts and SV1–SV12 exposures |
| ClusterSV + NMF | Clusters, footprints, CN-free counts and de novo signatures; optional timing |

[Parameters](../CONFIG_PARAMETERS.md) · [HPC deployment](hpc.md) · [Inputs and results](inputs-outputs.md) · [Method scope](workflow-status.md) · [Tests](validation.md) · [Performance](performance.md)
