# SV Workflows

Paired tumour–normal SV calling and selectable signature analysis for reproducible HPC research.

## Start Here

1. Follow the numbered [quick-start guide](../README.md): install tools, pull containers and prepare the exact alignment reference.
2. Copy a run configuration, prepare your samplesheet, and run preflight followed by a one-pair pilot.
3. Retrieve the reviewed caller VCFs, then choose COSMIC SV32 or ClusterSV analysis.

The [configuration guide](../CONFIG_PARAMETERS.md) lists defaults and user overrides. [HPC instructions](hpc.md) explain scheduler settings, concurrency, storage and resume.

## What You Can Run

| Task | Available workflow | Guide |
| --- | --- | --- |
| Call SVs from paired BAMs | Manta, DELLY, GRIDSS, SvABA, BRASS | [Quick start](../README.md#usage-instructions) |
| Fit reference SV signatures | COSMIC SV32 / SV1–SV12 | [Signature methods](signatures.md) |
| Cluster junctions without copy number | ClusterSV, footprints, timing annotations | [Signature methods](signatures.md) |
| Check setup with synthetic data | Offline tests and Nextflow toy workflow | [Setup test](../README.md#test-your-installation) |

Cue and JaBbA are research extensions; LUMPY is discontinued. See [readiness and limitations](workflow-status.md). No unavailable copy-number-dependent mechanism is fabricated.

This is a development release based on workflows previously used for paired tumour–normal WGS. The portable refactor still requires a real-pair equivalence pilot; successful workflow completion does not establish clinical accuracy. Public examples are synthetic; no patient data, cohort results or clinical metadata are distributed.

[Performance measurement](performance.md) · [Development and reproducibility](development.md) · [Source repository](https://github.com/jierui-cell/sv-workflows)
