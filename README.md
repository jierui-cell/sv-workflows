# *SV Workflows*

*SV Workflows*: A Nextflow workflow for paired tumour–normal structural-variant calling and SV signature analysis.

[![Tests](https://github.com/jierui-cell/sv-workflows/actions/workflows/tests.yml/badge.svg?branch=main)](https://github.com/jierui-cell/sv-workflows/actions/workflows/tests.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Nextflow](https://img.shields.io/badge/Nextflow-25.10.7-blue)](https://www.nextflow.io/)
[![Python](https://img.shields.io/badge/Python-3.11-blue)](requirements.txt)
[![Status](https://img.shields.io/badge/Status-research%20development-orange)](docs/workflow-status.md)

[Requirements](#requirements) · [Usage](#usage-instructions) · [Toy test](#test-your-installation) · [Signatures](#sv-signature-analysis) · [Results](#retrieve-your-results) · [Parameters](CONFIG_PARAMETERS.md)

## Overview

This open-source workflow runs multiple SV callers on matched tumour and normal BAMs, records settings and resource usage, and prepares caller-specific outputs for SV signature analysis. The same workflow and configuration format can be reused across cohorts and HPC systems.

**Key features:**

- Five core callers: **Manta, DELLY, GRIDSS, SvABA and BRASS**.
- Two analysis routes: **COSMIC SV32 refitting** or **ClusterSV clustering without copy-number data**.
- YAML defaults and user overrides, pinned software identities, input checks and resume support.
- Synthetic setup tests, caller QC, command logs and performance records.

```text
Matched tumour + normal BAMs → Input checks → Selected SV callers → VCFs + QC
                                                                    │
                                                    Choose analysis method
                                                      ├── COSMIC SV32
                                                      │     → SV1–SV12 exposures
                                                      └── ClusterSV
                                                            → Clusters, footprints
                                                              and timing annotations
```

**Research use only.** This public version has passed synthetic tests; validate it on a real matched tumour–normal pair before expanding to a cohort. ClusterSV without copy-number data (CN-free) reports clusters and annotations, rather than reference-signature exposures. See [readiness and limitations](docs/workflow-status.md).

## Requirements

### 0. Clone the repository

```bash
git clone https://github.com/jierui-cell/sv-workflows.git
cd sv-workflows
```

Run the commands below from this directory. **Keep BAMs and clinical metadata outside the repository**; the workflow reads them by path.

### 1. Install the launcher and workflow engine

You need Linux, **Python 3.11**, **Java 17+**, and **Nextflow 25.10.7**. Load your HPC's modules where available; names differ between sites. Install Nextflow using its [official instructions](https://www.nextflow.io/docs/latest/install.html), or load a module providing the required version.

```bash
python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
nextflow -version
python scripts/test_installation.py
```

The offline tests should finish with `OK`. Real caller tasks additionally need **samtools**, **bcftools**, and **Apptainer or Singularity** in the compute environment. Building the reference indexes below also requires **BWA**. BRASS has extra prerequisites; see [installation details](docs/setup.md).

### Test Your Installation

The bundled synthetic VCF contains a deletion, a duplication and an interchromosomal adjacency. On a workstation or inside a small compute allocation:

```bash
python scripts/workflow.py test --profile test
```

On Slurm, run that same local toy task **inside an allocated compute-node driver**:

```bash
mkdir -p logs
sbatch --account=YOUR_ACCOUNT --partition=YOUR_PARTITION \
  --cpus-per-task=1 --mem=4G --time=00:15:00 --output=logs/toy-%j.log \
  --wrap="python scripts/workflow.py test --profile test"
```

The expected output is `results/toy/PASS.txt` and **three retained events**. See [expected toy outputs](examples/toy_expected/README.md). This checks installation and conversion—not caller sensitivity or somatic accuracy.

For stronger installation validation, exercise the **actual caller entry point**, with real synthetic BAM preflight and stubbed caller tasks:

```bash
# Workstation or allocated compute node; requires samtools and Nextflow.
python scripts/test_main_stub.py --out results/main_stub_check01
```

This checks four caller branches plus cached resume, without running caller inference. After downloading Manta's image, add `--manta-image containers/images/manta-1.6.0--py27h9948957_6.sif` to also test its container configuration. Use a fresh output directory each time. CI runs the four-branch test and the original ClusterSV numerical control; the container-enabled five-branch control is checked separately on HPC.

### 2. Download the caller containers

Choose your site's runtime. Both commands download the same digest-pinned images into `containers/images/`.

**Option A — Apptainer**

```bash
python scripts/fetch_containers.py --engine apptainer --tools manta,delly,svaba,gridss,brass
```

**Option B — Singularity**

```bash
python scripts/fetch_containers.py --engine singularity --tools manta,delly,svaba,gridss,brass
```

Set `container_engine` accordingly. You can download only the tools you plan to run, for example `--tools manta,delly`. Pull images on a node where network access and temporary storage are permitted. Each caller has its own image; there is no single all-in-one container.

### 3. Prepare the alignment reference

Use the **exact FASTA used to align your BAMs**, not an arbitrary download with the same assembly name. Prepare indexes on a compute node:

```bash
samtools faidx /path/to/reference.fa
bwa index /path/to/reference.fa
samtools dict /path/to/reference.fa > /path/to/reference.dict
```

Your reference folder should look like this:

```text
reference/
├── reference.fa
├── reference.fa.fai       # All callers
├── reference.fa.amb       # SvABA / GRIDSS BWA indexes
├── reference.fa.ann
├── reference.fa.bwt
├── reference.fa.pac
├── reference.fa.sa
└── reference.dict         # GRIDSS
```

BRASS also requires a validated **chr-prefixed GRCh37 resource bundle**; the FASTA alone is insufficient. Follow [BRASS bundle setup](docs/setup.md#brass-bundle) before enabling it. GRCh38 compatibility is not currently established.

## Usage Instructions

Your working directory will look approximately like this:

```text
sv-workflows/
├── main.nf                  # SV calling
├── sv_signatures/           # Signature/feature analysis
├── modules/                 # Reusable caller modules
├── scripts/                 # Launcher, adapters and QC
├── config/
│   ├── defaults.yaml        # Versioned defaults; do not edit per cohort
│   └── private-cohort.yaml  # Your paths/settings; ignored by Git
├── examples/                # Synthetic inputs and copyable configs
├── containers/images/      # Downloaded images; ignored by Git
├── results/                 # Published outputs; ignored by Git
├── logs/                    # Settings, trace and failures; ignored by Git
└── work/                    # Nextflow resume cache; ignored by Git
```

### 1. Create your paired samplesheet

Choose an existing private folder outside this checkout, then copy the template (replace the destination below):

```bash
cp examples/samplesheet.tsv /path/to/private/paired_samples.tsv
```

Edit that file, **preserving its tabs and all nine column headers**. The template is already tab-separated; do not replace tabs with spaces or remove trailing empty metadata fields. Replace the synthetic row with your unique pair IDs and absolute BAM/index paths in the first five fields. Metadata values may be blank, but their headers must remain present. Missing purity/ploidy does not block the five core callers; BRASS explicitly records native-default values. See [input rules](docs/inputs-outputs.md).

### 2. Copy and edit the run configuration

```bash
cp examples/cohort.yaml config/private-cohort.yaml
```

Edit the samplesheet/reference paths, cohort label, scheduler settings and enabled callers. The template starts with **Manta and DELLY**; add `svaba,gridss,brass` once their resources and pilots are ready.

Set `genome_build` explicitly (`GRCh37` by default). Preflight checks primary chromosome lengths against that declaration as well as BAM/reference dictionary agreement; this is not a replacement for supplying the exact alignment FASTA. Optional region exclusions and their caller-specific behavior are explained in [parameters](CONFIG_PARAMETERS.md#optional-region-exclusions).

All scientific/resource defaults are in [config/defaults.yaml](config/defaults.yaml). Your private YAML only needs values you want to change. See [CONFIG_PARAMETERS.md](CONFIG_PARAMETERS.md) for names, defaults and tuning examples.

### 3. Check the configuration and inputs

Inspect resolved settings without submitting jobs:

```bash
python scripts/workflow.py callers --config config/private-cohort.yaml --dry-run
```

Then run BAM/reference preflight on compute nodes. For **Slurm**:

```bash
mkdir -p logs
sbatch --account=YOUR_ACCOUNT --partition=YOUR_PARTITION \
  --cpus-per-task=1 --mem=4G --time=4:00:00 --output=logs/preflight-%j.log \
  scripts/run.slurm.sh callers --config config/private-cohort.yaml --set callers=preflight
```

Replace submission flags **and** matching YAML account/partition with your site's values. Activate the environment and load required modules before submission. `--dry-run` alone does not inspect BAM contents.

Submit these commands **from the repository root**. The driver uses `SLURM_SUBMIT_DIR`, not the location of Slurm's copied script. To submit from elsewhere, export `SV_WORKFLOW_ROOT=/absolute/path/to/sv-workflows`; optionally set `SV_PYTHON` to the absolute launcher-environment Python executable.

### 4. Run a pilot, then expand

Run one pair first, replacing `PAIR_001` with an ID from your samplesheet:

```bash
sbatch --account=YOUR_ACCOUNT --partition=YOUR_PARTITION \
  --cpus-per-task=1 --mem=4G --time=48:00:00 --output=logs/pilot-%j.log \
  scripts/run.slurm.sh callers --config config/private-cohort.yaml \
  --set pair_ids=PAIR_001 --set outdir=results/pilot --set max_parallel_callers=1
```

After reviewing the pilot and [control tests](docs/validation.md), run the selected cohort in a new output directory:

```bash
sbatch --account=YOUR_ACCOUNT --partition=YOUR_PARTITION \
  --cpus-per-task=1 --mem=4G --time=7-00:00:00 --output=logs/cohort-%j.log \
  scripts/run.slurm.sh callers --config config/private-cohort.yaml \
  --set outdir=results/cohort_run01 --set max_parallel_callers=3
```

Concurrency is **per caller**, not global. Raise it to `8` only after a successful concurrent pilot and capacity checks. Driver time must cover the whole run. See [HPC instructions](docs/hpc.md) for PBS/LSF/SGE or an allocated local node; non-Slurm profiles are templates requiring site validation.

## SV Signature Analysis

Start from completed, reviewed caller-specific VCFs. Create a private five-column inventory following [examples/callsets.tsv](examples/callsets.tsv), then calculate source hashes:

```bash
mkdir -p manifests
python scripts/build_sv_manifest.py --inventory /path/to/private/callsets.tsv \
  --output manifests/sv_manifest.tsv
cp examples/signatures.yaml config/private-signatures.yaml
```

Edit that config and install the separate [signature environments/references](docs/signatures.md). Both routes currently use **GRCh37**; other assemblies need explicit implementation and validation changes.

| Choose this route | You receive | Interpretation |
| --- | --- | --- |
| COSMIC SV32 | 32-channel counts and SV1–SV12 exposures | Fit to reference signatures |
| ClusterSV | Clusters, footprints and timing annotations | Clustering and annotation without copy-number data |

Choose **one** analysis route:

```bash
# COSMIC: 32-channel matrix and SV1–SV12 exposures.
sbatch --account=YOUR_ACCOUNT --partition=YOUR_PARTITION \
  --cpus-per-task=1 --mem=4G --time=24:00:00 --output=logs/cosmic-%j.log \
  scripts/run.slurm.sh signatures --config config/private-signatures.yaml --set method=cosmic_sv32

# ClusterSV: CN-free clustering, footprints and timing annotations.
sbatch --account=YOUR_ACCOUNT --partition=YOUR_PARTITION \
  --cpus-per-task=1 --mem=4G --time=24:00:00 --output=logs/clustersv-%j.log \
  scripts/run.slurm.sh signatures --config config/private-signatures.yaml \
  --set method=pcawg_clustersv --set pcawg_stage=clustering --set outdir=results/clustersv_run01
```

The methods have different inclusion rules and outputs. Do not equate their event totals or interpret CN-free clusters as established complex-event mechanisms. Only COSMIC SV32 fits reference signatures; the ClusterSV route returns clustering and annotation results.

## Retrieve Your Results

Under your configured `outdir`:

| Result | Location |
| --- | --- |
| Input/reference checks | `preflight/PAIR/` |
| Manta somatic SVs | `manta/PAIR/manta/results/variants/somaticSV.vcf.gz` |
| DELLY somatic SVs | `delly/PAIR/output/somaticSV.vcf.gz` |
| SvABA somatic SVs | `svaba/PAIR/output/PAIR.svaba.somatic.sv.vcf` |
| GRIDSS somatic SVs | `gridss/PAIR/output/somatic/somatic.high_confidence.vcf.gz` and `somatic.full.vcf.gz` |
| GRIDSS raw joint calls | `gridss/PAIR/output/joint.vcf.gz` — **not a somatic callset** |
| BRASS annotated SVs | `brass/PAIR/output/callset/calls/*.annot.vcf.gz` |
| Runtime/memory | Manta: `manta/PAIR/benchmark.json`; others: `CALLER/PAIR/output/benchmark.json` |
| COSMIC matrices/exposures | `COHORT/CALLER/group_output/SV32.matrix.tsv` and `SV32.exposures.tsv` |
| ClusterSV clusters/footprints | `COHORT/CALLER/group_output/clusters/` and `cluster_summary.tsv` |
| Resolved settings and failures | `logs/WORKFLOW-TIMESTAMP/` in the repository |

Paths use `compact_publish: false`; signature smoke/pilot selections add grouping directories. See [output details](docs/inputs-outputs.md) and [performance interpretation](docs/performance.md).

## Configuration and Further Documentation

- [Configuration parameters](CONFIG_PARAMETERS.md)
- [Installation and BRASS resources](docs/setup.md)
- [HPC profiles, concurrency and resume](docs/hpc.md)
- [Input/output rules](docs/inputs-outputs.md)
- [Signature methods and references](docs/signatures.md)
- [Readiness and limitations](docs/workflow-status.md)
- [Validation and controls](docs/validation.md)
- [Development and reproducibility](docs/development.md)

## Common Issues

- **Missing column/invalid pair ID:** use the exact TSV headers, preserve blank metadata fields, and check selected IDs.
- **Missing image/index/bundle:** finish the relevant setup; never bypass validation markers or substitute an unrelated reference.
- **Command/module unavailable:** load tools in the compute environment or configure site module names in YAML.
- **Out of memory/no temporary space:** inspect task/native logs. RAM and scratch storage are different; page-cache-inclusive Slurm MaxRSS is not application RSS.
- **Restart after failure:** fix the cause and use `--resume SESSION_UUID` with original settings/work files. Changed scientific parameters require a new output directory.

## License

Our independent code is [MIT-licensed](LICENSE). Upstream tools, scientific sources and reference resources retain their own terms; see [third-party notices](THIRD_PARTY.md).

Only code, documentation and synthetic examples are published. **Do not commit real BAMs/VCFs, sample metadata, private configs, credentials or runtime logs.**
