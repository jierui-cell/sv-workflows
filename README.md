# SV Workflows

**From matched BAMs to structural variants and SV signatures.**

[![Tests](https://github.com/jierui-cell/sv-workflows/actions/workflows/tests.yml/badge.svg?branch=main)](https://github.com/jierui-cell/sv-workflows/actions/workflows/tests.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Nextflow](https://img.shields.io/badge/Nextflow-25.10.7-blue)](https://www.nextflow.io/)

[Setup](#requirements) · [Launch](#usage-instructions) · [Signatures](#sv-signature-analysis) · [Results](#retrieve-your-results) · [Parameters](CONFIG_PARAMETERS.md)

## Overview

Run **Manta, DELLY, GRIDSS, SvABA and BRASS** on matched tumour–normal samples, then fit **COSMIC SV32 signatures**. YAML controls callers, parameters, HPC modules and resources. Nextflow schedules the work and saves VCFs, QC, commands and runtime/memory records.

## Standard Folder Layout

~~~text
sv-workflows/
├── config/private-cohort.yaml   # Your cohort and HPC settings
├── inputs/samplesheet.tsv       # Matched tumour–normal pairs
├── inputs/bams/                 # BAMs and indexes, or links to existing files
├── resources/reference/reference.fa  # Shared alignment FASTA and indexes
├── resources/brass/bundle/      # Additional BRASS resources
├── containers/images/          # Downloaded caller images
├── results/                    # Caller and signature results
├── logs/                       # Settings, manifests, commands and run logs
└── work/                       # Nextflow task files for resume
~~~

Relative paths start at the repository root. Absolute paths also work.

## Requirements

Linux, Python 3.11, Java 17+, Nextflow 25.10.7, samtools, bcftools and Apptainer or Singularity. Caller software runs in pinned container images.

### 1. Install once

~~~bash
git clone https://github.com/jierui-cell/sv-workflows.git
cd sv-workflows
python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
cp examples/cohort.yaml config/private-cohort.yaml
curl -fL --retry 3 https://github.com/nextflow-io/nextflow/releases/download/v25.10.7/nextflow \
  -o .venv/bin/nextflow
chmod +x .venv/bin/nextflow
~~~

Set your HPC module names and container runtime in `config/private-cohort.yaml`, for example:

~~~yaml
modules: [java/20.0.1, samtools/1.20, bcftools/1.20]
container_engine: singularity
~~~

Use the module names supplied by your HPC, including a container-runtime module if needed. Use `modules: []` when tools are already installed. The launcher and compute jobs load this list automatically.

Initialize Nextflow once with internet access:

~~~bash
python scripts/hpc_environment.py -- .venv/bin/nextflow -version
~~~

### 2. Download caller images

~~~bash
python scripts/fetch_containers.py
~~~

This uses your cohort YAML and writes the pinned images to `containers/images/`.

### 3. Prepare the alignment reference

Place the FASTA used to align your BAMs at `resources/reference/reference.fa`, with its matching FASTA/BWA indexes and sequence dictionary. Existing files can be linked here. [Index commands](docs/setup.md#alignment-reference).

### 4. Set up BRASS's extra resources (once)

BRASS uses the **same human FASTA** as the other callers, plus the resources in its [reference guide](https://github.com/cancerit/BRASS/wiki). On an allocated compute node:

~~~bash
python scripts/setup_brass.py
~~~

The command downloads Sanger [BRASS resources](https://ftp.sanger.ac.uk/pub/cancer/dockstore/human/CNV_SV_ref_GRCh37d5_brass6+.tar.gz) and [VAGrENT annotations](https://ftp.sanger.ac.uk/pub/cancer/dockstore/human/VAGrENT_ref_GRCh37d5_ensembl_75.tar.gz) (~3.1 GB) and prepares `resources/brass/bundle/`. [Why preparation is needed](docs/setup.md#why-preparation-is-needed).

The BRASS setup supports compatible GRCh37 references with chromosome names such as `chr1` and `chr2`. For other reference builds, adapt [`scripts/setup_brass.py`](scripts/setup_brass.py) and its preparation helpers to use matching resources and chromosome names. [Setup details](docs/setup.md#brass-bundle).

## Usage Instructions

### 1. Supply your samples and settings

~~~bash
mkdir -p inputs/bams inputs/regions logs
cp examples/samplesheet.tsv inputs/samplesheet.tsv
~~~

Fill the samplesheet with unique pair IDs and BAM/index paths. Keep all nine headers; the first five columns are required. Sex, phenotype, purity and ploidy can be blank. [Samplesheet format](docs/inputs-outputs.md#paired-samplesheet).

In `config/private-cohort.yaml`, set your cohort label, output folder, account/partition and callers. The template selects all five callers with eight concurrent pairs per caller. Scientific defaults are in [config/defaults.yaml](config/defaults.yaml).

### 2. Launch the full cohort

From the repository root:

~~~bash
sbatch --account=YOUR_ACCOUNT --partition=YOUR_PARTITION \
  --cpus-per-task=1 --mem=4G --time=7-00:00:00 --output=logs/cohort-%j.log \
  scripts/run.slurm.sh callers --config config/private-cohort.yaml
~~~

The submission script selects `.venv`, and the launcher loads your YAML modules before starting Nextflow:

~~~text
sbatch → run.slurm.sh → workflow.py → load modules → Nextflow → compute jobs
~~~

Nextflow loads the same module list in each compute task, checks inputs and runs the selected callers. Set driver walltime to cover the cohort within your HPC limits. [Concurrency and other schedulers](docs/hpc.md).

### 3. Optional: preview or change settings

Preview resolved settings and the Nextflow command:

~~~bash
python scripts/workflow.py callers --config config/private-cohort.yaml --dry-run
~~~

Use `--set` for one-off overrides, such as `--set delly.options.map_qual=2`. For lasting changes, edit your cohort YAML. Give parameter experiments a new `outdir`. [Parameter reference](CONFIG_PARAMETERS.md).

### 4. Optional: restart an interrupted run

List previous runs using your configured modules:

~~~bash
python scripts/hpc_environment.py -- nextflow log
~~~

Copy the original run’s `SESSION ID`, then add it to the original submission:

~~~bash
sbatch --account=YOUR_ACCOUNT --partition=YOUR_PARTITION \
  --cpus-per-task=1 --mem=4G --time=7-00:00:00 --output=logs/resume-%j.log \
  scripts/run.slurm.sh callers --config config/private-cohort.yaml --resume SESSION_UUID
~~~

Keep the original inputs/settings, `work/` and `.nextflow/`; matching completed tasks are reused.

## Exome and Targeted Sequencing

Use the same samplesheet and launcher with **Manta**. Copy your cohort YAML to `config/private-capture.yaml` and update its input, cohort label and these settings:

~~~yaml
sequencing_type: WES
callers: manta
target_bed: null  # Optional capture-kit BED for WES.
outdir: results/exome_run01
~~~

For panels, use `sequencing_type: targeted` and supply `target_bed: inputs/regions/panel_targets.bed`. BED starts are 0-based and ends excluded; chromosome names match the reference.

Submit using the launch command above with the capture YAML. Manta receives `--exome` automatically. Preparing a supplied BED uses bgzip/tabix. COSMIC SV32 uses WGS callsets. [Assay support](docs/workflow-status.md#assay-support).

## SV Signature Analysis

Signatures automatically read the completed VCFs from the run configured in `config/private-cohort.yaml`. GRIDSS contributes its high-confidence somatic VCF.

~~~bash
cp examples/signatures.yaml config/private-signatures.yaml
~~~

Set the signature account/partition, complete the one-time [signature environment setup](docs/signatures.md#install-once), then submit:

~~~bash
sbatch --account=YOUR_ACCOUNT --partition=YOUR_PARTITION \
  --cpus-per-task=1 --mem=4G --time=24:00:00 --output=logs/cosmic-%j.log \
  scripts/run.slurm.sh signatures --config config/private-signatures.yaml
~~~

The launcher inherits the cohort modules, builds the hashed VCF manifest and starts COSMIC SV32 counts → SV1–SV12 fitting. Groups run independently by cohort/caller. [Input options and outputs](docs/signatures.md).

PCAWG/ClusterSV signature analysis is on hold.

## Retrieve Your Results

Under the caller `outdir`:

| Result | Location |
| --- | --- |
| Input checks | `preflight/PAIR/` |
| Manta | `manta/PAIR/manta/results/variants/somaticSV.vcf.gz` |
| DELLY | `delly/PAIR/output/somaticSV.vcf.gz` |
| SvABA | `svaba/PAIR/output/PAIR.svaba.somatic.sv.vcf` |
| GRIDSS | `gridss/PAIR/output/somatic/somatic.high_confidence.vcf.gz` |
| BRASS | `brass/PAIR/output/callset/calls/*.annot.vcf.gz` |
| Runtime / memory | Manta: `manta/PAIR/benchmark.json`; others: `CALLER/PAIR/output/benchmark.json` |

Signature results: `SIGNATURE_OUTDIR/COHORT/CALLER/group_output/SV32.matrix.tsv` and `SV32.exposures.tsv`. Settings and logs: `logs/WORKFLOW-TIMESTAMP/`. [Performance exports](docs/performance.md).

## Test Your Installation

Run optional checks on a workstation or allocated compute node. The Nextflow tests load the cohort YAML modules automatically:

~~~bash
python scripts/test_installation.py
python scripts/workflow.py test --profile test
python scripts/test_main_stub.py --out results/main_stub_check01
~~~

These check Python logic, toy SV conversion, caller task wiring and cached resume. [Test guide](docs/validation.md).

## Further Documentation

[Setup](docs/setup.md) · [HPC](docs/hpc.md) · [Parameters](CONFIG_PARAMETERS.md) · [Signatures](docs/signatures.md) · [Method scope](docs/workflow-status.md) · [Development](docs/development.md)

## License

Independent code is [MIT-licensed](LICENSE); caller licenses are listed in [upstream notices](THIRD_PARTY.md). Research use. Only source and synthetic examples are published.
