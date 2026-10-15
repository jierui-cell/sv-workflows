# SV Workflows

**From matched BAMs to structural variants and SV signatures—one framework, any cohort.**

[![Tests](https://github.com/jierui-cell/sv-workflows/actions/workflows/tests.yml/badge.svg?branch=main)](https://github.com/jierui-cell/sv-workflows/actions/workflows/tests.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Nextflow](https://img.shields.io/badge/Nextflow-25.10.7-blue)](https://www.nextflow.io/)
[![Python](https://img.shields.io/badge/Python-3.11-blue)](requirements.txt)

[Folders](#standard-folder-layout) · [Launch a cohort](#usage-instructions) · [Signatures](#sv-signature-analysis) · [Exome / panels](#exome-and-targeted-sequencing) · [Results](#retrieve-your-results) · [Parameters](CONFIG_PARAMETERS.md)

## Overview

Run **Manta, DELLY, GRIDSS, SvABA and BRASS** across a matched tumour–normal cohort. Nextflow handles scheduling, automatic input checks, parallel execution and resume. YAML gives you one place to select callers, tune parameters and size resources. Each run saves its settings, commands, QC and comparable runtime/memory records.

Use **COSMIC SV32 reference fitting** for SV signature analysis. **PCAWG/ClusterSV signature analysis is on hold**, including its feature/NMF extension. The implementation is preserved for later review; that route is not currently being run.

~~~text
Paired BAMs → Automatic checks → Five SV callers → VCFs + QC + benchmarks
                                                       │
                                                   COSMIC SV32
                                                   32-channel counts
                                                   → SV1–SV12 exposures
~~~

The framework is for research analysis. [Method scope](docs/workflow-status.md) describes reference builds, assay support and output interpretation.

## Standard Folder Layout

The defaults use this layout inside your checkout—no resource-path edits are needed when you follow it:

~~~text
sv-workflows/
├── inputs/
│   ├── samplesheet.tsv       # Your paired samples
│   ├── bams/                 # BAMs + indexes; existing files may be linked here
│   ├── regions/              # Optional targets/exclusions
│   └── callsets.tsv          # VCF inventory for signature analysis
├── resources/
│   ├── reference/reference.fa  # Alignment FASTA + matching indexes alongside it
│   └── brass/bundle/         # Extra BRASS annotations/masks; NOT a second human genome
├── containers/images/        # Created/populated by the image downloader
├── manifests/sv_manifest.tsv # Generated signature manifest
├── config/private-cohort.yaml
├── results/                  # Created by the workflow
├── logs/                     # Settings, setup checks and run logs
└── work/                     # Task files for resume
~~~

Relative paths in YAML and BAM paths in the samplesheet are resolved from the **repository root**, not the YAML file's folder. Absolute paths remain supported through YAML overrides. Inputs, resources, images and runtime outputs are Git-ignored; only source and synthetic examples are published.

## Requirements

### 1. Install once

Use Linux, Python 3.11, Java 17+ and Nextflow 25.10.7 from the [official launcher](https://www.nextflow.io/docs/latest/install.html) or your HPC modules. Compute tasks need samtools, bcftools and Apptainer or Singularity.

~~~bash
git clone https://github.com/jierui-cell/sv-workflows.git
cd sv-workflows
python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
~~~

### 2. Download caller images

~~~bash
python scripts/fetch_containers.py --engine apptainer --tools manta,delly,gridss,svaba,brass
~~~

Use `--engine singularity` for Singularity sites. Digest-pinned images go into `containers/images/`. Download only your selected callers if preferred, for example `--tools manta,delly`.

### 3. Prepare the alignment reference

Place the FASTA used to align your BAMs at **`resources/reference/reference.fa`**, with its matching indexes in the same folder. This filename is a convention, not a different reference assembly. Reuse the existing indexed reference or create indexes on a compute node with samtools and BWA:

~~~bash
mkdir -p resources/reference
samtools faidx resources/reference/reference.fa
bwa index resources/reference/reference.fa
samtools dict resources/reference/reference.fa > resources/reference/reference.dict
~~~

### 4. Set up BRASS's extra resources (once)

**All five callers use the same human FASTA above.** BRASS additionally needs annotation/filtering files in `resources/brass/bundle/`: depth/GC masks, gene/chromosome annotations and microbial/viral sequences. This is **not a second human reference**.

After preparing the shared reference, run this command from the checkout root on an allocated compute node with internet access:

~~~bash
python scripts/setup_brass.py
~~~

It downloads the two official Sanger archives (~3.1 GB), verifies their checksums, extracts the needed files, adapts chromosome labels to your shared reference and runs BRASS's native resource checks. It also downloads the pinned BRASS image if absent. The resulting folder already matches the YAML default; no path edits are needed. Repeat the command to verify/reuse an existing bundle.

**Using Singularity?** Run:

~~~bash
python scripts/setup_brass.py --engine singularity
~~~

Also set `container_engine: singularity` in your cohort YAML so caller jobs use Singularity too.

This setup supports **GRCh37 with `chr` chromosome names** (for example, `chr1`). Skip this step if you are not running BRASS. [Custom paths and other setup options](docs/setup.md#brass-bundle).

## Usage Instructions

### 1. Supply your samples and settings

~~~bash
mkdir -p inputs/bams inputs/regions logs
cp examples/samplesheet.tsv inputs/samplesheet.tsv
cp examples/cohort.yaml config/private-cohort.yaml
~~~

Place or link your BAMs and indexes in `inputs/bams/`. Edit `inputs/samplesheet.tsv` with unique pair IDs and the corresponding filenames; paths such as `inputs/bams/TOY_PAIR.tumor.bam` are supported. Absolute paths to existing BAMs also work, without copying them.

The TSV template has nine headers: `patient_id`, `tumor_bam`, `tumor_bai`, `normal_bam`, `normal_bai`, `sex`, `tumor_phenotype`, `purity`, `ploidy`. Fill the first five columns and keep all headers; metadata values can be blank. **Missing purity or ploidy does not block calling.**

In the YAML, set your dataset label and HPC account/partition, and choose callers/concurrency. The template selects all five callers and **eight concurrent pairs per caller**. It already uses the standard samplesheet/reference paths and inherits downloaded image filenames. Change a path only when using a different location. [config/defaults.yaml](config/defaults.yaml) defines the resource paths and scientific settings.

**Optional: check your setup before submission:**

~~~bash
python scripts/workflow.py callers --config config/private-cohort.yaml --check-files
~~~

This checks the selected pairs' BAMs/indexes, reference indexes, selected caller images, optional BEDs and BRASS bundle files. It reports all missing, unreadable or empty files, without running Nextflow or submitting jobs. A normal launch performs the same check automatically. Unselected callers' resources are not required.

### 2. Launch the full cohort

Activate your environment, load site modules, and submit from the repository root:

~~~bash
mkdir -p logs
sbatch --account=YOUR_ACCOUNT --partition=YOUR_PARTITION \
  --cpus-per-task=1 --mem=4G --time=7-00:00:00 --output=logs/cohort-%j.log \
  scripts/run.slurm.sh callers --config config/private-cohort.yaml
~~~

**No separate preflight or pilot submission is needed.** Input checks run automatically before each pair's callers start. The small driver schedules tasks on compute nodes; your scheduler assigns available resources. Set driver duration to cover your cohort within your site's walltime and maintenance limits.

Eight pairs per caller permits up to 40 caller tasks. The default submission window is 24 outstanding tasks; use a [site config](docs/hpc.md#capacity-and-concurrency) to raise it. Adjust `max_parallel_callers` to your available capacity. For PBS, LSF, SGE or an allocated local node, see [HPC instructions](docs/hpc.md).

### 3. Optional: preview or change settings

You can skip this section when using your existing settings. These commands help you check a configuration or plan a parameter experiment; they are separate from launching the cohort.

**Preview your configuration without running anything:**

~~~bash
python scripts/workflow.py callers --config config/private-cohort.yaml --dry-run
~~~

This prints the final settings and the Nextflow command. **`--dry-run` does not submit jobs or call SVs.**

**Preview different parameter values:**

~~~bash
python scripts/workflow.py callers --config config/private-cohort.yaml \
  --set delly.options.map_qual=2 --set max_parallel_callers=8 --dry-run
~~~

Here, `delly.options.map_qual=2` sets DELLY's minimum read mapping quality to 2, and `max_parallel_callers=8` allows up to eight pairs **per caller**, subject to scheduler limits. `--set` overrides a value for that invocation; it does not edit either YAML file. With `--dry-run`, this example is still only a preview.

To **run** that experiment, put the settings in your `config/private-cohort.yaml` and submit the command in Step 2. Add or edit these values inside the file's existing sections, keeping your other settings:

~~~yaml
outdir: results/mapq2_run01
max_parallel_callers: 8
delly:
  options:
    map_qual: 2
~~~

The new `outdir` keeps the experiment's results separate from your previous run. Edit your private YAML, not `config/defaults.yaml`.

### 4. Optional: restart an interrupted run

Use this after a run stops—for example, because of maintenance or a failed task. Fix the reported problem, then restart with the **same inputs, settings and output directory**. Nextflow reuses completed tasks whose cached inputs/settings still match and runs the remaining tasks.

From the repository root, find the stopped run's Nextflow session ID:

~~~bash
nextflow log
~~~

Copy its value from the `SESSION ID` column. This is a **Nextflow ID, not a Slurm job ID**. Resubmit your original command with `--resume` added:

~~~bash
sbatch --account=YOUR_ACCOUNT --partition=YOUR_PARTITION \
  --cpus-per-task=1 --mem=4G --time=7-00:00:00 --output=logs/resume-%j.log \
  scripts/run.slurm.sh callers --config config/private-cohort.yaml --resume SESSION_UUID
~~~

Replace `SESSION_UUID` with the ID you copied. Include any `--set` or `--site-config` options from the original run. Keep both `work/` (task files) and `.nextflow/` (execution history/cache) in the checkout. [Resume details](https://www.nextflow.io/docs/latest/cache-and-resume.html).

## Exome and Targeted Sequencing

This section is for **whole-exome sequencing (WES)** or a **targeted gene panel**, rather than whole-genome sequencing (WGS). Use the same nine-column matched tumour–normal samplesheet. The pipeline calls SVs from these BAMs using **Manta**; the other four caller adapters currently use WGS settings.

**1. Make a configuration for your exome or panel cohort.** Copy your already-filled cohort YAML so it retains the reference, image paths and HPC settings:

~~~bash
cp config/private-cohort.yaml config/private-capture.yaml
~~~

Update `input` to the exome/panel samplesheet, set a suitable `dataset_id`, and replace the following keys in that new file. For WES:

~~~yaml
sequencing_type: WES
callers: manta
target_bed: null  # Optional: replace with your capture-kit BED to restrict call regions.
outdir: results/exome_run01
~~~

For a targeted panel, use these values instead:

~~~yaml
sequencing_type: targeted
callers: manta
target_bed: inputs/regions/panel_targets.bed
outdir: results/panel_run01
~~~

`target_bed` lists the genomic intervals covered by your capture kit or panel, usually supplied by its manufacturer. It is optional for WES and required for targeted panels. Use BED coordinates: start is 0-based, end is excluded (for example, `1 100 200` covers 100 bases). Chromosome names must match the BAM/reference—for example, `1` and `chr1` are different names.

**2. Submit the new configuration:**

~~~bash
sbatch --account=YOUR_ACCOUNT --partition=YOUR_PARTITION \
  --cpus-per-task=1 --mem=4G --time=7-00:00:00 --output=logs/capture-%j.log \
  scripts/run.slurm.sh callers --config config/private-capture.yaml
~~~

For either assay, the pipeline automatically sets Manta's `--exome` option; you do not need to add native Manta flags. When a target BED is supplied, it validates/merges the intervals, subtracts any configured exclusions, and prepares Manta's indexed call-region file. Have `bgzip` and `tabix` on the compute-node PATH for this preparation.

You receive the same Manta VCF, QC and benchmark files described below. **This launches SV calling, not signature analysis.** The current COSMIC SV32 signature route accepts WGS callsets, not WES/panel callsets. [Assay support](docs/workflow-status.md#assay-support) · [Alternative two-file capture configuration](examples/capture.yaml).

## SV Signature Analysis

Copy [examples/callsets.tsv](examples/callsets.tsv) to `inputs/callsets.tsv`, fill it with your completed caller VCFs, then build the hashed manifest:

~~~bash
mkdir -p inputs manifests
cp examples/callsets.tsv inputs/callsets.tsv
python scripts/build_sv_manifest.py --inventory inputs/callsets.tsv \
  --output manifests/sv_manifest.tsv
cp examples/signatures.yaml config/private-signatures.yaml
~~~

The manifest and Python environment paths already match the standard layout. Set your HPC values, complete the one-time [signature setup](docs/signatures.md), and run COSMIC SV32. You can check its files with `python scripts/workflow.py signatures --config config/private-signatures.yaml --check-files`.

| Route | What it does | Main outputs |
| --- | --- | --- |
| `cosmic_sv32` | SV32 counts → COSMIC SV1–SV12 fit | Count matrix and sample exposures |

~~~bash
# COSMIC reference fitting
sbatch --account=YOUR_ACCOUNT --partition=YOUR_PARTITION \
  --cpus-per-task=1 --mem=4G --time=24:00:00 --output=logs/cosmic-%j.log \
  scripts/run.slurm.sh signatures --config config/private-signatures.yaml --set method=cosmic_sv32
~~~

Groups run independently by cohort and caller. The PCAWG/ClusterSV route is paused; selecting `method=pcawg_clustersv` reports the hold before submitting analysis tasks. Existing code and outputs remain preserved. [Signature guide and status](docs/signatures.md).

## Retrieve Your Results

Under your configured `outdir`:

| Result | Location |
| --- | --- |
| Input checks | `preflight/PAIR/` |
| Manta somatic calls | `manta/PAIR/manta/results/variants/somaticSV.vcf.gz` |
| DELLY somatic calls | `delly/PAIR/output/somaticSV.vcf.gz` |
| SvABA somatic calls | `svaba/PAIR/output/PAIR.svaba.somatic.sv.vcf` |
| GRIDSS somatic calls | `gridss/PAIR/output/somatic/somatic.high_confidence.vcf.gz`, `somatic.full.vcf.gz` |
| BRASS annotated calls | `brass/PAIR/output/callset/calls/*.annot.vcf.gz` |
| Runtime / memory | Manta: `manta/PAIR/benchmark.json`; others: `CALLER/PAIR/output/benchmark.json` |
| COSMIC counts / exposures | `COHORT/CALLER/group_output/SV32.matrix.tsv`, `SV32.exposures.tsv` |
| Settings / logs | `logs/WORKFLOW-TIMESTAMP/` |

Native outputs and QC remain alongside these files. GRIDSS's separate `joint.vcf.gz` contains joint calls before somatic filtering. [Output formats](docs/inputs-outputs.md) · [Performance exports](docs/performance.md).

## Test Your Installation

Optional setup checks:

~~~bash
python scripts/test_installation.py
# Workstation or allocated compute node: three synthetic SV events.
python scripts/workflow.py test --profile test
# Actual main.nf preflight + caller stubs + cached resume.
python scripts/test_main_stub.py --out results/main_stub_check01
~~~

CI checks installation, caller wiring, cached resume and publication. PCAWG/ClusterSV numerical workflows are paused. See [tests and scientific controls](docs/validation.md).

## Further Documentation

[Parameters](CONFIG_PARAMETERS.md) · [Setup](docs/setup.md) · [HPC](docs/hpc.md) · [Signatures](docs/signatures.md) · [Method scope](docs/workflow-status.md) · [Development](docs/development.md)

## License

Independent code is [MIT-licensed](LICENSE); upstream tools retain their own [licenses](THIRD_PARTY.md). Public examples are synthetic. Patient data, private configs and runtime outputs stay outside Git.
