# SV Workflows

**From matched BAMs to structural variants and SV signatures—one framework, any cohort.**

[![Tests](https://github.com/jierui-cell/sv-workflows/actions/workflows/tests.yml/badge.svg?branch=main)](https://github.com/jierui-cell/sv-workflows/actions/workflows/tests.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Nextflow](https://img.shields.io/badge/Nextflow-25.10.7-blue)](https://www.nextflow.io/)
[![Python](https://img.shields.io/badge/Python-3.11-blue)](requirements.txt)

[Launch a cohort](#usage-instructions) · [Signatures](#sv-signature-analysis) · [Exome / panels](#exome-and-targeted-sequencing) · [Results](#retrieve-your-results) · [Parameters](CONFIG_PARAMETERS.md)

## Overview

Run **Manta, DELLY, GRIDSS, SvABA and BRASS** across a matched tumour–normal cohort. Nextflow handles scheduling, automatic input checks, parallel execution and resume. YAML gives you one place to select callers, tune parameters and size resources. Each run saves its settings, commands, QC and comparable runtime/memory records.

Choose **COSMIC SV32 reference fitting**, or **ClusterSV + de novo NMF** to discover cohort-specific SV patterns without copy-number inputs.

~~~text
Paired BAMs → Automatic checks → Five SV callers → VCFs + QC + benchmarks
                                                       │
                                  ┌────────────────────┴─────────────────┐
                                  COSMIC SV32                            ClusterSV
                                  32-channel counts                      Clusters + footprints
                                  → SV1–SV12 exposures                   → CN-free 38-channel counts
                                                                         → De novo NMF signatures
~~~

The framework is for research analysis. [Method scope](docs/workflow-status.md) describes reference builds, assay support and output interpretation.

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

Use the FASTA used to align your BAMs. Reuse existing indexes or create them on a compute node with samtools and BWA:

~~~bash
samtools faidx /path/to/reference.fa
bwa index /path/to/reference.fa
samtools dict /path/to/reference.fa > /path/to/reference.dict
~~~

BRASS also needs a GRCh37 resource bundle: [prepare or reuse a bundle](docs/setup.md#brass-bundle). Enable the callers you have prepared; [setup](docs/setup.md) lists each caller's resources.

## Usage Instructions

### 1. Supply your samples and settings

~~~bash
cp examples/samplesheet.tsv /path/to/private/paired_samples.tsv
cp examples/cohort.yaml config/private-cohort.yaml
~~~

The TSV template has nine headers: `patient_id`, `tumor_bam`, `tumor_bai`, `normal_bam`, `normal_bai`, `sex`, `tumor_phenotype`, `purity`, `ploidy`. Fill the first five columns with a unique pair ID and absolute BAM/index paths. Keep all headers; metadata values can be blank. **Missing purity or ploidy does not block calling.**

In the YAML, set input/reference paths, dataset label, output directory and HPC account/partition. The template selects all five callers and **eight concurrent pairs per caller**. Remove callers you do not need. Scientific defaults are already defined in [config/defaults.yaml](config/defaults.yaml).

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

### 3. Tune or resume

Inspect or override settings without editing defaults:

~~~bash
python scripts/workflow.py callers --config config/private-cohort.yaml \
  --set delly.options.map_qual=2 --set max_parallel_callers=8 --dry-run
~~~

Resume an interrupted run by adding `--resume SESSION_UUID` to the original command. Keep its `work/` directory. For a parameter experiment, choose a new `outdir`.

## Exome and Targeted Sequencing

Use the same samplesheet, launcher, QC and benchmark outputs. **The implemented capture-assay route is Manta**, with native `--exome` settings for both WES and targeted sequencing.

~~~yaml
sequencing_type: WES  # WGS (default), WES or targeted
callers: manta
target_bed: /absolute/path/to/capture_targets.bed
outdir: results/exome_run01
~~~

`target_bed` is optional for WES and required for targeted panels. Use 0-based, half-open BED coordinates and reference-matched contig names. The workflow normalizes targets, applies optional exclusions, and builds Manta's indexed inclusion list; `bgzip` and `tabix` are needed on compute nodes when a BED is supplied. [Capture overlay](examples/capture.yaml).

The other four adapters currently use WGS settings. Genome-wide signature analysis is for WGS callsets; capture assays need a capture-aware signature model. [Assay support](docs/workflow-status.md#assay-support).

## SV Signature Analysis

Create a five-column VCF inventory following [examples/callsets.tsv](examples/callsets.tsv), then build the hashed manifest:

~~~bash
mkdir -p manifests
python scripts/build_sv_manifest.py --inventory /path/to/private/callsets.tsv \
  --output manifests/sv_manifest.tsv
cp examples/signatures.yaml config/private-signatures.yaml
~~~

Set the manifest and HPC values, complete the one-time [signature setup](docs/signatures.md), and choose a route:

| Route | What it does | Main outputs |
| --- | --- | --- |
| `cosmic_sv32` | SV32 counts → COSMIC SV1–SV12 fit | Count matrix and sample exposures |
| `pcawg_clustersv` | ClusterSV → CN-free features → de novo NMF | Clusters, footprints, 38-channel counts, signatures and activities; optional timing |

~~~bash
# COSMIC reference fitting
sbatch --account=YOUR_ACCOUNT --partition=YOUR_PARTITION \
  --cpus-per-task=1 --mem=4G --time=24:00:00 --output=logs/cosmic-%j.log \
  scripts/run.slurm.sh signatures --config config/private-signatures.yaml --set method=cosmic_sv32

# ClusterSV + de novo NMF
sbatch --account=YOUR_ACCOUNT --partition=YOUR_PARTITION \
  --cpus-per-task=1 --mem=4G --time=24:00:00 --output=logs/clustersv-%j.log \
  scripts/run.slurm.sh signatures --config config/private-signatures.yaml \
  --set method=pcawg_clustersv --set outdir=results/clustersv_run01
~~~

Groups run independently by cohort and caller. NMF selects a rank from repeated seeded fits; `--set pcawg_rank=3` selects three components. Use `pcawg_stage=features` or `clustering` to stop earlier. Timing is optional and requires no setup for the default de novo route. The CN-free 38-channel scheme and `cSV1…` labels describe cohort-specific patterns, not the published PCAWG 44-feature matrix or reference signatures. [Feature definitions and NMF](docs/signatures.md#cn-free-features-and-de-novo-nmf).

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
| ClusterSV counts | `COHORT/CALLER/group_output/ClusterSV38.matrix.tsv` |
| De novo signatures / activities | `COHORT/CALLER/group_output/ClusterSV.NMF.signatures.tsv`, `ClusterSV.NMF.activities.tsv` |
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

CI runs these checks plus ClusterSV and NMF controls. See [tests and scientific controls](docs/validation.md).

## Further Documentation

[Parameters](CONFIG_PARAMETERS.md) · [Setup](docs/setup.md) · [HPC](docs/hpc.md) · [Signatures](docs/signatures.md) · [Method scope](docs/workflow-status.md) · [Development](docs/development.md)

## License

Independent code is [MIT-licensed](LICENSE); upstream tools retain their own [licenses](THIRD_PARTY.md). Public examples are synthetic. Patient data, private configs and runtime outputs stay outside Git.
