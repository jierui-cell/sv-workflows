# Installation and Resources

Follow the [README setup](../README.md#requirements) once; reuse the tools and resources across cohorts.

## Standard Locations

Inputs go in `inputs/`, the shared alignment FASTA in `resources/reference/`, caller images in `containers/images/` and BRASS resources in `resources/brass/bundle/`. Paths can be changed in your cohort YAML.

## Launcher and Compute Tools

Put your site modules in `config/private-cohort.yaml`:

~~~yaml
modules: [java/20.0.1, samtools/1.20, bcftools/1.20]
container_engine: singularity
~~~

Use module names from your HPC, adding Python or container-runtime modules where needed; `modules: []` uses installed tools. The launcher loads this list before Nextflow starts, and Nextflow loads it inside compute tasks. Slurm submissions select `.venv` automatically.

Install Nextflow in `.venv/bin/nextflow` using the [README command](../README.md#1-install-once), then initialize it with Java from your YAML:

~~~bash
python scripts/hpc_environment.py -- .venv/bin/nextflow -version
~~~

This first start needs internet access. The same helper runs other setup commands with your configured modules.

## Caller Containers

~~~bash
python scripts/fetch_containers.py
~~~

The downloader reads `config/private-cohort.yaml`, loads its modules and uses its container engine. Add `--tools manta,delly` to select images or `--config` for another YAML. Image digests are pinned in [config/containers.tsv](../config/containers.tsv).

## Alignment Reference

Place the alignment FASTA and matching indexes at:

| Files | Used by |
| --- | --- |
| `resources/reference/reference.fa` and `.fa.fai` | All callers |
| `reference.fa.amb`, `.ann`, `.bwt`, `.pac`, `.sa` | SvABA and GRIDSS |
| `resources/reference/reference.dict` | GRIDSS and BRASS preparation |

Existing files can be linked here. To build indexes on a compute node, include BWA in your modules and run:

~~~bash
python scripts/hpc_environment.py -- samtools faidx resources/reference/reference.fa
python scripts/hpc_environment.py -- bwa index resources/reference/reference.fa
python scripts/hpc_environment.py -- samtools dict resources/reference/reference.fa > resources/reference/reference.dict
~~~

## BRASS Bundle

BRASS uses the shared human FASTA plus additional masks and annotations. Prepare them on a compute node:

~~~bash
python scripts/setup_brass.py
~~~

It reads the cohort YAML, loads its modules and prepares `resources/brass/bundle/`. Use `--config` for another YAML or `--engine singularity` to select that runtime explicitly.

| Sanger archive (~3.1 GB total) | Contents |
| --- | --- |
| [CNV_SV_ref_GRCh37d5_brass6+.tar.gz](https://ftp.sanger.ac.uk/pub/cancer/dockstore/human/CNV_SV_ref_GRCh37d5_brass6+.tar.gz) | Depth/GC masks, chromosome annotations and microbial/viral sequences |
| [VAGrENT_ref_GRCh37d5_ensembl_75.tar.gz](https://ftp.sanger.ac.uk/pub/cancer/dockstore/human/VAGrENT_ref_GRCh37d5_ensembl_75.tar.gz) | Transcript annotations and sequences |

### Why preparation is needed

The [BRASS reference guide](https://github.com/cancerit/BRASS/wiki) and [v6.3.4 inputs](https://github.com/cancerit/BRASS/blob/v6.3.4/perl/bin/brass.pl#L363-L380) specify these resources separately from the human FASTA. Sanger describes [bundles for different naming conventions](https://github.com/cancerit/dockstore-cgpwgs/wiki/Reference-archives#generating-alternative-reference-sets).

Our adapter verifies sequence identities against the included GRCh37d5 dictionary, changes labels such as `1` to `chr1` and rebuilds indexes. Coordinates stay unchanged. Other builds require changes to [setup](../scripts/setup_brass.py), [preparation](../scripts/prepare_brass_bundle.py), [adaptation](../scripts/adapt_brass_resources.py) and [validation](../scripts/validate_brass_cache.pl), with build-matched resources.

Run setup again to reuse a completed bundle or resume downloads. Keep `resources/brass/source/`, which supplies linked files. Logs are in `logs/brass-setup-TIMESTAMP/`. Reuse an existing prepared bundle for the same reference by setting `brass.bundle` in YAML.

## Optional Setup Checks

~~~bash
python scripts/workflow.py callers --config config/private-cohort.yaml --check-files
~~~

File checks also run at launch. [Installation tests](validation.md).
