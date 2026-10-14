# Installation and Resources

Prepare the tools and shared resources once, then reuse them across cohorts.

## Launcher and Compute Tools

Use Linux, Python 3.11, Java 17+ and Nextflow 25.10.7 through its [official launcher](https://www.nextflow.io/docs/latest/install.html) or site modules:

~~~bash
python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
~~~

Compute tasks use host `python3`, `samtools`, `bcftools` and Apptainer/Singularity. Reference indexing also uses `bwa`. BRASS uses `jq`, `patch` and standard Linux utilities. When supplying Manta region BEDs, also provide `bgzip` and `tabix`.

## Caller Containers

~~~bash
python scripts/fetch_containers.py --engine apptainer --tools manta,delly,gridss,svaba,brass
~~~

Use `--engine singularity` if appropriate. Images are pinned by OCI digest in `config/containers.tsv`; the downloader records a SIF SHA-256 receipt. It reuses matching downloads. Set image paths in your private YAML; the cohort template already contains the standard downloaded filenames.

## Alignment Reference

Use the FASTA used for alignment, plus its indexes:

| Resource | Used by |
| --- | --- |
| FASTA + `.fai` | All callers |
| BWA `.amb .ann .bwt .pac .sa` | SvABA, GRIDSS |
| `.dict` | GRIDSS and BRASS bundle preparation |

Reuse existing indexes, or build them on a compute node:

~~~bash
samtools faidx /path/to/reference.fa
bwa index /path/to/reference.fa
samtools dict /path/to/reference.fa > /path/to/reference.dict
~~~

Automatic preflight checks BAM/reference dictionaries, the declared build, sample names, sortedness, quickcheck and indexed access.

## BRASS Bundle

BRASS uses a GRCh37 bundle containing depth/GC masks, centromere/telomere and cytoband annotations, a VAGrENT transcript cache and microbial/viral resources. The bundle adapter matches sequence identities, adapts chromosome labels and checks the native readers. Coordinates remain unchanged.

Resource archive identities are listed in `config/brass-resources.tsv`. After downloading and extracting those publisher resources, run bundle preparation in a compute allocation:

~~~bash
python scripts/prepare_brass_bundle.py --help
~~~

For resources already adapted by this framework, use `--reuse-adapted` to prepare a new bundle without repeating adaptation. It checks recorded hashes, reruns native validation and creates a reference/image-bound receipt. Supply the reference, image and new output path shown by `--help`; fresh adaptation also takes `--source-dict`.

Set `brass.bundle` to the resulting directory. The current bundle route uses chr-prefixed GRCh37 resources. Enable BRASS once this shared resource is prepared; cohort execution then uses it automatically.

## Optional Setup Checks

~~~bash
python scripts/test_installation.py
# Workstation / allocated compute node:
python scripts/workflow.py test --profile test
python scripts/test_main_stub.py --out results/main_stub_check01
~~~

The toy test produces three events. The entry-point test exercises real synthetic BAM preflight and stubbed caller tasks with resume. [Test details](validation.md).
