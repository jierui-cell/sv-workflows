# Installation and Resources

Prepare the tools and shared resources once, then reuse them across cohorts.

## Standard Locations

Use `inputs/samplesheet.tsv` for pair metadata, `inputs/bams/` for BAMs/indexes, `resources/reference/reference.fa` for the alignment reference and `resources/brass/bundle/` for prepared BRASS resources. Images go in `containers/images/`. These paths are already in the YAML defaults; relative paths are anchored to the checkout root. [Folder layout](../README.md#standard-folder-layout).

Existing BAM/reference locations can be used through absolute YAML/TSV overrides or suitable symlinks. Keep linked references with their original index files as well as the configured names; the checker reports missing source indexes where an adapter follows the FASTA link.

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

Use `--engine singularity` if appropriate. Images are pinned by OCI digest in `config/containers.tsv`; the downloader records a SIF SHA-256 receipt. It reuses matching downloads. Its output filenames already match `config/defaults.yaml`, so no image-path edits are needed.

## Alignment Reference

Place the FASTA used for alignment at `resources/reference/reference.fa`, plus its indexes:

| Resource | Used by |
| --- | --- |
| `resources/reference/reference.fa` + `reference.fa.fai` | All callers |
| `reference.fa.amb`, `.ann`, `.bwt`, `.pac`, `.sa` alongside the FASTA | SvABA, GRIDSS |
| `resources/reference/reference.dict` | GRIDSS and BRASS bundle preparation |

Reuse existing indexes, or build them on a compute node:

~~~bash
samtools faidx resources/reference/reference.fa
bwa index resources/reference/reference.fa
samtools dict resources/reference/reference.fa > resources/reference/reference.dict
~~~

Automatic preflight checks BAM/reference dictionaries, the declared build, sample names, sortedness, quickcheck and indexed access.

## BRASS Bundle

BRASS uses a GRCh37 bundle containing depth/GC masks, centromere/telomere and cytoband annotations, a VAGrENT transcript cache and microbial/viral resources. The bundle adapter matches sequence identities, adapts chromosome labels and checks the native readers. Coordinates remain unchanged.

Resource archive identities are listed in `config/brass-resources.tsv`. Collect the extracted publisher files in `resources/brass/source/`, with the publisher reference's sequence dictionary at `resources/brass/source/reference.dict`. That dictionary must include the original sequence identities (MD5s), which the adapter compares with your alignment reference. Run preparation in a compute allocation:

~~~bash
python scripts/prepare_brass_bundle.py \
  --source resources/brass/source \
  --source-dict resources/brass/source/reference.dict \
  --reference resources/reference/reference.fa \
  --image containers/images/brass-v6.3.4.sif \
  --output resources/brass/bundle --engine apptainer
~~~

Use `--engine singularity` on Singularity sites. For resources already adapted by this framework, `--reuse-adapted` accepts the existing directory as `--source`, without `--source-dict`. It checks hashes, reruns native readers and creates a new validated bundle; keep `--output resources/brass/bundle` to use the default. It does not overwrite an existing output directory.

The default `brass.bundle` already points to that output. The current bundle route uses chr-prefixed GRCh37 resources. Native validation creates `COMPATIBILITY_VALIDATED`, `native-validation.json`, `bundle-receipt.json` and the runtime overlay along with the resource files. The setup checker verifies their presence/status before BRASS task submission.

## Optional Setup Checks

~~~bash
python scripts/workflow.py callers --config config/private-cohort.yaml --check-files
python scripts/test_installation.py
# Workstation / allocated compute node:
python scripts/workflow.py test --profile test
python scripts/test_main_stub.py --out results/main_stub_check01
~~~

`--check-files` is a read-only inventory check and works without Nextflow or caller execution. Real launches perform it automatically and retain `setup-check.json` in their logs. The toy test produces three events; the entry-point test exercises synthetic BAM preflight and cached caller stubs. [Test details](validation.md).
