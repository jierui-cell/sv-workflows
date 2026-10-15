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

**BRASS uses the same human FASTA as every other caller:** `reference` in your YAML. The additional `brass.bundle` folder holds annotations and filtering resources, not a different human assembly.

### One-command setup

Use the BRASS image installed in [Caller Containers](#caller-containers) and the shared FASTA, `.fai` and `.dict` prepared above. Run from the checkout root in a compute allocation with internet access, using Apptainer (or Singularity), `curl`, `patch` and the launcher's Python environment:

~~~bash
python scripts/setup_brass.py
~~~

With no explicit choice, setup uses the installed runtime: Apptainer first, otherwise Singularity. To select Singularity explicitly:

~~~bash
python scripts/setup_brass.py --engine singularity
~~~

If your cohort YAML changes resource paths or the container engine, use it here too:

~~~bash
python scripts/setup_brass.py --config config/private-cohort.yaml
~~~

Runtime selection follows `--engine` → `--set container_engine=...` → explicit YAML `container_engine` → automatic detection. A requested runtime is never silently replaced. To use Singularity for both setup and caller jobs, put `container_engine: singularity` in your private YAML and pass that YAML to both commands. The setup command does not change YAML defaults or your cohort configuration.

The command performs these steps automatically:

1. Checks the shared reference's sequence names, lengths and MD5 identities against the versioned GRCh37d5 dictionary, before downloading.
2. Downloads the two Sanger resource archives (~3.1 GB) into `resources/brass/downloads/`, verifies their pinned MD5 and SHA-256 checksums and resumes interrupted transfers.
3. Extracts only the BRASS and VAGrENT files used by this workflow into `resources/brass/source/`. It leaves unrelated ASCAT files and the optional normal panel unused.
4. Adapts chromosome labels, rebuilds resource indexes and checks the resources with BRASS/VAGrENT's native readers. **No coordinates are lifted over, BAMs realigned or human genomes replaced.**
5. Publishes the validated folder at `brass.bundle` (default `resources/brass/bundle/`) and records the outcome in `logs/brass-setup-TIMESTAMP/status.json` and `setup.log`.

The download URLs and checksums are versioned in [config/brass-resources.tsv](../config/brass-resources.tsv). These are Sanger's [published reference archives](https://github.com/cancerit/dockstore-cgpwgs/wiki/Reference-archives):

| Archive | Resources used by BRASS |
| --- | --- |
| [CNV_SV_ref_GRCh37d5_brass6+.tar.gz](https://ftp.sanger.ac.uk/pub/cancer/dockstore/human/CNV_SV_ref_GRCh37d5_brass6+.tar.gz) | High-depth masks, GC windows, centromere/telomere and cytoband annotations, viral/bacterial sequences |
| [VAGrENT_ref_GRCh37d5_ensembl_75.tar.gz](https://ftp.sanger.ac.uk/pub/cancer/dockstore/human/VAGrENT_ref_GRCh37d5_ensembl_75.tar.gz) | Transcript annotation cache and transcript sequences |

### Why preparation is needed

BRASS's [reference guide](https://github.com/cancerit/BRASS/wiki) and [v6.3.4 input specification](https://github.com/cancerit/BRASS/blob/v6.3.4/perl/bin/brass.pl#L363-L380) list these resources separately from the human FASTA. Sanger supplies them in the two archives above and describes [building bundles for different naming conventions](https://github.com/cancerit/dockstore-cgpwgs/wiki/Reference-archives#generating-alternative-reference-sets).

Our [adapter](../scripts/adapt_brass_resources.py) matches the archives to the shared reference: it verifies sequence identities, converts chromosome labels such as `1` to `chr1`, and rebuilds the resource indexes. This label adaptation is part of this workflow; genomic coordinates stay unchanged.

The small [GRCh37d5 dictionary](../config/brass-GRCh37d5.dict) is included as public reference metadata: chromosome names, lengths and sequence MD5s derived from the [1000 Genomes hs37d5 reference](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/). Its header records the source URL and compressed FASTA SHA-256. You do **not** need to download another human FASTA or find a separate dictionary yourself.

This route prepares **chr-prefixed GRCh37** resources (`chr1`–`chr22`, `chrX`, `chrY` for BRASS calling). For other reference builds, revise [setup_brass.py](../scripts/setup_brass.py) and its helpers ([bundle preparation](../scripts/prepare_brass_bundle.py), [resource adaptation](../scripts/adapt_brass_resources.py) and [native validation](../scripts/validate_brass_cache.pl)) to use build-matched archives, sequence dictionaries and chromosome names.

### Reuse and restart

Run the same setup command again to verify and reuse a completed bundle, without downloading or rebuilding it. Shared resources can be reused across cohorts with the same reference and pinned BRASS image. Keep `resources/brass/source/`: the bundle links unchanged transcript/microbial files there. After successful setup, the cached archives in `resources/brass/downloads/` are optional to retain.

Interrupted transfers retain `.partial` files. Completed extraction is cached with checksums, and the final bundle appears only after validation succeeds. If preparation fails, inspect the named step in `status.json` and its `setup.log`, then rerun the same command. Failed staged bundles remain under `.brass-bundle-*` beside the requested output for diagnosis; they are not used for calling. An existing completed bundle with changed inputs is reported as a mismatch: select a new `brass.bundle` directory for the new reference/image.

### Existing publisher files or an adapted bundle

If you already have the extracted publisher files, you can skip downloading. Put the `CNV_SV_ref/` and `VAGrENT_ref_GRCh37d5_ensembl_75/` directories under `resources/brass/source/`, then run the lower-level preparation command:

~~~bash
python scripts/prepare_brass_bundle.py \
  --source resources/brass/source \
  --source-dict config/brass-GRCh37d5.dict \
  --reference resources/reference/reference.fa \
  --image containers/images/brass-v6.3.4.sif \
  --output resources/brass/bundle --engine apptainer
~~~

Use `--engine singularity` on Singularity sites. For resources already adapted by this framework, `--reuse-adapted` accepts the existing directory as `--source`, without `--source-dict`. It checks hashes, reruns native readers and creates a new validated bundle. Use a new `--output` directory: the lower-level preparation command does not overwrite existing bundles.

Native validation creates `COMPATIBILITY_VALIDATED`, `native-validation.json`, `bundle-receipt.json` and the runtime overlay alongside the resources. The caller launcher's automatic setup check verifies required files and validation status before BRASS task submission. These checks are separate from the shared human reference check because they validate different resources.

## Optional Setup Checks

~~~bash
python scripts/workflow.py callers --config config/private-cohort.yaml --check-files
python scripts/test_installation.py
# Workstation / allocated compute node:
python scripts/workflow.py test --profile test
python scripts/test_main_stub.py --out results/main_stub_check01
~~~

`--check-files` is a read-only inventory check and works without Nextflow or caller execution. Real launches perform it automatically and retain `setup-check.json` in their logs. The toy test produces three events; the entry-point test exercises synthetic BAM preflight and cached caller stubs. [Test details](validation.md).
