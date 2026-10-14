# Installation and resources

Use Linux, Python 3.11 and Java 17 or newer. Load site modules if appropriate; versions/module names are site-specific. Install Nextflow **25.10.7**, then run `nextflow -version`. The launcher requests this version via `NXF_VER`; an unrelated system installation still needs checking.

```bash
python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python scripts/test_installation.py
# HPC: run the following toy task only inside an allocated compute session.
python scripts/workflow.py test --profile test
```

For a batch allocation, use the [Slurm toy-test example](../README.md#test-your-installation). Container downloads and WGS reference indexing are not needed for this synthetic test.

Real tasks require host `python3`, `samtools`, `bcftools` and the selected SIF runtime. Generating BWA reference indexes requires `bwa`; load/install it before the indexing commands below. BRASS additionally uses `jq`, `patch`, `realpath`, `sha256sum` and standard Linux utilities. Put the activated environment and modules in the compute-node environment, not only your login shell. Docker execution is not implemented: these adapters explicitly use the Singularity/Apptainer CLI.

## Caller containers

`config/containers.tsv` contains immutable OCI digests, not mutable image tags. Download once on a permitted network/compute node:

```bash
python scripts/fetch_containers.py --engine apptainer --tools manta,delly,svaba,gridss,brass
```

Each image gets a SHA-256 receipt. Existing images without a matching receipt require review; the helper never silently replaces them. Set the resulting absolute SIF paths in your private overlay. Container pulling on another machine may produce a different SIF byte hash from the same pinned OCI image; record both identities.

## Alignment reference

Use the **exact FASTA used to align the BAMs**, not merely a reference with the same assembly name. The reconstructed GRCh37_ICGC_standard_v2 reference was accepted for the original cohorts; it is not automatically valid for another dataset. Provide `.fai`; SvABA/GRIDSS also need `.amb`, `.ann`, `.bwt`, `.pac`, `.sa`; GRIDSS needs a sequence dictionary. Prepare indexes on a compute node, retaining the FASTA/index hashes:

```bash
samtools faidx /path/to/reference.fa
bwa index /path/to/reference.fa
samtools dict /path/to/reference.fa > /path/to/reference.dict
```

Preflight checks header contigs/lengths, sortedness metadata, BAM quickcheck, sample names, indexes and an indexed query. These checks do not prove donor identity or every alignment's reference compatibility.

## BRASS bundle

BRASS requires more than a FASTA: depth masks, GC windows, centromere/telomere and cytoband annotations, a serialized VAGrENT transcript cache, and microbial/viral resources. The current adapter is limited to a **validated chr-prefixed GRCh37 bundle**; arbitrary GRCh38 support is not claimed.

`config/brass-resources.tsv` records publisher archive MD5s; `scripts/adapt_brass_resources.py` checks source/target sequence dictionaries including sequence MD5 before renaming contigs. It does not lift coordinates. `scripts/prepare_brass_bundle.py` performs adaptation, native reader validation, the guarded chromosome-name overlay and a reference-bound receipt. Use it on an allocated node with downloaded/extracted publisher resources; see its `--help`. The marker is created only after native validation passes. Do not create the marker manually. Resource downloading/extraction remains a documented prerequisite, not an automatic cohort-specific hack.

If you already have resources produced by this adapter, `--reuse-adapted --source /path/to/existing/bundle` can create a **new** validated bundle without repeating adaptation. It verifies the prior adaptation's target-dictionary hash and every recorded adapted-resource hash, reruns native readers, and regenerates the image-derived overlay and reference/image-bound receipt. It does not modify the original bundle or accept its marker alone. See `--help` for the required reference, image and new output path; `--source-dict` is only required for fresh publisher adaptation.
