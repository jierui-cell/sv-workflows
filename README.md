# SV Workflows

Reproducible paired tumour–normal structural-variant calling, followed by selectable SV signature analysis. Built from workflows used on melanoma and pancreatic neuroendocrine tumour WGS; public examples contain synthetic data only.

This is a development release for research, not a clinically validated diagnostic pipeline. The production callsets are unchanged; the portable refactor needs a real-pair equivalence pilot before cohort deployment.

## Workflows and readiness

| Workflow | What it produces | Validation boundary |
| --- | --- | --- |
| Manta, DELLY, SvABA | Paired somatic SV VCFs and QC | Original workflows completed both cohorts; refactor equivalence pilot required |
| GRIDSS | Joint VCF plus native somatic high-confidence/full VCFs | No panel of normals or BSgenome; documented filtering limitation |
| BRASS | Annotated paired SV VCF/BEDPE | GRCh37, validated reference bundle required; missing purity/ploidy use explicitly recorded native defaults |
| Cue | Tumour/normal predictions and provisional paired candidates | Research-only, not established somatic truth |
| LUMPY | Historical adapter/testing code only | Discontinued; launcher rejects it |
| JaBbA | Copy-number/junction graph | Licensed solver required; toy solver passed, real pilot did not converge |
| COSMIC SV32 | 32-channel matrix and SV1–SV12 refitting | Pinned SigProfiler clustering control required |
| ClusterSV | CN-free clusters, footprints and timing annotations | Original clustering method; not the full mechanistic PCAWG classifier |
| Advanced 44-feature analysis | Matrix from supplied validated classifications | Separate newer-paper feature profile; cannot derive missing CN-dependent classes from clusters |

## Start here

Requirements: Linux, Python 3.11, Java 17+, Nextflow 25.10.7. Real caller tasks also need Apptainer or Singularity, samtools and bcftools. No GPU is required by the supported CPU workflows. Caller images and scientific source revisions are pinned; licensed software and patient data are not distributed.

```bash
git clone https://github.com/jierui-cell/sv-workflows.git
cd sv-workflows
python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python scripts/test_installation.py
python scripts/workflow.py test --profile test
```

The toy test checks schema/parameter validation, VCF-to-BEDPE conversion, resource-recording behavior and a real Nextflow task. It does not run a caller or measure SV sensitivity. See [validation](docs/validation.md) for container checks and biological controls.

## Configure once, reuse across cohorts

Versioned defaults are in [config/defaults.yaml](config/defaults.yaml), [config/signatures.yaml](config/signatures.yaml) and [config/replication.yaml](config/replication.yaml). Keep your paths and scheduler settings in a private YAML overlay. Precedence is defaults → repeated overlays in order → `--set dotted.key=value`. Unknown/duplicate YAML keys fail early. Never store private overlays or genomic data in Git.

```yaml
# config/private-cohort.yaml (ignored by Git)
input: /absolute/path/to/paired_samples.tsv
reference: /absolute/path/to/alignment_reference.fa
dataset_id: PRIVATE_COHORT
outdir: /absolute/path/to/results/run01
callers: manta,delly,svaba,gridss,brass
slurm_account: YOUR_ACCOUNT
slurm_partition: YOUR_PARTITION
manta_image: /absolute/path/to/manta-1.6.0--py27h9948957_6.sif
delly:
  image: /absolute/path/to/delly-2.6.0--h3752d28_0.sif
svaba:
  image: /absolute/path/to/svaba-1.2.0--h69ac913_1.sif
gridss:
  image: /absolute/path/to/gridss-2.13.2--h96c455f_6.sif
brass:
  image: /absolute/path/to/brass-v6.3.4.sif
  bundle: /absolute/path/to/validated-brass-bundle
```

Check the resolved configuration before submitting. Keep modules configurable with `modules_samtools`, `modules_bcftools`, `modules_calling`, or preload them in your site environment.

```bash
python scripts/workflow.py callers --config config/private-cohort.yaml --dry-run
mkdir -p logs
sbatch --account=YOUR_ACCOUNT --partition=YOUR_PARTITION --cpus-per-task=1 \
  --mem=4G --time=48:00:00 --output=logs/driver-%j.log \
  scripts/run.slurm.sh callers --config config/private-cohort.yaml --set callers=preflight
# After preflight and an equivalence pilot, use a NEW output directory:
sbatch --account=YOUR_ACCOUNT --partition=YOUR_PARTITION --cpus-per-task=1 \
  --mem=4G --time=7-00:00:00 --output=logs/driver-%j.log \
  scripts/run.slurm.sh callers --config config/private-cohort.yaml \
  --set max_parallel_callers=8 --set manta_options.minPassSomaticScore=30
```

The concurrency limit is per caller, not global: five callers × eight pairs can submit 40 tasks, additionally bounded by Nextflow's executor queue. Allocate within your account limits and measured disk capacity. Driver walltime must exceed the entire workflow duration.

## Signatures from completed VCFs

Both methods take the same six-column manifest with VCF checksums. There is no melanoma-specific path resolver.

```bash
python scripts/build_sv_manifest.py --inventory examples/callsets.tsv --output /tmp/toy-sv-manifest.tsv
# Use a private inventory for real data; preserve the recorded QC status.
python scripts/workflow.py signatures --config config/private-signatures.yaml \
  --set method=cosmic_sv32 --dry-run
python scripts/workflow.py signatures --config config/private-signatures.yaml \
  --set method=pcawg_clustersv --set pcawg_stage=clustering --dry-run
```

See [signature setup and scientific scope](docs/signatures.md), including the replication-timing tertile reference. Decomposition against the original PCAWG 16-signature matrix is deferred: that numerical matrix is not available here. No unavailable feature is fabricated or filled with zero.

## Documentation

- [Installation and reference resources](docs/setup.md)
- [HPC portability and storage](docs/hpc.md)
- [Input/output schemas and provenance](docs/inputs-outputs.md)
- [Scientific validation and toy tests](docs/validation.md)
- [Signature methods](docs/signatures.md)
- [Performance comparison](docs/performance.md)
- [Development and publishing](docs/development.md)

The same documentation can be published with GitHub Pages. No identifiable cohort metadata, real VCFs/BAMs, reports, model weights or solver installation files are included.

## License

Our independent workflow code is [MIT-licensed](LICENSE). Upstream callers, scientific sources and external resources remain subject to their original terms; see [third-party notices](THIRD_PARTY.md).
