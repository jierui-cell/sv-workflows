# SV Signature Analysis

Build SV32 counts from completed caller VCFs and fit COSMIC SV1–SV12 reference signatures. Results are grouped by cohort and caller.

## Use Your Caller Results

Copy the signature template and set your account/partition:

~~~bash
cp examples/signatures.yaml config/private-signatures.yaml
~~~

`caller_config: config/private-cohort.yaml` automatically supplies the completed callsets and HPC modules. No VCF list is needed. For a caller run launched with command-line overrides, use its saved `logs/callers-TIMESTAMP/resolved.yaml` instead.

| Caller | Selected output |
| --- | --- |
| Manta / DELLY | `somaticSV.vcf.gz` |
| SvABA | `PAIR.svaba.somatic.sv.vcf` |
| GRIDSS | `somatic.high_confidence.vcf.gz` |
| BRASS | `*.annot.vcf.gz` |

The launcher checks completion records and selects GRCh37 WGS callsets. Use `groups` and `pair_ids` for subsets. It saves the hashed input manifest and selection record in `logs/signatures-TIMESTAMP/`.

## Install Once

Create the pinned Python environment using the [README setup commands](../README.md#sv-signature-analysis). The workflow selects it automatically; no separate activation is needed.

## Launch COSMIC SV32

From the repository root:

~~~bash
mkdir -p logs
sbatch --account=YOUR_ACCOUNT --partition=YOUR_PARTITION \
  --cpus-per-task=1 --mem=4G --time=24:00:00 --output=logs/cosmic-%j.log \
  scripts/run.slurm.sh signatures --config config/private-signatures.yaml
~~~

~~~text
Somatic VCFs → BEDPE → SV clustering → 32-channel counts → COSMIC SV1–SV12 fit
~~~

The workflow applies caller-specific filters, retains intrachromosomal events ≥1 kb, uses SigProfiler clustering and fits COSMIC v3.6. Input checks and the numerical clustering control run automatically. Set `skip_fit: true` for counts only.

## Retrieve Results

Find results under `OUTDIR/COHORT/CALLER/group_output/`:

| File | Contents |
| --- | --- |
| `SV32.matrix.tsv` | Feature counts per sample |
| `SV32.exposures.tsv` | Reference-signature exposures |
| Conversion/source audits | Retained events and input hashes |
| Reconstruction QC | Fit quality and residuals |

[All signature settings](../CONFIG_PARAMETERS.md#signature-settings).

## External VCFs (Optional)

For callsets produced elsewhere, fill the [inventory template](../examples/callsets.tsv) and build a manifest:

~~~bash
mkdir -p inputs manifests
cp examples/callsets.tsv inputs/callsets.tsv
# Fill inputs/callsets.tsv with your VCF paths, then:
python scripts/build_sv_manifest.py --inventory inputs/callsets.tsv \
  --output manifests/sv_manifest.tsv
~~~

Set `manifest: manifests/sv_manifest.tsv` in the signature YAML. Leave it `null` to use the previous caller run automatically.

## PCAWG Status

PCAWG/ClusterSV clustering, feature construction and NMF are on hold. Code and previous outputs are preserved; selecting `method=pcawg_clustersv` stops before submitting analysis tasks.
