# SV Signature Analysis

**COSMIC SV32 is the active signature-analysis route.** It builds a 32-channel matrix from caller-specific WGS VCFs and fits COSMIC SV1–SV12 reference signatures.

**PCAWG/ClusterSV signature analysis is on hold.** Clustering, feature construction, de novo NMF and related analysis are not being run. The implementation and existing outputs are retained for a later review. Selecting `method=pcawg_clustersv` reports the hold before submitting analysis tasks.

## Use Your Caller Results

After SV calling finishes, copy the signature settings:

~~~bash
cp examples/signatures.yaml config/private-signatures.yaml
~~~

Set your account/partition and signature output directory. The default `caller_config: config/private-cohort.yaml` points to the YAML used for SV calling. The launcher reads that file's output directory, cohort label, samplesheet and caller selection, selects the completed SV VCFs, and creates the hashed manifest automatically.

| Caller | VCF used for signatures |
| --- | --- |
| Manta | `somaticSV.vcf.gz` |
| DELLY | `somaticSV.vcf.gz` |
| SvABA | `PAIR.svaba.somatic.sv.vcf` |
| GRIDSS | `somatic.high_confidence.vcf.gz` |
| BRASS | `*.annot.vcf.gz` |

Both native and compact publishing layouts are supported. Completion records and technical output QC identify the finished callsets. Signature conversion retains the existing caller-specific event filters. The automatic route uses GRCh37 WGS caller results and checks that build/assay settings agree.

For a different caller YAML, change `caller_config` in your signature YAML. If the caller run used command-line overrides, point it to the saved `logs/callers-TIMESTAMP/resolved.yaml` to use that run's exact settings. Use `groups` and `pair_ids` to select a subset of callers/samples. Each selected callset must be completed.

The generated `input_manifest.tsv` and `signature-inputs.json` are saved in `logs/signatures-TIMESTAMP/`. A content-addressed manifest under `.nextflow/signature-manifests/` keeps its path and modification time stable for unchanged inputs. The pipeline verifies VCF hashes and primary-contig lengths during conversion. Groups run independently by cohort and caller.

## Install Once

~~~bash
python3.11 -m venv sv_signatures/.venv
sv_signatures/.venv/bin/python -m pip install -r sv_signatures/requirements.lock.txt
~~~

COSMIC does not require R, ClusterSV, copy-number inputs or replication-timing reference construction.

You can verify input files and the environment executable with `python scripts/workflow.py signatures --config config/private-signatures.yaml --check-files`. The same check runs automatically at launch.

## External VCFs (Optional)

For VCFs produced outside this framework, use the five-column [inventory template](../examples/callsets.tsv):

~~~bash
mkdir -p inputs manifests
cp examples/callsets.tsv inputs/callsets.tsv
# Fill the inventory with your external VCF paths, then:
python scripts/build_sv_manifest.py --inventory inputs/callsets.tsv \
  --output manifests/sv_manifest.tsv
~~~

Set this in your signature YAML:

~~~yaml
manifest: manifests/sv_manifest.tsv
~~~

An explicit `manifest` selects these inputs instead of reading `caller_config`. This also supports a combined inventory from multiple cohorts. Leave `manifest: null` for automatic caller-result discovery.

## Launch COSMIC SV32

From the repository root with your site's driver allocation:

~~~bash
mkdir -p logs
sbatch --account=YOUR_ACCOUNT --partition=YOUR_PARTITION \
  --cpus-per-task=1 --mem=4G --time=24:00:00 --output=logs/cosmic-%j.log \
  scripts/run.slurm.sh signatures --config config/private-signatures.yaml --set method=cosmic_sv32
~~~

~~~text
Somatic VCFs → BEDPE → SV clustering → 32-channel counts → COSMIC SV1–SV12 fit
~~~

The workflow converts retained VCF events to BEDPE, applies SigProfiler clustering, builds SV32 counts and refits COSMIC v3.6 reference signatures. The intrachromosomal minimum is 1 kb. The locked pandas semantics and numerical clustering control are checked automatically.

## Retrieve Results

Results are grouped under `OUTDIR/COHORT/CALLER/group_output/`:

| File | Contents |
| --- | --- |
| `SV32.matrix.tsv` | 32-channel feature counts per sample |
| `SV32.exposures.tsv` | Fitted reference-signature exposures |
| Conversion/source audits | Retained events, source hashes and input decisions |
| Reconstruction QC | Fit quality and residuals |

Set `skip_fit: true` for counts only. Scientific/resource settings and overrides are listed in [the parameter guide](../CONFIG_PARAMETERS.md#signature-settings).

## PCAWG Status

The `pcawg_clustersv` route, its CN-free feature/NMF extension and associated numerical workflows are paused. They are excluded from active launch instructions and CI execution. Source, settings, tests and earlier outputs remain available; no files are deleted. We will revisit the method design before enabling it again.
