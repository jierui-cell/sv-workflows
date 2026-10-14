# Configuration Parameters

Start with the [caller template](examples/cohort.yaml) or [signature template](examples/signatures.yaml). Copy it to a private YAML file and specify the values you want to change; versioned defaults fill in the rest.

Resource paths are already defined for the [standard layout](README.md#standard-folder-layout). Relative YAML and BAM paths are repository-root-relative; absolute overrides remain supported. The image downloader writes exactly the filenames in the defaults.

## How Settings Are Resolved

```text
Versioned defaults → YAML overlays in supplied order → --set overrides
```

```bash
python scripts/workflow.py callers --config config/private-cohort.yaml \
  --set delly.options.map_qual=2 --set max_parallel_callers=8 --dry-run
```

This example is a **preview only**: `--dry-run` prints the final configuration and command, without submitting jobs. `delly.options.map_qual=2` changes DELLY's minimum read mapping quality to 2; `max_parallel_callers=8` permits eight simultaneous pairs per caller, subject to scheduler limits.

`--set` accepts dotted keys and YAML values. Overrides affect that invocation and do not edit your YAML files. For a persistent change, update your private YAML and submit it as usual. The launcher saves resolved settings for each real run; BAM checks happen automatically on compute nodes. [Preview, change and resume examples](README.md#3-optional-preview-or-change-settings).

## Inputs and Execution Settings

| Parameter | Default / meaning |
| --- | --- |
| `input` | `inputs/samplesheet.tsv`; your paired samplesheet |
| `reference` | `resources/reference/reference.fa`; place the alignment FASTA and matching indexes here |
| `genome_build` | `GRCh37`; checked against primary reference chromosome lengths |
| `sequencing_type` | `WGS`; `WES` / `targeted` use the implemented Manta capture route |
| `target_bed` | `null`; optional WES / required targeted inclusion BED |
| `exclude_bed` | `null`; optional validated 0-based BED, supported for Manta/DELLY/GRIDSS |
| `dataset_id` | `my_cohort`; set a distinct cohort label for benchmark identity |
| `outdir` | `results/my_run`; choose a new directory for each experiment |
| `callers` | `preflight`; comma-separated caller names without spaces |
| `pair_ids` | Empty: all pairs; optionally `PAIR_001,PAIR_002` |
| `container_engine` | `apptainer`, or `singularity` |
| `slurm_account`, `slurm_partition` | `null`; set site values when required |
| `scheduler_options` | Empty; additional site-owned scheduler flags |
| `modules_samtools`, `modules_bcftools`, `modules_calling` | `null`; optionally supply site module names, colon-separated for multiple modules |
| `max_parallel_callers` | `3` per caller process, not a global total |
| `max_parallel_preflight` | `3` simultaneous input checks |
| `compact_publish` | `false`; preserve native output hierarchy |

Choose the backend with `--profile`; use `--site-config` for scheduler-specific Nextflow settings. The Slurm driver fixes the Slurm profile. `--resume SESSION_UUID` resumes an existing session. Local real-data execution requires `--allow-local` on an already allocated node/workstation.

## Resource Paths and File Checks

| Resource | Default location |
| --- | --- |
| Manta image | `containers/images/manta-1.6.0--py27h9948957_6.sif` |
| DELLY image | `containers/images/delly-2.6.0--h3752d28_0.sif` |
| SvABA image | `containers/images/svaba-1.2.0--h69ac913_1.sif` |
| GRIDSS image | `containers/images/gridss-2.13.2--h96c455f_6.sif` |
| BRASS image | `containers/images/brass-v6.3.4.sif` |
| BRASS bundle | `resources/brass/bundle/` |

Run `python scripts/workflow.py callers --config config/private-cohort.yaml --check-files` to inspect input/resource presence without submission. It requires only selected callers' files and selected pairs' BAMs/indexes. `--dry-run` instead previews YAML/commands; the flags serve different purposes.

Real launches run file checks automatically before Nextflow starts, and save `setup-check.json` with the run logs. BAM structure/reference agreement is checked on compute nodes. The launcher saves an absolute-path samplesheet snapshot and uses a stable, content-addressed input copy for resume; your original samplesheet is unchanged.

## Caller Resources

| Caller | CPU key / default | Memory key / default | Task time key / default |
| --- | --- | --- | --- |
| Manta | `manta_cpus`: 8 | `manta_memory_gb`: 64 | `manta_time`: `24h` |
| DELLY | `delly.cpus`: 8 | `delly.memory_gb`: 64 | `delly.time`: `24h` |
| SvABA | `svaba.cpus`: 8 | `svaba.memory_gb`: 64 | `svaba.time`: `24h` |
| GRIDSS | `gridss.cpus`: 8 | `gridss.memory_gb`: 96 | `gridss.time`: `24h` |
| BRASS | `brass.cpus`: 8 | `brass.memory_gb`: 64 | `brass.time`: `24h` |

These allocations are configurable per caller. Manta's workflow memory is 56 GB within its 64-GB allocation. GRIDSS uses 40-GB and 4-GB heaps within its 96-GB allocation. The small driver schedules these separate compute tasks.

## Scientific Parameters

Exposed options are typed and recorded. Other options retain the pinned caller's native defaults, captured with its help/config records.

| Caller | Exposed examples and defaults |
| --- | --- |
| Manta | `manta_options.minCandidateVariantSize`: 8; `minSomaticScore`: 10; `minPassSomaticScore`: 30 |
| DELLY | `delly.options`: `map_qual` 1, `minclip` 25, `mad_cutoff` 9; `filter_options`: `quality` 300, `altaf` 0.03, `coverage` 10, `controlcont` 0 |
| SvABA | `svaba.options`: `mate_lookup_min` 3, `chunk_size` 25000, `max_reads` 50000, `max_coverage` 500 |
| GRIDSS | Native joint-call defaults; `somatic_filter: true`; explicit `somatic_filter_options` in the defaults file |
| BRASS | `min_group`: 2; `min_cn_change`: 0.3; `missing_statistics_policy`: `native_defaults`; assembly `GRCh37` |

See all values in [config/defaults.yaml](config/defaults.yaml). Example tuning overlay:

```yaml
manta_options:
  minPassSomaticScore: 40
delly:
  options:
    map_qual: 2
  filter_options:
    altaf: 0.05
```

Use a new output directory for each tuning experiment. GRIDSS's baseline has no panel of normals or BSgenome; use its filtered somatic output for downstream somatic analysis.

## Optional Region Exclusions

Set `exclude_bed` to an absolute, reference-matched 0-based half-open BED path, or leave it `null` to preserve the original unmasked baseline. Intervals are checked against FASTA contig names/lengths and merged; aliases are not inferred. The source and normalized region hashes are recorded.

- DELLY receives the normalized exclusions via `-x`.
- GRIDSS receives them via `--blacklist`.
- Manta receives a **complement inclusion list** via `--callRegions`, generated with `bgzip` and `tabix`. These two tools must be available when exclusions are enabled for Manta.

For SvABA or BRASS, use a separate run without this shared setting. Exclusions are a scientific parameter recorded with the outputs; their native interpretation differs by caller.

## Capture-Assay Settings

For WES or a targeted panel, copy your completed cohort YAML to a private capture configuration, retain its reference/image/HPC settings, and update the input samplesheet and dataset label. Set `callers: manta` and choose `sequencing_type: WES` or `targeted`.

`target_bed` is the capture kit/panel's genomic interval list. WES can use `target_bed: null`; targeted panels require a BED path. BED starts are 0-based and ends excluded; contig names must match the reference. Manta automatically receives `--exome`. For a supplied BED, the workflow merges intervals, subtracts exclusions and generates indexed `--callRegions`, using compute-node `bgzip`/`tabix`.

Submit that private configuration with the same Slurm driver used for WGS. This runs Manta SV calling, not signatures; COSMIC SV32 remains a WGS-only route. [Complete WES/panel example](README.md#exome-and-targeted-sequencing).

## Signature Settings

Use **a separate YAML**, based on [config/signatures.yaml](config/signatures.yaml).

| Parameter | Default / meaning |
| --- | --- |
| `manifest` | `manifests/sv_manifest.tsv`; generated six-column hashed VCF manifest |
| `genome_build` | `GRCh37`; active COSMIC SV32 route uses this assembly |
| `sequencing_type` | `WGS`; genome-wide signature routes use WGS callsets |
| `python` | `sv_signatures/.venv/bin/python`; separate pinned environment |
| `method` | `cosmic_sv32`; PCAWG/ClusterSV is on hold |
| `outdir` | Empty: method-specific default; use an explicit new path for experiments |
| `groups` | `all`; optionally comma-separated `cohort:caller` names |
| `pair_id`, `pair_ids` | Empty: all pairs; optional single-pair or multi-pair selection |
| `min_size_bp` | `1000`; COSMIC intrachromosomal minimum |
| `cosmic_version` | `3.6`; SV1–SV12 reference fitting |
| `skip_fit` | `false`; true constructs SV32 counts only |
| `account`, `partition` | `null`; signature scheduler keys use these names |
| `group_cpus`, `group_memory_gb`, `group_time` | `4`, `32`, `4h` |
| `max_parallel_groups` | `4` simultaneous cohort/caller groups |

PCAWG/ClusterSV signature analysis is on hold, including clustering, feature construction and NMF. Its settings remain in the defaults for future work, but selecting `method=pcawg_clustersv` stops before task submission. COSMIC settings and caller workflows are unchanged. See [status](docs/signatures.md#pcawg-status).

## Replication-Timing Reference

The replication-timing implementation and [configuration](config/replication.yaml) are retained with the paused PCAWG work. No timing-reference construction is needed for COSMIC SV32.

See [method scope](docs/workflow-status.md) for assay/build support and output interpretation.
