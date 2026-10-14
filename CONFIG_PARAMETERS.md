# Configuration Parameters

Start with the [caller template](examples/cohort.yaml) or [signature template](examples/signatures.yaml). Copy it to a private YAML file; do not edit versioned defaults for each cohort.

## How Settings Are Resolved

```text
Versioned defaults → YAML overlays in supplied order → --set overrides
```

```bash
python scripts/workflow.py callers --config config/private-cohort.yaml \
  --set delly.options.map_qual=2 --set max_parallel_callers=8 --dry-run
```

`--set` accepts dotted keys and YAML values. Unknown/duplicate keys and malformed mappings are rejected. Dry-run validates configuration, not BAM contents. Actual runs retain the resolved settings.

## Inputs and Execution Settings

| Parameter | Default / meaning |
| --- | --- |
| `input` | Example TSV; replace with your private paired samplesheet |
| `reference` | Placeholder; replace with the exact alignment FASTA |
| `genome_build` | `GRCh37`; checked against primary reference chromosome lengths |
| `exclude_bed` | `null`; optional validated 0-based BED, supported for Manta/DELLY/GRIDSS |
| `dataset_id` | Placeholder label; required for caller benchmark identity |
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

## Caller Resources

| Caller | CPU key / default | Memory key / default | Task time key / default |
| --- | --- | --- | --- |
| Manta | `manta_cpus`: 8 | `manta_memory_gb`: 64 | `manta_time`: `24h` |
| DELLY | `delly.cpus`: 8 | `delly.memory_gb`: 64 | `delly.time`: `24h` |
| SvABA | `svaba.cpus`: 8 | `svaba.memory_gb`: 64 | `svaba.time`: `24h` |
| GRIDSS | `gridss.cpus`: 8 | `gridss.memory_gb`: 96 | `gridss.time`: `24h` |
| BRASS | `brass.cpus`: 8 | `brass.memory_gb`: 64 | `brass.time`: `24h` |

These are starting allocations, not universal measured requirements. Manta's workflow memory defaults to 56 GB and must be below its allocation. GRIDSS heaps default to 40 GB + 4 GB and require allocation headroom. Child-task requests differ from the small driver allocation.

## Scientific Parameters

Only reviewed options are exposed. Others retain the pinned caller's native defaults and help/config records; adding an option requires a tested adapter change.

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

These are examples, not recommended thresholds. Changes require revalidation and a new output directory. The GRIDSS baseline has no panel of normals or BSgenome; raw joint calls are not equivalent to filtered somatic calls.

## Optional Region Exclusions

Set `exclude_bed` to an absolute, reference-matched 0-based half-open BED path, or leave it `null` to preserve the original unmasked baseline. Intervals are checked against FASTA contig names/lengths and merged; aliases are not inferred. The source and normalized region hashes are recorded.

- DELLY receives the normalized exclusions via `-x`.
- GRIDSS receives them via `--blacklist`.
- Manta receives a **complement inclusion list** via `--callRegions`, generated with `bgzip` and `tabix`. These two tools must be available when exclusions are enabled for Manta.

Do not enable this shared option for SvABA or BRASS: the launcher rejects unsupported combinations rather than silently ignoring the mask. Run those callers separately without it. A mask changes the scientific configuration; revalidate in a fresh output directory. Caller-specific interpretation may differ, so this is not a guarantee of identical callable regions across callers.

## Signature Settings

Use **a separate YAML**, based on [config/signatures.yaml](config/signatures.yaml).

| Parameter | Default / meaning |
| --- | --- |
| `manifest` | Placeholder; six-column hashed VCF manifest |
| `genome_build` | `GRCh37`; both supported analysis routes reject other assemblies |
| `python` | `sv_signatures/.venv/bin/python`; separate pinned environment |
| `method` | `cosmic_sv32` or `pcawg_clustersv` |
| `outdir` | Empty: method-specific default; use an explicit new path for experiments |
| `groups` | `all`; optionally comma-separated `cohort:caller` names |
| `pair_id`, `pair_ids` | Empty; single-pair smoke or multi-pair pilot, not both |
| `min_size_bp` | `1000`; COSMIC intrachromosomal minimum |
| `cosmic_version` | `3.6`; SV1–SV12 reference fitting |
| `skip_fit` | `false`; true constructs SV32 counts only |
| `pcawg_stage` | `clustering`; supported CN-free path |
| `pcawg_min_size_bp` | `1`; preserve short oriented junctions |
| `pcawg_fdr_first`, `pcawg_fdr_final`, `pcawg_footprint_fdr` | `0.01`, `0.05`, `0.01` |
| `pcawg_mask_policy` | `local_cen_telo_exclusion`; local endpoint filtering, or `primary_only` |
| `pcawg_timing_source` | `reference`; verified fixed genome-wide tertiles |
| `account`, `partition` | `null`; signature scheduler keys use these names |
| `group_cpus`, `group_memory_gb`, `group_time` | `4`, `32`, `4h` |
| `max_parallel_groups` | `4` simultaneous cohort/caller groups |

Use `pcawg_stage: clustering` for the ClusterSV route. It reports clustering and annotations, not inferred copy-number-dependent mechanisms. See [methods](docs/signatures.md).

## Replication-Timing Reference

Timing reference setup uses [config/replication.yaml](config/replication.yaml): 1-kb bins, full per-track bin coverage and fixed genome-wide tertiles. Copy [examples/timing.yaml](examples/timing.yaml) to `config/private-timing.yaml`, then edit the site settings. It requires a separate environment and is independent of cohort/caller.

See [readiness](docs/workflow-status.md) for the five callers and the two analysis routes.
